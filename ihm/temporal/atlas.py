"""Build source-separated temporal evidence; no cross-subject causal fusion."""
from pathlib import Path
import csv
import hashlib
import json
import re
import numpy as np
from . import spectral_estimate, cross_spectrum, finite_laplace, fit_predictor

LIMITATIONS=[
    'Native BioGears traces are source simulations, not human measurements or independent validation.',
    'Coherence measures linear association and does not identify causal direction or intervention response.',
    'The native 60-second baseline cannot identify slow endocrine, immune, renal or thermoregulatory dynamics.',
    'Hann segments remove their means; low-frequency drift, DC and finite-horizon leakage require separate interpretation.',
    'Reduced forecasts are autonomous empirical approximations; holdout error can exceed persistence and no clinical validity is claimed.',
    'Measured BIDMC record 01 is one critically ill hospital patient, not matched to the simulated patient; no direct amplitude calibration is performed.'
]


def _read(path, sampling_rate_hz=None):
    with path.open() as f:
        reader=csv.reader(f); header=[v.strip() for v in next(reader)]
        a=np.array([[float(v) for v in row] for row in reader if row])
    if a.ndim!=2 or a.shape[1]!=len(header) or not np.isfinite(a).all():
        raise ValueError('CSV must contain a rectangular finite numeric table')
    t=a[:,0]
    if sampling_rate_hz is not None:
        dt=1/float(sampling_rate_hz)
        regular=t[0]+np.arange(len(t))*dt
        if not np.isfinite(dt) or dt<=0 or np.max(abs(t-regular))>.00500001:
            raise ValueError('CSV timestamps disagree with declared source-header cadence')
        t=regular
    dt=np.median(np.diff(t)) if sampling_rate_hz is None else dt
    if dt<=0 or not np.allclose(np.diff(t),dt,rtol=1e-5,atol=1e-8):
        raise ValueError('CSV time must be increasing and regularly sampled in seconds')
    return header,t,a[:,1:],float(dt)


def _run(path, root, run_id, source_kind, units=None, sampling_rate_hz=None, *, slow=False, run_notes=None):
    header,t,y,dt=_read(path,sampling_rate_hz)
    if 'CTSresistance' in header:
        keep=[i for i,h in enumerate(header[1:]) if h!='CTSresistance']
        y=y[:,keep]; header=[header[0]]+[header[i+1] for i in keep]
    omitted=[]
    if slow:
        omitted=[h for h in header[1:] if h.startswith(('ArterialPressure(', 'TotalLungVolume('))]
        keep=[i for i,h in enumerate(header[1:]) if h not in omitted]
        y=y[:,keep];header=[header[0]]+[header[i+1] for i in keep]
    fs=1/dt; nperseg=min(len(t)//2,int(round((900 if slow else 10)*fs)))
    variables=[]
    for i,name in enumerate(header[1:]):
        match=re.fullmatch(r'(.+)\(([^()]*)\)',name)
        label=match[1] if match else name
        unit=(units or {}).get(name,match[2] if match else 'dimensionless')
        spec=spectral_estimate(y[:,i],fs,nperseg=nperseg)
        std=float(y[:,i].std()); threshold=max(1,abs(float(y[:,i].mean())))*1e-12
        peak=int(np.argmax(spec['psd'][1:]))+1
        variables.append(dict(id=label,label=label,unit=unit,psd_unit=f'({unit})^2/Hz',
            mean=float(y[:,i].mean()),variance=std**2,constant=std<=threshold,
            dominant_frequency_hz=spec['frequency_hz'][peak] if std>threshold else None,**spec))
    # Full synchronized pairs for selected physiological relationships only.
    ids=[v['id'] for v in variables]
    desired=[('ArterialPressure','TotalLungVolume'),('ArterialPressure','CardiacOutput'),
             ('TotalLungVolume','ArterialCarbonDioxidePressure'),('RESP','II'),('RESP','PLETH')]
    pairs=[]
    for xname,yname in desired:
        if xname in ids and yname in ids:
            ix,iy=ids.index(xname),ids.index(yname)
            pairs.append(dict(x=xname,y=yname,unit=f"{variables[ix]['unit']}*{variables[iy]['unit']}/Hz",
                 **cross_spectrum(y[:,ix],y[:,iy],fs,nperseg=nperseg)))
    lf=np.linspace(0,min(.05 if slow else 5,fs/2),101); sigma=[0.,1/3600,1/600,1/60] if slow else [0.,.1,1.]
    laplace=[]
    for i,var in enumerate(variables):
        z=np.array([finite_laplace(t,y[:,i]-y[:,i].mean(),sig+2j*np.pi*lf) for sig in sigma])
        laplace.append(dict(id=var['id'],unit=var['unit']+'*s',real=z.real.tolist(),imag=z.imag.tolist()))
    model=fit_predictor(y,dt,rank=min(8,y.shape[1]))
    n=model.diagnostics['train_samples']; forecast=model.forecast(y[n-1],len(y)-n)
    stride=max(1,int(np.ceil(len(t)/1200))); hi=np.arange(n,len(t),max(1,int(np.ceil((len(t)-n)/600))))
    prediction=dict(model=model.to_dict(),variable_ids=ids,units=[v['unit'] for v in variables],
        holdout_time_s=t[hi].tolist(),holdout_actual=y[hi].T.tolist(),holdout_predicted=forecast[hi-n].T.tolist(),
        first_holdout_time_s=float(t[n]),last_training_time_s=float(t[n-1]),
        display_decimation='stride only for plotted trajectories; fits and spectra use full-rate data')
    return dict(id=run_id,source_kind=source_kind,source=str(path.relative_to(root)),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),sample_rate_hz=fs,dt_s=dt,
        samples=len(t),duration_s=float(t[-1]-t[0]),sample_start_s=float(t[0]),sample_end_s=float(t[-1]),
        time_basis='seconds from source CSV; Laplace uses elapsed time from first sample',
        variables=variables,cross_spectra=pairs,spectral_interpretation='Observed finite-window Fourier PSD, not a causal transfer function',
        omitted_undersampled_channels=omitted,window_duration_s=nperseg/fs,run_notes=run_notes or [],
        trace=dict(time_s=t[::stride].tolist(),variable_ids=ids,values=y[::stride].T.tolist()),
        laplace=dict(sigma_per_s=sigma,frequency_hz=lf.tolist(),detrend='whole-record mean removed',
            quadrature='trapezoid; finite observed horizon only; no tail extrapolation',variables=laplace),
        predictor=prediction,limitations=LIMITATIONS+(['One-Hz output has Nyquist 0.5 Hz; instantaneous arterial pressure and lung volume are excluded because their pulse/breath oscillations can alias. No anti-alias filtering is established for other source output.', 'The 900-second window resolves 1/900 Hz; mean-removed drift is not evidence of stationary oscillation or a physiological cycle.'] if slow else []))


def build_atlas(root, native_csv=None, output=None):
    root=Path(root).resolve(); native_csv=Path(native_csv or root/'data/derived/physiology/biogears_native_run/native_multisystem.csv').resolve()
    runs=[_run(native_csv,root,'biogears_native_baseline','source_simulation')]
    measured=root/'data/raw/temporal/bidmc/bidmc_01_Signals.csv'
    provenance=measured.parent/'provenance.json'
    if measured.exists() and provenance.exists():
        hea=(measured.parent/'bidmc01.hea').read_text().splitlines()
        record=hea[0].split(); cadence=float(record[2]); count=int(record[3])
        units={line.split()[-1].rstrip(','):line.split()[2].split('/')[-1] for line in hea[1:6]}
        measured_run=_run(measured,root,'bidmc_01','human_measurement',units,cadence)
        if measured_run['samples']!=count: raise ValueError('waveform count differs from WFDB header')
        measured_run['time_basis']='125 Hz sample index reconstructed from verified WFDB header; CSV clock rounded to 0.01 s above 100 s, max error 0.004 s'
        measured_run['measurement_notes']='RESP uses header unit pm, an impedance respiratory surrogate, not calibrated lung volume; PLETH uses NU; ECG leads use mV.'
        measured_run['provenance']=json.loads(provenance.read_text())
        runs.append(measured_run)
    for directory in sorted((root/'data/derived/physiology').glob('native_hour_*')):
        path=directory/'native_multisystem.csv'
        if not path.exists():continue
        config_path=directory/'configuration.json'
        config=json.loads(config_path.read_text()) if config_path.exists() else {}
        if config.get('seconds')!=3600 or config.get('sample_hz')!=1:continue
        notes=['Long native source experiment, not a stationary healthy baseline. Cool-room 22 C / 0.5 clo rest cools substantially; initial GI water redistribution contributes to fluid shifts.', 'Thermal comparisons are distinct source conditions; inspect the attached run configuration.']
        run=_run(path,root,directory.name,'source_simulation',slow=True,run_notes=notes)
        run['configuration']=config
        if config_path.exists():run['configuration_sha256']=hashlib.sha256(config_path.read_bytes()).hexdigest()
        runs.append(run)
    causal=[]
    circuit_path=root/'data/derived/coupling/native-skin-circuit.json'
    if circuit_path.exists():
        from ihm.coupling.circuit import FluidCircuit
        artifact=json.loads(circuit_path.read_text());model=FluidCircuit.from_dict(artifact['model'])
        frequency=np.geomspace(1e-7,.1,121);sigma=[0.,1/3600]
        for sig in sigma:
            response=model.response(sig+2j*np.pi*frequency,'Aorta1')
            causal.append(dict(id='native_skin_lymph',kind='causal_frozen_descriptor_resolvent',boundary='Aorta1',sigma_per_s=sig,
                frequency_hz=frequency.tolist(),pressures_pa_per_pa={k:dict(real=v.real.tolist(),imag=v.imag.tolist()) for k,v in response['pressures_pa_per_pa'].items()},
                flows_m3_s_per_pa={k:dict(real=v.real.tolist(),imag=v.imag.tolist()) for k,v in response['flows_m3_s_per_pa'].items()},
                finite_poles_per_s=dict(real=model.finite_poles_per_s().real.tolist(),imag=model.finite_poles_per_s().imag.tolist()),
                source=artifact['model']['source'],artifact_path=str(circuit_path.relative_to(root)),artifact_sha256=hashlib.sha256(circuit_path.read_bytes()).hexdigest(),
                limitations=['Causal small-signal response of source-coefficient hydraulic equations, not observed Fourier power.', 'Native gates and coefficients frozen; external boundary flows fixed. Zero-frequency excluded because the descriptor has an intracellular storage integrator.', 'No empirical identification or clinical validation from short observed trajectories.']))
    inspected=[]
    for name in ['fields/uncertainty/spectral.py','fields/priors.py','forge/spectra.py','processes/base.py','processes/mechanical.py']:
        path=root.parent/'IBM-1/ibm'/name
        if path.exists(): inspected.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    result=dict(schema_version=1,runs=runs,causal_responses=causal,additional_source_trajectories=[dict(path='data/derived/csf/index.json',interpretation='Nonstationary prescribed-pressure responses; use source ODE and local poles rather than interpreting short transients as cycles.'),dict(path='data/derived/reproductive/index.json',interpretation='Nonperiodic prescribed hormonal inputs; the ten-day example cannot establish a menstrual-cycle spectrum.')],limitations=LIMITATIONS,ibm_reference=dict(
        inspected_sources=inspected,interpretation='Cyclic temporal Laplacian eigenvectors are Fourier modes. SpectralGaussian stores coefficient uncertainty; power density needs explicit cadence/window normalization. LTI transfer functions are separately evaluated on i*omega.'))
    out=Path(output or root/'data/derived/temporal/index.json'); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,separators=(',',':'),allow_nan=False))
    for run in runs:
        (out.parent/(run['id']+'.json')).write_text(json.dumps(run,separators=(',',':'),allow_nan=False))
    return result
