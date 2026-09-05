"""Evidence-bound display reduction of continuing native body experiments."""
import json
from pathlib import Path
import numpy as np
from ihm.native import _sha
from ihm.temporal import finite_laplace, spectral_estimate
from .respiration import BodyRespiration

ALIASES={'heart_rate_per_min':'HeartRate(1/min)', 'mean_arterial_pressure_mmhg':'MeanArterialPressure(mmHg)',
         'lung_volume_ml':'TotalLungVolume(mL)', 'arterial_co2_mmhg':'ArterialCarbonDioxidePressure(mmHg)',
         'oxygen_saturation':'OxygenSaturation','core_temperature_c':'CoreTemperature(degC)'}
SELECTED=('heart_rate_per_min','arterial_pressure_mmhg','lung_volume_ml','respiratory_request_cmh2o',
          'respiratory_applied_cmh2o','airway_flow_l_per_s','arterial_co2_mmhg','arterial_o2_mmhg',
          'stomach_carbohydrate_g','Aorta.Glucose.concentration_mg_per_dl',
          'insulin_synthesis_pmol_per_min','liver_glycogen_g','metabolic_rate_w')


def project_systemic(root, path):
    root, path=Path(root).resolve(),Path(path).resolve()
    data=json.loads(path.read_text())
    for source,expected in data['runtime_sources'].items():
        if _sha(root/source)!=expected:raise ValueError('Systemic source changed: '+source)
    t=np.array([frame['time_s'] for frame in data['frames']])
    dt=float(np.median(np.diff(t))) if len(t)>1 else float('nan')
    if len(t)<4 or not np.isfinite(t).all() or dt<=0 or not np.allclose(np.diff(t),dt,rtol=1e-8,atol=1e-8):
        raise ValueError('Projection requires a uniform increasing native record with at least four frames')
    anatomy_path=root/'data/derived/canonical/anatomy.json'
    respiration_path=root/'data/derived/canonical/respiration.json'
    anatomy=json.loads(anatomy_path.read_text())
    model=json.loads(respiration_path.read_text())
    if model['anatomy_sha256']!=_sha(anatomy_path):raise ValueError('Thorax anatomical source changed')
    # Coarse chemical observations cannot resolve respiration. Never interpolate a
    # breathing cycle into a six-hour nutrient record sampled every thirty seconds.
    dense=dt<=.1+1e-9
    body=BodyRespiration(model) if dense else None
    reference=data['frames'][0]['values']['lung_volume_ml']
    if reference is None or reference<=0:raise ValueError('Missing initial lung volume')
    frames=[]
    max_energy_error=0.
    for i,frame in enumerate(data['frames']):
        values=frame['values']
        rendered={'time_s':frame['time_s'],'entities':{},'systemic_values':values,
                  'physiology':{alias:values[key] for key,alias in ALIASES.items()}}
        if dense:
            volume=values['lung_volume_ml']
            if volume is None or not np.isfinite(volume) or volume<=0:raise ValueError('Missing native lung volume')
            if i:
                state=body.step(float(t[i]-t[i-1]), {'volume_change_m3':(volume-reference)*1e-6,
                   'lung_volume_ratios':{key:volume/reference for key in body.lung_ids}})
                rendered['entities']=state['entities']
                rendered['respiration']=state
                max_energy_error=max(max_energy_error,abs(state['audit']['energy_balance_residual_j']))
            else:
                rendered['respiration']={'skin_field':{**model['skin_field'],'displacement_m':[0,0,0]},
                                         'time_s':0,'displacement_m':[0,0,0],'diaphragm_descent_m':0}
        frames.append(rendered)
    duration=float(t[-1]-t[0])
    frequency=np.linspace(0,min(5.,.5/dt),101)
    sigma=np.array([0.,1/duration,4/duration])
    grid=np.array([s+2j*np.pi*f for s in sigma for f in frequency])
    variables=[]
    for name in SELECTED:
        if any(frame['values'].get(name) is None for frame in data['frames']):continue
        x=np.array([frame['values'][name] for frame in data['frames']])
        mean=float(x.mean());laplace=finite_laplace(t,x-mean,grid).reshape(len(sigma),len(frequency))
        variables.append(dict(name=name,unit=data['fields'][name]['unit'],removed_sample_mean=mean,
            psd=spectral_estimate(x,1/dt,nperseg=min(len(x),512)),
            laplace_real=laplace.real.tolist(),laplace_imag=laplace.imag.tolist()))
    sources={**data['runtime_sources'],str(path.relative_to(root)):_sha(path)}
    for source in ('ihm/assembly/systemic_projection.py','ihm/assembly/respiration.py','ihm/temporal/spectra.py',
                   'data/derived/canonical/anatomy.json','data/derived/canonical/respiration.json'):
        sources[source]=_sha(root/source)
    return dict(schema='ihm.systemic-display.v1',model_id='ihm-body',has_body_projection=dense,
                configuration=data['configuration'],frames=frames,fields=data['fields'],
                centroids_m={entity['id']:entity['centroid_m'] for entity in anatomy['entities']},
                clock=dict(start_s=float(t[0]),end_s=float(t[-1]),sample_interval_s=dt,native_step_s=.02,
                           coupling_mode='native observation → reduced thoracic projection' if dense else 'chemical observation; geometry at reference',
                           playback='finite computed record; no interpolated or invented cycles'),
                mechanism_edges=data['mechanism_edges'],actions=data['actions'],native_manifest=data['native_manifest'],
                runtime_sources=sources,summary=data['summary'],checks=data['checks'],
                spectra=dict(frequency_hz=frequency.tolist(),sigma_per_s=sigma.tolist(),variables=variables,
                             convention='s = sigma + 2 pi i f; finite trapezoidal integral of mean-subtracted signal',
                             unit='signal unit times seconds',interpretation='Finite-window descriptors, not physiological poles or identified causal transfer functions'),
                projection_audit=dict(maximum_energy_balance_residual_j=max_energy_error if dense else None,
                    lung_partition='equal fractional total-volume change assigned to all five lobes; uncalibrated projection prior' if dense else None,
                    whole_body_mechanical_feedback=False),
                limitations=data['limitations']+[
                    'Thoracic projection uses generic mode fractions and common fractional lobe expansion; individual lobe gas volumes are not resolved here.',
                    'Only the thoracic modes and skin field move in dense recordings. Remaining geometry stays at its reference pose.',
                    'Sparse nutrient records cannot resolve cardiac or respiratory cycles. Their Nyquist frequency is explicitly limited by output cadence.'])
