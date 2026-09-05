"""Retain finite-window spectra of actual regional body experiment observables."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.temporal import finite_laplace,spectral_estimate
from ihm.app.experiments import read_experiment


def spectrum(id,time,columns,units):
    t=np.array(time);dt=float(np.median(np.diff(t)))
    if not np.allclose(np.diff(t),dt,rtol=1e-9,atol=1e-12):raise ValueError('Regional spectral samples must be uniform')
    frequency=np.linspace(0,.5/dt,121);sigma=[0.,1/(t[-1]-t[0]),1.];variables=[];laplace=[]
    for key,values in columns.items():
        x=np.array(values);unit=units[key]
        z=np.array([finite_laplace(t,x-x.mean(),s+2j*np.pi*frequency) for s in sigma])
        variables.append({'id':key,'label':key,'unit':unit,'psd_unit':unit+'²/Hz',**spectral_estimate(x,1/dt,nperseg=len(x))})
        laplace.append({'id':key,'unit':unit+'·s','real':z.real.tolist(),'imag':z.imag.tolist(),'removed_sample_mean':float(x.mean())})
    return {'id':id,'source_kind':'Recorded regional body experiment','sample_rate_hz':1/dt,'duration_s':float(t[-1]-t[0]),
        'variables':variables,'laplace':{'frequency_hz':frequency.tolist(),'sigma_per_s':sigma,'variables':laplace},
        'limitations':['Single finite record; no confidence interval or inferred physiological poles.','Fourier/Laplace signal descriptors are distinct from the mechanistic IBM transfer or RC resolvent.'],
        'time_basis':'elapsed seconds from first recorded sample; trapezoidal finite Laplace, full-record mean removed, no tail extrapolation'}


def build():
    touch=read_experiment(ROOT,'forearm-touch');electric=read_experiment(ROOT,'skin-electric');runs=[]
    columns={k:[f[k] for f in touch['frames']] for k in ['indentation_m','reaction_n','elastic_energy_j','rapid_response','slow_response']}
    runs.append(spectrum('forearm-touch',[f['time_s'] for f in touch['frames']],columns,dict(zip(columns,['m','N','J','donor response','donor response']))))
    for name,experiment in electric['experiments'].items():
        frames=experiment['frames']
        columns={'Central membrane voltage':[f['membrane_voltage_V'][4] for f in frames],
            'Central apical voltage':[f['apical_voltage_V'][4] for f in frames],
            'Maximum absolute lateral field':[max(abs(v) for v in f['edge_field_V_m']) for f in frames]}
        runs.append(spectrum('skin-electric:'+name,[f['time_s'] for f in frames],columns,dict(zip(columns,['V','V','V/m']))))
    inputs=['data/derived/canonical/forearm-touch.json','data/derived/canonical/skin-electric.json']
    code=['scripts/build_body_experiment_spectra.py','ihm/temporal/spectra.py']
    result={'schema_version':1,'runs':runs,'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        'runtime_sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in code}}
    out=ROOT/'data/derived/canonical/regional-spectra.json';out.write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
    print('Recorded regional Fourier/Laplace descriptors:',len(runs));return result

if __name__=='__main__':build()
