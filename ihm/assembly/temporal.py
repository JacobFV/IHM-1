"""Finite-window spectra of one canonical body's actual simulated variables."""
import numpy as np
from ihm.temporal import finite_laplace, spectral_estimate
from .body import read_native,digest,write_json
from pathlib import Path

def build_body_spectra(native_directory,output,body_payload):
    directory=Path(native_directory);series=read_native(directory);t=np.array(series['time_s']);dt=float(np.median(np.diff(t)))
    if len(t)<4 or not np.allclose(np.diff(t),dt,rtol=1e-6,atol=1e-9):raise ValueError('Spectral density requires at least four uniformly spaced native samples')
    frequency=np.linspace(0,min(5.,.5/dt),121);sigma=[0.,1/(t[-1]-t[0]),.1,1.]
    grid=np.array([s+2j*np.pi*f for s in sigma for f in frequency]);variables=[]
    selected=[('physiology',name) for name in ['ArterialPressure(mmHg)','TotalLungVolume(mL)','BloodVolume(mL)','CoreTemperature(degC)','OxygenSaturation']]
    selected += [('compartments',name) for name in ['LeftHeartVolume(mL)','RightHeartVolume(mL)','LeftLungPulmonaryGasVolume(mL)','RightLungPulmonaryGasVolume(mL)','LymphVolume(mL)','SkinTissueExtracellularVolume(mL)']]
    for group,name in selected:
        x=np.array([r[name] for r in series[group]]);mean=float(x.mean());laplace=finite_laplace(t,x-mean,grid).reshape(len(sigma),len(frequency))
        variables.append({'name':name,'source_group':group,'removed_sample_mean':mean,'psd':spectral_estimate(x,1/dt,nperseg=min(len(t),512)),
                          'laplace_real':laplace.real.tolist(),'laplace_imag':laplace.imag.tolist()})
    if series['summary'].get('input_state_sha256') is not None or series['summary']['configuration']['patient']!=body_payload['profile']['native_patient']:raise ValueError('Canonical spectra require fresh initialization of the shared patient')
    if series['summary']['patient_sha256']!=body_payload['profile']['native_patient_sha256']:raise ValueError('Spectra patient differs from canonical body')
    result={'model_id':'ihm-body','canonical_sources':body_payload['sources'],'runtime_sources':body_payload['runtime_sources'],'native_directory':str(directory.resolve()),'native_summary_sha256':series['input_hashes']['summary.json'],'source_files':{name:series['input_hashes'][name] for name in ['native_multisystem.csv','body_compartments.csv']},'time_interval_s':[float(t[0]),float(t[-1])],
            'sample_interval_s':dt,'laplace':{'sigma_per_s':sigma,'frequency_hz':frequency.tolist(),'s_convention':'sigma + 2 pi i frequency_hz','detrend':'subtract full-record sample mean','integration':'trapezoidal; elapsed time from first sample','units':'signal unit times seconds'},
            'variables':variables,'interpretation':'Finite-window descriptors, not physiological poles, calibrated predictor transfer functions or evidence of causal coupling.'}
    write_json(output,result);return result
