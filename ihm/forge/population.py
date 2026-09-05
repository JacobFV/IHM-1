"""initialize shared state beliefs from measured population covariance.

Population dispersion is uncertainty about an unspecified individual, not assay
noise or a calibration of dynamical coefficients. No dynamic parameter is changed.
"""
import json
from pathlib import Path
import numpy as np
from ihm.forge.acquisition import sha256
from ihm.runtime.gaussian import covariance


def initialize_population(model, path):
    if model.time != 0 or model.seen:
        raise ValueError('population initialization requires a fresh unconditioned model')
    path=Path(path); prior=json.loads(path.read_text())
    if prior.get('basis')!='survey_weighted_empirical_population_state_prior' or prior.get('schema_version')!=1:
        raise ValueError('unsupported population prior schema')
    components=prior['components']; units=prior['units']; mean=np.asarray(prior['mean'],float)
    if len(set(components))!=len(components) or len(units)!=len(components) or mean.shape!=(len(components),) or not np.isfinite(mean).all():
        raise ValueError('invalid population prior components, units or mean')
    cov=covariance(prior['covariance'],len(components))
    if set(prior['fit_subject_ids']) & set(prior['holdout_subject_ids']):
        raise ValueError('population fitting and holdout subjects overlap')
    definitions={c['id']:c for c in model.provenance['components']}
    selected=[c for c in components if c in model.components]
    if not selected: raise ValueError('population prior does not overlap this materialization')
    model_ids=[model.components.index(c) for c in selected]; prior_ids=[components.index(c) for c in selected]
    for c,j in zip(selected,prior_ids):
        if definitions[c]['unit']!=units[j]: raise ValueError(f'population prior unit mismatch: {c}')
    new_mean=model.mean.copy(); new_cov=model.cov.copy()
    new_mean[model_ids]=mean[prior_ids]
    # The prior supplies no covariance to components it did not measure.
    new_cov[model_ids,:]=0.;new_cov[:,model_ids]=0.
    new_cov[np.ix_(model_ids,model_ids)]=cov[np.ix_(prior_ids,prior_ids)]
    covariance(new_cov,len(new_mean))
    model.mean=new_mean;model.cov=new_cov
    for c,i in zip(selected,model_ids):
        definition=definitions[c]
        definition['original_declared_prior']={k:definition[k] for k in ('mean','std','provenance')}
        definition['mean']=float(new_mean[i]);definition['std']=float(np.sqrt(new_cov[i,i]))
        definition['provenance']='measured_population_prior:'+prior['source']
    model.provenance['population_initialization']={
        'source':prior['source'],'path':str(path),'sha256':sha256(path),'population':prior['population'],
        'components':selected,'n_fit':prior['n_fit'],'n_internal_holdout':prior['n_internal_holdout'],
        'basis':prior['basis'],'weight':prior['weight'],'bindings':prior['bindings'],
        'parameter_calibration':False,'limitations':prior['limitations'],
        'detection_limit_handling':prior.get('detection_limit_handling'),
        'decoder':prior.get('decoder'),
        'forward_dynamics':'unchanged; population fitting does not validate process coefficients'}
    return model
