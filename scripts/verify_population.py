"""verify source-derived population beliefs and preserve the dynamics boundary."""
import json
from pathlib import Path
import tempfile
import numpy as np
from ihm import body,materialize,Request
from ihm.materialize.model import Model

path=Path('data/derived/population/nhanes-2017-2018/joint-population-prior.json')
prior=json.loads(path.read_text());r=body();request=Request(('metabolic.glucose','blood.hemoglobin'),'population-reference')
a=materialize(r,request);b=materialize(r,request,population_prior=path)
assert set(prior['fit_subject_ids']).isdisjoint(prior['holdout_subject_ids'])
assert len(prior['fit_subject_ids'])==prior['n_fit'] and prior['n_fit']>1000
for c in set(prior['components']) & set(b.components):
    i=b.components.index(c);j=prior['components'].index(c)
    assert b.mean[i]==prior['mean'][j]
    assert b.cov[i,i]==prior['covariance'][j][j]
assert np.array_equal(a.a,b.a) and np.array_equal(a.b,b.b) and np.array_equal(a.q,b.q)
assert b.provenance['population_initialization']['parameter_calibration'] is False
with tempfile.TemporaryDirectory() as root:
    target=Path(root)/'population.npz';b.save(target);loaded=Model.load(target)
    assert np.array_equal(b.mean,loaded.mean)
    bad=Path(root)/'bad.json';wrong=dict(prior);wrong['units']=['wrong']*len(prior['units']);bad.write_text(json.dumps(wrong))
    try: materialize(r,request,population_prior=bad)
    except ValueError: pass
    else: raise AssertionError('unit-mismatched population prior accepted')
print('verified: measured population initialization, fit/holdout separation, serialization, units; dynamics unchanged')

import pandas as pd
root=Path('data/derived/population/nhanes-2017-2018')
glu=pd.read_parquet(root/'GLU_J.parquet');demo=pd.read_parquet(root/'DEMO_J.parquet')
assert int((glu.WTSAF2YR==0).sum())==325
assert int((demo.WTMEC2YR==0).sum())==550
assert not np.any(glu.WTSAF2YR.to_numpy()==2.**-260)
linked=pd.read_parquet(root/'adult-linked-measurements.parquet')
assert set(linked['endocrine.insulin__detection_flag'].dropna().unique()) <= {0.,1.}
assert prior['detection_limit_handling']['fit_below_detection_counts']
for definition in b.provenance['components']:
    c=definition['id']
    if c in prior['components']:
        i=b.components.index(c)
        assert definition['mean']==b.mean[i]
        assert definition['std']==np.sqrt(b.cov[i,i])
        assert 'original_declared_prior' in definition
print('verified: XPORT zero weights, detection flags and initialized component metadata')
