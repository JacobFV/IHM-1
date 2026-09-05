"""full error covariance, cross-format identity, source manifests and input errors."""
import json
import tempfile
from pathlib import Path
import numpy as np
from ihm import body, materialize, Request
from ihm.forge.ingest import Evidence, SourceCard, load

r = body(); request = Request(('metabolic.glucose',), 's')
m = materialize(r, request); i = m.components.index('metabolic.glucose')
mean, variance = m.mean[i], m.cov[i,i]
evidence = [Evidence('lab', f'row-{j}', 's', 0, 'metabolic.glucose', 7., 1., 'measured') for j in range(2)]
noise = np.array([[1., .9], [.9, 1.]])
m.assimilate(evidence, error_covariance=noise)
precision = np.ones(2) @ np.linalg.solve(noise, np.ones(2))
assert np.isclose(m.cov[i,i], 1/(1/variance+precision))
independent = materialize(r, request); independent.assimilate(evidence)
assert m.cov[i,i] > independent.cov[i,i]
# Bad batches are atomic even when they are at a future timestamp.
bad = Evidence('lab','future','s',10,'metabolic.glucose',7.,1.,'measured')
before = (m.mean.copy(),m.cov.copy(),m.time)
try:
    m.assimilate([bad], error_covariance=np.array([[-1.]]))
except ValueError:
    pass
else:
    raise AssertionError('invalid observation covariance accepted')
assert np.array_equal(before[0],m.mean) and np.array_equal(before[1],m.cov) and before[2]==m.time
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    row = dict(id='record', subject='s', time=0., component='metabolic.glucose', value=90., unit='mg/dL', variance=81.)
    card = SourceCard('multiformat','synthetic',('metabolic.glucose',),'synthetic verification')
    js = root/'data.json'; js.write_text(json.dumps([row]))
    npz = root/'data.npz'; np.savez(npz, **{k: np.array([v]) for k,v in row.items()})
    first, second = load(js,card,r), load(npz,card,r)
    assert first[0].value == second[0].value and first[0].variance == second[0].variance
    assert first[0].key == second[0].key
    m = materialize(r,request); m.assimilate(first)
    try:
        m.assimilate(second)
    except ValueError:
        pass
    else:
        raise AssertionError('format conversion defeated deduplication')
    for column, invalid in (('unit','unknown'),('variance',0),('value',float('nan')),('component','missing')):
        js.write_text(json.dumps([{**row,column:invalid}]))
        try:
            load(js,card,r)
        except ValueError:
            pass
        else:
            raise AssertionError(f'invalid {column} accepted')
    js.write_text(json.dumps([{**row,'value':None}]))
    assert load(js,card,r) == []
try:
    r.anatomy['organ_partition'].memberships['heart']['heart'] = -1
except TypeError:
    pass
else:
    raise AssertionError('sealed anatomy is mutable')
print('verified: correlated likelihood, atomic rejection, JSON/NPZ equivalence, deduplication, bad inputs, missing values, immutable anatomy')
