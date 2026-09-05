"""Materialization, covariance conditioning and reproducible evidence binding."""
import json
import tempfile
from pathlib import Path
import numpy as np
from ihm.human import ImplicitHuman, PopulationBelief

p=PopulationBelief(['a','b'],['m','m'],[0,0],[[1,.5],[.5,1]],{'source':'analytic'})
q=p.condition({'a':(2.,1.,'m')})
assert np.allclose(q.mean,[1,.5])
assert np.allclose(q.cov,[[.5,.25],[.25,.875]])
assert np.allclose(p.mean,[0,0])
for data in ({'a':(1,-1,'m')},{'a':(1,1,'mm')},{'c':(1,1,'m')}):
    try:p.condition(data)
    except ValueError:pass
    else:raise AssertionError('invalid observation accepted')
root=Path(__file__).resolve().parents[1];human=ImplicitHuman.open(root)
assert human.describe()['coverage']['declared_components']==180
assert len(human.fields())>500
assert human.materialize('population').mean.shape==(23,)
assert human.materialize('skin-field').predict('18-29','female','LA-1')['mean_V_per_m']>0
c=human.materialize('skin-lymph'); assert c.step(.02)['time_s']==.02
r=human.describe()['temporal_runs'][0];m=human.materialize('temporal',run_id=r)
assert np.isfinite(m.forecast(m.center,5)).all()
with tempfile.TemporaryDirectory() as td:
    path=Path(td)/'human.json';human.save(path)
    assert ImplicitHuman.load(path,root).describe()==human.describe()
    data=json.loads(path.read_text());data['assets']['coverage']['sha256']='0'*64;path.write_text(json.dumps(data))
    try:ImplicitHuman.load(path,root)
    except ValueError:pass
    else:raise AssertionError('changed evidence accepted')
print('verified integrated materialization, real evidence availability, Gaussian conditioning, immutable input, source hash rejection')
first=human.materialize('skin-field');expected=first.evidence['fit']['parameters'][0]
first.evidence['fit']['parameters'][0]+=100
assert human.materialize('skin-field').evidence['fit']['parameters'][0]==expected
view=human.describe();view['assets']['coverage']['sha256']='tampered'
assert human.describe()['assets']['coverage']['sha256']!='tampered'
rank_one=PopulationBelief(['a','b'],['m','m'],[0,0],[[1,1],[1,1]],{})
tiny=rank_one.condition({'a':(1,1e-20,'m'),'b':(1,1e-20,'m')})
assert np.allclose(tiny.mean,[1,1],atol=1e-14)
assert np.allclose(tiny.cov,5e-21*np.ones((2,2)),rtol=1e-10,atol=1e-35)

assert "thermal" in human.describe()["materializations"]
thermal=human.materialize("thermal",profile="supine_blanket",seconds=2,dt=1).run()
assert thermal["configuration"]["posture"]=="lying"
assert len(thermal["node_temperature_C"])==3
assert thermal["audit"]["maximum_heat_balance_residual_W"]<1e-6

# Opening the evidence API works before canonical materialization has been built.
# A request for absent canonical data fails before importing optional runtime code.
with tempfile.TemporaryDirectory() as td:
    empty = ImplicitHuman.open(td)
    assert 'body' not in empty.describe()['materializations']
    try:
        empty.materialize('body')
    except ValueError as error:
        assert 'Required evidence unavailable: canonical_' in str(error), str(error)
    else:
        raise AssertionError('canonical body accepted an unbuilt workspace')

# Cached evidence must not hide a changed file from the root-loading body API.
with tempfile.TemporaryDirectory() as td:
    base = Path(td)
    canonical = base / 'data/derived/canonical'
    canonical.mkdir(parents=True)
    for name in ('anatomy', 'profile', 'brain', 'mechanics', 'body'):
        (canonical / (name + '.json')).write_text('{}')
    binding = ImplicitHuman.open(base)
    binding.save(base / 'human.json')
    cached = ImplicitHuman.load(base / 'human.json', base)
    (canonical / 'brain.json').write_text('{"changed":true}')
    try:
        cached.materialize('body')
    except ValueError as error:
        assert 'Evidence changed: canonical_brain' in str(error), str(error)
    else:
        raise AssertionError('changed cached canonical source accepted')
    try:
        cached.materialize('body', output_hz=10)
    except ValueError:
        pass
    else:
        raise AssertionError('implicit body simulation options accepted')

if 'body' in human.describe()['materializations']:
    canonical = human.materialize('body')
    assert canonical.__class__.__name__ == 'CanonicalBody'
    assert isinstance(canonical.describe(), dict)
    print('verified lazy canonical body materialization and current evidence identity')
