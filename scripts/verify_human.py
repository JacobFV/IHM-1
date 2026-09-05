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
