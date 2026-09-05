"""Numerical and source-boundary checks, independent of clinical validation."""
from pathlib import Path
import sys,json,hashlib,tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.thermal import total_resistance, heat_step_audit, run_thermal, load_source, bedding_boundary

rt,ret=total_resistance(2.78)
assert abs(rt-.4309)<1e-14 and abs(ret-rt/(.38*16.5))<1e-14
for bad in [0,-1,float('nan'),float('inf')]:
    try:total_resistance(bad)
    except ValueError:pass
    else:raise AssertionError('invalid resistance accepted')
# Two equal-capacity nodes, one conservative conductance, backward Euler.
c=np.array([2.,2.]);w=np.array([[0.,1.],[1.,0.]])
k=np.diag(w.sum(axis=1))-w; old=np.array([10.,0.]);new=np.linalg.solve(np.diag(c)+k,c*old)
a=heat_step_audit(c,w,np.zeros(2),np.zeros(2),np.zeros(2),old,new,1.)
assert a['heat_balance_residual_W']<1e-12 and a['internal_column_sum_max_W_K']<1e-12
# One-node finite-capacity cooling checks boundary sign and endpoint evaluation.
cooled=heat_step_audit([2],[[0]],[1],[0],[0],[10],[20/3],1)
assert cooled['heat_balance_residual_W']<1e-12 and cooled['boundary_heat_into_body_W']<0
# Internal directed circulatory cycle is conservative without symmetry.
cycle=np.array([[0,0,2],[2,0,0],[0,2,0]],float)
kcycle=np.diag(cycle.sum(axis=1))-cycle
old3=np.array([10.,20.,30.]);new3=np.linalg.solve(np.eye(3)+kcycle,old3)
assert heat_step_audit(np.ones(3),cycle,np.zeros(3),np.zeros(3),np.zeros(3),old3,new3,1)['heat_balance_residual_W']<1e-12
# Nonconservative flow matrix must be diagnosed, never silently balanced.
bad=w.copy();bad[0,1]=2
assert heat_step_audit(c,bad,np.zeros(2),np.zeros(2),np.zeros(2),old,new,1.)['internal_column_sum_max_W_K']>0
# Invalid clocks must fail before source import or model initialization.
for seconds,dt in [(True,1),(1,True),(1e-300,30),(3600,1e-10),(604801,1),(3600,.01),(1,1.01),(1,float('inf'))]:
    try:run_thermal(Path('/not-a-source-checkout'),seconds=seconds,dt=dt)
    except ValueError:pass
    else:raise AssertionError('invalid/unbounded thermal clock accepted')
try:heat_step_audit([1],[[0]],[1],[0],[1e308],[0],[1e308],1e-300)
except FloatingPointError:pass
else:raise AssertionError('overflowed heat ledger accepted')
for capacity,transfer,boundary in [([],[],[]),([1],[[-1]],[0]),([1],[[0]],[-1])]:
    try:heat_step_audit(capacity,transfer,boundary,[0],[0],[0],[0],1)
    except ValueError:pass
    else:raise AssertionError('invalid thermal coefficients accepted')
if '--artifacts' in sys.argv:
    root=Path(__file__).resolve().parents[1];index=json.loads((root/'data/derived/thermal/index.json').read_text())
    assert index['models']==index['runs']  # API descriptor compatibility
    assert index['node_count']==85 and len(index['body_regions'])==17
    package,_=load_source(root)
    import importlib
    th=importlib.import_module(package.__name__+'.thermoregulation'); original=th.dry_r
    original_wet=th.wet_r
    try:
        with bedding_boundary(package,2.78):
            raise RuntimeError('intentional boundary failure')
    except RuntimeError:pass
    assert th.dry_r is original and th.wet_r is original_wet
    with tempfile.TemporaryDirectory() as directory:
        temporary=Path(directory);raw=temporary/'data/raw/thermal';raw.mkdir(parents=True)
        (raw/'JOS-3').symlink_to(root/'data/raw/thermal/JOS-3',target_is_directory=True)
        metadata=json.loads((root/'data/raw/thermal/source.json').read_text());metadata['source_sha256']={}
        (raw/'source.json').write_text(json.dumps(metadata))
        try:load_source(temporary)
        except ValueError as error:assert 'inventory' in str(error)
        else:raise AssertionError('empty source hash inventory accepted')
    short=run_thermal(root,'supine_blanket',seconds=15,dt=15)
    assert short['source']['runtime_adapter_sha256']==hashlib.sha256((root/'ihm/native/thermal.py').read_bytes()).hexdigest()
    assert short['source']['adapter_sha256']==json.loads((root/'data/raw/thermal/source.json').read_text())['adapter_sha256']
    assert th.dry_r is original
    again=run_thermal(root,'supine_blanket',seconds=15,dt=15)
    assert short['node_temperature_C']==again['node_temperature_C']
    for case in index['runs']:
        a=json.loads((root/case['trajectory_path']).read_text())
        coefficients=json.loads((root/case['coefficient_path']).read_text())
        capacities=np.asarray(coefficients['capacity_J_K']);conductance=np.asarray(coefficients['conductance_W_K'])
        assert capacities.shape==(85,) and np.isfinite(capacities).all() and (capacities>0).all()
        assert conductance.shape==(85,85) and np.isfinite(conductance).all() and (conductance>=0).all()
        assert np.allclose(conductance,conductance.T,rtol=0,atol=1e-12)
        assert a['time_s'][-1]==3600 and all(np.isfinite([v for v in ch['values'] if v is not None]).all() for ch in a['channels'])
        assert a['audit']['maximum_heat_balance_residual_W']<1e-7
        assert a['audit']['maximum_internal_column_sum_W_K']<1e-10
        assert case['refinement']['dt1.875_vs_dt0.9375_final_max_node_temperature_C'] < .001
        assert case['refinement']['dt1.875_vs_dt0.9375_final_max_node_temperature_C'] < case['refinement']['dt3.75_vs_dt1.875_final_max_node_temperature_C']
print('Thermal source boundary, analytic heat closure and artifact checks PASS')
