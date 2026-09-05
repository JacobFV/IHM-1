"""Analytic RC descriptor, native hydraulic conservation and transfer checks."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.coupling.circuit import FluidCircuit, SKIN_NODES, SKIN_BOUNDARIES, build_skin_circuit


def fixture():
    def prop(v,u):return dict(si_value=v,si_unit=u,value=v,unit=u,status='finite')
    def node(name,p):return dict(id='fluid:'+name,name=name,family='fluid',properties={'Pressure':prop(p,'Pa'),'Volume':prop(1,'m^3')})
    def path(name,a,b,properties):return dict(id='fluid:'+name,name=name,family='fluid',source='fluid:'+a,target='fluid:'+b,properties=properties,gate_states={},circuits=['fluid:RC'])
    return dict(source={'sha256':'analytic_fixture'},nodes=[node('in',0),node('x',0),node('ground',0)],paths=[
        path('R','in','x',{'Resistance':prop(2e9,'Pa s/m^3'),'Flow':prop(0,'m^3/s')}),
        path('C','x','ground',{'Compliance':prop(3e-9,'m^3/Pa'),'Flow':prop(0,'m^3/s')})])


def main():
    m=FluidCircuit.from_native(fixture(),['in','x','ground'],['in','ground'],circuit_name='RC')
    pole=m.finite_poles_per_s();np.testing.assert_allclose(pole,[-1/6],rtol=1e-9)
    response=m.response(np.array([0,1j,2j]),'in')
    np.testing.assert_allclose(response['pressures_pa_per_pa']['x'],1/(1+6*np.array([0,1j,2j])),rtol=1e-10)
    dt=.01
    for _ in range(600): result=m.step(dt,{'in':100})
    np.testing.assert_allclose(result['pressures_pa']['x'],100*(1-(1+dt/6)**-600),rtol=1e-8)
    assert result['balance']['max_abs_free_node_residual_m3_s']<1e-18
    np.testing.assert_allclose(result['volume_deltas_m3']['x'],result['flows_m3_s']['C']*dt)
    restored=FluidCircuit.from_dict(json.loads(json.dumps(m.to_dict())))
    np.testing.assert_allclose(m.step(.1)['pressures_pa']['x'],restored.step(.1)['pressures_pa']['x'])
    gated=fixture();gated['paths'][0]['gate_states']={'Valve':'Closed'}
    closed=FluidCircuit.from_native(gated,['in','x','ground'],['in','ground'],circuit_name='RC')
    assert closed.step(.1,{'in':-100})['gate_violations']
    gated['paths'][0]['gate_states']={'Valve':'Open'}
    opened=FluidCircuit.from_native(gated,['in','x','ground'],['in','ground'],circuit_name='RC')
    assert opened.step(.1,{'in':100})['flows_m3_s']['R']==0
    sourced=fixture();sourced['paths'][0]['properties']['PressureSource']={'si_value':50.,'si_unit':'Pa'}
    driven=FluidCircuit.from_native(sourced,['in','x','ground'],['in','ground'],circuit_name='RC')
    np.testing.assert_allclose(driven.step(.1)['pressures_pa']['x'],50*.1/(6+.1))
    try:m.step(0)
    except ValueError:pass
    else:raise AssertionError('zero dt accepted')
    g=fixture();g['paths'][0]['properties']['Resistance']['si_unit']='K/W'
    try:FluidCircuit.from_native(g,['in','x','ground'],['in','ground'],circuit_name='RC')
    except ValueError:pass
    else:raise AssertionError('wrong physical dimension accepted')
    root=Path(__file__).resolve().parents[1];native=root/'data/derived/native-circuits/graph.json'
    if native.exists():
        graph=json.loads(native.read_text());m=FluidCircuit.from_native(graph,SKIN_NODES,SKIN_BOUNDARIES)
        r=m.step(.02)
        assert r['balance']['max_abs_free_node_residual_m3_s']<1e-13
        assert not r['gate_violations']
        assert all(np.isfinite(v) for v in r['pressures_pa'].values())
        assert max(m.finite_poles_per_s().real)<1e-10
        serialized=json.loads(json.dumps(m.to_dict(),allow_nan=False));restored=FluidCircuit.from_dict(serialized)
        a=m.step(.02);b=restored.step(.02)
        np.testing.assert_allclose(list(a['pressures_pa'].values()),list(b['pressures_pa'].values()))
        native_q={p['name']:p['properties'].get('Flow',{}).get('si_value') for p in m.native_paths}
        for name in ('SkinVToSkinE1','SkinToLymphValve','LymphToVenaCava'):
            np.testing.assert_allclose(r['flows_m3_s'][name],native_q[name],rtol=.005,atol=1e-14)
    print('circuit predictor: analytic RC pole/step/Laplace, KCL/storage, SI validation and native frozen one-step PASS')

if __name__=='__main__':main()
