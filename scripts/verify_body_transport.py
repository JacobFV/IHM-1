"""Executable conservation, intervention, refinement and restart regression checks."""
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.assembly.skin_transport import RegionalSkinTransport, albumin_flux_kg_s


def verify(root):
    initial = RegionalSkinTransport.from_sources(root)
    assert initial.to_dict()['native_blood_storage_connected'] is False
    assert abs(albumin_flux_kg_s(0., 40., 20., 1e-12) - 2e-11) < 1e-24
    assert albumin_flux_kg_s(-1., 40., 20., 1e-12) == 0.
    # SI dimensional checks: Pa/(Pa s/m3)=m3/s, ug/mL=.001 kg/m3.
    p=initial.parameters
    rates,protein,_,pressure=initial.rates()
    assert abs(rates['arterial_supply']*p['resistance_pa_s_m3']['arterial_supply']-(p['arterial_pa']-pressure['blood'])) < 1e-10
    assert 40 < initial.mass['blood']/initial.volumes['blood'] < 50
    for q in (-1e-15,0.,1e-15):
        assert abs(albumin_flux_kg_s(q,40.,20.,1e-12)-2e-11) < 2e-15
    assert initial.rates(venous_delta_pa=1e5)[0]['lymph_return']==0.
    results = {}
    for name, kw in [('baseline', {}), ('venous_pressure', {'venous_delta_pa': 600.}),
                     ('inflow_reduction', {'inflow_fraction': .25}),
                     ('lymph_obstruction', {'lymph_obstruction': 1.})]:
        model = RegionalSkinTransport.from_dict(initial.to_dict())
        for _ in range(120): model.step(5., **kw)
        state = model.to_dict()
        assert min(state['volumes_m3'].values()) > 0
        assert min(state['albumin_kg'].values()) >= 0
        assert abs(state['balance']['volume_residual_m3']) < 1e-20
        assert abs(state['balance']['albumin_residual_kg']) < 1e-18
        assert max(map(abs,state['balance']['node_volume_residual_m3'].values())) < 1e-18
        assert max(map(abs,state['balance']['node_albumin_residual_kg'].values())) < 1e-16
        assert state['integrated_flows_m3']['lymph_return'] >= 0
        results[name] = state
    assert results['lymph_obstruction']['volumes_m3']['interstitium'] > results['baseline']['volumes_m3']['interstitium']
    assert results['venous_pressure']['volumes_m3']['interstitium'] > results['baseline']['volumes_m3']['interstitium']
    assert results['inflow_reduction']['integrated_flows_m3']['arterial_supply'] < results['baseline']['integrated_flows_m3']['arterial_supply']
    a = RegionalSkinTransport.from_sources(root)
    a.step(30.)
    b = RegionalSkinTransport.from_dict(json.loads(json.dumps(a.to_dict())))
    assert a.step(40.) == b.step(40.)
    for field,value in [('time_s',-1.),('volumes_m3',{'blood':0.})]:
        bad=initial.to_dict();bad[field]=value
        try: RegionalSkinTransport.from_dict(bad)
        except ValueError: pass
        else: raise AssertionError('invalid restart accepted')
    refinements = []
    for dt in (1., .5, .25):
        m = RegionalSkinTransport.from_sources(root)
        for _ in range(round(60 / dt)): m.step(dt, venous_delta_pa=600.)
        refinements.append(m.to_dict()['volumes_m3']['interstitium'])
    assert abs(refinements[2]-refinements[1]) < abs(refinements[1]-refinements[0])
    stress = RegionalSkinTransport.from_sources(root)
    stress.step(3600., venous_delta_pa=5000., inflow_fraction=0., lymph_obstruction=1.)
    assert min(stress.volumes.values()) > 0 and min(stress.mass.values()) >= 0
    return {'scenarios': results, 'refinement_interstitial_volumes_m3': refinements,
            'refinement_difference_ratio':abs(refinements[1]-refinements[0])/abs(refinements[2]-refinements[1]),
            'dimensional_checks_passed': True,'restart_exact': True, 'stress_balance': stress.to_dict()['balance']}


if __name__ == '__main__':
    result = verify(Path(__file__).resolve().parents[1])
    print(json.dumps({k:v for k,v in result.items() if k != 'scenarios'}, indent=2))
