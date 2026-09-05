#!/usr/bin/env python3
"""Numerical constitutive and canonical assembly checks (no empirical claims)."""
import json
import sys
from pathlib import Path
import numpy as np
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
from ihm.assembly.mechanics import BodyMechanics, neo_hookean, tetra_force_energy


def main():
    # Objectivity, undeformed stress, and energy-gradient consistency.
    mu, lam = 1200., 3400.
    w, p = neo_hookean(np.eye(3), mu, lam)
    assert abs(w) < 1e-12 and np.linalg.norm(p) < 1e-12
    theta=.4
    r=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    f=np.array([[1.1,.1,0],[0,.95,.02],[0,0,1.03]])
    assert abs(neo_hookean(r@f,mu,lam)[0]-neo_hookean(f,mu,lam)[0]) < 1e-9
    x=np.array([[0.,0,0],[.1,0,0],[0,.2,0],[0,0,.3]])
    y=x@f.T
    forces,energy=tetra_force_energy(x,y,mu,lam)
    assert np.linalg.norm(forces.sum(axis=0)) < 1e-10
    assert np.linalg.norm(np.cross(y,forces).sum(axis=0)) < 1e-10
    eps=1e-7
    for i in range(4):
        for j in range(3):
            yp=y.copy();ym=y.copy();yp[i,j]+=eps;ym[i,j]-=eps
            numeric=-(tetra_force_energy(x,yp,mu,lam)[1]-tetra_force_energy(x,ym,mu,lam)[1])/(2*eps)
            assert abs(numeric-forces[i,j]) < 1e-5
    try: neo_hookean(np.diag([-1,1,1]),mu,lam)
    except ValueError: pass
    else: raise AssertionError('Inversion accepted')
    print('PASS neo-Hookean objectivity, tetra force/energy gradient and force/torque conservation')
    artifact=BASE/'data/derived/canonical/mechanics.json'
    if not artifact.exists():
        print('Canonical artifact unavailable; kernels only');return
    payload=json.loads(artifact.read_text());body=BodyMechanics.from_dict(payload)
    ids={e['id'] for e in payload['entities']}
    assert len(ids)==len(payload['entities'])
    assert len(payload['native_muscles'])==80
    source=json.loads((BASE/'data/derived/anatomy/opensim__Rajagopal__Rajagopal2016.json').read_text())
    source_by_name={m['name']:m for m in source['muscles']}
    for m in payload['muscles']:
        if 'source_name' in m:
            par=source_by_name[m['source_name']]['parameters']
            assert m['max_isometric_force_n']==float(par['max_isometric_force']['value'])
            assert m['pennation_angle_rad']==float(par['pennation_angle_at_optimal']['value'])
    graph={id:set() for id in ids}
    for link in payload['links']:
        graph[link['a']].add(link['b']);graph[link['b']].add(link['a'])
    visited=set();pending=[next(iter(ids))]
    while pending:
        id=pending.pop()
        if id in visited:continue
        visited.add(id);pending.extend(graph[id]-visited)
    assert visited==ids
    for e in payload['entities']:
        if e['role']=='fluid_cavity':assert e['material_volume_m3']==0 and e['constitutive']=='affine_boundary_carrier'
    assert all(n['entity_id'] in ids for m in payload['muscles'] for n in m['anchors'])
    assert all(m['canonical_entity_id'] in ids for m in payload['muscles'])
    baseline=body.step(.002,{})
    assert baseline['audit']['internal_force_residual_n'] < 1e-7
    assert baseline['audit']['internal_torque_residual_nm'] < 1e-7
    assert max(np.linalg.norm(e['translation_m']) for e in baseline['entities'].values()) < 1e-10
    target=next(m for m in payload['muscles'] if m.get('source_name')=='recfem_r')
    result=body.step(.01,{'activation':{target['id']:.2}})
    assert result['muscle_forces_n'][target['id']]>0
    assert max(np.linalg.norm(result['entities'][n['entity_id']]['translation_m']) for n in target['anchors'])>0
    assert np.linalg.norm(np.array(result['entities'][target['canonical_entity_id']]['deformation_gradient'])-np.eye(3))>0
    assert result['audit']['internal_force_residual_n'] < 1e-7
    assert result['audit']['internal_torque_residual_nm'] < 1e-6
    organ=next(e for e in payload['entities'] if 'liver' in e['name'].lower())
    fresh=BodyMechanics.from_dict(payload)
    result=fresh.step(.01,{'volume_ratios':{organ['id']:1.05}})
    assert abs(np.linalg.det(result['entities'][organ['id']]['deformation_gradient'])-1.05)<1e-10
    assert result['audit']['elastic_energy_j']>0
    pressure=BodyMechanics.from_dict(payload).step(.01,{'pressure_pa':{organ['id']:100.}})
    assert np.linalg.det(pressure['entities'][organ['id']]['deformation_gradient'])>1
    previous=(fresh.x.copy(),fresh.deformation.copy(),fresh.time,fresh.affine_work)
    for bad in [{'volume_ratios':{organ['id']:1.1},'external_forces_n':{organ['id']:[float('nan'),0,0]}},
                {'volume_ratios':{organ['id']:1.1},'pressure_pa':{organ['id']:1.}},
                {'volume_ratios':{organ['id']:True}}]:
        try:fresh.step(.02,bad)
        except ValueError:pass
        else:raise AssertionError('Invalid driver accepted')
        assert np.array_equal(fresh.x,previous[0]) and np.array_equal(fresh.deformation,previous[1])
        assert (fresh.time,fresh.affine_work)==previous[2:]
    try:fresh.step(True,{})
    except ValueError:pass
    else:raise AssertionError('Boolean dt accepted')
    # Time refinement tests dynamic response, independently of biological priors.
    end=[]
    for h in [.002,.001,.0005]:
        b=BodyMechanics.from_dict(payload)
        for _ in range(round(.01/h)):z=b.step(h,{'activation':{target['id']:.05}})
        end.append(np.array(z['entities'][target['anchors'][-1]['entity_id']]['translation_m']))
    assert np.linalg.norm(end[1]-end[2]) <= np.linalg.norm(end[0]-end[2])*1.2+1e-12
    # Sustained localized activation catches rotational/support divergence that
    # a single short perturbation cannot expose.
    stable=BodyMechanics.from_dict(payload)
    for _ in range(50):
        sustained=stable.step(.02,{'activation':{'recfem_r':.1}})
        assert max(np.linalg.norm(e['translation_m']) for e in sustained['entities'].values())<.5
        assert sustained['audit']['internal_force_residual_n']<1e-7
    assert np.linalg.norm((stable.mass[:,None]*stable.v).sum(axis=0))<1e-6
    assert abs(sum(e['mass_kg'] for e in payload['entities'])-77.1107029)<1e-8
    report={'status':'pass','entities':len(ids),'muscles':len(payload['muscles']),'native_muscles':80,'links':len(payload['links']),'time_refinement_displacement_difference_m':float(np.linalg.norm(end[1]-end[2])),'sustained_activation_duration_s':1.,'transactional_invalid_inputs':True,'constitutive_gradient_objectivity':True,'source_parameter_fidelity_80':True,'whole_graph_connected':True,'empirical_validation':False}
    (BASE/'data/derived/canonical/mechanics_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
