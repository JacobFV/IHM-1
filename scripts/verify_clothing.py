#!/usr/bin/env python3
"""Independent force, energy and Coulomb impulse checks for cloth mechanics."""
from pathlib import Path
import sys
import hashlib
import json
import subprocess
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.clothing import Cloth, PlaneContact


def main():
    audit={'schema_version':1,'parameter_status':'engineering priors; no material calibration'}
    x=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    cloth=Cloth(x,[[0,1,2]],areal_density_kg_m2=.2,edge_stiffness_n_m=12.)
    moved=x.copy();moved[1,0]=1.1
    f,energy=cloth.elastic_forces(moved)
    assert np.linalg.norm(f.sum(axis=0))<1e-12
    assert np.linalg.norm(np.cross(moved,f).sum(axis=0))<1e-12
    eps=1e-6; plus=moved.copy();minus=moved.copy();plus[1,0]+=eps;minus[1,0]-=eps
    derivative=(cloth.elastic_forces(plus)[1]-cloth.elastic_forces(minus)[1])/(2*eps)
    assert abs(f[1,0]+derivative)<1e-7
    assert energy>0 and np.isclose(cloth.mass_kg.sum(),.1)
    audit['elastic']={'force_gradient_error_n':float(abs(f[1,0]+derivative)),
                      'net_force_n':f.sum(axis=0).tolist(),
                      'net_torque_nm':np.cross(moved,f).sum(axis=0).tolist(),
                      'mass_kg':float(cloth.mass_kg.sum())}
    plane=PlaneContact(normal=[0,0,1],offset_m=0,friction_static=.6,friction_kinetic=.4)
    # Impulses are physically measured in N s. Sliding decelerates and obeys the cone.
    contact=plane.resolve(np.array([[0.,0.,-.001]]),np.array([[2.,0.,-1.]]),np.array([.2]),.01)
    assert contact['positions_m'][0,2]>=0
    assert np.allclose(contact['impulses_ns'][0],[-.08,0,.2])
    assert contact['friction_dissipation_j']>0
    assert np.allclose(contact['body_reaction_impulses_ns'],-contact['impulses_ns'])
    stick=plane.resolve(np.array([[0.,0.,0.]]),np.array([[.2,0.,-1.]]),np.array([.2]),.01)
    assert np.allclose(stick['velocities_m_s'],0)
    free=plane.resolve(np.array([[0.,0.,1.]]),np.array([[2.,0.,0.]]),np.array([.2]),.01)
    assert np.allclose(free['impulses_ns'],0)
    moving=PlaneContact(normal=[0,0,1],offset_m=0,friction_static=.6,friction_kinetic=.4,surface_velocity_m_s=[1,0,0])
    dragged=moving.resolve(np.array([[0.,0.,0.]]),np.array([[0.,0.,-1.]]),np.array([.2]),.01)
    kinetic_change=.1*np.sum(dragged['velocities_m_s']**2)-.1
    work_residual=kinetic_change-dragged['prescribed_surface_work_j']+dragged['friction_dissipation_j']+dragged['normal_impact_dissipation_j']
    assert abs(work_residual)<1e-14
    audit['contact']={'sliding_impulse_ns':contact['impulses_ns'].tolist(),
                      'static_final_velocity_m_s':stick['velocities_m_s'].tolist(),
                      'moving_surface_energy_residual_j':float(work_residual),
                      'friction_static':.6,'friction_kinetic':.4}
    # A loaded strip rests on its body support; tangential friction arrests it.
    x=x*.1;x[:,2]=.0001
    cloth=Cloth(x,[[0,1,2]],areal_density_kg_m2=.2,edge_stiffness_n_m=12.)
    cloth.velocity_m_s[:]=[.1,0,0]
    reaction=np.zeros(3);dissipation=0.;momentum_error=0.;energy_defect=0.
    for _ in range(400):
        result=cloth.step(.0005,gravity_m_s2=[0,0,-9.81],contact=plane)
        reaction+=result['body_reaction_impulse_ns'];dissipation+=result['friction_dissipation_j']
        momentum_error=max(momentum_error,float(np.linalg.norm(result['momentum_residual_ns'])))
        energy_defect+=result['numerical_energy_defect_j']
    assert np.max(abs(cloth.velocity_m_s))<1e-10
    assert np.min(cloth.position_m[:,2])>=-1e-12
    assert reaction[2]<0 and dissipation>0
    assert momentum_error<1e-14
    audit['supported_strip']={'duration_s':cloth.time_s,'dt_s':.0005,
                             'body_reaction_impulse_ns':reaction.tolist(),
                             'friction_dissipation_j':dissipation,
                             'max_momentum_residual_ns':momentum_error,
                             'numerical_energy_defect_j':energy_defect}
    before=(cloth.position_m.copy(),cloth.velocity_m_s.copy(),cloth.time_s)
    for dt in [True,float('nan'),0.,1.]:
        try:cloth.step(dt)
        except ValueError:pass
        else:raise AssertionError('invalid step accepted')
        assert np.array_equal(cloth.position_m,before[0]) and np.array_equal(cloth.velocity_m_s,before[1]) and cloth.time_s==before[2]
    convergence=[]
    for dt in [.0001,.00005,.000025]:
        c=Cloth(x,[[0,1,2]],areal_density_kg_m2=.2,edge_stiffness_n_m=12.)
        c.position_m[1,0]+=.001;initial=c.elastic_forces(c.position_m)[1];error=0.
        for _ in range(round(.1/dt)):
            r=c.step(dt,gravity_m_s2=[0,0,0]);error=max(error,abs(r['kinetic_energy_j']+r['elastic_energy_j']-initial)/initial)
        convergence.append({'dt_s':dt,'max_relative_total_energy_error':error})
    assert convergence[-1]['max_relative_total_energy_error']<.005
    assert all(b['max_relative_total_energy_error']<.55*a['max_relative_total_energy_error'] for a,b in zip(convergence,convergence[1:]))
    audit['elastic_time_refinement']=convergence
    root=Path(__file__).resolve().parents[1]
    out=root/'data/derived/clothing';out.mkdir(parents=True,exist_ok=True)
    # Retain the exact same geometric materialization that the browser computes,
    # including source and constructor identities. No server/API mutation needed.
    materialize=r'''
import fs from 'node:fs';import zlib from 'node:zlib';import crypto from 'node:crypto';
import {buildGarments} from './app/src/clothing.js';
const path='data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz';
const bytes=fs.readFileSync(path),sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const source={id:'body-bp3d-FJ2810',path,geometry_sha256:sha(bytes),license:'CC BY 4.0',attribution:'BodyParts3D, The Database Center for Life Science'};
const garments=buildGarments(JSON.parse(zlib.gunzipSync(bytes)),source);
fs.writeFileSync('data/derived/clothing/garments.json',JSON.stringify({schema_version:1,source,constructor_sha256:sha(fs.readFileSync('app/src/clothing.js')),garments}));
'''
    subprocess.run(['node','--input-type=module','-e',materialize],cwd=root,check=True,capture_output=True,text=True)
    generated=json.loads((out/'garments.json').read_text())
    garment_stats=[]
    for garment in generated['garments']:
        v=np.asarray(garment['positions']).reshape(-1,3);f=np.asarray(garment['indices']).reshape(-1,3)
        c=Cloth(v,f,areal_density_kg_m2=.18,edge_stiffness_n_m=12.)
        garment_stats.append({'id':garment['id'],'vertices':len(v),'triangles':len(f),
                              'assumed_areal_density_kg_m2':.18,'mass_at_assumed_density_kg':float(c.mass_kg.sum()),
                              'physical_contact_solved':False})
    audit['garments']=garment_stats
    audit['sources']={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in [
        'ihm/assembly/clothing.py','scripts/verify_clothing.py','app/src/clothing.js',
        'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz']}
    audit['limitations']=['No bending constitutive law or fabric calibration','Node-plane contact only; no triangle CCD or self-contact',
                          'Prescribed body reactions returned, not applied to whole-body owner','No deformable genital tucking',
                          'Garment display is source-conditioned geometry with kinematic respiratory motion, not the independent contact solve',
                          'Numerical energy defect includes finite-step and position-projection error; it is reported, not balanced away']
    (out/'verification.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('PASS elastic force/energy/torque, contact impulses, static/kinetic friction and supported dynamics')


if __name__=='__main__':main()
