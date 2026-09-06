"""Bounded hair/native port fixture; no native engine or anatomical assets."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.hair_dynamics import ElasticHairState
from ihm.assembly.hair_feedback import HairFeedback


class Registration:
    bodies=['fixture-body'];basis=np.eye(3);global_map=np.eye(4)
    def transforms(self,state):return {'fixture-body':np.asarray(state['bodies']['fixture-body']['transform_ground'])}


class Native:
    """Finite rigid fixture with explicit mass and isotropic rotational inertia."""
    def __init__(self):self.x=np.zeros(3);self.v=np.array([0.,0.,.001]);self.omega=np.zeros(3);self.rotation=np.eye(3);self.time=0.;self.mass=1.;self.inertia=.1;self.tokens=0
    def snapshot(self):
        t=np.eye(4);t[:3,3]=self.x;t[:3,:3]=self.rotation
        return {'time_s':self.time,'gravity_m_s2':[0,0,0],'bodies':{'fixture-body':{'transform_ground':t.tolist(),'origin_velocity_m_s':self.v.tolist(),'angular_velocity_rad_s':self.omega.tolist()}}}
    def checkpoint(self):self.tokens+=1;return self.snapshot()
    def restore(self,s):
        body=s['bodies']['fixture-body'];transform=np.asarray(body['transform_ground'])
        self.x=transform[:3,3].copy();self.rotation=transform[:3,:3].copy();self.v=np.asarray(body['origin_velocity_m_s']).copy();self.omega=np.asarray(body['angular_velocity_rad_s']).copy();self.time=s['time_s']
    def release(self,s):self.tokens-=1
    def advance(self,dt,loads,actuation):
        force=np.sum([p['force_n'] for p in loads],axis=0) if loads else np.zeros(3)
        torque=np.sum([np.cross(np.asarray(p['point_m'])-self.x,p['force_n']) for p in loads],axis=0) if loads else np.zeros(3)
        alpha=torque/self.inertia;self.rotation=Rotation.from_rotvec(dt*self.omega+.5*dt*dt*alpha).as_matrix()@self.rotation;self.omega+=dt*alpha
        a=force/self.mass;self.x+=dt*self.v+.5*dt*dt*a;self.v+=dt*a;self.time+=dt
        return self.snapshot()


def fixture():
    x=np.array([[-1.,-1.,0.],[1.,-1.,0.],[0.,1.,0.]])
    data={'centerlines_m':[[0,0,0],[.005,0,0],[.01,0,0],[.015,0,0]],'strand_offsets':[0,4],'radius_m':[40e-6],
          'tensile_modulus_pa':7.11e9,'bending_modulus_pa':5.7e9,'density_kg_m3':1312}
    hair=ElasticHairState(data,max_substep_s=.001)
    owner=HairFeedback(hair,x,np.array([[0,1,2]]),np.zeros(3,int),Registration(),[0],[[.25,.25,.5]],
                       identity={'source_id':'fixture-surface','status':'synthetic port fixture'},friction_static=0.,friction_kinetic=0.)
    return Native(),owner


def main():
    _,uncoupled=fixture();beam=uncoupled.hair;reference=beam.position_m.copy()
    beam.step(.001,roots_m=beam.root_m,root_tangents=beam.tangent,gravity_m_s2=[0,0,0])
    assert np.array_equal(beam.position_m,reference)
    beam_receipt=beam.step(.001,roots_m=beam.root_m,root_tangents=beam.tangent,gravity_m_s2=[0,-9.81,0])
    assert beam.position_m[-1,1]<0 and np.array_equal(beam.position_m[:2],reference[:2])
    assert np.linalg.norm(beam_receipt['hair_body_impulse_residual_ns'])<1e-18
    assert beam_receipt['hair_numerical_energy_defect_j']<=0
    before=beam.checkpoint();beam.max_substep_s=1e-9
    try:beam.step(.02,roots_m=beam.root_m,root_tangents=beam.tangent,gravity_m_s2=[0,0,0])
    except ValueError:pass
    else:raise AssertionError('Unbounded substep request accepted')
    assert np.array_equal(beam.position_m,before['position_m']) and beam.time_s==before['time_s']
    native,owner=fixture();hair=owner.hair
    initial_momentum=native.mass*native.v+(hair.mass_kg[:,None]*hair.velocity_m_s).sum(0)
    initial_energy=.5*native.mass*float(native.v@native.v)+.5*native.inertia*float(native.omega@native.omega)+hair.energy_j()
    initial_angular=np.cross(native.x,native.mass*native.v)+native.inertia*native.omega+np.cross(hair.position_m,hair.mass_kg[:,None]*hair.velocity_m_s).sum(0)
    result,receipt=owner.advance(native,.001,[],{},tolerance_m=1e-12)
    final_momentum=native.mass*native.v+(hair.mass_kg[:,None]*hair.velocity_m_s).sum(0)
    assert receipt['contact_count']>0 and np.linalg.norm(receipt['surface_contact_impulse_ns'])>0
    assert np.linalg.norm(final_momentum-initial_momentum)<1e-18
    assert np.linalg.norm(receipt['hair_body_impulse_residual_ns'])<1e-18
    final_angular=np.cross(native.x,native.mass*native.v)+native.inertia*native.omega+np.cross(hair.position_m,hair.mass_kg[:,None]*hair.velocity_m_s).sum(0)
    assert np.linalg.norm(final_angular-initial_angular)<1e-18
    assert np.linalg.norm(receipt['hair_body_angular_impulse_residual_nms'])<1e-18
    assert np.linalg.norm(native.omega)>0 and receipt['within_small_deflection']
    assert np.isfinite(receipt['interface_work_quadrature_defect_j'])
    final_energy=.5*native.mass*float(native.v@native.v)+.5*native.inertia*float(native.omega@native.omega)+hair.energy_j()
    assert abs((final_energy-initial_energy)-receipt['hair_plus_native_interface_work_j'])<1e-18
    assert hair.time_s==native.time and native.tokens==0
    assert json.loads(json.dumps(owner.frame()))['viewer_synchronized'] is False
    saved=owner.checkpoint();native_saved=native.snapshot()
    owner.advance(native,.001,[],{},tolerance_m=1e-12);expected=hair.position_m.copy();expected_velocity=hair.velocity_m_s.copy()
    owner.restore(saved);native.restore(native_saved);owner.advance(native,.001,[],{},tolerance_m=1e-12)
    assert np.array_equal(hair.position_m,expected) and np.array_equal(hair.velocity_m_s,expected_velocity)
    native,owner=fixture();before=owner.checkpoint();native_before=native.snapshot()
    try:owner.advance(native,.001,[],{},tolerance_m=1e-30,max_iterations=1)
    except RuntimeError:pass
    else:raise AssertionError('Nonconverged hair/body interval accepted')
    assert native.snapshot()==native_before and native.tokens==0
    assert np.array_equal(owner.hair.position_m,before['hair']['position_m']) and owner.hair.time_s==0
    # A loose surface tolerance must not commit roots too far from the native
    # endpoint for the next coupled interval to accept the same state.
    native,owner=fixture();native.mass=.001;native.v=np.array([0.,0.,.01])
    owner.advance(native,.01,[],{})
    owner.advance(native,.01,[],{})
    assert owner.last['root_position_residual_m']<=owner.last['root_tolerance_m']
    assert owner.last['clamp_position_residual_m']<=owner.last['root_tolerance_m']
    native,owner=fixture();native.time=.1
    try:owner.advance(native,.001,[],{})
    except ValueError:pass
    else:raise AssertionError('Independent hair/native clocks accepted')
    assert native.tokens==0 and owner.hair.time_s==0
    print('PASS moving source face contact, finite-body force/torque reaction, total linear/angular momentum, explicit energy receipt, exact checkpoint replay and atomic nonconvergence rollback')
    print({k:receipt[k] for k in ('contact_count','normal_impact_dissipation_j','hair_numerical_energy_defect_j','interface_work_quadrature_defect_j','boundary_position_residual_m')})


if __name__=='__main__':main()
