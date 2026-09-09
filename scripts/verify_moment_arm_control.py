#!/usr/bin/env python3
"""Numerical derivative and bounded torque allocation checks, not gait evidence."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.native.moment_arm_control import FittedMomentArms, polynomial_exponents, allocate_excitation
root=Path(__file__).resolve().parents[1]
assert polynomial_exponents(2,2)==[(0,0),(0,1),(0,2),(1,0),(1,1),(2,0)]
geometry=FittedMomentArms(root/'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled_FunctionBasedPathSet.xml')
count=0
for name,(coordinates,_,_) in geometry.paths.items():
    q={n:.07*(i+1) for i,n in enumerate(coordinates)}
    length,arms=geometry.length_and_moment_arms(name,q)
    for n,arm in arms.items():
        a=dict(q); b=dict(q);a[n]+=1e-6;b[n]-=1e-6
        derivative=(geometry.length_and_moment_arms(name,a)[0]-geometry.length_and_moment_arms(name,b)[0])/2e-6
        assert abs(derivative+arm)<1e-8,(name,n,derivative,arm)
        count+=1
u=allocate_excitation([[10,-10]],[5],baseline=0,regularization=1e-6)
assert np.all((u>=0)&(u<=1)) and abs(np.array([10,-10])@u-5)<1e-5
u=allocate_excitation([[10,-10]],[25],baseline=0,regularization=1e-6)
assert abs(u[0]-1)<1e-7 and u[1]<1e-7
print(f'PASS: {count} fitted analytic moment arms match finite differences; feasible and saturated allocation checked.')

# Inject bound roundoff and pass it through the real Python native-port validator.
from unittest.mock import patch
from types import SimpleNamespace
from threading import RLock
from ihm.native.mechanical_stream import NativeMechanicalStream
with patch('scipy.optimize.lsq_linear',return_value=SimpleNamespace(success=True,x=np.array([-1e-15,1+1e-15]))):
    values=allocate_excitation([[10,-10]],[5])
assert values.tolist()==[0.,1.]
stream=NativeMechanicalStream.__new__(NativeMechanicalStream)
stream.state={'muscles':{'a':{},'b':{}},'bodies':{}}
stream.lock=RLock()
stream._request=lambda line: stream.state
stream.advance(.01,actuation=dict(zip(('a','b'),values)))
for bad in (-1e-4,1.01,float('nan')):
    with patch('scipy.optimize.lsq_linear',return_value=SimpleNamespace(success=True,x=np.array([bad,.5]))):
        try: allocate_excitation([[10,-10]],[5])
        except RuntimeError: pass
        else: raise AssertionError('Infeasible optimizer output accepted')
print('PASS: tiny optimizer bound excursions clipped; native command validator accepts; material excursions rejected.')
# Native-only lumbar geometry adds a controlled joint without fitted polynomials.
from ihm.native.moment_arm_control import JointPosturalController
reference={'coordinates':{'pelvis_tilt':{'value':0.,'speed':0.},'pelvis_list':{'value':0.,'speed':0.},'lumbar_extension':{'value':0.,'speed':0.}},'muscles':{n:{'fiber_length_m':.12,'optimal_fiber_length_m':.12,'max_isometric_force_n':1000.} for n in ('native_extensor','native_flexor')}}
policy=JointPosturalController(reference,native_moment_arms={'native_extensor':{'lumbar_extension':.05},'native_flexor':{'lumbar_extension':-.05}},pelvis_gain=0,pelvis_damping=0,com_position_gain=0,com_velocity_gain=0)
import copy
perturbed=copy.deepcopy(reference);perturbed['coordinates']['lumbar_extension']['value']=-.2
command=policy.commands(perturbed)
assert command['native_extensor']>command['native_flexor']
assert policy.targets['lumbar_extension']==0 and policy.last_allocation['requested_torque_nm']['lumbar_extension']==20
policy.update_native_moment_arms({'native_extensor':{'lumbar_extension':.04},'native_flexor':{'lumbar_extension':-.04}})
assert policy.commands(perturbed)['native_extensor']>command['native_extensor']
for bad in ({'wrong_muscle':{'lumbar_extension':.1}},{'native_extensor':{'lumbar_extension':.1}},{'native_extensor':{'pelvis_ty':.1}}):
    try:policy.update_native_moment_arms(bad)
    except ValueError:pass
    else:raise AssertionError('Invalid native geometry accepted')
assert policy.identity()['unmodeled_muscles']==[]
print('PASS: native-only lumbar feedback allocation, refresh, and fail-fast catalog validation.')
from ihm.native.moment_arm_control import native_center_of_mass
body={'mass_kg':2.,'transform_ground':np.eye(4).tolist(),'mass_center_local_m':[0,1,0],'angular_velocity_rad_s':[0,0,2],'origin_velocity_m_s':[1,0,0]}
com,velocity=native_center_of_mass({'bodies':{'body':body}})
np.testing.assert_allclose(com,[0,1,0]);np.testing.assert_allclose(velocity,[-1,0,0])
sign_reference=copy.deepcopy(reference)
sign_reference['coordinates']['ankle_angle_r']={'value':0.,'speed':0.}
sign_reference['bodies']={'body':copy.deepcopy(body)}
sign_reference['bodies']['body']['angular_velocity_rad_s']=[0,0,0]
sign_reference['bodies']['body']['origin_velocity_m_s']=[0,0,0]
geometry={'native_extensor':{'ankle_angle_r':.05},'native_flexor':{'ankle_angle_r':-.05}}
policy=JointPosturalController(sign_reference,native_moment_arms=geometry,pelvis_gain=60,pelvis_damping=10,com_position_gain=200,com_velocity_gain=50)
state=copy.deepcopy(sign_reference);state['coordinates']['pelvis_tilt']['value']=.1
policy.commands(state)
assert abs(policy.last_allocation['requested_torque_nm']['ankle_angle_r']-3)<1e-12
state['coordinates']['pelvis_tilt']['value']=0
state['bodies']['body']['transform_ground'][0][3]=.1
state['bodies']['body']['origin_velocity_m_s'][0]=.2
policy.commands(state)
assert abs(policy.last_allocation['requested_torque_nm']['ankle_angle_r']+15)<1e-12
print('PASS: rigid body COM velocity includes angular offset; pelvis and COM ankle corrections have correct frame signs.')
policy=JointPosturalController(sign_reference,native_moment_arms=geometry,pelvis_gain=0,pelvis_damping=0,com_position_gain=200,com_velocity_gain=50,com_target_x_m=.1)
policy.commands(sign_reference)
assert abs(policy.last_allocation['requested_torque_nm']['ankle_angle_r']-10)<1e-12
assert policy.identity()['com_target_x_m']==.1
for bad in (float('nan'),float('inf'),True):
    try:JointPosturalController(sign_reference,native_moment_arms=geometry,com_target_x_m=bad)
    except ValueError:pass
    else:raise AssertionError('Invalid COM target accepted')
print('PASS: explicit COM target generates forward correction and rejects nonfinite/bool target.')
lateral_reference=copy.deepcopy(sign_reference)
for side in ('r','l'):lateral_reference['coordinates']['hip_adduction_'+side]={'value':0.,'speed':0.}
lateral_geometry={'native_extensor':{'hip_adduction_r':.05,'hip_adduction_l':-.05},'native_flexor':{'hip_adduction_r':-.05,'hip_adduction_l':.05}}
policy=JointPosturalController(lateral_reference,native_moment_arms=lateral_geometry,pelvis_gain=0,pelvis_damping=0,com_lateral_position_gain=200,com_lateral_velocity_gain=50)
state=copy.deepcopy(lateral_reference);state['bodies']['body']['transform_ground'][2][3]=.1;state['bodies']['body']['origin_velocity_m_s'][2]=.2
policy.commands(state)
assert abs(policy.last_allocation['requested_torque_nm']['hip_adduction_r']+15)<1e-12
assert abs(policy.last_allocation['requested_torque_nm']['hip_adduction_l']-15)<1e-12
print('PASS: positive lateral COM error gives empirically corrective right-negative/left-positive hip torque.')
nonzero_reference=copy.deepcopy(lateral_reference)
nonzero_reference['coordinates']['pelvis_tilt']['value']=-.138
nonzero_reference['coordinates']['pelvis_list']['value']=.02
policy=JointPosturalController(nonzero_reference,native_moment_arms=lateral_geometry)
policy.commands(nonzero_reference)
assert all(abs(v)<1e-12 for v in policy.last_allocation['requested_torque_nm'].values())
print('PASS: vestibular feedback references initial pose, including nonzero tilt/list.')
equilibrium={'native_extensor':.27,'native_flexor':.41}
policy=JointPosturalController(reference,native_moment_arms={'native_extensor':{'lumbar_extension':.05},'native_flexor':{'lumbar_extension':-.05}},equilibrium_excitations=equilibrium)
command=policy.commands(reference)
np.testing.assert_allclose([command[n] for n in equilibrium],list(equilibrium.values()),atol=1e-10)
for bad in ({'native_extensor':.2},{'native_extensor':float('nan'),'native_flexor':.3}):
    try:JointPosturalController(reference,equilibrium_excitations=bad)
    except ValueError:pass
    else:raise AssertionError('Invalid equilibrium catalog accepted')
print('PASS: zero-error equilibrium excitation vector preserved; malformed catalogs rejected.')
force_reference=copy.deepcopy(reference);force_reference['time_s']=0.
force_reference['muscles']['native_extensor'].update(tendon_force_n=100.,activation=.2)
force_reference['muscles']['native_flexor'].update(tendon_force_n=40.,activation=.1)
geometry={'native_extensor':{'lumbar_extension':.05},'native_flexor':{'lumbar_extension':-.05}}
policy=JointPosturalController(force_reference,native_moment_arms=geometry)
policy.commands(force_reference)
assert abs(policy.last_allocation['actual_muscle_torque_nm']['lumbar_extension']-3)<1e-12
observed=copy.deepcopy(force_reference);observed['time_s']=.03;observed['muscles']['native_extensor']['tendon_force_n']=160.
first_commands=policy.commands(observed)
assert abs(policy.last_allocation['actual_delta_torque_nm']['lumbar_extension']-3)<1e-12
assert policy.last_allocation['initial_reference_muscle_torque_nm']['lumbar_extension']==3.
assert policy.last_allocation['actual_activation']['native_extensor']==.2
policy.feedforward['lumbar_extension']=200.
second_commands=policy.commands(observed)
assert first_commands!=second_commands
assert abs(policy.last_allocation['actual_muscle_torque_nm']['lumbar_extension']-6)<1e-12
assert policy.last_allocation['excitations_at_upper_bound']==['native_extensor']
policy.update_native_moment_arms(geometry,time_s=.02);policy.commands(observed)
assert abs(policy.last_allocation['torque_observation_clock']['native_moment_arm_age_s']-.01)<1e-12
print('PASS: measured tendon torques stay independent of issued command; reference deltas, bounds, activation, and native geometry timestamps recorded.')
