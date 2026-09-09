#!/usr/bin/env python3
"""Empirical native matched-checkpoint ankle torque-request sign experiment."""
from pathlib import Path
import json,sys,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import JointPosturalController,native_center_of_mass
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(tempfile.mkdtemp(prefix='native_ankle_pulse_',dir=ROOT/'data/research/locomotion_control'))
native=NativeMechanicalStream(ROOT,OUT/'native',environment='upright',target_mass_kg=70,augmented_registration='data/derived/mechanics/initial_stance86/registration.json')
report={'schema':'ihm.native-ankle-pulse.v1','seconds':.05,'step_s':.01,'kp':100,'kd':15,'baseline':.02,'arms':{},'scope':'Native matched-checkpoint local ankle excitation response; not a stance stability test.'}
try:
    initial=native.snapshot();checkpoint=native.checkpoint()
    for label,torque in [('zero',0.),('positive',20.),('negative',-20.)]:
        native.restore(checkpoint)
        policy=JointPosturalController(initial,kp=100,kd=15,baseline=.02,pelvis_gain=0,pelvis_damping=0,com_position_gain=0,com_velocity_gain=0,feedforward_torques={'ankle_angle_r':torque,'ankle_angle_l':torque})
        frames=[]
        for i in range(5):
            command=policy.commands(native.state)
            state=native.advance(.01,actuation=command)
            com,velocity=native_center_of_mass(state)
            frames.append({'time_s':state['time_s'],'com_ground_m':com.tolist(),'com_velocity_ground_m_s':velocity.tolist(),'ankle_r':state['coordinates']['ankle_angle_r'],'ankle_l':state['coordinates']['ankle_angle_l'],'allocation':policy.last_allocation,'commands':command})
        report['arms'][label]={'requested_feedforward_per_ankle_nm':torque,'frames':frames}
    native.release(checkpoint)
finally:native.close()
zero=report['arms']['zero']['frames'][-1]
for arm in report['arms'].values():
    final=arm['frames'][-1]
    arm['delta_com_vx_from_zero_m_s']=final['com_velocity_ground_m_s'][0]-zero['com_velocity_ground_m_s'][0]
    arm['delta_ankle_r_from_zero_rad']=final['ankle_r']['value']-zero['ankle_r']['value']
(OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))
print(json.dumps({'output':str(OUT/'report.json'),'results':{label:{k:v for k,v in arm.items() if k!='frames'} for label,arm in report['arms'].items()}},indent=2))
