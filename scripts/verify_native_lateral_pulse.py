#!/usr/bin/env python3
"""Matched native lateral hip torque-request causal sign test, not stability."""
from pathlib import Path
import argparse,hashlib,json,sys,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import JointPosturalController,native_center_of_mass
ROOT=Path(__file__).resolve().parents[1]
def run(output=None):
    out=Path(output).resolve() if output else Path(tempfile.mkdtemp(prefix='native_lateral_pulse_',dir=ROOT/'data/research/locomotion_control'))
    out.mkdir(parents=True,exist_ok=True)
    regpath=ROOT/'data/derived/mechanics/resting_stance86/registration.json';raw=regpath.read_bytes();reg=json.loads(raw)
    native=NativeMechanicalStream(ROOT,out/'native',environment='upright',target_mass_kg=70,augmented_registration=str(regpath.relative_to(ROOT)))
    report={'schema':'ihm.native-lateral-pulse.v1','seconds':.15,'step_s':.01,'kp':150,'kd':25,'arms':{},
        'registration_sha256':hashlib.sha256(raw).hexdigest(),'scope':'Paired hip-adduction torque requests on resting stance86; equilibrium excitation baseline; COM and vestibular feedback disabled; native muscle actuation only; no stability claim.'}
    sourcepaths=[Path(__file__),ROOT/'ihm/native/moment_arm_control.py'];sources={p.name:p.read_bytes() for p in sourcepaths}
    report['source_sha256']={name:hashlib.sha256(data).hexdigest() for name,data in sources.items()}
    checkpoint=None
    try:
        initial=native.snapshot();checkpoint=native.checkpoint();lumbar=[n for n in initial['muscles'] if n.startswith('gait2392_')]
        for label,torque in [('zero',0.),('positive',20.),('negative',-20.)]:
            native.restore(checkpoint)
            arms=native.moment_arms(muscles=lumbar,coordinates=['lumbar_extension','lumbar_bending','lumbar_rotation'])['moment_arms_m'] if lumbar else {}
            policy=JointPosturalController(initial,kp=150,kd=25,pelvis_gain=0,pelvis_damping=0,
                com_position_gain=0,com_velocity_gain=0,com_lateral_position_gain=0,com_lateral_velocity_gain=0,
                equilibrium_excitations=reg['equilibrium_excitations'],native_moment_arms=arms,
                feedforward_torques={'hip_adduction_r':torque,'hip_adduction_l':-torque})
            frames=[]
            for i in range(15):
                if i and i%5==0 and lumbar:policy.update_native_moment_arms(native.moment_arms(muscles=lumbar,coordinates=['lumbar_extension','lumbar_bending','lumbar_rotation'])['moment_arms_m'])
                commands=policy.commands(native.state);state=native.advance(.01,actuation=commands)
                com,velocity=native_center_of_mass(state)
                frames.append({'time_s':state['time_s'],'com_ground_m':com.tolist(),'com_velocity_ground_m_s':velocity.tolist(),
                    'pelvis_list':state['coordinates']['pelvis_list'],'hip_adduction_r':state['coordinates']['hip_adduction_r'],
                    'hip_adduction_l':state['coordinates']['hip_adduction_l'],'allocation':policy.last_allocation,'commands':commands})
            report['arms'][label]={'requested_right_nm':torque,'requested_left_nm':-torque,'frames':frames}
            (out/'partial_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            print(label,'completed',flush=True)
    finally:
        if checkpoint is not None:native.release(checkpoint)
        native.close()
    for label,arm in report['arms'].items():
        arm['comparisons']={}
        for index in (4,14):
            zero=report['arms']['zero']['frames'][index];frame=arm['frames'][index]
            arm['comparisons'][str(round(frame['time_s'],2))]={
                'delta_com_vz_from_zero_m_s':frame['com_velocity_ground_m_s'][2]-zero['com_velocity_ground_m_s'][2],
                'delta_com_z_from_zero_m':frame['com_ground_m'][2]-zero['com_ground_m'][2],
                **{'delta_'+name+'_from_zero_rad':frame[name]['value']-zero[name]['value'] for name in ('pelvis_list','hip_adduction_r','hip_adduction_l')}}
    for name,data in sources.items():(out/name).write_bytes(data)
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'output':str(out/'report.json'),'results':{label:arm['comparisons'] for label,arm in report['arms'].items()}},indent=2),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output');args=parser.parse_args();run(args.output)
