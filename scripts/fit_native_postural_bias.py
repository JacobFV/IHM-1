#!/usr/bin/env python3
"""Identify bounded feedforward joint torques from native checkpoint perturbations.

This fits local initial acceleration cancellation, not gait or a prescribed
trajectory. Every observation is a muscle-excitation-driven native rollout.
"""
from pathlib import Path
import argparse,json,sys,tempfile,time
import numpy as np
from scipy.optimize import lsq_linear
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import JointPosturalController,native_center_of_mass
ROOT=Path(__file__).resolve().parents[1]
GROUPS=[('hip_flexion',1),('hip_adduction',-1),('hip_rotation',-1),('knee_angle',1),('ankle_angle',1),('lumbar_extension',0)]

def fit(output,seconds=.05,iterations=2):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    native=NativeMechanicalStream(ROOT,output/'native',environment='upright',target_mass_kg=70,augmented_registration='data/derived/mechanics/initial_stance86/registration.json')
    initial=native.snapshot();checkpoint=native.checkpoint()
    lumbar=[n for n in initial['muscles'] if n.startswith('gait2392_')]
    joints=['lumbar_extension','lumbar_bending','lumbar_rotation']
    initial_arms=native.moment_arms(muscles=lumbar,coordinates=joints)['moment_arms_m']
    evaluations=[];history=[]
    def torque_map(vector):
        result={}
        for (name,mirror),value in zip(GROUPS,vector):
            if mirror:
                result[name+'_r']=float(value);result[name+'_l']=float(mirror*value)
            else:result[name]=float(value)
        return result
    def evaluate(vector,label):
        native.restore(checkpoint)
        policy=JointPosturalController(initial,kp=150,kd=25,baseline=.02,pelvis_gain=0,pelvis_damping=0,com_position_gain=0,com_velocity_gain=0,native_moment_arms=initial_arms,feedforward_torques=torque_map(vector))
        frames=[];start=time.monotonic()
        for _ in range(round(seconds/.01)):
            policy.update_native_moment_arms(native.moment_arms(muscles=lumbar,coordinates=joints)['moment_arms_m'])
            commands=policy.commands(native.state)
            state=native.advance(.01,actuation=commands)
            frames.append({'time_s':state['time_s'],'coordinates':state['coordinates'],'allocation':policy.last_allocation,'commands':commands})
        com,velocity=native_center_of_mass(state)
        # All actuated limb/lumbar coordinate speeds remain separate, preserving
        # asymmetric responses even though requested bias uses bilateral symmetry.
        names=[n for name,mirror in GROUPS for n in ([name+'_r',name+'_l'] if mirror else [name])]+['lumbar_bending','lumbar_rotation']
        residual=[state['coordinates'][n]['speed'] for n in names]+[2*velocity[0],2*velocity[2],2*velocity[1]]
        record={'label':label,'bias_vector_nm':list(map(float,vector)),'feedforward_torques_nm':torque_map(vector),'residual_names':names+['2*com_vx','2*com_vz','2*com_vy'],'residual':list(map(float,residual)),'residual_norm':float(np.linalg.norm(residual)),'com_ground_m':com.tolist(),'com_velocity_ground_m_s':velocity.tolist(),'frames':frames,'wall_seconds':time.monotonic()-start}
        index=len(evaluations);path=output/f'trial_{index:03d}.json';path.write_text(json.dumps(record,allow_nan=False))
        summary={k:v for k,v in record.items() if k!='frames'};summary['artifact']=str(path.relative_to(ROOT));evaluations.append(summary)
        print(json.dumps({'trial':index,'label':label,'bias':record['bias_vector_nm'],'norm':record['residual_norm'],'wall_seconds':record['wall_seconds']}),flush=True)
        return np.array(residual)
    bias=np.zeros(6)
    try:
        residual=evaluate(bias,'zero_bias')
        for iteration in range(iterations):
            jacobian=np.zeros((len(residual),6));delta=5.
            for j in range(6):
                probe=bias.copy();probe[j]+=delta
                jacobian[:,j]=(evaluate(probe,f'iteration{iteration}_jacobian{j}')-residual)/delta
            regularization=.001
            matrix=np.vstack((jacobian,regularization*np.eye(6)))
            target=np.r_[-residual,np.zeros(6)]
            step=lsq_linear(matrix,target,bounds=(np.maximum(-50.,-150.-bias),np.minimum(50.,150.-bias)),tol=1e-8,max_iter=200).x
            candidate=np.clip(bias+step,-150,150);new_residual=evaluate(candidate,f'iteration{iteration}_candidate')
            accepted=np.linalg.norm(new_residual)<np.linalg.norm(residual)
            history.append({'iteration':iteration,'jacobian':jacobian.tolist(),'proposed_step_nm':step.tolist(),'accepted':bool(accepted),'before_norm':float(np.linalg.norm(residual)),'after_norm':float(np.linalg.norm(new_residual))})
            if accepted:bias,residual=candidate,new_residual
            else:
                candidate=np.clip(bias+.5*step,-150,150);new_residual=evaluate(candidate,f'iteration{iteration}_halfstep')
                if np.linalg.norm(new_residual)<np.linalg.norm(residual):bias,residual=candidate,new_residual
        report={'schema':'ihm.native-postural-bias-identification.v1','horizon_s':seconds,'step_s':.01,'iterations':iterations,'groups':GROUPS,'baseline_policy':{'kp':150,'kd':25,'baseline':.02,'COM_gains':0,'vestibular_gains':0},'feedforward_torques_nm':torque_map(bias),'final_residual_norm':float(np.linalg.norm(residual)),'zero_bias_residual_norm':evaluations[0]['residual_norm'],'evaluations':evaluations,'history':history,'scope':'Local native initial acceleration compensation only; no sustained stance or walking acceptance.','native_matched_checkpoint_each_trial':True,'external_support':False,'prescribed_motion':False}
        (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        native.release(checkpoint)
        print(json.dumps({'report':str(output/'report.json'),'feedforward_torques_nm':torque_map(bias),'final_norm':float(np.linalg.norm(residual)),'initial_norm':evaluations[0]['residual_norm']}),flush=True)
        return report
    finally:native.close()
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output');parser.add_argument('--seconds',type=float,default=.05);parser.add_argument('--iterations',type=int,default=2);args=parser.parse_args()
    if not .01<=args.seconds<=.2 or abs(args.seconds/.01-round(args.seconds/.01))>1e-9:parser.error('seconds must align to10ms and be.01..2')
    if not 1<=args.iterations<=5:parser.error('iterations must be1..5')
    output=Path(args.output) if args.output else Path(tempfile.mkdtemp(prefix='native_bias_fit_',dir=ROOT/'data/research/locomotion_control'))/'fit'
    fit(output,args.seconds,args.iterations)
