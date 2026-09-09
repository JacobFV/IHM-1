#!/usr/bin/env python3
"""Bounded slow donor joint tracking through native muscle commands only."""
from pathlib import Path
import argparse,hashlib,json,sys,tempfile,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.gait_tracking_probe import SlowJointReference
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.native.moment_arm_control import native_center_of_mass
ROOT=Path(__file__).resolve().parents[1]

def run(output,artifact,registration,*,seconds=6.,amplitude=.05,mass_kg=70.):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    reg=json.loads((ROOT/registration).read_text());model=NativeStanceLQR(artifact,model_sha256=reg['model_sha256'],dt_s=.01)
    if model.target_mass_kg is not None and abs(model.target_mass_kg-mass_kg)>1e-8:raise ValueError('Policy/native mass mismatch')
    reference=SlowJointReference(amplitude=amplitude,duration_s=seconds-1.)
    controlled_coordinates={p.strip('/').split('/')[-2] for p in model.state_names if p.startswith('/jointset/') and p.endswith('/value')}
    if set(reference.names)-controlled_coordinates:raise ValueError('Donor targets include coordinates outside independent LQR chart')
    source_paths=[Path(artifact),Path(__file__),ROOT/'ihm/native/gait_tracking_probe.py',ROOT/'ihm/native/gait_reference.py',ROOT/'ihm/native/stance_lqr.py',ROOT/registration]
    sources={}
    for i,path in enumerate(source_paths):
        data=path.read_bytes();sources[str(path.resolve().relative_to(ROOT))]=hashlib.sha256(data).hexdigest();(out/f'source_{i}_{path.name}').write_bytes(data)
    native=NativeMechanicalStream(ROOT,out/'plant',environment='upright',target_mass_kg=mass_kg,augmented_registration=registration,initial_pose=reg.get('initial_pose'))
    initial=native.snapshot();base_q={n:float(initial['coordinates'][n]['value']) for n in reference.names}
    # Artifact x0 is the controller reference; retain it exactly where available.
    for path,value in zip(model.state_names,model.x0):
        parts=path.strip('/').split('/')
        if parts[0]=='jointset' and parts[-1]=='value' and parts[-2] in base_q:base_q[parts[-2]]=float(value)
    frames=[];error=None;start=time.monotonic();peak_clip=0;peak_floor=0
    try:
        while native.state['time_s']<seconds-1e-9:
            before=native.snapshot();observed,target=reference.controller_observation(before,before['time_s'])
            commands,diag=model.commands(observed)
            assert native.state['coordinates']==before['coordinates'],'Reference helper mutated native state'
            state=native.advance(.01,actuation=commands)
            target_after=reference.sample(state['time_s']);com,velocity=native_center_of_mass(state)
            foot_forces={side:sum(max(0.,c['force_n'][1]) for c in state['contacts'] if c.get('body_frame') in ('calcn_'+side,'toes_'+side)) for side in ('r','l')}
            target_angles={n:base_q[n]+target_after['joint_offsets_rad'][n] for n in reference.names}
            errors={n:state['coordinates'][n]['value']-target_angles[n] for n in reference.names}
            peak_clip=max(peak_clip,diag['clipped_muscles']);peak_floor=max(peak_floor,diag.get('requested_below_minimum_activation_muscles') or 0)
            frames.append({'time_s':state['time_s'],'target':target_after,'target_angles_rad':target_angles,'tracking_error_rad':errors,'coordinates':state['coordinates'],'com_ground_m':com.tolist(),'com_velocity_ground_m_s':velocity.tolist(),'foot_vertical_force_n':foot_forces,'liftoff_under_5n':{s:foot_forces[s]<5 for s in ('r','l')},'commands':commands,'actual_activation':{n:m['activation'] for n,m in state['muscles'].items()},'diagnostics':diag})
            if len(frames)%50==0:print(json.dumps({'time_s':state['time_s'],'height':state['coordinates']['pelvis_ty']['value'],'tracking_rmse_rad':float(np.sqrt(np.mean(np.square(list(errors.values()))))),'clip':diag['clipped_muscles'],'below_floor':diag.get('requested_below_minimum_activation_muscles'),'feet_n':foot_forces}),flush=True)
            if state['coordinates']['pelvis_ty']['value']<.6 or any(abs(state['coordinates'][n]['value'])>1 for n in ('pelvis_tilt','pelvis_list')):break
    except Exception as exc:error=type(exc).__name__+': '+str(exc)
    finally:final=native.snapshot();native.close()
    active=[f for f in frames if f['time_s']>1.]
    actual_error=np.array([[f['tracking_error_rad'][n] for n in reference.names] for f in active])
    frozen_error=np.array([[f['target']['joint_offsets_rad'][n] for n in reference.names] for f in active])
    per_joint={n:{'tracking_rmse_rad':float(np.sqrt(np.mean(actual_error[:,i]**2))),'frozen_pose_rmse_rad':float(np.sqrt(np.mean(frozen_error[:,i]**2))),'maximum_target_offset_rad':float(np.max(np.abs(frozen_error[:,i]))),'actual_range_rad':float(np.ptp([f['coordinates'][n]['value'] for f in active]))} for i,n in enumerate(reference.names)} if active else {}
    initial_com,_=native_center_of_mass(initial);final_com,final_velocity=native_center_of_mass(final)
    report={'schema':'ihm.native-small-joint-tracking-probe.v1','registration':registration,'target_mass_kg':mass_kg,'amplitude':amplitude,'stance_leadin_s':1.,'seconds_requested':seconds,'elapsed_s':final['time_s'],'completed_horizon':final['time_s']>=seconds-1e-8,'error':error,'wall_seconds':time.monotonic()-start,'target_joint_count':len(reference.names),'source_sha256':sources,'donor_provenance':reference.reference.provenance,'target_tracking':per_joint,'tracking_rmse_rad':None if not active else float(np.sqrt(np.mean(actual_error**2))),'frozen_pose_rmse_rad':None if not active else float(np.sqrt(np.mean(frozen_error**2))),'peak_clipped_muscles':peak_clip,'peak_requested_below_native_activation_floor':peak_floor,'final_com_displacement_m':(final_com-initial_com).tolist(),'final_com_speed_m_s':float(np.linalg.norm(final_velocity)),'liftoff_frames_under_5n':{s:sum(f['liftoff_under_5n'][s] for f in frames) for s in ('r','l')},'walking_demonstrated':False,'brain_trained':False,'prescribed_motion':False,'external_forces_applied':False,'root_reference_changed':False,'source_CMC_excitations_used':False,'scope':'Small finite-recording joint-target tracking probe through muscle-only local LQR; no gait, periodicity, or walking acceptance.'}
    (out/'trajectory.json').write_text(json.dumps(frames,allow_nan=False));(out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'report':str(out/'report.json'),**{k:v for k,v in report.items() if k not in ('target_tracking','donor_provenance','source_sha256')}},indent=2),flush=True)
    return report
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output');parser.add_argument('--artifact',default='data/research/locomotion_control/linearization_eyxmbi_z/discrete_margin/linearization.npz');parser.add_argument('--registration',default='data/derived/mechanics/resting_stance98/registration.json');parser.add_argument('--seconds',type=float,default=6.);parser.add_argument('--amplitude',type=float,default=.05);parser.add_argument('--mass-kg',type=float,default=70.);args=parser.parse_args()
    if not 1.5<=args.seconds<=10 or abs(args.seconds/.01-round(args.seconds/.01))>1e-8:parser.error('Horizon must align to10ms in1.5..10s')
    out=Path(args.output) if args.output else Path(tempfile.mkdtemp(prefix='native_joint_tracking_',dir=ROOT/'data/research/locomotion_control'))/'probe'
    run(out,args.artifact,args.registration,seconds=args.seconds,amplitude=args.amplitude,mass_kg=args.mass_kg)
