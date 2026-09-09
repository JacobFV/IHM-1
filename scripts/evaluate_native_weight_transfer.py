#!/usr/bin/env python3
"""Test muscle-only interpolation between two native equilibrium states."""
import argparse,hashlib,json,time,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.native.balance_observation import observe_balance

def run(output,target,transition_s=5.,seconds=8.,fraction=1.,target_policy=None):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    bundle=ROOT/'data/models/engineering_stance_v1';reg=json.loads((bundle/'registration.json').read_text());policy=NativeStanceLQR(bundle/'linearization.npz',model_sha256=reg['model_sha256'],target_mass_kg=reg['target_mass_kg'],dt_s=.01)
    with np.load(target,allow_pickle=False) as t:
        if t['state_names'].tolist()!=policy.state_names or t['muscle_names'].tolist()!=policy.muscle_names or abs(float(t['target_mass_kg'])-policy.target_mass_kg)>1e-8:raise ValueError('Target state catalog or mass mismatch')
        if str(t['source_model_sha256'].item())!=policy.model_sha256:raise ValueError('Target source plant identity mismatch')
        xgoal=t['x_target'].copy();ugoal=t['u_target'].copy()
        if xgoal.shape!=policy.x0.shape or ugoal.shape!=policy.u0.shape or not np.isfinite(xgoal).all() or not np.isfinite(ugoal).all() or np.any((ugoal<0)|(ugoal>1)):raise ValueError('Invalid target state or excitation')
    goal_gain=policy.K
    if target_policy is not None:
        goal=NativeStanceLQR(target_policy,target_mass_kg=policy.target_mass_kg,dt_s=.01)
        if goal.state_names!=policy.state_names or goal.muscle_names!=policy.muscle_names or not np.array_equal(goal.x0,xgoal) or not np.array_equal(goal.u0,ugoal):raise ValueError('Target gain equilibrium mismatch')
        goal_gain=goal.K
    files={}
    for name,path in [('target.npz',Path(target)),('policy.npz',bundle/'linearization.npz'),('evaluator.py',Path(__file__))]:
        b=path.read_bytes();(out/name).write_bytes(b);files[name]=hashlib.sha256(b).hexdigest()
    if target_policy is not None:
        b=Path(target_policy).read_bytes();(out/'target_policy.npz').write_bytes(b);files['target_policy.npz']=hashlib.sha256(b).hexdigest()
    stream=NativeMechanicalStream(ROOT,out/'plant',environment='upright',target_mass_kg=policy.target_mass_kg,augmented_registration=str((bundle/'registration.json').relative_to(ROOT)),initial_pose=reg['initial_pose'])
    frames=[];error=None;started=time.monotonic();delta=xgoal-policy.x0;du=ugoal-policy.u0
    values={path.rsplit('/',1)[0]:i for i,path in enumerate(policy.state_names) if path.endswith('/value')}
    try:
        while stream.state['time_s']<seconds-1e-9:
            t=stream.state['time_s'];phase=np.clip((t-1)/transition_s,0,1);blend=fraction*(10*phase**3-15*phase**4+6*phase**5);rate=fraction*(30*phase**2-60*phase**3+30*phase**4)/transition_s
            reference=policy.x0+blend*delta
            for i,path in enumerate(policy.state_names):
                if path.endswith('/speed'):reference[i]+=rate*delta[values[path.rsplit('/',1)[0]]]
            gain=policy.K+blend*(goal_gain-policy.K)
            raw=policy.u0+blend*du-gain@(policy.state_vector(stream.state)-reference);commands=dict(zip(policy.muscle_names,map(float,np.clip(raw,0,1))))
            state=stream.advance(.01,actuation=commands);balance=observe_balance(state);forces=balance['per_foot_normal_force_n'];total=sum(forces.values());left=forces['l']/total if total>1 else None
            clearance={side:min(c['center_m'][1]-c['radius_m'] for c in state['contacts'] if c.get('body_frame') in ('calcn_'+side,'toes_'+side)) for side in ('r','l')}
            frame={'time_s':state['time_s'],'per_foot_clearance_m':clearance,'blend':blend,'left_load_fraction':left,'balance':balance,'coordinates':state['coordinates'],'muscles':state['muscles'],'commands':commands,'clipped_count':int(np.sum((raw<0)|(raw>1))),'native_floor_violation_count':int(np.sum(raw<policy.minimum_activation)),'state_reference_error_norm':float(np.linalg.norm(policy.state_vector(state)-reference))};frames.append(frame)
            if len(frames)%50==0:print(json.dumps({k:frame[k] for k in ['time_s','blend','left_load_fraction','clipped_count','state_reference_error_norm']}),flush=True)
            q=state['coordinates']
            if q['pelvis_ty']['value']<.6 or any(abs(q[k]['value'])>1 for k in ('pelvis_tilt','pelvis_list')):break
    except Exception as exc:error=type(exc).__name__+': '+str(exc)
    finally:
        (out/'final_snapshot.json').write_text(json.dumps(stream.snapshot()))
        stream.close()
    report={'elapsed_s':frames[-1]['time_s'] if frames else 0,'completed_horizon':bool(frames and frames[-1]['time_s']>=seconds-1e-9),'final_left_load_fraction':frames[-1]['left_load_fraction'] if frames else None,'final_com_speed_m_s':float(np.linalg.norm(frames[-1]['balance']['com_velocity_ground_m_s'])) if frames else None,'final_state_reference_error_norm':frames[-1]['state_reference_error_norm'] if frames else None,'max_clipped_count':max((f['clipped_count'] for f in frames),default=0),'fraction':fraction,'gain_scheduled':target_policy is not None,'transition_s':transition_s,'seconds_requested':seconds,'error':error,'wall_seconds':time.monotonic()-started,'source_sha256':files,'walking_demonstrated':False,'scope':'Desired full equilibrium states and excitations are interpolated only in feedback reference; plant coordinates freely integrate. No root forces or prescribed motion. Endpoint LQR gains optionally interpolate; no guarantee along interpolation.'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');(out/'trajectory.json').write_text(json.dumps(frames));print(json.dumps(report,indent=2));return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);p.add_argument('--target',default='data/derived/mechanics/patient_left80_stance98/lqr_target/target.npz');p.add_argument('--transition-s',type=float,default=5.);p.add_argument('--seconds',type=float,default=8.);p.add_argument('--target-policy');p.add_argument('--fraction',type=float,default=1.);a=p.parse_args()
    if not 0<a.transition_s<=20 or not 0<a.seconds<=30 or not 0<a.fraction<=1:p.error('Invalid bounded horizon or interpolation fraction')
    run(a.output,a.target,a.transition_s,a.seconds,a.fraction,a.target_policy)
