#!/usr/bin/env python3
"""Run finite muscle-only equilibrium waypoint programs; no prescribed motion."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.native.balance_observation import observe_balance

def run(program_path,output,*,keep_open=False,continuation=None):
    program=json.loads(Path(program_path).read_text());out=Path(output);out.mkdir(parents=True,exist_ok=False)
    bundle=ROOT/'data/models/engineering_stance_v1';reg=json.loads((bundle/'registration.json').read_text());base=NativeStanceLQR(bundle/'linearization.npz',model_sha256=reg['model_sha256'],target_mass_kg=reg['target_mass_kg'],dt_s=.01)
    sources={};segments=[];now=1. if continuation is None else continuation[0].state['time_s'];prev=base if continuation is None else continuation[1];initial_reference=prev
    def retain(path,name):
        raw=Path(path).read_bytes();(out/name).write_bytes(raw);sources[name]=hashlib.sha256(raw).hexdigest()
    retain(program_path,'program.json');retain(__file__,'evaluator.py');retain(bundle/'linearization.npz','base_policy.npz')
    for i,entry in enumerate(program['segments']):
        duration=float(entry['duration_s']);hold=float(entry.get('hold_s',0))
        if not 0<duration<=15 or not 0<=hold<=10:raise ValueError('Invalid segment timing')
        target=NativeStanceLQR(entry['policy'],target_mass_kg=base.target_mass_kg,dt_s=.01)
        if target.state_names!=base.state_names or target.muscle_names!=base.muscle_names:raise ValueError('Waypoint catalog mismatch')
        # Target model physical equivalence is documented by the owned endpoint manifest.
        proof=json.loads(Path(entry['target_manifest']).read_text())
        if proof.get('same_physical_model_except_initial_activation_defaults') is not True:raise ValueError('Endpoint changes physical model parameters')
        if proof.get('source_model_sha256')!=base.model_sha256 or proof.get('target_model_sha256')!=target.model_sha256:raise ValueError('Endpoint model proof identity mismatch')
        target_path=Path(entry['target_manifest']).parent/'target.npz'
        if hashlib.sha256(target_path.read_bytes()).hexdigest()!=proof['target_sha256']:raise ValueError('Endpoint target bytes changed')
        with np.load(target_path,allow_pickle=False) as target_state:
            if not np.array_equal(target_state['x0'],target.x0) or not np.array_equal(target_state['u0'],target.u0):raise ValueError('Endpoint policy and target state differ')
        retain(target_path,f'target_{i}.npz')
        retain(entry['policy'],f'policy_{i}.npz');retain(entry['target_manifest'],f'target_manifest_{i}.json')
        segments.append((now,now+duration,now+duration+hold,prev,target,entry['name']));now+=duration+hold;prev=target
    if not segments or now-(0. if continuation is None else continuation[0].state['time_s'])>40:raise ValueError('Finite nonempty program <=40s required')
    stream=continuation[0] if continuation is not None else NativeMechanicalStream(ROOT,out/'plant',environment='upright',target_mass_kg=base.target_mass_kg,augmented_registration=str((bundle/'registration.json').relative_to(ROOT)),initial_pose=reg['initial_pose'])
    frames=[];error=None;started=time.monotonic();value_indices={p.rsplit('/',1)[0]:i for i,p in enumerate(base.state_names) if p.endswith('/value')}
    try:
        while stream.state['time_s']<now-1e-9:
            t=stream.state['time_s'];reference=initial_reference.x0.copy();u=initial_reference.u0;K=initial_reference.K;name='initial_stance';blend=0.
            if t>=segments[0][0]:
                segment=next((s for s in segments if t<s[2]-1e-9),segments[-1]);start,end,stop,old,target,name=segment
                phase=float(np.clip((t-start)/(end-start),0,1));blend=10*phase**3-15*phase**4+6*phase**5;rate=(30*phase**2-60*phase**3+30*phase**4)/(end-start)
                delta=target.x0-old.x0;reference=old.x0+blend*delta;u=old.u0+blend*(target.u0-old.u0);K=old.K+blend*(target.K-old.K)
                for i,path in enumerate(base.state_names):
                    if path.endswith('/speed'):reference[i]+=rate*delta[value_indices[path.rsplit('/',1)[0]]]
            raw=u-K@(base.state_vector(stream.state)-reference);commands=dict(zip(base.muscle_names,map(float,np.clip(raw,0,1))));state=stream.advance(.01,actuation=commands);balance=observe_balance(state)
            clearance={side:min(c['center_m'][1]-c['radius_m'] for c in state['contacts'] if c.get('body_frame') in ('calcn_'+side,'toes_'+side)) for side in ('r','l')}
            centers={side:np.mean([c['center_m'] for c in state['contacts'] if c.get('body_frame') in ('calcn_'+side,'toes_'+side)],axis=0).tolist() for side in ('r','l')}
            f={'time_s':state['time_s'],'segment':name,'blend':blend,'balance':balance,'per_foot_clearance_m':clearance,'per_foot_centroid_ground_m':centers,'coordinates':state['coordinates'],'muscles':state['muscles'],'commands':commands,'clipped_count':int(sum((raw<0)|(raw>1))),'state_reference_error_norm':float(np.linalg.norm(base.state_vector(state)-reference))};frames.append(f)
            if len(frames)%100==0:print(json.dumps({k:f[k] for k in ['time_s','segment','blend','per_foot_clearance_m','clipped_count','state_reference_error_norm']}),flush=True)
            q=state['coordinates']
            if q['pelvis_ty']['value']<.6 or any(abs(q[k]['value'])>1 for k in ('pelvis_tilt','pelvis_list')):break
    except Exception as exc:error=type(exc).__name__+': '+str(exc)
    finally:
        (out/'final_snapshot.json').write_text(json.dumps(stream.snapshot()))
        if not keep_open:stream.close()
    final=frames[-1] if frames else None
    report={'elapsed_s':final['time_s'] if final else 0,'native_stream_identity':stream.identity,'continued_existing_native_process':continuation is not None,'completed_horizon':bool(final and final['time_s']>=now-1e-9),'horizon_s':now,'error':error,'final_com_speed_m_s':float(np.linalg.norm(final['balance']['com_velocity_ground_m_s'])) if final else None,'final_per_foot_clearance_m':final['per_foot_clearance_m'] if final else None,'final_per_foot_centroid_ground_m':final['per_foot_centroid_ground_m'] if final else None,'max_clipped_count':max((f['clipped_count'] for f in frames),default=0),'source_sha256':sources,'wall_seconds':time.monotonic()-started,'walking_demonstrated':False,'scope':'Finite desired equilibrium and gain sequence applied through native muscle excitation only. No prescribed coordinates, root forces, residual actuators, or external GRFs.'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');(out/'trajectory.json').write_text(json.dumps(frames));print(json.dumps(report,indent=2),flush=True);return (report,(stream,segments[-1][4])) if keep_open else report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--program',required=True);p.add_argument('--output',required=True);p.add_argument('--keep-open',action='store_true');a=p.parse_args()
    if not a.keep_open:run(a.program,a.output)
    else:
        receipt,live=run(a.program,a.output,keep_open=True)
        try:
            while True:
                print('READY_FOR_CONTINUATION',flush=True);line=sys.stdin.readline()
                if not line:break
                command=json.loads(line)
                if command.get('close'):break
                receipt,live=run(command['program'],command['output'],keep_open=True,continuation=live)
        finally:live[0].close()
