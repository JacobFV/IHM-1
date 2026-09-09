"""Matched nonlinear native full/sever stance check; sever retains same tonic drive."""
from pathlib import Path
import argparse,json,sys,hashlib,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from ihm.native.cortical_stance import load_cortical_stance,PersistentStanceCommands
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import native_center_of_mass
from ihm.native.postural_control import posture_metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--artifact',required=True);p.add_argument('--output',required=True)
    p.add_argument('--registration',default='data/derived/mechanics/resting_stance98/registration.json')
    p.add_argument('--seconds',type=float,default=2.)
    p.add_argument('--target-mass',type=float)
    p.add_argument('--perturb-x',type=float,default=5.)
    p.add_argument('--perturb-z',type=float,default=5.)
    args=p.parse_args()
    if not .02<=args.seconds<=10 or abs(args.seconds/.01-round(args.seconds/.01))>1e-9:raise ValueError('Bounded10ms horizon required')
    torch.set_num_threads(1);out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    registration=json.loads(Path(args.registration).read_text())
    policy,artifact=load_cortical_stance(args.artifact,model_sha256=registration['model_sha256'],dt_s=.01)
    mass=artifact['provenance'].get('target_mass_kg',70.)
    if args.target_mass is not None and abs(args.target_mass-mass)>1e-8:raise ValueError('Native mass differs from trained artifact')
    if not np.isfinite(args.perturb_x) or not np.isfinite(args.perturb_z):raise ValueError('Finite force required')
    if artifact['provenance']['native_activation_floor']!=.01:raise ValueError('Unexpected native activation floor')
    retained=out/'sources';retained.mkdir()
    sources={}
    for path in (Path(__file__),ROOT/'ihm/native/cortical_stance.py',Path(args.artifact),Path(args.registration)):
        raw=path.read_bytes();(retained/path.name).write_bytes(raw);sources[str(path.resolve())]=hashlib.sha256(raw).hexdigest()
    stream=NativeMechanicalStream(ROOT,out/'plant',environment='upright',target_mass_kg=mass,augmented_registration=args.registration,initial_pose=registration.get('initial_pose'))
    initial=stream.snapshot();token=stream.checkpoint();initial_com,_=native_center_of_mass(initial)
    records={};started=time.monotonic()
    try:
        for arm in ('full','sever'):
            if stream.closed:break
            stream.restore(token);adapter=PersistentStanceCommands(policy,dt_s=.01);frames=[];error=None
            try:
                while stream.state['time_s']<args.seconds-1e-9:
                    commands=adapter.commands(stream.state,sever=arm=='sever')
                    if arm=='sever':assert np.allclose([commands[m] for m in policy.muscle_names],policy.u0.numpy(),atol=1e-8,rtol=0)
                    forces=[]
                    if 1.-1e-9<=stream.state['time_s']<1.1-1e-9:
                        point=stream.body_point(body='pelvis',station_m=initial['bodies']['pelvis']['mass_center_local_m'])['point_source_m']
                        forces=[{'body':'pelvis','point_m':point,'force_n':[args.perturb_x,0.,args.perturb_z]}]
                    state=stream.advance(.01,actuation=commands,forces=forces)
                    com,velocity=native_center_of_mass(state)
                    frames.append({'time_s':state['time_s'],'com_m':com.tolist(),'com_velocity_m_s':velocity.tolist(),
                        'coordinates':state['coordinates'],'muscles':state['muscles'],'commands':commands,'force_pulse':bool(forces),'commands_at_physical_floor':sum(v<=.01+1e-10 for v in commands.values())})
                    if len(frames)%50==0:print(json.dumps({'arm':arm,'time_s':state['time_s'],'com_displacement_m':float(np.linalg.norm(com-initial_com))}),flush=True)
                    if state['coordinates']['pelvis_ty']['value']<.6:break
            except Exception as exc:error=type(exc).__name__+': '+str(exc)
            final=stream.snapshot();com,velocity=native_center_of_mass(final)
            report=posture_metrics(initial,final)
            report.update(arm=arm,error=error,completed_horizon=final['time_s']>=args.seconds-1e-9,
                elapsed_s=final['time_s'],peak_com_displacement_m=max([float(np.linalg.norm(np.array(f['com_m'])-initial_com)) for f in frames],default=0.),
                final_com_displacement_m=float(np.linalg.norm(com-initial_com)),final_com_speed_m_s=float(np.linalg.norm(velocity)))
            records[arm]={'report':report,'frames':frames}
            (out/'partial.json').write_text(json.dumps(records)+'\n')
            print(json.dumps(report,indent=2),flush=True)
    finally:
        if not stream.closed:stream.release(token)
        stream.close()
    for arm in ('full','sever'):
        if arm not in records:records[arm]={'report':{'arm':arm,'error':'Not run: preceding native failure closed the matched checkpoint owner','completed_horizon':False,'peak_com_displacement_m':None},'frames':[]}
    full=records['full']['report'];sever=records['sever']['report']
    report={'schema':'ihm.native-cortical-stance-eval.v1','artifact_sha256':artifact['artifact_sha256'],
        'model_sha256':registration['model_sha256'],'target_mass_kg':mass,
        'perturbation_force_n':[args.perturb_x,0.,args.perturb_z],'seconds_requested':args.seconds,
        'full_improves_peak_com':full['completed_horizon'] and full['error'] is None and sever['peak_com_displacement_m'] is not None and full['peak_com_displacement_m']<sever['peak_com_displacement_m'],
        'arms':{k:v['report'] for k,v in records.items()},'source_sha256':sources,'wall_seconds':time.monotonic()-started,
        'scope':['Exact matched native checkpoint; same muscle excitation baseline in both arms',
            f'Full uses persistent{policy.dyn.n}site trained IBM E/I; sever zeroes association but retains baseline',
            f'PelvisCOM pulse{args.perturb_x},0,{args.perturb_z}N at1..1.1s; no extra stabilizer, PD bypass, pose prescription or artificial support',
            f'{mass}kg native mechanics only; no full physiological exchange; no walking claim']}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/'trajectory.json').write_text(json.dumps(records)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
