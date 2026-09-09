#!/usr/bin/env python3
"""Nonlinear native verification of local engineering LQR; no brain claim."""
from pathlib import Path
import argparse,json,sys,hashlib,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.native.postural_control import posture_metrics
from ihm.native.moment_arm_control import native_center_of_mass
from ihm.native.balance_observation import observe_balance
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def evaluate(output,artifact,registration,seconds=8.,arm='feedback',perturb=(0.,0.),gain=1.,target_mass_kg=70.,command_delay_ticks=0):
    if type(command_delay_ticks) is not int or not 0<=command_delay_ticks<=20:raise ValueError('Command delay must be integer ticks in [0,20]')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    registration_identity=json.loads((ROOT/registration).read_text())
    model=NativeStanceLQR(Path(artifact),model_sha256=registration_identity['model_sha256'],dt_s=.01);sources={}
    if model.target_mass_kg is not None and not np.isclose(model.target_mass_kg,target_mass_kg,rtol=0,atol=1e-8):raise ValueError('LQR target mass differs from native trial')
    for path in (Path(artifact),Path(__file__),ROOT/'ihm/native/stance_lqr.py'):
        raw=path.read_bytes();(out/path.name).write_bytes(raw);sources[str(path.resolve().relative_to(ROOT))]=hashlib.sha256(raw).hexdigest()
    stream=NativeMechanicalStream(ROOT,out/'plant',environment='upright',target_mass_kg=target_mass_kg,augmented_registration=registration,initial_pose=registration_identity.get('initial_pose'))
    initial=stream.snapshot();frames=[];error=None;started=time.monotonic()
    command_queue=[dict(zip(model.muscle_names,map(float,model.u0))) for _ in range(command_delay_ticks)]
    try:
        while stream.state['time_s']<seconds-1e-9:
            commands,diag=model.commands(stream.state,gain=gain)
            if arm=='held':commands=dict(zip(model.muscle_names,map(float,model.u0)))
            elif arm=='zero':commands=dict.fromkeys(commands,0.)
            if arm!='feedback':diag={'motor_owner':arm,'inactive_feedback_diagnostics':diag}
            requested_commands=dict(commands)
            if command_delay_ticks:
                command_queue.append(commands);commands=command_queue.pop(0)
            diag['command_delay_ticks']=command_delay_ticks
            forces=[]
            if 1.-1e-9<=stream.state['time_s']<1.1-1e-9 and any(perturb):forces=[{'body':'pelvis','point_m':stream.body_point(body='pelvis',station_m=initial['bodies']['pelvis']['mass_center_local_m'])['point_source_m'],'force_n':[perturb[0],0.,perturb[1]]}]
            state=stream.advance(min(.01,seconds-stream.state['time_s']),actuation=commands,forces=forces)
            com,com_velocity=native_center_of_mass(state)
            frames.append({'balance':observe_balance(state),'muscles':state['muscles'],'com_ground_m':com.tolist(),'com_velocity_ground_m_s':com_velocity.tolist(),'time_s':state['time_s'],'coordinates':state['coordinates'],'commands':commands,'requested_commands':requested_commands,'diagnostics':diag,'external_forces':forces,'contacts':state['contacts']})
            q=state['coordinates']
            if len(frames)%50==0:print(json.dumps({'arm':arm,'t':state['time_s'],'height':q['pelvis_ty']['value'],'x':q['pelvis_tx']['value'],'z':q['pelvis_tz']['value'],'tilt':q['pelvis_tilt']['value'],'diagnostics':diag}),flush=True)
            if q['pelvis_ty']['value']<.6 or any(abs(q[k]['value'])>1 for k in ('pelvis_tilt','pelvis_list')):break
    except Exception as exc:error=type(exc).__name__+': '+str(exc)
    finally:final=stream.snapshot();stream.close()
    initial_com,_=native_center_of_mass(initial)
    final_com,final_velocity=native_center_of_mass(final)
    peak_com=max([float(np.linalg.norm(np.asarray(f['com_ground_m'])-initial_com)) for f in frames],default=0.)
    report=posture_metrics(initial,final);report.update(final_com_displacement_m=(final_com-initial_com).tolist(),final_com_speed_m_s=float(np.linalg.norm(final_velocity)),peak_com_displacement_m=peak_com,arm=arm,gain=gain,seconds_requested=seconds,completed_horizon=final['time_s']>=seconds-1e-9,error=error,wall_seconds=time.monotonic()-started,source_sha256=sources,walking_demonstrated=False,brain_trained=False,target_mass_kg=target_mass_kg,command_delay_ticks=command_delay_ticks,perturbation={'start_s':1.,'duration_s':.1,'force_n':[perturb[0],0.,perturb[1]],'point_basis':'current pelvis COM in source ground; recomputed each step'})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');(out/'trajectory.json').write_text(json.dumps(frames,allow_nan=False));return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--artifact',default='data/research/locomotion_control/linearization_ftxx3psl/linearization.npz');p.add_argument('--registration',default='data/derived/mechanics/resting_stance86_refined/registration.json');p.add_argument('--command-delay-ticks',type=int,default=0);p.add_argument('--target-mass-kg',type=float,default=70.);p.add_argument('--gain',type=float,default=1);p.add_argument('--seconds',type=float,default=8);p.add_argument('--arm',choices=['feedback','held','zero'],default='feedback');p.add_argument('--perturb-x-n',type=float,default=0);p.add_argument('--perturb-z-n',type=float,default=0);a=p.parse_args()
    if not 0<a.seconds<=30:p.error('horizon must be in(0,30]')
    print(json.dumps(evaluate(a.output,a.artifact,a.registration,a.seconds,a.arm,(a.perturb_x_n,a.perturb_z_n),a.gain,a.target_mass_kg,a.command_delay_ticks),indent=2))
