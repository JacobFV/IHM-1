#!/usr/bin/env python3
"""Measure engineering moment-arm muscle control on freely integrated plant."""
import argparse,json,sys,time,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import JointPosturalController,foot_support_center
from ihm.native.postural_control import posture_metrics

def evaluate(root,out,seconds=2.,augmented_registration=None,com_target_support=False,perturb_x_n=0.,perturb_z_n=0.,**gains):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    sources={}
    for source in (Path(__file__),Path(root)/'ihm/native/moment_arm_control.py'):
        raw=source.read_bytes();(out/source.name).write_bytes(raw);sources[str(source.relative_to(root))]=hashlib.sha256(raw).hexdigest()
    stream=NativeMechanicalStream(root,out/'plant',environment='upright',target_mass_kg=70,augmented_registration=augmented_registration)
    try:
        initial=stream.snapshot();lumbar=[n for n in initial['muscles'] if n.startswith('gait2392_')]
        native_arms=stream.moment_arms(muscles=lumbar,coordinates=['lumbar_extension','lumbar_bending','lumbar_rotation'])['moment_arms_m'] if lumbar else {}
        if com_target_support:
            support=foot_support_center(initial)
            if support is None:raise ValueError('No initial native foot support')
            gains['com_target_x_m']=float(support[0])
        policy=JointPosturalController(initial,native_moment_arms=native_arms,**gains);frames=[];error=None;start=time.monotonic()
    except BaseException:
        stream.close()
        raise
    try:
        while stream.state['time_s']<seconds-1e-9:
            if lumbar and len(frames)%5==0:policy.update_native_moment_arms(stream.moment_arms(muscles=lumbar,coordinates=['lumbar_extension','lumbar_bending','lumbar_rotation'])['moment_arms_m'])
            commands=policy.commands(stream.state)
            dt=min(.01,seconds-stream.state['time_s'])
            forces=[]
            if 1.-1e-9<=stream.state['time_s']<1.1-1e-9 and (perturb_x_n or perturb_z_n):
                forces=[{'body':'pelvis','point_m':stream.body_point(body='pelvis',station_m=initial['bodies']['pelvis']['mass_center_local_m'])['point_source_m'],'force_n':[perturb_x_n,0.,perturb_z_n]}]
            state=stream.advance(dt,actuation=commands,forces=forces)
            frames.append({'time_s':state['time_s'],'coordinates':state['coordinates'],'commands':commands,'allocation':policy.last_allocation,'contacts':state['contacts'],'external_forces':forces})
            q=state['coordinates']
            if len(frames)%20==0:print(json.dumps({'time_s':state['time_s'],'height_m':q['pelvis_ty']['value'],'pelvis_tilt_rad':q['pelvis_tilt']['value'],'lumbar_extension_rad':q['lumbar_extension']['value'],'pelvis_x_m':q['pelvis_tx']['value'],'com':policy.last_allocation.get('com_feedback')}),flush=True)
            if q['pelvis_ty']['value']<.6 or any(abs(q[k]['value'])>1 for k in ('pelvis_tilt','pelvis_list')):break
    except Exception as exc:error=type(exc).__name__+': '+str(exc)
    finally:
        final=stream.snapshot();stream.close()
    report=posture_metrics(initial,final);report.update(perturbation={'start_s':1.,'duration_s':.1,'force_ground_n':[perturb_x_n,0.,perturb_z_n],'station':'pelvis local mass center'},source_sha256=sources,completed_horizon=final['time_s']>=seconds-1e-9,fall_detected=final['coordinates']['pelvis_ty']['value']<.6 or any(abs(final['coordinates'][k]['value'])>1 for k in ('pelvis_tilt','pelvis_list')),lumbar_extension_rad=final['coordinates']['lumbar_extension']['value'],error=error,wall_seconds=time.monotonic()-start,seconds_requested=seconds,walking_demonstrated=False,brain_trained=False,policy=policy.identity())
    (out/'report.json').write_text(json.dumps(report,indent=2));(out/'trajectory.json').write_text(json.dumps(frames,allow_nan=False))
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--augmented-registration');p.add_argument('--perturb-x-n',type=float,default=0);p.add_argument('--perturb-z-n',type=float,default=0);p.add_argument('--feedforward-json');p.add_argument('--equilibrium-json');p.add_argument('--seconds',type=float,default=2);p.add_argument('--kp',type=float,default=100);p.add_argument('--kd',type=float,default=15);p.add_argument('--com-target-support',action='store_true');p.add_argument('--com-lateral-position-gain',type=float,default=0);p.add_argument('--com-lateral-velocity-gain',type=float,default=0);p.add_argument('--com-position-gain',type=float,default=0);p.add_argument('--com-velocity-gain',type=float,default=0);p.add_argument('--pelvis-gain',type=float,default=60);p.add_argument('--pelvis-damping',type=float,default=10);a=p.parse_args()
    if not 0<a.seconds<=30:p.error('invalid horizon')
    equilibrium=None if a.equilibrium_json is None else json.loads(Path(a.equilibrium_json).read_text())
    feedforward={} if a.feedforward_json is None else json.loads(Path(a.feedforward_json).read_text())
    print(json.dumps(evaluate(Path(__file__).resolve().parents[1],a.output,a.seconds,perturb_x_n=a.perturb_x_n,perturb_z_n=a.perturb_z_n,equilibrium_excitations=equilibrium,feedforward_torques=feedforward,augmented_registration=a.augmented_registration,com_target_support=a.com_target_support,kp=a.kp,kd=a.kd,pelvis_gain=a.pelvis_gain,pelvis_damping=a.pelvis_damping,com_position_gain=a.com_position_gain,com_velocity_gain=a.com_velocity_gain,com_lateral_position_gain=a.com_lateral_position_gain,com_lateral_velocity_gain=a.com_lateral_velocity_gain),indent=2))
