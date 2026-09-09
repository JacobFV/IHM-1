#!/usr/bin/env python3
"""Evaluate muscle-only feedback against zero and constant excitation controls."""
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.postural_control import PosturalController,PosturalConfig,posture_metrics

def evaluate(root,out,seconds=1.,config=PosturalConfig(),baselines=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    report={'schema':'ihm.native-postural-evaluation.v1','seconds_requested':seconds,'arms':{},'walking_demonstrated':False,'brain_trained':False}
    for arm in ('feedback','constant','zero'):
        stream=NativeMechanicalStream(root,out/arm,environment='upright',target_mass_kg=70)
        initial=stream.snapshot();policy=PosturalController(initial,config,baselines)
        frames=[];error=None;start=time.monotonic()
        try:
            while stream.state['time_s']<seconds-1e-9:
                commands=policy.commands(stream.state) if arm=='feedback' else {k:(0 if arm=='zero' else policy.baselines.get(k,config.baseline)) for k in initial['muscles']}
                state=stream.advance(min(.01,seconds-stream.state['time_s']),actuation=commands)
                frames.append({'time_s':state['time_s'],'coordinates':state['coordinates'],'contact_forces':state.get('contacts',[]),'commands':commands})
                if state['coordinates']['pelvis_ty']['value']<.6 or any(abs(state['coordinates'][k]['value'])>1 for k in ('pelvis_tilt','pelvis_list')):break
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        finally:
            final=stream.snapshot();stream.close()
        metrics=posture_metrics(initial,final)
        metrics.update(completed_horizon=metrics['elapsed_s']>=seconds-1e-8,fall_detected=metrics['pelvis_height_m']<.6 or abs(metrics['pelvis_tilt_rad'])>1 or abs(metrics['pelvis_list_rad'])>1,error=error,wall_seconds=time.monotonic()-start)
        report['arms'][arm]=metrics
        (out/arm/'trajectory.json').write_text(json.dumps(frames,allow_nan=False))
        (out/arm/'policy.json').write_text(json.dumps(policy.identity(),indent=2))
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--seconds',type=float,default=1);p.add_argument('--baseline',type=float,default=.03);p.add_argument('--length-gain',type=float,default=2);p.add_argument('--velocity-gain',type=float,default=.1);a=p.parse_args()
    if not 0<a.seconds<=30:p.error('seconds must be in (0,30]')
    print(json.dumps(evaluate(Path(__file__).resolve().parents[1],a.output,a.seconds,PosturalConfig(a.baseline,a.length_gain,a.velocity_gain)),indent=2))
