#!/usr/bin/env python3
"""Bounded derivative-free engineering controller search on native dynamics.

No IBM checkpoint is read or written. Every candidate uses the same opaque
native initial state, and only muscle excitations act on the articulated body.
"""
import argparse,json,sys,time,hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.postural_control import PosturalConfig,PosturalController,posture_metrics

GROUPS=('vas','glmax','glmed','soleus','tibant','hamstring','hipflexor','other')
def group(name):
    if name.startswith(('semimem','semiten','bflh','bfsh')):return 'hamstring'
    if name.startswith(('psoas','iliacus','recfem')):return 'hipflexor'
    for g in GROUPS[:5]:
        if name.startswith(g):return g
    return 'other'
def decode(x,names):
    config=PosturalConfig(.03,float(np.exp(x[8])),float(np.exp(x[9])))
    baseline={n:float(1/(1+np.exp(-x[GROUPS.index(group(n))]))) for n in names}
    return config,baseline

def train(root,out,trials=32,seconds=1.2,seed=7):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(seed);stream=NativeMechanicalStream(root,out/'plant',environment='upright',target_mass_kg=70)
    initial=stream.snapshot();checkpoint=stream.checkpoint();names=list(initial['muscles'])
    best=np.r_[np.full(8,-3.5),np.log(2),np.log(.1)];best_score=-np.inf;history=[]
    try:
        for trial in range(trials):
            if trial==0:x=best.copy()
            else:x=best+rng.normal(0,max(.2,1.2*(1-trial/trials)),10)
            x[:8]=np.clip(x[:8],-6,1);x[8]=np.clip(x[8],-4,np.log(50));x[9]=np.clip(x[9],-6,np.log(3))
            config,baselines=decode(x,names);stream.restore(checkpoint);policy=PosturalController(initial,config,baselines)
            score=0.;error=None;frames=[];started=time.monotonic()
            try:
                while stream.state['time_s']<seconds-1e-9:
                    state=stream.advance(min(.01,seconds-stream.state['time_s']),actuation=policy.commands(stream.state));q=state['coordinates']
                    height=q['pelvis_ty']['value'];tilt=q['pelvis_tilt']['value'];lean=q['pelvis_list']['value']
                    score+=.01*(1+height-2*tilt*tilt-2*lean*lean-.1*q['pelvis_tx']['speed']**2)
                    frames.append({'time_s':state['time_s'],'coordinates':q})
                    if height<.6 or abs(tilt)>1 or abs(lean)>1:break
            except Exception as exc:error=str(exc);score-=2
            record={'trial':trial,'score':score,'metrics':posture_metrics(initial,stream.state),'error':error,'parameters':x.tolist(),'wall_seconds':time.monotonic()-started}
            history.append(record)
            if score>best_score:
                best_score=score;best=x.copy();(out/'best_policy.json').write_text(json.dumps(policy.identity(),indent=2));(out/'best_trajectory.json').write_text(json.dumps(frames))
            (out/'history.json').write_text(json.dumps(history,indent=2,allow_nan=False));print(json.dumps(record),flush=True)
    finally:
        try:
            if not stream.closed:stream.release(checkpoint)
        finally:stream.close()
    receipt={'schema':'ihm.native-posture-search.v1','seed':seed,'trials':trials,'horizon_s':seconds,'best_score':best_score,'best_parameters':best.tolist(),'brain_trained':False,'walking_demonstrated':False,'objective':'time integral of survival + pelvis height - squared pelvis tilt/list - forward speed penalty; early termination below 0.6m or rotation beyond 1 rad','source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'ihm/native/postural_control.py',root/'scripts/train_native_posture.py')}}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--trials',type=int,default=32);p.add_argument('--seconds',type=float,default=1.2);p.add_argument('--seed',type=int,default=7);a=p.parse_args()
    if not 1<=a.trials<=1000 or not 0<a.seconds<=30:p.error('invalid bounded training horizon/trials')
    train(Path(__file__).resolve().parents[1],a.output,a.trials,a.seconds,a.seed)
