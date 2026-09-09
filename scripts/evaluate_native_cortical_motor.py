#!/usr/bin/env python3
"""Matched native held-out targets through persistent actual IBM E/I cortex."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.cortical_motor import load_cortical_policy,PersistentCorticalCommands
from ihm.native.mechanical_stream import NativeMechanicalStream

def evaluate(root,out,artifact,seconds=.3):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1)
    source_paths=[Path(__file__),root/'ihm/native/cortical_motor.py',root/'ihm/native/mechanical_stream.py']
    sources={str(p.relative_to(root)):p.read_bytes() for p in source_paths}
    policy,metadata=load_cortical_policy(artifact)
    stream=NativeMechanicalStream(root,out/'plant',environment='free',target_mass_kg=70)
    initial=stream.snapshot();token=stream.checkpoint();records=[]
    try:
        for offset in (-.12,.12,.22):
            target=initial['coordinates']['ankle_angle_r']['value']+offset
            for arm in ('full','sever','no_controller'):
                stream.restore(token);adapter=PersistentCorticalCommands(policy,dt_s=.02);trajectory=[]
                for step in range(round(seconds/.02)):
                    commands=({name:0. for name in initial['muscles']} if arm=='no_controller' else adapter.commands(stream.state,target,sever=arm=='sever'))
                    state=stream.advance(.02,actuation=commands);q=state['coordinates']['ankle_angle_r']
                    trajectory.append({'time_s':state['time_s'],'angle_rad':q['value'],'speed_rad_s':q['speed'],
                        'tibant_r':commands['tibant_r'],'soleus_r':commands['soleus_r']})
                record={'arm':arm,'offset_rad':offset,'target_rad':target,'tracking_mse_rad2':float(np.mean([(r['angle_rad']-target)**2 for r in trajectory])),
                    'final_error_rad':trajectory[-1]['angle_rad']-target,'trajectory':trajectory}
                records.append(record);(out/'partial_results.json').write_text(json.dumps(records,indent=2)+'\n')
                print(json.dumps({k:v for k,v in record.items() if k!='trajectory'}),flush=True)
    finally:stream.release(token);stream.close()
    means={arm:float(np.mean([r['tracking_mse_rad2'] for r in records if r['arm']==arm])) for arm in ('full','sever','no_controller')}
    comparisons=[{'offset_rad':offset,'full_mse':next(r['tracking_mse_rad2'] for r in records if r['offset_rad']==offset and r['arm']=='full'),
        'sever_mse':next(r['tracking_mse_rad2'] for r in records if r['offset_rad']==offset and r['arm']=='sever')} for offset in (-.12,.12,.22)]
    receipt={'schema':'ihm.native-cortical-motor-eval.v1','artifact_path':str(Path(artifact).resolve()),'artifact_sha256':metadata['artifact_sha256'],
        'source_identity':metadata['source_identity'],'sites':policy.dyn.n,'metrics':means,'paired_comparisons':comparisons,'records':records,
        'kernel_contribution_demonstrated':all(r['full_mse']<r['sever_mse'] for r in comparisons),
        'source_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in sources.items()},
        'limits':['Actual IBM E/I equations with persistent state and disjoint sensory/motor ports; copied embedding and motor decoder trained',
            'Synthetic engineering PD imitation with privileged native joint feedback; not anatomical afferent/cord loop',
            'Held-out native ankle targets in free falling articulated body, not standing or walking',
            'Cortical state reinitialized only between matched arms; no resets within rollouts',
            'Severed association keeps decoder and E/I dynamics; no-controller zeros all muscle commands']}
    retained=out/'sources';retained.mkdir()
    for path,raw in sources.items():(retained/Path(path).name).write_bytes(raw)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(means),flush=True)
    return receipt
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);p.add_argument('--artifact',default='data/runtime/motor-learning/cortical-20260908/cortical_motor.pt');p.add_argument('--seconds',type=float,default=.3);a=p.parse_args()
    if not .02<=a.seconds<=2 or abs(a.seconds/.02-round(a.seconds/.02))>1e-8:p.error('Integral20ms steps, horizon .02..2s required')
    evaluate(Path(__file__).resolve().parents[1],a.output,a.artifact,a.seconds)
