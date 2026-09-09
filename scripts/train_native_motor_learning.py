#!/usr/bin/env python3
"""Train copied IBM kernel sites against ankle servo and evaluate native ablations."""
import argparse,json,sys,hashlib,time
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.motor_learning import KernelAnklePolicy,source_embedding
from ihm.native.mechanical_stream import NativeMechanicalStream

def run(root,out,checkpoint,epochs=600,seconds=.4,seed=17):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1);torch.manual_seed(seed)
    embedding,provenance=source_embedding(checkpoint)
    policy=KernelAnklePolicy(embedding);untrained=KernelAnklePolicy(embedding)
    # Offline synthetic error/velocity coverage, not native demonstration data.
    g=torch.Generator().manual_seed(seed)
    x=torch.stack((torch.rand(4096,generator=g)*.8-.4,torch.rand(4096,generator=g)*4-2),-1)
    differential=3*x[:,0]-.15*x[:,1]
    labels=torch.stack((differential.clamp(0,1),(-differential).clamp(0,1)),-1)
    opt=torch.optim.Adam(policy.parameters(),lr=.003);losses=[]
    for epoch in range(epochs):
        loss=(policy(x)-labels).square().mean();opt.zero_grad();loss.backward();opt.step()
        if epoch%50==0 or epoch==epochs-1:losses.append({'epoch':epoch,'mse':float(loss.detach())})
    artifact={'schema':'ihm.ibm-native-ankle-kernel.v1','embed':policy.embed.detach(),
              'provenance':provenance,'teacher':{'position_gain':3.,'velocity_gain':.15},
              'materialization':'8 uniform IBM embedding sites; normalized tanh kernel; fixed signed ports; no E/I integration',
              'training':'synthetic error/velocity PD imitation; embedding only','seed':seed}
    torch.save(artifact,out/'motor_kernel.pt')
    stream=NativeMechanicalStream(root,out/'plant',environment='free',target_mass_kg=70)
    initial=stream.snapshot();token=stream.checkpoint();records=[]
    try:
        q0=initial['coordinates']['ankle_angle_r']['value']
        # Held out target values. Each arm restores identical complete native state.
        for offset in (-.12,.12,.22):
            target=q0+offset
            for arm in ('full','sever','untrained','no_controller','teacher'):
                stream.restore(token);trajectory=[]
                for step in range(round(seconds/.01)):
                    q=stream.state['coordinates']['ankle_angle_r']
                    if arm=='no_controller':cmd={n:0. for n in initial['muscles']}
                    elif arm=='teacher':
                        d=3*(target-q['value'])-.15*q['speed'];cmd={n:0. for n in initial['muscles']};cmd['tibant_r']=float(np.clip(d,0,1));cmd['soleus_r']=float(np.clip(-d,0,1))
                    else:cmd=(untrained if arm=='untrained' else policy).commands(stream.state,target,sever=arm=='sever')
                    state=stream.advance(.01,actuation=cmd);q=state['coordinates']['ankle_angle_r']
                    trajectory.append({'time_s':state['time_s'],'angle_rad':q['value'],'speed_rad_s':q['speed'],'target_rad':target,'tibant_r':cmd['tibant_r'],'soleus_r':cmd['soleus_r']})
                mse=float(np.mean([(v['angle_rad']-target)**2 for v in trajectory]))
                record={'arm':arm,'offset_rad':offset,'target_rad':target,'tracking_mse_rad2':mse,'final_error_rad':trajectory[-1]['angle_rad']-target,'trajectory':trajectory}
                records.append(record);print(json.dumps({k:v for k,v in record.items() if k!='trajectory'}),flush=True)
    finally:stream.release(token);stream.close()
    means={arm:float(np.mean([r['tracking_mse_rad2'] for r in records if r['arm']==arm])) for arm in ('full','sever','untrained','no_controller','teacher')}
    receipt={'schema':'ihm.ibm-native-motor-learning-eval.v1','training_losses':losses,'provenance':provenance,'metrics':means,'records':records,
       'kernel_contribution_demonstrated':means['full']<means['sever'] and means['full']<means['no_controller'],
       'training_improves_untrained':means['full']<means['untrained'],
       'limits':['Isolated ankle target tracking in free falling native articulated body, no contact support','Synthetic PD imitation, no demonstrations or reward learning','Reduced 8-site kernel materialization, not IBM cortical E/I dynamics','No walking or general human brain fidelity claim','Original IBM checkpoint not changed'],
       'artifact_sha256':hashlib.sha256((out/'motor_kernel.pt').read_bytes()).hexdigest(),
       'source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'ihm/native/motor_learning.py',root/'scripts/train_native_motor_learning.py')}}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);p.add_argument('--checkpoint',default=str(Path.home()/'Documents/IBM-1/ckpt/ibm1_implicit.pt'));p.add_argument('--epochs',type=int,default=600);p.add_argument('--seconds',type=float,default=.4);a=p.parse_args()
    if not 1<=a.epochs<=10000 or not .01<=a.seconds<=2:p.error('Invalid bounded epochs or horizon')
    run(Path(__file__).resolve().parents[1],a.output,a.checkpoint,a.epochs,a.seconds)
