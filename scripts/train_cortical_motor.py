"""Train copied actual IBM E/I materialization on synthetic ankle PD examples."""
import argparse
import json
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import torch
from ihm.assembly.ibm_controller import IBMImplicitController
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog
from ihm.native.cortical_motor import CorticalAnklePolicy


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--epochs',type=int,default=300)
    ap.add_argument('--out',default='data/runtime/motor-learning/cortical-20260908')
    args=ap.parse_args();torch.set_num_threads(1);torch.manual_seed(42)
    adapter=IBMImplicitController.from_root(ROOT,muscle_catalog=native_muscle_catalog(ROOT),sites=128)
    policy=CorticalAnklePolicy(adapter.loop.dyn,adapter.loop.sense_idx,adapter.loop.motor_idx)
    initial={k:v.detach().clone() for k,v in policy.state_dict().items()}
    opt=torch.optim.Adam([policy.dyn.embed,*policy.decoder.parameters()],lr=.005)
    start=time.monotonic()
    for epoch in range(args.epochs):
        x=torch.rand(32,2)*2-1;x[:,0]*=.18;x[:,1]*=.6
        target=(3*x[:,0]-.22*x[:,1]).clamp(-.45,.45)
        y=torch.stack((target.clamp_min(0),(-target).clamp_min(0)),-1)
        out,state=policy.advance(x,policy.state(len(x)),ticks=100)
        loss=(out-y).square().mean()
        opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(policy.parameters(),1.);opt.step()
        if epoch%25==0:print(json.dumps({'epoch':epoch,'mse':float(loss.detach()),'elapsed_s':time.monotonic()-start}),flush=True)
    torch.manual_seed(981);x=torch.rand(256,2)*2-1;x[:,0]*=.18;x[:,1]*=.6
    target=(3*x[:,0]-.22*x[:,1]).clamp(-.45,.45)
    y=torch.stack((target.clamp_min(0),(-target).clamp_min(0)),-1)
    with torch.no_grad():
        full,_=policy.advance(x,policy.state(len(x)),ticks=100)
        sever,_=policy.advance(x,policy.state(len(x)),ticks=100,sever=True)
    metrics={'full_mse':float((full-y).square().mean()),'sever_mse':float((sever-y).square().mean()),
        'full_signed_corr':float(torch.corrcoef(torch.stack((full[:,0]-full[:,1],target)))[0,1]),
        'sever_command_sd':float(sever.std(0).max()),'wall_s':time.monotonic()-start,
        'epochs':args.epochs,'teacher':'clip(3*position_error-.22*angular_velocity,-.45,.45)',
        'evidence':'Heldout synthetic observations only; native tracking not yet measured',
        'ports':'fixed signed sensory stimulation; trained bias-free precentral readout',
        'dynamic_model':'actual IBM CorticalDynamics128sites;100ms initial-value trials,1ms integration'}
    out=ROOT/args.out;out.mkdir(parents=True,exist_ok=True)
    (out/'pretrain_video_loop.py').write_bytes(adapter.source_bytes['pretrain_video_loop.py'])
    (out/'cortical_motor.py').write_bytes((ROOT/'ihm/native/cortical_motor.py').read_bytes())
    torch.save({'schema':'ihm.ibm-cortical-ankle.v1','state_dict':policy.state_dict(),
        'initial_state_dict':initial,'source_identity':adapter.identity,'metrics':metrics},out/'cortical_motor.pt')
    (out/'training.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics,indent=2),flush=True)

if __name__=='__main__':main()
