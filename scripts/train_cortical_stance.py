"""Bounded E/I stance imitation pilot on retained perturbed native snapshots."""
import argparse,hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.native.cortical_stance import CorticalStancePolicy
from ihm.assembly.ibm_controller import IBMImplicitController
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--teacher',default='data/research/locomotion_control/linearization_3q5uh3o3/discrete/linearization.npz')
    ap.add_argument('--trajectory',default='data/research/locomotion_control/lqr98_discrete_push/trajectory.json')
    ap.add_argument('--out',default='data/runtime/motor-learning/cortical-stance-pilot-20260908')
    ap.add_argument('--epochs',type=int,default=100)
    ap.add_argument('--max-time',type=float,default=2.)
    ap.add_argument('--nominal-weight',type=float,default=1.)
    ap.add_argument('--row-regularization',type=float,default=0.)
    ap.add_argument('--variable-horizon',action='store_true')
    ap.add_argument('--reference-normalization',action='store_true');args=ap.parse_args()
    torch.set_num_threads(1);torch.manual_seed(47)
    teacher=NativeStanceLQR(args.teacher,dt_s=.01)
    frames=json.loads(Path(args.trajectory).read_text())
    frames=[f for f in frames if f['time_s']<=args.max_time+1e-9]
    if len(frames)<100 or not any(f.get('external_forces') for f in frames):raise ValueError('Retained native perturbation frames required')
    vectors=np.array([teacher.state_vector(f) for f in frames]);raw=teacher.u0-(vectors-teacher.x0)@teacher.K.T
    targets=np.clip(raw,.01,1.)
    # Complete native state order retained; the frozen random cortical projection
    # does not contain K or a teacher command shortcut.
    scales=[]
    for name in teacher.state_names:
        scales.append(.05 if name.endswith('/speed') else .01 if '/jointset/' in name else .005 if name.endswith('/activation') else .0001)
    adapter=IBMImplicitController.from_root(ROOT,muscle_catalog=native_muscle_catalog(ROOT),sites=256)
    policy=CorticalStancePolicy(adapter.loop.dyn,adapter.loop.sense_idx,adapter.loop.motor_idx,
        state_names=teacher.state_names,muscle_names=teacher.muscle_names,x0=teacher.x0,
        u0=np.maximum(teacher.u0,.01),state_scale=scales,reference_normalization=args.reference_normalization)
    initial={k:v.detach().clone() for k,v in policy.state_dict().items()}
    indices=np.arange(len(vectors));test=indices[indices%5==4];train=indices[indices%5!=4]
    X=torch.tensor(vectors,dtype=torch.float32);Y=torch.tensor(targets,dtype=torch.float32)
    # Error weights prevent inactive muscles and tiny nominal perturbations from
    # making an unchanged baseline look like successful feedback imitation.
    loss_scale=torch.tensor(np.maximum(np.std(targets[train],axis=0),.002),dtype=torch.float32)
    optimizer=torch.optim.Adam([policy.dyn.embed,*policy.decoder.parameters()],lr=.003)
    start=time.monotonic()
    for epoch in range(args.epochs):
        choose=np.random.default_rng(epoch).choice(train,24)
        x=torch.cat((X[choose],policy.x0[None].expand(8,-1)))
        y=torch.cat((Y[choose],policy.u0[None].expand(8,-1)))
        ticks=[20,40,80,160][epoch%4] if args.variable_horizon else 80
        prediction,_=policy.advance(x,policy.state(len(x)),ticks=ticks)
        residual=((prediction-y)/loss_scale)**2
        loss=residual[:24].mean()+args.nominal_weight*residual[24:].mean()
        if args.row_regularization:loss=loss+args.row_regularization*policy.dyn.edge_weights().sum(-1).var()
        optimizer.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(policy.parameters(),1.);optimizer.step()
        if epoch%20==0:print(json.dumps({'epoch':epoch,'normalized_loss':float(loss.detach()),'wall_s':time.monotonic()-start}),flush=True)
    with torch.no_grad():
        full,_=policy.advance(X[test],policy.state(len(test)),ticks=80)
        sever,_=policy.advance(X[test],policy.state(len(test)),ticks=80,sever=True)
        equilibrium,_=policy.advance(policy.x0[None],policy.state(),ticks=80)
        saved=policy.dyn.embed.clone();policy.dyn.embed.copy_(initial['dyn.embed'])
        reset,_=policy.advance(X[test],policy.state(len(test)),ticks=80);policy.dyn.embed.copy_(saved)
        errors={name:float((value-Y[test]).square().mean()) for name,value in [('full',full),('sever',sever),('reset_embedding',reset)]}
        equilibrium_error=float((equilibrium-policy.u0).abs().max())
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    provenance={'teacher_sha256':hashlib.sha256(Path(args.teacher).read_bytes()).hexdigest(),
        'trajectory_sha256':hashlib.sha256(Path(args.trajectory).read_bytes()).hexdigest(),
        'model_sha256':teacher.model_sha256,'dt_s':teacher.dt_s,'native_activation_floor':.01,
        'brain_checkpoint_sha256':adapter.identity['checkpoint_sha256'],
        'cortical_source_sha256':adapter.identity['implementation_sha256'],
        'reference_normalization':args.reference_normalization,
        'reference_basis':'Persistent zero-input counterfactual E/I state subtracts intrinsic activity; same kernel; no native command bypass' if args.reference_normalization else 'none',
        'teacher':'sampled-data native LQR; effective commands clipped to physical0.01activation floor',
        'dataset':f'Actual perturbed native snapshots through{args.max_time}s; interleaved within-trajectory holdout; no independent trial generalization',
        'state_names':teacher.state_names,'muscle_names':teacher.muscle_names}
    report={'errors':errors,'nominal_max_command_error':equilibrium_error,'epochs':args.epochs,
        'nominal_weight':args.nominal_weight,'row_regularization':args.row_regularization,'variable_horizon':args.variable_horizon,
        'training_samples':len(train),'heldout_samples':len(test),'sites':policy.dyn.n,'wall_s':time.monotonic()-start,
        'core_imitation_improved':errors['full']<errors['sever'],
        'native_stance_demonstrated':False,'provenance':provenance,
        'limitations':['80ms initial-value imitation trials; deployed cortex must preserve history',
            'No nonlinear native full/sever stance trial of this policy yet',
            'Teacher itself not robust to sustained pushes; nominal error matters near activation floor',
            'Intrinsic cortical propagation delay is separate from imitation MSE']}
    torch.save({'schema':'ihm.ibm-cortical-stance.v1','state_dict':policy.state_dict(),
        'initial_state_dict':initial,'provenance':provenance,'report':report},out/'cortical_stance.pt')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/'pretrain_video_loop.py').write_bytes(adapter.source_bytes['pretrain_video_loop.py'])
    (out/'cortical_stance.py').write_bytes((ROOT/'ihm/native/cortical_stance.py').read_bytes())
    (out/'train_cortical_stance.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({k:v for k,v in report.items() if k!='provenance'},indent=2),flush=True)

if __name__=='__main__':main()
