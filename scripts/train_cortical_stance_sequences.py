"""Patient-bound temporal E/I imitation from four native perturbation directions.

Negative diagonal is held out. Loss is on actual successive10ms teacher commands
through persistent cortical history, not static endpoints or direct PD bypass.
"""
import argparse,hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.native.cortical_stance import CorticalStancePolicy,load_cortical_stance
from ihm.assembly.ibm_controller import IBMImplicitController
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data',default='data/runtime/motor-learning/patient-stance-axis-data-20260908')
    ap.add_argument('--out',default='data/runtime/motor-learning/cortical-stance-patient-sequences-20260908')
    ap.add_argument('--epochs',type=int,default=160);ap.add_argument('--sites',type=int,default=512)
    ap.add_argument('--encoder-kind',choices=['random_projection','signed_identity'],default='random_projection')
    ap.add_argument('--teacher');ap.add_argument('--initial');ap.add_argument('--extra-trajectory',action='append',default=[])
    # Which IBM association kernel initialises the cortex. Default keeps the original
    # fused implicit checkpoint this artifact family was trained from.
    ap.add_argument('--checkpoint',help='Retained IBM kernel to initialise dyn.embed; default is the donor fused implicit checkpoint')
    ap.add_argument('--lr',type=float,default=.0015);ap.add_argument('--reserve-weight',action='store_true')
    args=ap.parse_args();torch.set_num_threads(1);torch.manual_seed(71)
    source_raw=Path(__file__).read_bytes();policy_raw=(ROOT/'ihm/native/cortical_stance.py').read_bytes()
    bundle=ROOT/'data/models/engineering_stance_v1';manifest=json.loads((bundle/'manifest.json').read_text())
    teacher_path=Path(args.teacher) if args.teacher else bundle/'linearization.npz';teacher_raw=teacher_path.read_bytes()
    teacher=NativeStanceLQR(teacher_path,model_sha256=manifest['native_model_sha256'],dt_s=.01,target_mass_kg=manifest['target_mass_kg'])
    data=Path(args.data);receipt=json.loads((data/'report.json').read_text())
    if receipt['model_sha256']!=teacher.model_sha256 or receipt['artifact_sha256']!=hashlib.sha256((bundle/'linearization.npz').read_bytes()).hexdigest():
        raise ValueError('Training state owner/teacher differs from promoted patient')
    if receipt['target_mass_kg']!=manifest['target_mass_kg'] or not receipt['all_four_complete']:raise ValueError('Incomplete/wrong patient dataset')
    collector_teacher=NativeStanceLQR(bundle/'linearization.npz',model_sha256=teacher.model_sha256,dt_s=.01,target_mass_kg=manifest['target_mass_kg'])
    if teacher.state_names!=collector_teacher.state_names or teacher.muscle_names!=collector_teacher.muscle_names or not np.array_equal(teacher.x0,collector_teacher.x0) or not np.array_equal(teacher.u0,collector_teacher.u0):raise ValueError('Relabeling changes native baseline/state ownership')
    xs=[];ys=[];data_hashes={}
    for arm in ('positive_x','negative_x','positive_z','negative_z'):
        path=data/(arm+'.json');raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
        if digest!=receipt['arms'][arm]['frames_sha256']:raise ValueError('Changed native trajectory')
        frames=json.loads(raw);vectors=np.array([teacher.state_vector(frame) for frame in frames])
        if not all(abs(frame['time_s']-i*.01)<1e-8 for i,frame in enumerate(frames)):raise ValueError('Missing native sequence sample')
        raw_commands=teacher.u0-(vectors-teacher.x0)@teacher.K.T
        xs.append(torch.tensor(vectors,dtype=torch.float32));ys.append(torch.tensor(np.clip(raw_commands,.01,1),dtype=torch.float32))
        data_hashes[arm]=digest
    for number,extra in enumerate(args.extra_trajectory):
        path=Path(extra);raw=path.read_bytes();extra_report=json.loads(path.with_name('report.json').read_text())
        if extra_report.get('target_mass_kg')!=manifest['target_mass_kg'] or hashlib.sha256(teacher_raw).hexdigest() not in extra_report.get('source_sha256',{}).values():raise ValueError('Extra native trajectory teacher/mass differs')
        frames=[f for f in json.loads(raw) if f['time_s']<=3.+1e-9]
        vectors=np.array([teacher.state_vector(frame) for frame in frames])
        raw_commands=teacher.u0-(vectors-teacher.x0)@teacher.K.T
        xs.append(torch.tensor(vectors,dtype=torch.float32));ys.append(torch.tensor(np.clip(raw_commands,.01,1),dtype=torch.float32))
        data_hashes['extra_'+str(number)]=hashlib.sha256(raw).hexdigest()
    scales=[.05 if name.endswith('/speed') else .01 if '/jointset/' in name else .005 if name.endswith('/activation') else .0001 for name in teacher.state_names]
    core=IBMImplicitController.from_root(ROOT,muscle_catalog=native_muscle_catalog(ROOT),sites=args.sites,
        checkpoint_path=Path(args.checkpoint).resolve() if args.checkpoint else None)
    sensory,motor=core.loop.sense_idx,core.loop.motor_idx
    if args.encoder_kind=='signed_identity':
        width=len(teacher.state_names);motors=2*len(teacher.muscle_names)
        if args.sites<2*width+motors:raise ValueError('Insufficient sites for full signed sensory and motor ports')
        sensory=torch.arange(2*width);motor=torch.arange(2*width,2*width+motors)
    policy=CorticalStancePolicy(core.loop.dyn,sensory,motor,
        state_names=teacher.state_names,muscle_names=teacher.muscle_names,x0=teacher.x0,
        u0=np.maximum(teacher.u0,.01),state_scale=scales,reference_normalization=True,encoder_kind=args.encoder_kind)
    initial={k:v.detach().clone() for k,v in policy.state_dict().items()};lineage=None
    if args.initial:
        policy,previous=load_cortical_stance(args.initial,model_sha256=teacher.model_sha256,dt_s=.01)
        if policy.dyn.n!=args.sites or policy.muscle_names!=teacher.muscle_names or policy.state_names!=teacher.state_names or not torch.equal(policy.x0,torch.tensor(teacher.x0,dtype=torch.float32)) or not torch.equal(policy.u0,torch.tensor(teacher.u0,dtype=torch.float32)):raise ValueError('Fine-tune controller native binding differs')
        initial=previous['initial_state_dict'];lineage=previous['artifact_sha256']
    loss_scale=torch.stack(ys).flatten(0,1).std(0).clamp_min(.0002)
    if args.reserve_weight:loss_scale=torch.minimum(loss_scale,(policy.u0-.01).clamp_min(.00002))
    optimizer=torch.optim.Adam([policy.dyn.embed,*policy.decoder.parameters()],lr=args.lr)
    generator=np.random.default_rng(71);started=time.monotonic();window,burn,batch=30,10,8
    for epoch in range(args.epochs):
        arms=generator.integers(0,len(xs),size=batch)
        # Half the windows include force onset; the rest cover nominal/recovery.
        starts=np.r_[generator.integers(75,110,size=batch//2),generator.integers(0,270,size=batch//2)]
        inputs=torch.stack([xs[a][s:s+window] for a,s in zip(arms,starts)],1)
        targets=torch.stack([ys[a][s:s+window] for a,s in zip(arms,starts)],1)
        with torch.no_grad():
            _,single=policy.advance(policy.x0[None],policy.state(),ticks=500)
            state=tuple(value.expand(batch,-1).clone() for value in single)
        weights=policy.dyn.edge_weights();loss=0.
        for index in range(window):
            prediction,state=policy.advance(inputs[index],state,ticks=10,weights=weights)
            if index==burn-1:state=tuple(value.detach() for value in state)
            if index>=burn:loss=loss+(((prediction-targets[index])/loss_scale)**2).mean()/(window-burn)
        optimizer.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(policy.parameters(),1.);optimizer.step()
        if epoch%20==0:print(json.dumps({'epoch':epoch,'normalized_temporal_loss':float(loss.detach()),'wall_s':time.monotonic()-started}),flush=True)
    # Full3s temporal replay of training axes is a fit diagnostic only.
    metrics={}
    with torch.no_grad():
        weights=policy.dyn.edge_weights()
        for sever in (False,True):
            predictions=[]
            for inputs in xs:
                state=policy.state();values=[]
                for vector in inputs:
                    value,state=policy.advance(vector[None],state,ticks=10,weights=weights,sever=sever)
                    values.append(value[0])
                predictions.append(torch.stack(values))
            errors=(torch.stack(predictions)-torch.stack(ys)).square()
            metrics['sever' if sever else 'full']=float(errors.mean())
        _,rest=policy.advance(policy.x0[None],policy.state(),ticks=1000,weights=weights)
        nominal,_=policy.advance(policy.x0[None],rest,ticks=10,weights=weights)
        assert torch.equal(nominal[0],policy.u0)
    provenance={'teacher_sha256':hashlib.sha256(teacher_raw).hexdigest(),'collection_teacher_sha256':receipt['artifact_sha256'],'initial_artifact_sha256':lineage,'trajectory_sha256':data_hashes,
        'model_sha256':teacher.model_sha256,'target_mass_kg':manifest['target_mass_kg'],'dt_s':.01,
        'encoder_kind':policy.encoder_kind,
        'port_allocation':'Engineering full signed state ports, disjoint motor ports; no anatomical region claim' if policy.encoder_kind=='signed_identity' else 'Original materialized region strips',
        'native_activation_floor':.01,'brain_checkpoint_sha256':core.identity['checkpoint_sha256'],
        'brain_checkpoint_source':args.checkpoint or 'donor IBM-1 ckpt/ibm1_implicit.pt',
        'cortical_source_sha256':core.identity['implementation_sha256'],'reference_normalization':True,
        'reference_basis':'Persistent same-kernel zero-input counterfactual subtracts intrinsic neural activity; no native control bypass',
        'teacher':'Exact patient margin-aware native LQR; current-frame commands recomputed; physical0.01floor',
        'dataset':f'{len(xs)} complete actualpatient3s trajectories; four axis pushes plus declaredextra; negative diagonal held out entirely',
        'state_names':teacher.state_names,'muscle_names':teacher.muscle_names}
    report={'epochs':args.epochs,'sites':args.sites,'learning_rate':args.lr,'reserve_weight':args.reserve_weight,'state_width':len(teacher.state_names),'motor_width':len(teacher.muscle_names),
        'sensory_sites':len(policy.sensory_sites),'motor_sites':len(policy.motor_sites),
        'temporal_training_axis_mse':metrics,'zero_error_baseline_exact':True,'wall_s':time.monotonic()-started,
        'training_window_s':.3,'burn_in_s':.1,'scored_contiguous_s':.2,
        'native_stance_demonstrated':False,'scope':'Temporal imitationfit diagnostic; actual native and independentnegative-diagonal evaluation still required'}
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    torch.save({'schema':'ihm.ibm-cortical-stance.v1','state_dict':policy.state_dict(),
        'initial_state_dict':initial,'provenance':provenance,'report':report},out/'cortical_stance.pt')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    for name,raw in [('pretrain_video_loop.py',core.source_bytes['pretrain_video_loop.py']),('cortical_stance.py',policy_raw),('train_cortical_stance_sequences.py',source_raw),('teacher.npz',teacher_raw)]:
        (out/name).write_bytes(raw)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
