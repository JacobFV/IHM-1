#!/usr/bin/env python3
"""Decoder-only ridge refit using persistent paired-reference IBM E/I features."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from ihm.native.cortical_stance import load_cortical_stance
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--artifact',default='data/runtime/motor-learning/cortical-stance-patient-slow-20260908/cortical_stance.pt',help='Source cortical artifact; dimensions and encoder kind are loaded from its provenance')
 p.add_argument('--output',default='data/runtime/motor-learning/cortical-stance-patient-readout-20260908')
 args=p.parse_args();torch.set_num_threads(1)
 source=(ROOT/args.artifact).resolve()
 teacher_path=ROOT/'data/research/locomotion_control/delay_design_e5s_lepu/candidate_state1.0_R1000.0.npz'
 paths=[ROOT/'data/runtime/motor-learning/patient-stance-axis-data-20260908'/f'{name}.json' for name in ('positive_x','negative_x','positive_z','negative_z')]
 paths.append(ROOT/'data/research/locomotion_control/patient_slow_lqr_delay4_com_push/trajectory.json')
 out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
 policy,artifact=load_cortical_stance(source,dt_s=.01)
 teacher=NativeStanceLQR(teacher_path,model_sha256=artifact['provenance']['model_sha256'],dt_s=.01,target_mass_kg=MECHANICAL_TARGET_MASS_KG)
 if not policy.reference_normalization:raise ValueError('Paired-reference cortical normalization required')
 assert policy.muscle_names==teacher.muscle_names and policy.state_names==teacher.state_names
 architecture={'sites':policy.dyn.n,'encoder_kind':policy.encoder_kind,
  'sensory_sites':policy.sensory_sites.tolist(),'motor_sites':policy.motor_sites.tolist(),
  'encoder_shape':list(policy.encoder.shape),'decoder_shape':list(policy.decoder.weight.shape),
  'reference_normalization':policy.reference_normalization}
 initial={k:v.clone() for k,v in policy.state_dict().items()}
 features=[];targets=[];predictions=[];severed=[];traces=[];started=time.monotonic()
 with torch.no_grad():
  weights=policy.dyn.edge_weights().detach()
  for path in paths:
   frames=[f for f in json.loads(path.read_text()) if f['time_s']<=3.+1e-9]
   state=policy.state();sever_state=policy.state();begin=len(features)
   for frame in frames:
    vector=policy.snapshot_vector(frame)
    prediction,state=policy.advance(vector,state,ticks=10,weights=weights)
    sever,sever_state=policy.advance(vector,sever_state,ticks=10,weights=weights,sever=True)
    motor=state[1][:,policy.motor_sites]-state[5][:,policy.motor_sites]
    anchored=10*(motor-motor[:,:1])
    feature=anchored-anchored.mean(-1,keepdim=True)
    features.append(feature[0].numpy().astype(np.float64))
    target=np.clip(teacher.u0-teacher.K@(teacher.state_vector(frame)-teacher.x0),.01,1.)
    targets.append(target);predictions.append(prediction[0].numpy());severed.append(sever[0].numpy())
   traces.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),'start_index':begin,'end_index':len(features),'frames':len(frames)})
   print(json.dumps({'features_collected':len(features),'trace':str(path)}),flush=True)
 F=np.asarray(features);Y=np.asarray(targets);baseline=policy.u0.numpy().astype(np.float64)
 old=np.asarray(predictions);S=np.asarray(severed)
 assert np.array_equal(S,np.broadcast_to(policy.u0.numpy(),S.shape)), 'Sever must equal baseline exactly'
 ratio=(Y-baseline)/.2
 latent=np.arctanh(np.clip(ratio,-1+1e-6,1-1e-6))
 held=np.arange(len(F))%5==4;train=~held
 U,s,Vt=np.linalg.svd(F[train],full_matrices=False)
 rank=int(np.sum(s>s[0]*max(F[train].shape)*np.finfo(float).eps))
 sweeps=[];candidates=[]
 for relative in (1e-10,1e-8,1e-6,1e-4,1e-3,1e-2):
  ridge=float(relative*s[0]**2)
  W=(Vt.T*(s/(s*s+ridge)))@(U.T@latent[train])
  # Match actual float32 decoder and exact anchor computation for comparison.
  W=W.T.astype(np.float32)
  W=W-W.mean(axis=1,keepdims=True)
  prediction=np.clip(baseline+.2*np.tanh(F@W.T),.01,1.)
  item={'relative_ridge':relative,'ridge':ridge,'train_mse':float(np.mean((prediction[train]-Y[train])**2)),
   'validation_mse':float(np.mean((prediction[held]-Y[held])**2)),'all_mse':float(np.mean((prediction-Y)**2)),
   'decoder_norm':float(np.linalg.norm(W))}
  sweeps.append(item);candidates.append(W)
 choice=int(np.argmin([r['validation_mse'] for r in sweeps]))
 with torch.no_grad():policy.decoder.weight.copy_(torch.from_numpy(candidates[choice]))
 actual=[]
 with torch.no_grad():
  for trace,path in zip(traces,paths):
   state=policy.state()
   for frame in [f for f in json.loads(path.read_text()) if f['time_s']<=3.+1e-9]:
    output,state=policy.advance(policy.snapshot_vector(frame),state,ticks=10,weights=weights)
    actual.append(output[0].numpy())
 actual=np.array(actual)
 assert all(torch.equal(value,policy.state_dict()[key]) for key,value in initial.items() if key!='decoder.weight')
 report={'schema':'ihm.cortical-stance-readout-refit.v1','change':'decoder-only refit; embedding and all other tensors exactly unchanged',
  'samples':len(F),'feature_width':F.shape[1],'sites':policy.dyn.n,'paired_reference':True,'architecture':architecture,
  'rank':rank,'singular_max':float(s[0]),'singular_min':float(s[-1]),
  'condition_number':float(s[0]/s[-1]),'effective_condition_number':float(s[0]/s[rank-1]),
  'sweep':sweeps,'selected':sweeps[choice],'baseline_mse':float(np.mean((baseline-Y)**2)),
  'original_persistent_mse':float(np.mean((old-Y)**2)),'refit_actual_persistent_mse':float(np.mean((actual-Y)**2)),
  'sever_exact_baseline':True,'latent_target_clipped_values':int(np.count_nonzero(abs(ratio)>=1-1e-6)),
  'traces':traces,'parent_artifact_sha256':sha(source),'teacher_sha256':sha(teacher_path),
  'limits':['Only training-distribution replay; no native test or stability claim',
   'Interleaved within-trace validation selects ridge; not independent generalization',
   'Four axis traces are pre-step snapshots; positive slowteacher trace is retained post-step snapshots',
   'No negative diagonal held-out trajectory accessed; no teacher bypass in runtime'],
  'wall_seconds':time.monotonic()-started}
 result=copy.deepcopy(artifact);result.pop('artifact_sha256',None)
 result['state_dict']={k:v.detach().clone() for k,v in policy.state_dict().items()}
 result['parent_report']=result.get('report');result['report']=report
 result['provenance']['decoder_refit']={'kind':'decoder-only-refit','parent_artifact_sha256':sha(source),
  'parent_artifact_path':str(source),'architecture':architecture,
  'teacher_sha256':sha(teacher_path),'traces':traces,'script_sha256':sha(__file__),
  'feature_basis':'10 * mean-free (motor response minus paired zero-input reference), evaluated through exact 10ms persistent cortical replay',
  'selected_relative_ridge':sweeps[choice]['relative_ridge']}
 torch.save(result,out/'cortical_stance.pt')
 for path in (source.parent/'pretrain_video_loop.py',ROOT/'ihm/native/cortical_stance.py',Path(__file__)):
  (out/path.name).write_bytes(path.read_bytes())
 report['artifact_sha256']=sha(out/'cortical_stance.pt')
 (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 np.savez_compressed(out/'fit_arrays.npz',features=F,targets=Y,original=old,refit=actual,validation_mask=held)
 loaded,_=load_cortical_stance(out/'cortical_stance.pt',model_sha256=teacher.model_sha256,dt_s=.01)
 assert torch.equal(loaded.dyn.embed,initial['dyn.embed'])
 print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
