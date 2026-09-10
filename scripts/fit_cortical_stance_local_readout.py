#!/usr/bin/env python3
"""Fit expanded motor decoder to actual frozen E/I local steady-state responses."""
import argparse,copy,hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from scipy import sparse
from scipy.sparse.linalg import splu
from ihm.native.cortical_stance import CorticalStancePolicy,load_cortical_stance
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--artifact',default='data/runtime/motor-learning/cortical-stance-patient-fullstate-20260908/cortical_stance.pt')
 ap.add_argument('--output',default='data/runtime/motor-learning/cortical-stance-patient-fullstate-localreadout-20260908')
 args=ap.parse_args();torch.set_num_threads(1);started=time.monotonic()
 source=ROOT/args.artifact;out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
 policy,artifact=load_cortical_stance(source,dt_s=.01)
 teacher_path=ROOT/'data/research/locomotion_control/delay_design_e5s_lepu/candidate_state1.0_R1000.0.npz'
 teacher=NativeStanceLQR(teacher_path,model_sha256=artifact['provenance']['model_sha256'],dt_s=.01,target_mass_kg=MECHANICAL_TARGET_MASS_KG)
 assert policy.encoder_kind=='signed_identity' and policy.reference_normalization
 assert policy.state_names==teacher.state_names and policy.muscle_names==teacher.muscle_names
 initial={k:v.detach().clone() for k,v in policy.state_dict().items()}
 with torch.no_grad():weights=policy.dyn.edge_weights().detach().double()
 policy.double();dyn=policy.dyn;n=dyn.n;nx=len(policy.x0);h=.001
 state=tuple(v.double() for v in dyn.init_state(1,'cpu'));zero=torch.zeros(1,n,dtype=torch.float64)
 with torch.no_grad():
  for _ in range(10000):state=dyn.step(state,zero,h,weights)
  nextstate=dyn.step(state,zero,h,weights)
 rest_delta=max(float((a-b).abs().max()) for a,b in zip(nextstate,state))
 v,r,a,gi=[x[0].numpy() for x in state];I=sparse.eye(n,format='csr');Z=sparse.csr_matrix((n,n))
 scalar=lambda x:float(x.detach());tm=dyn.tau_m;ta=dyn.tau_a
 W=sparse.coo_matrix((weights.numpy().ravel(),(np.repeat(np.arange(n),dyn.idx.shape[1]),dyn.idx.numpy().ravel())),shape=(n,n)).tocsr()
 rate=dyn.r_max/(1+np.exp(-(v-dyn.v_half)/dyn.slope));rp=rate*(1-rate/dyn.r_max)/dyn.slope
 J=sparse.bmat([[sparse.diags(1-h/tm*(1+np.maximum(gi,0))),h/tm*20/dyn.r_max*(scalar(dyn.w_ee)*I+scalar(dyn.w_assoc)*W),-h/tm*I,sparse.diags(-h/tm*(v-dyn.e_rev)*(gi>0))],
  [sparse.diags(h/.005*rp),(1-h/.005)*I,Z,Z],[sparse.diags(h/ta*scalar(dyn.a_gain)*rp),Z,(1-h/ta)*I,Z],
  [Z,h/.008*scalar(dyn.w_ei)/dyn.r_max*I,Z,(1-h/.008)*I]],format='csr')
 G=np.zeros((4*n,nx));G[policy.sensory_sites.numpy()]=h/tm*12*policy.encoder.detach().numpy()
 system=(sparse.eye(4*n)-J).tocsc();H=splu(system).solve(G)
 solve_error=float(np.linalg.norm(system@H-G)/np.linalg.norm(G))
 sensory=set(policy.sensory_sites.tolist());motors=torch.tensor([i for i in range(n) if i not in sensory],dtype=torch.int64)
 response=H[n+motors.numpy()];response=response-response.mean(axis=0,keepdims=True)
 features=(2*response).T;target=(-teacher.K*policy.state_scale.numpy()[None,:]).T
 U,s,Vt=np.linalg.svd(features,full_matrices=False);tolerance=s[0]*max(features.shape)*np.finfo(float).eps
 rank=int(np.sum(s>tolerance));projection=U[:,:rank]@U[:,:rank].T@target
 null_residual=float(np.linalg.norm(target-projection)/np.linalg.norm(target))
 sweep=[];candidates=[]
 for relative in (1e-12,1e-10,1e-8,1e-6,1e-4):
  ridge=relative*s[0]**2
  decoder=((Vt.T*(s/(s*s+ridge)))@(U.T@target)).T.astype(np.float32)
  decoder-=decoder.mean(axis=1,keepdims=True)
  fit=features@decoder.T
  sweep.append({'relative_ridge':relative,'ridge':float(ridge),'gain_relative_error':float(np.linalg.norm(fit-target)/np.linalg.norm(target)),
   'gain_mse':float(np.mean((fit-target)**2)),'decoder_norm':float(np.linalg.norm(decoder))})
  candidates.append(decoder)
 selected=int(np.argmin([r['gain_relative_error'] for r in sweep]));decoder=candidates[selected]
 # Architecture change is explicit: all nonsensory sites become readout sites.
 policy.float();policy.motor_sites=motors;policy.decoder=torch.nn.Linear(len(motors),len(teacher.muscle_names),bias=False)
 with torch.no_grad():policy.decoder.weight.copy_(torch.from_numpy(decoder))
 assert all(torch.equal(value,policy.state_dict()[key]) for key,value in initial.items() if key not in ('motor_sites','decoder.weight'))
 report={'schema':'ihm.cortical-local-readout-fit.v1','change':'Decoder-only steady-state gain refit with declared expansion to every nonsensory motor site',
  'sites':n,'sensory_count':len(sensory),'old_motor_count':len(initial['motor_sites']),'motor_count':len(motors),
  'state_width':nx,'H_motor_rank':rank,'H_singular_max':float(s[0]),'H_singular_min':float(s[-1]),
  'H_condition_number':float(s[0]/s[-1]),'teacher_nullspace_relative_residual':null_residual,
  'steady_solve_relative_residual':solve_error,'rest_tick_max_abs_delta':rest_delta,'sweep':sweep,'selected':sweep[selected],
  'parent_artifact_sha256':sha(source),'teacher_sha256':sha(teacher_path),'script_sha256':sha(__file__),
  'nondecoder_core_tensors_unchanged':True,'scope':['K supplies offline local target only; deployed output is actual E/I motor rates through bias-free decoder',
   'Steady-state gain matching does not remove intrinsic neural latency or establish closed-loop stability',
   'Normalized native directions are linearized at nominal; clipping and tanh limit finite excursions',
   'All existing embedding, equation parameters and signed sensory encoding retained exactly; motor indices explicitly expanded']}
 result=copy.deepcopy(artifact);result.pop('artifact_sha256',None);result['parent_report']=result.get('report');result['report']=report
 result['state_dict']={k:v.detach().clone() for k,v in policy.state_dict().items()}
 result['provenance']['local_decoder_refit']={'parent_artifact_sha256':sha(source),'teacher_sha256':sha(teacher_path),
  'method':'(I-J)^-1G on actual frozen 1ms E/I Jacobian; normalized native fullstate signed input; mean-free motor-rate decoder',
  'motor_sites':motors.tolist(),'selected_relative_ridge':sweep[selected]['relative_ridge'],'script_sha256':sha(__file__)}
 for p in (source.parent/'pretrain_video_loop.py',ROOT/'ihm/native/cortical_stance.py',Path(__file__)):(out/p.name).write_bytes(p.read_bytes())
 torch.save(result,out/'cortical_stance.pt')
 loaded,_=load_cortical_stance(out/'cortical_stance.pt',dt_s=.01,model_sha256=teacher.model_sha256)
 # Float32 deployed replay on finite held-constant normalized perturbations. This
 # is a core replay only and intentionally has no native mechanical feedback.
 rng=np.random.default_rng(47);directions=rng.standard_normal((4,nx));directions/=np.linalg.norm(directions,axis=1,keepdims=True)
 amplitude=.01;X=loaded.x0[None]+torch.tensor(amplitude*directions,dtype=torch.float32)*loaded.state_scale[None]
 runtime=loaded.state(4);sever=loaded.state(4)
 with torch.no_grad():
  actualweights=loaded.dyn.edge_weights().detach()
  for _ in range(300):
   prediction,runtime=loaded.advance(X,runtime,ticks=10,weights=actualweights)
   severout,sever=loaded.advance(X,sever,ticks=10,weights=actualweights,sever=True)
   assert torch.equal(severout,loaded.u0[None].expand_as(severout))
 actual=(prediction-loaded.u0).numpy();wanted=amplitude*directions@target
 report['float32_constant_input_replay']={'duration_s':3.,'normalized_amplitude':amplitude,'directions':4,
  'correction_relative_error':float(np.linalg.norm(actual-wanted)/np.linalg.norm(wanted)),
  'correction_mse':float(np.mean((actual-wanted)**2)),'sever_exact_baseline':True,
  'target_correction_norm':float(np.linalg.norm(wanted)),'actual_correction_norm':float(np.linalg.norm(actual))}
 report['artifact_sha256']=sha(out/'cortical_stance.pt');report['wall_seconds']=time.monotonic()-started
 (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 np.savez_compressed(out/'local_fit_arrays.npz',singular_values=s,features=features,target=target,decoder=decoder,
  normalized_test_directions=directions,replay_correction=actual,replay_target=wanted)
 print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
