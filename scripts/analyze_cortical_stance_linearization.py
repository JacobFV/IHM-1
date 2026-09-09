#!/usr/bin/env python3
"""Local frozen-equation cortical/native augmented spectral diagnostic; no native run."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, eigs, ArpackNoConvergence
from ihm.native.cortical_stance import load_cortical_stance


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--artifact',default='data/runtime/motor-learning/cortical-stance-patient-slow-20260908/cortical_stance.pt')
 parser.add_argument('--plant',default='data/research/locomotion_control/delay_design_e5s_lepu/candidate_state1.0_R1000.0.npz')
 parser.add_argument('--output',default='data/runtime/motor-learning/cortical-stance-patient-slow-linearization-20260908')
 parser.add_argument('--rest-seconds',type=float,default=10.)
 parser.add_argument('--maxiter',type=int,default=400)
 args=parser.parse_args();torch.set_num_threads(1)
 if not 1<=args.rest_seconds<=30 or not 10<=args.maxiter<=2000:raise ValueError('Bounded rest/maxiter required')
 artifact_path=ROOT/args.artifact;plant_path=ROOT/args.plant
 output=ROOT/args.output;output.mkdir(parents=True,exist_ok=False)
 start=time.monotonic();policy,artifact=load_cortical_stance(artifact_path,dt_s=.01)
 if not policy.reference_normalization:raise ValueError('Paired-reference dynamics required')
 with np.load(plant_path,allow_pickle=False) as data:
  Ad=data['Ad'].copy();Bd=data['Bd'].copy();x0=data['x0'].copy();u0=data['u0'].copy()
  if str(data['model_sha256'].item())!=artifact['provenance']['model_sha256']:raise ValueError('Plant identity mismatch')
  if not np.isclose(float(data['target_mass_kg']),77.6122029,rtol=0,atol=1e-9):raise ValueError('Patient mass mismatch')
  if float(data['dt_s'])!=.01:raise ValueError('10ms plant required')
  if data['state_names'].tolist()!=policy.state_names or data['muscle_names'].tolist()!=policy.muscle_names:raise ValueError('Exact catalogs required')
  if not np.allclose(x0,policy.x0.numpy(),rtol=1e-6,atol=1e-7) or not np.allclose(u0,policy.u0.numpy(),rtol=1e-6,atol=1e-7):raise ValueError('Nominal baseline mismatch')
 # Freeze association weights in the artifact's actual dtype first, then evaluate equations in
 # float64 to avoid finite-difference roundoff masquerading as a zero derivative.
 source_dtype=str(policy.dyn.embed.dtype)
 with torch.no_grad():frozen_weights=policy.dyn.edge_weights().detach().double()
 policy=policy.double();dyn=policy.dyn;n=dyn.n;h=.001;state=dyn.init_state(1,'cpu')
 state=tuple(v.double() for v in state);zero=torch.zeros(1,n,dtype=torch.float64)
 with torch.no_grad():
  for _ in range(round(args.rest_seconds/h)):state=dyn.step(state,zero,h,frozen_weights)
  next_state=dyn.step(state,zero,h,frozen_weights)
 rest=np.concatenate([v[0].numpy() for v in state]);rest_delta=np.concatenate([(a-b)[0].numpy() for a,b in zip(next_state,state)])
 v,r,a,gi=[x[0].numpy() for x in state]
 eye=sparse.eye(n,format='csr');Z=sparse.csr_matrix((n,n));tm=dyn.tau_m;ta=dyn.tau_a
 scalar=lambda p:float(p.detach())
 row=np.repeat(np.arange(n),dyn.idx.shape[1]);col=dyn.idx.numpy().ravel()
 W=sparse.coo_matrix((frozen_weights.numpy().ravel(),(row,col)),shape=(n,n)).tocsr()
 rate=dyn.r_max/(1+np.exp(-(v-dyn.v_half)/dyn.slope));rate_derivative=rate*(1-rate/dyn.r_max)/dyn.slope
 J=sparse.bmat([
  [sparse.diags(1-h/tm*(1+np.maximum(gi,0))),h/tm*20/dyn.r_max*(scalar(dyn.w_ee)*eye+scalar(dyn.w_assoc)*W),-h/tm*eye,sparse.diags(-h/tm*(v-dyn.e_rev)*(gi>0))],
  [sparse.diags(h/.005*rate_derivative),(1-h/.005)*eye,Z,Z],
  [sparse.diags(h/ta*scalar(dyn.a_gain)*rate_derivative),Z,(1-h/ta)*eye,Z],
  [Z,h/.008*scalar(dyn.w_ei)/dyn.r_max*eye,Z,(1-h/.008)*eye]],format='csr')
 # Native x -> drive -> one neural tick. No teacher K enters this operator.
 E=np.zeros((n,len(x0)));E[policy.sensory_sites.numpy()]=12*policy.encoder.detach().numpy()/policy.state_scale.numpy()[None,:]
 D=np.zeros((4*n,len(x0)));D[:n]=h/tm*E
 decoder=policy.decoder.weight.detach().numpy();decoder=decoder-decoder.mean(axis=1,keepdims=True)
 C=np.zeros((len(u0),4*n));C[:,n+policy.motor_sites.numpy()]=2*decoder
 # Preserve actual anchored evaluation algebra, including mean roundoff.
 C[:,n+int(policy.motor_sites[0])]-=2*decoder.sum(axis=1)
 def neural(z,x):
  for _ in range(10):z=J@z+D@x
  return z
 def multiply(vector):
  x=vector[:len(x0)];z=vector[len(x0):];znew=neural(z,x)
  return np.concatenate((Ad@x+Bd@(C@znew),znew))
 size=len(x0)+4*n;operator=LinearOperator((size,size),matvec=multiply,dtype=np.float64)
 # Compare analytic one-tick and ten-tick JVPs against the retained equations.
 rng=np.random.default_rng(47);direction=rng.standard_normal(4*n);direction/=np.linalg.norm(direction)
 xdirection=rng.standard_normal(len(x0));xdirection/=np.linalg.norm(xdirection)
 def evaluate(z,x,ticks):
  s=tuple(torch.from_numpy(z[i*n:(i+1)*n][None].copy()) for i in range(4))
  drive=torch.from_numpy((E@x)[None])
  with torch.no_grad():
   for _ in range(ticks):s=dyn.step(s,drive,h,frozen_weights)
  return np.concatenate([q[0].numpy() for q in s])
 sanity=[]
 for ticks in (1,10):
  analytic=direction.copy()
  for _ in range(ticks):analytic=J@analytic+D@xdirection
  eps=1e-6
  finite=(evaluate(rest+eps*direction,eps*xdirection,ticks)-evaluate(rest-eps*direction,-eps*xdirection,ticks))/(2*eps)
  sanity.append({'ticks':ticks,'relative_error':float(np.linalg.norm(finite-analytic)/np.linalg.norm(analytic)),
    'max_absolute_error':float(abs(finite-analytic).max())})
 converged=True;message=None
 try:values,vectors=eigs(operator,k=8,which='LM',tol=1e-8,maxiter=args.maxiter,v0=rng.standard_normal(size),ncv=40)
 except ArpackNoConvergence as exc:
  converged=False;message=str(exc);values=exc.eigenvalues;vectors=exc.eigenvectors
 modes=[]
 for i,value in enumerate(values):
  vector=vectors[:,i];residual=np.linalg.norm(multiply(vector)-value*vector)/np.linalg.norm(vector)
  modes.append({'real':float(value.real),'imag':float(value.imag),'magnitude':float(abs(value)),
   'growth_per_s':float(np.log(abs(value))/.01),'frequency_hz':float(np.angle(value)/(2*np.pi*.01)),
   'relative_eigen_residual':float(residual),'native_state_norm_fraction':float(np.linalg.norm(vector[:len(x0)])/np.linalg.norm(vector))})
 modes.sort(key=lambda item:item['magnitude'],reverse=True)
 report={'schema':'ihm.cortical-native-local-linearization.v1','artifact_sha256':digest(artifact_path),
  'plant_sha256':digest(plant_path),'source_sha256':digest(__file__),'sites':n,'encoder_kind':policy.encoder_kind,
  'source_dtype':source_dtype,'diagnostic_dtype':'torch.float64',
  'augmented_dimension':size,'dt_s':.01,'neural_tick_s':h,'rest_seconds':args.rest_seconds,
  'rest_one_tick_max_abs_delta':float(abs(rest_delta).max()),'finite_difference_jvp':sanity,
  'arpack_converged':converged,'arpack_message':message,'maxiter':args.maxiter,'dominant_modes':modes,
  'unstable_mode_found':any(m['magnitude']>1+1e-7 for m in modes),
  'baseline_at_activation_clip_count':int(np.count_nonzero((u0<=.01+1e-8)|(u0>=1-1e-8))),
  'scope':['Local near-stationary paired-reference difference dynamics; four neural states per site',
   'Actual frozen association weights, encoder, decoder; no teacher K or hidden control enters augmented operator',
   'One 10ms neural update precedes native Ad/Bd update, matching current cortical command adapter ordering',
   'Unclipped derivative at nominal state; activation-boundary derivatives are directional and nonlinear excursions invalidate this model',
   'Native Ad/Bd already approximate local mechanics; no native rollout or standing/walking claim',
   'Finite differences use actual equations in float64 with weights frozen from source artifact dtype '+source_dtype,
   'Ten-tick constant-rest Jacobian approximation valid only to measured rest residual; reference common-mode perturbations excluded'],
  'wall_seconds':time.monotonic()-start}
 (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 (output/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
 print(json.dumps(report,indent=2),flush=True)
 if any(s['relative_error']>1e-4 for s in sanity):raise SystemExit('JVP sanity failed')

if __name__=='__main__':main()
