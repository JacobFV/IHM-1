#!/usr/bin/env python3
"""Build an isolated exact-native finite-difference linearization and LQR gain.

The native initialization is snapshotted from native_mechanical_stream.cpp;
no active runtime executable or shared build pointer is changed.
"""
from pathlib import Path
import argparse, hashlib, json, os, re, shlex, shutil, subprocess, sys, tempfile
import numpy as np
from scipy.linalg import solve_continuous_are
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import SOURCE_FILES

def main():
 p=argparse.ArgumentParser();p.add_argument('--registration',default='data/derived/mechanics/resting_stance86_refined/registration.json');p.add_argument('--output');p.add_argument('--target-mass-kg',type=float,default=70.0);p.add_argument('--reference-state-json');a=p.parse_args()
 if not np.isfinite(a.target_mass_kg) or a.target_mass_kg<=0:raise ValueError('Positive finite target mass required')
 out=Path(a.output).resolve() if a.output else Path(tempfile.mkdtemp(prefix='linearization_',dir=ROOT/'data/research/locomotion_control'))
 out.mkdir(parents=True,exist_ok=True);inputs=out/'inputs';inputs.mkdir()
 reg=json.loads((ROOT/a.registration).read_text());original=ROOT/'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
 pose=reg.get('initial_pose')
 if pose is not None:
  if not isinstance(pose,dict) or not pose or any(not isinstance(name,str) or re.fullmatch(r'[A-Za-z0-9_]+',name) is None or isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) for name,value in pose.items()):raise ValueError('Invalid explicit registration initial pose')
  (inputs/'initial_pose.txt').write_text('IHM_INITIAL_POSE_V1 '+str(len(pose))+'\n'+''.join(name+' '+str(value)+'\n' for name,value in pose.items()))
 for name in SOURCE_FILES:shutil.copyfile(ROOT/reg['model_path'] if name=='subject_walk_scaled.osim' else original/name,inputs/name)
 reference_sha256=None
 if a.reference_state_json:
  reference_bytes=Path(a.reference_state_json).read_bytes();reference=json.loads(reference_bytes);variables=reference.get('state_variables')
  if 'source_model_sha256' in reference and reference['source_model_sha256']!=reg['model_sha256']:raise ValueError('Reference state source model mismatch')
  if 'source_registration_sha256' in reference and reference['source_registration_sha256']!=hashlib.sha256((ROOT/a.registration).read_bytes()).hexdigest():raise ValueError('Reference state registration mismatch')
  if 'target_mass_kg' in reference and (isinstance(reference['target_mass_kg'],bool) or not np.isclose(reference['target_mass_kg'],a.target_mass_kg,rtol=0,atol=1e-9)):raise ValueError('Reference state mass mismatch')
  if not isinstance(variables,dict) or not variables or any(not isinstance(name,str) or not name.startswith('/') or any(c.isspace() for c in name) or isinstance(value,bool) or not isinstance(value,(float,int)) or not np.isfinite(value) for name,value in variables.items()):raise ValueError('Reference JSON requires finite named state_variables')
  reference_sha256=hashlib.sha256(reference_bytes).hexdigest()
  (inputs/'reference_state.json').write_bytes(reference_bytes)
  (inputs/'reference_state.txt').write_text('IHM_REFERENCE_STATE_V1 '+str(len(variables))+'\n'+''.join(name+' '+str(value)+'\n' for name,value in variables.items()))
 native=(ROOT/'scripts/native_mechanical_stream.cpp').read_text();prefix=native[:native.index(' const char* mass_mode=')]
 source=out/'linearize.cpp';source.write_text(prefix+(ROOT/'scripts/native_stance_linearization_tail.inc').read_text())
 runtime=ROOT/'data/runtime/opensim';flags=['-std=c++20','-O0','-DSWIG_PYTHON','-I'+str(ROOT/'scripts')]
 for d in ['install/opensim/include','install/opensim/include/OpenSim','install/simbody/include/simbody']:flags+=['-isystem',str(runtime/d)]
 previous=shlex.split((runtime/'dynamics-adapter/build/CMakeFiles/native_opensim_dynamics.dir/link.txt').read_text());libs=[v for v in previous if v.endswith('.so') or '.so.' in v or v.startswith(('-Wl,','-l'))]
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 with (out/'compile.log').open('w') as log:subprocess.run(['c++',*flags,str(source),'-o',str(out/'linearize'),*libs,'-ldl'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 with (out/'native.log').open('w') as log:subprocess.run([str(out/'linearize'),str(inputs),str(out),'upright',str(a.target_mass_kg)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 raw=json.loads((out/'linearization.json').read_text());X=np.array([c['dx'] for c in raw['state_columns']]).T;F=np.array([c['df'] for c in raw['state_columns']]).T
 if a.reference_state_json and 'equilibrium_excitations' in reference:
  excitation=reference['equilibrium_excitations']
  if set(excitation)!=set(raw['muscle_names']) or any(abs(excitation[name]-value)>1e-12 for name,value in zip(raw['muscle_names'],raw['u0'])):raise ValueError('Native default excitation differs from reference equilibrium input')
 A=np.linalg.solve(X.T,F.T).T;B=np.array(raw['input_columns']).T;n=len(A);m=B.shape[1]
 q=np.ones(n)
 for i,name in enumerate(raw['state_names']):
  if name.endswith('/value'):q[i]=1000 if '/pelvis_' in name else 100
  elif name.endswith('/speed'):q[i]=10
  elif name.endswith('/fiber_length'):q[i]=100
 Q=np.diag(q);R=np.eye(m)*10
 report={'schema':'ihm.native-stance-linearization.v1','registration':a.registration,'state_count':n,'muscle_count':m,'tangent_condition':float(np.linalg.cond(X)),'base_derivative_norm':float(np.linalg.norm(raw['f0'])),'open_loop_max_real_eigenvalue':float(np.max(np.linalg.eigvals(A).real)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'limitations':['Local linearization at initialized native state; no nonlinear balance or walking acceptance.','Excitations require clipping, which invalidates unconstrained LQR guarantees.']}
 try:
  P=solve_continuous_are(A,B,Q,R);K=np.linalg.solve(R,B.T@P);closed=np.linalg.eigvals(A-B@K)
  report.update(lqr_solved=True,closed_loop_max_real_eigenvalue=float(np.max(closed.real)),care_residual_relative=float(np.linalg.norm(A.T@P+P@A-P@B@np.linalg.solve(R,B.T)@P+Q)/np.linalg.norm(Q)))
 except Exception as e:K=np.zeros((m,n));report.update(lqr_solved=False,lqr_error=str(e))
 report['minimum_activation_basis']='Native Millard2012EquilibriumMuscle/Thelen2003Muscle.getMinimumActivation()'
 report['minimum_activation_range']=[min(raw['minimum_activation']),max(raw['minimum_activation'])]
 report['target_mass_kg']=a.target_mass_kg
 report['initial_pose']=pose
 report['reference_state_json']=a.reference_state_json
 report['reference_state_sha256']=reference_sha256
 np.savez_compressed(out/'linearization.npz',A=A,B=B,K=K,Q=Q,R=R,x0=raw['x0'],u0=raw['u0'],minimum_activation=raw['minimum_activation'],state_names=raw['state_names'],muscle_names=raw['muscle_names'],target_mass_kg=a.target_mass_kg,reference_state_sha256=reference_sha256 or '')
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'output':str(out),**report},indent=2))
if __name__=='__main__':main()
