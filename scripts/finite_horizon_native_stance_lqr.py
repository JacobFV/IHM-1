#!/usr/bin/env python3
"""Finite-horizon Riccati fallback for neutral-mode DARE numerical failures.

Same native ZOH dynamics and margin-aware stage costs. This is a finite-horizon
receding feedback design, not a claim to have solved the algebraic Riccati eq.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from scipy.linalg import expm

def main():
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('--horizon-s',type=float,default=30.);p.add_argument('--dt',type=float,default=.01);a=p.parse_args()
 if not np.isfinite(a.dt) or a.dt<=0 or not np.isfinite(a.horizon_s) or a.horizon_s<=0:raise ValueError('Finite positive timestep and horizon required')
 source=Path(a.source);z=np.load(source,allow_pickle=False);data={k:z[k] for k in z.files};A,B=data['A'],data['B'];n,m=B.shape;E=expm(np.block([[A,B],[np.zeros((m,n+m))]])*a.dt);Ad,Bd=E[:n,:n],E[:n,n:]
 Q=data['Q']*a.dt;R=np.diag(np.diag(data['R']*a.dt)*(.05/np.maximum(data['u0']-data['minimum_activation'],.005))**2)
 P=Q.copy();steps=int(np.ceil(a.horizon_s/a.dt));samples=[];previous=None
 for step in range(1,steps+1):
  K=np.linalg.solve(R+Bd.T@P@Bd,Bd.T@P@Ad);F=Ad-Bd@K;P=Q+K.T@R@K+F.T@P@F;P=(P+P.T)*.5
  if step%500==0 or step==steps:
   row={'steps':step,'spectral_radius':float(max(abs(np.linalg.eigvals(F)))),'gain_max':float(abs(K).max()),'gain_change_relative':None if previous is None else float(np.linalg.norm(K-previous)/np.linalg.norm(K))};samples.append(row);previous=K.copy();print(json.dumps(row),flush=True)
 root=Path(__file__).resolve().parents[1];report=json.loads(source.with_name('report.json').read_text());reg=json.loads((root/report['registration']).read_text());retained=hashlib.sha256((source.parent/'inputs/subject_walk_scaled.osim').read_bytes()).hexdigest()
 if retained!=reg['model_sha256']:raise ValueError('Retained model mismatch')
 report.update(schema='ihm.native-stance-finite-horizon-lqr.v1',design='finite-horizon receding feedback, terminal cost Q',horizon_steps=steps,dt_s=a.dt,margin_aware=True,discrete_gain_spectral_radius=samples[-1]['spectral_radius'],samples=samples,model_sha256=retained,source_artifact_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
 out=source.parent/'finite_horizon_margin';out.mkdir(exist_ok=False);data.update(K=K,Ad=Ad,Bd=Bd,Q_discrete=Q,R_discrete=R,dt_s=np.array(a.dt),model_sha256=np.array(retained),horizon_steps=np.array(steps));np.savez_compressed(out/'linearization.npz',**data);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('OUTPUT '+str(out),flush=True)
if __name__=='__main__':main()
