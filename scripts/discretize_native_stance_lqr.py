#!/usr/bin/env python3
"""Design sample-and-hold muscle excitation feedback from retained native A/B."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
from scipy.linalg import expm, solve_discrete_are

def main():
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('--dt',type=float,default=.01);p.add_argument('--margin-aware',action='store_true');a=p.parse_args()
 source=Path(a.source);original=np.load(source,allow_pickle=False);data={k:original[k] for k in original.files};A,B=data['A'],data['B'];n,m=B.shape
 if not np.isfinite(a.dt) or a.dt<=0:raise ValueError('Positive finite sampling interval required')
 transition=expm(np.block([[A,B],[np.zeros((m,n+m))]])*a.dt);Ad,Bd=transition[:n,:n],transition[:n,n:]
 Q=data['Q']*a.dt;R=data['R']*a.dt
 if a.margin_aware:
  if 'minimum_activation' not in data:raise ValueError('Margin-aware design requires native minimum activation values; relinearize this old artifact')
  R=np.diag(np.diag(R)*(.05/np.maximum(data['u0']-data['minimum_activation'],.005))**2)
 P=solve_discrete_are(Ad,Bd,Q,R);K=np.linalg.solve(R+Bd.T@P@Bd,Bd.T@P@Ad)
 report=json.loads(source.with_name('report.json').read_text());root=Path(__file__).resolve().parents[1];registration=json.loads((root/report['registration']).read_text())
 retained_model_sha256=hashlib.sha256((source.parent/'inputs/subject_walk_scaled.osim').read_bytes()).hexdigest()
 if retained_model_sha256!=registration['model_sha256']:raise ValueError('Retained linearization model differs from current registration')
 rho=lambda matrix:float(np.max(np.abs(np.linalg.eigvals(matrix))))
 report['continuous_lqr_solved']=report.get('lqr_solved',False)
 if 'lqr_error' in report:report['continuous_lqr_error']=report.pop('lqr_error')
 report['lqr_solved']=True
 report['continuous_design_max_real_eigenvalue']=report.pop('closed_loop_max_real_eigenvalue',None)
 report.update(schema='ihm.native-stance-discrete-lqr.v1',dt_s=a.dt,margin_aware=a.margin_aware,model_sha256=retained_model_sha256,source_artifact_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),continuous_gain_sampled_spectral_radius=rho(Ad-Bd@data['K']),discrete_gain_spectral_radius=rho(Ad-Bd@K),discrete_gain_max=float(abs(K).max()),dare_residual_relative=float(np.linalg.norm(Ad.T@P@Ad-P-Ad.T@P@Bd@K+Q)/np.linalg.norm(Q)))
 if not report['continuous_lqr_solved']:report['continuous_gain_sampled_spectral_radius']=None
 out=source.parent/('discrete_margin' if a.margin_aware else 'discrete');out.mkdir(exist_ok=True)
 data.update(K=K,Ad=Ad,Bd=Bd,Q_discrete=Q,R_discrete=R,dt_s=np.array(a.dt),model_sha256=np.array(registration['model_sha256']))
 np.savez_compressed(out/'linearization.npz',**data);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'output':str(out),**report},indent=2))
if __name__=='__main__':main()
