#!/usr/bin/env python3
"""Offline exact discrete delay-margin analysis; no native acceptance implied."""
from pathlib import Path
import argparse,json,tempfile,hashlib
import numpy as np
from scipy.linalg import solve_discrete_are


def delayed_matrix(Ad,Bd,K,ticks):
 n,m=Bd.shape
 if ticks==0:return Ad-Bd@K
 out=np.zeros((n+m*ticks,n+m*ticks));out[:n,:n]=Ad;out[:n,n:n+m]=Bd
 for j in range(ticks-1):out[n+j*m:n+(j+1)*m,n+(j+1)*m:n+(j+2)*m]=np.eye(m)
 out[-m:,:n]=-K
 return out


def main():
 p=argparse.ArgumentParser();p.add_argument('--source',default='data/models/engineering_stance_v1/linearization.npz');p.add_argument('--ticks',type=int,default=4);a=p.parse_args()
 source=Path(a.source);z=np.load(source,allow_pickle=False);data={k:z[k] for k in z.files};Ad,Bd=data['Ad'],data['Bd'];baseQ,baseR=data['Q_discrete'],data['R_discrete'];names=data['state_names'];root=Path(__file__).resolve().parents[1]
 out=Path(tempfile.mkdtemp(prefix='delay_design_',dir=root/'data/research/locomotion_control'));rows=[]
 for state_weight in (1.,.1,0.):
  Q=baseQ.copy()
  for i,name in enumerate(names):
   if name.startswith('/forceset/'):Q[i,i]*=state_weight
  for scale in (1.,3.,10.,30.,100.,300.,1000.,3000.,10000.,100000.,1000000.):
   R=baseR*scale
   try:
    P=solve_discrete_are(Ad,Bd,Q,R);K=np.linalg.solve(R+Bd.T@P@Bd,Bd.T@P@Ad)
    rho=float(np.max(abs(np.linalg.eigvals(delayed_matrix(Ad,Bd,K,a.ticks)))))
    row={'muscle_state_cost_scale':state_weight,'R_scale':scale,'delay_ticks':a.ticks,'spectral_radius':rho,'candidate_linear_delay_stable':rho<=1+1e-8}
    if row['candidate_linear_delay_stable']:
     path=out/f'candidate_state{state_weight}_R{scale}.npz';np.savez_compressed(path,**dict(data,K=K,Q_discrete=Q,R_discrete=R,design_delay_ticks=a.ticks));row['artifact']=str(path.relative_to(root))
   except Exception as e:row={'muscle_state_cost_scale':state_weight,'R_scale':scale,'error':str(e)}
   rows.append(row);print(json.dumps(row),flush=True)
   (out/'report.json').write_text(json.dumps({'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'dt_s':float(data['dt_s']),'rows':rows,'limitations':['Linear time invariant pure input delay only; no nonlinear acceptance or nerve-loop equivalence implied.']},indent=2)+'\n')
 print('OUTPUT '+str(out),flush=True)
if __name__=='__main__':main()
