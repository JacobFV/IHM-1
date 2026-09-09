#!/usr/bin/env python3
"""Audit exact linear pure-input-delay feedback, not neural/native acceptance."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np

def analyze(path, ticks=4):
    if type(ticks) is not int or ticks<1:raise ValueError('Positive integer delay ticks required')
    with np.load(path,allow_pickle=False) as d:A,B,K=d['Ad'],d['Bd'],d['K'];dt=float(d['dt_s'])
    n,m=B.shape;F=A-B@K
    base=np.zeros((n+ticks*m,n+ticks*m));base[:n,:n]=A;base[:n,n:n+m]=B
    for j in range(ticks-1):base[n+j*m:n+(j+1)*m,n+(j+1)*m:n+(j+2)*m]=np.eye(m)
    def modes(H):
        e=np.linalg.eigvals(H)
        return {'spectral_radius':float(max(abs(e))),'modes_above_1_plus_1e8':int(sum(abs(e)>1+1e-8))}
    result={}
    for name,G in [('K',K),('K_closed_loop_power',K@np.linalg.matrix_power(F,ticks)),('K_open_loop_power',K@np.linalg.matrix_power(A,ticks))]:
        H=base.copy();H[-m:,:n]=-G;result[name]=modes(H)
    H=base.copy();H[-m:,:n]=-K@np.linalg.matrix_power(A,ticks)
    for j in range(ticks):H[-m:,n+j*m:n+(j+1)*m]=-K@np.linalg.matrix_power(A,ticks-1-j)@B
    result['queued_history_predictor']=modes(H);result['undelayed']=modes(F)
    return {'artifact_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'delay_ticks':ticks,'dt_s':dt,'results':result,'scope':'Exact linear pure input delay, centered state/input, queue oldest to newest. Not a nonlinear or neural circuit guarantee.'}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--artifact',default='data/models/engineering_stance_v1/linearization.npz');p.add_argument('--ticks',type=int,default=4);p.add_argument('--output');a=p.parse_args();r=analyze(a.artifact,a.ticks);raw=json.dumps(r,indent=2)+'\n'
    if a.output:Path(a.output).write_text(raw)
    print(raw)
