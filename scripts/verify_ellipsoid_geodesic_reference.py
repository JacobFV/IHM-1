"""Held native endpoints; independent stationary-route virtual-work checks."""
from pathlib import Path
import json,time
import numpy as np
from ellipsoid_geodesic_reference import solve_path
ROOT=Path(__file__).resolve().parents[1]
rows=[json.loads(s) for s in (ROOT/'data/derived/arm26-wrap-candidate-po_0h6uq/observed.jsonl').read_text().splitlines()]
start=time.monotonic();reports=[]
for side in ('r','l'):
    row=next(r for r in rows if r['case']==f'local_arm_add_{side}_0:arm26_BIClong_{side}' and r.get('wrapped'))
    p,q=np.array(row['p1']),np.array(row['p2']);seed=np.array(row['r1']);base=solve_path(p,q,seed,row['wrap_length']);errors=[]
    for endpoint in range(2):
        for axis in range(3):
            delta=np.zeros((2,3));delta[endpoint,axis]=1e-7;points=np.array([p,q])
            plus=solve_path(*(points+delta),seed,row['wrap_length']);minus=solve_path(*(points-delta),seed,row['wrap_length'])
            fd=(plus['length_m']-minus['length_m'])/2e-7;errors.append(abs(fd-base['endpoint_gradients'][endpoint,axis]))
    assert max(errors)<1e-7
    reports.append({'side':side,'length_m':base['length_m'],'arc_m':base['arc_m'],'maximum_endpoint_gradient_error':max(errors),'shooting_residual':base['shooting_residual_max'],'seed_shift_m':base['native_seed_contact_shift_m']})
print(json.dumps({'passed':True,'wall_s':time.monotonic()-start,'native_seed_cases':reports},indent=2))
