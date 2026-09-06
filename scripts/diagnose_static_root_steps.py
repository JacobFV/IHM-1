"""Native-free comparison of bounded solvers on a retained local Jacobian."""
from pathlib import Path
import argparse,json,sys
import numpy as np
from scipy.optimize import least_squares
from static_pose_journal import digest


def compare(directory):
    directory=Path(directory);identity=json.loads((directory/'cache_identity.json').read_text())
    records=[]
    for line in (directory/'candidates.jsonl').read_text().splitlines():
        record=json.loads(line);assert digest(record['payload'])==record['sha256'];records.append(record['payload'])
    # TRF's first interior-adjusted point follows the explicit seed observation.
    base=records[1];q=np.array(base['values']);a=np.array(base['entry']['native']['udot']);columns={}
    for record in records[2:]:
        delta=np.array(record['values'])-q;changed=np.flatnonzero(delta)
        if len(changed)==1 and int(changed[0]) not in columns:
            i=int(changed[0]);columns[i]=(np.array(record['entry']['native']['udot'])-a)/delta[i]
        if len(columns)==len(q):break
    if len(columns)!=len(q):raise ValueError('No complete first Jacobian in journal')
    matrix=np.column_stack([columns[i] for i in range(len(q))]);bounds=np.array(identity['bounds'])
    results={}
    for origin_name,origin in [('requested',np.array(records[0]['values'])),('interior',q)]:
        for method in ['trf','dogbox']:
            history=[];offset=a+matrix@(origin-q)
            result=least_squares(lambda d:offset+matrix@d,np.zeros(len(q)),jac=lambda d:matrix,
                bounds=(bounds[:,0]-origin,bounds[:,1]-origin),method=method,x_scale=.03,
                ftol=None,xtol=None,gtol=1e-10,max_nfev=6,callback=lambda d:history.append(d.tolist()))
            results[origin_name+'_'+method]={'residual_norm':float(np.linalg.norm(result.fun)),'evaluations':result.nfev,
                 'maximum_q_change':float(np.max(np.abs(result.x))),'first_step_max':None if not history else float(np.max(np.abs(history[0])))}
    return dict(native_run=False,scope='Frozen first-Jacobian linear surrogate only; no actual changed-pose physics or material-domain verification',
        source_cache_identity=digest(identity),initial_residual_norm=float(np.linalg.norm(a)),comparison=results)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('directory',type=Path);args=parser.parse_args()
    print(json.dumps(compare(args.directory),indent=2))
