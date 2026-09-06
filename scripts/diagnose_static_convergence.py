"""Read-only retained-Jacobian convergence audit, never a native solver."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from bounded_static_root import constrained_local_step, balanced_backtracked_trial
from static_pose_journal import pose_key


def linear_diagnostic(J,a,A,c):
    Z=null_space(A)
    particular=np.linalg.lstsq(A,-c,rcond=None)[0]
    B=J@Z
    b=a+J@particular
    y,_,rank,singular=np.linalg.lstsq(B,-b,rcond=None)
    return dict(projected_half_cost_gradient_norm=float(np.linalg.norm(Z.T@J.T@a)),
                linearized_equality_residual_floor_norm=float(np.linalg.norm(b+B@y)),
                tangent_rank=int(rank),tangent_dimension=int(Z.shape[1]),
                tangent_singular_values=singular.tolist(),
                unconstrained_tangent_step_max=float(np.max(np.abs(particular+Z@y))))


def audit(path):
    identity=json.loads((path/'cache_identity.json').read_text())
    cache={}
    for line in (path/'candidates.jsonl').read_text().splitlines():
        row=json.loads(line);payload=row['payload']
        # Journal's serialization/hash contract is independently checked below.
        from static_pose_journal import digest
        if digest(payload)!=row['sha256']:raise ValueError('Journal record hash mismatch')
        cache[pose_key(payload['values'])]=payload['entry']
    names=identity['coordinate_order'];bounds=np.array(identity['bounds'])
    accepted=[json.loads(x) for x in (path/'optimizer_iterates.jsonl').read_text().splitlines()]
    rows=[]
    for index,record in enumerate(accepted):
        q=np.array(record['values']);entry=cache[pose_key(q)];a=np.array(entry['native']['udot']);c=np.array(entry['support_constraints']);J=[];A=[]
        for k,(_,hi) in enumerate(bounds):
            h=1e-5 if q[k]+1e-5<=hi else -1e-5;s=q.copy();s[k]+=h
            probe=cache.get(pose_key(s))
            if probe is None:break
            J.append((np.array(probe['native']['udot'])-a)/h)
            A.append((np.array(probe['support_constraints'])-c)/h)
        if len(J)!=len(q):continue
        J=np.ascontiguousarray(np.array(J).T);A=np.ascontiguousarray(np.array(A).T)
        d=linear_diagnostic(J,a,A,c)
        step,_=constrained_local_step(J,a,q,bounds,A,c)
        d['smaller_box_resolves']=[]
        for radius in (.015,.0075):
            small,_=constrained_local_step(J,a,q,bounds,A,c,radius=radius)
            d['smaller_box_resolves'].append(dict(radius=radius,predicted_cost=float(np.linalg.norm(a+J@small)**2),
                scaled_original_predicted_cost=float(np.linalg.norm(a+J@step*(radius/.03))**2),
                direction_difference_norm=float(np.linalg.norm(small-step*(radius/.03)))))
        d.update(accepted_index=index,cost=entry['cost'],
            bounded_predicted_cost=float(np.linalg.norm(a+J@step)**2),
            local_box_active_coordinates=[names[k] for k in range(len(q)) if abs(abs(step[k])-.03)<1e-7],
            source_bound_active_coordinates=[names[k] for k in range(len(q)) if min(abs(q[k]+step[k]-bounds[k,0]),abs(q[k]+step[k]-bounds[k,1]))<1e-7])
        if index+1<len(accepted):
            nextq=np.array(accepted[index+1]['values']);nextentry=cache[pose_key(nextq)]
            pred=a+J@(nextq-q)
            trace=[]
            def cached_trial(candidate):
                observed=cache.get(pose_key(candidate))
                if observed is None:raise ValueError('Exact rejected-trial cache replay missing; no nearest-neighbor substitution')
                support_max=max(map(abs,observed['support_constraints']))
                gauge_max=max(map(abs,observed['gauge_residual']))
                trace.append(dict(cost=observed['cost'],maximum_coordinate_step=float(np.max(np.abs(candidate-q))),
                    support_max=support_max,gauge_max=gauge_max,
                    support_filter_pass=bool(max(support_max,gauge_max)<=1e-4),objective_decreases=bool(observed['cost']<entry['cost']),
                    source_bounds_pass=bool(np.all(candidate>=bounds[:,0]) and np.all(candidate<=bounds[:,1])),
                    largest_accelerations=sorted(zip(observed['native']['mobility_coordinate_names'],observed['native']['udot']),key=lambda item:-abs(item[1]))[:5],
                    coordinate_sha256=hashlib.sha256(json.dumps(candidate.tolist()).encode()).hexdigest()))
                return observed
            replay_q,replay_entry=balanced_backtracked_trial(q,step,bounds,entry['cost'],cached_trial,A,
                [names.index(n) for n in ('pelvis_tx','pelvis_tilt','pelvis_rotation')])
            if replay_q is None or not np.array_equal(replay_q,nextq):raise ValueError('Cached trial replay did not reproduce accepted q exactly')
            d['exact_trial_replay']=trace
            d['actual_next_cost']=nextentry['cost']
            d['actual_next_maximum_coordinate_step']=float(np.max(np.abs(nextq-q)))
            d['actual_next_predicted_cost']=float(pred@pred)
            d['actual_next_prediction_error_norm']=float(np.linalg.norm(np.array(nextentry['native']['udot'])-pred))
            active=lambda e:{v['name'] for v in e['native']['contacts'] if np.linalg.norm(v['force_n'])>0}
            d['active_body_contacts_added']=sorted(active(nextentry)-active(entry))
            d['active_body_contacts_removed']=sorted(active(entry)-active(nextentry))
        rows.append(d)
    return dict(scope='Offline local linearization of archived actual evaluations; not a global minimum or physical convergence certificate',
        coordinate_metric='Native q units (radians/metres), matching solver unit acceleration scaling; singular values and gradient depend on this declared metric',
        contact_resolution='Body aggregate contact presence only; per-quadrature active sets were not emitted and cannot be reconstructed from these receipts',
        source_sha256={str(path/f):hashlib.sha256((path/f).read_bytes()).hexdigest() for f in ('cache_identity.json','candidates.jsonl','optimizer_iterates.jsonl')},
        diagnostics=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('receipt',type=Path);args=p.parse_args()
    print(json.dumps(audit(args.receipt),indent=2,allow_nan=False))
