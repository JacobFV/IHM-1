"""Offline source-bounded descent and complete passive-component inventory."""
import argparse,json,hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import linprog
from bounded_static_root import constrained_local_step
from static_pose_journal import pose_key,digest
from build_supine_initial_state import passive_expression


def tangent_descent(gradient,A,q,bounds,radius):
    lo=np.maximum(bounds[:,0]-q,-radius);hi=np.minimum(bounds[:,1]-q,radius)
    result=linprog(gradient,A_eq=A,b_eq=np.zeros(len(A)),bounds=list(zip(lo,hi)),method='highs')
    if not result.success:raise ValueError('Source-bounded tangent LP failed')
    return dict(half_cost_directional_derivative=float(gradient@result.x),linear_support_tangent_error=float(np.max(np.abs(A@result.x))),step=result.x.tolist())


def minimax_step(J,a,A,c,q,bounds,radius):
    n=len(q);lo=np.maximum(bounds[:,0]-q,-radius);hi=np.minimum(bounds[:,1]-q,radius)
    objective=np.r_[np.zeros(n),1.]
    inequalities=np.vstack([np.c_[J,-np.ones(len(a))],np.c_[-J,-np.ones(len(a))]])
    result=linprog(objective,A_ub=inequalities,b_ub=np.r_[-a,a],A_eq=np.c_[A,np.zeros(len(A))],b_eq=-c,
        bounds=list(zip(lo,hi))+[(0,None)],method='highs')
    if not result.success:raise ValueError('Bounded minimax diagnostic failed')
    return result.x[:-1]


def audit(path):
    identity=json.loads((path/'cache_identity.json').read_text());names=identity['coordinate_order'];bounds=np.array(identity['bounds']);cache={}
    for line in (path/'candidates.jsonl').read_text().splitlines():
        row=json.loads(line);p=row['payload']
        if digest(p)!=row['sha256']:raise ValueError('Changed cached response')
        cache[pose_key(p['values'])]=p['entry']
    records=[json.loads(s) for s in (path/'optimizer_iterates.jsonl').read_text().splitlines()]
    latest=json.loads((path/'last_optimizer_iterate.json').read_text());n=latest['native'];results=[]
    for record in records:
        q=np.array(record['values']);entry=cache[pose_key(q)];a=np.array(entry['native']['udot']);c=np.array(entry['support_constraints']);J=np.zeros((len(a),len(q)));A=np.zeros((3,len(q)))
        for k,(_,hi) in enumerate(bounds):
            h=1e-5 if q[k]+1e-5<=hi else -1e-5;s=q.copy();s[k]+=h;probe=cache.get(pose_key(s))
            if probe is None:break
            J[:,k]=(np.array(probe['native']['udot'])-a)/h;A[:,k]=(np.array(probe['support_constraints'])-c)/h
        else:
            hip=names.index('hip_rotation_r');row=np.zeros(len(q));row[hip]=1
            free,_=constrained_local_step(J,a,q,bounds,A,c,radius=.00375)
            fixed,_=constrained_local_step(J,a,q,bounds,np.vstack([A,row]),np.r_[c,0.],radius=.00375)
            tangent=tangent_descent(J.T@a,A,q,bounds,.00375)
            minimax=minimax_step(J,a,A,c,q,bounds,.00375)
            results.append(dict(cost=entry['cost'],hip_q=q[hip],source_bound_qp_predicted_cost=float(np.linalg.norm(a+J@free)**2),
                current_maximum_acceleration=float(np.max(np.abs(a))),
                qp_predicted_maximum_acceleration=float(np.max(np.abs(a+J@free))),
                minimax_predicted_maximum_acceleration=float(np.max(np.abs(a+J@minimax))),
                minimax_predicted_cost=float(np.linalg.norm(a+J@minimax)**2),
                minimax_step=minimax.tolist(),
                hip_fixed_qp_predicted_cost=float(np.linalg.norm(a+J@fixed)**2),hip_step=float(free[hip]),
                source_bound_tangent=tangent,largest_qp_steps=sorted(zip(names,free.tolist()),key=lambda v:-abs(v[1]))[:8]))
    model=path/'native/assembled_model.osim';root=ET.parse(model).getroot();passive=[]
    for force in root.iter('ExpressionBasedCoordinateForce'):
        name=force.findtext('coordinate');q=n['coordinates'][name]['value'];expr=force.findtext('expression')
        passive.append(dict(component=force.get('name'),coordinate=name,q=q,expression=expr,zero_speed_torque_nm=passive_expression(expr,q)))
    conventions={}
    for joint in root.iter('CustomJoint'):
        axes=[dict(coordinate=ax.findtext('coordinates'),axis=ax.findtext('axis')) for ax in joint.findall('./SpatialTransform/TransformAxis') if ax.findtext('coordinates')]
        if any(v['coordinate'].startswith(('arm_','lumbar_','hip_')) for v in axes):conventions[joint.get('name')]=axes
    torques=sorted(zip(n['mobility_coordinate_names'],n['constrained_zero_acceleration_residual_mobility_force']),key=lambda v:-abs(v[1]))
    return dict(scope='Offline actual-cache linearization, not nonlinear feasible-step verification or physical acceptance',
        latest_cost=latest['cost'],dominant_generalized_residuals=torques,passive_coordinate_forces=passive,
        coordinate_axes_in_joint_frames=conventions,complete_jacobians=results,
        correction='20 passive expression components exist outside ForceSet; ForceSet-only inventory is incomplete. No discrete hip ROM stop exists, but a retained smooth hip passive torque is active.',
        source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (model,path/'candidates.jsonl',path/'last_optimizer_iterate.json',Path(__file__))})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('receipt',type=Path);a=p.parse_args();print(json.dumps(audit(a.receipt),indent=2,allow_nan=False))
