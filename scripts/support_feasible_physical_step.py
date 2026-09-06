"""Physical acceleration metric with explicit native support feasibility."""
import numpy as np
from scipy.optimize import minimize,LinearConstraint,Bounds
from verify_native_physical_metric import audit_metric


def native_metric(response):
    audit=audit_metric(response)
    if not audit['passed']:raise ValueError('Native inverse-mass/sign/constraint gate failed')
    return -np.asarray(response['inverse_mass_matrix'])@np.asarray(response['constrained_residual'])


def physical_step(acceleration,jacobian,q,bounds,support,support_jacobian,radius):
    a=np.asarray(acceleration,float);B=np.asarray(jacobian,float);q=np.asarray(q,float);bounds=np.asarray(bounds,float)
    A=np.asarray(support_jacobian,float);s=np.asarray(support,float)
    if not 0<radius<=.03:raise ValueError('Invalid native-coordinate trust box')
    lo=np.maximum(bounds[:,0]-q,-radius);hi=np.minimum(bounds[:,1]-q,radius)
    # Column preconditioning only: no residual weighting can suppress support.
    norms=np.linalg.norm(B,axis=0)
    if np.any(norms==0):raise ValueError('Unobservable physical step column')
    cost_scale=max(1.,np.linalg.norm(a));c=cost_scale/norms;C=B*c/cost_scale;b=a/cost_scale
    row=np.linalg.norm(A*c,axis=1)
    if np.any(row==0):raise ValueError('Unresponsive support constraint')
    E=A*c/row[:,None];target=-s/row
    fit=minimize(lambda z:.5*np.sum((b+C@z)**2),np.zeros(len(q)),jac=lambda z:C.T@(b+C@z),
        method='SLSQP',bounds=Bounds(lo/c,hi/c),constraints=[LinearConstraint(E,target,target)],
        options=dict(maxiter=150,ftol=1e-12))
    d=fit.x*c
    if not fit.success or np.max(np.abs(s+A@d))>1e-8 or np.any(d<lo-1e-10) or np.any(d>hi+1e-10):raise ValueError('Support-feasible physical QP failed: '+fit.message)
    singular=np.linalg.svd(C,compute_uv=False)
    return d,dict(predicted_acceleration_cost=float(.5*np.sum((a+B@d)**2)),initial_acceleration_cost=float(.5*a@a),
        maximum_linear_support_error=float(np.max(np.abs(s+A@d))),singular_values=singular.tolist(),
        rank=int(np.linalg.matrix_rank(C)),maximum_coordinate_step=float(np.max(abs(d))),
        active_source_or_box_columns=np.flatnonzero((abs(d-lo)<1e-9)|(abs(d-hi)<1e-9)).tolist(),
        coordinate_preconditioner=c.tolist(),radius=radius)


def feasible_trial(current,candidate,predicted_reduction):
    # Both dictionaries must retain the original physical-gate quantities.
    feasible=max(map(abs,candidate['support_constraints']+candidate['gauge_residual']))<=1e-4
    actual=.5*(np.sum(np.asarray(current['udot'])**2)-np.sum(np.asarray(candidate['udot'])**2))
    ratio=actual/predicted_reduction if predicted_reduction>0 else None
    return dict(feasible=feasible,actual_reduction=float(actual),reduction_ratio=ratio,
        accepted=bool(feasible and predicted_reduction>0 and actual>0 and ratio>=.1))
