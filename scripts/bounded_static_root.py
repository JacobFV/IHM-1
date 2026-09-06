"""Explicit local Newton boxes; no trust radius inferred from absolute q."""
import numpy as np
from scipy.optimize import lsq_linear


def interior_origin(q,bounds,*,preserve=False):
    q=np.asarray(q,float);bounds=np.asarray(bounds,float)
    if not np.all(np.isfinite(q)) or np.any(q<bounds[:,0]) or np.any(q>bounds[:,1]):
        raise ValueError('Initial root pose outside source bounds')
    if preserve:return q.copy()
    margin=np.minimum(1e-12,np.diff(bounds,axis=1)[:,0]/4)
    return np.clip(q,bounds[:,0]+margin,bounds[:,1]-margin)


def local_step(jacobian,acceleration,q,bounds,radius=.03):
    q=np.asarray(q,float);bounds=np.asarray(bounds,float)
    if not 0<radius<=.03:raise ValueError('Invalid local coordinate radius')
    lo=np.maximum(bounds[:,0]-q,-radius);hi=np.minimum(bounds[:,1]-q,radius)
    result=lsq_linear(jacobian,-np.asarray(acceleration,float),bounds=(lo,hi),method='bvls',tol=1e-10,max_iter=100)
    if not result.success or not np.all(np.isfinite(result.x)):raise ValueError('Bounded linear root step failed')
    if np.any(result.x<lo-1e-12) or np.any(result.x>hi+1e-12):raise ValueError('Local linear step violated bounds')
    return np.clip(result.x,lo,hi)


class StaticDomainRejection(ValueError):
    """Exact copied-state skin/bed domain rejection; never an accepted response."""


def backtracked_trial(q,direction,bounds,current_cost,evaluate):
    for backtrack in range(6):
        candidate=np.clip(q+direction*(.5**backtrack),bounds[:,0],bounds[:,1])
        try:trial=evaluate(candidate)
        except StaticDomainRejection:continue
        if trial['cost']<current_cost:return candidate,trial
    return None,None


def constrained_local_step(jacobian,acceleration,q,bounds,support_jacobian,support,radius=.03):
    from scipy.optimize import minimize,LinearConstraint,Bounds
    q=np.asarray(q,float);bounds=np.asarray(bounds,float);A=np.asarray(support_jacobian,float);support=np.asarray(support,float)
    if not 0<radius<=.03:raise ValueError('Invalid constrained local radius')
    lo=np.maximum(bounds[:,0]-q,-radius);hi=np.minimum(bounds[:,1]-q,radius)
    # Constant scaling conditions the numerical QP; it does not change its minimizer.
    scale=max(1.,float(np.linalg.norm(acceleration)));B=np.asarray(jacobian,float)/scale;a=np.asarray(acceleration,float)/scale
    row_scale=np.linalg.norm(A,axis=1)
    if np.any(row_scale<=1e-12):raise ValueError('Rank-deficient support equality row')
    C=A/row_scale[:,None];target=-support/row_scale
    result=minimize(lambda d:.5*float((a+B@d)@(a+B@d)),np.zeros(len(q)),jac=lambda d:B.T@(a+B@d),
        method='SLSQP',bounds=Bounds(lo,hi),constraints=[LinearConstraint(C,target,target)],
        options={'maxiter':150,'ftol':1e-12})
    residual=support+A@result.x
    if not result.success or np.max(np.abs(residual))>1e-8 or np.any(result.x<lo-1e-10) or np.any(result.x>hi+1e-10):
        raise ValueError('Bounded support-constrained linear step failed: '+str(result.message))
    return result.x,dict(iterations=int(result.nit),linear_support_max=float(np.max(np.abs(residual))),
        predicted_acceleration_norm=float(np.linalg.norm(np.asarray(acceleration)+np.asarray(jacobian)@result.x)),
        maximum_coordinate_step=float(np.max(np.abs(result.x))))


def balanced_backtracked_trial(q,direction,bounds,current_cost,evaluate,support_jacobian,root_indices,*,radius=.03,maximum_backtracks=6):
    if not 0<radius<=.03 or not 1<=maximum_backtracks<=6:raise ValueError("Invalid support trial box or attempt cap")
    root_block=np.asarray(support_jacobian)[:,root_indices]
    if not np.isfinite(np.linalg.cond(root_block)) or np.linalg.cond(root_block)>1e10:
        raise ValueError('Ill-conditioned support-coordinate correction')
    lo=np.maximum(bounds[:,0],q-radius);hi=np.minimum(bounds[:,1],q+radius)
    for backtrack in range(maximum_backtracks):
        candidate=np.clip(q+direction*(.5**backtrack),lo,hi)
        try:
            for correction in range(4):
                trial=evaluate(candidate)
                if max(map(abs,trial['support_constraints']+trial['gauge_residual']))<=1e-4:
                    if trial['cost']<current_cost:return candidate,trial
                    break
                if correction==3:break
                change=np.linalg.solve(root_block,-np.asarray(trial['support_constraints']))
                corrected=candidate.copy();corrected[root_indices]+=change;corrected=np.clip(corrected,lo,hi)
                if np.array_equal(corrected,candidate):break
                candidate=corrected
        except StaticDomainRejection:continue
    return None,None


def recomputed_box_trial(jacobian,acceleration,q,bounds,support_jacobian,support,current_cost,evaluate,root_indices,radius):
    """At most six local QPs/four evaluations each; no physical tolerance changes."""
    trace=[]
    for attempt in range(6):
        trial_radius=radius*(.5**attempt)
        direction,diagnostic=constrained_local_step(jacobian,acceleration,q,bounds,support_jacobian,support,radius=trial_radius)
        candidate,entry=balanced_backtracked_trial(q,direction,bounds,current_cost,evaluate,support_jacobian,root_indices,
            radius=trial_radius,maximum_backtracks=1)
        trace.append(dict(radius=trial_radius,accepted=candidate is not None,**diagnostic))
        if candidate is not None:return candidate,entry,trial_radius,trace
    return None,None,radius*(.5**5),trace
