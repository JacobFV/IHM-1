"""Explicit local Newton boxes; no trust radius inferred from absolute q."""
import numpy as np
from scipy.optimize import lsq_linear


def interior_origin(q,bounds):
    q=np.asarray(q,float);bounds=np.asarray(bounds,float)
    if not np.all(np.isfinite(q)) or np.any(q<bounds[:,0]) or np.any(q>bounds[:,1]):
        raise ValueError('Initial root pose outside source bounds')
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
