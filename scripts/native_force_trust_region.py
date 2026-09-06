"""Bounded generalized-force root steps; numerical merit is not stored energy."""
import numpy as np
from scipy.optimize import lsq_linear


def sensitivity_scales(jacobian, bounds):
    """Two-sided equilibration from source spans and actual force sensitivities."""
    J=np.asarray(jacobian,float);bounds=np.asarray(bounds,float);span=bounds[:,1]-bounds[:,0]
    if not np.all(np.isfinite(J)) or np.any(span<=0):raise ValueError('Invalid source Jacobian/bounds')
    row=np.linalg.norm(J*span,axis=1)
    # A truly insensitive equation is retained in the merit in its native unit.
    # This does not declare it solved; physical acceptance separately checks it.
    insensitive=row==0;row[insensitive]=1.
    column=np.linalg.norm(J*span/row[:,None],axis=0)
    if np.any(column==0):raise ValueError('Unobservable free coordinate: choose an explicit physical gauge')
    coordinate=span/column
    row=np.linalg.norm(J*coordinate,axis=1);row[insensitive]=1.
    return coordinate,row,dict(insensitive_rows=np.flatnonzero(insensitive).tolist(),
        basis='Source q spans equilibrated by native force Jacobian column sensitivity; force row units N or Nm, q units m or rad; frozen within a chunk')


def bounded_step(J,residual,q,bounds,coordinate_scale,residual_scale,radius):
    J=np.asarray(J,float);r=np.asarray(residual,float);q=np.asarray(q,float);bounds=np.asarray(bounds,float)
    c=np.asarray(coordinate_scale);s=np.asarray(residual_scale)
    if not 0<radius<=.1 or np.any(c<=0) or np.any(s<=0):raise ValueError('Invalid trust region/scales')
    # Preserve all XML bounds and the previous diagnostic's 0.03 native-unit cap.
    lo=np.maximum((bounds[:,0]-q)/c,-np.minimum(radius,.03/c))
    hi=np.minimum((bounds[:,1]-q)/c,np.minimum(radius,.03/c))
    A=J*c/s[:,None];b=r/s
    fit=lsq_linear(A,-b,bounds=(lo,hi),method='bvls',tol=1e-11,max_iter=150)
    if not fit.success:raise ValueError('Bounded force step failed: '+fit.message)
    step=c*fit.x;prediction=.5*(b@b-(b+A@fit.x)@(b+A@fit.x))
    singular=np.linalg.svd(A,compute_uv=False);threshold=np.finfo(float).eps*max(A.shape)*singular[0]
    rank=int(np.sum(singular>threshold));floor=np.linalg.norm(b-A@np.linalg.lstsq(A,b,rcond=None)[0])
    return step,dict(predicted_reduction=float(prediction),singular_values=singular.tolist(),rank=rank,
        linearized_unbounded_residual_floor=float(floor),scaled_gradient=(A.T@b).tolist(),
        source_or_box_active_columns=np.flatnonzero((abs(fit.x-lo)<1e-9)|(abs(fit.x-hi)<1e-9)).tolist(),
        maximum_coordinate_step=float(np.max(abs(step))),radius=radius)


def assess_trial(residual,trial,J,actual_step,scale):
    r=np.asarray(residual)/scale;t=np.asarray(trial)/scale
    predicted=r+(np.asarray(J)@actual_step)/scale
    prediction=.5*float(r@r-predicted@predicted);actual=.5*float(r@r-t@t)
    ratio=actual/prediction if prediction>0 else None
    return dict(predicted_reduction=prediction,actual_reduction=actual,reduction_ratio=ratio,
        accepted=bool(prediction>0 and actual>0 and ratio>=.1),
        model_error_norm=float(np.linalg.norm(t-predicted)))
