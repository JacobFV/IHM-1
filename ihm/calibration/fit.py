"""Dimensionless residual/Jacobian diagnostics for bounded least squares."""
import numpy as np
from scipy.optimize import least_squares


def bounded_fit(predict, observed, sigma, initial, bounds, *, parameter_scale=None):
    """Fit residuals (prediction-observation)/sigma in scaled parameter coordinates.

    sigma is an explicitly supplied observation error scale, never guessed from
    population or repeated-scan SEM. Covariance assumes those errors known and
    independent; rank deficiency or active bounds suppress its usual local form.
    """
    y=np.asarray(observed,float); sigma=np.asarray(sigma,float)
    p=np.asarray(initial,float); lower,upper=(np.asarray(v,float) for v in bounds)
    scale=np.asarray(parameter_scale if parameter_scale is not None else np.maximum(abs(p),1),float)
    if y.ndim!=1 or len(y)<1 or sigma.shape!=y.shape or p.ndim!=1 or len(p)<1 or lower.shape!=p.shape or upper.shape!=p.shape or scale.shape!=p.shape:
        raise ValueError('incompatible data, bounds or parameter shapes')
    if not all(np.isfinite(v).all() for v in (y,sigma,p,lower,upper,scale)) or np.any(sigma<=0) or np.any(scale<=0) or np.any(lower>=upper) or np.any(p<lower) or np.any(p>upper):
        raise ValueError('finite bounded parameters and positive uncertainty/scales required')
    def residual(q):
        prediction=np.asarray(predict(q*scale),float)
        if prediction.shape!=y.shape or not np.isfinite(prediction).all():
            raise ValueError('predictor must return matching finite observation vector')
        return (prediction-y)/sigma
    result=least_squares(residual,p/scale,bounds=(lower/scale,upper/scale),jac='3-point',xtol=1e-12,gtol=1e-12,ftol=1e-12,max_nfev=5000)
    if not result.success: raise RuntimeError('bounded fitting failed: '+result.message)
    _,sing,vh=np.linalg.svd(result.jac,full_matrices=True)
    threshold=(sing[0] if len(sing) else 0)*max(result.jac.shape)*1e-8
    rank=int(np.sum(sing>threshold)); full=rank==len(p)
    active=np.flatnonzero(result.active_mask).tolist()
    cov=None
    if full and not active:
        # Inverse in dimensionless coordinates, transformed back to parameters.
        inverse=(vh.T[:,:len(sing)]/sing**2)@vh[:len(sing)]
        cov=(scale[:,None]*inverse*scale[None,:]).tolist()
    return dict(parameters=(result.x*scale).tolist(),parameter_scale=scale.tolist(),
        lower_bounds=lower.tolist(),upper_bounds=upper.tolist(),rank=rank,parameters_count=len(p),identifiable=full,
        singular_values=sing.tolist(),rank_relative_tolerance=1e-8,
        condition_number=float(sing[0]/sing[-1]) if full else None,
        null_directions_scaled=vh[rank:].tolist(),active_bounds=active,covariance=cov,
        covariance_assumption='local linear Gaussian; supplied independent known observation sigma; absent at deficient rank or active bounds',
        residuals_standardized=result.fun.tolist(),weighted_sum_squares=float(result.fun@result.fun),
        residual_degrees_of_freedom=len(y)-rank,evaluations=result.nfev,success=bool(result.success))
