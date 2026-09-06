"""Training-only standardization/SVD and stable discrete autonomous predictors."""
from dataclasses import dataclass
import numpy as np


@dataclass
class ReducedPredictor:
    center: np.ndarray
    scale: np.ndarray
    basis: np.ndarray
    transition: np.ndarray
    dt_s: float
    diagnostics: dict

    def forecast(self, initial, steps):
        initial=np.asarray(initial,float)
        if initial.shape!=self.center.shape or not np.isfinite(initial).all() or not isinstance(steps,int) or steps<0:
            raise ValueError('invalid forecast state or step count')
        z=((initial-self.center)/self.scale)@self.basis
        out=np.empty((steps,len(self.center)))
        for i in range(steps):
            z=z@self.transition
            out[i]=self.center+(z@self.basis.T)*self.scale
        return out

    def to_dict(self):
        return dict(schema_version=1,kind='stable_reduced_autonomous_linear',dt_s=self.dt_s,
                    center=self.center.tolist(),scale=self.scale.tolist(),basis=self.basis.tolist(),
                    transition=self.transition.tolist(),diagnostics=self.diagnostics)

    @classmethod
    def from_dict(cls,d):
        if d.get('schema_version')!=1 or d.get('kind')!='stable_reduced_autonomous_linear':
            raise ValueError('unsupported predictor schema')
        c,sc,b,a=(np.asarray(d[k],float) for k in ('center','scale','basis','transition'))
        if c.ndim!=1 or sc.shape!=c.shape or b.ndim!=2 or b.shape[0]!=len(c) or a.shape!=(b.shape[1],b.shape[1]) or not np.all(sc>0) or not all(np.isfinite(v).all() for v in (c,sc,b,a)) or not np.isfinite(d['dt_s']) or d['dt_s']<=0:
            raise ValueError('invalid predictor arrays or cadence')
        if np.linalg.svd(a,compute_uv=False)[0]>1+1e-9:
            raise ValueError('predictor must be contractive')
        return cls(c,sc,b,a,float(d['dt_s']),d['diagnostics'])


def fit_predictor(values, dt_s, *, rank=8, train_fraction=.75, ridge=1e-5, max_gain=.9999):
    """Fit reduced VAR(1), clipping singular values to bound every-step gain.

    Holdout is an uninterrupted forecast initialized at the last training sample.
    Rank, center, scale and all coefficients use training samples only.
    No mechanism or intervention response is identified by this empirical fit.
    Relative constant detection preserves unit scaling while variance remains
    representable in float64; underflowed variance cannot resolve dynamics.
    """
    y=np.asarray(values,float)
    if y.ndim!=2 or len(y)<12 or y.shape[1]<1 or not np.isfinite(y).all():
        raise ValueError('finite samples-by-variables matrix with at least 12 rows required')
    if not np.isfinite(dt_s) or dt_s<=0 or not .1<=train_fraction<=.9 or not isinstance(rank,int) or rank<1 or not np.isfinite(ridge) or ridge<0 or not 0<max_gain<1:
        raise ValueError('invalid predictor hyperparameters')
    n=int(len(y)*train_fraction)
    if n<4: raise ValueError('insufficient training samples')
    train=y[:n]; center=train.mean(0); scale=train.std(0)
    constant=scale<=np.max(abs(train),axis=0)*1e-12
    scale=np.where(constant,1.,scale)
    x=(train-center)/scale
    _,sing,vt=np.linalg.svd(x,full_matrices=False)
    tol=sing[0]*max(x.shape)*np.finfo(float).eps if sing[0]>0 else 0
    data_rank=int(np.sum(sing>tol)); r=min(rank,max(1,data_rank))
    basis=vt[:r].T; z=x@basis
    singular=np.linalg.svd(z[:-1],compute_uv=False)
    a=np.linalg.lstsq(np.vstack([z[:-1],np.sqrt(ridge)*np.eye(r)]),np.vstack([z[1:],np.zeros((r,r))]),rcond=None)[0]
    u,g,v=np.linalg.svd(a); stable=u@np.diag(np.minimum(g,max_gain))@v
    eig=np.linalg.eigvals(stable)
    model=ReducedPredictor(center,scale,basis,stable,float(dt_s),{})
    hold=model.forecast(train[-1],len(y)-n)
    train_hat=center+((z[:-1]@stable)@basis.T)*scale
    rmse=lambda a,b: np.sqrt(np.mean((a-b)**2,axis=0)).tolist()
    model.diagnostics=dict(train_samples=n,holdout_samples=len(y)-n,train_fraction=train_fraction,
        holdout_protocol='chronological uninterrupted forecast; initial state last training sample',
        data_rank=data_rank,retained_rank=r,singular_values=sing.tolist(),
        reduced_design_condition=float(singular[0]/singular[-1]) if singular[-1]>tol else None,
        ridge=ridge,max_gain=max_gain,unconstrained_max_singular_value=float(g[0]),
        spectral_radius=float(max(abs(eig))),eigenvalues_real=eig.real.tolist(),eigenvalues_imag=eig.imag.tolist(),
        train_one_step_rmse=rmse(train[1:],train_hat),holdout_rmse=rmse(y[n:],hold),
        holdout_persistence_rmse=rmse(y[n:],np.broadcast_to(train[-1],y[n:].shape)),
        holdout_normalized_rmse=(np.sqrt(np.mean((y[n:]-hold)**2,axis=0))/scale).tolist(),
        constant_training_variables=np.flatnonzero(constant).tolist())
    return model
