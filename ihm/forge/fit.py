"""fit process pressure, not an unrelated latent predictor."""
from dataclasses import dataclass
import numpy as np
from ihm.processes import Process


@dataclass(frozen=True)
class AffineFit:
    inputs: tuple[str, ...]
    output: str
    weights: tuple[float, ...]
    bias: float
    residual_variance: float
    source: str
    training_subjects: tuple[str, ...]

    def process(self, ident, topology, noise=0.):
        # Derivative regression residual variance is not a continuous diffusion
        # intensity. Caller supplies the latter explicitly in state^2 / second.
        return Process(ident, self.inputs, self.output, topology, self.weights,
                       self.bias, noise,
                       f'fitted:{self.source};subjects={",".join(self.training_subjects)}')


def fit_affine(x, derivative, inputs, output, source, training_subjects):
    x = np.asarray(x, float); y = np.asarray(derivative, float)
    if x.ndim != 2 or x.shape[1] != len(inputs) or y.shape != (len(x),):
        raise ValueError('regression shape mismatch')
    if len(x) <= len(inputs)+1 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('finite overdetermined training data required')
    if not source or not training_subjects or not all(training_subjects):
        raise ValueError('training source and subject provenance required')
    design = np.column_stack((x, np.ones(len(x))))
    weights, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank != design.shape[1]:
        raise ValueError('process parameters are not identifiable from this design')
    residual = y - design @ weights
    return AffineFit(tuple(inputs), output, tuple(weights[:-1]), float(weights[-1]),
                     float(residual @ residual / (len(x)-rank)), source, tuple(training_subjects))


def evaluate(fit, x, derivative, subjects):
    if not subjects or set(subjects) & set(fit.training_subjects):
        raise ValueError('evaluation must use held-out subjects')
    x = np.asarray(x, float); y = np.asarray(derivative, float)
    if x.ndim != 2 or x.shape[1] != len(fit.inputs) or y.shape != (len(x),) or not len(x):
        raise ValueError('evaluation shape mismatch')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('evaluation values must be finite')
    residual = y - (x @ fit.weights + fit.bias)
    return {'rmse': float(np.sqrt(np.mean(residual**2))), 'n': len(y),
            'subjects': list(subjects), 'source': fit.source}
