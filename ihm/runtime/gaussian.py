"""exact transitions for affine stochastic differential equations."""
import numpy as np
from scipy.linalg import expm


def transition(a, b, q, dt):
    n = len(b)
    # Doubling avoids overflow in Van Loan's -A block for long stable intervals.
    scale = max(1, int(np.ceil(np.linalg.norm(a, ord=np.inf)*dt / 8)))
    squarings = int(np.ceil(np.log2(scale)))
    h = dt / 2**squarings
    block = np.zeros((2*n, 2*n))
    block[:n, :n] = a; block[:n, n:] = q; block[n:, n:] = -a.T
    e = expm(block*h); f = e[:n, :n]; noise = e[:n, n:] @ f.T
    affine = np.zeros((n+1, n+1)); affine[:n, :n] = a; affine[:n, n] = b
    offset = expm(affine*h)[:n, n]
    for _ in range(squarings):
        noise = f @ noise @ f.T + noise
        offset = f @ offset + offset
        f = f @ f
    return f, offset, (noise + noise.T)/2


def covariance(value, n, positive=False):
    value = np.asarray(value, float)
    if n < 1 or value.shape != (n, n) or not np.isfinite(value).all():
        raise ValueError('invalid covariance shape or values')
    diagonal = np.diag(value)
    if np.any(diagonal < 0) or (positive and np.any(diagonal == 0)):
        raise ValueError('negative variance or singular measurement covariance')
    # Check correlation scale: an absolute tolerance is invalid when voltages,
    # ion concentrations and cell counts coexist in one covariance matrix.
    scale = np.sqrt(np.where(diagonal > 0, diagonal, 1.))
    normalized = value / scale[:, None] / scale[None, :]
    zero = diagonal == 0
    if np.any(value[zero, :] != 0) or np.any(value[:, zero] != 0):
        raise ValueError('zero variance state cannot have nonzero covariance')
    if not np.allclose(normalized, normalized.T, rtol=1e-12, atol=1e-14):
        raise ValueError('covariance must be symmetric')
    eig = np.linalg.eigvalsh((normalized+normalized.T)/2)
    tolerance = 64*np.finfo(float).eps*n*max(1., float(eig.max()))
    if eig.min() < -tolerance or (positive and eig.min() <= 0):
        raise ValueError('covariance must be positive definite' if positive else 'covariance must be positive semidefinite')
    return value
