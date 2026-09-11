"""One smooth space warp from atlas space to the scaffold's ground, continuous across joints.

    W(x) = y + d(y),   y = G x   (G: the global binding similarity, atlas -> ground)
    d(y) = sum_i w_i phi(|y - c_i|) + A^T [y; 1],   phi(r) = -r,   P^T w = 0

phi = -r is the 3D biharmonic (thin-plate) kernel.  It is conditionally positive definite of
order 1, so on the constraint P^T w = 0 the quadratic form w^T K w is non-negative and is the
spline's bending energy.  With w = 0 and A = 0 the warp IS the global map, bit for bit: the
skin bundle and the whole-skin enclosure both reproduce the binding-map column through it.

The warp is one field over all of space.  It has no per-segment pieces, so it has no seams; it
can fold, and `jacobian` is what the no-folding gate reads.
"""
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist

def phi(r):
    return -r

class Warp:
    def __init__(self, base, centres, weights, affine, meta=None):
        self.base = np.asarray(base, float)
        self.centres = np.asarray(centres, float).reshape(-1, 3)
        self.weights = np.asarray(weights, float).reshape(-1, 3)
        self.affine = np.asarray(affine, float).reshape(4, 3)
        self.meta = dict(meta or {})
        if len(self.centres) != len(self.weights): raise ValueError('centres and weights differ in length')
        if np.linalg.det(self.base[:3, :3]) <= 0: raise ValueError('base map is not orientation-preserving')

    @classmethod
    def zero(cls, base):
        return cls(base, np.zeros((0, 3)), np.zeros((0, 3)), np.zeros((4, 3)), dict(kind='zero displacement'))

    @classmethod
    def load(cls, path):
        z = np.load(Path(path), allow_pickle=False)
        meta = {k[5:]: z[k].item() if z[k].shape == () else z[k].tolist() for k in z.files if k.startswith('meta_')}
        return cls(z['base'], z['centres'], z['weights'], z['affine'], meta)

    def save(self, path):
        extra = {'meta_' + k: np.asarray(v) for k, v in self.meta.items()}
        np.savez(Path(path), base=self.base, centres=self.centres, weights=self.weights, affine=self.affine, **extra)

    def ground(self, x):
        """The global map alone -- the same expression the skin bundle has always used."""
        G = self.base
        return np.asarray(x, float) @ G[:3, :3].T + G[:3, 3]

    def displacement(self, y, block=4096):
        y = np.asarray(y, float)
        out = y @ self.affine[:3] + self.affine[3]
        if len(self.centres):
            for s in range(0, len(y), block):
                out[s:s + block] += phi(cdist(y[s:s + block], self.centres)) @ self.weights
        return out

    def apply(self, x):
        y = self.ground(x)
        return y + self.displacement(y)

    def displacement_gradient(self, y, block=2048):
        """D[n, k, j] = d d_k / d y_j at each y.  phi'(r) = -1, so d phi / d y = -(y - c)/r."""
        y = np.asarray(y, float)
        D = np.broadcast_to(self.affine[:3].T, (len(y), 3, 3)).copy()
        for s in range(0, len(y), block):
            yb = y[s:s + block]
            diff = yb[:, None, :] - self.centres[None, :, :]          # (b, N, 3)
            r = np.linalg.norm(diff, axis=2)
            if np.any(r == 0): raise ValueError('gradient requested exactly at a spline centre')
            unit = diff / r[:, :, None]
            D[s:s + block] += -np.einsum('bnj,nk->bkj', unit, self.weights)
        return D

    def jacobian(self, x):
        """dW/dx: (I + Dd(y)) G_lin.  Its determinant is what gate 4 reads."""
        y = self.ground(x)
        return (np.eye(3)[None] + self.displacement_gradient(y)) @ self.base[:3, :3]

    def bending_energy(self):
        """sum over the three displacement components of w^T K w (>= 0 on P^T w = 0), units m."""
        if not len(self.centres): return 0.0
        K = phi(cdist(self.centres, self.centres))
        return float(np.einsum('ik,ij,jk->', self.weights, K, self.weights))

def _distances(Y, C, c2):
    """|Y - C| through BLAS, because scipy's cdist is single-threaded and the flow calls this
    2^N times per map."""
    d2 = (Y * Y).sum(1)[:, None] + c2[None, :] - 2.0 * (Y @ C.T)
    np.maximum(d2, 0.0, out=d2)
    return np.sqrt(d2, out=d2)

class FlowWarp:
    """W(x) = flow_1(G x): the time-1 flow of a STATIONARY velocity field, not a displacement.

        v(y) = sum_i w_i phi(|y - c_i|) + A^T [y; 1],   phi(r) = -r     (the same field family)
        psi_h(y) = y + h v(y),   h = 1 / 2^N,   flow_1 = psi_h squared N times

    A smooth velocity field's flow is a diffeomorphism, so det J > 0 holds BY CONSTRUCTION: the
    flow's Jacobian is the product of the per-step Jacobians I + h dv/dy, and every factor is
    positive-definite once h is small enough that the per-step displacement is small.  Folding is
    then a numerical question (how many squarings), not something a fit can produce.

    ON THE SQUARING.  Scaling and squaring is usually done on a GRID, where each squaring composes
    a stored displacement field with itself by interpolation.  This field is parametric and can be
    evaluated at any point, so squaring psi_h N times is exactly 2^N applications of psi_h, with no
    grid and no interpolation anywhere.  That is what makes the inverse exact rather than
    interpolation-limited: the inverse of the step y -> y + h v(y) is solved for its own fixed
    point to 1e-14, so forward composed with inverse returns a point to ~1e-13 m rather than to a
    grid's error.  With w = 0 and A = 0 the flow is the identity and the warp IS the global map,
    bit for bit, exactly as the spline's zero is.
    """
    def __init__(self, base, centres, weights, affine, steps, meta=None):
        self.base = np.asarray(base, float)
        self.centres = np.asarray(centres, float).reshape(-1, 3)
        self.weights = np.asarray(weights, float).reshape(-1, 3)
        self.affine = np.asarray(affine, float).reshape(4, 3)
        self.steps = int(steps)
        self.meta = dict(meta or {})
        if len(self.centres) != len(self.weights): raise ValueError('centres and weights differ in length')
        if np.linalg.det(self.base[:3, :3]) <= 0: raise ValueError('base map is not orientation-preserving')
        if self.steps < 1 or (self.steps & (self.steps - 1)): raise ValueError('steps must be 2^N, the squaring count')
        self._c2 = (self.centres ** 2).sum(1) if len(self.centres) else np.zeros(0)

    @classmethod
    def zero(cls, base, steps=1):
        return cls(base, np.zeros((0, 3)), np.zeros((0, 3)), np.zeros((4, 3)), steps, dict(kind='flow', zero='identity flow'))

    @classmethod
    def load(cls, path):
        z = np.load(Path(path), allow_pickle=False)
        meta = {k[5:]: z[k].item() if z[k].shape == () else z[k].tolist() for k in z.files if k.startswith('meta_')}
        return cls(z['base'], z['centres'], z['weights'], z['affine'], int(z['steps']), meta)

    def save(self, path):
        extra = {'meta_' + k: np.asarray(v) for k, v in self.meta.items()}
        np.savez(Path(path), kind=np.asarray('flow'), base=self.base, centres=self.centres,
                 weights=self.weights, affine=self.affine, steps=np.asarray(self.steps), **extra)

    def ground(self, x):
        G = self.base
        return np.asarray(x, float) @ G[:3, :3].T + G[:3, 3]

    def velocity(self, y, block=8192):
        y = np.asarray(y, float)
        out = y @ self.affine[:3] + self.affine[3]
        if len(self.centres):
            for s in range(0, len(y), block):
                r = _distances(y[s:s + block], self.centres, self._c2)
                np.negative(r, out=r)                       # phi(r) = -r
                out[s:s + block] += r @ self.weights
        return out

    def velocity_gradient(self, y, block=2048):
        """D[n, k, j] = dv_k / dy_j.  phi'(r) = -1, so d phi / d y = -(y - c) / r."""
        y = np.asarray(y, float)
        D = np.broadcast_to(self.affine[:3].T, (len(y), 3, 3)).copy()
        if not len(self.centres): return D
        for s in range(0, len(y), block):
            yb = y[s:s + block]
            diff = yb[:, None, :] - self.centres[None, :, :]
            r = np.linalg.norm(diff, axis=2)
            if np.any(r == 0): raise ValueError('velocity gradient requested exactly at a centre')
            diff /= r[:, :, None]
            D[s:s + block] += -np.einsum('bnj,nk->bkj', diff, self.weights)
        return D

    def flow(self, y, inverse=False, tol=1e-14, max_iter=64):
        """2^N applications of psi_h (or of its exact inverse, solved for its own fixed point)."""
        h = 1.0 / self.steps
        y = np.array(y, dtype=float, copy=True)
        if not inverse:
            for _ in range(self.steps): y += h * self.velocity(y)
            return y
        for _ in range(self.steps):
            p = y.copy()
            for _ in range(max_iter):
                q = y - h * self.velocity(p)
                delta = float(np.abs(q - p).max()); p = q
                if delta < tol: break
            else:
                raise RuntimeError('inverse step did not converge; the per-step displacement is too large')
            y = p
        return y

    def apply(self, x):
        return self.flow(self.ground(x))

    def unapply(self, y):
        """back to atlas space: the inverse flow, then the inverse of G."""
        return np.linalg.solve(self.base[:3, :3], (self.flow(y, inverse=True) - self.base[:3, 3]).T).T

    def jacobian_determinant(self, x, block=2048):
        """det dW/dx as the product of the per-step determinants times det G.

        Also returns the smallest single-step determinant seen, which is the numerical statement
        that the squaring is fine enough: every factor positive means the product cannot be
        negative, whatever the field does."""
        y = self.ground(x); h = 1.0 / self.steps
        det = np.full(len(y), float(np.linalg.det(self.base[:3, :3])))
        worst = np.inf
        eye = np.eye(3)[None]
        for _ in range(self.steps):
            step_det = np.linalg.det(eye + h * self.velocity_gradient(y, block=block))
            worst = min(worst, float(step_det.min()))
            det *= step_det
            y += h * self.velocity(y)
        return det, worst

    def maximum_step_displacement(self, y):
        return float(np.linalg.norm(self.velocity(y), axis=1).max() / self.steps)

    def bending_energy(self):
        """the VELOCITY field's w^T K w, the quantity the regularisation penalises."""
        if not len(self.centres): return 0.0
        K = phi(cdist(self.centres, self.centres))
        return float(np.einsum('ik,ij,jk->', self.weights, K, self.weights))

def load_warp(path):
    """A warp file is a spline (v1) or a flow (v2); the file says which."""
    z = np.load(Path(path), allow_pickle=False)
    kind = str(z['kind']) if 'kind' in z.files else 'spline'
    return FlowWarp.load(path) if kind == 'flow' else Warp.load(path)
