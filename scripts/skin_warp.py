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
