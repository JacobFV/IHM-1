"""Quasistatic prescribed-displacement solves on DeformableRegion, for seating soft tissue on a
rigid bed.

DeformableRegion (mechanics_backend.py) minimises compressible neo-Hookean energy
    W = mu/2 (I1 - 3) - mu ln J + lambda/2 (ln J)^2
over constant-strain tetrahedra with L-BFGS-B, and its boundary conditions are fixed nodes and
one planar indenter. A breast base conformed to a curved chest wall needs arbitrary PRESCRIBED
displacements per degree of freedom. L-BFGS-B already takes per-DOF bounds, so a prescribed DOF
is lo = hi = its value; everything else is the parent's energy and gradient, unchanged. The
displacement is applied in load steps, each warm-started from the last, so the first iterate
never carries the full base jump into elements that would invert.

This law is also FEBio's "neo-Hookean" material, so an FEBio run on the same mesh and boundary
conditions is an independent implementation of the same discretisation.
"""
import time
import numpy as np
from scipy.optimize import minimize
from .mechanics_backend import DeformableRegion, tetra_box   # noqa: F401  (tetra_box re-exported for callers)


def lame(young_pa, nu):
    return young_pa / (2 * (1 + nu)), young_pa * nu / ((1 + nu) * (1 - 2 * nu))


def tet_volumes(positions, tets):
    x = np.asarray(positions, float)
    return np.linalg.det(np.swapaxes(x[tets[:, 1:]] - x[tets[:, 0, None]], 1, 2)) / 6


class PrescribedRegion(DeformableRegion):
    def solve_prescribed(self, mask, displacement, *, load_steps=8, maxiter=20000, gtol=1e-10):
        """mask (N,3) bool: which DOFs are prescribed; displacement (N,3): their final values (m).
        Returns positions, displacement, minimum J, volume ratio, free-DOF residual relative to the
        prescribed reaction, and iterations. Raises if any tet leaves J > 0.2."""
        mask = np.asarray(mask, bool); target = np.asarray(displacement, float)
        if mask.shape != self.reference.shape or target.shape != self.reference.shape:
            raise ValueError("mask and displacement must be (N,3)")
        scale = float(np.min(np.ptp(self.reference, axis=0)))
        energy_scale = float(np.mean(self.mu)) * scale ** 3
        u = np.zeros_like(self.reference); iterations = 0; began = time.perf_counter()
        for k in range(1, load_steps + 1):
            value = target * (k / load_steps)
            lo = np.full_like(self.reference, -np.inf); hi = np.full_like(self.reference, np.inf)
            lo[mask] = value[mask]; hi[mask] = value[mask]
            u0 = u.copy(); u0[mask] = value[mask]
            def objective(z):
                y = self.reference + (z * scale).reshape(-1, 3)
                e, g = self.energy_gradient(y, trial=True)
                return e / energy_scale, g.ravel() * scale / energy_scale
            r = minimize(objective, u0.ravel() / scale, jac=True, method="L-BFGS-B",
                         bounds=list(zip(lo.ravel() / scale, hi.ravel() / scale)),
                         options={"ftol": 1e-16, "gtol": gtol, "maxiter": maxiter, "maxls": 50, "maxcor": 30})
            u = r.x.reshape(-1, 3) * scale; iterations += int(r.nit)
        y = self.reference + u
        j = np.linalg.det(self.deformation(y))
        if j.min() <= 0.2: raise RuntimeError(f"left the constitutive domain: min J {j.min():.3f} <= 0.2")
        _, g = self.energy_gradient(y)
        free, fixed = g[~mask], g[mask]
        reaction = float(np.max(np.abs(fixed))) if fixed.size else 1.0
        self.positions = y.copy()
        return dict(positions=y, displacement=u, minimum_jacobian=float(j.min()),
                    volume_ratio=float(tet_volumes(y, self.tets).sum() / tet_volumes(self.reference, self.tets).sum()),
                    residual_relative=float(np.max(np.abs(free)) / reaction) if free.size else 0.0,
                    iterations=iterations, message=str(r.message), wall_seconds=time.perf_counter() - began)
