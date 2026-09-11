"""Frictionless sliding of a soft solid on a rigid bed, on DeformableRegion's energy.

The judge (docs/BODY_PARAMETERS.md, "The sliding base") holds every base node that penetrates the
bed ON it along the bed's normal and leaves it free tangentially, and forbids the other base nodes
from entering the bed. Each base node gets a local frame whose first axis is the bed's normal at
its closest bed point; in that frame the constraint is a BOUND on one coordinate: held nodes fix it
(lo = hi), unilateral nodes bound it below, tangential coordinates are free. The frames are
re-linearised at the new closest points as the nodes slide (SeatingDriver-side, see
scripts/seat_breast_sliding.py), until the held nodes' normal gap converges.

The solve is projected Newton rather than L-BFGS-B: at nu = 0.49 L-BFGS-B did not finish one
12,513-node breast in over 20 minutes. Energy and gradient are the parent's, unchanged; the Hessian
is the analytic one of the same compressible neo-Hookean law
    W = mu/2 (I1 - 3) - mu ln J + lambda/2 (ln J)^2,
    dP = mu dF + lambda (F^-T : dF) F^-T - (lambda ln J - mu) F^-T dF^T F^-T,
projected to positive semi-definite per element (eigenvalues of the 9x9 material tangent clamped
at zero) so every Newton step is a descent step. Bounds are handled by the active-set rule for box
constraints (a bound component is active if it sits on its bound with the gradient pushing into
it), with an Armijo line search over steps projected onto the box. Every accepted iterate has all
J > 0 because the line search rejects any step that inverts an element.
"""
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from .mechanics_backend import DeformableRegion


class SlidingRegion(DeformableRegion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        inv = self.inverse                                                    # [t, a-1, j]
        self.dN = np.concatenate((-inv.sum(axis=1, keepdims=True), inv), axis=1)   # (T,4,3) shape-function gradients
        dofs = (3 * self.tets[:, :, None] + np.arange(3)).reshape(-1, 12)     # element dof 3a+i
        self._rows = np.repeat(dofs, 12, axis=1).ravel(); self._cols = np.tile(dofs, (1, 12)).ravel()
        edges = np.linalg.norm(self.reference[self.tets[:, [1, 2, 3, 2, 3, 3]]] - self.reference[self.tets[:, [0, 0, 0, 1, 1, 2]]], axis=2)
        self.edge_m = float(edges.mean())

    def energy(self, y):
        if np.linalg.det(self.deformation(y)).min() <= 0: return np.inf
        return self.energy_gradient(y)[0]

    def hessian(self, y, *, project=True):
        f = self.deformation(y); A = np.linalg.inv(f); lnj = np.log(np.linalg.det(f))
        c = self.lam * lnj - self.mu
        g = np.swapaxes(A, 1, 2).reshape(-1, 9)                               # vec(F^-T), index 3i+j
        M = np.einsum('tki,tjl->tijlk', A, A).reshape(-1, 9, 9)               # vec(F^-T dF^T F^-T) = M vec(dF)
        H = self.mu[:, None, None] * np.eye(9) + self.lam[:, None, None] * g[:, :, None] * g[:, None, :] - c[:, None, None] * M
        if project:
            w, V = np.linalg.eigh(H)
            H = np.einsum('tab,tb,tcb->tac', V, np.maximum(w, 0.0), V)
        K = np.einsum('taj,tijln,tbn->taibl', self.dN, H.reshape(-1, 3, 3, 3, 3), self.dN, optimize=True)
        K *= self.volumes[:, None, None, None, None]
        n = 3 * len(self.reference)
        return sp.csr_matrix((K.reshape(-1), (self._rows, self._cols)), shape=(n, n))

    def solve_sliding(self, frames, lo, hi, *, start=None, rtol=1e-7, max_newton=300, log=None):
        """frames (N,3,3): columns are each node's local axes, u = frames @ v. lo, hi (N,3): bounds
        on v (lo == hi fixes a component). start (N,3): initial displacement u (default: current).
        Returns the displacement, iterations, the projected-force residual and its tolerance."""
        R = np.asarray(frames, float); lo = np.asarray(lo, float); hi = np.asarray(hi, float)
        N = len(self.reference)
        node = np.arange(N)[:, None, None]
        rows = np.broadcast_to(3 * node + np.arange(3)[None, :, None], (N, 3, 3))   # u index 3n+i
        cols = np.broadcast_to(3 * node + np.arange(3)[None, None, :], (N, 3, 3))   # v index 3n+k
        Q = sp.csr_matrix((R.reshape(-1), (rows.reshape(-1), cols.reshape(-1))), shape=(3 * N, 3 * N))   # u = Q v
        u0 = (self.positions - self.reference) if start is None else np.asarray(start, float)
        v = np.clip(np.einsum('nik,ni->nk', R, u0), lo, hi)
        j0 = np.linalg.det(self.deformation(self.reference + np.einsum('nik,nk->ni', R, v)))
        if j0.min() <= 0:
            # clipping the start onto changed bounds moved constrained nodes past their neighbours;
            # the caller's load step is too large. Say so, rather than failing inside the energy.
            raise ValueError(f"bounded start inverts {int((j0 <= 0).sum())} elements; reduce the load step")
        tol = rtol * float(np.mean(self.mu)) * self.edge_m ** 2
        fixed = lo == hi; began = time.perf_counter(); res = np.inf
        for it in range(max_newton):
            y = self.reference + np.einsum('nik,nk->ni', R, v)
            E, gu = self.energy_gradient(y)
            gv = np.einsum('nik,ni->nk', R, gu)
            active = fixed | ((v <= lo) & (gv > 0)) | ((v >= hi) & (gv < 0))
            res = float(np.abs(np.where(active, 0.0, gv)).max())
            if log: log(f"    newton {it:3d}  E {E:.6e}  projected force {res:.3e} (tol {tol:.1e})  active {int(active.sum())}")
            if res <= tol: break
            Kv = (Q.T @ self.hessian(y) @ Q).tocsr()
            free = np.flatnonzero(~active.ravel())
            Af = Kv[free][:, free]
            Af = Af + sp.identity(len(free), format='csr') * (1e-10 * float(Af.diagonal().mean()))
            d = np.zeros(3 * N); d[free] = spsolve(Af.tocsc(), -gv.ravel()[free])
            d = d.reshape(N, 3); alpha = 1.0
            for _ in range(40):
                vt = np.clip(v + alpha * d, lo, hi); Et = self.energy(self.reference + np.einsum('nik,nk->ni', R, vt))
                if Et <= E + 1e-4 * float(np.sum(gv * (vt - v))) + 1e-13 * abs(E): break
                alpha *= 0.5
            else:
                raise RuntimeError(f"line search failed at Newton iteration {it} (projected force {res:.3e})")
            v = vt
        u = np.einsum('nik,nk->ni', R, v)
        self.positions = self.reference + u
        return dict(displacement=u, local=v, iterations=it + 1, residual_n=res, tolerance_n=tol,
                    converged=res <= tol, wall_seconds=time.perf_counter() - began,
                    minimum_jacobian=float(np.linalg.det(self.deformation(self.positions)).min()))


def tangent_frames(n, hint=(0.0, 1.0, 0.0)):
    """(M,3,3) frames with columns [n, t1, t2]. The hint axis sets t1 = hint x n, so a bed that is
    a cylinder about the hint axis gets t2 exactly along it (which lets symmetry pins be bounds)."""
    n = n / np.linalg.norm(n, axis=1, keepdims=True)
    h = np.broadcast_to(np.asarray(hint, float), n.shape).copy()
    h[np.abs(np.einsum('ij,ij->i', n, h)) > 0.9] = (1.0, 0.0, 0.0)
    t1 = np.cross(h, n); t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    return np.stack([n, t1, np.cross(n, t1)], axis=2)


def seat_on_bed(region, base, closest, *, load_steps=8, pins=(), gap_tol_m=5e-5, max_outer=25,
                hint=(0.0, 1.0, 0.0), rtol=1e-7, jump_limit_m=0.01, assoc_tol_m=1e-4, move=None,
                rigid_m=None, log=None):
    """The sliding base. base: node indices on the surface facing the bed. closest(points) ->
    (c, n): closest bed points and the bed's outward unit normals there. A base node BEHIND the bed
    in the registered position is HELD: its signed normal gap is ramped to zero over load_steps and
    then held there, both ways, while it slides freely along the bed. Every other base node is
    UNILATERAL: its gap may not go negative. Each load step re-linearises -- new closest points,
    new frames, new bounds -- until the held gap error and any unilateral penetration are both
    <= gap_tol_m. pins: (node, global axis) held at zero, only where a local axis is parallel to it."""
    X = region.reference; N = len(X); base = np.asarray(base)
    c0, n0 = closest(X[base]); gap0 = np.einsum('ij,ij->i', n0, X[base] - c0)   # signed; < 0 behind the bed
    held = gap0 < 0
    u = np.zeros_like(X); record = []
    assoc = [c0.copy(), n0.copy()]          # persistent association, snapshotted per accepted step

    def associate(points):
        """Update each base node's bed patch, rejecting a teleport. The bed is three overlapping
        pectoralis parts: a node sliding tangentially can find a different sheet nearest, and against
        that sheet its gap flips sign -- tens of millimetres of apparent penetration through a
        constraint that forbids any. Real sliding moves the association ~0.02 mm per pass."""
        c_new, n_new = closest(points)
        # An association is rejected if it teleports, or if the bed function reports none (NaN):
        # either way the node keeps the patch it had.
        teleport = ~np.isfinite(c_new).all(1) | ~np.isfinite(n_new).all(1) | \
                   (np.linalg.norm(np.nan_to_num(c_new - assoc[0]), axis=1) > jump_limit_m)
        c_new[teleport] = assoc[0][teleport]; n_new[teleport] = assoc[1][teleport]
        moved = float(np.linalg.norm(c_new - assoc[0], axis=1).max())
        assoc[0], assoc[1] = c_new, n_new
        return c_new, n_new, int(teleport.sum()), moved

    # How far each held node is asked to travel along the bed normal. By default it closes its own
    # gap exactly; `move` supplies a different per-node distance -- the DECLARED modelling choice of
    # a smoothed depth field, where closing each gap exactly would fold the tissue through itself.
    travel = -gap0[held] if move is None else np.asarray(move, float)[held]
    # CONTROL R (docs/BODY_PARAMETERS.md, 4a58efd): drive the SAME held set through the SAME stepping
    # with a rigid translation of the whole base. A rigid motion preserves every Jacobian exactly, so
    # nothing can invert for any reason of physics or mesh quality, and a whole-body translation by
    # rigid_m satisfies every held constraint at zero strain energy. It separates "this configuration
    # is infeasible" from "the stepping is wrong".
    rigid = None if rigid_m is None else np.asarray(rigid_m, float)

    def advance(fraction, u_start):
        """Re-linearise and solve at this fraction of the held nodes' travel."""
        u_local = u_start
        for outer in range(max_outer):
            c, n, kept, _ = associate((X + u_local)[base])
            step = travel if rigid is None else n[held] @ rigid      # rigid: the normal part of one vector
            target = gap0[held] + fraction * step
            frames = np.tile(np.eye(3), (N, 1, 1)); frames[base] = tangent_frames(n, hint)
            lo = np.full(X.shape, -np.inf); hi = np.full(X.shape, np.inf)
            on_plane = np.einsum('ij,ij->i', n, c - X[base])                     # n.u that puts the node on the tangent plane
            lo[base[held], 0] = hi[base[held], 0] = on_plane[held] + target
            lo[base[~held], 0] = on_plane[~held]
            for node, axis in pins:
                col = np.flatnonzero(np.abs(frames[node][axis, :]) > 1 - 1e-9)
                if len(col) != 1: raise ValueError(f"pin on node {node}, axis {axis}, is not along one of its frame axes")
                lo[node, col[0]] = hi[node, col[0]] = 0.0
            r = region.solve_sliding(frames, lo, hi, start=u_local, rtol=rtol); u_local = r['displacement']
            c, n, kept, moved = associate((X + u_local)[base]); gap = np.einsum('ij,ij->i', n, (X + u_local)[base] - c)
            step = travel if rigid is None else n[held] @ rigid
            target = gap0[held] + fraction * step
            held_err = float(np.abs(gap[held] - target).max()) if held.any() else 0.0
            pen = float(max(0.0, -gap[~held].min())) if (~held).any() else 0.0
            if log: log(f"  fraction {fraction:.4f} pass {outer}: Newton {r['iterations']}, held gap error {held_err*1e3:.4f} mm, "
                        f"unilateral penetration {pen*1e3:.4f} mm, min J {r['minimum_jacobian']:.3f}, association moved {moved*1e3:.4f} mm", flush=True)
            settled = moved <= assoc_tol_m and pen <= gap_tol_m
            if settled and (fraction < 1.0 - 1e-12 or held_err <= gap_tol_m):
                return u_local, dict(passes=outer + 1, held_gap_error_m=held_err, unilateral_penetration_m=pen,
                                     newton_last=r['iterations'], min_J=r['minimum_jacobian'], converged=bool(r['converged'])), gap, c, n
        raise RuntimeError("re-linearisation did not converge")

    # Adaptive load stepping: closing tens of millimetres of penetration in equal steps moves held
    # nodes past their free neighbours and inverts surface tets before the first Newton iteration.
    # On failure the increment is halved and retried, then allowed to grow back.
    fraction, ds, cutbacks = 0.0, 1.0 / load_steps, 0
    while fraction < 1.0 - 1e-12:
        trial = min(1.0, fraction + ds)
        saved = [assoc[0].copy(), assoc[1].copy()]      # a failed step must not leave a stale association
        try:
            u_new, info, gap, c, n = advance(trial, u)
        except (ValueError, RuntimeError) as failure:
            assoc[0], assoc[1] = saved
            ds *= 0.5; cutbacks += 1
            if log: log(f"  cut back to {ds:.5f} at fraction {trial:.4f}: {failure}", flush=True)
            if ds < 1e-4: raise RuntimeError(f"load stepping stalled at fraction {fraction:.4f}: {failure}")
            continue
        u = u_new; fraction = trial; record.append(dict(fraction=fraction, **info))
        ds = min(ds * 1.5, 1.0 / load_steps)
    return dict(displacement=u, gap_m=gap, held=held, initial_gap_m=gap0, closest_m=c, normal=n, steps=record,
                cutbacks=cutbacks, minimum_jacobian=record[-1]['min_J'], converged=record[-1]['converged'])
