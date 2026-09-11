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

    def solve_sliding(self, frames, lo, hi, *, start=None, rtol=1e-7, max_newton=300, project=True, log=None):
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
        v_prev = np.einsum('nik,ni->nk', R, u0)
        v = np.clip(v_prev, lo, hi)
        delta = v - v_prev                       # what the changed bounds ask of the constrained DOFs
        clip0 = np.abs(delta).copy()             # kept: delta is reassigned inside the active-set loop
        # The bound change enters through the STIFFNESS, not by clipping alone. Clipping moves the
        # constrained coordinates and leaves every other one at the previous solution, which is
        # infeasible by construction: that start could not carry even a RIGID TRANSLATION, which
        # preserves every Jacobian exactly (Control R, docs/BODY_PARAMETERS.md). One linear elastic
        # response at the previous state carries the free coordinates along with the constrained
        # ones -- the standard incremental prescribed-displacement step, using the Hessian and the
        # active set this solver already has.
        if np.any(delta):
            y_prev = self.reference + np.einsum('nik,nk->ni', R, v_prev)
            if np.linalg.det(self.deformation(y_prev)).min() > 0:
                Kv = (Q.T @ self.hessian(y_prev) @ Q).tocsr()
                # ACTIVE SET at the start, not a clip after it. The response alone reproduces a rigid
                # translation exactly (0.0000 mm, zero inversions); clipping it back afterwards threw
                # 79 one-sided nodes up to 7.34 mm off and inverted 185 elements, because a bound that
                # the response violates has to be carried BY the solve -- with its neighbours moving
                # too -- not imposed on the answer afterwards.
                active = (lo == hi) | (np.abs(delta) > 0)
                for _ in range(8):
                    v = np.where(active, np.clip(v_prev + delta, lo, hi), v_prev)
                    flat = active.ravel()
                    free = np.flatnonzero(~flat); cols = np.flatnonzero(flat)
                    if not len(free) or not len(cols): break
                    Aff = Kv[free][:, free]
                    Aff = Aff + sp.identity(len(free), format='csr') * (1e-10 * float(Aff.diagonal().mean()))
                    step = (v - v_prev).ravel()
                    response = np.zeros(3 * N)
                    response[free] = spsolve(Aff.tocsc(), -(Kv[free][:, cols] @ step[cols]))
                    trial = v_prev + step.reshape(N, 3) + response.reshape(N, 3)
                    violated = (trial < lo - 1e-15) | (trial > hi + 1e-15)
                    if not violated.any():
                        v = trial; break
                    active = active | violated                       # carry the new bounds in the next solve
                    delta = np.where(active, np.clip(trial, lo, hi) - v_prev, 0.0)
                else:
                    v = np.clip(trial, lo, hi)
        j0 = np.linalg.det(self.deformation(self.reference + np.einsum('nik,nk->ni', R, v)))
        if j0.min() <= 0:
            # State what was observed. The previous wording named a cause -- "reduce the load step" --
            # that was false, and it cost this line five suspects' worth of investigation.
            # The clip magnitude is part of the observation, not commentary: a start that inverts with
            # delta at zero cannot have been made by the bounds, and one that inverts with delta
            # finite at a ZERO load increment was made by the re-association and not by the step.
            # TWO numbers, because they are two different things and one sentence covering both was
            # wrong: by this line v has been overwritten by the active-set loop, so |v - v_prev| is
            # the WHOLE start displacement including the elastic response carried on free nodes --
            # it counted 5,660 nodes when only 4,097 are bounded at all, which is what exposed it.
            # clip0 is the bound-driven part alone.
            tot = np.abs(v - v_prev)
            raise ValueError(f"the start of this increment has {int((j0 <= 0).sum())} inverted "
                             f"elements; the BOUNDS moved {int((clip0.max(1) > 1e-3).sum())} nodes by "
                             f"more than 1 mm, at most {1000 * float(clip0.max()):.4f} mm, and the "
                             f"start including the elastic response moved at most "
                             f"{1000 * float(tot.max()):.4f} mm")
        tol = rtol * float(np.mean(self.mu)) * self.edge_m ** 2
        fixed = lo == hi; began = time.perf_counter(); res = np.inf
        for it in range(max_newton):
            y = self.reference + np.einsum('nik,nk->ni', R, v)
            E, gu = self.energy_gradient(y)
            gv = np.einsum('nik,ni->nk', R, gu)
            active = fixed | ((v <= lo) & (gv > 0)) | ((v >= hi) & (gv < 0))
            res = float(np.abs(np.where(active, 0.0, gv)).max())
            if log: log(f"    newton {it:3d}  E {E:.6e}  projected force {res:.3e} (tol {tol:.1e})  "
                        f"active {int(active.sum())}  min J {float(np.linalg.det(self.deformation(y)).min()):.4f}")
            if res <= tol: break
            Kv = (Q.T @ self.hessian(y, project=project) @ Q).tocsr()
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
                    bound_displacement_m=float(clip0.max()),
                    bound_nodes_over_1mm=int((clip0.max(1) > 1e-3).sum()),
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
                rigid_m=None, freeze_frames=False, project=True, solver_log=None, stop_fraction=1.0,
                association='persistent', lost_bed='hold', facet_m=1e-3, drive_full=False,
                prescribe=None, phases=(1.0,), on_stall=None, on_step=None, log=None):
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
    assoc = [c0.copy(), n0.copy(), np.zeros(len(base), bool), X[base].copy()]   # ..., positions when it was last taken
    lost_history = []; refusals = []          # persistent association, snapshotted per accepted step

    def associate(points):
        """Update each base node's bed patch, rejecting a teleport. The bed is three overlapping
        pectoralis parts: a node sliding tangentially can find a different sheet nearest, and against
        that sheet its gap flips sign -- tens of millimetres of apparent penetration through a
        constraint that forbids any. Real sliding moves the association ~0.02 mm per pass."""
        # association: 'persistent' rejects a teleport or an invalid update and keeps the previous one
        # (the original behaviour); 'none' never updates after the first; 'all' takes every valid update.
        # MIXING the two -- freezing 46.6% of constraints while updating the rest -- is what made the
        # constraint set inconsistent and deformed the body (gate S).
        # lost_bed names what happens to a node whose ray no longer meets the bed within the reach:
        # 'hold' keeps its last valid association, 'release' drops its constraint for that step. Declared,
        # not chosen silently, and counted per association.
        if association in ('none', 'adaptive') and lost_history:
            # 'adaptive' never updates DURING a step: the update is attempted after the step, and a
            # step whose update would exceed any node's bar is rejected and shrunk (gate V). The
            # quantity driving the step size is then association motion -- the thing that actually
            # breaks -- rather than min J.
            return assoc[0], assoc[1], 0, 0.0
        c_new, n_new = closest(points)
        if association == 'allornothing':
            # THRESHOLD, derived rather than chosen: the step's own node motion plus the bed's facet
            # scale. A closest point on a locally smooth bed cannot outrun the node that owns it by
            # more than the surface's own resolution; anything further has changed FEATURE, which is
            # what a teleport is. An update is taken WHOLE or refused WHOLE, because mixing updated
            # and stale constraints is what deformed the body (gate S).
            # PER-NODE bar. Under a rigid drive every node moves the same distance and a single
            # scalar sufficed; under a non-rigid one it does not, and a global bar is then too tight
            # for the fast nodes or too loose for the slow ones. The refuse-or-take DECISION stays
            # global -- that is what keeps every constraint in one epoch -- and only the bar is
            # per-node. With equal motion this reduces exactly to the previous rule, which is what
            # the T, U and W regressions check.
            bar = np.linalg.norm(points - assoc[3], axis=1) + facet_m
            move = np.linalg.norm(np.nan_to_num(c_new - assoc[0]), axis=1)
            invalid = ~np.isfinite(c_new).all(1) | ~np.isfinite(n_new).all(1)
            over = int((move > bar).sum()) + int(invalid.sum())
            lost_history.append(int(invalid.sum()))
            refusals.append(bool(over))
            if over:                                   # refuse the whole update; every constraint stays in one epoch
                return assoc[0], assoc[1], over, 0.0
            assoc[0], assoc[1], assoc[3] = c_new, n_new, points.copy()
            return c_new, n_new, 0, float(move.max())
        invalid = ~np.isfinite(c_new).all(1) | ~np.isfinite(n_new).all(1)
        rejected = invalid if association == 'all' else (
            invalid | (np.linalg.norm(np.nan_to_num(c_new - assoc[0]), axis=1) > jump_limit_m))
        c_new[rejected] = assoc[0][rejected]; n_new[rejected] = assoc[1][rejected]
        moved = float(np.linalg.norm(c_new - assoc[0], axis=1).max())
        assoc[0], assoc[1] = c_new, n_new
        assoc[2] = invalid if lost_bed == 'release' else np.zeros(len(invalid), bool)
        lost_history.append(int(invalid.sum()))
        teleport = rejected
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

    def advance(fraction, u_start, theta, aim):
        """Re-linearise and solve, closing `theta` of the distance remaining to fraction `aim`.

        THE REPAIR (gate BB). The bound used to be set ABSOLUTELY to `on_plane + gap0 +
        fraction*travel`, where `on_plane = n.(c - X[base])` is the displacement that puts the node
        on its association plane. The load fraction scales `travel` and does NOT scale `on_plane`,
        so once the association had drifted, the head of the next increment demanded that drift
        instantly -- 7.19 mm over 26 nodes at a step of size ZERO, four inverted elements, identical
        when repeated. Shrinking the step could not touch it, which is what made the inversion count
        invariant across a 256x range and what every gate from R to AA was unknowingly fighting.

        Now the bound is interpolated from where the node IS to where the constraint wants it, by
        theta -- the fraction of the REMAINING drive this increment consumes. Every term, the
        association term included, is under the increment parameter. At theta = 0 the bound is the
        node's own current position, so a step of size zero moves nothing, exactly. At theta = 1 the
        bound is the constraint itself, so the endpoint is unchanged: this reschedules the approach,
        it does not relax where the drive lands.

        With freeze_frames, the association is taken ONCE at the head of the step and the constraint
        directions are held fixed while the solve runs, so a drift inside a step cannot be the
        association's doing; re-association happens only between steps (gate R'', 04f32a0).
        """
        u_local = u_start
        passes = 1 if freeze_frames else max_outer
        for outer in range(passes):
            c, n, kept, _ = associate((X + u_local)[base])
            step = travel if rigid is None else n[held] @ rigid      # rigid: the normal part of one vector
            frames = np.tile(np.eye(3), (N, 1, 1)); frames[base] = tangent_frames(n, hint)
            lo = np.full(X.shape, -np.inf); hi = np.full(X.shape, np.inf)
            on_plane = np.einsum('ij,ij->i', n, c - X[base])                     # n.u that puts the node on the tangent plane
            if drive_full:
                # GATE W: the base is driven as a whole, all three components, so a TANGENTIAL motion
                # is actually imposed. With a normal-only constraint a tangential translation is driven
                # by nothing -- u = 0 already satisfies n.u = 0 -- and the association would never be
                # asked to follow the bed, which is the regime every earlier control missed.
                full = np.einsum('nik,ni->nk', frames[base], np.broadcast_to(fraction * rigid, (len(base), 3)))
                lo[base] = hi[base] = full
            else:
                # v0 is the node's CURRENT normal displacement, in the coordinate the bound is written
                # in. It uses the SAME einsum over the SAME frames that solve_sliding uses to form
                # v_prev, not the mathematically equal n . u: at theta = 0 the bound must land
                # bit-identically on the node's own coordinate, or gate BB reads float dust instead
                # of the exact zero it is asking for.
                v0 = np.einsum('nik,ni->nk', frames[base], u_local[base])[:, 0]
                want = on_plane[held] + gap0[held] + aim * step          # fully seated at the aim
                lo[base[held], 0] = hi[base[held], 0] = v0[held] + theta * (want - v0[held])
                # The unilateral bound carries the association term too, so it gets the same
                # treatment -- but only where the drifted association has left the node in
                # violation. Where the node already satisfies it, the true constraint stands;
                # relaxing it there would invent a push that the contact does not ask for.
                free_side = on_plane[~held]
                lo[base[~held], 0] = np.minimum(free_side, v0[~held] + theta * (free_side - v0[~held]))
            if assoc[2].any():                          # released: no constraint for this step
                lo[base[assoc[2]]] = -np.inf; hi[base[assoc[2]]] = np.inf
            if prescribe is not None:
                # GATE X: an arbitrary per-node prescribed displacement, scaled by the load fraction.
                # A homogeneous field is exact only if the WHOLE boundary is compatible with it --
                # driving the base alone would let the free surface relax and there would be no
                # closed-form answer to compare against.
                idx, values = prescribe
                lo[idx] = hi[idx] = np.einsum('nik,ni->nk', frames[idx], fraction * values)
            for node, axis in pins:
                col = np.flatnonzero(np.abs(frames[node][axis, :]) > 1 - 1e-9)
                if len(col) != 1: raise ValueError(f"pin on node {node}, axis {axis}, is not along one of its frame axes")
                lo[node, col[0]] = hi[node, col[0]] = 0.0
            j_entry = float(np.linalg.det(region.deformation(X + u_local)).min())
            r = region.solve_sliding(frames, lo, hi, start=u_local, rtol=rtol, project=project,
                                     log=solver_log); u_local = r['displacement']
            c, n, kept, moved = associate((X + u_local)[base]); gap = np.einsum('ij,ij->i', n, (X + u_local)[base] - c)
            step = travel if rigid is None else n[held] @ rigid
            # Measured against the AIM, not against a nominal mid-drive schedule the repair no longer
            # tracks: the increment now closes theta of what remains, so "distance from
            # gap0 + fraction*travel" would be a number the drive is not trying to hit. At the aim,
            # where the acceptance test below actually reads it, the two coincide.
            to_aim = np.abs(gap[held] - (gap0[held] + aim * step)) if held.any() else np.zeros(1)
            held_err = float(to_aim.max())
            # THE MEDIAN TOO. Reporting only the max made 'distance closed' unreadable: the drive
            # ran to fraction 0.3115 while the max sat at 13-15 mm against 14.18 mm at the start,
            # so the load fraction was measuring load APPLIED and not distance CLOSED, and there
            # was no way to tell from the log whether one node or the whole sheet was stuck.
            held_mid = float(np.median(to_aim))
            pen = float(max(0.0, -gap[~held].min())) if (~held).any() else 0.0
            if log: log(f"  fraction {fraction:.4f} pass {outer}: Newton {r['iterations']}{'' if r['converged'] else ' UNCONVERGED'}, held gap to the aim {held_err*1e3:.4f} mm max / {held_mid*1e3:.4f} median, "
                        f"unilateral penetration {pen*1e3:.4f} mm, min J {r['minimum_jacobian']:.3f}, "
                        + (f"association moved {moved*1e3:.4f} mm" if association != 'adaptive'
                           else "association held for the step (adaptive updates after it)"), flush=True)
            settled = moved <= assoc_tol_m and pen <= gap_tol_m
            if freeze_frames or (settled and (abs(fraction - aim) > 1e-12 or held_err <= gap_tol_m)):
                return u_local, dict(passes=outer + 1, held_gap_error_m=held_err, unilateral_penetration_m=pen,
                                     bound_displacement_m=r['bound_displacement_m'],
                                     bound_nodes_over_1mm=r['bound_nodes_over_1mm'],
                                     newton_last=r['iterations'], min_J=r['minimum_jacobian'],
                                     min_J_entry=j_entry, association_moved_m=moved,
                                     converged=bool(r['converged'])), gap, c, n
        raise RuntimeError("re-linearisation did not converge")

    # Adaptive load stepping: closing tens of millimetres of penetration in equal steps moves held
    # nodes past their free neighbours and inverts surface tets before the first Newton iteration.
    # On failure the increment is halved and retried, then allowed to grow back.
    # phases: the schedule of target fractions. (1.0,) drives forward; (1.0, 0.0) drives forward and
    # back for the reversibility gate, in ONE call so the association history stays continuous.
    fraction, ds, cutbacks = 0.0, 1.0 / load_steps, 0
    targets = [min(t, stop_fraction) for t in phases]
    target = targets.pop(0)
    while True:
        if abs(fraction - target) <= 1e-12:
            if not targets: break
            target = targets.pop(0); ds = 1.0 / load_steps
            if abs(fraction - target) <= 1e-12: continue
        direction = 1.0 if target > fraction else -1.0
        trial = fraction + direction * min(ds, abs(target - fraction))
        # theta: the share of the REMAINING drive this increment consumes. It is what now carries the
        # association term, so a cut-back shrinks EVERY term of the prescribed displacement.
        theta = (trial - fraction) / (target - fraction)
        saved = [assoc[0].copy(), assoc[1].copy()]      # a failed step must not leave a stale association
        try:
            u_new, info, gap, c, n = advance(trial, u, theta, target)
        except (ValueError, RuntimeError) as failure:
            assoc[0], assoc[1] = saved
            ds *= 0.5; cutbacks += 1
            if log: log(f"  cut back to {ds:.5f} at fraction {trial:.4f}: {failure}", flush=True)
            if ds < 1e-4:
                # A ZERO-SIZED STEP, offered to the caller before the stall is raised. advance() at
                # the fraction already accepted asks for no new travel at all, so anything it does is
                # the RE-ASSOCIATION acting on an already-deformed state rather than the increment.
                # It is the only way to separate the two, since every cut-back keeps re-associating.
                if on_stall is not None: on_stall(fraction, advance, u, target)
                raise RuntimeError(f"load stepping stalled at fraction {fraction:.4f}: {failure}")
            continue
        if association == 'adaptive':
            points = (X + u_new)[base]
            c_new, n_new = closest(points)
            bar = np.linalg.norm(points - assoc[3], axis=1) + facet_m
            invalid = ~np.isfinite(c_new).all(1) | ~np.isfinite(n_new).all(1)
            over = int(((np.linalg.norm(np.nan_to_num(c_new - assoc[0]), axis=1) > bar) | invalid).sum())
            refusals.append(bool(over))
            if over:                                  # the STEP is rejected, not the update
                assoc[0], assoc[1] = saved
                ds *= 0.5; cutbacks += 1
                if log: log(f"  step rejected at fraction {trial:.4f}: {over} associations would exceed "
                            f"their own bar; shrinking to {ds:.5f}", flush=True)
                if ds < 1e-6: raise RuntimeError(f"association-driven stepping stalled at {fraction:.4f}")
                continue
            moved_now = float(np.linalg.norm(c_new - assoc[0], axis=1).max())
            assoc[0], assoc[1], assoc[3] = c_new, n_new, points.copy()
            lost_history.append(int(invalid.sum()))
            info = dict(info, association_moved_m=moved_now)
        u = u_new; fraction = trial; record.append(dict(fraction=fraction, **info))
        # GATE BB is asked here, after every ACCEPTED step, not only at a stall: if the repair
        # works there may be no stall at all, and a gate that only fires on failure would never
        # run on the code it is meant to certify.
        if on_step is not None: on_step(fraction, advance, u, target)
        ds = min(ds * 1.5, 1.0 / load_steps)
    return dict(displacement=u, gap_m=gap, held=held, initial_gap_m=gap0, closest_m=c, normal=n, steps=record,
                cutbacks=cutbacks, lost_bed_per_association=lost_history, association=association,
                refusals=refusals, refusal_rate=(float(np.mean(refusals)) if refusals else None),
                lost_bed_rule=lost_bed, minimum_jacobian=record[-1]['min_J'], converged=record[-1]['converged'])
