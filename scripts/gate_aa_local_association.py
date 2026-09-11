"""GATE AA: a LOCAL-SEARCH association for the anatomical bed.

Every rule tried so far -- per-node limit, all-or-nothing, derived threshold, adaptive stepping --
detects or tolerates a jump AFTER a global query has made it. One node of 3,123 jumps more than 1 mm
under an arbitrarily small motion, and no step size fixes a discontinuity. This constrains the QUERY
instead: the new association must lie within a bounded geodesic neighbourhood of the previous one,
grown over the bed's face adjacency, so a jump between sheets is impossible by construction.

The neighbourhood radius is the same measured pair the threshold used -- the node's own motion this
step plus the bed's facet scale -- so the method stays parameter-free: an association cannot need to
travel further than its node moved plus the surface's own resolution.

REPORTED, not hidden: a local search can hold the WRONG sheet when the true contact genuinely moves
to another one. The count of nodes whose local best is worse than the global closest point, and by
how much, is the price of the construction and is printed every time it runs.
"""
import sys, importlib.util, time
sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
CAP_S = float(sys.argv[sys.argv.index("--cap") + 1]) if "--cap" in sys.argv else 1500.0
RADIUS_SCALE = float(sys.argv[sys.argv.index("--radius-scale") + 1]) if "--radius-scale" in sys.argv else 1.0
import numpy as np
from scipy.spatial import cKDTree
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed
from ihm.assembly.prescribed_deformation import lame

ROOT = "/home/brandonin/Documents/IHM-1"
D = f"{ROOT}/data/derived/female-breast-sliding-v1/s1159/left"
spec = importlib.util.spec_from_file_location("seat", f"{ROOT}/scripts/seat_breast_sliding.py")
sys.argv = ["x"]; S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
from ihm.anatomy.normal_shooting import vertex_normals


class LocalAssociation:
    """Ray association restricted to faces reachable from the previous one within a bounded radius."""

    def __init__(self, V, F, dirs, facet_m, radius_scale=1.0):
        self.V, self.F, self.D = V, F, dirs / np.maximum(np.linalg.norm(dirs, axis=1, keepdims=True), 1e-30)
        # ONE normal field, defined exactly as bed_rays defines it: area-weighted face normals
        # accumulated to vertices. The first call used to take its normal from bed_rays and every
        # later call from vertex_normals(), so a node could be handed a different normal at ZERO
        # motion -- a second way for this object not to be a function of position.
        fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
        self.vn = np.zeros_like(V)
        for k in range(3): np.add.at(self.vn, F[:, k], fn)
        self.vn /= np.maximum(np.linalg.norm(self.vn, axis=1, keepdims=True), 1e-30)
        self.index = S.build_bed_index(V, F)
        # The neighbourhood is grown over FACE ADJACENCY, not by centroid distance: this bed's faces
        # run 1.6 mm median but 28.9 mm max, so a radius test on centroids excludes the very face the
        # association sits on. Same two-tier hazard as the ray index, in a new place.
        import collections
        edges = collections.defaultdict(list)
        for fi, (a, b, c) in enumerate(F):
            for e in ((a, b), (b, c), (c, a)): edges[(min(e), max(e))].append(fi)
        self.adj = collections.defaultdict(list)
        for fs in edges.values():
            if len(fs) == 2:
                self.adj[fs[0]].append(fs[1]); self.adj[fs[1]].append(fs[0])
        self.edge_m = float(np.median(np.linalg.norm(V[F[:, 1]] - V[F[:, 0]], axis=1)))
        self.facet_m, self.scale = facet_m, radius_scale
        self.face = None; self.point = None; self.normal = None; self.anchor = None
        self.worse_than_global = []          # the price of the construction, counted every update

    def _normal(self, c, f, d):
        """The bed normal at hit point c on face f: barycentric interpolation of self.vn, oriented
        anteriorly. Identical arithmetic to bed_rays, so the seeding call and every local update
        return the same normal for the same hit."""
        tri = self.V[self.F[f]]
        w = np.array([np.linalg.norm(np.cross(tri[1] - c, tri[2] - c)),
                      np.linalg.norm(np.cross(tri[2] - c, tri[0] - c)),
                      np.linalg.norm(np.cross(tri[0] - c, tri[1] - c))])
        w /= max(w.sum(), 1e-30)
        n = w @ self.vn[self.F[f]]
        n /= max(np.linalg.norm(n), 1e-30)
        return -n if n @ d > 0 else n

    def _hit(self, p, d, faces):
        v0 = self.V[self.F[faces][:, 0]]
        e1 = self.V[self.F[faces][:, 1]] - v0; e2 = self.V[self.F[faces][:, 2]] - v0
        best_t, best_f, best_sign = np.inf, -1, 1.0
        for sign in (1.0, -1.0):
            h = np.cross(sign * d, e2); a = np.einsum('fj,fj->f', e1, h)
            ok = np.abs(a) > 1e-16
            inv = 1.0 / np.where(ok, a, 1.0)
            s = p - v0
            u = inv * np.einsum('fj,fj->f', s, h)
            q = np.cross(s, e1)
            v = inv * (q @ (sign * d))
            t = inv * np.einsum('fj,fj->f', e2, q)
            ok &= (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9) & (t > 1e-12)
            if ok.any():
                j = int(np.argmin(np.where(ok, t, np.inf)))
                if t[j] < best_t: best_t, best_f, best_sign = t[j], faces[j], sign
        # the SIGNED distance: the ray is cast both ways and which way it hit is part of the answer.
        # Returning |t| and guessing the side from a vertex put associations on the wrong sheet and
        # inverted 43 elements in the first increment -- a bug of this script, not of local search.
        return best_sign * best_t, best_f

    def __call__(self, points):
        points = np.ascontiguousarray(np.asarray(points, float))
        if self.face is None:                                   # first call: the global query, once
            # THE SEED IS THE FACE THE RAY QUERY RETURNED. It used to be the face with the nearest
            # CENTROID, which on a bed running 1.6 mm median against 28.9 mm max is routinely not
            # the face the association sits on, so the neighbourhood grew from the wrong place and
            # the next call disagreed with this one at identical positions -- 141 of 3,123 nodes, up
            # to 6.5 mm. That disagreement, not the stepping, was gate AA's 42 inverted elements.
            t_back, f_back = S.ray_hits_indexed(points, self.D, self.V, self.F, self.index, S.RAY_REACH_M)
            t_front, f_front = S.ray_hits_indexed(points, -self.D, self.V, self.F, self.index, S.RAY_REACH_M)
            front = t_front < t_back
            t = np.where(front, t_front, t_back); f = np.where(front, f_front, f_back)
            hit = np.isfinite(t)
            c = points + np.where(front, -1.0, 1.0)[:, None] * np.where(hit, t, 0.0)[:, None] * self.D
            # KNOWN ANSWER for the duplication: this must reproduce bed_rays to the last bit, or the
            # face I just took is not the face that query used and the seed is wrong again.
            association, _ = S.bed_rays(self.V, self.F, self.D, S.RAY_REACH_M)
            cg, ng = association(points)
            ok = np.isfinite(cg).all(1)
            if not np.array_equal(ok, hit):
                raise AssertionError(f"seed reconstruction disagrees on {int((ok != hit).sum())} hits")
            dp = float(np.abs(c[ok] - cg[ok]).max()) if ok.any() else 0.0
            nl = np.stack([self._normal(c[i], int(f[i]), self.D[i]) for i in np.flatnonzero(ok)]) \
                 if ok.any() else np.zeros((0, 3))
            dn = float(np.abs(nl - ng[ok]).max()) if ok.any() else 0.0
            if dp > 1e-12 or dn > 1e-9:
                raise AssertionError(f"seed reconstruction differs from bed_rays: point {dp:.3e} m, "
                                     f"normal {dn:.3e}")
            print(f"    seed: face taken from the ray query, reproducing bed_rays to "
                  f"{dp:.2e} m and {dn:.2e} in the normal, over {int(ok.sum())} of {len(points)} hits",
                  flush=True)
            self.point = np.where(hit[:, None], c, points)
            # the seed's normal comes from _normal, not from bed_rays' own copy of the same
            # arithmetic, so the seeding and updating paths are bit-identical rather than merely
            # equal to 3e-14 -- there is no reason to leave a difference that can be removed.
            self.normal = (-self.D).copy()
            self.normal[np.flatnonzero(ok)] = nl
            self.face = np.where(hit, f, -1)
            self.anchor = points.copy()
            return self.point.copy(), self.normal.copy()
        motion = np.linalg.norm(points - self.anchor, axis=1)
        radius = self.scale * (motion + self.facet_m)
        worse = 0; excess = []
        for i in range(len(points)):
            if radius[i] <= 0: continue                         # discrimination arm: never updates
            if self.face[i] < 0: continue                       # no bed at the seed: hold
            rings = int(np.ceil(radius[i] / max(self.edge_m, 1e-6))) + 1
            seen = {int(self.face[i])}; frontier = [int(self.face[i])]
            for _ in range(min(rings, 8)):
                nxt = [g for f0 in frontier for g in self.adj[f0] if g not in seen]
                seen.update(nxt); frontier = nxt
                if not frontier: break
            faces = np.fromiter(seen, int, len(seen))
            t, f = self._hit(points[i], self.D[i], faces)
            if f < 0: continue                                  # nothing in the neighbourhood: hold
            c = points[i] + t * self.D[i]
            if not np.isfinite(c).all(): continue
            self.point[i] = c; self.face[i] = f
            self.normal[i] = self._normal(c, int(f), self.D[i])
        self.anchor = points.copy()
        return self.point.copy(), self.normal.copy()


def global_comparison(local, points):
    """How often is the local best worse than the global closest point, and by how much?"""
    association, _ = S.bed_rays(local.V, local.F, local.D, S.RAY_REACH_M)
    g, _ = association(points)
    ok = np.isfinite(g).all(1)
    dl = np.linalg.norm(local.point[ok] - points[ok], axis=1)
    dg = np.linalg.norm(g[ok] - points[ok], axis=1)
    worse = dl > dg + 1e-9
    return int(worse.sum()), int(ok.sum()), (1000 * (dl - dg)[worse] if worse.any() else np.zeros(0))


P = np.load(f"{D}/prepared.npz")
X, T, base = P["X"], P["T"], P["base"]
move = np.load(f"{D}/smoothed_move.npy")
mu, lam = lame(1000.0, 0.49)
region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
local = LocalAssociation(P["bedV"], P["bedF"], P["ray_directions"], S.GAP_TOL_M, RADIUS_SCALE)


def idempotence_gate(local, points, label):
    """A query meant to be a function of position must return the same answer called twice at the
    same position. This is NOT a known answer: a known answer tests the value, idempotence tests
    whether the instrument is a function at all. It costs one extra call, and it is the only thing
    that catches a non-idempotent query consumed inside a cut-back loop -- which is what made gate
    AA's inversion count invariant to the step size and my reasoning from that count worthless."""
    c1, n1 = local(points)
    g1 = local.face.copy()
    c2, n2 = local(points)
    dc = float(np.abs(c2 - c1).max()); dn = float(np.abs(n2 - n1).max())
    moved = int((local.face != g1).sum())
    far = int((np.linalg.norm(c2 - c1, axis=1) > 1e-3).sum())
    ok = moved == 0 and dc <= 1e-12 and dn <= 1e-9
    print(f"  IDEMPOTENCE at {label}: faces changed {moved} of {len(points)}, point moved at most "
          f"{1000*dc:.6f} mm ({far} over 1 mm), normal at most {dn:.2e} -- "
          f"{'PASS' if ok else 'FAIL'}", flush=True)
    return ok


print("GATE AA idempotence, run BEFORE any AA number is quoted:")
p0 = X[base].copy()
gate = idempotence_gate(local, p0, "the initial positions")
gate &= idempotence_gate(local, p0 + 5e-4 * local.D, "0.5 mm along each ray")
gate &= idempotence_gate(local, p0, "the initial positions again")
if not gate:
    print("\nGATE AA: NOT RUN -- the association is not a function of position, so no number it "
          "produces is evidence about local search.")
    raise SystemExit(1)

start = time.time()
class Cap(RuntimeError): pass
def logger(m, flush=True):
    print(f"    [{time.time()-start:6.0f}s] {m}", flush=True)
    if time.time() - start > CAP_S: raise Cap(f"wall-clock cap {CAP_S:.0f}s")

print(f"GATE AA: local-search association, radius = {RADIUS_SCALE} x (node motion + "
      f"{1000*S.GAP_TOL_M:.0f} mm facet), on the anatomical bed's seating drive")
r, err = None, None
stalled = {}


def on_stall(fraction, advance, u):
    """At the stall: measure the price where the body has actually MOVED, then ask whether a step of
    size ZERO also fails. The price is taken first, because advance() re-associates and would
    overwrite the state being measured."""
    stalled["fraction"] = fraction
    pts = (X + u)[base]
    stalled["price"] = global_comparison(local, pts)
    print(f"\n  at the stall, fraction {fraction:.4f}, the body having moved "
          f"{1000 * np.abs(u).max():.2f} mm at most:", flush=True)
    try:
        advance(fraction, u)
        stalled["zero"] = "a ZERO-sized step SUCCEEDS -- the failure needs a finite increment"
    except Exception as e:
        stalled["zero"] = f"a ZERO-sized step ALSO FAILS -- {type(e).__name__}: {e}"
    print(f"  ZERO-INCREMENT PROBE: {stalled['zero']}", flush=True)


try:
    r = seat_on_bed(region, base, local, load_steps=8, gap_tol_m=S.GAP_TOL_M, jump_limit_m=1e9,
                    assoc_tol_m=S.ASSOC_TOL_M, move=move, freeze_frames=True, association="persistent",
                    lost_bed="hold", facet_m=S.GAP_TOL_M, on_stall=on_stall, log=logger)
except Cap as e: err = str(e)
except Exception as e: err = f"{type(e).__name__}: {e}"
elapsed = time.time() - start
# The price is only measurable where the run actually moved the body: comparing at the initial
# positions returns "0 worse" by construction, since the first association IS the global one.
measurable = r is not None or "price" in stalled
if measurable:
    worse, total, excess = (global_comparison(local, (X + r["displacement"])[base])
                            if r is not None else stalled["price"])
if r is not None:
    print(f"\nGATE AA: COMPLETED to 1.0 in {elapsed:.0f}s, {len(r['steps'])} increments, "
          f"cutbacks {r['cutbacks']}, min J {r['steps'][-1]['min_J']:.4f}")
else:
    print(f"\nGATE AA: did not complete in {elapsed:.0f}s -- {err}")
print(f"  against gate V on the same drive: 0.02% of the drive in 705 s")
if measurable:
    where = ("the completed drive" if r is not None
             else f"the stalled state at fraction {stalled['fraction']:.4f}")
    print(f"  THE PRICE, measured at {where}: local best worse than the global closest point for "
          f"{worse} of {total} nodes"
          + (f", by a median {np.median(excess):.3f} mm and at most {excess.max():.3f} mm" if len(excess) else ""))
    if "zero" in stalled: print(f"  ZERO-INCREMENT PROBE: {stalled['zero']}")
else:
    print("  THE PRICE: not measurable -- the run did not complete, and at the starting positions the")
    print("  local association IS the global one, so a comparison there would read 0 by construction.")
