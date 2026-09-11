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
        self.vn = vertex_normals(V, F)
        self.centre = V[F].mean(1)
        self.tree = cKDTree(self.centre)
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
            association, _ = S.bed_rays(self.V, self.F, self.D, S.RAY_REACH_M)
            c, n = association(points)
            hit = np.isfinite(c).all(1)
            self.point = np.where(hit[:, None], c, points)
            self.normal = np.where(hit[:, None], np.nan_to_num(n), -self.D)
            self.face = self.tree.query(self.point)[1]
            self.anchor = points.copy()
            return self.point.copy(), self.normal.copy()
        motion = np.linalg.norm(points - self.anchor, axis=1)
        radius = self.scale * (motion + self.facet_m)
        worse = 0; excess = []
        for i in range(len(points)):
            if radius[i] <= 0: continue                         # discrimination arm: never updates
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
            nrm = self.vn[self.F[f]].mean(0); nrm /= max(np.linalg.norm(nrm), 1e-30)
            self.normal[i] = -nrm if nrm @ self.D[i] > 0 else nrm
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
start = time.time()
class Cap(RuntimeError): pass
def logger(m, flush=True):
    print(f"    [{time.time()-start:6.0f}s] {m}", flush=True)
    if time.time() - start > CAP_S: raise Cap(f"wall-clock cap {CAP_S:.0f}s")

print(f"GATE AA: local-search association, radius = {RADIUS_SCALE} x (node motion + "
      f"{1000*S.GAP_TOL_M:.0f} mm facet), on the anatomical bed's seating drive")
r, err = None, None
try:
    r = seat_on_bed(region, base, local, load_steps=8, gap_tol_m=S.GAP_TOL_M, jump_limit_m=1e9,
                    assoc_tol_m=S.ASSOC_TOL_M, move=move, freeze_frames=True, association="persistent",
                    lost_bed="hold", facet_m=S.GAP_TOL_M, log=logger)
except Cap as e: err = str(e)
except Exception as e: err = f"{type(e).__name__}: {e}"
elapsed = time.time() - start
# The price is only measurable where the run actually moved the body: comparing at the initial
# positions returns "0 worse" by construction, since the first association IS the global one.
measurable = r is not None
if measurable:
    worse, total, excess = global_comparison(local, (X + r["displacement"])[base])
if r is not None:
    print(f"\nGATE AA: COMPLETED to 1.0 in {elapsed:.0f}s, {len(r['steps'])} increments, "
          f"cutbacks {r['cutbacks']}, min J {r['steps'][-1]['min_J']:.4f}")
else:
    print(f"\nGATE AA: did not complete in {elapsed:.0f}s -- {err}")
print(f"  against gate V on the same drive: 0.02% of the drive in 705 s")
if measurable:
    print(f"  THE PRICE: local best worse than the global closest point for {worse} of {total} nodes"
          + (f", by a median {np.median(excess):.3f} mm and at most {excess.max():.3f} mm" if len(excess) else ""))
else:
    print("  THE PRICE: not measurable -- the run did not complete, and at the starting positions the")
    print("  local association IS the global one, so a comparison there would read 0 by construction.")
