"""A skin carrier continuous across joints: one smooth space warp, fitted on the bones.

The gates are pre-registered in docs/SEGMENT_CONTACT_SURFACES.md (commit 93d3d57, "A skin
carrier continuous across joints").  This script does two things, in two stages, so the warp is
fixed before three of its four gates are computed:

  --stage fit    known-answer controls, correspondences, the regularisation rule, the warp, gate 1
  --stage score  gate 4 and the reported quantities from the warp; gates 2 and 3 read from the
                 instruments that own them (measure_skin_enclosure_whole.py --warp, and the
                 skin_minus_bone_minimum_y_m of a bundle built by build_skin_contact_meshes.py --warp)

THE WARP.  W(x) = y + d(y), y = G x, with G the global binding similarity (exactly the map the
0.888 column is measured through) and d a regularised 3D thin-plate spline (scripts/skin_warp.py).

CORRESPONDENCES (fixed before any fit).  For every segment, area-weighted samples a on the atlas
bone group (200 per segment; torso 600, pelvis 300).  Source: G a.  Target: the nearest point ON
the scaffold's own bone mesh (exact point-to-triangle, at binding.json's reference pose) to the
per-segment similarity's image M_seg a (registration.json, fit_segment_registration.py).  The
farthest 10% of each segment's correspondences by |target - M_seg a| are dropped -- the same 10%
trim the per-segment fit itself uses, for atlas parts the scaffold mesh does not carry.

THE REGULARISATION RULE (fixed before gates 2-4 were computed; it reads the bone
correspondences and nothing else).  Smoothing lambda on the fixed grid {0} U logspace(-8, 0, 17)
(metres, the kernel's unit).  5-fold cross-validation over the correspondences (folds from
default_rng(7)); held-out error = RMS |W(source) - target| on the held-out fold.  lambda is the
LARGEST grid value whose cross-validated RMS is within 1% of the minimum -- ties go to the
smoother warp.  The final warp is fitted on all correspondences at that lambda.  Nothing about
the skin, the enclosure, the heel or folding enters the choice, and the rule is not revisited
after they are seen.
"""
import argparse, gzip, hashlib, importlib.util, json, sys, time
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import trimesh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from skin_warp import Warp, FlowWarp, load_warp, phi

FAMILY = "spline"          # rebound in __main__; "flow" is the v2 instrument (3b1526c)
# The flow's own numerical settings.  The squaring count is not free: it is whatever makes the
# per-step displacement <= STEP_TARGET_M, which is the numerical statement of "cannot fold".
STEP_TARGET_M, OUTER, OUTER_TOL = 1e-3, 6, 1e-4
OUT = ROOT / "data/derived/skin-warp-v1"
REGISTRATION = ROOT / "data/derived/anatomy-segment-registration/registration.json"
BUNDLE = "data/derived/segment-contact-meshes/skin-warp-v1"
REFERENCE_BUNDLE = ROOT / "data/derived/segment-contact-meshes/skin-binding"

# ---- fixed before the fit
N_CORR, N_CORR_TORSO, N_CORR_PELVIS = 200, 600, 300
# There is no correspondence trim.  A 10% trim was imported into the v1 pre-registration by analogy
# with the per-segment fit's ICP outlier rejection, never measured, and the sweep in b034206 found
# its cost strictly monotone: 0% 1.666 mm, 2% 1.729, 5% 1.847, 10% 2.138 to-surface recovery at a
# fixed fitted count.  It is deleted rather than set to zero so there is nothing here to retune.
# CAVEAT carried wherever this is quoted: the control's truth is a known field on this same body, so
# its largest separations are hard but GENUINE; on a real subject-to-scaffold registration they may
# be WRONG correspondences, where a trim could be protective, and no ground truth exists there to
# tell.  The only evidence available says remove it; none says keep it.
CONTROL_FITTED_FRACTION = 0.90   # count-matching for the historical control arms ONLY, not a trim
FOLDS, LAMBDAS, TIE = 5, [0.0] + [float(x) for x in np.logspace(-8, 0, 17)], 0.01
# ---- the gates, as pre-registered in 93d3d57
GATE1_MARGIN_M, ENCLOSURE_MEAN, ENCLOSURE_FOOT, HEEL_RANGE_M = 1e-3, 0.95, 0.95, (-0.025, -0.005)
FOOT = ("calcn_l", "calcn_r", "toes_l", "toes_r")
# ---- fit_segment_registration.py's own sampling, replayed so gate 1 is scored on its samples
N_FIT, N_FIT_TORSO, N_EVAL, FIT_TRIM = 8000, 20000, 8000, 0.10

def say(*a): print(*a, flush=True)
def load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def apply(M, p): return p @ M[:3, :3].T + M[:3, 3]
def sim_matrix(s, R, t):
    M = np.eye(4); M[:3, :3] = s * R; M[:3, 3] = t; return M

B = load("btfe", "scripts/build_tissue_force_elements.py")
bscm = load("bscm", "scripts/build_skin_contact_meshes.py")
bind = B.load_module("bind_anatomy", ROOT / "scripts/bind_anatomy_to_segments.py")
render = B.load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")

# ------------------------------------------------ geometry, as fit_segment_registration.py reads it
entities = json.loads(B.ANATOMY.read_text())["entities"]
gpath = {e["id"]: ROOT / e["reference_geometry"]["path"] for e in entities if e.get("reference_geometry")}
groups = {}
for e in entities:
    if e["role"] == "rigid_bone": groups.setdefault(bind.named_segment(e["name"]), []).append(e["id"])
segments = sorted(groups)

def atlas_mesh(seg):
    V, F, off = [], [], 0
    for i in groups[seg]:
        g = json.loads(gzip.decompress(gpath[i].read_bytes()))
        v = np.asarray(g["positions"], float).reshape(-1, 3); f = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        V.append(v); F.append(f + off); off += len(v)
    return np.concatenate(V), np.concatenate(F)

from ihm.spatial.vtk import surface as read_surface
GEOM = ROOT / "data/raw/anatomy/opensim-models/source/Geometry"
def osim_meshes():
    root = ET.parse(B.MODEL).getroot().find("Model"); out = {}
    for body in root.iter("Body"):
        V, F, off = [], [], 0
        for mesh in body.iter("Mesh"):
            fac = np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            p, f = read_surface(GEOM / mesh.findtext("mesh_file")); V.append(p * fac); F.append(f + off); off += len(p)
        if V: out[body.get("name")] = (np.concatenate(V), np.concatenate(F))
    return out

def sample(V, F, n, rng):
    tri = V[F]; a = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    k = rng.choice(len(F), n, p=a / a.sum()); r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = tri[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2]

def sample_faces(V, F, n, rng):
    tri = V[F]; a = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    k = rng.choice(len(F), n, p=a / a.sum()); r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = tri[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2], k

def residual_points(S, tgt):
    """fit_segment_registration.residual, on already-mapped points (same arithmetic)."""
    d1 = cKDTree(tgt).query(S)[0]; d2 = cKDTree(S).query(tgt)[0]; d = np.concatenate([d1, d2])
    return float(np.sqrt((d ** 2).mean())), float(np.sqrt((np.sort(d)[: int(len(d) * (1 - FIT_TRIM))] ** 2).mean()))

def closest_on_triangles(P, A, B, C):
    """Closest point on triangle ABC to P, row by row (Ericson, Real-Time Collision Detection 5.1.5).
    Written out because trimesh.triangles.closest_point returned points up to 1e-6 m off a face for
    points lying on it (C4)."""
    ab, ac, ap = B - A, C - A, P - A; bp, cp = P - B, P - C
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1); d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1); d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    out = np.empty_like(P); done = np.zeros(len(P), bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        def put(m, X):
            nonlocal done
            m = m & ~done; out[m] = X[m]; done = done | m
        put((d1 <= 0) & (d2 <= 0), A)
        put((d3 >= 0) & (d4 <= d3), B)
        put((vc <= 0) & (d1 >= 0) & (d3 <= 0), A + (d1 / (d1 - d3))[:, None] * ab)
        put((d6 >= 0) & (d5 <= d6), C)
        put((vb <= 0) & (d2 >= 0) & (d6 <= 0), A + (d2 / (d2 - d6))[:, None] * ac)
        put((va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0), B + ((d4 - d3) / ((d4 - d3) + (d5 - d6)))[:, None] * (C - B))
        s = va + vb + vc
        put(np.ones(len(P), bool), A + (vb / s)[:, None] * ab + (vc / s)[:, None] * ac)
    return out

def projector(V, F, rng, n=200000, return_face=False):
    """Exact nearest point on a triangle mesh.  The nearest of a dense surface sample is a point ON
    the surface at distance u, so the answer is within u; every triangle whose bounding sphere comes
    within u of the query is tested exactly.  (A k-nearest-sample candidate heuristic missed an 8 mm2
    face in C4; this search cannot, by construction.)"""
    pts, _ = sample_faces(V, F, n, rng); tree = cKDTree(pts); tri = V[F]
    cen = tri.mean(1); rad = np.linalg.norm(tri - cen[:, None], axis=2).max(1); ctree = cKDTree(cen); rmax = rad.max()
    def project(Q):
        u = tree.query(Q)[0]
        lists = ctree.query_ball_point(Q, u + rmax + 1e-12)
        qi = np.concatenate([np.full(len(l), i) for i, l in enumerate(lists)]); ti = np.concatenate([np.asarray(l, int) for l in lists])
        keep = np.linalg.norm(Q[qi] - cen[ti], axis=1) <= u[qi] + rad[ti] + 1e-12; qi, ti = qi[keep], ti[keep]
        Cc = closest_on_triangles(Q[qi], tri[ti, 0], tri[ti, 1], tri[ti, 2]); dc = np.linalg.norm(Cc - Q[qi], axis=1)
        order = np.lexsort((dc, qi)); first = order[np.unique(qi[order], return_index=True)[1]]
        if len(first) != len(Q): raise RuntimeError('a query found no candidate triangle')
        if return_face: return Cc[first], dc[first], u, ti[first]
        return Cc[first], dc[first], u
    return project

# ------------------------------------------------ the spline
def tps_system(S):
    K = phi(cdist(S, S)); P = np.hstack([S, np.ones((len(S), 1))])
    Q, R = np.linalg.qr(P, mode="complete"); Q1, Q2 = Q[:, :4], Q[:, 4:]
    ev, V = np.linalg.eigh(Q2.T @ K @ Q2)
    return K, Q1, Q2, R[:4], ev, V

def tps_solve(system, D, lam):
    """(K + lam I) w + P a = D,  P^T w = 0;  w = Q2 gamma."""
    K, Q1, Q2, R, ev, V = system
    w = Q2 @ (V @ ((V.T @ (Q2.T @ D)) / (ev + lam)[:, None]))
    a = np.linalg.solve(R, Q1.T @ (D - K @ w - lam * w))
    return w, a

def tps_predict(Y, S, w, a):
    return phi(cdist(Y, S)) @ w + Y @ a[:3] + a[3]

def local_preconditioner(S, lam, k=48, chunk=4096):
    """Beatson-Cherrie-Mouat local cardinal functions: column i of an approximate inverse from the
    local spline problem on i's k nearest centres.  Plain CG does not converge on a 1 mm-spaced
    centre set (20,000 iterations and still 1e-3); this is the standard fix for polyharmonic
    splines and it changes only the SPEED of the solve, never the system being solved."""
    from scipy.sparse import csc_matrix
    N = len(S); k = min(k, N); idx = cKDTree(S).query(S, k=k)[1]
    rows, cols, vals = [], [], []
    for s in range(0, N, chunk):
        J = idx[s:s + chunk]; b = len(J); Sl = S[J]                       # (b, k, 3)
        Kl = phi(np.linalg.norm(Sl[:, :, None, :] - Sl[:, None, :, :], axis=3)) + lam * np.eye(k)[None]
        Pl = np.concatenate([Sl, np.ones((b, k, 1))], axis=2)             # (b, k, 4)
        A = np.zeros((b, k + 4, k + 4))
        A[:, :k, :k] = Kl; A[:, :k, k:] = Pl; A[:, k:, :k] = np.transpose(Pl, (0, 2, 1))
        rhs = np.zeros((b, k + 4, 1)); rhs[:, 0, 0] = 1.0                  # nearest neighbour is self
        z = np.linalg.solve(A, rhs)[:, :k, 0]
        rows.append(J.ravel()); cols.append(np.repeat(np.arange(s, s + b), k)); vals.append(z.ravel())
    M = csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return (M + M.T) * 0.5

def tps_solve_matrix_free(S, D, lam, tol=1e-12, maxiter=20000, block=2048, report=None, precondition=True, k=48,
                          log_every=25, budget_s=7200.0):
    """The SAME system as tps_solve, solved without ever forming K.

    With anchors there are 38,821 centres, so K is 12 GB and the nullspace QR is another 12 GB --
    more than this machine has free.  The equations are unchanged: (K + lam I) w + P a = D with
    P^T w = 0.  Projecting onto the complement of span(P) (where phi = -r is positive definite)
    turns the first equation into a symmetric positive-definite system in w, which conjugate
    gradients solves with only streamed kernel products.  Validated against tps_solve on the
    bone-only system, where it must return the same weights.
    """
    from skin_warp import _distances
    N = len(S); P = np.hstack([S, np.ones((N, 1))]); PtP = np.linalg.inv(P.T @ P); c2 = (S * S).sum(1)
    def proj(X): return X - P @ (PtP @ (P.T @ X))
    def Kmul(X):
        out = np.empty_like(X)
        for s in range(0, N, block):
            r = _distances(S[s:s + block], S, c2); np.negative(r, out=r)
            out[s:s + block] = r @ X
        return out
    def A(X): return proj(Kmul(X) + lam * X)
    M = local_preconditioner(S, lam, k=k) if precondition else None
    def apply_M(X): return proj(M @ X) if M is not None else X
    b = proj(D); w = np.zeros_like(D); r = b - A(w); z = apply_M(r); p = z.copy(); rz = (r * z).sum(0)
    scale = np.maximum(np.sqrt((b * b).sum(0)), 1e-300); used = 0
    start = time.time(); stopped = None
    for it in range(1, maxiter + 1):
        Ap = A(p); alpha = rz / np.maximum((p * Ap).sum(0), 1e-300)
        w += alpha * p; r -= alpha * Ap
        used = it; rel = float(np.max(np.sqrt((r * r).sum(0)) / scale)); elapsed = time.time() - start
        # A solve that prints nothing for hours cannot be told from a wedged one.
        if log_every and (it % log_every == 0 or it == 1):
            say(f"    CG {it:6d}  relative residual {rel:.3e}  {elapsed:7.0f}s")
        if rel < tol: break
        if budget_s is not None and elapsed > budget_s:
            stopped = "wall-clock budget"
            say(f"    CG stopped on its {budget_s:.0f}s budget at iteration {it}, relative residual {rel:.3e}")
            break
        z = apply_M(r); rz_new = (r * z).sum(0)
        p = z + (rz_new / np.where(np.abs(rz) < 1e-300, 1e-300, rz)) * p; rz = rz_new
    residual = float(np.max(np.sqrt((r * r).sum(0)) / scale))
    a = PtP @ (P.T @ (D - Kmul(w) - lam * w))
    if report is not None:
        report.update(iterations=used, relative_residual=residual, seconds=time.time() - start, stopped_on=stopped)
    return w, a

# ------------------------------------------------ anchors (v3 instrument, 437e30e)
ANCHOR_MIN_M = 0.020        # fixed in the pre-registration; this body's median skin depth is 11.0 mm
ANCHOR_MEASURE = "correspondence"   # rebound in __main__: 437e30e's measure, or c0e88a2's bone surface

def exterior_vertices(geometry):
    F = np.asarray(geometry["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    return np.unique(F[ext])

def atlas_bone_surface():
    """Every atlas bone group as one surface, for the corrected 'no bone under it' measure."""
    V, F, off = [], [], 0
    for seg in segments:
        v, f = atlas_mesh(seg); V.append(v); F.append(f + off); off += len(v)
    return np.concatenate(V), np.concatenate(F)

def anchor_distances(measure, canonical, used, atlas_samples=None, block=2048):
    """How far a skin vertex is from bone, in atlas space, by one of the two measures.

    'correspondence' (437e30e): distance to the nearest of the 4,410 sampled bone correspondences.
    Those are SPARSE samples, so this reads ~27 mm even where bone is directly beneath the skin.
    'surface' (c0e88a2): exact distance to the nearest point on any atlas bone surface, which is
    what "no bone under it" means.  Same 20 mm threshold in both."""
    if measure == "correspondence":
        return cKDTree(atlas_samples).query(canonical[used])[0]
    Vb, Fb = atlas_bone_surface()
    project = projector(Vb, Fb, np.random.default_rng(31), n=400000)
    return np.concatenate([project(canonical[used[s:s + block]])[1] for s in range(0, len(used), block)])

def build_anchors(G, Mseg, canonical, used, distance):
    """An anchor for every exterior skin vertex with no bone under it.

    The target is the vertex's image under ITS OWN segment's per-segment similarity, which is the
    only local statement available where there is no bone beneath the skin.  Anchors are an
    assumption, not a measurement, so they never enter the choice of lambda."""
    pick = used[distance > ANCHOR_MIN_M]
    sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes()))
    names = [s["id"] for s in sb["segments"]]
    segment = np.asarray(sb["weights"], np.float32)[pick].argmax(1)
    source = apply(G, canonical[pick]); target = np.empty_like(source)
    for i, name in enumerate(names):
        m = segment == i
        if m.any(): target[m] = apply(Mseg[name], canonical[pick][m])
    counts = {names[i]: int((segment == i).sum()) for i in np.unique(segment)}
    return dict(vertices=pick, source=source, target=target, counts=counts, distance=distance,
                exterior_vertices=int(len(used)), names=names, segment=segment)

def anchor_atlas_samples(Mseg, rest, om):
    """The kept bone correspondence samples in atlas space -- the 437e30e measure's reference."""
    out = []
    for seg in segments:
        i = segments.index(seg); Va, Fa = atlas_mesh(seg); Vo, Fo = om[seg]; Vo_g = apply(rest[seg], Vo)
        n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
        a = sample(Va, Fa, n, np.random.default_rng(10000 + i))
        projector(Vo_g, Fo, np.random.default_rng(20000 + i))(apply(Mseg[seg], a))
        out.append(a)
    return np.concatenate(out)

def stage_anchors():
    """Report the anchor set the rule produces, before anything is fitted with it."""
    reg = json.loads(REGISTRATION.read_text()); Tb, frames_b, binding = bscm.binding_registration()
    G = np.asarray(reg["global_atlas_to_ground"], float)
    Mseg = {s: np.asarray(reg["segments"][s]["atlas_to_ground"], float) for s in reg["segments"]}
    rest = render.OsimModel(B.MODEL).forward(reg["reference_pose_rad"]); om = osim_meshes()
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    geometry = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    canonical = np.asarray(geometry["positions"], float).reshape(-1, 3)
    used = exterior_vertices(geometry); atlas = anchor_atlas_samples(Mseg, rest, om)
    out = {}
    for measure in ("correspondence", "surface"):
        d = anchor_distances(measure, canonical, used, atlas_samples=atlas)
        A = build_anchors(G, Mseg, canonical, used, d)
        disp = np.linalg.norm(A["target"] - A["source"], axis=1)
        label = "437e30e, to the nearest of 4,410 bone SAMPLES" if measure == "correspondence" else \
                "c0e88a2, to the nearest point on any atlas BONE SURFACE"
        say(f"\n== {measure} ({label})")
        say(f"  skin-to-bone distance over {len(used):,} exterior vertices: median {1e3*np.median(d):.1f} mm, "
            f"quartiles {1e3*np.quantile(d,.25):.1f}/{1e3*np.quantile(d,.75):.1f}, max {1e3*d.max():.1f} mm")
        say(f"  anchors (> {1e3*ANCHOR_MIN_M:.0f} mm): {len(A['vertices']):,} of {len(used):,} "
            f"({100*len(A['vertices'])/len(used):.1f}%)  -> {len(atlas) + len(A['vertices']):,} centres in total")
        if len(disp): say(f"  anchor targets move {1e3*np.median(disp):.1f} mm from the global map (median), {1e3*disp.max():.1f} max")
        for name in sorted(A["counts"], key=lambda k: -A["counts"][k]):
            say(f"   {name:10s} {A['counts'][name]:6d}")
        out[measure] = dict(count=int(len(A["vertices"])), by_segment=A["counts"], exterior_vertices=int(len(used)),
                            distance_median_m=float(np.median(d)), distance_max_m=float(d.max()))
    (ROOT / "data/derived/skin-warp-v3/anchor_measures.json").write_text(json.dumps(out, indent=2) + "\n")

# ------------------------------------------------ the flow (v2 instrument, 3b1526c)
def make_flow(base, S, w, a, probe):
    """the field with a squaring count fixed by its own size: per-step displacement <= 1 mm."""
    speed = FlowWarp(base, S, w, a, 1)
    points = S if probe is None else np.vstack([S, probe])
    vmax = float(np.linalg.norm(speed.velocity(points), axis=1).max()) if len(points) else 0.0
    n = 4 if vmax <= 0 else int(max(4, np.ceil(np.log2(max(vmax / STEP_TARGET_M, 1.0)))))
    return FlowWarp(base, S, w, a, 1 << min(n, 10)), vmax

def fit_flow(system, S, T, lam, base, probe):
    """A velocity field whose FLOW hits the correspondences.

    The flow of v is not v, so the field cannot be read off the displacements: each correction is
    solved by the same regularised system against the residual of the ACTUAL flow, and the
    corrections are accumulated into one stationary field.  Same lambda, same system, same
    correspondences as the spline -- only what is being fitted has changed."""
    w, a, hist = np.zeros((len(S), 3)), np.zeros((4, 3)), []
    for _ in range(OUTER):
        F, _ = make_flow(base, S, w, a, probe)
        R = T - F.flow(S)
        hist.append(float(np.abs(R).max()))
        if hist[-1] < OUTER_TOL: break
        dw, da = tps_solve(system, R, lam)
        w, a = w + dw, a + da
    F, vmax = make_flow(base, S, w, a, probe)
    return F, hist, vmax

# ================================================= stage fit
def stage_fit():
    t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    reg = json.loads(REGISTRATION.read_text())
    Tb, frames_b, binding = bscm.binding_registration()
    G = np.asarray(reg["global_atlas_to_ground"], float)
    controls = {}
    say("== known-answer controls ==")
    ok = bool(np.array_equal(G, Tb)); controls["registration_global_is_binding_map"] = ok
    say(f"C1 registration.json global map == the skin bundle's binding map, bit for bit: {ok}")
    if not ok: sys.exit("C1 failed: the warp would not be built on the 0.888 map")
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    canonical = np.asarray(g["positions"], float).reshape(-1, 3)
    Z = FlowWarp.zero(Tb, steps=64) if FAMILY == "flow" else Warp.zero(Tb)
    ok = bool(np.array_equal(Z.apply(canonical), canonical @ Tb[:3, :3].T + Tb[:3, 3])); controls["zero_warp_bitwise"] = ok
    say(f"C2 zero {'flow (64 squaring steps of the identity)' if FAMILY == 'flow' else 'displacement'} reproduces the binding "
        f"map on all {len(canonical):,} skin vertices, bit for bit: {ok}")
    if not ok: sys.exit("C2 failed")
    Z.save(OUT / "zero.npz")
    if FAMILY == "flow":
        # C7, the flow's own known answer: a CONSTANT velocity field flows to a pure translation.
        c = np.array([0.013, -0.021, 0.007]); A = np.zeros((4, 3)); A[3] = c
        const = FlowWarp(Tb, np.zeros((0, 3)), np.zeros((0, 3)), A, 64)
        Y = apply(Tb, canonical[np.random.default_rng(3).choice(len(canonical), 5000, replace=False)])
        err = float(np.abs(const.flow(Y) - (Y + c)).max())
        back = float(np.abs(const.flow(const.flow(Y), inverse=True) - Y).max())
        # 1e-12 m is the double-precision accumulation bound for 2^N additions over metre-scale
        # ground coordinates, not a tolerance fitted to the answer; the exact result is a translation.
        controls.update(constant_velocity_translation_max_m=err, constant_velocity_roundtrip_max_m=back)
        ok = err < 1e-12 and back < 1e-12
        say(f"C7 constant velocity flows to a pure translation: {err:.1e} m (< 1e-12, roundoff); its own inverse returns "
            f"{back:.1e} m -> {'PASS' if ok else 'FAIL'}")
        if not ok: sys.exit("C7 failed")
    if segments != sorted(reg["segments"]) or reg["reference_pose_rad"] != binding["reference_pose_rad"]:
        sys.exit("C3 precondition: segments or reference pose differ from registration.json")
    rest = render.OsimModel(B.MODEL).forward(reg["reference_pose_rad"])
    om = osim_meshes()
    Mseg = {s: np.asarray(reg["segments"][s]["atlas_to_ground"], float) for s in segments}

    # C3: replay the per-segment fit's own evaluation samples; its stored RMS must come back
    ev_samples, worst = {}, 0.0
    for seg in segments:
        r = np.random.default_rng(segments.index(seg) + 1)
        Va, Fa = atlas_mesh(seg); Vo, Fo = om[seg]; Vo_g = apply(rest[seg], Vo)
        n = N_FIT_TORSO if seg == "torso" else N_FIT
        sample(Va, Fa, n, r); sample(Vo_g, Fo, n, r)            # the fit's own draws, discarded
        es, et = sample(Va, Fa, N_EVAL, r), sample(Vo_g, Fo, N_EVAL, r)
        ev_samples[seg] = (es, et)
        rs = residual_points(apply(Mseg[seg], es), et)[0]; rg = residual_points(apply(G, es), et)[0]
        worst = max(worst, abs(rs - reg["segments"][seg]["rms_nearest_surface_segment_m"]),
                    abs(rg - reg["segments"][seg]["rms_nearest_surface_global_m"]))
    ok = worst < 1e-9; controls["gate1_instrument_replay_max_abs_m"] = worst
    say(f"C3 gate 1's instrument reproduces registration.json's stored RMS (global and per-segment, 22 segments): "
        f"max |diff| {worst:.2e} m (< 1e-9) -> {'PASS' if ok else 'FAIL'}")
    if not ok: sys.exit("C3 failed")

    # C4: the projector on points that lie on the mesh, and never further than the nearest sample
    Vo, Fo = om["calcn_l"]; Vo_g = apply(rest["calcn_l"], Vo)
    proj = projector(Vo_g, Fo, np.random.default_rng(99))
    onm = sample(Vo_g, Fo, 2000, np.random.default_rng(98))
    _, d_on, _ = proj(onm)
    off = onm + np.random.default_rng(97).normal(scale=0.01, size=onm.shape)
    _, d_off, d_samp = proj(off)
    ok = d_on.max() < 1e-9 and bool(np.all(d_off <= d_samp + 1e-15))
    controls["projector_on_surface_max_m"] = float(d_on.max()); controls["projector_never_beyond_nearest_sample"] = bool(np.all(d_off <= d_samp + 1e-15))
    say(f"C4 projector: points on the mesh return {d_on.max():.1e} m (< 1e-9); off-surface never beyond the nearest "
        f"sample: {controls['projector_never_beyond_nearest_sample']} -> {'PASS' if ok else 'FAIL'}")
    if not ok: sys.exit("C4 failed")

    # C5: the spline -- interpolates at lambda 0, reproduces an affine field with zero bending,
    # its kernel is conditionally positive definite; C6: its Jacobian against finite differences
    rs_ = np.random.default_rng(5); S5 = rs_.normal(scale=0.3, size=(300, 3)); D5 = rs_.normal(scale=0.01, size=(300, 3))
    sy5 = tps_system(S5); w5, a5 = tps_solve(sy5, D5, 0.0); interp = np.abs(tps_predict(S5, S5, w5, a5) - D5).max()
    Aff = rs_.normal(scale=0.05, size=(4, 3)); Da = np.hstack([S5, np.ones((300, 1))]) @ Aff
    wa, aa = tps_solve(sy5, Da, 1e-3); Wa = Warp(np.eye(4), S5, wa, aa)
    Y5 = rs_.normal(scale=0.3, size=(200, 3)); aff_err = np.abs(Wa.displacement(Y5) - (np.hstack([Y5, np.ones((200, 1))]) @ Aff)).max()
    Wr = Warp(Tb, S5, w5, a5); X5 = np.linalg.solve(Tb[:3, :3], (Y5 - Tb[:3, 3]).T).T; h = 1e-6
    J = Wr.jacobian(X5); Jfd = np.stack([(Wr.apply(X5 + h * e) - Wr.apply(X5 - h * e)) / (2 * h) for e in np.eye(3)], axis=2)
    jac_err = np.abs(J - Jfd).max() / np.abs(J).max()
    ok = interp < 1e-8 and aff_err < 1e-9 and Wa.bending_energy() < 1e-12 and sy5[4].min() > 0 and jac_err < 1e-6
    controls.update(tps_interpolation_max_m=float(interp), tps_affine_max_m=float(aff_err), tps_affine_bending=Wa.bending_energy(),
                    kernel_min_eigenvalue=float(sy5[4].min()), jacobian_vs_finite_difference_rel=float(jac_err))
    say(f"C5 spline: lambda-0 interpolation {interp:.1e} m (< 1e-8); affine field {aff_err:.1e} m (< 1e-9), bending "
        f"{Wa.bending_energy():.1e} (< 1e-12); kernel min eigenvalue on P-perp {sy5[4].min():.2e} (> 0)")
    say(f"C6 Jacobian vs central differences: {jac_err:.1e} relative (< 1e-6) -> {'PASS' if ok else 'FAIL'}")
    if not ok: sys.exit("C5/C6 failed")

    # ---------------------------------------------- correspondences
    say(f"\n== correspondences: atlas bone group -> nearest point on the scaffold bone mesh to M_seg a  [{time.time()-t0:.0f}s] ==")
    S, T, L, per, atlas_kept = [], [], [], {}, []
    for seg in segments:
        i = segments.index(seg); Va, Fa = atlas_mesh(seg); Vo, Fo = om[seg]; Vo_g = apply(rest[seg], Vo)
        n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
        a = sample(Va, Fa, n, np.random.default_rng(10000 + i))
        p = apply(Mseg[seg], a); t, d, _ = projector(Vo_g, Fo, np.random.default_rng(20000 + i))(p)
        s = apply(G, a)                      # every sample is kept: the trim is gone (8caf15b)
        S.append(s); T.append(t); L += [seg] * len(a); atlas_kept.append(a)
        per[seg] = dict(sampled=n, kept=int(len(a)), projection_median_mm=1e3 * float(np.median(d)),
                        projection_max_mm=1e3 * float(d.max()),
                        displacement_from_global_median_mm=1e3 * float(np.median(np.linalg.norm(t - s, axis=1))),
                        displacement_from_global_max_mm=1e3 * float(np.linalg.norm(t - s, axis=1).max()))
        say(f"  {seg:10s} kept {len(a):4d}/{n}  |t - M a| median {per[seg]['projection_median_mm']:5.2f} mm   "
            f"|t - G a| median {per[seg]['displacement_from_global_median_mm']:5.1f} max {per[seg]['displacement_from_global_max_mm']:5.1f} mm")
    S, T, L = np.concatenate(S), np.concatenate(T), np.asarray(L); D = T - S
    atlas_kept = np.concatenate(atlas_kept)
    say(f"  {len(S)} correspondences")

    # ---------------------------------------------- the regularisation rule
    say(f"\n== lambda by {FOLDS}-fold cross-validation on the bone correspondences  [{time.time()-t0:.0f}s] ==")
    probe = apply(G, canonical[np.random.default_rng(11).choice(len(canonical), 20000, replace=False)])
    perm = np.random.default_rng(7).permutation(len(S)); folds = np.array_split(perm, FOLDS)
    sse = np.zeros(len(LAMBDAS))
    for f in folds:
        tr = np.setdiff1d(perm, f); sy = tps_system(S[tr])
        for li, lam in enumerate(LAMBDAS):
            if FAMILY == "flow":
                F, _, _ = fit_flow(sy, S[tr], T[tr], lam, G, probe)
                sse[li] += ((F.flow(S[f]) - T[f]) ** 2).sum()
            else:
                w, a = tps_solve(sy, D[tr], lam)
                sse[li] += ((S[f] + tps_predict(S[f], S[tr], w, a) - T[f]) ** 2).sum()
        say(f"  fold done [{time.time()-t0:.0f}s]")
    cv = np.sqrt(sse / len(S)); best = cv.min()
    lam = max(l for l, e in zip(LAMBDAS, cv) if e <= best * (1 + TIE))
    for l, e in zip(LAMBDAS, cv): say(f"  lambda {l:9.2e}  CV RMS {1e3*e:7.3f} mm" + ("   <- chosen" if l == lam else ""))
    say(f"  rule: largest lambda within {100*TIE:.0f}% of the minimum CV RMS ({1e3*best:.3f} mm) -> lambda = {lam:.2e}"
        + ("   (the grid's edge)" if lam in (LAMBDAS[0], LAMBDAS[-1]) else ""))
    sy = tps_system(S)
    rule = "largest lambda within 1% of the minimum 5-fold CV RMS on the bone correspondences"
    flow_report = anchor_report = None
    if FAMILY == "flow":
        W, hist, vmax = fit_flow(sy, S, T, lam, G, probe)
        fit_res = np.linalg.norm(W.flow(S) - T, axis=1)
        step_max = W.maximum_step_displacement(np.vstack([S, apply(G, canonical)]))
        flow_report = dict(steps=W.steps, squarings=int(np.log2(W.steps)), max_velocity_m=vmax,
                           max_step_displacement_m=step_max, step_budget_m=STEP_TARGET_M,
                           outer_max_residual_m=hist, outer_iterations=len(hist))
        W.meta = dict(kind="flow", lam=lam, correspondences=int(len(S)), steps=W.steps, rule=rule)
        say(f"  flow: {W.steps} steps (2^{int(np.log2(W.steps))} squarings), max |v| {1e3*vmax:.1f} mm, "
            f"max per-step displacement {1e3*step_max:.3f} mm (budget {1e3*STEP_TARGET_M:.1f} mm)")
        say(f"  outer corrections, worst |flow - target| per pass: " + ", ".join(f"{1e3*h:.2f}" for h in hist) + " mm")
    elif FAMILY in ("anchored", "anchored-surface"):
        # C9: with the anchor set EMPTY this pipeline must BE the 96e5f1b spline.
        w0, a0 = tps_solve(sy, D, lam); W0 = Warp(G, S, w0, a0)
        v1 = Warp.load(ROOT / "data/derived/skin-warp-v1/warp.npz")
        d9 = float(np.abs(W0.apply(canonical) - v1.apply(canonical)).max())
        controls["anchor_empty_reproduces_v1_max_m"] = d9
        say(f"C9 anchor set empty reproduces the 96e5f1b spline on all {len(canonical):,} skin vertices: "
            f"{d9:.2e} m -> {'PASS' if d9 == 0.0 else 'FAIL'}")
        if d9 != 0.0: sys.exit("C9 failed")
        # C10: the anchored system is 38,821 centres -- 12 GB as a dense factorisation, which this
        # machine does not have free -- so it is solved matrix-free.  Same equations; prove it here,
        # where the direct answer is available.
        rep10 = {}; w1, a1 = tps_solve_matrix_free(S, D, lam, report=rep10)
        d10 = float(np.abs(Warp(G, S, w1, a1).apply(canonical) - W0.apply(canonical)).max())
        controls.update(matrix_free_vs_direct_max_m=d10, matrix_free_iterations=rep10["iterations"],
                        matrix_free_relative_residual=rep10["relative_residual"])
        say(f"C10 matrix-free PCG vs the direct solve, same bone system: {d10:.2e} m over the whole skin "
            f"({rep10['iterations']} iterations, relative residual {rep10['relative_residual']:.1e})")
        used_v = exterior_vertices(g)
        A = build_anchors(G, Mseg, canonical, used_v,
                          anchor_distances(ANCHOR_MEASURE, canonical, used_v, atlas_samples=atlas_kept))
        disp = np.linalg.norm(A["target"] - A["source"], axis=1)
        say(f"\n== anchors: skin vertices more than {1e3*ANCHOR_MIN_M:.0f} mm (atlas) from "
            f"{'any sampled bone correspondence' if ANCHOR_MEASURE == 'correspondence' else 'the nearest point on any atlas bone surface'} ==")
        say(f"  {len(A['source']):,} of {A['exterior_vertices']:,} exterior vertices "
            f"({100*len(A['source'])/A['exterior_vertices']:.1f}%); {len(S) + len(A['source']):,} centres in total")
        say(f"  their targets move {1e3*np.median(disp):.1f} mm from the global map (median), {1e3*disp.max():.1f} max")
        for name in sorted(A["counts"], key=lambda k: -A["counts"][k]):
            say(f"   {name:10s} {A['counts'][name]:6d}")
        anchor_report = dict(threshold_m=ANCHOR_MIN_M, measure=ANCHOR_MEASURE, count=int(len(A["source"])), by_segment=A["counts"],
                             exterior_vertices=A["exterior_vertices"],
                             displacement_median_m=float(np.median(disp)), displacement_max_m=float(disp.max()))
        S_all, T_all = np.vstack([S, A["source"]]), np.vstack([T, A["target"]])
        say(f"\n== the anchored fit: {len(S_all):,} centres, lambda {lam:.2e}, matrix-free  [{time.time()-t0:.0f}s] ==")
        rep = {}; w, a = tps_solve_matrix_free(S_all, T_all - S_all, lam, report=rep)
        anchor_report.update(solver="projected PCG with local cardinal-function preconditioner",
                             iterations=rep["iterations"], relative_residual=rep["relative_residual"])
        say(f"  solved in {rep['iterations']} iterations, relative residual {rep['relative_residual']:.1e}  [{time.time()-t0:.0f}s]")
        W = Warp(G, S_all, w, a, dict(kind="thin-plate spline with anchors", lam=lam,
                                      correspondences=int(len(S_all)), anchors=int(len(A["source"])), rule=rule))
        fit_res = np.linalg.norm(W.apply(np.linalg.solve(G[:3, :3], (S_all - G[:3, 3]).T).T) - T_all, axis=1)
        bone_res = fit_res[:len(S)]
        say(f"  residual at the BONE correspondences RMS {1e3*np.sqrt((bone_res**2).mean()):.3f} mm, max {1e3*bone_res.max():.2f} mm")
    else:
        w, a = tps_solve(sy, D, lam)
        W = Warp(G, S, w, a, dict(kind="thin-plate spline displacement on the binding map", lam=lam,
                                  correspondences=int(len(S)), rule=rule))
        fit_res = np.linalg.norm(W.apply(np.linalg.solve(G[:3, :3], (S - G[:3, 3]).T).T) - T, axis=1)
    W.save(OUT / "warp.npz")
    say(f"  fitted: residual at the correspondences RMS {1e3*np.sqrt((fit_res**2).mean()):.3f} mm, max {1e3*fit_res.max():.2f} mm; "
        f"bending energy {W.bending_energy():.4e}")

    # ---------------------------------------------- gate 1
    say(f"\n== GATE 1: warped atlas bone group vs the per-segment similarity (RMS nearest-surface, fit's own eval samples) ==")
    say(f"{'segment':10s} {'global':>8s} {'per-seg':>8s} {'warped':>8s} {'limit':>8s}   mm")
    g1 = {}
    for seg in segments:
        es, et = ev_samples[seg]
        rw = residual_points(W.apply(es), et)[0]; rs = reg["segments"][seg]["rms_nearest_surface_segment_m"]
        g1[seg] = dict(global_m=reg["segments"][seg]["rms_nearest_surface_global_m"], per_segment_m=rs, warped_m=rw,
                       limit_m=rs + GATE1_MARGIN_M, pass_=bool(rw <= rs + GATE1_MARGIN_M))
        say(f"{seg:10s} {1e3*g1[seg]['global_m']:8.2f} {1e3*rs:8.2f} {1e3*rw:8.2f} {1e3*(rs+GATE1_MARGIN_M):8.2f}   "
            + ("PASS" if g1[seg]["pass_"] else "FAIL"))
    gate1 = all(v["pass_"] for v in g1.values())
    say(f"GATE 1 -> {'PASS' if gate1 else 'FAIL'} ({sum(v['pass_'] for v in g1.values())}/{len(g1)} segments)")
    if FAMILY == "flow":
        # C8, required because a flow is not a spline: forward then inverse must return EVERY skin
        # vertex to itself.  The inverse is the exact inverse of each step, not a second field.
        say(f"\n== C8: forward then inverse on all {len(canonical):,} skin vertices  [{time.time()-t0:.0f}s] ==")
        there = W.apply(canonical)
        back = W.flow(there, inverse=True)
        err = float(np.abs(back - W.ground(canonical)).max())
        controls["forward_inverse_max_m"] = err
        ok = err < 1e-9
        say(f"C8 max |inverse(forward(x)) - x| = {err:.2e} m (< 1e-9) -> {'PASS' if ok else 'FAIL'}   [{time.time()-t0:.0f}s]")
        if not ok: sys.exit("C8 failed")
    (OUT / "fit.json").write_text(json.dumps(dict(schema="ihm.skin-warp-fit.v1", family=FAMILY, flow=flow_report,
        anchors=anchor_report, controls=controls,
        correspondences=per, n_correspondences=int(len(S)), cv=dict(lambdas=LAMBDAS, rms_m=cv.tolist(), chosen=lam, tie=TIE, folds=FOLDS),
        bending_energy=W.bending_energy(), gate1=dict(passed=gate1, margin_m=GATE1_MARGIN_M, segments=g1),
        provenance=dict(script=str(Path(__file__).relative_to(ROOT)), script_sha256=sha(__file__),
                        inputs={str(Path(p).relative_to(ROOT)): sha(p) for p in (REGISTRATION, B.ANATOMY, B.BINDING, B.MODEL)})),
        indent=2) + "\n")
    say(f"wrote {(OUT/'warp.npz').relative_to(ROOT)}, {(OUT/'fit.json').relative_to(ROOT)}   [{time.time()-t0:.0f}s]")

# ================================================= stage score
def stage_score():
    fit = json.loads((OUT / "fit.json").read_text()); W = load_warp(OUT / "warp.npz")
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    used = np.unique(F[ext])
    # gate 4: Jacobian determinant at every exterior skin vertex; triangle orientation before/after
    det0 = np.linalg.det(W.base[:3, :3])
    if hasattr(W, "steps"):
        # the flow's determinant is the product of its per-step determinants; the smallest single
        # step is the numerical statement that the squaring is fine enough to forbid a fold
        det, worst_step = W.jacobian_determinant(V[used])
    else:
        det = np.concatenate([np.linalg.det(W.jacobian(V[used[s:s + 2048]])) for s in range(0, len(used), 2048)])
        worst_step = None
    Y0, Y1 = np.zeros_like(V), np.zeros_like(V); Y0[used] = W.ground(V[used]); Y1[used] = W.apply(V[used])
    def normals(Y, f): t = Y[f]; return np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    dots = (normals(Y0, F[ext]) * normals(Y1, F[ext])).sum(1)
    # the caps of the whole-skin surface gate 2 measures: invented surface, reported beside
    sv, sf = V[used], np.searchsorted(used, F[ext])
    cv_, cf_ = bscm.repair(sv, sf); cvv, cff, _ = bscm.cap_boundaries(cv_, cf_)
    cap = np.any(cff >= len(cv_), axis=1)
    cdots = (normals(W.ground(cvv), cff[cap]) * normals(W.apply(cvv), cff[cap])).sum(1)
    folds_v, folds_t = int((det <= 0).sum()), int((dots < 0).sum())
    gate4 = folds_v == 0 and folds_t == 0
    # reported: per-segment area change, over the bundle's own partition
    sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes())); names = [s["id"] for s in sb["segments"]]
    Wt = np.asarray(sb["weights"], np.float32); owner = ((Wt[F[:, 0]] + Wt[F[:, 1]] + Wt[F[:, 2]]) / 3).argmax(1)[ext]
    a0 = np.linalg.norm(normals(Y0, F[ext]), axis=1) / 2; a1 = np.linalg.norm(normals(Y1, F[ext]), axis=1) / 2
    area = {n: dict(before_m2=float(a0[owner == i].sum()), after_m2=float(a1[owner == i].sum()),
                    ratio=float(a1[owner == i].sum() / a0[owner == i].sum())) for i, n in enumerate(names) if (owner == i).any()}
    # zero-warp controls on the two instruments that own gates 2 and 3
    ez = json.loads((OUT / "enclosure_zero.json").read_text()); ew = json.loads((OUT / "enclosure_warp.json").read_text())
    zero_enc = max(abs(ez["warped"][s] - ez["binding_map"][s]) for s in ez["segments"])
    ref = {r["body"]: r for r in json.loads((REFERENCE_BUNDLE / "manifest.json").read_text())["records"]}
    zb = {r["body"]: r for r in json.loads((OUT / "bundle-zero/manifest.json").read_text())["records"]}
    zero_heel = {s: (zb[s]["skin_minus_bone_minimum_y_m"], ref[s]["skin_minus_bone_minimum_y_m"]) for s in FOOT}
    zero_ok = zero_enc == 0 and all(a == b for a, b in zero_heel.values()) and round(ez["mean"]["binding_map"], 3) == 0.888 \
        and round(ez["mean"]["own_bones_ceiling"], 3) == 0.997
    wb = {r["body"]: r for r in json.loads((ROOT / BUNDLE / "manifest.json").read_text())["records"]}
    heel = {s: wb[s]["skin_minus_bone_minimum_y_m"] for s in ("calcn_l", "calcn_r")}
    gate2 = ew["mean"]["warped"] >= ENCLOSURE_MEAN and all(ew["warped"][s] >= ENCLOSURE_FOOT for s in FOOT)
    gate3 = all(HEEL_RANGE_M[0] <= v <= HEEL_RANGE_M[1] for v in heel.values())
    gate1 = fit["gate1"]["passed"]
    say("== zero-warp controls (the same instruments, displacement set to zero) ==")
    say(f"  enclosure: max |zero-warp - binding map| over 22 segments {zero_enc:.3f}; binding mean {ez['mean']['binding_map']:.3f} "
        f"(0.888), ceiling {ez['mean']['own_bones_ceiling']:.3f} (0.997)")
    for s, (a, b) in zero_heel.items(): say(f"  {s:8s} skin - bone minimum y: zero-warp {1e3*a:+.1f} mm, skin-binding bundle {1e3*b:+.1f} mm")
    say(f"  -> {'PASS' if zero_ok else 'FAIL'}")
    g1w = max(fit["gate1"]["segments"].items(), key=lambda kv: kv[1]["warped_m"] - kv[1]["limit_m"])
    say("\n== the four gates (93d3d57) ==")
    say(f"  1 bones: worst segment {g1w[0]} {1e3*g1w[1]['warped_m']:.2f} mm vs limit {1e3*g1w[1]['limit_m']:.2f} mm -> {'PASS' if gate1 else 'FAIL'}")
    say(f"  2 enclosure: mean {ew['mean']['warped']:.3f} (>= {ENCLOSURE_MEAN}); " + ", ".join(f"{s} {ew['warped'][s]:.3f}" for s in FOOT)
        + f" (each >= {ENCLOSURE_FOOT}) -> {'PASS' if gate2 else 'FAIL'}")
    say(f"  3 heel: " + ", ".join(f"{s} {1e3*v:+.1f} mm" for s, v in heel.items()) + f" (in [{1e3*HEEL_RANGE_M[0]:.0f}, {1e3*HEEL_RANGE_M[1]:.0f}] mm) -> {'PASS' if gate3 else 'FAIL'}")
    say(f"  4 folding: det J <= 0 at {folds_v} of {len(det):,} skin vertices (min det {det.min():.4f}, global map {det0:.4f}); "
        f"{folds_t} of {len(dots):,} skin triangles inverted -> {'PASS' if gate4 else 'FAIL'}")
    if worst_step is not None:
        f_ = fit.get("flow") or {}
        say(f"     flow: {W.steps} steps (2^{int(np.log2(W.steps))} squarings), smallest single-step det "
            f"{worst_step:.6f} (> 0 is what forbids the fold), max per-step displacement "
            f"{1e3*(f_.get('max_step_displacement_m') or 0):.3f} mm")
    say(f"     caps (invented surface, not gated): {int((cdots < 0).sum())} of {len(cdots)} inverted")
    # where the folds are (reported after the verdict, not gated): by partition segment, and how far from
    # the nearest spline centre -- a fold pulled into the skin by the bones under it sits close to them
    vseg = np.asarray(sb["weights"], np.float32).argmax(1)
    bad_v = used[det <= 0]; bad_t = np.where(dots < 0)[0]
    near = cKDTree(W.centres).query(W.ground(V[bad_v]))[0] if len(bad_v) else np.zeros(0)
    fold_where = dict(vertices={names[i]: int((vseg[bad_v] == i).sum()) for i in np.unique(vseg[bad_v])},
                      triangles={names[i]: int((owner[bad_t] == i).sum()) for i in np.unique(owner[bad_t])},
                      vertex_to_nearest_centre_mm=[1e3 * float(q) for q in np.quantile(near, (0, .5, 1))] if len(near) else None)
    say(f"     folded vertices by segment {fold_where['vertices']}; inverted triangles by segment {fold_where['triangles']}")
    if len(near): say(f"     folded vertex to nearest spline centre min/median/max {near.min()*1e3:.1f}/{np.median(near)*1e3:.1f}/{near.max()*1e3:.1f} mm")
    say(f"\n  reported: bending energy {W.bending_energy():.4e}; skin area {a0.sum():.4f} -> {a1.sum():.4f} m2")
    for n, v in area.items(): say(f"     {n:10s} area x{v['ratio']:.3f}")
    say("  enclosure per segment (binding -> warped): " + ", ".join(f"{s} {ew['binding_map'][s]:.3f}->{ew['warped'][s]:.3f}" for s in ew["segments"]))
    verdict = all((gate1, gate2, gate3, gate4))
    say(f"\nVERDICT (all four gates): {'PASS' if verdict else 'FAIL'}")
    (OUT / "verdict.json").write_text(json.dumps(dict(schema="ihm.skin-warp-verdict.v1", verdict=verdict, zero_warp_controls=dict(
        passed=zero_ok, enclosure_max_abs=zero_enc, heel={s: dict(zero_warp_m=a, skin_binding_m=b) for s, (a, b) in zero_heel.items()},
        binding_mean=ez["mean"]["binding_map"], ceiling_mean=ez["mean"]["own_bones_ceiling"]),
        gate1=dict(passed=gate1, worst=g1w[0]), gate2=dict(passed=gate2, mean=ew["mean"]["warped"], segments=ew["warped"]),
        gate3=dict(passed=gate3, skin_minus_bone_minimum_y_m=heel, range_m=HEEL_RANGE_M),
        gate4=dict(passed=gate4, vertices=len(det), nonpositive_det=folds_v, min_det=float(det.min()), global_det=float(det0),
                   triangles=len(dots), inverted=folds_t, caps=len(cdots), caps_inverted=int((cdots < 0).sum()),
                   steps=getattr(W, "steps", None), min_step_det=worst_step, folds=fold_where),
        reported=dict(bending_energy=W.bending_energy(), area=area, area_total_before_m2=float(a0.sum()), area_total_after_m2=float(a1.sum()))),
        indent=2) + "\n")

# ================================================= stage diagnose (after the verdict; reported, never used to refit)
def inside_mask(skin_vertices, skin_faces, bone_points, samples=500, seed=0):
    """bscm.enclosure with the same draws, returning which tested points are inside, not only the mean."""
    generator = np.random.default_rng(seed); points = np.asarray(bone_points, float)
    if len(points) > samples: points = points[generator.choice(len(points), samples, replace=False)]
    triangles = np.asarray(skin_vertices, float)[np.asarray(skin_faces)]
    a = triangles[:, 0]; edge1 = triangles[:, 1] - a; edge2 = triangles[:, 2] - a
    direction = generator.normal(size=3); direction /= np.linalg.norm(direction); counts = []
    for sign in (1., -1.):
        d = sign * direction; pvec = np.cross(d, edge2); det = (edge1 * pvec).sum(axis=1)
        parallel = np.abs(det) < 1e-14; inverse = np.where(parallel, 0., 1. / np.where(parallel, 1., det))
        hits = np.zeros(len(points), dtype=np.int64)
        for start in range(0, len(points), 256):
            block = points[start:start + 256]; tvec = block[:, None, :] - a[None, :, :]
            u = (tvec * pvec[None, :, :]).sum(axis=2) * inverse[None, :]; qvec = np.cross(tvec, edge1[None, :, :])
            v = (qvec * d).sum(axis=2) * inverse[None, :]; t = (qvec * edge2[None, :, :]).sum(axis=2) * inverse[None, :]
            hits[start:start + 256] = ((~parallel[None, :]) & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 1e-9)).sum(axis=1)
        counts.append(hits % 2 == 1)
    return points, counts[0] & counts[1]

def stage_diagnose():
    """Where the bone points that are NOT inside the skin sit, in their own segment's frame
    (OpenSim: x anterior/distal along the foot, y superior, z right), and how far outside."""
    W = load_warp(OUT / "warp.npz"); Tb, frames_b, _ = bscm.binding_registration()
    ew = json.loads((OUT / "enclosure_warp.json").read_text())
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    used, inv = np.unique(F[ext], return_inverse=True); sv, sf = V[used], inv.reshape(-1, 3)
    if bscm.simtk_precondition(sv, sf) is not None:
        cv_, cf_ = bscm.repair(sv, sf); sv, sf, _ = bscm.cap_boundaries(cv_, cf_)
    osim = bscm.bone_clouds(); out = {}
    for seg in ("toes_l", "toes_r", "calcn_l", "calcn_r", "pelvis", "hand_l", "hand_r"):
        M = frames_b[seg]; bones = osim[seg] @ M[:3, :3].T + M[:3, 3]; lb = osim[seg]
        say(f"== {seg}   bone extent in its frame: x {1e3*lb[:,0].min():+.0f}..{1e3*lb[:,0].max():+.0f}  "
            f"y {1e3*lb[:,1].min():+.0f}..{1e3*lb[:,1].max():+.0f}  z {1e3*lb[:,2].min():+.0f}..{1e3*lb[:,2].max():+.0f} mm")
        out[seg] = {}
        for name, S in (("binding", W.ground(sv)), ("warped", W.apply(sv))):
            pts, ins = inside_mask(S, sf, bones)
            check = ew["binding_map" if name == "binding" else "warped"][seg]
            if float(ins.mean()) != check: sys.exit(f"{seg} {name}: mask mean {ins.mean()} != enclosure() {check}; diagnosis void")
            o = pts[~ins]; local = (o - M[:3, 3]) @ M[:3, :3]
            dist = cKDTree(S).query(o)[0] if len(o) else np.zeros(0)
            rec = dict(inside=float(ins.mean()), outside=int((~ins).sum()), tested=int(len(pts)))
            if len(o):
                rec.update(nearest_skin_vertex_mm=dict(median=1e3 * float(np.median(dist)), max=1e3 * float(dist.max())),
                           outside_x_mm=[1e3 * float(q) for q in np.quantile(local[:, 0], (0, .5, 1))],
                           outside_y_mm=[1e3 * float(q) for q in np.quantile(local[:, 1], (0, .5, 1))],
                           outside_z_mm=[1e3 * float(q) for q in np.quantile(local[:, 2], (0, .5, 1))])
            out[seg][name] = rec
            say(f"   {name:8s} inside {rec['inside']:.3f} (= enclosure()); {rec['outside']:3d}/{rec['tested']} outside"
                + ("" if not len(o) else
                   f"; nearest skin vertex median {rec['nearest_skin_vertex_mm']['median']:.1f} max {rec['nearest_skin_vertex_mm']['max']:.1f} mm; "
                   f"outside x {rec['outside_x_mm'][0]:+.0f}/{rec['outside_x_mm'][1]:+.0f}/{rec['outside_x_mm'][2]:+.0f}  "
                   f"y {rec['outside_y_mm'][0]:+.0f}/{rec['outside_y_mm'][1]:+.0f}/{rec['outside_y_mm'][2]:+.0f}  "
                   f"z {rec['outside_z_mm'][0]:+.0f}/{rec['outside_z_mm'][1]:+.0f}/{rec['outside_z_mm'][2]:+.0f} mm (min/median/max)"))
    (OUT / "diagnose.json").write_text(json.dumps(out, indent=2) + "\n")

def choose_lambda(S, T):
    """The line's own rule, factored out so the recovery control uses the identical one."""
    D = T - S
    perm = np.random.default_rng(7).permutation(len(S)); fold_idx = np.array_split(perm, FOLDS)
    sse = np.zeros(len(LAMBDAS))
    for f in fold_idx:
        tr = np.setdiff1d(perm, f); sy = tps_system(S[tr])
        for li, lam in enumerate(LAMBDAS):
            w, a = tps_solve(sy, D[tr], lam)
            sse[li] += ((S[f] + tps_predict(S[f], S[tr], w, a) - T[f]) ** 2).sum()
    cv = np.sqrt(sse / len(S)); best = cv.min()
    return max(l for l, e in zip(LAMBDAS, cv) if e <= best * (1 + TIE)), cv

# Two magnitudes, because the d^2/R bias grows with the displacement: the chest wall's, so the two
# lines are comparable, and this line's own, which is what its gates are actually read at.
RECOVERY_SCALES = (("chest-wall scale", 0.0015), ("this line's scale", 0.020))

def stage_recovery():
    """The known answer this line never had (d626af3).

    Every gate here compares one fit to another fit.  This displaces THIS BODY'S OWN bone groups by
    a field we chose -- a thin-plate spline displacement, the same family the warp is fitted in --
    and asks whether the pipeline returns it.  Recovery is reported against the TRUTH, pointwise
    and to the surface, beside the residual at the correspondences, which is agreement with the
    targets and is the quantity that has been flattering this line."""
    t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    meshes = {seg: atlas_mesh(seg) for seg in segments}
    allV, allF = atlas_bone_surface()
    centres = sample(allV, allF, 200, np.random.default_rng(102))
    w0 = np.random.default_rng(101).normal(size=(len(centres), 3))
    P = np.hstack([centres, np.ones((len(centres), 1))]); Q, _ = np.linalg.qr(P)
    w0 -= Q @ (Q.T @ w0)                      # P^T w = 0: a pure spline displacement, no affine part
    probe = sample(allV, allF, 20000, np.random.default_rng(103))
    unit = float(np.linalg.norm(Warp(np.eye(4), centres, w0, np.zeros((4, 3))).displacement(probe), axis=1).mean())
    out = {}
    for label, magnitude in RECOVERY_SCALES:
        truth = Warp(np.eye(4), centres, w0 * (magnitude / unit), np.zeros((4, 3)))
        moved = float(np.linalg.norm(truth.displacement(probe), axis=1).mean())
        say(f"\n== recovery control, {label}: known field moves this body's bones {1e3*moved:.2f} mm on average "
            f"[{time.time()-t0:.0f}s]")
        S, T, ES, TRUTH_ES, seg_of_es = [], [], [], [], []
        for seg in segments:
            i = segments.index(seg); Va, Fa = meshes[seg]
            Vd = truth.apply(Va)                      # the displaced bone plays the scaffold's part
            n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
            a = sample(Va, Fa, n, np.random.default_rng(10000 + i))
            # this line's correspondence: the nearest point ON the displaced surface, from the
            # starting map's image -- the identity here, as the global similarity is in the real one
            t, d, _ = projector(Vd, Fa, np.random.default_rng(20000 + i))(a)
            keep = d <= np.quantile(d, 1 - 0.10)   # RETIRED control: its 10% is frozen as a literal so
            S.append(a[keep]); T.append(t[keep])   # the numbers already on the record stay reproducible
            e = sample(Va, Fa, 2000, np.random.default_rng(30000 + i))
            ES.append(e); TRUTH_ES.append(truth.apply(e)); seg_of_es += [seg] * len(e)
        S, T = np.concatenate(S), np.concatenate(T); ES, TRUTH_ES = np.concatenate(ES), np.concatenate(TRUTH_ES)
        seg_of_es = np.asarray(seg_of_es)
        # the correspondence as an instrument, before anything is fitted with it
        terr = np.linalg.norm(T - truth.apply(S), axis=1)
        say(f"   target error (nearest-point target vs the known truth): mean {1e3*terr.mean():.3f} mm, "
            f"p90 {1e3*np.quantile(terr, .9):.3f}, max {1e3*terr.max():.3f}")
        lam, cv = choose_lambda(S, T)
        say(f"   lambda by the same rule: {lam:.2e}  (CV min {1e3*cv.min():.3f} mm)  [{time.time()-t0:.0f}s]")
        w, a = tps_solve(tps_system(S), T - S, lam)
        W = Warp(np.eye(4), S, w, a)
        fit_res = np.linalg.norm(W.apply(S) - T, axis=1)
        got = W.apply(ES); point = np.linalg.norm(got - TRUTH_ES, axis=1)
        surface = np.zeros(len(ES))
        for seg in segments:
            m = seg_of_es == seg
            Va, Fa = meshes[seg]
            surface[m] = projector(truth.apply(Va), Fa, np.random.default_rng(40000 + segments.index(seg)))(got[m])[1]
        rec = dict(known_field_mean_displacement_m=moved, lam=lam,
                   target_error_mean_m=float(terr.mean()), target_error_p90_m=float(np.quantile(terr, .9)),
                   target_error_max_m=float(terr.max()),
                   residual_at_correspondences_rms_m=float(np.sqrt((fit_res ** 2).mean())),
                   recovery_pointwise_rms_m=float(np.sqrt((point ** 2).mean())),
                   recovery_pointwise_p90_m=float(np.quantile(point, .9)), recovery_pointwise_max_m=float(point.max()),
                   recovery_to_surface_rms_m=float(np.sqrt((surface ** 2).mean())),
                   recovery_to_surface_p90_m=float(np.quantile(surface, .9)))
        say(f"   residual at the correspondences (agreement with TARGETS): {1e3*rec['residual_at_correspondences_rms_m']:.3f} mm RMS")
        say(f"   RECOVERY vs the TRUTH, pointwise: {1e3*rec['recovery_pointwise_rms_m']:.3f} mm RMS, "
            f"p90 {1e3*rec['recovery_pointwise_p90_m']:.3f}, max {1e3*rec['recovery_pointwise_max_m']:.3f}")
        say(f"   RECOVERY to the surface: {1e3*rec['recovery_to_surface_rms_m']:.3f} mm RMS, "
            f"p90 {1e3*rec['recovery_to_surface_p90_m']:.3f}   [{time.time()-t0:.0f}s]")
        out[label] = rec
    (OUT / "recovery.json").write_text(json.dumps(dict(schema="ihm.skin-warp-recovery.v1", scales=out,
        basis="a known thin-plate-spline displacement of this body's own bone groups, fitted by this line's own "
              "correspondences (nearest point on the displaced surface) and lambda rule; recovery measured against "
              "the known field, not against the targets"), indent=2) + "\n")

# The control has to sit at the separation the pipeline operates at, and separation is what is LEFT
# after the per-segment similarity.  A field with 200 centres over the whole body is locally almost a
# similarity on any single bone, so the starting map absorbs it and the separation comes out at
# 0.09 mm -- ten times below the pipeline's 0.78-6.88 mm.  The field therefore needs bone-scale
# spatial frequency, and its amplitude is calibrated to put the separation inside the measured range.
FIELD_CENTRES = 2000
SEPARATION_TARGETS = (0.001, 0.003, 0.007)   # metres, spanning the pipeline's own 0.78-6.88 mm

# Normal shooting (fec11d2/5239b2a), adopted explicitly with its settings declared here rather than
# defaulted: the module's cap defaults to 20 mm, which is marginal against the 12.32 mm separations
# measured at amplitude 3, and its agreement filter defaults to OFF.
SHOOT_CAP_M = 0.050          # 4x the largest measured per-segment separation (12.32 mm)
SHOOT_AGREEMENT = 0.0        # ON: drop a target whose normal opposes the source's -- the thin-sheet far wall
SHOOT_RETURN_TOL_M = 1e-3    # the module's own return test, unchanged

def correspond(mode, Va, Fa, Vd, M0, n, seed_a, seed_b, keep_target):
    """One correspondence, two rules, everything downstream identical.

    Returns the ATLAS-space source points, their targets on the displaced surface, the point the
    rule started from, and a census of what it dropped and why."""
    if mode == "nearest":
        a = sample(Va, Fa, n, np.random.default_rng(seed_a))
        start = apply(M0, a)
        tgt, d, _ = projector(Vd, Fa, np.random.default_rng(seed_b))(start)
        if len(a) > keep_target:             # no trim; uniform thinning only, for count-matching
            idx = np.random.default_rng(seed_b + 5).choice(len(a), keep_target, replace=False)
            a, tgt, start = a[idx], tgt[idx], start[idx]
        return a, tgt, start, dict(rule="nearest point on the surface, no trim",
                                   sampled=int(n), kept=int(len(a)))
    if mode == "nearest_untrimmed":
        # 2x2 cell A (269f1bb): the nearest arm with the trim REMOVED and nothing else changed --
        # same sampler, same targets, same fitted count, thinned uniformly instead of by separation.
        a = sample(Va, Fa, n, np.random.default_rng(seed_a))
        start = apply(M0, a)
        tgt, _, _ = projector(Vd, Fa, np.random.default_rng(seed_b))(start)
        if len(a) > keep_target:
            idx = np.random.default_rng(seed_b + 5).choice(len(a), keep_target, replace=False)
            a, tgt, start = a[idx], tgt[idx], start[idx]
        return a, tgt, start, dict(rule="nearest point, NO separation trim, uniform thinning (cell A)",
                                   sampled=int(n), survived=int(n), kept=int(len(a)))
    if mode in ("hybrid", "hybrid_trimmed"):
        # 380df28: the shot's target where it survives, the nearest rule's where it does not, over ONE
        # population -- so coverage is 100% by construction and only the targets differ.
        from ihm.anatomy.normal_shooting import shoot_pairs, surface_samples
        src_V = apply(M0, Va)
        samples = surface_samples(src_V, Fa, n, seed_a)
        P = samples[0]
        r = shoot_pairs(src_V, Fa, Vd, Fa, samples=samples, cap_m=SHOOT_CAP_M,
                        return_tol_m=SHOOT_RETURN_TOL_M, min_normal_agreement=SHOOT_AGREEMENT)
        keep = r["keep"]; tgt = np.empty_like(P)
        if keep.any(): tgt[keep] = r["target"]
        if (~keep).any():
            tgt[~keep] = projector(Vd, Fa, np.random.default_rng(seed_b))(P[~keep])[0]
        if mode == "hybrid_trimmed":           # 2x2 cell B: hybrid targets WITH the nearest arm's trim
            sep = np.linalg.norm(tgt - P, axis=1); m = sep <= np.quantile(sep, 1 - 0.10)
            P, tgt, keep = P[m], tgt[m], keep[m]
        if len(P) > keep_target:               # uniform thinning to the fitted count; coverage unchanged
            idx = np.random.default_rng(seed_b + 3).choice(len(P), keep_target, replace=False)
            P, tgt, keep = P[idx], tgt[idx], keep[idx]
        return apply(np.linalg.inv(M0), P), tgt, P, dict(
            rule=("hybrid targets WITH the separation trim (cell B)" if mode == "hybrid_trimmed" else
                  "hybrid: normal shooting where the shot is kept, nearest point where it is dropped"),
            sampled=int(n), survived=int(len(P)), kept=int(len(P)),
            from_shot=int(keep.sum()), from_nearest=int((~keep).sum()))
    from ihm.anatomy.normal_shooting import shoot_pairs   # modes "shooting" and "nearest_subset"
    src_V = apply(M0, Va)                      # shoot from the starting map's image, in the target's space
    r = shoot_pairs(src_V, Fa, Vd, Fa, n=3 * n, cap_m=SHOOT_CAP_M, return_tol_m=SHOOT_RETURN_TOL_M,
                    seed=seed_a, min_normal_agreement=SHOOT_AGREEMENT)
    P, Q = r["source"], r["target"]
    survived = int(len(P))                     # BEFORE the subsample: the rule's own keep rate
    if len(P) > keep_target:                   # match the other arm's fitted count, so the fit is comparable
        idx = np.random.default_rng(seed_b).choice(len(P), keep_target, replace=False); P, Q = P[idx], Q[idx]
    a = apply(np.linalg.inv(M0), P)
    census = dict(rule=f"normal shooting, return test {1e3*SHOOT_RETURN_TOL_M:.1f} mm, cap {1e3*SHOOT_CAP_M:.0f} mm, "
                       f"min_normal_agreement {SHOOT_AGREEMENT}", sampled=int(r["sampled"]),
                  survived=survived, kept=int(len(P)),
                  no_hit=int(r["no_hit"]), no_return=int(r["no_return"]),
                  return_too_far=int(r["return_too_far"]), normal_disagreed=int(r["normal_disagreed"]))
    if mode == "shooting": return a, Q, P, census
    # Control C (0213dfc): the SAME points and the SAME starting points, but the nearest rule's
    # targets -- correspondence rule held constant against the full-coverage arm, only coverage varies.
    tgt, _, _ = projector(Vd, Fa, np.random.default_rng(seed_b + 7))(P)
    census["rule"] = "nearest point, restricted to the population normal shooting keeps (Control C)"
    return a, tgt, P, census

def stage_recovery_separation(mode="nearest"):
    """The recovery control at the separation this pipeline actually operates at (73dfc35).

    POST-HOC: this replaces a stop rule that was written at 20 mm, the amplitude of the skin's
    displacement, and fired on it.  The correspondence never sees 20 mm.  In the real pipeline the
    target is the nearest point on the scaffold FROM THE PER-SEGMENT SIMILARITY'S IMAGE, and those
    separations measure 0.78 mm (radius) to 6.88 mm (torso).  So the starting map here is the
    best-fit similarity from the atlas bone to the displaced bone -- what M_seg is -- and what is
    scored is the few millimetres left after it.  The 20 mm threshold is NOT carried across at any
    weight: a threshold set at one separation says nothing at another, which is the error being
    corrected.  Pointwise and to-surface are reported separately, because a random spline field
    slides each surface along itself and no surface method can see that component."""
    t0 = time.time(); OUT.mkdir(parents=True, exist_ok=True)
    # Known answer for the decomposition itself (ae34bd2): for an ISOTROPIC field the tangential
    # share is pi/4 by mean and sqrt(2/3) by rms, and the normal share is 1/2.  The measured
    # fraction is reported against this, so a field that is not isotropic declares itself.
    rr = np.random.default_rng(7).normal(size=(400000, 3)); nn = np.array([0., 0., 1.])
    vn = np.abs(rr @ nn); vt = np.linalg.norm(rr - (rr @ nn)[:, None] * nn, axis=1); vm = np.linalg.norm(rr, axis=1)
    say(f"decomposition known answer, isotropic: mean|v_t|/mean|v| {vt.mean()/vm.mean():.4f} (pi/4 {np.pi/4:.4f}), "
        f"rms {np.sqrt((vt**2).mean())/np.sqrt((vm**2).mean()):.4f} (sqrt(2/3) {np.sqrt(2/3):.4f}), "
        f"mean|v_n|/mean|v| {vn.mean()/vm.mean():.4f} (0.5)")
    meshes = {seg: atlas_mesh(seg) for seg in segments}
    allV, allF = atlas_bone_surface()
    centres = sample(allV, allF, FIELD_CENTRES, np.random.default_rng(102))
    w0 = np.random.default_rng(101).normal(size=(len(centres), 3))
    P = np.hstack([centres, np.ones((len(centres), 1))]); Q, _ = np.linalg.qr(P)
    w0 -= Q @ (Q.T @ w0)
    probe = sample(allV, allF, 20000, np.random.default_rng(103))
    unit = float(np.linalg.norm(Warp(np.eye(4), centres, w0, np.zeros((4, 3))).displacement(probe), axis=1).mean())
    # calibration: what a unit field leaves after the starting map, so the amplitudes can be set to
    # land on the pipeline's separations.  Linear in the amplitude, so one pass fixes all three.
    unit_field = Warp(np.eye(4), centres, w0 / unit, np.zeros((4, 3)))
    left = []
    for seg in segments:
        Va, Fa = meshes[seg]; Vd = unit_field.apply(Va)
        s_, R_, t_ = bind.umeyama(Va, Vd)
        a = sample(Va, Fa, 200, np.random.default_rng(50000 + segments.index(seg)))
        left.append(np.linalg.norm(unit_field.apply(a) - apply(sim_matrix(s_, R_, t_), a), axis=1))
    per_metre = float(np.median(np.concatenate(left)))
    say(f"calibration: a field of {FIELD_CENTRES} centres leaves {1e3*per_metre:.3f} mm after the starting map "
        f"per 1 m of mean displacement  [{time.time()-t0:.0f}s]")
    out = {}
    for target in SEPARATION_TARGETS:
        magnitude = target / per_metre
        truth = Warp(np.eye(4), centres, w0 * (magnitude / unit), np.zeros((4, 3)))
        moved = float(np.linalg.norm(truth.displacement(probe), axis=1).mean())
        say(f"\n-- aiming at {1e3*target:.1f} mm of separation")
        say(f"\n== known field moves this body's bones {1e3*moved:.2f} mm on average  [{time.time()-t0:.0f}s]")
        S, T, SEP, TERR, RESID, TOTAL, ES, TRUTH_ES, seg_es, per_seg = [], [], [], [], [], [], [], [], [], {}
        TN, TT, LAB = [], [], []   # the error vector split at each point: normal, tangential (the floor), segment
        CENSUS, ON = [], []        # what the correspondence dropped and why; targets' distance to the surface
        for seg in segments:
            i = segments.index(seg); Va, Fa = meshes[seg]; Vd = truth.apply(Va)
            s_, R_, t_ = bind.umeyama(Va, Vd)                 # the starting map: what M_seg is
            M0 = sim_matrix(s_, R_, t_)
            n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
            a, tgt, start, census = correspond(mode, Va, Fa, Vd, M0, n, 10000 + i, 20000 + i,
                                               int(round(n * CONTROL_FITTED_FRACTION)))
            CENSUS.append((seg, census))
            true_a = truth.apply(a)
            # The target's own face on the displaced surface, for the decomposition.  A target that
            # is not ON that surface declares itself here as a non-zero distance.
            _, d_on, _, face = projector(Vd, Fa, np.random.default_rng(60000 + i), return_face=True)(tgt)
            ON.append(d_on)
            # Split the error against the surface it is measured on.  The tangential part is
            # invisible to ANY surface-based correspondence -- a floor, not a defect of this one --
            # and only the normal part may be called target error.
            tri = Vd[Fa[face]]
            nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
            err = tgt - true_a
            en = np.einsum('ij,ij->i', err, nrm)
            TN.append(np.abs(en)); TT.append(np.linalg.norm(err - en[:, None] * nrm, axis=1))
            LAB += [seg] * len(a)
            sep = np.linalg.norm(tgt - start, axis=1)
            S.append(a); T.append(tgt); SEP.append(sep)
            TERR.append(np.linalg.norm(err, axis=1))
            RESID.append(np.linalg.norm(true_a - start, axis=1))
            TOTAL.append(np.linalg.norm(true_a - a, axis=1))
            per_seg[seg] = float(np.median(sep))
            e = sample(Va, Fa, 2000, np.random.default_rng(30000 + i))
            ES.append(e); TRUTH_ES.append(truth.apply(e)); seg_es += [seg] * len(e)
        S, T = np.concatenate(S), np.concatenate(T); SEP = np.concatenate(SEP)
        TERR, RESID, TOTAL = np.concatenate(TERR), np.concatenate(RESID), np.concatenate(TOTAL)
        ES, TRUTH_ES, seg_es = np.concatenate(ES), np.concatenate(TRUTH_ES), np.asarray(seg_es)
        ON = np.concatenate(ON)
        sampled = sum(c["sampled"] for _, c in CENSUS); kept = sum(c["kept"] for _, c in CENSUS)
        say(f"   correspondence: {CENSUS[0][1]['rule']}")
        if mode == "hybrid":
            shot = sum(c["from_shot"] for _, c in CENSUS); near_n = sum(c["from_nearest"] for _, c in CENSUS)
            say(f"   fitted {kept:,} points, coverage 100% of quota by construction; targets from the shot "
                f"{shot:,} ({100*shot/max(shot+near_n,1):.1f}%), from the nearest rule {near_n:,}")
            say("   per-segment coverage (fitted / quota) and where each segment's targets came from:")
            for s, c in CENSUS:
                q = int(round((N_CORR_TORSO if s == "torso" else N_CORR_PELVIS if s == "pelvis" else N_CORR) * CONTROL_FITTED_FRACTION))
                say(f"      {s:10s} {c['kept']:4d} / {q:4d}  {100*c['kept']/q:5.1f}%   shot {c['from_shot']:4d}, nearest {c['from_nearest']:4d}")
        say(f"   kept {kept:,} of {sampled:,} sampled ({100*kept/sampled:.1f}%)"
            + ("" if mode == "nearest" else
               "; dropped: no hit " + str(sum(c.get("no_hit", 0) for _, c in CENSUS))
               + ", no return " + str(sum(c.get("no_return", 0) for _, c in CENSUS))
               + ", return too far " + str(sum(c.get("return_too_far", 0) for _, c in CENSUS))
               + ", normal disagreed " + str(sum(c.get("normal_disagreed", 0) for _, c in CENSUS))))
        say(f"   targets lie on the surface to {1e3*ON.max():.2e} mm (max)")
        say(f"   separation after the starting map: median {1e3*np.median(SEP):.2f} mm, "
            f"range over segments {1e3*min(per_seg.values()):.2f}-{1e3*max(per_seg.values()):.2f} mm "
            f"(the pipeline's own: 0.78-6.88)")
        say(f"   residual displacement the correspondence must supply: median {1e3*np.median(RESID):.2f} mm "
            f"(total field displacement {1e3*np.median(TOTAL):.2f} mm)")
        TN, TT, LAB = np.concatenate(TN), np.concatenate(TT), np.asarray(LAB)
        frac = float(TT.mean() / TERR.mean())
        # Thin-sheet hits are a DIFFERENT phenomenon from the separation scaling: the nearest point
        # sits on the far wall of a rib or scapula, so the error is a sheet thickness and does not
        # scale with separation.  Reported as a count of affected correspondences, never as a max.
        thin = TN > max(0.010, 5 * float(np.median(RESID)))
        thin_by_seg = {s: int((thin & (LAB == s)).sum()) for s in segments if (thin & (LAB == s)).any()}
        say(f"   thin-sheet hits (normal error > {1e3*max(0.010, 5*float(np.median(RESID))):.1f} mm, the far wall of a "
            f"sheet): {int(thin.sum())} of {len(TN)} correspondences ({100*thin.mean():.2f}%) "
            + (", ".join(f"{s} {n}" for s, n in sorted(thin_by_seg.items(), key=lambda kv: -kv[1])[:6]) or "none"))
        say(f"   error vs the known truth, whole vector: mean {1e3*TERR.mean():.3f} mm "
            f"({100*TERR.mean()/RESID.mean():.1f}% of the residual displacement) -- NOT the correspondence's error")
        say(f"      tangential (the floor, invisible to any surface method): mean {1e3*TT.mean():.3f} mm, "
            f"fraction {frac:.4f}  (isotropic prediction pi/4 = 0.7854)")
        say(f"      NORMAL, the only part a correspondence can be blamed for: mean {1e3*TN.mean():.3f} mm, "
            f"p90 {1e3*np.quantile(TN,.9):.3f}, max {1e3*TN.max():.3f}  "
            f"=  {100*TN.mean()/RESID.mean():.1f}% of the residual displacement")
        lam, cv = choose_lambda(S, T)
        w, a_ = tps_solve(tps_system(S), T - S, lam)
        W = Warp(np.eye(4), S, w, a_)
        fit_res = np.linalg.norm(W.apply(S) - T, axis=1)
        got = W.apply(ES); point = np.linalg.norm(got - TRUTH_ES, axis=1); surf = np.zeros(len(ES))
        for seg in segments:
            m = seg_es == seg; Va, Fa = meshes[seg]
            surf[m] = projector(truth.apply(Va), Fa, np.random.default_rng(40000 + segments.index(seg)))(got[m])[1]
        say(f"   lambda {lam:.2e}; residual at the correspondences {1e3*np.sqrt((fit_res**2).mean()):.3f} mm RMS")
        say(f"   RECOVERY pointwise vs truth: {1e3*np.sqrt((point**2).mean()):.3f} mm RMS, p90 {1e3*np.quantile(point,.9):.3f}"
            f"   (includes the tangential slide no surface method can see)")
        say(f"   RECOVERY to the surface:     {1e3*np.sqrt((surf**2).mean()):.3f} mm RMS, p90 {1e3*np.quantile(surf,.9):.3f}"
            f"   [{time.time()-t0:.0f}s]")
        # DIAGNOSTIC, no threshold attaches and nothing is rescored on it (d4ce65e): agreement with
        # the targets is not a small version of distance from the truth.  The retired control put
        # this at 3.4x (0.149 -> 0.510 mm) and 7.7x (0.437 -> 3.345); what it reads at the
        # operating point is what licenses or forbids quoting a correspondence residual as accuracy.
        understatement = float(np.sqrt((surf ** 2).mean()) / max(np.sqrt((fit_res ** 2).mean()), 1e-30))
        say(f"   residual at correspondences {1e3*np.sqrt((fit_res**2).mean()):.3f} mm vs recovery to surface "
            f"{1e3*np.sqrt((surf**2).mean()):.3f} mm  ->  understatement {understatement:.1f}x  (diagnostic only)")
        out[f"{1e3*moved:.1f}mm"] = dict(
            field_mean_displacement_m=moved, separation_median_m=float(np.median(SEP)),
            separation_by_segment_m=per_seg, residual_displacement_median_m=float(np.median(RESID)),
            total_displacement_median_m=float(np.median(TOTAL)),
            whole_error_mean_m=float(TERR.mean()), whole_error_p90_m=float(np.quantile(TERR, .9)),
            whole_error_max_m=float(TERR.max()),
            whole_error_over_residual=float(TERR.mean() / RESID.mean()),
            tangential_floor_mean_m=float(TT.mean()), tangential_fraction=frac, isotropic_fraction=float(np.pi / 4),
            target_error_normal_mean_m=float(TN.mean()), target_error_normal_p90_m=float(np.quantile(TN, .9)),
            target_error_normal_max_m=float(TN.max()),
            target_error_normal_over_residual=float(TN.mean() / RESID.mean()), lam=lam,
            residual_at_correspondences_rms_m=float(np.sqrt((fit_res ** 2).mean())),
            recovery_pointwise_rms_m=float(np.sqrt((point ** 2).mean())),
            recovery_to_surface_rms_m=float(np.sqrt((surf ** 2).mean())),
            recovery_to_surface_p90_m=float(np.quantile(surf, .9)),
            understatement_to_surface_over_residual=understatement,
            correspondence=CENSUS[0][1]["rule"], sampled=sampled, kept=kept,
            kept_fraction=float(kept / sampled), census_by_segment={s: c for s, c in CENSUS},
            targets_on_surface_max_m=float(ON.max()),
            thin_sheet_hits=int(thin.sum()), thin_sheet_by_segment=thin_by_seg)
    name = {"nearest": "recovery_separation.json", "shooting": "recovery_separation_shooting.json",
            "nearest_subset": "recovery_separation_controlC.json", "hybrid": "recovery_separation_hybrid.json"}[mode]
    (OUT / name).write_text(json.dumps(dict(
        schema="ihm.skin-warp-recovery-separation.v1", status="post-hoc replacement of the 20 mm stop rule (73dfc35)",
        correspondence_mode=mode, carried_over_threshold=None,
        shooting_settings=(None if mode == "nearest" else
                           dict(cap_m=SHOOT_CAP_M, cap_basis="4x the largest measured per-segment separation (12.32 mm)",
                                min_normal_agreement=SHOOT_AGREEMENT, return_tol_m=SHOOT_RETURN_TOL_M,
                                module="ihm/anatomy/normal_shooting.py")),
        scales=out), indent=2) + "\n")

def known_field_setup():
    """Rebuild the control's known field and its calibration, with the same seeds and constants."""
    meshes = {seg: atlas_mesh(seg) for seg in segments}
    allV, allF = atlas_bone_surface()
    centres = sample(allV, allF, FIELD_CENTRES, np.random.default_rng(102))
    w0 = np.random.default_rng(101).normal(size=(len(centres), 3))
    P = np.hstack([centres, np.ones((len(centres), 1))]); Q, _ = np.linalg.qr(P)
    w0 -= Q @ (Q.T @ w0)
    probe = sample(allV, allF, 20000, np.random.default_rng(103))
    unit = float(np.linalg.norm(Warp(np.eye(4), centres, w0, np.zeros((4, 3))).displacement(probe), axis=1).mean())
    unit_field = Warp(np.eye(4), centres, w0 / unit, np.zeros((4, 3)))
    left = []
    for seg in segments:
        Va, Fa = meshes[seg]; Vd = unit_field.apply(Va)
        s_, R_, t_ = bind.umeyama(Va, Vd)
        a = sample(Va, Fa, 200, np.random.default_rng(50000 + segments.index(seg)))
        left.append(np.linalg.norm(unit_field.apply(a) - apply(sim_matrix(s_, R_, t_), a), axis=1))
    return meshes, centres, w0, unit, float(np.median(np.concatenate(left)))

def stage_subset():
    """Score the NEAREST rule on exactly the points normal shooting kept (ded65bf).

    1.933 mm against 0.444 mm compares two rules on two populations: shooting declines to answer
    for the far-wall hits that inflate the nearest rule's mean, so its kept set may be enriched for
    points the nearest rule also handles well.  Here both rules are scored on the SAME points.
    The reconstruction must reproduce the shooting arm's own published number or it is void."""
    t0 = time.time(); meshes, centres, w0, unit, per_metre = known_field_setup()
    out = {}
    for aim, label in ((0.003, "operating amplitude"), (0.007, "largest amplitude")):
        truth = Warp(np.eye(4), centres, w0 * ((aim / per_metre) / unit), np.zeros((4, 3)))
        near_n, shoot_n, per_seg, tot_s, tot_k = [], [], {}, 0, 0
        for seg in segments:
            i = segments.index(seg); Va, Fa = meshes[seg]; Vd = truth.apply(Va)
            s_, R_, t_ = bind.umeyama(Va, Vd); M0 = sim_matrix(s_, R_, t_)
            n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
            a, tgt_shoot, start, census = correspond("shooting", Va, Fa, Vd, M0, n, 10000 + i, 20000 + i,
                                                     int(round(n * CONTROL_FITTED_FRACTION)))
            # the rule's OWN keep rate, before the subsample that matches the other arm's count
            survived = census.get("survived", census["kept"])
            per_seg[seg] = (survived, census["sampled"]); tot_k += survived; tot_s += census["sampled"]
            if not len(a): continue
            # the NEAREST rule, from the same starting points, on the same population
            tgt_near, _, _ = projector(Vd, Fa, np.random.default_rng(20000 + i))(start)
            true_a = truth.apply(a)
            _, _, _, face = projector(Vd, Fa, np.random.default_rng(60000 + i), return_face=True)(tgt_shoot)
            tri = Vd[Fa[face]]
            nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
            for tgt, acc in ((tgt_near, near_n), (tgt_shoot, shoot_n)):
                err = tgt - true_a
                acc.append(np.abs(np.einsum('ij,ij->i', err, nrm)))
        near_n, shoot_n = np.concatenate(near_n), np.concatenate(shoot_n)
        say(f"\n== {label} (aim {1e3*aim:.0f} mm): {tot_k:,} of {tot_s:,} kept ({100*tot_k/tot_s:.1f}%)  [{time.time()-t0:.0f}s]")
        say(f"   on the SAME {len(near_n):,} points, normal error: nearest {1e3*near_n.mean():.3f} mm, "
            f"shooting {1e3*shoot_n.mean():.3f} mm  ->  {near_n.mean()/max(shoot_n.mean(),1e-30):.2f}x")
        say(f"   (the shooting number must match this amplitude's published one, or the subset is not the same set)")
        say(f"   kept fraction per segment:")
        for s in sorted(per_seg, key=lambda k: per_seg[k][0] / max(per_seg[k][1], 1)):
            k_, s_ = per_seg[s]
            say(f"      {s:10s} {k_:5d} / {s_:5d}  {100*k_/max(s_,1):5.1f}%")
        out[label] = dict(aim_m=aim, kept=tot_k, sampled=tot_s, kept_fraction=float(tot_k / tot_s),
                          points=int(len(near_n)), nearest_normal_mean_m=float(near_n.mean()),
                          shooting_normal_mean_m=float(shoot_n.mean()),
                          ratio=float(near_n.mean() / max(shoot_n.mean(), 1e-30)),
                          kept_by_segment={s: dict(kept=k_, sampled=s_, fraction=k_ / max(s_, 1)) for s, (k_, s_) in per_seg.items()})
    (OUT / "subset_control.json").write_text(json.dumps(dict(
        schema="ihm.skin-warp-subset-control.v1",
        basis="both correspondence rules scored on the population normal shooting keeps, from the same starting points",
        scales=out), indent=2) + "\n")

TRIM_SWEEP = (0.0, 0.02, 0.05, 0.10)

def stage_trim_sweep():
    """What the borrowed 10% costs, swept (b034206).

    The trim entered the v1 pre-registration by ANALOGY with the per-segment fit's ICP trim, never
    measured.  Here it is varied with everything else fixed: nearest targets, amplitude 2, lambda
    1e-3, and the fitted COUNT held constant at the trimmed arm's, so only WHICH points are kept
    changes, not how many.

    CAVEAT, which limits what this can conclude: this control's truth is a known field applied to
    this same body, so its largest separations are hard but GENUINE.  On a real subject-to-scaffold
    registration the largest separations may be WRONG correspondences rather than merely hard ones,
    and a trim that costs accuracy here could be protective there.  No ground truth exists on real
    data -- which is why this control exists -- so the only evidence available says remove the trim
    and no evidence says keep it.  That is weaker than the numbers look."""
    t0 = time.time(); meshes, centres, w0, unit, per_metre = known_field_setup()
    aim, LAM = 0.003, 1e-3
    truth = Warp(np.eye(4), centres, w0 * ((aim / per_metre) / unit), np.zeros((4, 3)))
    rows = []
    for trim in TRIM_SWEEP:
        S, T, TN, ES, TRUTH_ES, seg_es = [], [], [], [], [], []
        for seg in segments:
            i = segments.index(seg); Va, Fa = meshes[seg]; Vd = truth.apply(Va)
            s_, R_, t_ = bind.umeyama(Va, Vd); M0 = sim_matrix(s_, R_, t_)
            n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
            quota = int(round(n * CONTROL_FITTED_FRACTION))
            a = sample(Va, Fa, n, np.random.default_rng(10000 + i)); start = apply(M0, a)
            tgt, d, _ = projector(Vd, Fa, np.random.default_rng(20000 + i))(start)
            if trim > 0:
                m = d <= np.quantile(d, 1 - trim); a, tgt, start = a[m], tgt[m], start[m]
            if len(a) > quota:                       # count held constant across the sweep
                idx = np.random.default_rng(20005 + i).choice(len(a), quota, replace=False)
                a, tgt, start = a[idx], tgt[idx], start[idx]
            _, _, _, face = projector(Vd, Fa, np.random.default_rng(60000 + i), return_face=True)(tgt)
            tri = Vd[Fa[face]]
            nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
            err = tgt - truth.apply(a)
            TN.append(np.abs(np.einsum('ij,ij->i', err, nrm)))
            S.append(a); T.append(tgt)
            e = sample(Va, Fa, 2000, np.random.default_rng(30000 + i))
            ES.append(e); TRUTH_ES.append(truth.apply(e)); seg_es += [seg] * len(e)
        S, T, TN = np.concatenate(S), np.concatenate(T), np.concatenate(TN)
        ES, TRUTH_ES, seg_es = np.concatenate(ES), np.concatenate(TRUTH_ES), np.asarray(seg_es)
        w, a_ = tps_solve(tps_system(S), T - S, LAM)
        got = Warp(np.eye(4), S, w, a_).apply(ES); surf = np.zeros(len(ES))
        for seg in segments:
            m = seg_es == seg; Va, Fa = meshes[seg]
            surf[m] = projector(truth.apply(Va), Fa, np.random.default_rng(40000 + segments.index(seg)))(got[m])[1]
        rec = float(np.sqrt((surf ** 2).mean()))
        rows.append(dict(trim=trim, points=int(len(S)), to_surface_rms_m=rec,
                         target_error_normal_mean_m=float(TN.mean())))
        say(f"   trim {100*trim:5.1f}%  points {len(S):,}  to-surface {1e3*rec:.3f} mm  "
            f"normal target error {1e3*TN.mean():.3f} mm  [{time.time()-t0:.0f}s]")
    best = min(rows, key=lambda r: r["to_surface_rms_m"])
    monotone = all(rows[i]["to_surface_rms_m"] <= rows[i + 1]["to_surface_rms_m"] for i in range(len(rows) - 1))
    say(f"\nbest at trim {100*best['trim']:.0f}%; monotone in the trim: {monotone}")
    say("caveat: this control's largest separations are genuine; on real data they may be wrong "
        "correspondences instead, and no ground truth exists there to tell.")
    (OUT / "trim_sweep.json").write_text(json.dumps(dict(schema="ihm.skin-warp-trim-sweep.v1", amplitude_aim_m=aim,
        lam=LAM, rows=rows, monotone=monotone, best_trim=best["trim"],
        caveat="the control's largest separations are hard but genuine; on real data they may be wrong "
               "correspondences, where a trim could be protective. No ground truth exists on real data."), indent=2) + "\n")

def stage_twobytwo():
    """The 2x2 that separates the trim from the targets (269f1bb), at the operating amplitude.

                        trimmed 10%   untrimmed
      nearest targets   2.138 (have)  cell A
      hybrid targets    cell B        1.686 (have)

    Within a row only the trim changes, so each row isolates it exactly.  Across a column the
    sampler also differs (surface_samples vs the area-weighted sampler), which is stated rather
    than controlled.  lambda is held at 1e-3, what both existing arms' own rule chose here."""
    t0 = time.time(); meshes, centres, w0, unit, per_metre = known_field_setup()
    aim, LAM = 0.003, 1e-3
    truth = Warp(np.eye(4), centres, w0 * ((aim / per_metre) / unit), np.zeros((4, 3)))
    out = {}
    for mode, label in (("nearest_untrimmed", "cell A: nearest targets, untrimmed"),
                        ("hybrid_trimmed", "cell B: hybrid targets, trimmed")):
        S, T, TN, ES, TRUTH_ES, seg_es = [], [], [], [], [], []
        for seg in segments:
            i = segments.index(seg); Va, Fa = meshes[seg]; Vd = truth.apply(Va)
            s_, R_, t_ = bind.umeyama(Va, Vd); M0 = sim_matrix(s_, R_, t_)
            n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
            a, tgt, start, _ = correspond(mode, Va, Fa, Vd, M0, n, 10000 + i, 20000 + i,
                                          int(round(n * CONTROL_FITTED_FRACTION)))
            _, _, _, face = projector(Vd, Fa, np.random.default_rng(60000 + i), return_face=True)(tgt)
            tri = Vd[Fa[face]]
            nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
            err = tgt - truth.apply(a)
            TN.append(np.abs(np.einsum('ij,ij->i', err, nrm)))
            S.append(a); T.append(tgt)
            e = sample(Va, Fa, 2000, np.random.default_rng(30000 + i))
            ES.append(e); TRUTH_ES.append(truth.apply(e)); seg_es += [seg] * len(e)
        S, T, TN = np.concatenate(S), np.concatenate(T), np.concatenate(TN)
        ES, TRUTH_ES, seg_es = np.concatenate(ES), np.concatenate(TRUTH_ES), np.asarray(seg_es)
        w, a_ = tps_solve(tps_system(S), T - S, LAM)
        got = Warp(np.eye(4), S, w, a_).apply(ES); surf = np.zeros(len(ES))
        for seg in segments:
            m = seg_es == seg; Va, Fa = meshes[seg]
            surf[m] = projector(truth.apply(Va), Fa, np.random.default_rng(40000 + segments.index(seg)))(got[m])[1]
        rec = float(np.sqrt((surf ** 2).mean()))
        out[mode] = dict(label=label, to_surface_rms_m=rec, target_error_normal_mean_m=float(TN.mean()),
                         target_error_normal_p90_m=float(np.quantile(TN, .9)), points=int(len(S)))
        say(f"{label}: to-surface {1e3*rec:.3f} mm | normal target error mean {1e3*TN.mean():.3f} mm, "
            f"p90 {1e3*np.quantile(TN,.9):.3f} mm  [{time.time()-t0:.0f}s]")
    say(f"\n{'':22s} {'trimmed 10%':>13s} {'untrimmed':>13s}")
    say(f"{'nearest targets':22s} {'2.138':>13s} {1e3*out['nearest_untrimmed']['to_surface_rms_m']:13.3f}")
    say(f"{'hybrid targets':22s} {1e3*out['hybrid_trimmed']['to_surface_rms_m']:13.3f} {'1.686':>13s}")
    say("(rows isolate the trim exactly; columns also differ in sampler, which is stated not controlled)")
    say(f"nearest targets, trim effect: 2.138 -> {1e3*out['nearest_untrimmed']['to_surface_rms_m']:.3f} mm")
    say(f"hybrid targets,  trim effect: 1.686 -> {1e3*out['hybrid_trimmed']['to_surface_rms_m']:.3f} mm")
    say(f"cell A target error {1e3*out['nearest_untrimmed']['target_error_normal_mean_m']:.3f} mm against the "
        f"trimmed arm's 1.933 mm -- if far above, the trim was doing quality work as well as coverage work")
    (OUT / "twobytwo.json").write_text(json.dumps(dict(schema="ihm.skin-warp-2x2.v1", amplitude_aim_m=aim, lam=LAM,
        have=dict(nearest_trimmed_to_surface_m=0.002138, hybrid_untrimmed_to_surface_m=0.001686),
        cells=out, caveat="rows isolate the trim; columns also differ in sampler"), indent=2) + "\n")

def stage_per_segment():
    """Per-segment to-surface recovery, hybrid against full-coverage nearest (80afd28).

    Tests the mechanism offered for the hybrid's win: that it improves the bulk of points where the
    shot survives, while the torso -- 6.8% shot share -- uses nearest targets anyway and should
    therefore barely move.  If the torso improves as much as the well-covered segments, the
    mechanism is wrong.  Both arms are fitted at lambda 1e-3, which is what each arm's own rule
    chose at this amplitude, so no cross-validation is repeated and only the targets differ."""
    t0 = time.time(); meshes, centres, w0, unit, per_metre = known_field_setup()
    aim, LAM = 0.003, 1e-3
    truth = Warp(np.eye(4), centres, w0 * ((aim / per_metre) / unit), np.zeros((4, 3)))
    share, out = {}, {}
    for mode in ("nearest", "hybrid"):
        S, T, ES, TRUTH_ES, seg_es = [], [], [], [], []
        for seg in segments:
            i = segments.index(seg); Va, Fa = meshes[seg]; Vd = truth.apply(Va)
            s_, R_, t_ = bind.umeyama(Va, Vd); M0 = sim_matrix(s_, R_, t_)
            n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
            a, tgt, start, census = correspond(mode, Va, Fa, Vd, M0, n, 10000 + i, 20000 + i,
                                               int(round(n * CONTROL_FITTED_FRACTION)))
            if mode == "hybrid": share[seg] = census["from_shot"] / max(census["kept"], 1)
            S.append(a); T.append(tgt)
            e = sample(Va, Fa, 2000, np.random.default_rng(30000 + i))
            ES.append(e); TRUTH_ES.append(truth.apply(e)); seg_es += [seg] * len(e)
        S, T = np.concatenate(S), np.concatenate(T)
        ES, TRUTH_ES, seg_es = np.concatenate(ES), np.concatenate(TRUTH_ES), np.asarray(seg_es)
        w, a_ = tps_solve(tps_system(S), T - S, LAM)
        got = Warp(np.eye(4), S, w, a_).apply(ES)
        rms = {}
        for seg in segments:
            m = seg_es == seg; Va, Fa = meshes[seg]
            d = projector(truth.apply(Va), Fa, np.random.default_rng(40000 + segments.index(seg)))(got[m])[1]
            rms[seg] = float(np.sqrt((d ** 2).mean()))
        out[mode] = rms
        say(f"   {mode} fitted and scored  [{time.time()-t0:.0f}s]")
    say(f"\n{'segment':10s} {'shot share':>11s} {'nearest':>9s} {'hybrid':>9s} {'change':>9s}")
    rows = []
    for seg in sorted(segments, key=lambda s: -share[s]):
        n_, h_ = out["nearest"][seg], out["hybrid"][seg]
        rows.append((seg, share[seg], n_, h_))
        say(f"{seg:10s} {100*share[seg]:10.1f}% {1e3*n_:9.3f} {1e3*h_:9.3f} {100*(h_-n_)/n_:8.1f}%")
    hi = [r for r in rows if r[1] >= 0.5]; lo = [r for r in rows if r[1] < 0.5]
    for label, group in (("shot share >= 50%", hi), ("shot share < 50%", lo)):
        if group:
            say(f"{label}: {len(group)} segments, mean change "
                f"{100*np.mean([(h-n)/n for _, _, n, h in group]):.1f}%")
    (OUT / "per_segment_recovery.json").write_text(json.dumps(dict(
        schema="ihm.skin-warp-per-segment.v1", amplitude_aim_m=aim, lam=LAM,
        basis="hybrid vs full-coverage nearest, identical fit path, lambda held at what both arms' own rule chose",
        rows=[dict(segment=s, shot_share=sh, nearest_to_surface_m=n_, hybrid_to_surface_m=h_,
                   change=(h_ - n_) / n_) for s, sh, n_, h_ in rows]), indent=2) + "\n")

def stage_toe_geometry():
    """Does the toe skin carry more material than its warped bone frame can hold? (8f78bf2)

    The mechanism under test: if the scaffold's forefoot differs from this specimen's by ~21%, skin
    carried rigidly through the warp has the wrong amount of material for the frame, the excess
    buckles (gate 4's folds) and the bone pushes through where the skin pulled away (gate 2's
    shortfall).  It predicts compression NEAR THE SCALE FACTOR and folds where local area
    compression is greatest.  Note the direction: this warp maps the atlas onto the SCAFFOLD, whose
    foot is the larger, so the segment totals are expected to EXPAND; only a local measure can test
    the mechanism, and that is what is computed here."""
    W = load_warp(OUT / "warp.npz")
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    used = np.unique(F[ext])
    Y0, Y1 = np.zeros_like(V), np.zeros_like(V)
    Y0[used], Y1[used] = W.ground(V[used]), W.apply(V[used])
    def nrm(Y, f): t = Y[f]; return np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    n0, n1 = nrm(Y0, F[ext]), nrm(Y1, F[ext])
    a0, a1 = np.linalg.norm(n0, axis=1) / 2, np.linalg.norm(n1, axis=1) / 2
    ratio = np.where(a0 > 0, a1 / np.maximum(a0, 1e-30), np.nan)
    inverted = (n0 * n1).sum(1) < 0
    sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes())); names = [s["id"] for s in sb["segments"]]
    Wt = np.asarray(sb["weights"], np.float32)
    owner = ((Wt[F[:, 0]] + Wt[F[:, 1]] + Wt[F[:, 2]]) / 3).argmax(1)[ext]
    say(f"{'segment':10s} {'area total':>10s} {'median tri':>11s} {'p10 tri':>8s} {'frac<1':>7s} "
        f"{'principal extents (after/before)':>34s}  {'folds':>6s} {'fold median ratio':>18s}")
    rows = []
    for i, name in enumerate(names):
        m = owner == i
        if not m.any(): continue
        vm = np.unique(F[ext][m])
        P0, P1 = Y0[vm], Y1[vm]
        ext0 = np.linalg.svd(P0 - P0.mean(0), compute_uv=False) / max(len(P0) - 1, 1) ** 0.5
        ext1 = np.linalg.svd(P1 - P1.mean(0), compute_uv=False) / max(len(P1) - 1, 1) ** 0.5
        pr = ext1 / np.maximum(ext0, 1e-30)
        f_m = inverted & m
        rows.append(dict(segment=name, area_total=float(a1[m].sum() / a0[m].sum()),
                         median_triangle=float(np.nanmedian(ratio[m])), p10=float(np.nanquantile(ratio[m], .1)),
                         fraction_below_one=float(np.nanmean(ratio[m] < 1)), principal=[float(x) for x in pr],
                         folds=int(f_m.sum()),
                         fold_median_ratio=(float(np.nanmedian(ratio[f_m])) if f_m.any() else None)))
        r = rows[-1]
        say(f"{name:10s} {r['area_total']:10.3f} {r['median_triangle']:11.3f} {r['p10']:8.3f} "
            f"{r['fraction_below_one']:7.3f} {pr[0]:10.3f} {pr[1]:10.3f} {pr[2]:10.3f}  {r['folds']:6d} "
            + ("      -" if r['fold_median_ratio'] is None else f"{r['fold_median_ratio']:18.3f}"))
    fold_all = ratio[inverted]; rest = ratio[~inverted]
    say(f"\nfolded triangles: {int(inverted.sum())}, median area ratio {np.nanmedian(fold_all):.3f}; "
        f"all others {np.nanmedian(rest):.3f}")
    say(f"folds with local COMPRESSION (ratio < 1): {int(np.nansum(fold_all < 1))} of {int(inverted.sum())}")
    (OUT / "toe_geometry.json").write_text(json.dumps(dict(schema="ihm.skin-warp-toe-geometry.v1",
        basis="per-triangle area change and per-segment principal extents under the fitted warp; "
              "the warp maps atlas -> scaffold, whose foot is larger, so expansion is the expected direction",
        rows=rows, folded_median_ratio=float(np.nanmedian(fold_all)), others_median_ratio=float(np.nanmedian(rest)),
        folds_compressive=int(np.nansum(fold_all < 1)), folds_total=int(inverted.sum())), indent=2) + "\n")

def stage_seam():
    """Do folds sit at the MTP seam for a reason BEYOND being where compression is? (43a63e4)

    A bare 'folds cluster at the seam' is circular if the seam is simply the most compressed place,
    so folded triangles are compared against unfolded ones MATCHED on local area ratio.  If matched
    controls do not exist -- if folds occupy a compression range unfolded triangles never reach --
    that is reported rather than papered over, and it answers the question in the other direction."""
    CALIPER = 0.02
    W = load_warp(OUT / "warp.npz")
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    used = np.unique(F[ext])
    Y0, Y1 = np.zeros_like(V), np.zeros_like(V)
    Y0[used], Y1[used] = W.ground(V[used]), W.apply(V[used])
    tri = F[ext]
    def nrm(Y): t = Y[tri]; return np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    n0, n1 = nrm(Y0), nrm(Y1)
    a0, a1 = np.linalg.norm(n0, axis=1) / 2, np.linalg.norm(n1, axis=1) / 2
    ratio = a1 / np.maximum(a0, 1e-30); inverted = (n0 * n1).sum(1) < 0
    sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes())); names = [s["id"] for s in sb["segments"]]
    Wt = np.asarray(sb["weights"], np.float32)
    owner = ((Wt[F[:, 0]] + Wt[F[:, 1]] + Wt[F[:, 2]]) / 3).argmax(1)[ext]
    out = {}
    for side in ("l", "r"):
        ti, ci = names.index("toes_" + side), names.index("calcn_" + side)
        toes, calcn = owner == ti, owner == ci
        if not toes.any(): continue
        # the seam: vertices carried by BOTH partitions, i.e. the MTP boundary on the skin
        seam = np.intersect1d(np.unique(tri[toes]), np.unique(tri[calcn]))
        if not len(seam): say(f"toes_{side}: no shared boundary with calcn_{side}"); continue
        cen = Y0[tri[toes]].mean(1)
        d = cKDTree(Y0[seam]).query(cen)[0]
        r_, inv_ = ratio[toes], inverted[toes]
        say(f"\n== toes_{side}: {len(cen):,} triangles, {int(inv_.sum())} folded, seam of {len(seam)} vertices")
        say(f"   RAW: folded {1e3*d[inv_].mean():.1f} mm from the seam, unfolded {1e3*d[~inv_].mean():.1f} mm")
        say(f"   compression: folded ratio {np.median(r_[inv_]):.3f} (min {r_[inv_].min():.3f}, max {r_[inv_].max():.3f}); "
            f"unfolded {np.median(r_[~inv_]):.3f}")
        # matched on local area ratio
        lo, hi = r_[inv_].min() - CALIPER, r_[inv_].max() + CALIPER
        pool = (~inv_) & (r_ >= lo) & (r_ <= hi)
        say(f"   unfolded triangles inside the folded compression range [{lo:.3f}, {hi:.3f}]: {int(pool.sum())}")
        pairs, matched_f, matched_u = 0, [], []
        for k in np.flatnonzero(inv_):
            c = (~inv_) & (np.abs(r_ - r_[k]) <= CALIPER)
            if c.any(): pairs += 1; matched_f.append(d[k]); matched_u.append(d[c].mean())
        if pairs:
            mf, mu = np.asarray(matched_f), np.asarray(matched_u)
            diff = mf - mu
            t = float(diff.mean() / (diff.std(ddof=1) / np.sqrt(len(diff)))) if len(diff) > 1 and diff.std(ddof=1) > 0 else float('nan')
            say(f"   MATCHED (caliper {CALIPER} on area ratio): {pairs} of {int(inv_.sum())} folded triangles matched")
            say(f"      folded {1e3*mf.mean():.1f} mm vs equally-compressed unfolded {1e3*mu.mean():.1f} mm "
                f"-> difference {1e3*diff.mean():+.1f} mm, paired t = {t:+.2f}")
        else:
            say(f"   MATCHED: NO unfolded triangle is within {CALIPER} of any folded one's area ratio -- "
                f"the folds occupy a compression range the rest of the toe skin never reaches, so the "
                f"matched test cannot be run and compression alone separates them")
        out[f"toes_{side}"] = dict(triangles=int(len(cen)), folded=int(inv_.sum()), seam_vertices=int(len(seam)),
            raw_folded_mm=1e3 * float(d[inv_].mean()), raw_unfolded_mm=1e3 * float(d[~inv_].mean()),
            folded_ratio_median=float(np.median(r_[inv_])), folded_ratio_min=float(r_[inv_].min()),
            folded_ratio_max=float(r_[inv_].max()), unfolded_ratio_median=float(np.median(r_[~inv_])),
            pool_in_range=int(pool.sum()), matched_pairs=int(pairs),
            matched_folded_mm=(1e3 * float(np.mean(matched_f)) if pairs else None),
            matched_unfolded_mm=(1e3 * float(np.mean(matched_u)) if pairs else None))
    (OUT / "seam_test.json").write_text(json.dumps(dict(schema="ihm.skin-warp-seam.v1", caliper=CALIPER,
        basis="folded vs unfolded toe triangles, matched on local area ratio, compared by distance to the "
              "toes/calcn partition seam (the MTP boundary on the skin)", sides=out), indent=2) + "\n")

def stage_contact_folds():
    """Do the folded triangles reach the contact layer the engine uses? (a173957)

    Option (a) -- accept the folds and gate the body anyway -- is viable only if they are not in the
    contact set, or carry negligible area and no load.  A fold reverses the surface locally, so a
    triangle there has an outward normal pointing inward and any force it carried would push the
    wrong way.  Checked here: whether the bundle the engine loads still contains them, their share
    of the contact area, and where they sit relative to the plantar band that meets the floor."""
    W = load_warp(OUT / "warp.npz")
    Tb, frames_b, _ = bscm.binding_registration()
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    source = W.apply(V)                                   # exactly what the bundle is built from
    used = np.unique(F[ext])
    Y0 = np.zeros_like(V); Y0[used] = W.ground(V[used])
    def nrm(Y, f): t = Y[f]; return np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    inverted = (nrm(Y0, F[ext]) * nrm(source, F[ext])).sum(1) < 0
    sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes())); names = [s["id"] for s in sb["segments"]]
    Wt = np.asarray(sb["weights"], np.float32)
    owner = ((Wt[F[:, 0]] + Wt[F[:, 1]] + Wt[F[:, 2]]) / 3).argmax(1)[ext]
    manifest = json.loads((ROOT / BUNDLE / "manifest.json").read_text())
    rec = {r["body"]: r for r in manifest["records"]}
    say(f"folded triangles: {int(inverted.sum())} of {len(ext):,} exterior")
    out = {}
    for i, name in enumerate(names):
        m = owner == i
        if not (m & inverted).any(): continue
        sel = np.flatnonzero(m)                            # this segment's exterior triangles, build()'s order
        fold_local = inverted[m]
        tri_g = source[F[ext][m]]
        area = np.linalg.norm(np.cross(tri_g[:, 1] - tri_g[:, 0], tri_g[:, 2] - tri_g[:, 0]), axis=1) / 2
        r = rec.get(name)
        in_bundle = r is not None and r["exterior_triangles"] == len(sel)
        # into the segment's own frame, where the bundle stores it and where the floor is met
        M = np.linalg.inv(frames_b[name])
        loc = tri_g.reshape(-1, 3) @ M[:3, :3].T + M[:3, 3]
        cen_y = loc.reshape(-1, 3, 3).mean(1)[:, 1]
        piece_min = float(cen_y.min()) if not len(loc) else float((loc[:, 1]).min())
        h = cen_y[fold_local] - piece_min                  # height above the piece's lowest point
        plantar = int((h <= 0.010).sum())
        say(f"\n{name}: {int(fold_local.sum())} folded of {len(sel):,} triangles in this piece")
        say(f"   still in the bundle the engine loads: {in_bundle} "
            f"(manifest records {r['exterior_triangles'] if r else 'no record'} exterior triangles, build selected {len(sel)})")
        say(f"   folded area {1e6*area[fold_local].sum():.1f} mm2 of {1e6*area.sum():.0f} mm2 in this piece "
            f"({100*area[fold_local].sum()/area.sum():.3f}%); of the bundle's {1e4*r['cut_surface_area_m2']:.0f} cm2 cut area "
            f"that is {100*area[fold_local].sum()/r['cut_surface_area_m2']:.3f}%" if r else "")
        say(f"   height above the piece's lowest point: median {1e3*np.median(h):.1f} mm, min {1e3*h.min():.1f} mm; "
            f"{plantar} of {int(fold_local.sum())} lie within the lowest 10 mm (the band that meets the floor)")
        out[name] = dict(folded=int(fold_local.sum()), piece_triangles=int(len(sel)), in_bundle=bool(in_bundle),
                         folded_area_m2=float(area[fold_local].sum()), piece_area_m2=float(area.sum()),
                         folded_area_fraction=float(area[fold_local].sum() / area.sum()),
                         height_above_lowest_median_m=float(np.median(h)), height_above_lowest_min_m=float(h.min()),
                         within_plantar_10mm=plantar)
    (OUT / "contact_folds.json").write_text(json.dumps(dict(schema="ihm.skin-warp-contact-folds.v1",
        bundle=BUNDLE, basis="folded triangles located in the bundle the engine loads, by area share and by "
        "height above the piece's lowest point in its own segment frame", segments=out), indent=2) + "\n")

def stage_crawl_clearance(trajectory):
    """SUPPLEMENTARY: how close the folded triangles come to the floor over a crawl (3070d00).

    THIS IS NOT THE COMMISSIONED MEASUREMENT.  a0653f5 asked for clearance over a crawl that
    RESPECTS the declared joint limits.  No such trajectory exists in this repo: all 67 stored
    trajectories are outside the declared ranges, and the only two crawls plantarflex the ankle to
    -1.073 rad against a declared +-0.873.  So this is measured against an INADMISSIBLE trajectory
    and it does not close the branch.  Clearance here means the standing-pose conclusion survives a
    pose that over-plantarflexes by 11.5 degrees; it does not mean an admissible crawl was tested."""
    spec = importlib.util.spec_from_file_location("crawl", ROOT / "scripts/crawl.py")
    crawl = importlib.util.module_from_spec(spec); spec.loader.exec_module(crawl)
    ranges = crawl.declared_ranges(B.MODEL)
    t = json.loads((ROOT / trajectory).read_text()); frames = t["frames"]
    def pose_of(f):
        d = f.get("joints") or f.get("coordinates") or {}
        return {k: (v["value"] if isinstance(v, dict) else v) for k, v in d.items()}
    worst, where = 0.0, None; outside = set()
    for f in frames:
        for k, v in pose_of(f).items():
            if k in ranges and isinstance(v, (int, float)):
                lo, hi = ranges[k]; e = max(lo - v, v - hi, 0.0)
                if e > 1e-9: outside.add(k)
                if e > worst: worst, where = e, k
    say(f"trajectory {trajectory}: {len(frames)} frames")
    say(f"AUDIT: outside the declared ranges on {len(outside)} of {len(ranges)} coordinates; worst "
        f"{worst:.4f} rad ({np.degrees(worst):.1f} deg) on {where}  -> INADMISSIBLE, supplementary only")
    W = load_warp(OUT / "warp.npz"); Tb, frames_b, _ = bscm.binding_registration()
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    source = W.apply(V); used = np.unique(F[ext])
    Y0 = np.zeros_like(V); Y0[used] = W.ground(V[used])
    def nrm(Y, f): tt = Y[f]; return np.cross(tt[:, 1] - tt[:, 0], tt[:, 2] - tt[:, 0])
    inverted = (nrm(Y0, F[ext]) * nrm(source, F[ext])).sum(1) < 0
    sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes())); names = [s["id"] for s in sb["segments"]]
    Wt = np.asarray(sb["weights"], np.float32)
    owner = ((Wt[F[:, 0]] + Wt[F[:, 1]] + Wt[F[:, 2]]) / 3).argmax(1)[ext]
    fold_local, piece_local = {}, {}
    for i, name in enumerate(names):
        m = (owner == i) & inverted
        if not m.any(): continue
        M = np.linalg.inv(frames_b[name])
        fold_local[name] = source[np.unique(F[ext][m])] @ M[:3, :3].T + M[:3, 3]
        piece_local[name] = source[np.unique(F[ext][owner == i])] @ M[:3, :3].T + M[:3, 3]
    model = render.OsimModel(B.MODEL)
    best = dict(clearance=np.inf); series = []
    for n, f in enumerate(frames):
        tr = model.forward(pose_of(f))
        lo_fold, lo_piece = np.inf, np.inf
        for name, P in fold_local.items():
            Mt = tr[name]
            lo_fold = min(lo_fold, float((P @ Mt[:3, :3].T + Mt[:3, 3])[:, 1].min()))
            Q = piece_local[name]
            lo_piece = min(lo_piece, float((Q @ Mt[:3, :3].T + Mt[:3, 3])[:, 1].min()))
        series.append(lo_fold)
        if lo_fold < best["clearance"]:
            best = dict(clearance=lo_fold, frame=n, time_s=f.get("time_s"), piece_lowest=lo_piece)
    s = np.asarray(series)
    say(f"\nfolded triangles vs the floor (y = 0) over {len(frames)} frames:")
    say(f"   minimum clearance {1e3*best['clearance']:+.1f} mm at frame {best['frame']} (t = {best['time_s']}s); "
        f"the skin of those segments reaches {1e3*best['piece_lowest']:+.1f} mm at that frame")
    say(f"   median over frames {1e3*np.median(s):+.1f} mm, 5th percentile {1e3*np.quantile(s,.05):+.1f} mm")
    say(f"   frames with a folded triangle at or below the floor: {int((s <= 0).sum())} of {len(s)}")
    (OUT / "crawl_clearance.json").write_text(json.dumps(dict(schema="ihm.skin-warp-crawl-clearance.v1",
        status="SUPPLEMENTARY -- measured against an INADMISSIBLE trajectory; not the measurement "
               "commissioned in a0653f5, and it does not close that branch",
        trajectory=str(trajectory), frames=len(frames),
        audit=dict(coordinates_outside=len(outside), coordinates_total=len(ranges),
                   worst_excursion_rad=worst, worst_coordinate=where),
        minimum_clearance_m=float(best["clearance"]), minimum_frame=int(best["frame"]),
        minimum_time_s=best["time_s"], piece_lowest_at_minimum_m=float(best["piece_lowest"]),
        median_clearance_m=float(np.median(s)), frames_at_or_below_floor=int((s <= 0).sum())), indent=2) + "\n")

def stage_slivers():
    """The triangles gate 4 caught: slivers in the SOURCE mesh, or ordinary triangles?

    det J > 0 everywhere is LOCAL invertibility, at a point.  A triangle of finite size can still
    come out with a reversed normal under a map that is locally orientation-preserving at each of
    its three corners, and the thinner the triangle the smaller the rotation needed.  So the
    question is not whether the flow folded (it cannot) but whether these triangles are degenerate
    to begin with, which their aspect ratio against the whole mesh's answers."""
    W = load_warp(OUT / "warp.npz")
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    used = np.unique(F[ext]); loc = np.searchsorted(used, F[ext])
    Y0, Y1 = W.ground(V[used]), W.apply(V[used])
    def nrm(Y): t = Y[loc]; return np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    n0, n1 = nrm(Y0), nrm(Y1); dots = (n0 * n1).sum(1)
    tri = V[F[ext]]
    edges = np.stack([np.linalg.norm(tri[:, 1] - tri[:, 0], axis=1), np.linalg.norm(tri[:, 2] - tri[:, 1], axis=1),
                      np.linalg.norm(tri[:, 0] - tri[:, 2], axis=1)], axis=1)
    area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    inradius = np.where(area > 0, 2 * area / edges.sum(1), 0.0)
    aspect = np.where(inradius > 0, edges.max(1) / (2 * inradius), np.inf)   # 1 = equilateral
    sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes())); names = [s["id"] for s in sb["segments"]]
    Wt = np.asarray(sb["weights"], np.float32); owner = ((Wt[F[:, 0]] + Wt[F[:, 1]] + Wt[F[:, 2]]) / 3).argmax(1)[ext]
    bad = np.where(dots < 0)[0]
    say(f"exterior triangles {len(ext):,}; aspect ratio (longest edge / 2*inradius, 1 = equilateral): "
        f"median {np.median(aspect):.2f}, 99% {np.quantile(aspect, .99):.2f}, 99.99% {np.quantile(aspect, .9999):.2f}, max {aspect.max():.1f}")
    out = []
    for i in bad:
        pct = 100.0 * float((aspect < aspect[i]).mean())
        rank = int((aspect > aspect[i]).sum()) + 1
        out.append(dict(triangle=int(ext[i]), segment=names[owner[i]], aspect_ratio=float(aspect[i]),
                        percentile=pct, rank_worst=rank, area_mm2=1e6 * float(area[i]),
                        edges_mm=[1e3 * float(x) for x in edges[i]], normal_dot=float(dots[i]),
                        area_ratio_after=float(np.linalg.norm(n1[i]) / max(np.linalg.norm(n0[i]), 1e-30))))
        say(f"  triangle {ext[i]} on {names[owner[i]]}: aspect {aspect[i]:.1f} (worse than {pct:.4f}% of the mesh; "
            f"rank {rank} of {len(ext):,}), area {1e6*area[i]:.4f} mm2, edges {edges[i][0]*1e3:.2f}/{edges[i][1]*1e3:.2f}/"
            f"{edges[i][2]*1e3:.2f} mm, normal dot {dots[i]:.3e}")
    (OUT / "slivers.json").write_text(json.dumps(dict(inverted=out, exterior_triangles=int(len(ext)),
        aspect_median=float(np.median(aspect)), aspect_p99=float(np.quantile(aspect, .99)),
        aspect_p9999=float(np.quantile(aspect, .9999)), aspect_max=float(aspect.max())), indent=2) + "\n")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("fit", "score", "diagnose", "slivers", "anchors", "recovery",
                                        "recovery-separation", "recovery-shooting", "subset",
                                        "recovery-controlC", "recovery-hybrid", "per-segment", "twobytwo", "trim-sweep", "toe-geometry", "seam", "contact-folds", "crawl-clearance"), required=True)
    ap.add_argument("--family", choices=("spline", "flow", "anchored", "anchored-surface"), default="spline",
                    help="spline: the v1 thin-plate warp (data/derived/skin-warp-v1). "
                         "flow: the v2 stationary-velocity-field flow, which cannot fold (skin-warp-v2). "
                         "anchored: the v3 spline with anchors for skin that has no bone under it (skin-warp-v3).")
    ap.add_argument("--out", default=None, help="override the family's output directory")
    ap.add_argument("--bundle", default=None, help="override the family's contact-mesh bundle path")
    ap.add_argument("--trajectory", default="data/derived/crawl-best/trajectory.json",
                    help="crawl-clearance only: the trajectory to measure against; it is audited against the "
                         "declared joint ranges and the verdict is recorded with the number")
    args = ap.parse_args()
    # module-level rebinding, so every stage reads the family's own directories
    FAMILY = args.family
    if FAMILY == "flow":
        OUT = ROOT / "data/derived/skin-warp-v2"
        BUNDLE = "data/derived/segment-contact-meshes/skin-warp-v2"
    elif FAMILY == "anchored":
        OUT = ROOT / "data/derived/skin-warp-v3"
        BUNDLE = "data/derived/segment-contact-meshes/skin-warp-v3"
    elif FAMILY == "anchored-surface":
        # c0e88a2: the same instrument, with "no bone under it" measured to the bone SURFACE
        ANCHOR_MEASURE = "surface"
        OUT = ROOT / "data/derived/skin-warp-v4"
        BUNDLE = "data/derived/segment-contact-meshes/skin-warp-v4"
    if args.out: OUT = ROOT / args.out
    if args.bundle: BUNDLE = args.bundle
    {"fit": stage_fit, "score": stage_score, "diagnose": stage_diagnose, "slivers": stage_slivers,
     "anchors": stage_anchors, "recovery": stage_recovery,
     "recovery-separation": stage_recovery_separation,
     "recovery-shooting": lambda: stage_recovery_separation("shooting"),
     "subset": stage_subset,
     "recovery-controlC": lambda: stage_recovery_separation("nearest_subset"),
     "recovery-hybrid": lambda: stage_recovery_separation("hybrid"),
     "per-segment": stage_per_segment,
     "twobytwo": stage_twobytwo,
     "trim-sweep": stage_trim_sweep,
     "toe-geometry": stage_toe_geometry,
     "seam": stage_seam,
     "contact-folds": stage_contact_folds,
     "crawl-clearance": lambda: stage_crawl_clearance(args.trajectory)}[args.stage]()
