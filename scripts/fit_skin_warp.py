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
from skin_warp import Warp, phi

OUT = ROOT / "data/derived/skin-warp-v1"
REGISTRATION = ROOT / "data/derived/anatomy-segment-registration/registration.json"
BUNDLE = "data/derived/segment-contact-meshes/skin-warp-v1"
REFERENCE_BUNDLE = ROOT / "data/derived/segment-contact-meshes/skin-binding"

# ---- fixed before the fit
N_CORR, N_CORR_TORSO, N_CORR_PELVIS, CORR_TRIM = 200, 600, 300, 0.10
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

def projector(V, F, rng, n=200000):
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
    Z = Warp.zero(Tb)
    ok = bool(np.array_equal(Z.apply(canonical), canonical @ Tb[:3, :3].T + Tb[:3, 3])); controls["zero_warp_bitwise"] = ok
    say(f"C2 zero displacement reproduces the binding map on all {len(canonical):,} skin vertices, bit for bit: {ok}")
    if not ok: sys.exit("C2 failed")
    Z.save(OUT / "zero.npz")
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
    S, T, L, per = [], [], [], {}
    for seg in segments:
        i = segments.index(seg); Va, Fa = atlas_mesh(seg); Vo, Fo = om[seg]; Vo_g = apply(rest[seg], Vo)
        n = N_CORR_TORSO if seg == "torso" else N_CORR_PELVIS if seg == "pelvis" else N_CORR
        a = sample(Va, Fa, n, np.random.default_rng(10000 + i))
        p = apply(Mseg[seg], a); t, d, _ = projector(Vo_g, Fo, np.random.default_rng(20000 + i))(p)
        keep = d <= np.quantile(d, 1 - CORR_TRIM)
        s = apply(G, a)
        S.append(s[keep]); T.append(t[keep]); L += [seg] * int(keep.sum())
        per[seg] = dict(sampled=n, kept=int(keep.sum()), projection_median_mm=1e3 * float(np.median(d)),
                        projection_kept_max_mm=1e3 * float(d[keep].max()),
                        displacement_from_global_median_mm=1e3 * float(np.median(np.linalg.norm(t[keep] - s[keep], axis=1))),
                        displacement_from_global_max_mm=1e3 * float(np.linalg.norm(t[keep] - s[keep], axis=1).max()))
        say(f"  {seg:10s} kept {keep.sum():4d}/{n}  |t - M a| median {per[seg]['projection_median_mm']:5.2f} mm   "
            f"|t - G a| median {per[seg]['displacement_from_global_median_mm']:5.1f} max {per[seg]['displacement_from_global_max_mm']:5.1f} mm")
    S, T, L = np.concatenate(S), np.concatenate(T), np.asarray(L); D = T - S
    say(f"  {len(S)} correspondences")

    # ---------------------------------------------- the regularisation rule
    say(f"\n== lambda by {FOLDS}-fold cross-validation on the bone correspondences  [{time.time()-t0:.0f}s] ==")
    perm = np.random.default_rng(7).permutation(len(S)); folds = np.array_split(perm, FOLDS)
    sse = np.zeros(len(LAMBDAS))
    for f in folds:
        tr = np.setdiff1d(perm, f); sy = tps_system(S[tr])
        for li, lam in enumerate(LAMBDAS):
            w, a = tps_solve(sy, D[tr], lam)
            sse[li] += ((S[f] + tps_predict(S[f], S[tr], w, a) - T[f]) ** 2).sum()
        say(f"  fold done [{time.time()-t0:.0f}s]")
    cv = np.sqrt(sse / len(S)); best = cv.min()
    lam = max(l for l, e in zip(LAMBDAS, cv) if e <= best * (1 + TIE))
    for l, e in zip(LAMBDAS, cv): say(f"  lambda {l:9.2e}  CV RMS {1e3*e:7.3f} mm" + ("   <- chosen" if l == lam else ""))
    say(f"  rule: largest lambda within {100*TIE:.0f}% of the minimum CV RMS ({1e3*best:.3f} mm) -> lambda = {lam:.2e}"
        + ("   (the grid's edge)" if lam in (LAMBDAS[0], LAMBDAS[-1]) else ""))
    sy = tps_system(S); w, a = tps_solve(sy, D, lam)
    W = Warp(G, S, w, a, dict(kind="thin-plate spline displacement on the binding map", lam=lam, correspondences=int(len(S)),
                              rule="largest lambda within 1% of the minimum 5-fold CV RMS on the bone correspondences"))
    W.save(OUT / "warp.npz")
    fit_res = np.linalg.norm(W.apply(np.linalg.solve(G[:3, :3], (S - G[:3, 3]).T).T) - T, axis=1)
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
    (OUT / "fit.json").write_text(json.dumps(dict(schema="ihm.skin-warp-fit.v1", controls=controls,
        correspondences=per, n_correspondences=int(len(S)), cv=dict(lambdas=LAMBDAS, rms_m=cv.tolist(), chosen=lam, tie=TIE, folds=FOLDS),
        bending_energy=W.bending_energy(), gate1=dict(passed=gate1, margin_m=GATE1_MARGIN_M, segments=g1),
        provenance=dict(script=str(Path(__file__).relative_to(ROOT)), script_sha256=sha(__file__),
                        inputs={str(Path(p).relative_to(ROOT)): sha(p) for p in (REGISTRATION, B.ANATOMY, B.BINDING, B.MODEL)})),
        indent=2) + "\n")
    say(f"wrote {(OUT/'warp.npz').relative_to(ROOT)}, {(OUT/'fit.json').relative_to(ROOT)}   [{time.time()-t0:.0f}s]")

# ================================================= stage score
def stage_score():
    fit = json.loads((OUT / "fit.json").read_text()); W = Warp.load(OUT / "warp.npz")
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    used = np.unique(F[ext])
    # gate 4: Jacobian determinant at every exterior skin vertex; triangle orientation before/after
    det = np.concatenate([np.linalg.det(W.jacobian(V[used[s:s + 2048]])) for s in range(0, len(used), 2048)])
    det0 = np.linalg.det(W.base[:3, :3])
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
                   triangles=len(dots), inverted=folds_t, caps=len(cdots), caps_inverted=int((cdots < 0).sum())),
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
    W = Warp.load(OUT / "warp.npz"); Tb, frames_b, _ = bscm.binding_registration()
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

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", choices=("fit", "score", "diagnose"), required=True)
    stage = ap.parse_args().stage
    {"fit": stage_fit, "score": stage_score, "diagnose": stage_diagnose}[stage]()
