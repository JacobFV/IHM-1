#!/usr/bin/env python3
"""register ONE female torso -- TotalSegmentator subject s0790 -- into this body's
canonical atlas frame, by the bones the two bodies share.

WHY. This body is male; a female breast and a female torso envelope come from s0790,
the one TotalSegmentator subject that passes every gate in
scripts/build_female_torso_from_totalsegmentator.py (chest whole, breast whole, 39 of
39 registration labels, laterality from the data). Its meshes are in that CT's own
world frame. This puts them in this body's frame.

HOW. Correspondences are bone to bone: 24 ribs, both clavicles, T1-T12, and the
sternum (manubrium + body + xiphoid). Each side is reduced to its AREA-WEIGHTED
surface centroid, so marching-cubes vertices (uniform) and BodyParts3D vertices (not
uniform) are measured the same way. ONE similarity by Umeyama on the 39 centroid
pairs, then an optional trimmed symmetric ICP between the concatenated bone surfaces.

DUPLICATES, found before fitting and excluded. This body carries two copies of four
of these bones, a BodyParts3D entity and a Z-Anatomy entity of the same bone:
  fifth / eleventh / twelfth thoracic vertebra  vs  'vertebra t5' / 't11' / 't12'
  'manubrium'                                   vs  'manubrium of sternum'
same extents within a few mm, centroids 3.6-8.4 mm apart, 1.5-3.0 mm median nearest
vertex. Only the BodyParts3D entities (id body-bp3d-*) are used, and the script refuses
to run if any entity would enter the fit twice.

DECLARED BEFORE ANY RUN, so no threshold is chosen after seeing a result:
  LOBO_RATIO_MAX   leave-one-bone-out median error <= 0.5 x the null (the median
                   distance between neighbouring rib centroids in this body)
  CONTAIN_MIN      this body's ribs and sternum inside the mapped female trunk >= 0.95
  BREAST_RIB_MAX   share of this body's rib vertices inside a mapped breast <= 0.01
  REFINE RULE      ICP-refined map is used only if it LOWERS the trimmed surface RMS
                   AND does not raise the centroid RMS by more than 5 mm; otherwise the
                   centroid fit is used.

WHY body_trunc AND NOT skin FOR CONTAINMENT. TotalSegmentator's `skin` label is a
shell with thickness; marching cubes of it gives two nested surfaces, and ray parity
through two surfaces reports every interior point as outside. The trunk mask's surface
is one closed surface. skin.obj is still mapped and written, as the deliverable.

licence: TotalSegmentator dataset CC-BY-4.0 (Zenodo 10047292); breasts and body
subtask models Apache-2.0. cite Wasserthal et al., Radiology: AI 2023,
doi:10.1148/ryai.230024.
"""
import gzip, hashlib, importlib.util, json, subprocess, sys, time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
SUBJECT = "s0790"
SRC = ROOT / "data/derived/female-torso-totalsegmentator-v1" / SUBJECT
RAW = ROOT / "data/raw/anatomy/totalsegmentator" / SUBJECT
OUT = ROOT / "data/derived/female-torso-registered-v1"
ORD = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth"]
TRIM, ITERS, TOL, SAMPLES_PER_BONE = 0.10, 300, 1e-10, 600
LOBO_RATIO_MAX, CONTAIN_MIN, BREAST_RIB_MAX = 0.5, 0.95, 0.01
REFINE_CENTROID_SLACK_M = 0.005
MALE_SKIN_SAMPLES = 600_000

spec = importlib.util.spec_from_file_location("bscm", ROOT / "scripts/build_skin_contact_meshes.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
t0 = time.time()
def say(*a): print(*a, flush=True)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

# ------------------------------------------------------------------ geometry helpers
def read_obj(p):
    V, F = [], []
    for line in open(p):
        if line.startswith("v "): V.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "): F.append([int(x.split("/")[0]) - 1 for x in line.split()[1:4]])
    return np.asarray(V, float), np.asarray(F, np.int64)

def write_obj(p, V, F, header):
    with open(p, "w") as h:
        h.write(f"# {header}\n")
        for v in V: h.write("v %.9g %.9g %.9g\n" % tuple(v))
        for f in F: h.write("f %d %d %d\n" % (f[0] + 1, f[1] + 1, f[2] + 1))

def entity_mesh(e):
    g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
    return np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3)

def areas(V, F):
    t = V[F]; return np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2

def area_centroid(V, F):
    a = areas(V, F); return (a[:, None] * V[F].mean(1)).sum(0) / a.sum()

def sample(V, F, n, rng):
    tri = V[F]; a = areas(V, F); k = rng.choice(len(F), n, p=a / a.sum())
    r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = tri[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2]

def umeyama(a, b):
    ca, cb = a.mean(0), b.mean(0); x, y = a - ca, b - cb
    u, s, vt = np.linalg.svd(x.T @ y)
    d = np.diag([1.0, 1.0, np.sign(np.linalg.det(vt.T @ u.T))])
    r = vt.T @ d @ u.T; scale = float((s @ np.diag(d)).sum() / (x ** 2).sum())
    return scale, r, cb - scale * r @ ca

def sim(s, R, t): M = np.eye(4); M[:3, :3] = s * R; M[:3, 3] = t; return M
def apply(M, p): return p @ M[:3, :3].T + M[:3, 3]
def scale_of(M): return float(np.cbrt(np.linalg.det(M[:3, :3])))

def icp(src, tgt, M0):
    tt = cKDTree(tgt); M = M0.copy()
    for k in range(1, ITERS + 1):
        S = apply(M, src); d1, j1 = tt.query(S); d2, j2 = cKDTree(S).query(tgt)
        k1 = d1 <= np.quantile(d1, 1 - TRIM); k2 = d2 <= np.quantile(d2, 1 - TRIM)
        A = np.vstack([src[k1], src[j2[k2]]]); Bt = np.vstack([tgt[j1[k1]], tgt[k2]])
        Mn = sim(*umeyama(A, Bt))
        if np.abs(Mn - M).max() < TOL: return Mn, k
        M = Mn
    return M, ITERS

def surface_rms(M, src, tgt):
    S = apply(M, src); d = np.concatenate([cKDTree(tgt).query(S)[0], cKDTree(S).query(tgt)[0]])
    return float(np.sqrt((np.sort(d)[: int(len(d) * (1 - TRIM))] ** 2).mean()))

def parity_inside(V, F, points, seed=0):
    """per-point version of build_skin_contact_meshes.enclosure(): the same Moller-Trumbore
    two-opposite-ray parity, the same seeded ray direction, returning the flags instead of
    their mean. checked against enclosure() on identical points below."""
    rng = np.random.default_rng(seed); points = np.asarray(points, float)
    tri = np.asarray(V, float)[np.asarray(F)]; a = tri[:, 0]; e1 = tri[:, 1] - a; e2 = tri[:, 2] - a
    direction = rng.normal(size=3); direction /= np.linalg.norm(direction); counts = []
    for sign in (1., -1.):
        d = sign * direction; pvec = np.cross(d, e2); det = (e1 * pvec).sum(axis=1)
        par = np.abs(det) < 1e-14; invd = np.where(par, 0., 1. / np.where(par, 1., det))
        hits = np.zeros(len(points), np.int64)
        for s in range(0, len(points), 256):
            blk = points[s:s + 256]; tv = blk[:, None, :] - a[None]
            u = (tv * pvec[None]).sum(2) * invd[None]; qv = np.cross(tv, e1[None])
            v = (qv * d).sum(2) * invd[None]; t = (qv * e2[None]).sum(2) * invd[None]
            hits[s:s + 256] = ((~par[None]) & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 1e-9)).sum(1)
        counts.append(hits % 2 == 1)
    return counts[0] & counts[1]

def mask_surface(path, affine_mm, step):
    import nibabel as nib
    from skimage import measure
    m = np.asanyarray(nib.load(str(path)).dataobj) > 0
    v, f, _, _ = measure.marching_cubes(np.pad(m.astype(np.uint8), 1), 0.5, step_size=step)
    v = v - 1
    return (v @ affine_mm[:3, :3].T + affine_mm[:3, 3]) / 1000.0, f.astype(np.int64)

def main():
    import nibabel as nib
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    # -------------------------------------------------------- correspondences
    pairs = {}
    for side in ("left", "right"):
        for i in range(12): pairs[f"rib_{side}_{i + 1}"] = [f"{side} {ORD[i]} rib"]
        pairs[f"clavicula_{side}"] = [f"{side} clavicle"]
    for i in range(12): pairs[f"vertebrae_T{i + 1}"] = [f"{ORD[i]} thoracic vertebra"]
    pairs["sternum"] = ["manubrium", "body of sternum", "xiphoid process"]
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    by_name = {}
    for e in ents:
        if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-"):
            by_name.setdefault(e["name"], []).append(e)
    used_ids = []
    labels = sorted(pairs)
    src_c, dst_c, src_surf, dst_surf, body_mesh = [], [], [], [], {}
    for lab in labels:
        Vs, Fs = read_obj(SRC / "meshes" / f"{lab}.obj")
        parts = []
        for nm in pairs[lab]:
            cand = by_name.get(nm, [])
            if len(cand) != 1: raise SystemExit(f"{lab}: expected exactly one BodyParts3D entity named {nm!r}, found {len(cand)}")
            used_ids.append(cand[0]["id"]); parts.append(entity_mesh(cand[0]))
        off = np.cumsum([0] + [len(v) for v, _ in parts[:-1]])
        Vb = np.vstack([v for v, _ in parts]); Fb = np.vstack([f + o for (_, f), o in zip(parts, off)])
        body_mesh[lab] = (Vb, Fb)
        src_c.append(area_centroid(Vs, Fs)); dst_c.append(area_centroid(Vb, Fb))
        src_surf.append(sample(Vs, Fs, SAMPLES_PER_BONE, rng)); dst_surf.append(sample(Vb, Fb, SAMPLES_PER_BONE, rng))
    if len(used_ids) != len(set(used_ids)): raise SystemExit("an entity would enter the fit twice")
    src_c, dst_c = np.array(src_c), np.array(dst_c)
    src_all, dst_all = np.vstack(src_surf), np.vstack(dst_surf)
    say(f"{len(labels)} bone correspondences, {len(used_ids)} body entities (all BodyParts3D, none twice)")

    # -------------------------------------------------------- gate a
    Rt = Rotation.from_rotvec(np.deg2rad(11) * np.array([1, -2, 3]) / np.sqrt(14)).as_matrix()
    Mt = sim(1.07, Rt, np.array([0.12, -0.31, 0.05]))
    err_a = float(np.abs(sim(*umeyama(src_c, apply(Mt, src_c))) - Mt).max())
    say(f"GATE a: recovery of a known similarity on s0790's own centroids, transform error {err_a:.2e} (< 1e-9) -> "
        + ("PASS" if err_a < 1e-9 else "FAIL"))
    if err_a >= 1e-9: sys.exit(1)

    # -------------------------------------------------------- fit
    Mc = sim(*umeyama(src_c, dst_c))
    res_c = np.linalg.norm(apply(Mc, src_c) - dst_c, axis=1)
    Mr, it = icp(src_all, dst_all, Mc)
    res_r = np.linalg.norm(apply(Mr, src_c) - dst_c, axis=1)
    surf_c, surf_r = surface_rms(Mc, src_all, dst_all), surface_rms(Mr, src_all, dst_all)
    rms = lambda r: float(np.sqrt((r ** 2).mean()))
    say(f"centroid fit : scale {scale_of(Mc):.4f}  centroid RMS {1000*rms(res_c):.1f} mm  max {1000*res_c.max():.1f} mm "
        f"({labels[int(res_c.argmax())]})  trimmed surface RMS {1000*surf_c:.1f} mm")
    say(f"ICP refined  : scale {scale_of(Mr):.4f}  centroid RMS {1000*rms(res_r):.1f} mm  max {1000*res_r.max():.1f} mm "
        f"({labels[int(res_r.argmax())]})  trimmed surface RMS {1000*surf_r:.1f} mm  [{it} iterations]")
    use_refined = surf_r < surf_c and rms(res_r) <= rms(res_c) + REFINE_CENTROID_SLACK_M
    M = Mr if use_refined else Mc
    say(f"REFINE RULE (declared): {'refined map USED' if use_refined else 'centroid fit USED'}")

    # -------------------------------------------------------- gate b
    held = []
    for i in range(len(labels)):
        keep = np.arange(len(labels)) != i
        Mi = sim(*umeyama(src_c[keep], dst_c[keep]))
        held.append(float(np.linalg.norm(apply(Mi, src_c[i:i + 1])[0] - dst_c[i])))
    held = np.array(held)
    neigh = [np.linalg.norm(dst_c[labels.index(f"rib_{s}_{k}")] - dst_c[labels.index(f"rib_{s}_{k + 1}")])
             for s in ("left", "right") for k in range(1, 12)]
    null = float(np.median(neigh)); ratio = float(np.median(held)) / null
    say(f"GATE b: leave-one-bone-out (centroid fit) median {1000*np.median(held):.1f} mm, max {1000*held.max():.1f} mm "
        f"({labels[int(held.argmax())]}); null = median neighbouring-rib centroid distance {1000*null:.1f} mm; "
        f"ratio {ratio:.3f} (<= {LOBO_RATIO_MAX}) -> " + ("PASS" if ratio <= LOBO_RATIO_MAX else "FAIL"))
    if ratio > LOBO_RATIO_MAX: sys.exit(1)

    # -------------------------------------------------------- gate c
    mc = apply(M, src_c); idx = labels.index
    c1 = np.linalg.norm(mc[idx("clavicula_left")] - dst_c[idx("clavicula_left")]) < np.linalg.norm(mc[idx("clavicula_left")] - dst_c[idx("clavicula_right")])
    c2 = np.linalg.norm(mc[idx("rib_left_6")] - dst_c[idx("rib_left_6")]) < np.linalg.norm(mc[idx("rib_left_6")] - dst_c[idx("rib_right_6")])
    say(f"GATE c: laterality after mapping -- s0790 left clavicle nearer this body's left clavicle: {c1}; "
        f"rib_left_6 nearer left sixth rib: {c2} -> " + ("PASS" if c1 and c2 else "FAIL"))
    if not (c1 and c2): sys.exit(1)

    # -------------------------------------------------------- map the meshes
    ct = nib.load(str(RAW / "ct.nii.gz")); A = ct.affine
    written = {}
    for nm in ("breast_left", "breast_right", "skin"):
        V, F = read_obj(SRC / "meshes" / f"{nm}.obj"); Vm = apply(M, V)
        write_obj(OUT / f"{nm}.obj", Vm, F, f"{SUBJECT} {nm}, registered into the canonical atlas frame, metres")
        written[nm] = (Vm, F)
    Vt, Ft = mask_surface(SRC / "body" / "body_trunc.nii.gz", A, step=2); Vt = apply(M, Vt)
    write_obj(OUT / "body_trunc_surface.obj", Vt, Ft, f"{SUBJECT} body_trunc surface (marching cubes, step 2), canonical frame, metres")
    say(f"mapped meshes written   [{time.time()-t0:.0f}s]")

    # -------------------------------------------------------- gate d
    import trimesh
    say(f"female trunk surface: {len(Ft):,} triangles, watertight {trimesh.Trimesh(Vt, Ft, process=False).is_watertight}")
    ribs_l = np.vstack([body_mesh[f"rib_left_{k}"][0] for k in range(1, 13)])
    ribs_r = np.vstack([body_mesh[f"rib_right_{k}"][0] for k in range(1, 13)])
    stern = body_mesh["sternum"][0]
    groups = {"left ribs": ribs_l, "right ribs": ribs_r, "sternum": stern}
    cont = {g: B.enclosure(Vt, Ft, p, samples=1500) for g, p in groups.items()}
    ok_d1 = all(v >= CONTAIN_MIN for v in cont.values())
    say("GATE d1: this body's bones inside the MAPPED female trunk -- "
        + ", ".join(f"{g} {v:.3f}" for g, v in cont.items()) + f" (each >= {CONTAIN_MIN}) -> " + ("PASS" if ok_d1 else "FAIL"))
    ribs_all = np.vstack([ribs_l, ribs_r])
    in_breast = {s: B.enclosure(*written[f"breast_{s}"], ribs_all, samples=3000) for s in ("left", "right")}
    ok_d2 = all(v <= BREAST_RIB_MAX for v in in_breast.values())
    say(f"GATE d2: this body's rib vertices inside a mapped breast -- left {in_breast['left']:.4f}, right {in_breast['right']:.4f} "
        f"(each <= {BREAST_RIB_MAX}) -> " + ("PASS" if ok_d2 else "FAIL"))
    if not (ok_d1 and ok_d2):
        # a run stopped at gate d still writes the fitted transforms: they are exactly
        # what diagnosing the stop needs, and the first run discarded them
        (OUT / "manifest.json").write_text(json.dumps(dict(stopped_at="gate d", centroid_fit_transform=Mc.tolist(), icp_refined_transform=Mr.tolist(), transform=M.tolist(), used=("icp_refined" if use_refined else "centroid_fit"), containment=cont, rib_in_breast=in_breast), indent=2) + "\n")
        sys.exit(1)

    # -------------------------------------------------------- e: the envelope
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / B.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    used, inv = np.unique(F[ext], return_inverse=True); sv, sf = V[used], inv.reshape(-1, 3)
    if B.simtk_precondition(sv, sf) is not None:
        cv, cf = B.repair(sv, sf); sv, sf, _ = B.cap_boundaries(cv, cf)
    dense = sample(sv, sf, MALE_SKIN_SAMPLES, rng); tree = cKDTree(dense)
    spacing = float(np.sqrt(areas(sv, sf).sum() / MALE_SKIN_SAMPLES))
    probe = written["breast_left"][0][rng.choice(len(written["breast_left"][0]), 1500, replace=False)]
    pm, em = float(parity_inside(sv, sf, probe).mean()), B.enclosure(sv, sf, probe, samples=len(probe))
    say(f"internal check: per-point parity {pm:.6f} vs enclosure() {em:.6f} on the same 1500 points -> "
        + ("identical" if pm == em else "DIFFERENT -- the envelope numbers below are not trustworthy"))
    def signed_mm(P):
        d = tree.query(P)[0]; inside = parity_inside(sv, sf, P)
        return 1000 * np.where(inside, -d, d)          # positive = OUTSIDE the male skin
    env = {}
    for nm in ("breast_left", "breast_right", "body_trunc"):
        P = written[nm][0] if nm != "body_trunc" else Vt
        P = P[rng.choice(len(P), min(len(P), 4000), replace=False)]
        s = signed_mm(P)
        env[nm] = dict(median_mm=float(np.median(s)), p95_mm=float(np.percentile(s, 95)), max_mm=float(s.max()),
                       share_of_surface_outside=float((s > 0).mean()))
        say(f"envelope {nm:12s} signed distance to the male skin (+ = outside): median {env[nm]['median_mm']:+.1f} mm, "
            f"95th {env[nm]['p95_mm']:+.1f} mm, max {env[nm]['max_mm']:+.1f} mm; surface outside {env[nm]['share_of_surface_outside']:.3f}")
    br = (np.asanyarray(nib.load(str(SRC / "breasts" / "breast.nii.gz")).dataobj) > 0) & \
         (np.asanyarray(nib.load(str(SRC / "body" / "body.nii.gz")).dataobj) > 0)
    ijk = np.argwhere(br); vox_ml = abs(np.linalg.det(A[:3, :3])) / 1000.0
    pick = ijk[rng.choice(len(ijk), min(len(ijk), 20000), replace=False)]
    P = apply(M, (pick @ A[:3, :3].T + A[:3, 3]) / 1000.0)
    out_share = float((~parity_inside(sv, sf, P)).mean())
    env["breast_volume"] = dict(total_ml=float(len(ijk) * vox_ml), share_outside_male_skin=out_share,
                                ml_outside=float(out_share * len(ijk) * vox_ml))
    say(f"envelope breast volume: {env['breast_volume']['total_ml']:.0f} mL (breast clipped to body), "
        f"{100*out_share:.1f}% ({env['breast_volume']['ml_outside']:.0f} mL) outside the male skin; "
        f"male-skin sample spacing {1000*spacing:.2f} mm (distance resolution)")

    # -------------------------------------------------------- manifest
    try: git = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    except Exception: git = "unknown"
    inputs = [SRC / "meshes" / f"{l}.obj" for l in labels] + [SRC / "meshes" / f"{n}.obj" for n in ("breast_left", "breast_right", "skin")] \
        + [SRC / "body" / "body_trunc.nii.gz", SRC / "body" / "body.nii.gz", SRC / "breasts" / "breast.nii.gz", RAW / "ct.nii.gz",
           ROOT / "data/derived/canonical/anatomy.json", ROOT / "data/derived/canonical/mechanics.json", ROOT / B.EVIDENCE]
    (OUT / "manifest.json").write_text(json.dumps(dict(
        schema="ihm.female-torso-registered.v1", subject=SUBJECT,
        frame_from=f"TotalSegmentator {SUBJECT} CT world frame (RAS), metres", frame_to="canonical atlas frame (bodyparts3d-display-m), metres",
        transform=M.tolist(), scale=scale_of(M), used=("icp_refined" if use_refined else "centroid_fit"),
        centroid_fit=dict(transform=Mc.tolist(), scale=scale_of(Mc), centroid_rms_m=rms(res_c), centroid_max_m=float(res_c.max()), surface_trimmed_rms_m=surf_c),
        icp_refined=dict(transform=Mr.tolist(), scale=scale_of(Mr), iterations=it, centroid_rms_m=rms(res_r), centroid_max_m=float(res_r.max()), surface_trimmed_rms_m=surf_r),
        correspondences={l: pairs[l] for l in labels}, body_entities=used_ids,
        duplicates_excluded=["vertebra t5", "vertebra t11", "vertebra t12", "manubrium of sternum"],
        gates=dict(a_recovery_error=err_a, b_lobo_median_m=float(np.median(held)), b_lobo_max_m=float(held.max()), b_null_m=null, b_ratio=ratio,
                   c_laterality=bool(c1 and c2), d_containment=cont, d_rib_in_breast=in_breast),
        envelope=env, thresholds=dict(lobo_ratio_max=LOBO_RATIO_MAX, contain_min=CONTAIN_MIN, breast_rib_max=BREAST_RIB_MAX,
                                      refine_centroid_slack_m=REFINE_CENTROID_SLACK_M),
        licence="dataset CC-BY-4.0 (Zenodo 10047292); breasts and body subtask models Apache-2.0",
        citation="Wasserthal et al., Radiology: Artificial Intelligence 2023, doi:10.1148/ryai.230024",
        caveat="one clinical subject's anatomy as a segmentation model drew it; not the skeleton's subject nor the skin's",
        provenance=dict(git_sha=git, script=str(Path(__file__).resolve().relative_to(ROOT)), script_sha256=sha(__file__),
                        inputs={str(p.relative_to(ROOT)): sha(p) for p in inputs})), indent=2) + "\n")
    say(f"wrote {(OUT / 'manifest.json').relative_to(ROOT)}   [{time.time()-t0:.0f}s]")

if __name__ == "__main__": main()
