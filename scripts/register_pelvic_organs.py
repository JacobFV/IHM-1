"""register UT-EndoMRI uterus and ovaries into this body by pelvic bone landmarks.

gates fixed before any fitting: docs/BODY_PARAMETERS.md, "Uterus and ovaries: the pelvic
registration, and its gates, fixed before any fitting".

WHY ONE-WAY. The T2 is 5 mm slices over ~160-210 mm, so the MRI's femurs, and possibly its
iliac crests, are CUT. A cut bone's centroid is dragged toward the scan (s1067's cut ribs did
exactly that). So ONE similarity is fitted by one-way trimmed ICP: points on the MRI's PARTIAL
bone surfaces (marching cubes of TotalSegmentator total_mr labels, world metres) find their
nearest points on this body's COMPLETE matching bone -- per label, hip to hip, sacrum to
sacrum -- and never the reverse, so missing bone pulls on nothing.

INITIALISATION, stated: Umeyama on three centroids -- hip_left, hip_right, sacrum -- MRI vs
this body. Those centroids are cut-biased, which is acceptable for an initialiser only; the
femurs, the most cut, are excluded from it. The known-answer gate exercises this same path.

GATES:
  known answer  truncate this body's own pelvis as a scan would (keep a 160 mm superior-
                inferior band: iliac crests and femurs cut), move it by a known similarity, and
                the one-way fit recovers it: translation within 2 mm (at the pelvic centroid),
                rotation within 2 deg, scale within 1%
  laterality    the MRI's left hip bone maps nearer this body's left hip bone
  containment   the mapped uterus lies inside this body's pelvic ring (convex hull of hip bones
                and sacrum) -- read as >= 0.99 of area samples of its surface inside the hull --
                and overlaps no bone: <= 1% of its interior points inside any bone mesh
REPORTED, NOT JUDGED: overlap with this male body's prostate, seminal vesicles and bladder
(expected: no female pelvis variant exists), the trimmed one-way residual, and the ovaries.

CAVEAT on every artefact: endometriosis cohort, pathology-selected, NOT a typical-anatomy
reference; non-commercial research use; cite Liang et al.
"""
import argparse, gzip, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree, Delaunay
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
OUT_DEFAULT = ROOT / "data/derived/ut-endomri-registered-v1"
CAVEAT = ("endometriosis cohort, pathology-selected, NOT a typical-anatomy reference; "
          "non-commercial research use; cite Liang et al.")
BONES = {"hip_left": "left hip bone", "hip_right": "right hip bone", "sacrum": "sacrum",
         "femur_left": "left femur", "femur_right": "right femur"}
INIT = ("hip_left", "hip_right", "sacrum")
MALE = ("prostate", "left seminal vesicle", "right seminal vesicle", "urinary bladder")
TRIM, ITERS, TOL, N_SRC = 0.10, 400, 1e-10, 4000
KA_T_MM, KA_R_DEG, KA_S = 2.0, 2.0, 0.01
RING_MIN, BONE_MAX = 0.99, 0.01

def load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
CW = load("cwo", "scripts/measure_female_chest_wall_offset.py")   # per-point parity inside()

def say(*a): print(*a, flush=True)
def sim(s, R, t): M = np.eye(4); M[:3, :3] = s * R; M[:3, 3] = t; return M
def apply(M, p): return p @ M[:3, :3].T + M[:3, 3]
def scale_of(M): return float(np.cbrt(np.linalg.det(M[:3, :3])))
def umeyama(a, b):
    ca, cb = a.mean(0), b.mean(0); x, y = a - ca, b - cb
    u, s, vt = np.linalg.svd(x.T @ y); d = np.diag([1.0, 1.0, np.sign(np.linalg.det(vt.T @ u.T))])
    r = vt.T @ d @ u.T; return float((s @ np.diag(d)).sum() / (x ** 2).sum()), r, cb - (s @ np.diag(d)).sum() / (x ** 2).sum() * r @ ca
def areas(V, F): t = V[F]; return np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
def area_samples(V, F, n, seed=0):
    rng = np.random.default_rng(seed); t = V[F]; a = areas(V, F); k = rng.choice(len(F), n, p=a / a.sum())
    r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = t[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2]
def read_obj(p):
    V, F = [], []
    for l in open(p):
        if l.startswith("v "): V.append([float(x) for x in l.split()[1:4]])
        elif l.startswith("f "): F.append([int(x.split("/")[0]) - 1 for x in l.split()[1:4]])
    return np.asarray(V, float), np.asarray(F, np.int64)
def write_obj(p, V, F, header):
    with open(p, "w") as h:
        h.write(f"# {header}\n# {CAVEAT}\n")
        for v in V: h.write("v %.9g %.9g %.9g\n" % tuple(v))
        for f in F: h.write("f %d %d %d\n" % (f[0] + 1, f[1] + 1, f[2] + 1))

def body_meshes():
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    by = {}
    for e in ents:
        if e["id"].startswith("body-bp3d-") and e.get("reference_geometry"): by.setdefault(e["name"], []).append(e)
    out = {}
    for n in list(BONES.values()) + list(MALE):
        c = by.get(n, [])
        if len(c) != 1: raise SystemExit(f"expected exactly one BodyParts3D entity named {n!r}, found {len(c)}")
        g = json.loads(gzip.decompress((ROOT / c[0]["reference_geometry"]["path"]).read_bytes()))
        out[n] = (np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3))
    return out

def mr_surfaces(t2, seg):
    import nibabel as nib
    from skimage import measure
    A = nib.load(str(t2)).affine; out = {}
    for lab in BONES:
        p = Path(seg) / f"{lab}.nii.gz"
        if not p.exists(): continue
        m = np.asanyarray(nib.load(str(p)).dataobj) > 0
        if m.sum() < 50: continue
        v, f, _, _ = measure.marching_cubes(np.pad(m.astype(np.uint8), 1), 0.5); v = v - 1
        out[lab] = ((v @ A[:3, :3].T + A[:3, 3]) / 1000.0, f.astype(np.int64))
    return out

def fit(src, body, rng_seed=0):
    """src: {label: (V, F)} partial; body: {label: (V, F)} complete, same labels. one-way ICP."""
    labs = [l for l in src if l in body]
    for l in INIT:
        if l not in labs: raise SystemExit(f"initialiser needs {l}; the scan has no usable {l}")
    ca = np.array([area_samples(*src[l], 3000).mean(0) for l in INIT])
    cb = np.array([area_samples(*body[l], 3000).mean(0) for l in INIT])
    M = sim(*umeyama(ca, cb))
    S = {l: area_samples(*src[l], N_SRC, rng_seed + i) for i, l in enumerate(labs)}
    T = {l: cKDTree(area_samples(*body[l], 40000, 100 + i)) for i, l in enumerate(labs)}
    for it in range(1, ITERS + 1):
        A_, B_, D_ = [], [], []
        for l in labs:
            P = apply(M, S[l]); d, j = T[l].query(P)
            A_.append(S[l]); B_.append(T[l].data[j]); D_.append(d)
        A_, B_, D_ = np.vstack(A_), np.vstack(B_), np.concatenate(D_)
        keep = D_ <= np.quantile(D_, 1 - TRIM)
        Mn = sim(*umeyama(A_[keep], B_[keep]))
        if np.abs(Mn - M).max() < TOL: M = Mn; break
        M = Mn
    d = np.concatenate([T[l].query(apply(M, S[l]))[0] for l in labs])
    rms = float(np.sqrt((np.sort(d)[: int(len(d) * (1 - TRIM))] ** 2).mean()))
    return M, it, rms, labs

def known_answer(body):
    """truncate this body's own pelvis to a 160 mm SI band (crests and femurs cut), move it by a
    known similarity, and require the one-way fit to recover it."""
    B = {lab: body[n] for lab, n in BONES.items()}
    allv = np.vstack([v for v, _ in B.values()]); c = allv.mean(0)
    hips = np.vstack([B["hip_left"][0], B["hip_right"][0]])
    y_mid = float(np.median(hips[:, 1])); lo, hi = y_mid - 0.08, y_mid + 0.08   # +y superior in this frame
    cut = {}
    for lab, (V, F) in B.items():
        keepf = np.all((V[F][:, :, 1] >= lo) & (V[F][:, :, 1] <= hi), axis=1)
        if keepf.sum() < 50: continue
        cut[lab] = (V, F[keepf])
    frac = {lab: round(float(areas(V, F).sum() / areas(*B[lab]).sum()), 3) for lab, (V, F) in cut.items()}
    R = Rotation.from_rotvec(np.deg2rad(9) * np.array([1, -2, 1.5]) / np.linalg.norm([1, -2, 1.5])).as_matrix()
    Mt = sim(1.06, R, np.array([0.04, -0.03, 0.05]))                  # body -> "scan"
    scan = {lab: (apply(Mt, V), F) for lab, (V, F) in cut.items()}
    M, it, rms, _ = fit(scan, B)
    Minv = np.linalg.inv(Mt)
    terr = 1000 * float(np.linalg.norm(apply(M, apply(Mt, c[None]))[0] - c))
    Rf = M[:3, :3] / scale_of(M); Ri = Minv[:3, :3] / scale_of(Minv)
    rerr = float(np.degrees(np.arccos(np.clip((np.trace(Rf.T @ Ri) - 1) / 2, -1, 1))))
    serr = abs(scale_of(M) / scale_of(Minv) - 1)
    ok = terr <= KA_T_MM and rerr <= KA_R_DEG and serr <= KA_S
    say(f"  kept surface share after the 160 mm cut: {frac}")
    say(f"GATE known answer: translation {terr:.3f} mm (<= {KA_T_MM}), rotation {rerr:.3f} deg (<= {KA_R_DEG}), "
        f"scale {100*serr:.3f}% (<= {100*KA_S:.0f}%), {it} iterations, trimmed residual {1000*rms:.2f} mm -> {'PASS' if ok else 'FAIL'}")
    return ok, dict(translation_mm=terr, rotation_deg=rerr, scale_rel=serr, iterations=it, kept_share=frac)

def interior_points(V, F, n, seed=0):
    rng = np.random.default_rng(seed); lo, hi = V.min(0), V.max(0); pts = []
    while sum(len(p) for p in pts) < n:
        P = lo + (hi - lo) * rng.random((4 * n, 3)); pts.append(P[CW.inside(V, F, P)])
    return np.vstack(pts)[:n]

def register(sid, t2, seg, organs, out, body):
    say(f"\n=== {sid} ===")
    src = mr_surfaces(t2, seg)
    B = {lab: body[n] for lab, n in BONES.items()}
    M, it, rms, labs = fit(src, B)
    say(f"  labels fitted: {labs} | scale {scale_of(M):.4f} | {it} iterations | trimmed one-way residual {1000*rms:.2f} mm")
    cl = {l: apply(M, area_samples(*src[l], 3000)).mean(0) for l in ("hip_left", "hip_right")}
    bl = {l: area_samples(*B[l], 3000).mean(0) for l in ("hip_left", "hip_right")}
    lat = bool(np.linalg.norm(cl["hip_left"] - bl["hip_left"]) < np.linalg.norm(cl["hip_left"] - bl["hip_right"]))
    say(f"GATE laterality: MRI left hip maps nearer this body's left hip: {lat} -> {'PASS' if lat else 'FAIL'}")
    ring = Delaunay(np.vstack([body["left hip bone"][0], body["right hip bone"][0], body["sacrum"][0]]))
    rec = dict(subject=sid, transform=M.tolist(), scale=scale_of(M), iterations=it, residual_trimmed_mm=1000 * rms,
               labels_fitted=labs, gate_laterality=lat, caveat=CAVEAT)
    d = Path(out) / sid; d.mkdir(parents=True, exist_ok=True)
    Vu, Fu = read_obj(Path(organs) / "uterus.obj"); Vu = apply(M, Vu)
    write_obj(d / "uterus.obj", Vu, Fu, f"{sid} uterus, registered into the canonical frame, metres")
    ring_share = float((ring.find_simplex(area_samples(Vu, Fu, 20000)) >= 0).mean())
    inner = interior_points(Vu, Fu, 3000)
    bone_share = {n: float(CW.inside(*body[n], inner).mean()) for n in BONES.values()}
    male_share = {n: float(CW.inside(*body[n], inner).mean()) for n in MALE}
    contain = ring_share >= RING_MIN and max(bone_share.values()) <= BONE_MAX
    say(f"GATE containment: uterus surface inside the pelvic ring {ring_share:.4f} (>= {RING_MIN}); "
        f"interior inside a bone max {max(bone_share.values()):.4f} (<= {BONE_MAX}) -> {'PASS' if contain else 'FAIL'}")
    say(f"  reported: uterus interior inside male organs {({k: round(v, 3) for k, v in male_share.items()})}")
    ov = []
    for k in (1, 2):
        p = Path(organs) / f"ovary_{k}.obj"
        if not p.exists(): continue
        Vo, Fo = read_obj(p); Vo = apply(M, Vo)
        write_obj(d / f"ovary_{k}.obj", Vo, Fo, f"{sid} ovary piece {k}, registered into the canonical frame, metres")
        ov.append(dict(piece=k, inside_ring=float((ring.find_simplex(area_samples(Vo, Fo, 5000)) >= 0).mean())))
    say(f"  reported: ovaries {ov}")
    rec.update(gate_containment=contain, uterus_in_ring=ring_share, uterus_in_bone=bone_share, uterus_in_male_organs=male_share,
               ovaries=ov, passes=bool(lat and contain))
    return rec

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", nargs=4, action="append", metavar=("SUBJECT", "T2", "SEGDIR", "ORGANDIR"), default=[])
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    a = ap.parse_args()
    body = body_meshes()
    say("== known answer ==")
    ok, ka = known_answer(body)
    report = dict(schema="ihm.ut-endomri-registered.v1", caveat=CAVEAT, known_answer=ka, subjects={})
    a.out.mkdir(parents=True, exist_ok=True)
    # MANIFEST-MERGE. register_pelvic_organs_batch.py calls this once per subject into one output
    # directory, and every call used to REWRITE manifest.json with only its own subjects -- after
    # five subjects the shared manifest listed one. seed from what is already there, so both
    # writes below keep every subject registered before this call.
    prior = a.out / "manifest.json"
    if prior.exists():
        report["subjects"].update(json.loads(prior.read_text()).get("subjects", {}))
    if not ok:
        (a.out / "manifest.json").write_text(json.dumps(report, indent=2) + "\n"); sys.exit("known-answer gate failed; nothing registered")
    for sid, t2, seg, org in a.case:
        report["subjects"][sid] = register(sid, t2, seg, org, a.out, body)
        (a.out / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    n = sum(r["passes"] for r in report["subjects"].values())
    say(f"\n{n}/{len(report['subjects'])} subjects pass laterality and containment")

if __name__ == "__main__": main()
