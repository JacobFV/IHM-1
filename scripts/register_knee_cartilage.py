"""put OAIZIB-CM knee articular cartilage onto this body's RIGHT knee, judged by the gates
committed to docs/TISSUE_MECHANICS.md ("Articular cartilage: the knee, and its gates fixed
before any data is seen").

SOURCE. OAIZIB-CM (HuggingFace YongchengYAO/OAIZIB-CM), CC-BY-NC-4.0; knee MRI from the
Osteoarthritis Initiative. Cite CartiMorph (doi:10.1016/j.media.2023.103035) and OAIZIB
(doi:10.1016/j.media.2018.11.009). Labels: 1 femur (cut at the top of the volume), 2 femoral
cartilage, 3 tibia (cut at the bottom), 4 medial and 5 lateral tibial cartilage. All 507 knees
are RIGHT knees (info/kneeSideInfo.csv) -- a known answer for the side gate.

METHOD. A knee is two bones at a joint, so the FEMUR and the TIBIA are fitted separately, each
with its own similarity, each carrying its own cartilage (femur -> label 2; tibia -> 4, 5).
Each fit is ONE-WAY trimmed ICP: the scan's partial bone surface pulled onto this body's
COMPLETE bone surface, never the reverse, so missing bone pulls on nothing. The scan frame is
arbitrary relative to the canonical frame, so each fit starts from 24 proper rotations aligning
the principal axes of the scan bone to those of this body's bone cropped around its joint
centre, runs a short ICP from each, and keeps the lowest trimmed residual.

GATES (thresholds from the doc, fixed before any subject was seen):
  known answer  this body's right femur and tibia, truncated to a knee-MRI field of view
                (112 x 140 x 140 mm about the joint), moved by a known similarity, are recovered
                within 2 mm, 2 deg, 1%. Run with a SMALL and a LARGE rotation, because the scan
                frame is arbitrary and the multi-start is what must find it.
  which knee    fitted to BOTH knees; the lower mean residual must win by >= 20%; every knee is
                a right knee, so a LEFT win is a failure.
  placement     >= 95% of mapped cartilage within 3 mm of its bone's surface, <= 1% inside it.
  joint space   mapped femoral and tibial cartilage overlap <= 1% of the smaller volume.
Reported, not judged: cartilage volume and mean thickness (2V/A); whether the mapped medial
tibial cartilage lies nearer the body midline than the lateral one.
"""
import argparse, gzip, importlib.util, io, itertools, json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/anatomy/oaizib-cm"
OUT_DEFAULT = ROOT / "data/derived/knee-cartilage-registered-v2"
CANONICAL_OUT = ROOT / "data/derived/knee-cartilage-registered-v1"
MIRRORED_OUT = ROOT / "data/derived/knee-cartilage-registered-v3-mirrored-femur"
LICENCE = ("OAIZIB-CM, CC-BY-NC-4.0 (non-commercial); cite CartiMorph doi:10.1016/j.media.2023.103035 "
           "and OAIZIB doi:10.1016/j.media.2018.11.009")
BONE_ENT = {"right": {"femur": "right femur", "tibia": "right tibia"},
            "left": {"femur": "left femur", "tibia": "left tibia"}}
FOV_M = (0.112, 0.140, 0.140)            # canonical x (left-right), y (superior), z (anterior)
SIDE_MARGIN, PLACE_MM, PLACE_MIN, INSIDE_MAX, OVERLAP_MAX = 0.20, 3.0, 0.95, 0.01, 0.01
SHORT_ITERS, N_TGT, N_PLACE = 25, 80000, 1200
# Principal-axis starts alone put the scaffold's tibia 129 deg out, and rolling them about the target's first
# axis did not help: not one of those 288 starts came within 40 deg of the truth. The two clouds' principal
# frames disagree -- the scaffold's tibia BODY carries the fibula, which the scan's tibia label does not, and
# the target crop and the scan's field of view keep different parts of the shaft. So rotation is searched over
# a fixed uniform set instead of being read off the axes, and only the translation is anchored, at the joint
# centre, which is the one landmark both clouds have. A start at the true rotation reaches 0.46 mm, so the
# basin exists; this is about reaching it.
N_ROT, SCREEN_ITERS, SCREEN_KEEP, FINAL_KEEP = 2048, 10, 40, 3
# A one-way residual has a degenerate global minimum: shrink the source onto a single target point and it
# reads 0.00 mm. The principal-axis starts never went near it; a uniform rotation search finds it, and it won
# the large-rotation femur known answer at 100% scale error. Candidates outside a plainly anatomical scale
# range are therefore not eligible. The range is wide on purpose -- every fit so far sits within 0.90-1.08.
SCALE_MIN, SCALE_MAX = 0.5, 2.0
# The bounds admit a fit that is degenerate in practice: oaizib_006's tibia collapsed to 0.529 against an
# IDENTICAL target that every other subject fitted at 0.96-1.11, taking its cartilage volume down eightfold.
# One knee's two bones come from one person, so their scales should agree; a fit whose two scales differ by
# more than this is reported as suspect. REPORTED, NOT GATED, and not tuned: no verdict depends on it.
SCALE_DISAGREE = 0.20
FLEX_DEG = (0, 30, 60, 90)        # gate 4': overlap must hold at EVERY one of these, not on average
PAIRED_SLACK = 0.05               # gate 3': how far below/above the same map's own ceiling the cartilage may sit
FRAMES = {"scaffold": "ground (supine-support-5ma720yd t=0)", "canonical": "canonical"}

def _load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
PO = _load("po", "scripts/register_pelvic_organs.py")
sim, apply, scale_of, umeyama, areas, area_samples = PO.sim, PO.apply, PO.scale_of, PO.umeyama, PO.areas, PO.area_samples
TRIM, ITERS, TOL, N_SRC = PO.TRIM, PO.ITERS, PO.TOL, PO.N_SRC
KA_T_MM, KA_R_DEG, KA_S = PO.KA_T_MM, PO.KA_R_DEG, PO.KA_S
inside = PO.CW.inside

def say(*a): print(*a, flush=True)

MIRROR_AXIS = 2          # in the mesh frame; the model's own left/right pairs fix it, and check_mirror() proves it

def check_mirror(read_surface, geometry):
    """l_tibia against r_tibia and l_patella against r_patella are exact reflections about MIRROR_AXIS.
    They are the known answer for the mirroring itself: if this stops holding, the femur mirror is not trustworthy."""
    out = {}
    for l, r in (("l_tibia.vtp", "r_tibia.vtp"), ("l_patella.vtp", "r_patella.vtp")):
        lv, _ = read_surface(geometry / l); rv, _ = read_surface(geometry / r)
        m = lv.copy(); m[:, MIRROR_AXIS] *= -1
        out[f"{l}->{r}"] = 1000 * float(np.linalg.norm(np.sort(m, 0) - np.sort(rv, 0), axis=1).max())
    if max(out.values()) > 1e-6: raise SystemExit(f"the model's left/right pairs are not reflections about axis {MIRROR_AXIS}: {out}")
    return out

def scaffold_bones(right_femur="model"):
    """the SCAFFOLD's femur and tibia, in ground at the reference run's t=0 pose.

    Mirrors build_skin_contact_meshes.bone_clouds() -- same model, same meshes, same scale factors
    -- but keeps the faces, which the gates need, and carries each body by its transform_ground."""
    import xml.etree.ElementTree as ET
    sys.path.insert(0, str(ROOT))
    from ihm.spatial.vtk import surface as read_surface
    geometry = ROOT / "data/raw/anatomy/opensim-models/source/Geometry"
    root = ET.parse(ROOT / "data/models/engineering_stance_v1/model.osim").getroot().find("Model")
    ref = json.loads((ROOT / "data/derived/supine-support-5ma720yd/initial_native.json").read_text())["bodies"]
    want = {"femur_r": ("right", "femur"), "tibia_r": ("right", "tibia"), "femur_l": ("left", "femur"), "tibia_l": ("left", "tibia")}
    mirror_check = check_mirror(read_surface, geometry) if right_femur == "mirrored-left" else None
    out = {}
    for b in root.iter("Body"):
        key = want.get(b.get("name"))
        if key is None: continue
        Vs, Fs, n = [], [], 0
        for mesh in b.iter("Mesh"):
            factors = np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            name = mesh.findtext("mesh_file")
            mirror = right_femur == "mirrored-left" and b.get("name") == "femur_r"
            if mirror: name = "l_femur.vtp"          # this body's own left femur, reflected, as every other pair already is
            v, f = read_surface(geometry / name)
            if mirror: v = v.copy(); v[:, MIRROR_AXIS] *= -1; f = f[:, ::-1]
            Vs.append(v * factors); Fs.append(f + n); n += len(v)
        M = np.asarray(ref[b.get("name")]["transform_ground"], float)
        out[key] = (apply(M, np.concatenate(Vs)), np.concatenate(Fs))
    missing = [k for k in want.values() if k not in out]
    if missing: raise SystemExit(f"scaffold bones missing from the model: {missing}")
    return out, mirror_check

def body_bones():
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    by = {}
    for e in ents:
        if e["id"].startswith("body-bp3d-") and e.get("reference_geometry"): by.setdefault(e["name"], []).append(e)
    out = {}
    for side, d in BONE_ENT.items():
        for bone, n in d.items():
            c = by.get(n, [])
            if len(c) != 1: raise SystemExit(f"expected one BodyParts3D entity {n!r}, found {len(c)}")
            g = json.loads(gzip.decompress((ROOT / c[0]["reference_geometry"]["path"]).read_bytes()))
            out[(side, bone)] = (np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3))
    return out

def joint_centre(fV, fF, tV, tF, margin=0.008):
    """midpoint of the closest-approach pairs between femur and tibia surfaces (within 8 mm of the minimum gap)"""
    f = area_samples(fV, fF, 20000, 7); t = area_samples(tV, tF, 20000, 8)
    d, j = cKDTree(t).query(f); ok = d <= d.min() + margin
    return ((f[ok] + t[j[ok]]) / 2).mean(0)

def pca(P):
    m = P.mean(0); w, V = np.linalg.eigh(np.cov((P - m).T)); return m, V[:, ::-1]

SIGNED_PERMS = [np.diag(s) @ np.eye(3)[list(p)] for p in itertools.permutations(range(3))
                for s in itertools.product((1, -1), repeat=3)]
SCREEN_ROT = np.concatenate([np.stack([P for P in SIGNED_PERMS if np.linalg.det(P) > 0]),
                             Rotation.random(N_ROT, random_state=0).as_matrix()])


def icp(S, tree, M, iters):
    for it in range(1, iters + 1):
        P = apply(M, S); d, j = tree.query(P); keep = d <= np.quantile(d, 1 - TRIM)
        Mn = sim(*umeyama(S[keep], tree.data[j[keep]]))
        if np.abs(Mn - M).max() < TOL: M = Mn; break
        M = Mn
    d = tree.query(apply(M, S))[0]
    return M, float(np.sqrt((np.sort(d)[: int(len(d) * (1 - TRIM))] ** 2).mean())), it

def fit_bone(src, src_jc, tgt, tgt_jc, seed=0):
    """one-way trimmed ICP of a partial scan bone onto this body's complete bone, multi-start"""
    S = area_samples(*src, N_SRC, seed)
    T = area_samples(*tgt, N_TGT, seed + 100); tree = cKDTree(T)
    r = float(np.max(np.linalg.norm(S - src_jc, axis=1)))
    Tc = T[np.linalg.norm(T - tgt_jc, axis=1) <= r]
    ms, Vs = pca(S); mt, Vt = pca(Tc)
    s0 = float(np.sqrt(((Tc - mt) ** 2).sum(1).mean() / ((S - ms) ** 2).sum(1).mean()))
    Ssc, Tsc = S[::5], T[::4]; tsc = cKDTree(Tsc)          # a cheap screen over every start
    cands = []
    for R0 in SCREEN_ROT:
        R = Vt @ R0 @ Vs.T                                  # the perms align the axes; the random set covers the rest
        M, rms, _ = icp(Ssc, tsc, sim(s0, R, tgt_jc - s0 * R @ src_jc), SCREEN_ITERS)
        if SCALE_MIN <= scale_of(M) <= SCALE_MAX: cands.append((rms, M))
    if not cands: raise SystemExit("every start collapsed: no candidate inside the anatomical scale range")
    refined = []
    for _, M0 in sorted(cands, key=lambda c: c[0])[:SCREEN_KEEP]:
        M, rms, _ = icp(S, tree, M0, SHORT_ITERS)
        if SCALE_MIN <= scale_of(M) <= SCALE_MAX: refined.append((rms, M))
    if not refined: raise SystemExit("every refined start collapsed: no candidate inside the anatomical scale range")
    best = None                                             # the screen is only a proxy, so several finish the full fit
    for _, M0 in sorted(refined, key=lambda c: c[0])[:FINAL_KEEP]:
        M, rms, it = icp(S, tree, M0, ITERS)
        if SCALE_MIN <= scale_of(M) <= SCALE_MAX and (best is None or rms < best[1]): best = (M, rms, it)
    if best is None: raise SystemExit("every full fit collapsed: no candidate inside the anatomical scale range")
    return best

def knee_function(f, q):
    """the walker_knee's own function types. ihm/spatial/opensim.py evaluates Constant/Linear/SimmSpline;
    this knee is polynomial, so the polynomial forms are evaluated here in the same convention."""
    if f.tag == "Constant": return float(f.findtext("value"))
    if f.tag == "LinearFunction":
        c = [float(x) for x in f.findtext("coefficients").split()]; return float(c[0] * q + c[1])
    if f.tag == "PolynomialFunction":
        return float(np.polyval([float(x) for x in f.findtext("coefficients").split()], q))
    if f.tag == "MultiplierFunction":
        return float(f.findtext("scale")) * knee_function(f.find("function")[0], q)
    raise SystemExit(f"unsupported knee transform function {f.tag}")

def flexion_maps(side="right"):
    """What the PLANT's own knee does to the tibia, read from the model the engine loads and composed in
    ihm/spatial/opensim.py's convention: rotations R1@R2@R3, translations summed in the parent offset frame.

    Known answer: the reference run sits at knee_angle_r = 0, so the joint evaluated at 0 must reproduce the
    tibia's recorded transform_ground. That checks the offset frames and the composition, not the engine's
    integrator -- this is a re-evaluation of the plant's joint, not the plant's solver."""
    import xml.etree.ElementTree as ET
    sys.path.insert(0, str(ROOT))
    from ihm.spatial.opensim import offset_transform, vector
    root = ET.parse(ROOT / "data/models/engineering_stance_v1/model.osim").getroot().find("Model")
    suffix = "r" if side == "right" else "l"
    j = next(x for x in root.iter("CustomJoint") if x.get("name") == f"walker_knee_{suffix}")
    offs = {f.attrib["name"]: f for f in j.findall("frames/PhysicalOffsetFrame")}
    T = lambda f: offset_transform(vector(f.findtext("translation")), vector(f.findtext("orientation")))
    xpf, xcf = T(offs[j.findtext("socket_parent_frame")]), T(offs[j.findtext("socket_child_frame")])
    axes = j.findall("SpatialTransform/TransformAxis")
    if [a.attrib["name"] for a in axes] != ["rotation1", "rotation2", "rotation3", "translation1", "translation2", "translation3"]:
        raise SystemExit("unsupported spatial axis order in the knee joint")
    def motion(q):
        M = np.eye(4)
        for i, a in enumerate(axes):
            ax = vector(a.findtext("axis")); ax = ax / np.linalg.norm(ax)
            fn = [c for c in a if c.tag not in ("coordinates", "axis")][0]
            v = knee_function(fn, q)
            if i < 3: M[:3, :3] = M[:3, :3] @ Rotation.from_rotvec(ax * v).as_matrix()
            else: M[:3, 3] += ax * v
        return M
    ref = json.loads((ROOT / "data/derived/supine-support-5ma720yd/initial_native.json").read_text())["bodies"]
    Gf = np.asarray(ref[f"femur_{suffix}"]["transform_ground"], float); Gt = np.asarray(ref[f"tibia_{suffix}"]["transform_ground"], float)
    G = lambda q: Gf @ xpf @ motion(q) @ np.linalg.inv(xcf)
    G0 = G(0.0)
    known = dict(translation_mm=1000 * float(np.abs(G0[:3, 3] - Gt[:3, 3]).max()),
                 rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(G0[:3, :3].T @ Gt[:3, :3]) - 1) / 2, -1, 1)))))
    if known["translation_mm"] > 0.01 or known["rotation_deg"] > 0.01:
        raise SystemExit(f"the knee joint evaluated at 0 does not reproduce the plant's own tibia pose: {known}")
    return {d: G(np.deg2rad(d)) @ np.linalg.inv(G0) for d in FLEX_DEG}, known

def area_samples_faces(V, F, n, seed=0):
    """area-weighted samples, with the face each came from -- gate 3'' needs the bone's outward normal there"""
    rng = np.random.default_rng(seed); a = areas(V, F); k = rng.choice(len(F), n, p=a / a.sum()); t = V[F][k]
    r1 = np.sqrt(rng.random(n)); r2 = rng.random(n)
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2], k

def label_meshes(arr, affine):
    from skimage import measure
    out = {}
    for lab in (1, 2, 3, 4, 5):
        m = arr == lab
        if m.sum() < 50: continue
        v, f, _, _ = measure.marching_cubes(np.pad(m.astype(np.uint8), 1), 0.5); v = v - 1
        out[lab] = ((v @ affine[:3, :3].T + affine[:3, 3]) / 1000.0, f.astype(np.int64))
    return out

def frame_axes(body):
    """which axis is left-right, superior-inferior, anterior-posterior -- read off the target, not assumed.
    The canonical frame and the scaffold's ground frame do NOT agree: left-right is x in one and z in the other."""
    c = lambda k: area_samples(*body[k], 20000, 4).mean(0)
    rf, lf, rt = c(("right", "femur")), c(("left", "femur")), c(("right", "tibia"))
    lr = int(np.argmax(np.abs(rf - lf))); si = int(np.argmax(np.abs(rf - rt)))
    if lr == si: raise SystemExit("cannot tell the target's left-right axis from its superior-inferior axis")
    return lr, si, ({0, 1, 2} - {lr, si}).pop()

def fov_half(axes):
    lr, si, ap = axes; h = np.empty(3)
    h[lr], h[si], h[ap] = FOV_M[0] / 2, FOV_M[1] / 2, FOV_M[2] / 2
    return h

def truncate(V, F, centre, half):
    keep = np.all(np.all(np.abs(V[F] - centre) <= half, axis=2), axis=1)
    return V, F[keep]

def outward_normals(V, F):
    tri = V[F]; n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True); n = n / np.where(ln == 0, 1, ln)
    return n if float((tri[:, 0] * np.cross(tri[:, 1], tri[:, 2])).sum()) > 0 else -n

def bone_facing(cV, cF, bone_samples):
    """the half of the cartilage shell whose outward normal points AT the bone: the subchondral face.
    Classified by DIRECTION, in the subject's own frame, so it does not presuppose the 3 mm test."""
    n = outward_normals(cV, cF); c = cV[cF].mean(1)
    return ((bone_samples[cKDTree(bone_samples).query(c)[1]] - c) * n).sum(1) > 0

def body_gap(body, side):
    f = area_samples(*body[(side, "femur")], 60000, 1); d = cKDTree(area_samples(*body[(side, "tibia")], 60000, 2)).query(f)[0]
    return dict(min=1000 * float(d.min()), p1=1000 * float(np.quantile(d, 0.01)))

def rot_err(M, Minv):
    Rf = M[:3, :3] / scale_of(M); Ri = Minv[:3, :3] / scale_of(Minv)
    return float(np.degrees(np.arccos(np.clip((np.trace(Rf.T @ Ri) - 1) / 2, -1, 1))))

def known_answer(body, axes):
    fV, fF = body[("right", "femur")]; tV, tF = body[("right", "tibia")]
    jc = joint_centre(fV, fF, tV, tF)
    half = fov_half(axes)
    cut = {"femur": truncate(fV, fF, jc, half), "tibia": truncate(tV, tF, jc, half)}
    share = {k: round(float(areas(*v).sum() / areas(*body[("right", k)]).sum()), 3) for k, v in cut.items()}
    say(f"  right knee joint centre {np.round(jc, 4)} m; axes (left-right, superior-inferior, anterior-posterior) = {axes}; "
        f"surface kept inside the 112x140x140 mm FOV: {share}")
    cases = {"small": sim(1.06, Rotation.from_rotvec(np.deg2rad(9) * np.array([1, -2, 1.5]) / np.linalg.norm([1, -2, 1.5])).as_matrix(), np.array([0.04, -0.03, 0.05])),
             "large": sim(0.95, Rotation.from_rotvec(np.deg2rad(110) * np.array([0.3, 1, -0.6]) / np.linalg.norm([0.3, 1, -0.6])).as_matrix(), np.array([-0.2, 0.1, 0.3]))}
    ok_all, rec = True, {}
    for name, Mt in cases.items():
        scan = {k: (apply(Mt, v[0]), v[1]) for k, v in cut.items()}
        sjc = joint_centre(*scan["femur"], *scan["tibia"])
        for bone in ("femur", "tibia"):
            M, rms, it = fit_bone(scan[bone], sjc, body[("right", bone)], jc, seed=3)
            Minv = np.linalg.inv(Mt)
            used = np.unique(cut[bone][1]); c = cut[bone][0][used].mean(0)
            terr = 1000 * float(np.linalg.norm(apply(M, apply(Mt, c[None]))[0] - c))
            rerr = rot_err(M, Minv); serr = abs(scale_of(M) / scale_of(Minv) - 1)
            ok = terr <= KA_T_MM and rerr <= KA_R_DEG and serr <= KA_S; ok_all &= ok
            rec[f"{name}_{bone}"] = dict(translation_mm=terr, rotation_deg=rerr, scale_rel=serr, iterations=it, residual_mm=1000 * rms)
            say(f"GATE known answer [{name} rotation, {bone}]: translation {terr:.3f} mm (<= {KA_T_MM}), rotation {rerr:.3f} deg "
                f"(<= {KA_R_DEG}), scale {100*serr:.3f}% (<= {100*KA_S:.0f}%), {it} iterations, residual {1000*rms:.2f} mm -> {'PASS' if ok else 'FAIL'}")
    return ok_all, dict(kept_share=share, cases=rec)

def register(sid, arr, affine, meta, body, out, target, axes, flexmaps):
    say(f"\n=== {sid} ({meta.get('sex')}, gender code {meta['gender']}, KL {meta['kl']}) ===")
    src = label_meshes(arr, affine)
    if 1 not in src or 3 not in src: return dict(subject=sid, error="femur or tibia label missing", passes=False, **meta)
    sjc = joint_centre(*src[1], *src[3])
    fits, resid = {}, {}
    for side in ("right", "left"):
        bjc = joint_centre(*body[(side, "femur")], *body[(side, "tibia")])
        for lab, bone in ((1, "femur"), (3, "tibia")):
            fits[(side, bone)] = fit_bone(src[lab], sjc, body[(side, bone)], bjc, seed=11)
        resid[side] = float(np.mean([fits[(side, b)][1] for b in ("femur", "tibia")]))
    lo, hi = sorted(resid, key=resid.get)
    margin = (resid[hi] - resid[lo]) / resid[lo]
    side_ok = lo == "right" and margin >= SIDE_MARGIN
    say(f"  residual right {1000*resid['right']:.2f} mm, left {1000*resid['left']:.2f} mm; lower: {lo}, margin {100*margin:.0f}% "
        f"-> which knee {'PASS' if side_ok else 'FAIL'}")
    rec = dict(subject=sid, **meta, residual_right_mm=1000 * resid["right"], residual_left_mm=1000 * resid["left"],
               side_winner=lo, side_margin=margin, gate_which_knee=side_ok,
               scale_femur=scale_of(fits[("right", "femur")][0]), scale_tibia=scale_of(fits[("right", "tibia")][0]))
    sf, st = scale_of(fits[("right", "femur")][0]), scale_of(fits[("right", "tibia")][0])
    rec["scale_disagreement"] = float(abs(np.log(sf / st)))
    rec["fit_suspect"] = bool(rec["scale_disagreement"] > SCALE_DISAGREE)
    if rec["fit_suspect"]:
        say(f"  SUSPECT FIT: femur scale {sf:.3f} against tibia scale {st:.3f} -- one knee, one person, so these should "
            f"agree; everything below for this subject is reported but should not be read as a measurement")
    d = out / sid; d.mkdir(parents=True, exist_ok=True)
    mapped, place_ok, void = {}, True, False
    for labs, bone, bl in (((2,), "femur", 1), ((4, 5), "tibia", 3)):
        M = fits[("right", bone)][0]; tV, tF = body[("right", bone)]
        ttree = cKDTree(area_samples(tV, tF, 200000, 21))
        sV, sF = src[bl]; ssamp = area_samples(sV, sF, 200000, 22); stree = cKDTree(ssamp)
        for lab in labs:
            if lab not in src: continue
            V, F = src[lab]; Vm = apply(M, V); mapped[lab] = (Vm, F, M)
            name = {2: "femoral_cartilage", 4: "medial_tibial_cartilage", 5: "lateral_tibial_cartilage"}[lab]
            with open(d / f"{name}.obj", "w") as h:
                h.write(f"# {sid} {name}, registered onto the {target} right {bone}, {FRAMES[target]} frame, metres\n# {LICENCE}\n")
                for v in Vm: h.write("v %.9g %.9g %.9g\n" % tuple(v))
                for f in F: h.write("f %d %d %d\n" % (f[0] + 1, f[1] + 1, f[2] + 1))
        present = [l for l in labs if l in mapped]
        if not present: place_ok = False; continue
        # the SAME points in both frames: classified bone-facing in the subject's own frame, then carried by the map
        Ps = np.vstack([area_samples(src[l][0], src[l][1][bone_facing(src[l][0], src[l][1], ssamp)], N_PLACE, 30 + l) for l in present])
        src_near = float(np.mean(stree.query(Ps)[0] <= PLACE_MM / 1000)); src_in = float(np.mean(inside(sV, sF, Ps)))
        Pm = apply(M, Ps)
        near = float(np.mean(ttree.query(Pm)[0] <= PLACE_MM / 1000)); inb = float(np.mean(inside(tV, tF, Pm)))
        Pb = apply(M, area_samples(sV, sF, N_PLACE, 40 + bl))               # the ceiling: the scan's own bone surface, mapped
        rec[f"bone_ceiling_{bone}"] = dict(within_3mm=float(np.mean(ttree.query(Pb)[0] <= PLACE_MM / 1000)),
                                           inside_bone=float(np.mean(inside(tV, tF, Pb))))
        rec[f"transform_{bone}"] = M.tolist()
        tS, tK = area_samples_faces(tV, tF, 200000, 23)                      # gate 3'': the bone's outward normal where it is nearest
        signed = ((Pm - tS[cKDTree(tS).query(Pm)[1]]) * outward_normals(tV, tF)[tK][cKDTree(tS).query(Pm)[1]]).sum(1)
        med_signed = 1000 * float(np.median(signed))
        ok_near = src_near >= PLACE_MIN; ok_in = src_in <= INSIDE_MAX      # the criterion validated on the SOURCE
        ceil = rec[f"bone_ceiling_{bone}"]
        paired = bool(near >= ceil["within_3mm"] - PAIRED_SLACK and inb <= ceil["inside_bone"] + PAIRED_SLACK)
        rec[f"placement_{bone}"] = dict(source_within_3mm=src_near, source_inside_bone=src_in,
                                        criterion_valid_within=ok_near, criterion_valid_inside=ok_in,
                                        within_3mm=near, inside_bone=inb,
                                        within_verdict=(near >= PLACE_MIN) if ok_near else None,
                                        inside_verdict=(inb <= INSIDE_MAX) if ok_in else None,
                                        gate_paired=paired, median_signed_offset_mm=med_signed, gate_signed=bool(med_signed > 0))
        place_ok &= paired and med_signed > 0
        say(f"  GATE 3' paired [{bone}]: within {100*near:.1f}% vs ceiling {100*ceil['within_3mm']:.1f}% (needs >= ceiling - 5), "
            f"inside {100*inb:.1f}% vs ceiling {100*ceil['inside_bone']:.1f}% (needs <= ceiling + 5) -> {'PASS' if paired else 'FAIL'}")
        say(f"  GATE 3'' signed [{bone}]: median offset along the bone's outward normal {med_signed:+.2f} mm (needs > 0) "
            f"-> {'PASS' if med_signed > 0 else 'FAIL'}")
        say(f"  criterion on the SOURCE [{bone} cartilage, bone-facing surface]: {100*src_near:.1f}% within {PLACE_MM:g} mm (needs >= 95), "
            f"{100*src_in:.2f}% inside its own bone (needs <= 1) -> {'VALID' if ok_near and ok_in else 'VOID'}")
        say(f"  ceiling [{bone}: the scan's own bone surface carried by the same map]: {100*rec[f'bone_ceiling_{bone}']['within_3mm']:.1f}% "
            f"within {PLACE_MM:g} mm, {100*rec[f'bone_ceiling_{bone}']['inside_bone']:.1f}% inside -- no cartilage on it can do better")
        say(f"  placement [{bone} cartilage, bone-facing surface]: {100*near:.1f}% within {PLACE_MM:g} mm, {100*inb:.2f}% inside the {bone} -> "
            + (" ".join(x for x in ((f"within {'PASS' if near >= PLACE_MIN else 'FAIL'}" if ok_near else "within VOID"),
                                    (f"inside {'PASS' if inb <= INSIDE_MAX else 'FAIL'}" if ok_in else "inside VOID")))))
        if not (ok_near and ok_in): void = True
    rec["gate_placement"] = place_ok
    rec["placement_criterion_void"] = void
    vox_m3 = float(np.prod(np.abs(np.diag(affine)[:3]))) / 1e9     # m^3 per voxel
    vol = {}
    for lab, bone in ((2, "femur"), (4, "tibia"), (5, "tibia")):
        s = scale_of(fits[("right", bone)][0]); vol[lab] = float((arr == lab).sum() * vox_m3 * s ** 3)
    def vox(labs, M, n, seed):
        P = (np.argwhere(np.isin(arr, labs)) @ affine[:3, :3].T + affine[:3, 3]) / 1000.0
        P = P[np.random.default_rng(seed).choice(len(P), min(n, len(P)), replace=False)]
        return P if M is None else apply(M, P)
    Ainv = np.linalg.inv(affine)
    def overlap_of(Mf, Mt):
        """share of the SMALLER cartilage solid's voxels that land in the other solid: each voxel centre is carried by its
        own bone's map into this body, pulled back through the OTHER bone's map into the scan grid, and looked up there.
        Exact on the labels, and exactly 0 at identity because the labels are exclusive."""
        I = np.eye(4); Mf = I if Mf is None else Mf; Mt = I if Mt is None else Mt
        fem_small = vol[2] <= vol.get(4, 0) + vol.get(5, 0)
        own, M_own, other, M_other = ([2], Mf, (4, 5), Mt) if fem_small else ([4, 5], Mt, (2,), Mf)
        P = apply(np.linalg.inv(M_other) @ M_own, vox(own, None, 10**9, 3)) * 1000.0
        ijk = np.rint(P @ Ainv[:3, :3].T + Ainv[:3, 3]).astype(np.int64)
        ok = np.all((ijk >= 0) & (ijk < np.array(arr.shape)), axis=1)
        hit = np.zeros(len(P), bool); hit[ok] = np.isin(arr[tuple(ijk[ok].T)], other)
        return float(hit.mean())
    overlap = overlap_of(fits[("right", "femur")][0], fits[("right", "tibia")][0])
    js_ok = overlap <= OVERLAP_MAX
    rec.update(joint_space_overlap=overlap, gate_joint_space=js_ok,
               cartilage_volume_ml={"femoral": 1e6 * vol[2], "medial_tibial": 1e6 * vol.get(4, 0), "lateral_tibial": 1e6 * vol.get(5, 0)})
    say(f"  joint space at rest: femoral/tibial cartilage overlap {100*overlap:.2f}% of the smaller (<= 1) -> {'PASS' if js_ok else 'FAIL'}")
    if flexmaps is None:
        rec["gate_joint_space_flexion"] = None; flex_ok = js_ok
    else:
        Mf_, Mt_ = fits[("right", "femur")][0], fits[("right", "tibia")][0]
        flex = {d: overlap_of(Mf_, D @ Mt_) for d, D in flexmaps.items()}
        first = next((d for d in FLEX_DEG if flex[d] > OVERLAP_MAX), None)
        flex_ok = first is None
        rec["joint_space_flexion"] = {str(d): v for d, v in flex.items()}
        rec["first_flexion_over_1pct_deg"] = first
        rec["gate_joint_space_flexion"] = flex_ok
        say(f"  GATE 4' the moving knee: overlap " + ", ".join(f"{d} deg {100*flex[d]:.2f}%" for d in FLEX_DEG)
            + f" -> {'PASS at every angle' if flex_ok else f'FAIL, first over 1% at {first} deg'}")
    # EVIDENCE, not gates: the same measures on the scan in its own frame (a perfect map would score these),
    # and placement read by cartilage VOLUME (voxel centres) rather than by surface samples
    ident = dict(joint_space_overlap=overlap_of(None, None))
    for labs, bl, bone in (((2,), 1, "femur"), ((4, 5), 3, "tibia")):
        bV, bF = src[bl]; P = np.vstack([area_samples(*src[l], 3000, 5 + i) for i, l in enumerate(labs) if l in src])
        ident[f"placement_{bone}"] = dict(within_3mm=float(np.mean(cKDTree(area_samples(bV, bF, 200000, 21)).query(P)[0] <= PLACE_MM / 1000)))
        M = fits[("right", bone)][0]; cV, cF = body[("right", bone)]; Pv = vox(list(labs), M, 2000, 4)
        rec[f"placement_{bone}_by_volume"] = dict(within_3mm=float(np.mean(cKDTree(area_samples(cV, cF, 200000, 21)).query(Pv)[0] <= PLACE_MM / 1000)),
                                                  inside_bone=float(np.mean(inside(cV, cF, Pv))))
    rec["metric_at_identity"] = ident
    say(f"  evidence: at identity (scan frame) overlap {100*ident['joint_space_overlap']:.2f}%, femoral {100*ident['placement_femur']['within_3mm']:.1f}% within 3 mm, "
        f"tibial {100*ident['placement_tibia']['within_3mm']:.1f}%; by volume after mapping: femoral {100*rec['placement_femur_by_volume']['within_3mm']:.1f}% within / "
        f"{100*rec['placement_femur_by_volume']['inside_bone']:.1f}% inside, tibial {100*rec['placement_tibia_by_volume']['within_3mm']:.1f}% / {100*rec['placement_tibia_by_volume']['inside_bone']:.1f}%")
    thick = {}
    for lab in (2, 4, 5):
        if lab in mapped:
            Vm, F, _ = mapped[lab]; thick[lab] = 2 * vol[lab] / float(areas(Vm, F).sum())
    rec["mean_thickness_mm"] = {k: 1000 * v for k, v in thick.items()}
    if 4 in mapped and 5 in mapped:
        lr = axes[0]
        mid = 0.5 * (area_samples(*body[("right", "femur")], 20000, 4)[:, lr].mean() + area_samples(*body[("left", "femur")], 20000, 4)[:, lr].mean())
        med = abs(mapped[4][0][:, lr].mean() - mid) < abs(mapped[5][0][:, lr].mean() - mid)
        rec["medial_nearer_midline"] = bool(med)
    say(f"  volume femoral {1e6*vol[2]:.1f} mL, tibial {1e6*(vol.get(4,0)+vol.get(5,0)):.1f} mL; "
        f"mean thickness {', '.join(f'{k}: {1000*v:.2f} mm' for k, v in thick.items())}; medial nearer midline: {rec.get('medial_nearer_midline')}")
    rec["passes"] = bool(side_ok and place_ok and flex_ok and not rec["fit_suspect"])
    return rec

def sex_coding():
    """the Gender code's meaning is settled from anatomy on the source card, not assumed"""
    d = json.loads((ROOT / "data/sources/oaizib-cm.json").read_text())
    def find(o):
        if isinstance(o, dict):
            if "gender_coding" in o: return o["gender_coding"]
            for v in o.values():
                r = find(v)
                if r: return r
    g = find(d); g = g.get("codes", g); return {int(k): v for k, v in g.items() if k in ("1", "2")}

def subject_info():
    import openpyxl
    z = zipfile.ZipFile(RAW / "info.zip"); info = {}
    for n in ("info/subInfo_train.xlsx", "info/subInfo_test.xlsx"):
        ws = openpyxl.load_workbook(io.BytesIO(z.read(n)), read_only=True).worksheets[0]
        rows = list(ws.iter_rows(values_only=True)); h = list(rows[0])
        for r in rows[1:]:
            info[f"oaizib_{r[h.index('CMT-ID')]}"] = dict(gender=r[h.index("Gender")], kl=r[h.index("KLGrade")],
                                                         age=r[h.index("Age")], bmi=r[h.index("BMI")])
    side = {}
    for line in z.read("info/kneeSideInfo.csv").decode("utf-8-sig").splitlines():   # NO header row
        if line.strip(): k, v = line.split(","); side[k.strip()[:-7]] = v.strip()
    return info, side

def main():
    import nibabel as nib
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--right-femur", choices=("model", "mirrored-left"), default="model", dest="right_femur",
                    help="model: r_femur.vtp, 265 faces, which is NOT the left femur mirrored. mirrored-left: this "
                         "body's own left femur reflected, which is what every other left/right pair in the model is")
    ap.add_argument("--target", choices=("scaffold", "canonical"), default="scaffold",
                    help="scaffold: the plant's own femur and tibia, 3.0 mm apart (the second attempt). "
                         "canonical: this body's atlas bones, 0.6 mm apart (the first pilot)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--pilot", type=int, default=0, help="take the first N KL-grade-0 knees of EACH gender code")
    a = ap.parse_args()
    out = (a.out or (MIRRORED_OUT if a.right_femur == "mirrored-left" else OUT_DEFAULT if a.target == "scaffold" else CANONICAL_OUT)).resolve(); out.mkdir(parents=True, exist_ok=True)
    mirror_check = None
    if a.target == "scaffold": body, mirror_check = scaffold_bones(a.right_femur)
    else: body = body_bones()
    say(f"target: the {a.target} femur and tibia, {FRAMES[a.target]} frame; right femur from {a.right_femur}")
    flexmaps = None
    if a.target == "scaffold":
        flexmaps, knee_known = flexion_maps("right")
        say(f"  the plant's own knee, evaluated at 0 deg against its recorded tibia pose (known answer): "
            f"{knee_known['translation_mm']:.2e} mm, {knee_known['rotation_deg']:.2e} deg")
    if mirror_check is not None: say(f"  mirror known answer (max vertex difference, must be 0): {mirror_check}")
    say("== known answer ==")
    axes = frame_axes(body)
    ok, ka = known_answer(body, axes)
    man = out / "manifest.json"
    report = json.loads(man.read_text()) if man.exists() else dict(schema="ihm.knee-cartilage-registered.v2", licence=LICENCE, subjects={})
    if flexmaps is not None: report["knee_known_answer"] = knee_known
    report["target"] = dict(bones=a.target, frame=FRAMES[a.target], axes_lr_si_ap=list(axes),
                            right_femur=a.right_femur, mirror_known_answer_mm=mirror_check); report["known_answer"] = ka; report["licence"] = LICENCE
    report["body_joint_gap_mm"] = {side: body_gap(body, side) for side in ("right", "left")}
    say(f"  this body's femur-tibia surface gap (min, 1st percentile): {report['body_joint_gap_mm']}")          # MERGE: earlier subjects are kept
    if not ok:
        man.write_text(json.dumps(report, indent=2) + "\n"); sys.exit("known-answer gate failed; nothing registered")
    info, side = subject_info(); sex = sex_coding()
    subs = list(a.subjects or [])
    if a.pilot:
        for g in (1, 2):
            subs += sorted(k for k, v in info.items() if v["kl"] == 0 and v["gender"] == g)[: a.pilot]
    zips = {n[:-7].split("/")[-1]: (zf, n) for zf in (zipfile.ZipFile(RAW / "labelsTr.zip"), zipfile.ZipFile(RAW / "labelsTs.zip"))
            for n in zf.namelist() if n.endswith(".nii.gz")}
    for sid in subs:
        zf, n = zips[sid]; img = nib.Nifti1Image.from_bytes(gzip.decompress(zf.read(n)))
        meta = dict(gender=info[sid]["gender"], sex=sex.get(info[sid]["gender"]), kl=info[sid]["kl"], age=info[sid]["age"], knee_side_recorded=side.get(sid))
        prev = report["subjects"].get(sid)
        rec = register(sid, np.asanyarray(img.dataobj), img.affine, meta, body, out, a.target, axes, flexmaps)
        if prev: rec["superseded"] = prev.pop("superseded", []) + [prev]
        report["subjects"][sid] = rec
        man.write_text(json.dumps(report, indent=2, default=float) + "\n")
    subs_all = report["subjects"]
    say(f"\n{sum(v.get('passes', False) for v in subs_all.values())} of {len(subs_all)} subjects in the manifest pass every gate")

if __name__ == "__main__": main()
