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
OUT_DEFAULT = ROOT / "data/derived/knee-cartilage-registered-v1"
LICENCE = ("OAIZIB-CM, CC-BY-NC-4.0 (non-commercial); cite CartiMorph doi:10.1016/j.media.2023.103035 "
           "and OAIZIB doi:10.1016/j.media.2018.11.009")
BONE_ENT = {"right": {"femur": "right femur", "tibia": "right tibia"},
            "left": {"femur": "left femur", "tibia": "left tibia"}}
FOV_M = (0.112, 0.140, 0.140)            # canonical x (left-right), y (superior), z (anterior)
SIDE_MARGIN, PLACE_MM, PLACE_MIN, INSIDE_MAX, OVERLAP_MAX = 0.20, 3.0, 0.95, 0.01, 0.01
SHORT_ITERS, N_TGT = 25, 80000

def _load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
PO = _load("po", "scripts/register_pelvic_organs.py")
sim, apply, scale_of, umeyama, areas, area_samples = PO.sim, PO.apply, PO.scale_of, PO.umeyama, PO.areas, PO.area_samples
TRIM, ITERS, TOL, N_SRC = PO.TRIM, PO.ITERS, PO.TOL, PO.N_SRC
KA_T_MM, KA_R_DEG, KA_S = PO.KA_T_MM, PO.KA_R_DEG, PO.KA_S
inside = PO.CW.inside

def say(*a): print(*a, flush=True)

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
    best = None
    for P in SIGNED_PERMS:
        R = Vt @ P @ Vs.T
        if np.linalg.det(R) < 0: continue
        M, rms, _ = icp(S, tree, sim(s0, R, mt - s0 * R @ ms), SHORT_ITERS)
        if best is None or rms < best[1]: best = (M, rms)
    M, rms, it = icp(S, tree, best[0], ITERS)
    return M, rms, it

def label_meshes(arr, affine):
    from skimage import measure
    out = {}
    for lab in (1, 2, 3, 4, 5):
        m = arr == lab
        if m.sum() < 50: continue
        v, f, _, _ = measure.marching_cubes(np.pad(m.astype(np.uint8), 1), 0.5); v = v - 1
        out[lab] = ((v @ affine[:3, :3].T + affine[:3, 3]) / 1000.0, f.astype(np.int64))
    return out

def truncate(V, F, centre):
    half = np.array(FOV_M) / 2
    keep = np.all(np.all(np.abs(V[F] - centre) <= half, axis=2), axis=1)
    return V, F[keep]

def body_gap(body, side):
    f = area_samples(*body[(side, "femur")], 60000, 1); d = cKDTree(area_samples(*body[(side, "tibia")], 60000, 2)).query(f)[0]
    return dict(min=1000 * float(d.min()), p1=1000 * float(np.quantile(d, 0.01)))

def rot_err(M, Minv):
    Rf = M[:3, :3] / scale_of(M); Ri = Minv[:3, :3] / scale_of(Minv)
    return float(np.degrees(np.arccos(np.clip((np.trace(Rf.T @ Ri) - 1) / 2, -1, 1))))

def known_answer(body):
    fV, fF = body[("right", "femur")]; tV, tF = body[("right", "tibia")]
    jc = joint_centre(fV, fF, tV, tF)
    cut = {"femur": truncate(fV, fF, jc), "tibia": truncate(tV, tF, jc)}
    share = {k: round(float(areas(*v).sum() / areas(*body[("right", k)]).sum()), 3) for k, v in cut.items()}
    say(f"  right knee joint centre {np.round(jc, 4)} m; surface kept inside the 112x140x140 mm FOV: {share}")
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

def register(sid, arr, affine, meta, body, out):
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
    d = out / sid; d.mkdir(parents=True, exist_ok=True)
    mapped, place_ok = {}, True
    for labs, bone in (((2,), "femur"), ((4, 5), "tibia")):
        M = fits[("right", bone)][0]; bV, bF = body[("right", bone)]
        btree = cKDTree(area_samples(bV, bF, 200000, 21))
        for lab in labs:
            if lab not in src: continue
            V, F = src[lab]; Vm = apply(M, V); mapped[lab] = (Vm, F, M)
            name = {2: "femoral_cartilage", 4: "medial_tibial_cartilage", 5: "lateral_tibial_cartilage"}[lab]
            with open(d / f"{name}.obj", "w") as h:
                h.write(f"# {sid} {name}, registered onto this body's right {bone}, canonical frame, metres\n# {LICENCE}\n")
                for v in Vm: h.write("v %.9g %.9g %.9g\n" % tuple(v))
                for f in F: h.write("f %d %d %d\n" % (f[0] + 1, f[1] + 1, f[2] + 1))
        cart = [mapped[l] for l in labs if l in mapped]
        if not cart: place_ok = False; continue
        P = np.vstack([area_samples(Vm, F, 3000, 5 + i) for i, (Vm, F, _) in enumerate(cart)])
        near = float(np.mean(btree.query(P)[0] <= PLACE_MM / 1000))
        inb = float(np.mean(inside(bV, bF, P)))
        ok = near >= PLACE_MIN and inb <= INSIDE_MAX; place_ok &= ok
        rec[f"placement_{bone}"] = dict(within_3mm=near, inside_bone=inb, ok=ok)
        say(f"  placement [{bone} cartilage]: {100*near:.1f}% within {PLACE_MM:g} mm (>= 95), {100*inb:.2f}% inside the {bone} (<= 1) "
            f"-> {'PASS' if ok else 'FAIL'}")
    rec["gate_placement"] = place_ok
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
    say(f"  joint space: femoral/tibial cartilage overlap {100*overlap:.2f}% of the smaller (<= 1) -> {'PASS' if js_ok else 'FAIL'}")
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
        mid = 0.5 * (body[("right", "femur")][0][:, 0].mean() + body[("left", "femur")][0][:, 0].mean())
        med = abs(mapped[4][0][:, 0].mean() - mid) < abs(mapped[5][0][:, 0].mean() - mid)
        rec["medial_nearer_midline"] = bool(med)
    say(f"  volume femoral {1e6*vol[2]:.1f} mL, tibial {1e6*(vol.get(4,0)+vol.get(5,0)):.1f} mL; "
        f"mean thickness {', '.join(f'{k}: {1000*v:.2f} mm' for k, v in thick.items())}; medial nearer midline: {rec.get('medial_nearer_midline')}")
    rec["passes"] = bool(side_ok and place_ok and js_ok)
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
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--pilot", type=int, default=0, help="take the first N KL-grade-0 knees of EACH gender code")
    a = ap.parse_args()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    body = body_bones()
    say("== known answer ==")
    ok, ka = known_answer(body)
    man = out / "manifest.json"
    report = json.loads(man.read_text()) if man.exists() else dict(schema="ihm.knee-cartilage-registered.v1", licence=LICENCE, subjects={})
    report["known_answer"] = ka; report["licence"] = LICENCE
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
        rec = register(sid, np.asanyarray(img.dataobj), img.affine, meta, body, out)
        if prev: rec["superseded"] = prev.pop("superseded", []) + [prev]
        report["subjects"][sid] = rec
        man.write_text(json.dumps(report, indent=2, default=float) + "\n")
    subs_all = report["subjects"]
    say(f"\n{sum(v.get('passes', False) for v in subs_all.values())} of {len(subs_all)} subjects in the manifest pass every gate")

if __name__ == "__main__": main()
