"""female torso geometry from TotalSegmentator CTs: breast, skin envelope, and the
bones it will be registered by.

WHY. This body is male: 0 female-specific entities, 0 mammary glands in either sex
(docs/BODY_PARAMETERS.md). A female breast cannot simply be dropped inside it,
because the skin is male too -- a female torso needs its own envelope. Both come
from the same female CT: the `breasts` subtask (one class, 'breast') and the
`body` subtask (body, body_trunc, body_extremities, skin), both Apache-2.0 models,
run on CTs from the CC-BY-4.0 TotalSegmentator v2.0.1 dataset.

WHAT IT IS NOT. Each output is ONE clinical subject's anatomy as a segmentation
model drew it. 'breast' is a single undifferentiated soft-tissue label: no gland,
duct, nipple or areola. And the subject is not the skeleton's subject nor the
skin's; registering it into this body is a separate step (it needs the per-segment
registration, judged by containment -- see docs/SEGMENT_CONTACT_SURFACES.md).

GATES, each with an answer this script does not compute:
  bones-in-body      the SAME CT's own bones lie inside its own body mask (>= 0.99)
  breast-in-body     the breast lies inside the body mask grown by ONE voxel (>= 0.99).
                     AMENDED AFTER A RESULT, and said so: as first written this was the
                     strict body mask, and s0790 failed it at 0.9825. every one of its
                     6,127 outside voxels lies exactly 1 voxel (1.5 mm) from the body
                     mask, in 3,598 pieces of at most 18 voxels -- the breasts and body
                     models drawing the skin surface one voxel apart, not misplaced
                     anatomy. the strict value is still reported beside it.
  breast-off-ribs    the breast does not overlap any rib label (< 0.01)
  laterality         midline = the sternum's x; TotalSegmentator's OWN
                     clavicula_left must fall on the side this script calls left.
                     catches a flipped axis convention before it names a breast.
  volume             per-side breast volume within a COARSE sanity bound
                     (BREAST_ML_BOUND). a typo catcher, not a measurement.
  chest-whole        all 24 ribs, the sternum and both clavicles present, and NONE
                     touches the first or last slice along the superior-inferior
                     axis. the archive's study_type names the exam and the member is
                     a crop: s1218 'neck-thorax-abdomen-pelvis' is a 68 mm slab.
  breast-whole       the breast mask touches no face of the volume. a breast cut
                     by the field of view is a wrong breast, not a smaller one.
"""
import argparse, hashlib, json, subprocess
from pathlib import Path
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/anatomy/totalsegmentator"
OUT = ROOT / "data/derived/female-torso-totalsegmentator-v1"
TS = ROOT / ".venv-totalseg/bin/TotalSegmentator"
REGISTRATION_LABELS = ([f"rib_{s}_{i}" for s in ("left", "right") for i in range(1, 13)]
                       + ["sternum", "clavicula_left", "clavicula_right"]
                       + [f"vertebrae_T{i}" for i in range(1, 13)])
# STRAY FRAGMENTS. TotalSegmentator occasionally labels a speck far from the bone it names:
# s0970's "rib_left_4" is a 3,864-voxel rib plus 10- and 1-voxel specks 343 mm below it at the
# bottom of the scan; its "rib_right_12" carries a 25-voxel speck 202 mm away that drags the
# rib's centroid 10.0 mm. Those specks made s0970 fail chest-whole. But secondary pieces are
# not all strays: s0790's T2-T4 and s1159's T7 carry 77-224-voxel pieces 7.5-11.2 mm from the
# vertebra (inside its own 30-40 mm span), moving its centroid at most 1.2 mm -- bone split by a
# sub-voxel gap. So the rule is DISTANCE, not "largest piece": drop a piece only if its closest
# approach to the label's largest piece exceeds STRAY_GAP_MM. Across all eight subjects and 39
# labels, every piece KEPT lies within 11.2 mm of its bone and every piece DROPPED lies 73.5 mm
# or more from it (s0897 rib_right_8: 73.5-117.8 mm; s0970: 138-372 mm), so 50 mm sits in the
# empty interval 11.2-73.5 mm. (A first version of this comment said 11-201 mm, from six pieces;
# the full census narrowed it.) Every dropped piece is recorded in the manifest.
STRAY_GAP_MM = 50.0

# NOT a measured range: a bound wide enough that any real adult breast passes and a
# unit or laterality error (mm^3 read as mL, both sides in one) does not.
BREAST_ML_BOUND = (30.0, 3000.0)

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run_subtask(ct, task, dest, device):
    if dest.exists() and any(dest.glob("*.nii.gz")): return
    dest.mkdir(parents=True, exist_ok=True)
    cmd = [str(TS), "-i", str(ct), "-o", str(dest), "-ta", task, "--device", device]
    print("  $", " ".join(cmd[1:]), flush=True)
    subprocess.run(cmd, check=True)

def mesh(mask, affine):
    from skimage import measure
    padded = np.pad(mask.astype(np.uint8), 1)
    v, f, _, _ = measure.marching_cubes(padded, 0.5)
    v = v - 1
    world = v @ affine[:3, :3].T + affine[:3, 3]
    return world / 1000.0, f            # NIfTI world is mm; this repository is metres

def write_obj(path, v, f, header):
    with open(path, "w") as h:
        h.write(f"# {header}\n")
        for x in v: h.write("v %.7g %.7g %.7g\n" % tuple(x))
        for t in f: h.write("f %d %d %d\n" % (t[0] + 1, t[1] + 1, t[2] + 1))

def main():
    import nibabel as nib
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--device", default="gpu", choices=("gpu", "cpu"))
    ap.add_argument("--subjects", nargs="*", default=None,
                    help="run only these subjects, taking their metadata from the archive's meta.csv "
                         "(for re-gating a subject the current fetch manifest no longer lists)")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    OUT_ = a.out
    fetched = json.loads((RAW / "manifest.json").read_text())
    if a.subjects:
        import csv, io
        meta = {r["image_id"]: r for r in csv.DictReader(io.open(RAW / "meta.csv", encoding="utf-8-sig"), delimiter=";")}
        fetched = dict(fetched, subjects={sid: dict(meta=meta[sid]) for sid in a.subjects})
    OUT_.mkdir(parents=True, exist_ok=True)
    report = dict(schema="ihm.female-torso-totalsegmentator.v1", source=fetched["citation"],
                  dataset_licence=fetched["licence"], model_licence="Apache-2.0 (breasts, body subtasks)",
                  selection_rule=fetched["selection_rule"], breast_ml_bound=BREAST_ML_BOUND,
                  breast_ml_bound_basis="sanity bound, not a measurement", subjects={})
    for sid, meta in fetched["subjects"].items():
        print(f"\n=== {sid}  age {meta['meta']['age']}  {meta['meta']['study_type']}", flush=True)
        ct = RAW / sid / "ct.nii.gz"
        if not ct.exists(): raise SystemExit(f"{ct} missing: the fetch did not deliver this subject's CT")
        run_subtask(ct, "breasts", OUT_ / sid / "breasts", a.device)
        run_subtask(ct, "body", OUT_ / sid / "body", a.device)
        img = nib.load(str(ct)); affine = img.affine; voxel_ml = abs(np.linalg.det(affine[:3, :3])) / 1000.0
        load = lambda p: np.asanyarray(nib.load(str(p)).dataobj) > 0
        body = load(OUT_ / sid / "body" / "body.nii.gz")
        breast = load(OUT_ / sid / "breasts" / "breast.nii.gz")
        seg = RAW / sid / "segmentations"
        from scipy.spatial import cKDTree
        world_mm = lambda ijk: ijk @ affine[:3, :3].T + affine[:3, 3]
        strays = {}
        def clean(n, m):
            if not m.any(): return m
            lab, k = ndimage.label(m, structure=np.ones((3, 3, 3), bool))
            if k < 2: return m
            sizes = np.bincount(lab.ravel())[1:]; main = int(np.argmax(sizes)) + 1
            tree = cKDTree(world_mm(np.argwhere(lab == main)))
            for j in range(1, k + 1):
                if j == main: continue
                gap = float(tree.query(world_mm(np.argwhere(lab == j)))[0].min())
                if gap > STRAY_GAP_MM:
                    m = m & (lab != j)
                    strays.setdefault(n, []).append(dict(voxels=int(sizes[j - 1]), gap_to_main_mm=round(gap, 1)))
            return m
        bones = {n: clean(n, load(seg / f"{n}.nii.gz")) for n in REGISTRATION_LABELS if (seg / f"{n}.nii.gz").exists()}
        present = {n: m for n, m in bones.items() if m.any()}
        missing = sorted(set(REGISTRATION_LABELS) - set(present))
        rec = dict(meta=meta["meta"], ct_sha256=sha(ct), voxel_ml=voxel_ml, stray_pieces_dropped=strays,
                   registration_labels_present=len(present), registration_labels_missing=missing)
        # --- gates
        all_bone = np.zeros_like(body)
        for m in present.values(): all_bone |= m
        rec["gate_bones_in_body"] = float(body[all_bone].mean()) if all_bone.any() else None
        body_1 = ndimage.binary_dilation(body, iterations=1)
        rec["gate_breast_in_body_strict"] = float(body[breast].mean()) if breast.any() else None
        rec["gate_breast_in_body"] = float(body_1[breast].mean()) if breast.any() else None
        # the MESH never pierces the skin: the breast is clipped to the body mask, and the
        # volume clipped away is reported so the operation is visible
        rec["breast_ml_clipped_to_body"] = float((breast & ~body).sum() * voxel_ml)
        breast = breast & body
        ribs = np.zeros_like(body)
        for n, m in present.items():
            if n.startswith("rib_"): ribs |= m
        rec["gate_breast_on_ribs"] = float(ribs[breast].mean()) if breast.any() else None
        codes = nib.aff2axcodes(affine); si = next(i for i, c in enumerate(codes) if c in "SI")
        def touches_si(m):
            ix = np.argwhere(m)[:, si]; return bool(ix.min() == 0 or ix.max() == m.shape[si] - 1)
        chest = [f"rib_{x}_{i}" for x in ("left", "right") for i in range(1, 13)] + ["sternum", "clavicula_left", "clavicula_right"]
        rec["gate_chest_whole"] = bool(all(n in present for n in chest) and not any(touches_si(present[n]) for n in chest))
        def touches_face(m):
            ijk_ = np.argwhere(m); return bool((ijk_.min(0) == 0).any() or (ijk_.max(0) == np.array(m.shape) - 1).any())
        rec["gate_breast_whole"] = bool(breast.any() and not touches_face(breast))
        rec["volume_extent_mm"] = [float(n * z) for n, z in zip(body.shape, img.header.get_zooms()[:3])]
        def world_x(m):
            ijk = np.argwhere(m); return (ijk @ affine[:3, :3].T + affine[:3, 3])[:, 0]
        if "sternum" not in present: raise SystemExit(f"{sid}: no sternum label, so no midline; refusing to name sides")
        midline = float(np.median(world_x(present["sternum"])))
        # do not trust a remembered axis convention: decide it
        # from the data: whichever side TotalSegmentator's clavicula_left falls on IS left.
        left_sign = np.sign(np.median(world_x(present["clavicula_left"])) - midline) if "clavicula_left" in present else None
        right_sign = np.sign(np.median(world_x(present["clavicula_right"])) - midline) if "clavicula_right" in present else None
        rec["gate_laterality"] = bool(left_sign is not None and right_sign is not None and left_sign == -right_sign != 0)
        if not rec["gate_laterality"]: raise SystemExit(f"{sid}: clavicles do not straddle the sternum; laterality undefined")
        ijk = np.argwhere(breast); bx = (ijk @ affine[:3, :3].T + affine[:3, 3])[:, 0]
        on_left = np.sign(bx - midline) == left_sign
        rec["breast_ml"] = dict(left=float(on_left.sum() * voxel_ml), right=float((~on_left).sum() * voxel_ml))
        rec["gate_volume"] = all(BREAST_ML_BOUND[0] <= v <= BREAST_ML_BOUND[1] for v in rec["breast_ml"].values())
        print(f"  bones in body {rec['gate_bones_in_body']:.4f} | breast in body {rec['gate_breast_in_body']:.4f} "
              f"(strict {rec['gate_breast_in_body_strict']:.4f}, {rec['breast_ml_clipped_to_body']:.1f} mL clipped) | "
              f"breast on ribs {rec['gate_breast_on_ribs']:.4f} | laterality {rec['gate_laterality']} | "
              f"breast L {rec['breast_ml']['left']:.0f} mL R {rec['breast_ml']['right']:.0f} mL | "
              f"registration labels {len(present)}/{len(REGISTRATION_LABELS)} | "
              f"chest whole {rec['gate_chest_whole']} | breast whole {rec['gate_breast_whole']} | "
              f"extent {' x '.join('%.0f' % e for e in rec['volume_extent_mm'])} mm", flush=True)
        rec["passes"] = bool(rec["gate_bones_in_body"] is not None and rec["gate_bones_in_body"] >= .99
                             and rec["gate_breast_in_body"] >= .99 and rec["gate_breast_on_ribs"] < .01
                             and rec["gate_laterality"] and rec["gate_volume"]
                             and rec["gate_chest_whole"] and rec["gate_breast_whole"])
        # --- meshes, in the CT's own world frame (metres); registration is a later step
        d = OUT_ / sid / "meshes"; d.mkdir(parents=True, exist_ok=True)
        for side, sel in (("left", on_left), ("right", ~on_left)):
            m = np.zeros_like(breast); m[tuple(ijk[sel].T)] = True
            v, f = mesh(m, affine); write_obj(d / f"breast_{side}.obj", v, f, f"{sid} breast_{side}, CT world frame, metres")
        skin = OUT_ / sid / "body" / "skin.nii.gz"
        if skin.exists():
            v, f = mesh(load(skin), affine); write_obj(d / "skin.obj", v, f, f"{sid} skin, CT world frame, metres")
        for n, m in present.items():
            v, f = mesh(m, affine); write_obj(d / f"{n}.obj", v, f, f"{sid} {n}, CT world frame, metres")
        report["subjects"][sid] = rec
        (OUT_ / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    ok = [s for s, r in report["subjects"].items() if r["passes"]]
    print(f"\n{len(ok)}/{len(report['subjects'])} subjects pass every gate: {ok}")

if __name__ == "__main__": main()
