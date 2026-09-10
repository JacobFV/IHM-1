"""choose TotalSegmentator subjects by what their scan ACTUALLY contains, with no download.

the archive's `study_type` names the clinical exam, and the member is a CROP of it:
s1218, 'ct neck-thorax-abdomen-pelvis', is 45 slices at 1.5 mm, a 68 mm slab. a
breast spans ~150-200 mm, so selecting by study_type selects slabs.

what a crop contains is visible in the zip's central directory, already cached by
scripts/fetch_totalsegmentator_subjects.py: an all-zero label mask compresses to a
few hundred bytes and a present structure to far more. so a subject whose 24 ribs,
sternum and both clavicles are all above the threshold has the chest in its scan.

GATE (known answer): on a subject already on disk, the size rule must classify
every one of its 39 registration labels exactly as the labels' own voxel counts do.

LIMIT, and why a second gate exists downstream: the female-wide size histogram has
no clean gap at the threshold, so a rib present only as a clipped fragment passes
here. coverage is PROVEN per subject by build_female_torso_from_totalsegmentator.py
(no chest label may touch the first or last slice; the breast may touch no face).
"""
import argparse, csv, io, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/anatomy/totalsegmentator"
LABELS = ([f"rib_{s}_{i}" for s in ("left", "right") for i in range(1, 13)]
          + ["sternum", "clavicula_left", "clavicula_right"] + [f"vertebrae_T{i}" for i in range(1, 13)])
CHEST = [f"rib_{s}_{i}" for s in ("left", "right") for i in range(1, 13)] + ["sternum", "clavicula_left", "clavicula_right"]

def main():
    import nibabel as nib
    ap = argparse.ArgumentParser(); ap.add_argument("--calibrate-on", default="s1218")
    ap.add_argument("--target-age", type=float, default=40.0); a = ap.parse_args()
    cd = json.loads((RAW / "central_directory.json").read_text())["members"]
    meta = {r["image_id"]: r for r in csv.DictReader(io.open(RAW / "meta.csv", encoding="utf-8-sig"), delimiter=";")}
    sid = a.calibrate_on
    truth = {l: int((np.asanyarray(nib.load(str(RAW / sid / "segmentations" / f"{l}.nii.gz")).dataobj) > 0).sum()) for l in LABELS}
    csz = {l: cd[f"{sid}/segmentations/{l}.nii.gz"]["compressed"] for l in LABELS}
    empty = [csz[l] for l in LABELS if truth[l] == 0]; present = [csz[l] for l in LABELS if truth[l] > 0]
    if not empty or not present: raise SystemExit(f"{sid} has no empty or no present label; it cannot calibrate a threshold")
    thr = (max(empty) + min(present)) / 2
    agree = all((csz[l] > thr) == (truth[l] > 0) for l in LABELS)
    print(f"calibration on {sid}: empty {min(empty)}..{max(empty)} B, present {min(present)}..{max(present)} B, threshold {thr} B")
    print("GATE: size rule matches voxel counts on all 39 labels:", "PASS" if agree else "FAIL")
    if not agree: raise SystemExit("gate failed; not selecting anything")
    rows = []
    for s, r in meta.items():
        if r["gender"] != "f" or r["age"] in ("", "nan") or r["pathology"] != "no_pathology": continue
        if not all(f"{s}/segmentations/{l}.nii.gz" in cd for l in CHEST): continue
        if all(cd[f"{s}/segmentations/{l}.nii.gz"]["compressed"] > thr for l in CHEST):
            rows.append((s, float(r["age"])))
    rows.sort(key=lambda x: abs(x[1] - a.target_age))
    print(f"{len(rows)} female no_pathology subjects have all {len(CHEST)} chest labels above threshold")
    out = RAW / "coverage_selection.json"
    out.write_text(json.dumps(dict(threshold_bytes=thr, calibrated_on=sid, gate_pass=agree,
        rule="female; no_pathology; all 24 ribs + sternum + both clavicles above the empty-mask size; ages closest to %g" % a.target_age,
        candidates=[s for s, _ in rows]), indent=2) + "\n")
    print("wrote", out.relative_to(ROOT))

if __name__ == "__main__": main()
