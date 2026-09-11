#!/usr/bin/env python3
"""What is actually in the UT-EndoMRI cohort, in numbers.

The endometriosis caveat travels with every artefact derived from this dataset, and until now
it travelled as a sentence: *"an endometriosis cohort, NOT a typical-anatomy reference."* True,
and unfalsifiable as written. This turns it into figures a reader can check, so that anyone
using these organs knows how far from typical they are and in which direction.

It also exists because the extraction's own gates cannot answer this. `extract_ut_endomri.py`
says so in its header -- `volume_in_bound` is a COARSE typo catcher at 10-1000 mL, not a
plausibility check -- but the manifest did not, and on this cohort that gate passes 91 of 91
including a 765.7 mL uterus. A gate that cannot fail has not passed, and this script is where
the anatomical reading goes instead of being smuggled into a gate that was never meant to carry
it.

REFERENCE RANGES, stated before the data are read so they cannot be fitted to it. These are
textbook adult non-gravid figures and they are the comparison, not a gate:
  uterus  ~30-120 mL
  ovary   ~3-10 mL per side, and an ovary above ~20 mL in this cohort is endometrioma-scale
Nothing here PASSES or FAILS on them. They are the axis the cohort is described against.

KNOWN ANSWER, checked before any cohort figure is printed: the volumes are recomputed from the
manifest's own per-piece lists and must reproduce the recorded totals. If a total and its pieces
disagree the manifest is inconsistent and no summary of it means anything.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
UTERUS_NORMAL = (30.0, 120.0)
OVARY_NORMAL = (3.0, 10.0)
ENDOMETRIOMA_SCALE = 20.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest",
                    default="data/derived/ut-endomri-organs-v2-world/manifest.json")
    ap.add_argument("--out", default="data/derived/ut-endomri-organs-v2-world/cohort.json")
    a = ap.parse_args()

    m = json.loads((ROOT / a.manifest).read_text())
    S = m["subjects"]
    print(f"{a.manifest}\n  {m['caveat']}\n")
    print(f"  md5 verified at extraction: {m.get('md5_verified', 'NOT RECORDED -- older manifest')}"
          f"   (independently re-verified 2026-09-11: PASS)")

    # ---- known answer: the ovary pieces must sum to what is recorded --------------
    bad = []
    for sid, v in S.items():
        ov = v["organs"].get("ov")
        if not ov or not ov.get("pieces_ml"):
            continue
        if ov.get("count") != len(ov["pieces_ml"]):
            bad.append((sid, ov.get("count"), len(ov["pieces_ml"])))
    print(f"\nKNOWN ANSWER: every ovary `count` equals the length of its `pieces_ml`.")
    print(f"  mismatches: {len(bad)}  {'PASS' if not bad else 'FAIL ' + str(bad[:4])}")
    if bad:
        sys.exit("the manifest is internally inconsistent; no summary of it is meaningful")

    ut = np.array([v["organs"]["ut"]["volume_ml"] for v in S.values()
                   if "ut" in v["organs"] and v["organs"]["ut"].get("volume_ml") is not None])
    ovs = [v["organs"]["ov"] for v in S.values() if "ov" in v["organs"]]
    pieces = np.array([p for o in ovs for p in (o.get("pieces_ml") or [])])
    counts = [o.get("count", 0) for o in ovs]
    dice = [d for v in S.values() if "ut" in v["organs"]
            for d in v["organs"]["ut"]["inter_rater_dice"].values()]

    lo, hi = UTERUS_NORMAL
    print(f"\n{len(S)} subjects\n")
    print(f"UTERUS -- labelled in {len(ut)}, absent in {len(S) - len(ut)}")
    print(f"  median {np.median(ut):6.1f} mL   IQR {np.percentile(ut, 25):.1f}-{np.percentile(ut, 75):.1f}"
          f"   range {ut.min():.1f}-{ut.max():.1f}")
    print(f"  against a {lo:.0f}-{hi:.0f} mL adult non-gravid reference:")
    print(f"    above {hi:.0f} mL: {(ut > hi).sum():3d}  ({(ut > hi).mean():.0%} of those labelled)")
    print(f"    above {2*hi:.0f} mL: {(ut > 2*hi).sum():3d}  ({(ut > 2*hi).mean():.0%})")
    print(f"    below {lo:.0f} mL: {(ut < lo).sum():3d}  ({(ut < lo).mean():.0%})")
    print(f"  the largest is {ut.max():.1f} mL -- {ut.max()/hi:.1f}x the top of the reference range")

    olo, ohi = OVARY_NORMAL
    print(f"\nOVARY -- labelled in {len(ovs)} subjects, absent in {len(S) - len(ovs)}")
    print(f"  pieces per subject 0/1/2/3+: {counts.count(0)}/{counts.count(1)}/{counts.count(2)}/"
          f"{sum(c >= 3 for c in counts)}")
    print(f"  {counts.count(1)} of {len(ovs)} ({counts.count(1)/len(ovs):.0%}) have a SINGLE ovary "
          f"segmented. The dataset's own README names surgical resection as a reason.")
    print(f"  {len(pieces)} pieces: median {np.median(pieces):5.2f} mL   range {pieces.min():.2f}-{pieces.max():.2f}")
    print(f"  against a {olo:.0f}-{ohi:.0f} mL reference: above {ohi:.0f} mL "
          f"{(pieces > ohi).sum()} ({(pieces > ohi).mean():.0%}); "
          f"above {ENDOMETRIOMA_SCALE:.0f} mL (endometrioma-scale) "
          f"{(pieces > ENDOMETRIOMA_SCALE).sum()} ({(pieces > ENDOMETRIOMA_SCALE).mean():.0%})")

    if dice:
        d = np.array(dice)
        print(f"\nINTER-RATER AGREEMENT on the uterus (D1 only -- D2 has one rater)")
        print(f"  {len(d)} rater pairs: median Dice {np.median(d):.3f}   range {d.min():.3f}-{d.max():.3f}")
        print(f"  below 0.70: {(d < 0.70).sum()} pairs ({(d < 0.70).mean():.0%})")
        print(f"  the worst pair agrees at {d.min():.3f}, which is most of a disagreement about "
              f"where the organ is.")

    flagged = sum(1 for v in S.values() if v["flags"])
    print(f"\nSUBJECTS CARRYING EXTRACTION FLAGS: {flagged}/{len(S)} ({flagged/len(S):.0%})")

    print(f"\nWHAT THE CAVEAT MEANS, in numbers rather than as a sentence:")
    print(f"  {(ut > hi).mean():.0%} of labelled uteri are above the adult non-gravid reference and the")
    print(f"  largest is {ut.max()/hi:.1f}x its top; {(pieces > ENDOMETRIOMA_SCALE).mean():.0%} of ovary pieces are "
          f"endometrioma-scale; and")
    print(f"  {counts.count(1)/len(ovs):.0%} of subjects with an ovary label have only one. Any organ taken from")
    print(f"  here is pathology-selected, and a mean over this cohort is not a typical anatomy.")

    p = ROOT / a.out
    p.write_text(json.dumps({
        "caveat": m["caveat"],
        "reference_ranges_stated_in_advance": {
            "uterus_ml": list(UTERUS_NORMAL), "ovary_ml": list(OVARY_NORMAL),
            "endometrioma_scale_ml": ENDOMETRIOMA_SCALE,
            "note": "textbook adult non-gravid figures; a comparison axis, never a gate"},
        "n_subjects": len(S),
        "uterus": {"n_labelled": int(len(ut)), "median_ml": float(np.median(ut)),
                   "iqr_ml": [float(np.percentile(ut, 25)), float(np.percentile(ut, 75))],
                   "range_ml": [float(ut.min()), float(ut.max())],
                   "frac_above_reference": float((ut > hi).mean())},
        "ovary": {"n_subjects_labelled": len(ovs), "n_pieces": int(len(pieces)),
                  "median_ml": float(np.median(pieces)),
                  "range_ml": [float(pieces.min()), float(pieces.max())],
                  "frac_endometrioma_scale": float((pieces > ENDOMETRIOMA_SCALE).mean()),
                  "frac_subjects_unilateral": float(counts.count(1) / len(ovs))},
        "inter_rater_dice_uterus": ({"n_pairs": len(dice), "median": float(np.median(dice)),
                                     "min": float(min(dice)),
                                     "frac_below_0.70": float((np.array(dice) < 0.70).mean())}
                                    if dice else None),
        "subjects_with_flags": flagged}, indent=2) + "\n")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
