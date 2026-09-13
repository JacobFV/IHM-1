#!/usr/bin/env python3
"""Ribcage against ribcage: is gate d's failure anatomical?

Four attempts to compare this body's chest with the female trunks failed, every one of them on
SELECTION rather than on the measure: bone against skin, arms in the band, a lateral clip that cut
the torso at its widest, and a nearest-bone classifier that reached only 85.8% at the humerus
against a pre-set 90%. AP/ML itself passed its known answers every time.

This compares the one structure both bodies have unambiguously: **the ribcage**. No skin, no soft
tissue, no arms, no clipping. `data/derived/female-torso-totalsegmentator-v1/<subject>/meshes/`
carries 25 rib and sternum meshes per subject -- the same labels the registration built its 34 bone
correspondences from.

FRAMES. The female meshes are in scan frame and this body is in the atlas frame, so the raw axes
are not comparable. Applying the registration's own transform maps the female ribs into the atlas
frame, where +z is anterior for both. That transform is a SIMILARITY -- rotation, uniform scale,
translation -- so it makes the axes agree while leaving AP/ML, a ratio of two lengths, untouched by
the scale it chose.

KNOWN ANSWERS, all fixed before any comparison:
  1. AP/ML is invariant to a uniform scale, to 1e-9.
  2. The registration records this body's ribs as 97.3-99.9% contained inside the mapped female
     trunk, so the mapped female ribcage must be COMPARABLE IN SIZE to this body's -- within 25%
     on each axis. A gross mismatch would mean the mapping or the meshes are not what I think they
     are. This can fail, and it is the check that matters.
"""
import json, gzip, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
SRC = ROOT / "data/derived/female-torso-totalsegmentator-v1"


def read_obj(p):
    return np.asarray([ln.split()[1:4] for ln in open(p) if ln.startswith("v ")], float)


def ap_ml(P):
    ap = np.percentile(P[:, 2], 95) - np.percentile(P[:, 2], 5)
    ml = np.percentile(P[:, 0], 95) - np.percentile(P[:, 0], 5)
    return ap / ml, ap * 1e3, ml * 1e3


def main():
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    mine = []
    for e in ents:
        n = e.get("name", "").lower()
        if e.get("reference_geometry") and ("rib" in n or "sternum" in n or "manubrium" in n
                                            or "xiphoid" in n):
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            mine.append(np.asarray(g["positions"], float).reshape(-1, 3))
    MINE = np.vstack(mine)
    r_mine, ap_mine, ml_mine = ap_ml(MINE)
    print(f"this body's ribcage: {len(MINE):,} vertices  AP {ap_mine:.1f}  ML {ml_mine:.1f}  "
          f"AP/ML {r_mine:.3f}")

    print("\nKNOWN ANSWER 1: AP/ML invariant to a uniform scale")
    chk, _, _ = ap_ml(MINE * 1.5)
    ok1 = abs(chk - r_mine) < 1e-9
    print(f"  scaled 1.5x -> {chk:.6f} against {r_mine:.6f}   {'PASS' if ok1 else 'FAIL'}")
    if not ok1:
        sys.exit("AP/ML is not scale-free")

    print(f"\n{'subject':10s} {'ribs':>8s} {'AP mm':>7s} {'ML mm':>7s} {'AP/ML':>7s} {'size vs this body':>18s}")
    vals, ok2 = [], True
    for sid, d in REG.items():
        md = SRC / sid / "meshes"
        if not md.exists():
            print(f"{sid:10s} no meshes"); continue
        M = np.array(json.loads((ROOT / "data/derived" / d / "manifest.json").read_text())["transform"])
        pts = []
        for f in sorted(md.glob("rib_*.obj")) + sorted(md.glob("sternum*.obj")):
            V = read_obj(f)
            pts.append((np.c_[V, np.ones(len(V))] @ M.T)[:, :3])
        if not pts:
            print(f"{sid:10s} no rib meshes"); continue
        P = np.vstack(pts)
        rr, a, m = ap_ml(P); vals.append(rr)
        sz = f"AP {a/ap_mine:.2f}x ML {m/ml_mine:.2f}x"
        bad = not (0.75 < a / ap_mine < 1.25 and 0.75 < m / ml_mine < 1.25)
        ok2 &= not bad
        print(f"{sid:10s} {len(P):8,} {a:7.1f} {m:7.1f} {rr:7.3f} {sz:>18s}{'  <- OUT OF RANGE' if bad else ''}")

    print(f"\nKNOWN ANSWER 2: the mapped female ribcage must be within 25% of this body's on each "
          f"axis\n  (this body's ribs are 97.3-99.9% contained inside the mapped trunk)   "
          f"{'PASS' if ok2 else 'FAIL'}")
    if not ok2:
        sys.exit("the mapped ribcages are not comparable in size; the mapping or the meshes are "
                 "not what this assumes")

    v = np.array(vals)
    print(f"\n  female ribcages AP/ML: mean {v.mean():.3f}, range {v.min():.3f}-{v.max():.3f}")
    print(f"  this body's ribcage  : {r_mine:.3f}  ->  {100*(r_mine/v.mean()-1):+.1f}% relative")
    diff = abs(r_mine / v.mean() - 1)
    print("\n  VERDICT: " + (
        f"the ribcage proportions MATCH within {100*diff:.0f}%. A single uniform scale can fit "
        f"both,\n  so gate d is NOT anatomical in this sense and its mechanism remains unknown."
        if diff < 0.10 else
        f"the ribcage proportions DIFFER by {100*(r_mine/v.mean()-1):+.0f}%. A similarity transform "
        f"has one\n  scale and cannot match a chest that is relatively deeper on one body than the "
        f"other --\n  gate d is UNREACHABLE by this registration rather than merely unmet."))


if __name__ == "__main__":
    main()
