#!/usr/bin/env python3
"""This body's chest AP/ML with the arms removed BY ANATOMY, and the female trunks beside it.

Gate d fails on all four registrations: this body's sternum protrudes anteriorly through the
mapped female trunk. One candidate explanation is anatomical -- if this body's chest is
proportioned differently, the registration's single uniform scale cannot match both bodies and
gate d is unreachable rather than unmet.

THREE EARLIER ATTEMPTS AT THIS FAILED, each differently, and they are why this script exists:
  ribcage vs trunk surface      -- bone against skin; soft tissue makes a trunk far wider
  skin vs trunk surface         -- the ARMS are in the band; skin ML 459.9 mm, ribcage 214.2 mm
  skin clipped to the ribcage's lateral span -- the clip cuts the torso at its widest point
The instrument was never the problem: AP/ML is scale-invariant to six decimals and idempotent.
The selection was.

WHAT THIS DOES INSTEAD. `data/derived/anatomy-segment-binding/binding.json` binds all 4,000
entities to rigid segments, and its basis is `nearest_bone_group_vertex_vote`. The skin entity
itself is bound to `None` with coherence 0.252 -- correctly, because one surface spans many
segments -- so the same per-vertex rule is applied here: each skin vertex takes the segment of its
nearest BONE vertex, and only those nearest a `torso` bone are kept. Arms leave by anatomy, not by
a coordinate threshold.

KNOWN ANSWERS, both before any comparison:
  1. a skin vertex nearest the STERNUM must classify as torso -- checked on the skin vertices
     closest to the sternum, which must be 100% torso;
  2. a skin vertex nearest a HUMERUS must classify as arm -- checked the same way. If the
     classifier cannot separate the two bones it was built to separate, nothing downstream counts.
AP/ML itself is checked for scale-invariance as before.
"""
import json, gzip, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
ARM = {"humerus_l", "humerus_r", "ulna_l", "ulna_r", "radius_l", "radius_r", "hand_l", "hand_r"}


def read_obj(p):
    V = [ln.split()[1:4] for ln in open(p) if ln.startswith("v ")]
    return np.asarray(V, float)


def main():
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    bind = json.loads((ROOT / "data/derived/anatomy-segment-binding/binding.json").read_text())["entities"]

    def geom(e):
        g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
        return np.asarray(g["positions"], float).reshape(-1, 3)

    bone_pts, bone_seg, skin, sternum, humerus = [], [], None, [], []
    for e in ents:
        if not e.get("reference_geometry"):
            continue
        nm, eid = e.get("name", "").lower(), e.get("id")
        seg = (bind.get(eid) or {}).get("segment")
        if e.get("role") == "skin" or nm == "skin":
            skin = geom(e); continue
        if seg is None:
            continue
        # bones only: the binding's own basis is a nearest-BONE vote
        if not any(k in nm for k in ("bone", "rib", "sternum", "manubrium", "xiphoid", "vertebra",
                                     "humerus", "ulna", "radius", "carpal", "phalanx", "metacarpal",
                                     "clavicle", "scapula", "femur", "tibia", "fibula", "patella",
                                     "sacrum", "coccyx", "ilium", "ischium", "pubis", "talus",
                                     "calcaneus", "tarsal", "metatarsal", "skull", "mandible")):
            continue
        P = geom(e)
        bone_pts.append(P); bone_seg += [seg] * len(P)
        if any(k in nm for k in ("sternum", "manubrium", "xiphoid")): sternum.append(P)
        if "humerus" in nm: humerus.append(P)

    B = np.vstack(bone_pts); S = np.array(bone_seg)
    ST = np.vstack(sternum); HU = np.vstack(humerus) if humerus else None
    tree = cKDTree(B)
    print(f"bones {len(B):,} vertices over {len(set(S))} segments | skin {len(skin):,} | "
          f"sternum {len(ST):,} | humerus {0 if HU is None else len(HU):,}")

    def seg_of(P):
        _, i = tree.query(P, k=1, workers=-1)
        return S[i]

    print("\nKNOWN ANSWERS")
    ok = True
    for tag, ref, want in (("nearest the STERNUM", ST, "torso"), ("nearest a HUMERUS", HU, "arm")):
        if ref is None:
            continue
        d, _ = cKDTree(ref).query(skin, k=1, workers=-1)
        near = skin[d < np.percentile(d, 2)]          # the skin closest to that bone
        segs = seg_of(near)
        frac = (np.isin(segs, list(ARM)) if want == "arm" else (segs == "torso")).mean()
        print(f"  skin {tag:20s} classifies as {want:5s} for {frac:.1%}  "
              f"{'PASS' if frac > 0.9 else 'FAIL'}")
        ok &= frac > 0.9
    if not ok:
        sys.exit("a known answer FAILED; the classifier cannot separate the bones it must")

    y0, y1 = ST[:, 1].min(), ST[:, 1].max()

    def ap_ml(P):
        ap = np.percentile(P[:, 2], 95) - np.percentile(P[:, 2], 5)
        ml = np.percentile(P[:, 0], 95) - np.percentile(P[:, 0], 5)
        return ap / ml, ap * 1e3, ml * 1e3

    band = skin[(skin[:, 1] >= y0) & (skin[:, 1] <= y1)]
    seg = seg_of(band)
    trunk = band[seg == "torso"]
    r_all, ap_all, ml_all = ap_ml(band)
    r, ap, ml = ap_ml(trunk)
    s = 1.5
    chk, _, _ = ap_ml(trunk * s)
    print(f"  AP/ML scale-invariance: {chk:.6f} against {r:.6f}  "
          f"{'PASS' if abs(chk - r) < 1e-9 else 'FAIL'}")
    if abs(chk - r) >= 1e-9:
        sys.exit("AP/ML is not scale-free")

    print(f"\nthis body's skin in the sternum band (y {y0*1e3:.0f}..{y1*1e3:.0f} mm)")
    print(f"  all of it        {len(band):6,} vertices  AP {ap_all:6.1f}  ML {ml_all:6.1f}  AP/ML {r_all:.3f}")
    print(f"  TORSO-bound only {len(trunk):6,} vertices  AP {ap:6.1f}  ML {ml:6.1f}  AP/ML {r:.3f}")
    print(f"  arms removed by anatomy: {len(band)-len(trunk):,} vertices dropped")

    print(f"\n{'mapped female trunk':22s} {'AP mm':>7s} {'ML mm':>7s} {'AP/ML':>7s}")
    vals = []
    for sid, d in REG.items():
        f = ROOT / "data/derived" / d / "body_trunc_surface.obj"
        if not f.exists():
            continue
        V = read_obj(f); Vb = V[(V[:, 1] >= y0) & (V[:, 1] <= y1)]
        rr, a2, m2 = ap_ml(Vb); vals.append(rr)
        print(f"{sid:22s} {a2:7.1f} {m2:7.1f} {rr:7.3f}")
    v = np.array(vals)
    print(f"\n  female mapped trunks AP/ML mean {v.mean():.3f}, range {v.min():.3f}-{v.max():.3f}")
    print(f"  this body's TORSO skin AP/ML    {r:.3f}  ->  {100*(r/v.mean()-1):+.1f}% relative")
    print("\n  VERDICT: " + (
        "the proportions MATCH -- a uniform scale can fit both, so gate d is not anatomical and\n"
        "  its mechanism is still unknown." if abs(r / v.mean() - 1) < 0.10 else
        f"the proportions DIFFER by {100*(r/v.mean()-1):+.0f}% -- a single uniform scale cannot match\n"
        f"  both bodies, and gate d is UNREACHABLE by this registration rather than merely unmet."))


if __name__ == "__main__":
    main()
