"""Where the registered ovaries land, which the pelvic gates only ever reported.

87 subjects are registered; 66 carry an ovary mesh (57 one piece, 9 two). The gates judged the
uterus and reported the ovaries as a share inside the pelvic ring, which says nothing about
whether they sit where ovaries sit: lateral, against the pelvic side wall, in the ovarian fossa,
further from the midline than the uterus.

KNOWN ANSWERS, checked before any ovary number is read:
  1. this body's left hip bone lies on the left of its right hip bone -- that is what fixes the
     lateral axis and its sign, and the registration's own laterality gate is stated in those terms;
  2. the uterus is a midline organ: the median |lateral offset| over all 87 registered uteri must
     be small, and smaller than the ovaries'. If the uterus reads lateral, the frame is wrong and
     nothing below means anything;
  3. the two pieces of a two-piece subject fall on OPPOSITE sides. That is anatomy, not a fit:
     a woman has one ovary per side, so a rule that puts both on one side has mapped something
     wrong.

Caveat carried from the source: UT-EndoMRI is an endometriosis cohort, in which ovaries are
displaced, adherent and enlarged by endometriomas; this measures where the labels land, not where
a healthy ovary sits.
"""
import importlib.util, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
def load(name, path):
    s = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
R = load("rpo", ROOT / "scripts/register_pelvic_organs.py")
Q = load("q", ROOT / "scripts/score_uterus_lowest_quartile.py")
REG = ROOT / "data/derived/ut-endomri-registered-v1"

def main():
    body = R.body_meshes(); up, up_sign = Q.vertical()
    left, right = body["left hip bone"][0], body["right hip bone"][0]
    d = left.mean(0) - right.mean(0); d[up] = 0
    lat = int(np.argmax(np.abs(d))); lat_sign = float(np.sign(d[lat]))
    assert lat != up and abs(d[lat]) > 0.05, f"hip-to-hip offset {d} gives no clear lateral axis"
    print(f"known answer 1: lateral axis {lat}, +{lat_sign:.0f} toward this body's LEFT "
          f"(hip centroids {1e3*abs(d[lat]):.0f} mm apart)")
    hips = np.vstack([left, right]); hip_tree = cKDTree(hips)
    floor = (up_sign * hips[:, up]).min()
    mid = (left.mean(0)[lat] + right.mean(0)[lat]) / 2

    man = json.loads((REG / "manifest.json").read_text())["subjects"]
    ut_off, rows = [], []
    for sid, rec in sorted(man.items()):
        u = REG / sid / "uterus.obj"
        if u.exists():
            V, _ = R.read_obj(u); ut_off.append(abs(1e3 * (V.mean(0)[lat] - mid)))
        for ov in rec.get("ovaries", []):
            p = REG / sid / f"ovary_{ov['piece']}.obj"
            if not p.exists(): continue
            V, F = R.read_obj(p); c = V.mean(0)
            rows.append(dict(subject=sid, piece=ov["piece"],
                             lateral_mm=1e3 * lat_sign * (c[lat] - mid),
                             wall_mm=1e3 * hip_tree.query(c[None])[0][0],
                             height_mm=1e3 * (up_sign * c[up] - floor),
                             inside=ov["inside_ring"]))
    ut = float(np.median(ut_off))
    ov_lat = np.array([abs(r["lateral_mm"]) for r in rows])
    print(f"known answer 2: uterus median |lateral offset| {ut:.1f} mm over {len(ut_off)} subjects, "
          f"ovaries {np.median(ov_lat):.1f} mm over {len(rows)} pieces -> "
          f"{'PASS' if ut < np.median(ov_lat) else 'FAIL: the frame is wrong'}")

    pairs = {}
    for r in rows: pairs.setdefault(r["subject"], []).append(r)
    two = {k: v for k, v in pairs.items() if len(v) == 2}
    opp = [k for k, v in two.items() if v[0]["lateral_mm"] * v[1]["lateral_mm"] < 0]
    print(f"known answer 3: of {len(two)} two-piece subjects, {len(opp)} put their pieces on "
          f"opposite sides -> {'PASS' if len(opp) == len(two) else 'FAIL'}"
          + ("" if len(opp) == len(two) else f"  same-side: {sorted(set(two) - set(opp))}"))

    print(f"\n{len(rows)} ovary pieces over {len(pairs)} subjects")
    for name, key in (("lateral offset from the midline", "lateral_mm"),
                      ("distance to the nearest hip bone", "wall_mm"),
                      ("height above the ring floor", "height_mm")):
        v = np.array([r[key] for r in rows])
        print(f"  {name:34s} median {np.median(v):6.1f} mm   10-90% {np.percentile(v,10):6.1f} to {np.percentile(v,90):6.1f}")
    inside = np.array([r["inside"] for r in rows])
    print(f"  inside the pelvic ring                 median {np.median(inside):.3f}, "
          f"{int((inside >= 0.99).sum())} of {len(inside)} pieces at >= 0.99")
    side = np.array([r["lateral_mm"] for r in rows])
    print(f"  {int((side > 0).sum())} pieces on this body's left, {int((side < 0).sum())} on its right")
    (REG / "ovaries.json").write_text(json.dumps(dict(
        caveat="UT-EndoMRI is an endometriosis cohort: ovaries are displaced, adherent and enlarged by endometriomas",
        lateral_axis=lat, lateral_sign_left=lat_sign, uterus_median_abs_lateral_mm=ut, pieces=rows), indent=2) + "\n")
    print(f"\n{REG / 'ovaries.json'}")

if __name__ == "__main__": main()
