"""is the female-torso sternum failure a LABELLING mismatch or a registration failure?

register_female_torso.py stopped at gate d: a quarter of this body's sternum lies outside
the mapped female trunk (0.767), and the sternum is the largest-residual bone in every fit
(23.8 mm centroid, 18.1 after ICP). the registration's own reading was that one similarity
cannot carry the anterior chest wall between people. that is an inference.

a cheaper explanation has to be ruled out first. this body's sternum is three entities --
manubrium, body of sternum, xiphoid process -- and the correspondence took their UNION.
TotalSegmentator's 'sternum' is one label drawn on CT, where the xiphoid is often
cartilaginous and faint. if the two 'sternums' cover different bone, the centroid pair is
biased along the sternum's length by construction.

the test: map TotalSegmentator's sternum through the fitted transform and ask, for each
part of this body's sternum, what share of its vertices lie within NEAR_MM of it.
  manubrium and body near, xiphoid far  -> the label omits the xiphoid; correspondence bias
  every part equally far                -> not a labelling mismatch; the thorax fit is the cause

GATE 1 (known answer): with the registration's OWN centroid definition, the saved ICP
transform must reproduce the sternum residual register_female_torso.py reported, 18.1 mm.
GATE 2: the mapped TotalSegmentator sternum must lie in this body's sternal region at all --
its centroid within 30 mm of this body's sternum centroid.

CENTROIDS AND SAMPLES ARE AREA-WEIGHTED, and the first run of this script shows why. It took
raw VERTEX means, and BodyParts3D meshes are very unevenly tessellated: the definition alone
moves this body's sternum centroid 14.3 mm (TotalSegmentator's marching-cubes mesh, which is
uniform, 1.5 mm). By vertex means the refined residual read 33.2 mm and gate 2 FAILED; by the
registration's own area-weighted definition it is 18.1 mm, exactly what the registration
reported. The threshold was not changed; the measurement was corrected to the definition the
number it is compared against was computed with.
"""
import gzip, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
REG = ROOT / "data/derived/female-torso-registered-v1/manifest.json"
TS_STERNUM = ROOT / "data/derived/female-torso-totalsegmentator-v1/s0790/meshes/sternum.obj"
ANATOMY = ROOT / "data/derived/canonical/anatomy.json"
PARTS = ("manubrium", "body of sternum", "xiphoid process")
NEAR_MM = 5.0

def read_obj(p):
    V, F = [], []
    for l in open(p):
        if l.startswith("v "): V.append([float(x) for x in l.split()[1:4]])
        elif l.startswith("f "): F.append([int(x.split("/")[0]) - 1 for x in l.split()[1:4]])
    return np.asarray(V), np.asarray(F, np.int64)

def areas(V, F):
    t = V[F]; return np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2

def area_centroid(V, F):
    a = areas(V, F); return (a[:, None] * V[F].mean(1)).sum(0) / a.sum()

def area_samples(V, F, n, seed=0):
    """uniform on SURFACE AREA, so a densely tessellated patch does not outvote a sparse one"""
    rng = np.random.default_rng(seed); t = V[F]; a = areas(V, F); k = rng.choice(len(F), n, p=a / a.sum())
    r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = t[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2]

def main():
    m = json.loads(REG.read_text())
    if "transform" not in m:
        raise SystemExit(f"{REG} carries no transform (keys {list(m)}); rerun register_female_torso.py first")
    T = np.asarray(m["transform"], float)
    apply = lambda M, P: P @ np.asarray(M)[:3, :3].T + np.asarray(M)[:3, 3]
    Vs, Fs = read_obj(TS_STERNUM)
    ents = {e["name"]: e for e in json.loads(ANATOMY.read_text())["entities"]
            if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-")}
    part = {}
    for n in PARTS:
        g = json.loads(gzip.decompress((ROOT / ents[n]["reference_geometry"]["path"]).read_bytes()))
        part[n] = (np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3))
    off = np.cumsum([0] + [len(v) for v, _ in list(part.values())[:-1]])
    Vb = np.vstack([v for v, _ in part.values()]); Fb = np.vstack([f + o for (_, f), o in zip(part.values(), off)])
    cb = area_centroid(Vb, Fb)
    print(f"transform used: {m.get('used', '?')}")
    g1 = 1000 * np.linalg.norm(apply(m["icp_refined_transform"], area_centroid(Vs, Fs)[None])[0] - cb)
    print(f"GATE 1: refined sternum residual by the registration's definition {g1:.1f} mm (registration reported 18.1) -> "
          + ("PASS" if abs(g1 - 18.1) < 0.1 else "FAIL -- this is not the registration's measurement"))
    if abs(g1 - 18.1) >= 0.1: sys.exit(1)
    ts = area_samples(apply(T, Vs), Fs, 20000)
    offc = 1000 * np.linalg.norm(apply(T, area_centroid(Vs, Fs)[None])[0] - cb)
    print(f"GATE 2: mapped TS sternum centroid {offc:.1f} mm from this body's sternum centroid (must be < 30) -> "
          + ("PASS" if offc < 30 else "FAIL -- the transform does not put the sternum near the sternum"))
    if offc >= 30: sys.exit(1)
    body_s = area_samples(Vb, Fb, 20000)
    axis = np.linalg.eigh(np.cov((body_s - cb).T))[1][:, -1]          # the sternum's own long axis
    proj = lambda P: (P - cb) @ axis
    print(f"\nalong this body's sternal axis (mm, 0 = its area centroid):")
    print(f"  TotalSegmentator sternum, mapped   {1000*proj(ts).min():+7.1f} .. {1000*proj(ts).max():+7.1f}  "
          f"(length {1000*np.ptp(proj(ts)):.0f} mm)")
    tree = cKDTree(ts)
    print(f"\n{'part':18s} {'extent along axis (mm)':>26s} {'median to TS':>13s} {'within %g mm' % NEAR_MM:>12s}")
    rows = {}
    for n, (V, F) in part.items():
        P = area_samples(V, F, 4000); d = 1000 * tree.query(P)[0]; pr = 1000 * proj(P)
        rows[n] = float(np.mean(d <= NEAR_MM))
        print(f"{n:18s} {pr.min():+10.1f} .. {pr.max():+7.1f} {np.median(d):10.1f} mm {rows[n]:12.3f}")
    near = [n for n in PARTS if rows[n] >= 0.5]
    far = [n for n in PARTS if rows[n] < 0.5]
    print(f"\nnear the mapped TS sternum: {near or 'none'};  far from it: {far or 'none'}")
    if near and far == ["xiphoid process"]:
        print("-> the TotalSegmentator label omits the xiphoid; the union correspondence was biased by construction")
    elif not near:
        print("-> no part is near: not a labelling mismatch; the whole-thorax fit is misplacing the sternum")
    else:
        print("-> mixed: neither explanation is clean; see the per-part distances")

if __name__ == "__main__": main()
