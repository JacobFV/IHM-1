"""how much soft tissue lies between this body's skin and what is under it, region by region.

the skin contact model gives every patch ONE layer thickness (skin_layers: from the three
skin-layer entities). a real body has a few millimetres over the shin and the scalp and
centimetres over the buttocks, and the programme requires contact to be mediated by fat,
muscle and skin, never bone. this body's own geometry already implies the map: for every
point on its exterior skin, the distance to the nearest MUSCLE or BONE surface beneath it --
everything in between (subcutaneous tissue, fat the atlas does not segment) is the layer.

KNOWN ANSWER (anatomy, not set by this script): the anterior shin and the scalp are thin, the
buttocks are thick. the gate: median depth over the anterior shin and over the scalp each
below the median over the buttocks by at least a factor of 2. if not, the measurement is wrong.

nearest-surface distance, not depth along the inward normal: on convex regions the two agree;
in folds (axilla, groin) it underestimates. reported, not hidden.
"""
import gzip, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("bscm", ROOT / "scripts/build_skin_contact_meshes.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
OUT = ROOT / "data/derived/soft-tissue-depth-v1"

def area_samples(V, F, n, seed=0):
    rng = np.random.default_rng(seed); t = V[F]; a = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    if a.sum() <= 0: return np.empty((0, 3))
    k = rng.choice(len(F), n, p=a / a.sum()); r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = t[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2]

def mesh(e):
    g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
    return np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3)

def main():
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    SV, SF = mesh(skin)
    ext = np.asarray(json.loads((ROOT / B.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
    skin_pts = np.unique(SF[ext].ravel())
    rng = np.random.default_rng(0); skin_pts = SV[rng.choice(skin_pts, min(len(skin_pts), 40000), replace=False)]
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    deep = [e for e in ents if e["role"] in ("rigid_bone", "muscle") and e.get("reference_geometry") and e["id"].startswith("body-bp3d-")]
    P, owner = [], []
    for i, e in enumerate(deep):
        V, F = mesh(e); n = 600 if e["role"] == "rigid_bone" else 300
        s = area_samples(V, F, n, seed=i); P.append(s); owner += [i] * len(s)
    P = np.vstack(P); owner = np.asarray(owner)
    print(f"{len(skin_pts)} exterior skin points; {len(deep)} bone and muscle entities, {len(P)} surface samples", flush=True)
    d, j = cKDTree(P).query(skin_pts)
    near = [deep[k]["name"] for k in owner[j]]
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "depth.npz", skin_points_m=skin_pts, depth_m=d, nearest=np.array(near))
    mm = 1000 * d
    print(f"soft-tissue depth over the whole exterior skin: median {np.median(mm):.1f} mm, 10-90% {np.percentile(mm,10):.1f}-{np.percentile(mm,90):.1f} mm")
    def region(pred):
        sel = np.array([pred(n) for n in near]); return mm[sel], int(sel.sum())
    regions = {
        "anterior shin (nearest structure a tibia)": lambda n: n in ("right tibia", "left tibia"),
        "scalp (nearest structure a skull bone)": lambda n: any(k in n for k in ("frontal bone", "parietal bone", "occipital bone")),
        "buttock (nearest structure gluteus maximus)": lambda n: "gluteus maximus" in n,
        "thigh (nearest structure a quadriceps muscle)": lambda n: any(k in n for k in ("vastus", "rectus femoris")),
        "sternum (nearest structure the sternum)": lambda n: any(k in n for k in ("manubrium", "body of sternum", "xiphoid")),
    }
    res = {}
    for name, pred in regions.items():
        v, n = region(pred)
        res[name] = float(np.median(v)) if n else None
        print(f"  {name:48s} " + (f"n {n:5d}  median {np.median(v):5.1f} mm  10-90% {np.percentile(v,10):.1f}-{np.percentile(v,90):.1f}" if n else "no points"))
    shin, scalp, butt = (res[k] for k in list(regions)[:3])
    ok = all(x is not None for x in (shin, scalp, butt)) and shin * 2 <= butt and scalp * 2 <= butt
    print(f"\nGATE (anatomy): shin and scalp each <= half the buttock median -> {'PASS' if ok else 'FAIL'}"
          f"  (shin {shin}, scalp {scalp}, buttock {butt})")
    (OUT / "summary.json").write_text(json.dumps(dict(regions_median_mm=res, gate_pass=ok,
        method="nearest bone/muscle surface distance from exterior skin; underestimates in folds"), indent=2) + "\n")

if __name__ == "__main__": main()
