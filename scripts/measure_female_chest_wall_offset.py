"""which way is the female chest wall off, and by how much?

register_female_torso.py stopped at gate d with two failures:
  * this body's sternum partly OUTSIDE the mapped female trunk (manubrium 0.910 inside,
    body 0.787, xiphoid 0.545 -- graded along the bone, worst at the bottom)
  * this body's ribs 2-7 partly INSIDE the mapped breasts (4-22% of each rib's vertices)

if both come from one cause -- this body's anterior chest wall sitting ANTERIOR to the
mapped female chest wall -- then the sternum vertices outside the trunk lie in front of it,
and the rib vertices inside the breasts are the anterior rib ends, sunk to some depth. a
few millimetres is surface noise; centimetres is the whole-thorax similarity failing to
carry the anterior chest wall from one person to another.

the first diagnostic reported a median distance of each rib to the breast surface over ALL
its vertices (58-106 mm), which answers nothing: a rib curves round the back. this
classifies every vertex inside/outside first and measures depth only for the right ones.
"""
import gzip, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "data/derived/female-torso-registered-v1"
ANATOMY = ROOT / "data/derived/canonical/anatomy.json"
ORD = ("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth")

def read_obj(p):
    V, F = [], []
    for l in open(p):
        if l.startswith("v "): V.append([float(x) for x in l.split()[1:4]])
        elif l.startswith("f "): F.append([int(x.split("/")[0]) - 1 for x in l.split()[1:4]])
    return np.asarray(V), np.asarray(F)

def inside(V, F, points, seed=0):
    """per-point ray parity, two opposite rays, both odd = inside: the same test as
    build_skin_contact_meshes.enclosure(), without averaging it away."""
    rng = np.random.default_rng(seed); P = np.asarray(points, float)
    tri = V[F]; a = tri[:, 0]; e1 = tri[:, 1] - a; e2 = tri[:, 2] - a
    d0 = rng.normal(size=3); d0 /= np.linalg.norm(d0); odd = []
    for sgn in (1.0, -1.0):
        d = sgn * d0; pv = np.cross(d, e2); det = (e1 * pv).sum(1)
        par = np.abs(det) < 1e-14; iv = np.where(par, 0.0, 1.0 / np.where(par, 1.0, det))
        hits = np.zeros(len(P), np.int64)
        for s in range(0, len(P), 128):
            blk = P[s:s + 128]; tv = blk[:, None] - a[None]
            u = (tv * pv[None]).sum(2) * iv[None]; qv = np.cross(tv, e1[None])
            v = (qv * d).sum(2) * iv[None]; t = (qv * e2[None]).sum(2) * iv[None]
            hits[s:s + 128] = ((~par[None]) & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 1e-9)).sum(1)
        odd.append(hits % 2 == 1)
    return odd[0] & odd[1]

def main():
    ents = {e["name"]: e for e in json.loads(ANATOMY.read_text())["entities"] if e["role"] == "rigid_bone"}
    def verts(n, k, seed=0):
        g = json.loads(gzip.decompress((ROOT / ents[n]["reference_geometry"]["path"]).read_bytes()))
        P = np.asarray(g["positions"], float).reshape(-1, 3)
        return P[np.random.default_rng(seed).choice(len(P), min(len(P), k), replace=False)]
    # control: points on a known sphere must classify exactly
    ph = np.random.default_rng(1).normal(size=(800, 3)); ph /= np.linalg.norm(ph, axis=1, keepdims=True)
    import trimesh
    ico = trimesh.creation.icosphere(subdivisions=4, radius=0.1)
    ctl_in = inside(ico.vertices, ico.faces, 0.05 * ph).mean(); ctl_out = inside(ico.vertices, ico.faces, 0.5 * ph).mean()
    print(f"CONTROL: points at r=0.05 inside a 0.1 m sphere {ctl_in:.3f} (must be 1.000); at r=0.5 {ctl_out:.3f} (must be 0.000)")
    if not (ctl_in == 1.0 and ctl_out == 0.0): sys.exit("per-point parity is wrong; stopping")

    tV, tF = read_obj(R / "body_trunc_surface.obj"); ttree = cKDTree(tV)
    print("\n1. this body's sternum vertices OUTSIDE the mapped female trunk: how far out, and which way")
    print("   (AP axis +z anterior in this frame; positive = in front of the female chest wall)")
    for part in ("manubrium", "body of sternum", "xiphoid process"):
        P = verts(part, 1500); out = ~inside(tV, tF, P)
        if not out.any(): print(f"   {part:18s} none outside"); continue
        dist, j = ttree.query(P[out]); dz = (P[out] - tV[j])[:, 2]
        print(f"   {part:18s} {out.mean():.3f} outside | by {1000*np.median(dist):5.1f} mm median, {1000*np.percentile(dist, 95):5.1f} mm 95th | "
              f"anterior of the nearest trunk point: {np.mean(dz > 0):.2f} of them")

    print("\n2. this body's rib vertices INSIDE a mapped breast: how deep")
    for side in ("left", "right"):
        bV, bF = read_obj(R / f"breast_{side}.obj"); btree = cKDTree(bV)
        for i, o in enumerate(ORD, 1):
            P = verts(f"{side} {o} rib", 1500); ins = inside(bV, bF, P)
            if ins.sum() < 3: continue
            depth = 1000 * btree.query(P[ins])[0]
            print(f"   {side:5s} rib {i:2d}: {int(ins.sum()):4d} of {len(P)} inside, depth below the breast surface "
                  f"median {np.median(depth):5.1f} mm, max {depth.max():5.1f} mm")

if __name__ == "__main__": main()
