"""Where a registered uterus leaves the pelvic ring -- evidence after a verdict, never a gate.

For each subject named, the part of the uterus surface outside the ring (same bones, sampler and
Delaunay hull as register_pelvic_organs.py) is located relative to the uterus's own centroid:
along this body's vertical (calcanei below the frontal bone) and its anterior (the hip bones'
centroid lies anterior of the sacrum's), both derived from the anatomy and asserted.
Usage: locate_uterus_outside_ring.py RUN/SUBJECT [RUN/SUBJECT ...]
  e.g. ut-endomri-registered-v1/D1-035
"""
import importlib.util, sys
from pathlib import Path
import numpy as np
from scipy.spatial import Delaunay, cKDTree

ROOT = Path(__file__).resolve().parents[1]
def load(name, path):
    s = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
R = load("rpo", ROOT / "scripts/register_pelvic_organs.py")
Q = load("q", ROOT / "scripts/score_uterus_lowest_quartile.py")

def main(targets):
    body = R.body_meshes(); up_axis, up_sign = Q.vertical()
    hips = np.vstack([body["left hip bone"][0], body["right hip bone"][0]]); sacrum = body["sacrum"][0]
    fwd = hips.mean(0) - sacrum.mean(0); fwd[up_axis] = 0
    ant_axis = int(np.argmax(np.abs(fwd))); ant_sign = float(np.sign(fwd[ant_axis]))
    assert ant_axis != up_axis and abs(fwd[ant_axis]) > 0.02, f"hip-to-sacrum offset {fwd} gives no clear anterior"
    H = np.vstack([hips, sacrum]); ring = Delaunay(H); tree = cKDTree(H)
    floor = (up_sign * H[:, up_axis]).min()
    print(f"vertical: axis {up_axis} {up_sign:+.0f}; anterior: axis {ant_axis} {ant_sign:+.0f} (hips {1e3*abs(fwd[ant_axis]):.0f} mm anterior of the sacrum)")
    for t in targets:
        V, F = R.read_obj(ROOT / "data/derived" / t / "uterus.obj")
        P = R.area_samples(V, F, 20000); out = ring.find_simplex(P) < 0
        h = up_sign * P[:, up_axis] - floor
        print(f"{t}: {out.mean():.4f} of the surface outside the ring; uterus {1e3*h.min():.0f}-{1e3*h.max():.0f} mm above the ring's floor")
        if not out.any(): continue
        d = P[out].mean(0) - P.mean(0)
        dist = tree.query(P[out])[0]
        print(f"   outside part: {1e3*up_sign*d[up_axis]:+.0f} mm up and {1e3*ant_sign*d[ant_axis]:+.0f} mm anterior of the uterus centroid, "
              f"median height {1e3*np.median(h[out]):.0f} mm; {1e3*np.median(dist):.1f} mm (median) / {1e3*dist.max():.1f} mm (max) from the nearest ring-bone vertex")

if __name__ == "__main__": main(sys.argv[1:])
