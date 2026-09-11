"""Is every paired bone in the scaffold its opposite, mirrored? One of them is not.

The knee cartilage line failed a 3 mm placement gate on a target of 265 faces with 18.5 mm edges,
and the ceiling said the fit was fine (docs/TISSUE_MECHANICS.md). Chasing that found something the
registration cannot see: `l_femur.vtp` is 908 faces, `r_femur.vtp` is 265, and mirrored they differ
by a median 7.58 mm -- while `l_tibia`/`r_tibia` and `l_patella`/`r_patella` agree to 0.00 mm on
every vertex.

The per-segment fit is blind to it. Its residuals read femur_l 3.88 mm and femur_r 3.94: an ICP
against the atlas fits each side about as well whatever surface is there, so a wrong right femur
costs nothing it measures. Mirror agreement is the invariant that catches it, and nothing in this
repository tested it until now.

KNOWN ANSWER, checked before any pair is judged: a mesh mirrored against ITSELF through its own
mid-plane is not the test (a bone is not symmetric); the test is the left file against the right
file. Its zero case is the pairs that ARE exact -- if l_tibia against r_tibia does not read 0.000 mm
the comparison is wrong and no other row means anything.

Mirrored on the model's own left-right axis, taken from the pair's own centroid separation rather
than assumed, with the two centroids aligned before comparing: this asks whether the SURFACES are
reflections, not whether they sit in mirrored positions.
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.spatial.vtk import surface as read_surface

MODEL = ROOT / "data/models/engineering_stance_v1/model.osim"
GEOMETRY = ROOT / "data/raw/anatomy/opensim-models/source/Geometry"

def pairs():
    """Every mesh the model references that has an opposite-side counterpart."""
    root = ET.parse(MODEL).getroot().find("Model")
    files = {}
    for body in root.iter("Body"):
        for mesh in body.iter("Mesh"):
            files.setdefault(body.get("name"), []).append(mesh.findtext("mesh_file"))
    out = []
    for body, names in sorted(files.items()):
        if not body.endswith("_r"): continue
        left = body[:-2] + "_l"
        if left not in files: continue
        for right_file, left_file in zip(names, files[left]):
            out.append((body, left, right_file, left_file))
    return out

def mirrored_distance(right, left):
    """Median and max distance from each right vertex to the mirrored left surface."""
    axis = int(np.argmax(np.abs(left.mean(0) - right.mean(0))))
    flipped = left.copy()
    flipped[:, axis] *= -1
    flipped = flipped - flipped.mean(0) + right.mean(0)
    d = cKDTree(flipped).query(right)[0]
    back = cKDTree(right).query(flipped)[0]
    both = np.concatenate([d, back])
    return axis, float(np.median(both)), float(both.max())

def main():
    rows = []
    for body, left_body, right_file, left_file in pairs():
        R, FR = read_surface(GEOMETRY / right_file)
        L, FL = read_surface(GEOMETRY / left_file)
        axis, median, worst = mirrored_distance(R, L)
        rows.append(dict(body=body, right=right_file, left=left_file, right_faces=len(FR),
                         left_faces=len(FL), axis=axis, median_mm=1e3 * median, max_mm=1e3 * worst))
    exact = [r for r in rows if r["max_mm"] < 1e-6]
    assert exact, "no pair is an exact mirror: the comparison itself is wrong, so nothing below holds"
    print(f"known answer: {len(exact)} of {len(rows)} pairs are exact mirrors to 0.000 mm "
          f"({', '.join(sorted(r['body'][:-2] for r in exact))})\n")
    print(f"{'body':10s} {'right mesh':18s} {'faces R/L':>12s} {'median':>9s} {'max':>9s}")
    for r in sorted(rows, key=lambda r: -r["median_mm"]):
        flag = "  <- not a mirror" if r["median_mm"] > 0.5 else ""
        print(f"{r['body']:10s} {r['right']:18s} {r['right_faces']:5d}/{r['left_faces']:<6d} "
              f"{r['median_mm']:7.3f} mm {r['max_mm']:7.3f} mm{flag}")
    odd = [r for r in rows if r["median_mm"] > 0.5]
    print(f"\n{len(rows) - len(odd)} of {len(rows)} meshes are their opposite mirrored; {len(odd)} are not.")
    for r in odd:
        print(f"  {r['body']}: {r['right']} ({r['right_faces']} faces) against {r['left']} "
              f"({r['left_faces']}) differs by a median {r['median_mm']:.2f} mm, max {r['max_mm']:.2f} mm")
    if odd:
        print("\nA pair that is not a mirror is not a left-right difference this body measured; it is\n"
              "two different surfaces standing in for one bone, and every bilateral measurement made\n"
              "on them inherits it. The per-segment registration does not see it: femur_l 3.88 mm\n"
              "against femur_r 3.94 mm.")

if __name__ == "__main__": main()
