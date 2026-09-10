#!/usr/bin/env python3
"""Quick overview renders of this body: skin, skeleton, and the skin over the skeleton.

Work in progress, for looking at. Everything is the canonical anatomy in its own frame
(+y up, +z anterior); no plant, no pose, no contact. Written to data/derived/renders/.
"""
import gzip, importlib.util, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rft", ROOT / "scripts/render_female_torso.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
OUT = ROOT / "data/derived/renders"
SKIN, BONE, MUSCLE = (0.80, 0.66, 0.58), (0.91, 0.88, 0.80), (0.70, 0.30, 0.28)

def meshes(pred, limit=None):
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    out = []
    for e in ents:
        if not (e["id"].startswith("body-bp3d-") and e.get("reference_geometry")): continue
        if not pred(e["name"], e["role"]): continue
        g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
        out.append((np.asarray(g["positions"], float).reshape(-1, 3),
                    np.asarray(g["indices"], np.int64).reshape(-1, 3)))
        if limit and len(out) >= limit: break
    return out

def skin_mesh():
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    e = next(x for x in mech["entities"] if x["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / "data/research/engineered_skin_territories/materialization.json").read_text())
                     ["contact_eligible_triangle_ids"], np.int64)
    return V, F[ext]

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    skin = skin_mesh()
    bones = meshes(lambda n, r: r == "rigid_bone")
    muscles = meshes(lambda n, r: r == "muscle")
    print(f"skin {len(skin[1]):,} triangles; {len(bones)} bones; {len(muscles)} muscles", flush=True)
    centre = (skin[0].min(0) + skin[0].max(0)) / 2
    extent = float(np.linalg.norm(skin[0].max(0) - skin[0].min(0)))
    def shot(scene, direction, title, sub):
        im = R.render(scene, direction, W=520, H=1000, ss=2, fov=20, centre=centre, extent=extent)
        return R.label(im, title, sub)
    R.strip([shot([(*skin, SKIN)], (0, 0, 1), "skin", "the exterior surface, 109k triangles"),
             shot([(*skin, SKIN)], (0.72, 0, 0.69), "skin", "three-quarter"),
             shot([(V, F, BONE) for V, F in bones], (0, 0, 1), "skeleton", f"{len(bones)} bones"),
             shot([(V, F, MUSCLE) for V, F in muscles], (0, 0, 1), "muscles", f"{len(muscles)} entities"),
             ]).save(OUT / "body_overview.png")
    print(OUT / "body_overview.png")

if __name__ == "__main__": main()
