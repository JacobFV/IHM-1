#!/usr/bin/env python3
"""Render the female torsos: TotalSegmentator CT subjects registered into this body's frame.

Three figures, written to data/derived/renders/ (not in git):
  female_torso_s0790.png      s0790's trunk surface -- front, three-quarter, side
  female_torso_subjects.png   the four registered subjects' trunk surfaces, front
  female_breast_seating.png   s0790's breast tissue over THIS body's rib cage, sternum and
                              pectoralis major: the seating problem (float and rib overlap)

Source: anonymised CT, TotalSegmentator (CC BY 4.0). The trunk is what each scan's field of view
kept, so it is cut at the neck and hips. The meshes are in the canonical atlas frame, metres:
+y up, +z anterior (scripts/locate_uterus_outside_ring.py derives and asserts both).
Rendering reuses scripts/render_body_3d.py's z-buffered rasteriser.
"""
import gzip, importlib.util, json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rb3", ROOT / "scripts/render_body_3d.py")
RB = importlib.util.module_from_spec(spec); spec.loader.exec_module(RB)
OUT = ROOT / "data/derived/renders"
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
SKIN, BREAST, BONE, MUSCLE = (0.80, 0.66, 0.58), (0.93, 0.60, 0.62), (0.91, 0.88, 0.80), (0.70, 0.30, 0.28)

def read_obj(p):
    V, F = [], []
    for line in open(p):
        if line.startswith("v "): V.append(line.split()[1:4])
        elif line.startswith("f "): F.append([int(t.split("/")[0]) - 1 for t in line.split()[1:4]])
    return np.asarray(V, float), np.asarray(F, np.int64)

def canonical(names):
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    out = []
    for e in ents:
        if e["id"].startswith("body-bp3d-") and e.get("reference_geometry") and names(e["name"], e["role"]):
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            out.append((np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3)))
    return out

def render(meshes, direction, W=640, H=900, ss=2, fov=22.0, centre=None, extent=None):
    """meshes: [(V, F, rgb)]. Two-sided Lambert with a key and a fill light, supersampled ss x ss."""
    Vs = np.vstack([m[0] for m in meshes])
    c = (Vs.min(0) + Vs.max(0)) / 2 if centre is None else centre
    ext = np.linalg.norm(Vs.max(0) - Vs.min(0)) if extent is None else extent
    d = np.asarray(direction, float); d /= np.linalg.norm(d)
    eye = c + d * (ext / 2 / np.tan(np.radians(fov / 2)) * 1.08)
    cam = RB.Camera(eye, c, up=(0, 1, 0), fov_deg=fov, aspect=W / H)
    V, F, C, N, off = [], [], [], [], 0
    for Vm, Fm, rgb in meshes:
        V.append(Vm); F.append(Fm + off); C.append(np.tile(rgb, (len(Vm), 1))); N.append(RB.vertex_normals(Vm, Fm)); off += len(Vm)
    V, F, C, N = np.vstack(V), np.vstack(F), np.vstack(C), np.vstack(N)
    Wi, Hi = W * ss, H * ss
    yy = np.linspace(0, 1, Hi)[:, None, None]
    img = (np.array([0.16, 0.17, 0.19]) * (1 - yy) + np.array([0.09, 0.09, 0.10]) * yy) * np.ones((Hi, Wi, 3))
    sx, sy, z, _ = cam.project(V, Wi, Hi)
    r = RB.rasterise(sx, sy, z, F, Wi, Hi, backface=False)
    if r is not None:
        pix, tri, bary, _ = r
        n = (N[F[tri]] * bary[..., None]).sum(1); n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
        view = -cam.forward
        n[(n @ view) < 0] *= -1                                   # two-sided: open CT cuts show their inside
        key = view * 0.55 + cam.up * 0.55 - cam.right * 0.45; key /= np.linalg.norm(key)
        fill = view * 0.6 + cam.right * 0.6; fill /= np.linalg.norm(fill)
        col = (C[F[tri]] * bary[..., None]).sum(1)
        half = key + view; half /= np.linalg.norm(half)
        light = 0.22 + 0.68 * np.clip(n @ key, 0, None) + 0.18 * np.clip(n @ fill, 0, None)
        spec_ = 0.10 * np.clip(n @ half, 0, None) ** 24
        img.reshape(-1, 3)[pix] = np.clip(col * light[:, None] + spec_[:, None], 0, 1)
    img = img.reshape(H, ss, W, ss, 3).mean((1, 3))
    return Image.fromarray((img * 255).astype(np.uint8))

def label(im, text, sub=None):
    d = ImageDraw.Draw(im)
    try: f = ImageFont.truetype("DejaVuSans.ttf", 22); fs = ImageFont.truetype("DejaVuSans.ttf", 15)
    except OSError: f = fs = ImageFont.load_default()
    d.text((18, 14), text, fill=(235, 235, 235), font=f)
    if sub: d.text((18, 44), sub, fill=(170, 170, 170), font=fs)
    return im

def strip(ims):
    W = sum(i.width for i in ims); H = max(i.height for i in ims)
    out = Image.new("RGB", (W, H)); x = 0
    for i in ims: out.paste(i, (x, 0)); x += i.width
    return out

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Display only: the trunk is a marching-cubes surface of a CT label, terraced at the voxel
    # size. Taubin smoothing (volume-preserving, no shrinkage) removes the terraces for viewing;
    # the meshes on disk and everything measured from them are untouched.
    import trimesh
    def smooth(VF):
        m = trimesh.Trimesh(*VF, process=False); trimesh.smoothing.filter_taubin(m, iterations=30)
        return np.asarray(m.vertices), np.asarray(m.faces)
    trunk = {s: smooth(read_obj(ROOT / "data/derived" / d / "body_trunc_surface.obj")) for s, d in REG.items()}
    s0790 = [(*trunk["s0790"], SKIN)]
    views = [((0, 0, 1), "front"), ((0.72, 0, 0.69), "three-quarter"), ((1, 0, 0), "side")]
    ims = [label(render(s0790, d), f"s0790 - {name}", "CT trunk surface, registered to this body (display-smoothed)") for d, name in views]
    strip(ims).save(OUT / "female_torso_s0790.png")
    ims = [label(render([(*trunk[s], SKIN)], (0, 0, 1)), s, "front") for s in REG]
    strip(ims).save(OUT / "female_torso_subjects.png")
    # Rigid bones only: a name match on "costal" also pulls in the intercostal muscle sheets and
    # vessels, which cover the chest and hide what this figure is for.
    ribs = canonical(lambda n, r: r == "rigid_bone" and (n.endswith(" rib") or any(k in n for k in ("sternum", "manubrium", "xiphoid"))))
    pect = canonical(lambda n, r: r == "muscle" and "pectoralis major" in n)
    assert len(ribs) == 27 and len(pect) == 6, (len(ribs), len(pect))
    breasts = [read_obj(ROOT / "data/derived" / REG["s0790"] / f"breast_{s}.obj") for s in ("left", "right")]
    scene = [(V, F, BONE) for V, F in ribs] + [(V, F, MUSCLE) for V, F in pect] + [(V, F, BREAST) for V, F in breasts]
    c = np.vstack([b[0] for b in breasts] + [p[0] for p in pect]).mean(0)
    ims = [label(render(scene, d, W=760, H=760, centre=c, extent=0.46), f"seating - {name}",
                 "s0790 breast tissue (pink) on this body's pectoralis (red) and ribs")
           for d, name in [((0, 0, 1), "front"), ((1, 0, 0.35), "side")]]
    strip(ims).save(OUT / "female_breast_seating.png")
    for p in ("female_torso_s0790.png", "female_torso_subjects.png", "female_breast_seating.png"): print(OUT / p)

if __name__ == "__main__": main()
