#!/usr/bin/env python3
"""Close-ups of the breast and the toes, at a size worth opening.

Two sheets into data/derived/renders/ (gitignored):

  closeup_breast.png   the CT breast tissue of all four female subjects against THIS body's ribs,
                       sternum and pectoralis major. It is UNSEATED: the tissue passes through the
                       chest wall rather than resting on it. The solver meant to seat it has a
                       repaired stepping and a sound Newton solve, and still does not close the gap.
  closeup_toes.png     the forefoot, skin and phalanges, at the same camera.

The body itself has NO breast: a search of all 4,000 entities returns 0 named breast, 0 nipple and
0 mammary gland, and only four surface patches called "<side> mammary region". Everything pink in
the breast sheet is CT tissue from a female subject registered into this body's frame, not anatomy
this body owns.

Rendering reuses scripts/render_body_3d.py's z-buffered rasteriser, as render_female_torso.py does.
Note that rasterise() returns a FLAT pixel index, so assignment is img.reshape(-1, 3)[pix].
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
SKIN, BONE, MUSCLE, BREAST = (0.80, 0.66, 0.58), (0.91, 0.88, 0.80), (0.70, 0.30, 0.28), (0.93, 0.60, 0.62)


def read_obj(p):
    V, F = [], []
    for line in open(p):
        if line.startswith("v "): V.append(line.split()[1:4])
        elif line.startswith("f "): F.append([int(t.split("/")[0]) - 1 for t in line.split()[1:4]])
    return np.asarray(V, float), np.asarray(F, np.int64)


def entities(match):
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    out = []
    for e in ents:
        if e.get("reference_geometry") and match(e.get("name", ""), e.get("role", "")):
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            out.append((np.asarray(g["positions"], float).reshape(-1, 3),
                        np.asarray(g["indices"], np.int64).reshape(-1, 3)))
    return out


def render(meshes, direction, centre, extent, W=900, H=1000, ss=2, fov=24.0):
    d = np.asarray(direction, float); d /= np.linalg.norm(d)
    eye = centre + d * (extent / 2 / np.tan(np.radians(fov / 2)) * 1.08)
    cam = RB.Camera(eye, centre, up=(0, 1, 0), fov_deg=fov, aspect=W / H)
    V, F, C, N, off = [], [], [], [], 0
    for Vm, Fm, rgb in meshes:
        V.append(Vm); F.append(Fm + off); C.append(np.tile(rgb, (len(Vm), 1)))
        N.append(RB.vertex_normals(Vm, Fm)); off += len(Vm)
    V, F, C, N = np.vstack(V), np.vstack(F), np.vstack(C), np.vstack(N)
    Wi, Hi = W * ss, H * ss
    yy = np.linspace(0, 1, Hi)[:, None, None]
    img = (np.array([0.16, 0.17, 0.19]) * (1 - yy) + np.array([0.09, 0.09, 0.10]) * yy) * np.ones((Hi, Wi, 3))
    sx, sy, z, _ = cam.project(V, Wi, Hi)
    r = RB.rasterise(sx, sy, z, F, Wi, Hi, backface=False)
    if r is not None:
        pix, tri, bary, _ = r
        n = (N[F[tri]] * bary[..., None]).sum(1)
        n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
        view = -cam.forward
        n[(n @ view) < 0] *= -1
        key = view * 0.55 + cam.up * 0.55 - cam.right * 0.45; key /= np.linalg.norm(key)
        fill = view * 0.6 + cam.right * 0.6; fill /= np.linalg.norm(fill)
        col = (C[F[tri]] * bary[..., None]).sum(1)
        light = 0.22 + 0.68 * np.clip(n @ key, 0, None) + 0.18 * np.clip(n @ fill, 0, None)
        img.reshape(-1, 3)[pix] = np.clip(col * light[:, None], 0, 1)
    return Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8)).resize((W, H), Image.LANCZOS)


def label(im, title, sub):
    d = ImageDraw.Draw(im)
    try:
        ft = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 26)
        fs = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 15)
    except OSError:
        ft = fs = ImageFont.load_default()
    d.rectangle([0, 0, im.width, 34 + 20 * len(sub)], fill=(14, 15, 18))
    d.text((18, 6), title, fill=(238, 238, 238), font=ft)
    for i, line in enumerate(sub):
        d.text((18, 40 + 19 * i), line, fill=(168, 172, 180), font=fs)
    return im


def strip(ims, path):
    W = sum(i.width for i in ims)
    s = Image.new("RGB", (W, ims[0].height), (20, 21, 24))
    x = 0
    for i in ims:
        s.paste(i, (x, 0)); x += i.width
    s.save(path)
    print(f"wrote {path}  {s.size[0]}x{s.size[1]}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ribs = entities(lambda n, r: "rib" in n.lower() or "sternum" in n.lower())
    pect = entities(lambda n, r: "pectoralis major" in n.lower())
    skin = entities(lambda n, r: r == "skin" or n == "skin")

    # ---- breast: all four subjects' CT tissue against this body's chest wall
    panels = []
    for sid, d in REG.items():
        p = ROOT / "data/derived" / d
        try:
            br = [read_obj(p / f"breast_{s}.obj") for s in ("left", "right")]
        except FileNotFoundError:
            print(f"  {sid}: no breast meshes, skipped"); continue
        scene = [(V, F, BONE) for V, F in ribs] + [(V, F, MUSCLE) for V, F in pect] \
              + [(V, F, BREAST) for V, F in br]
        Vb = np.vstack([b[0] for b in br])
        c = (Vb.min(0) + Vb.max(0)) / 2
        ext = float(np.linalg.norm(Vb.max(0) - Vb.min(0))) * 1.25
        panels.append(label(render(scene, (0.25, 0.10, 0.96), c, ext), f"breast — {sid}",
                            ["pink: CT tissue registered into this body's frame.",
                             "red: this body's pectoralis major. cream: its ribs and sternum.",
                             "the tissue passes THROUGH the chest wall — it is not seated."]))
    if panels:
        strip(panels, OUT / "closeup_breast.png")

    # ---- toes
    toes = entities(lambda n, r: "phalanx" in n.lower() and "toe" in n.lower() and "left" in n.lower())
    Vt = np.vstack([v for v, _ in toes])
    c_toe = (Vt.min(0) + Vt.max(0)) / 2
    ext_toe = float(np.linalg.norm(Vt.max(0) - Vt.min(0))) * 1.8
    dirn = (0.55, 0.35, 0.76)
    skin_m = [(v, f, SKIN) for v, f in skin]
    toe_m = [(v, f, BONE) for v, f in toes]
    strip([label(render(skin_m, dirn, c_toe, ext_toe), "toes — skin",
                 ["the exterior surface over the left forefoot"]),
           label(render(toe_m, dirn, c_toe, ext_toe), "toes — the phalanges alone",
                 ["same camera. gate 2 measures 0.854 / 0.876 of toe-bone vertices",
                  "inside the skin against a 0.95 bar: about 1 in 7 lies outside,",
                  "and gate 4's 73 folded triangles are almost all in this skin"]),
           label(render(skin_m + toe_m, dirn, c_toe, ext_toe), "toes — bone through skin",
                 ["the same pair drawn together"])],
          OUT / "closeup_toes.png")


if __name__ == "__main__":
    main()
