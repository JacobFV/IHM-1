#!/usr/bin/env python3
"""Close-ups of the places today's gates failed: the toes, and the chest where a breast is not.

Three panels, written to data/derived/renders/closeups.png (not in git):

  toes, skin        the exterior surface over the forefoot
  toes, skin + bone the same view with the 50 toe phalanges drawn THROUGH the skin, which is
                    what gate 2 measures: 0.854/0.876 of toe-bone vertices inside the skin
                    against a 0.95 bar, so roughly one vertex in seven is outside it
  chest             the mammary region of THIS body. It carries no breast: a search of all
                    4,000 entities returns 0 for "breast", 0 for "nipple" and 0 for a mammary
                    GLAND -- only four surface patches named "<side> mammary region". The only
                    breast geometry in the repository is CT tissue from four female subjects,
                    which has never been seated (see data/derived/renders/female_breast_seating.png).

Rendering reuses scripts/render_body_3d.py's z-buffered rasteriser, as render_female_torso.py does.
"""
import gzip, importlib.util, json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rb3", ROOT / "scripts/render_body_3d.py")
RB = importlib.util.module_from_spec(spec); spec.loader.exec_module(RB)
OUT = ROOT / "data/derived/renders"
SKIN, BONE, REGION = (0.80, 0.66, 0.58), (0.91, 0.88, 0.80), (0.93, 0.60, 0.62)


def entities(match):
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    out = []
    for e in ents:
        if e.get("reference_geometry") and match(e.get("name", ""), e.get("role", "")):
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            out.append((np.asarray(g["positions"], float).reshape(-1, 3),
                        np.asarray(g["indices"], np.int64).reshape(-1, 3)))
    return out


def render(meshes, direction, centre, extent, W=760, H=900, ss=2, fov=24.0):
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
        # rasterise returns a FLAT pixel index (flat[win]), not (x, y) pairs
        img.reshape(-1, 3)[pix] = np.clip(col * light[:, None], 0, 1)
    im = Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8)).resize((W, H), Image.LANCZOS)
    return im


def label(im, title, sub):
    d = ImageDraw.Draw(im)
    try:
        ft = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
        fs = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
    except OSError:
        ft = fs = ImageFont.load_default()
    d.text((16, 14), title, fill=(235, 235, 235), font=ft)
    for i, line in enumerate(sub):
        d.text((16, 44 + 18 * i), line, fill=(165, 168, 175), font=fs)
    return im


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    skin = entities(lambda n, r: r == "skin" or n == "skin")
    if not skin:
        skin = entities(lambda n, r: "skin" in n.lower() and "layer" not in n.lower())
    toes = entities(lambda n, r: "phalanx" in n.lower() and "toe" in n.lower() and "left" in n.lower())
    mamm = entities(lambda n, r: "mammary region" in n.lower())
    print(f"skin surfaces {len(skin)}, left toe phalanges {len(toes)}, mammary regions {len(mamm)}")

    Vt = np.vstack([v for v, _ in toes])
    c_toe = (Vt.min(0) + Vt.max(0)) / 2
    ext_toe = float(np.linalg.norm(Vt.max(0) - Vt.min(0))) * 1.9
    Vm = np.vstack([v for v, _ in mamm])
    c_ch = (Vm.min(0) + Vm.max(0)) / 2
    ext_ch = float(np.linalg.norm(Vm.max(0) - Vm.min(0))) * 1.15

    skin_m = [(v, f, SKIN) for v, f in skin]
    toe_m = [(v, f, BONE) for v, f in toes]
    mam_m = [(v, f, REGION) for v, f in mamm]
    dirn = (0.55, 0.35, 0.76)

    panels = [
        label(render(skin_m, dirn, c_toe, ext_toe), "toes - skin",
              ["the exterior surface over the left forefoot"]),
        label(render(toe_m, dirn, c_toe, ext_toe), "toes - the phalanges alone",
              ["same camera. gate 2 measures 0.854 / 0.876 of toe-bone vertices",
               "inside the skin against a 0.95 bar: about 1 in 7 lies outside it,",
               "and gate 4's 73 folded triangles are almost all in this skin"]),
        label(render(skin_m + mam_m, (0.0, 0.15, 1.0), c_ch, ext_ch), "chest - where a breast is not",
              ["pink: the four 'mammary region' surface patches, all this body has.",
               "0 entities named breast, 0 nipple, 0 mammary gland, of 4,000.",
               "the only breast geometry is unseated CT tissue from 4 subjects."]),
    ]
    W = sum(p.width for p in panels)
    sheet = Image.new("RGB", (W, panels[0].height), (20, 21, 24))
    x = 0
    for p in panels:
        sheet.paste(p, (x, 0)); x += p.width
    path = OUT / "closeups.png"
    sheet.save(path)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
