#!/usr/bin/env python3
"""Two observations made by eye, turned into measurements.

1. "The toes look glued together — not individual protrusions."
   The decisive test is binary and needs no threshold: **in a real foot the space between two toes
   is OUTSIDE the body.** Air. So take the midpoint between adjacent toes' proximal-phalanx
   centroids and ask whether it is inside the skin. Outside -> there is a cleft. Inside -> the
   surface bridges the gap and the forefoot is a mitten.

2. "s1067 looks the most realistic of the four breasts."
   The render's pink sits in front of the pectoralis for s1067 and behind it for the others, which
   is a guess about a measurable thing: how much of each subject's breast tissue lies POSTERIOR to
   this body's chest wall, i.e. penetrating it. If the eye is tracking that, s1067 penetrates least.
   If it is not, the camera is flattering one subject and the impression is about the view.

THE FIRST VERSION USED RAY-PARITY CONTAINMENT AND ITS KNOWN ANSWER FAILED -- correctly. The skin's
own centroid read OUTSIDE. The detector was not buggy: the canonical skin has **1,512 boundary
edges**, so it is an OPEN surface, and a ray through a hole flips the parity. A containment test
does not apply to this mesh at all. Recorded because the failure is the useful part: the precondition
was never checked before the test was written.

Both questions are therefore answered by DISTANCE, which needs no containment:

  toes    -- the minimum distance from the straight segment between two adjacent toe-bone centroids
             to the skin surface, against the same distance from a point in the MIDDLE of a single
             toe. A real cleft brings skin close to the inter-toe line; a mitten leaves it as deep in
             flesh as the middle of a toe is. The intra-toe figure is the control, so the comparison
             is relative and no threshold is chosen.
  breast  -- the atlas frame is +y up, +z anterior (asserted in locate_uterus_outside_ring.py). For
             each breast vertex take the nearest chest-wall vertex; the breast vertex is BEHIND the
             wall when its z is smaller. That fraction is what the eye reads as "buried".
"""
import gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}


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
            out.append((e.get("name", ""),
                        np.asarray(g["positions"], float).reshape(-1, 3),
                        np.asarray(g["indices"], np.int64).reshape(-1, 3)))
    return out


def nearest(points, V, chunk=4000):
    """Distance from each point to the nearest vertex of V, in metres."""
    P = np.asarray(points, float)
    best = np.full(len(P), np.inf)
    for s in range(0, len(V), chunk):
        d = np.linalg.norm(P[:, None, :] - V[None, s:s + chunk, :], axis=2)
        best = np.minimum(best, d.min(1))
    return best


def boundary_edges(F):
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]], 0)
    e.sort(1)
    _, c = np.unique(e, axis=0, return_counts=True)
    return int((c == 1).sum())


def main():
    skin = entities(lambda n, r: r == "skin" or n == "skin")
    if not skin:
        sys.exit("no skin entity")
    _, SV, SF = skin[0]
    nb = boundary_edges(SF)
    print(f"skin: {len(SV):,} vertices, {len(SF):,} triangles, {nb:,} boundary edges "
          f"({'OPEN -- containment by parity does not apply' if nb else 'closed'})")

    order = ["great", "second", "third", "fourth", "fifth"]
    for side in ("left", "right"):
        cents = {}
        for name, V, F in entities(lambda n, r: "proximal phalanx" in n.lower()
                                   and "toe" in n.lower() and side in n.lower()):
            for k in order:
                if k in name.lower():
                    cents[k] = V.mean(0)
        have = [k for k in order if k in cents]
        if len(have) < 2:
            print(f"\n{side} forefoot: only {len(have)} proximal toe phalanges found, skipped")
            continue
        # control: how far is the skin from a point in the MIDDLE of a single toe?
        intra = nearest(np.array([cents[k] for k in have]), SV) * 1000
        print(f"\n{side} forefoot -- distance from the skin surface, in mm")
        print(f"  CONTROL, inside a single toe: median {np.median(intra):5.2f}, "
              f"range {intra.min():.2f}-{intra.max():.2f}")
        inter = []
        for x, y in zip(have, have[1:]):
            seg = cents[x] + (cents[y] - cents[x]) * np.linspace(0.25, 0.75, 9)[:, None]
            d = nearest(seg, SV).min() * 1000
            gap = float(np.linalg.norm(cents[x] - cents[y])) * 1000
            inter.append(d)
            print(f"  between {x:6s} and {y:6s} (centres {gap:5.1f} mm apart): {d:5.2f}")
        inter = np.array(inter)
        ratio = float(np.median(inter) / np.median(intra))
        print(f"  inter-toe median {np.median(inter):.2f} against the intra-toe control "
              f"{np.median(intra):.2f}  ->  ratio {ratio:.2f}")
        print("  a real cleft brings the skin CLOSER between the toes than inside one (ratio < 1).")
        print("  a mitten leaves the gap as deep in flesh as a toe's middle (ratio near or above 1).")

    # -------------------------------------------------------------- breast
    wall = entities(lambda n, r: "pectoralis major" in n.lower() or "rib" in n.lower()
                    or "sternum" in n.lower())
    WV = np.vstack([v for _, v, _ in wall])
    print(f"\nchest wall: {len(WV):,} vertices (pectoralis major + ribs + sternum)")
    print("breast vertices BEHIND the nearest chest-wall point (+z is anterior) -- the 'buried' fraction:")
    rows = []
    rng = np.random.default_rng(0)
    for sid, d in REG.items():
        p = ROOT / "data/derived" / d
        try:
            br = [read_obj(p / f"breast_{s}.obj") for s in ("left", "right")]
        except FileNotFoundError:
            continue
        Vb = np.vstack([v for v, _ in br])
        idx = rng.choice(len(Vb), size=min(1500, len(Vb)), replace=False)
        P = Vb[idx]
        # nearest wall vertex per sampled breast vertex, then compare anterior coordinate
        best_i = np.zeros(len(P), int)
        best_d = np.full(len(P), np.inf)
        for s0 in range(0, len(WV), 4000):
            dd = np.linalg.norm(P[:, None, :] - WV[None, s0:s0 + 4000, :], axis=2)
            j = dd.argmin(1); v = dd.min(1)
            upd = v < best_d
            best_d[upd] = v[upd]; best_i[upd] = j[upd] + s0
        behind = float((P[:, 2] < WV[best_i][:, 2]).mean())
        rows.append((behind, sid, len(Vb)))
    for behind, sid, n in sorted(rows):
        print(f"  {sid}  {behind:6.1%} behind the wall   ({n:,} vertices total)")
    if rows:
        print(f"\n  least buried: {sorted(rows)[0][1]}")


if __name__ == "__main__":
    main()
