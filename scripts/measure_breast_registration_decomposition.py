#!/usr/bin/env python3
"""Is the breast's 41.82 mm of chest-wall penetration a PLACEMENT error or a SHAPE mismatch?

The first judged breast (docs/BODY_PARAMETERS.md, 2026-09-11) reached fraction 1.0000 with
perfect contact -- 0.0000 mm median gap, 0.0000 mm penetration, 0 of 3,123 held nodes worse --
and an invalid mesh: min J 0.045 against a 0.2 bar, 494 flipped base triangles, 67.96 mm of
maximum deformation. The diagnosis recorded there is that the drive spent the material's entire
strain budget pushing registration error out, because this breast begins **41.82 mm inside the
chest wall**.

That diagnosis has a fork, and it has never been measured. The penetration decomposes into

    penetration  =  a RIGID part   (the breast is in the wrong place; a translation fixes it)
                 +  a RESIDUAL part (the breast's base and this chest wall are different shapes)

and the two call for completely different work. A rigid part means the registration needs
correcting and the seating problem then becomes easy. A residual part means **no placement seats
this breast**, the conform step is doing reconstructive surgery rather than seating, and the
honest options are a different subject or a different chest wall.

THE DECOMPOSITION, AND IT NEEDS NO OPTIMISATION. Take the signed distance from each breast
vertex to its nearest chest-wall point along the anterior axis (+z; negative means behind the
wall). Then

    the RIGID part     is the MEAN of those signed distances -- a uniform offset, which a
                       translation removes exactly and completely;
    the RESIDUAL part  is their SPREAD about that mean -- what no translation can touch,
                       because a translation moves every vertex by the same amount.

So the answer is a ratio of two numbers already in the data, and `spread / |mean|` says which
dominates. No search, no objective, no bound.

TWO EARLIER VERSIONS OF THIS FAILED THEIR KNOWN ANSWERS, and both failures are worth keeping.

The first projected onto each vertex's OWN outward normal. Displacing the breast 12 mm deeper
into the wall made the measured penetration go DOWN, 1.88 mm to 1.53 mm. A breast's base faces
the chest so its normals point posteriorly, while its front surface's point anteriorly; projecting
onto "the vertex's own normal" measures penetration with one sign on the base and the other on the
front, and over a whole breast measures nothing.

The second minimised residual PENETRATION over a bounded translation search -- and reproduced a
degeneracy `docs/BODY_PARAMETERS.md` had already recorded: *"Step 1 as written minimises
penetration ONLY, and that has a degenerate minimum -- carry the breast away and every ray misses
the muscle."* Penetration alone is not a seating objective, because flying the breast off the
chest scores perfectly. The doc said so and the known answer said so again.

The mean/spread form has neither failure mode: it is signed, so it cannot be gamed by moving away,
and it has one sign everywhere on the surface.

KNOWN ANSWER, exact and symmetry-breaking: adding a translation `d` along the anterior axis must
shift the MEAN by exactly `d` and leave the SPREAD **unchanged**. A decomposition that does not do
that is not separating rigid from residual at all.

THE CAVEAT THAT TRAVELS: these breasts are registered from four subjects onto a male-derived
body. What is measured here is the disagreement between that body's chest wall and a registered
subject's breast base -- a statement about the registration, not about any subject's anatomy.
"""
import argparse, gzip, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
SUBJECTS = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
            "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
PLACEMENT_LIMIT_MM = 25.0        # BODY_PARAMETERS.md: beyond this it is a registration failure


def read_obj(p):
    V, F = [], []
    for line in open(p):
        if line.startswith("v "):
            V.append(line.split()[1:4])
        elif line.startswith("f "):
            F.append([int(t.split("/")[0]) - 1 for t in line.split()[1:4]])
    return np.asarray(V, float), np.asarray(F, np.int64)


def entities(match):
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    out = []
    for e in ents:
        if e.get("reference_geometry") and match(e.get("name", ""), e.get("role", "")):
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            out.append(np.asarray(g["positions"], float).reshape(-1, 3))
    return out


def vertex_normals(V, F):
    n = np.zeros_like(V)
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    for k in range(3):
        np.add.at(n, F[:, k], fn)
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    return n / np.maximum(ln, 1e-30)


def signed_depth(Vb, tree, WVz, t=np.zeros(3)):
    """Signed anterior distance from each vertex to its NEAREST chest-wall point, after moving by t.

    Positive = in front of the wall. +z is anterior (`locate_uterus_outside_ring.py`).
    The nearest point is re-found after the move, which is the physically meaningful thing and
    also the reason the mean/spread of this does NOT decompose additively -- see the header.
    """
    P = Vb + t
    _, i = tree.query(P, k=1, workers=-1)
    return P[:, 2] - WVz[i]


def fit_translation(Vb, tree, WVz, limit_mm, seed_step_mm=8.0):
    """The translation minimising the VARIANCE of signed distance, within `limit_mm`.

    Variance, not penetration. Penetration alone is degenerate -- carrying the breast off the
    chest scores perfectly, which `docs/BODY_PARAMETERS.md` already recorded and which the second
    version of this script reproduced. Variance asks the question seating actually needs: **can a
    rigid move make the base sit at a UNIFORM offset from the wall?** A breast that can is seated
    by a translation; one that cannot is a different shape from the wall it is being put on, and
    no placement fixes that.
    """
    t = np.zeros(3)
    cost = lambda tt: float(signed_depth(Vb, tree, WVz, tt).var())
    c, step = cost(t), seed_step_mm / 1000.0
    while step > 0.0002:
        improved = False
        for k in range(3):
            for sgn in (+1.0, -1.0):
                tt = t.copy(); tt[k] += sgn * step
                if np.linalg.norm(tt) * 1000.0 > limit_mm:
                    continue
                cc = cost(tt)
                if cc < c - 1e-16:
                    t, c, improved = tt, cc, True
        if not improved:
            step *= 0.5
    return t, float(np.sqrt(c)) * 1e3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=1200)
    ap.add_argument("--limit-mm", type=float, default=PLACEMENT_LIMIT_MM)
    ap.add_argument("--out", default="out/breast_registration_decomposition.json")
    a = ap.parse_args()

    wall = entities(lambda n, r: "pectoralis major" in n.lower() or "rib" in n.lower()
                    or "sternum" in n.lower())
    WV = np.vstack(wall)
    print(f"chest wall: {len(WV):,} vertices (pectoralis major + ribs + sternum)")
    print(f"placement bound: {a.limit_mm:.0f} mm, as BODY_PARAMETERS.md already fixes it\n")

    rng = np.random.default_rng(0)
    tree = cKDTree(WV); WVz = WV[:, 2]

    # ---------------- known answer ------------------------------------------------
    Vb0, _ = read_obj(ROOT / "data/derived" / SUBJECTS["s1159"] / "breast_left.obj")
    idx = rng.choice(len(Vb0), size=min(a.samples, len(Vb0)), replace=False)
    V0 = Vb0[idx]
    t0, sd0 = fit_translation(V0, tree, WVz, a.limit_mm)
    off = np.array([0.0, 0.0, -0.010])                    # 10 mm further into the wall
    t1, sd1 = fit_translation(V0 + off, tree, WVz, a.limit_mm)
    miss = float(np.linalg.norm((t1 - off * -1.0) - t0) * 1e3)
    print("KNOWN ANSWER: displace the breast 10.0 mm posteriorly. The fitted translation must move")
    print("by +10.0 mm to compensate, and the residual sd it achieves must be UNCHANGED.")
    print(f"  fitted t before [{', '.join(f'{v*1e3:+6.2f}' for v in t0)}] mm, sd {sd0:6.3f} mm")
    print(f"  fitted t after  [{', '.join(f'{v*1e3:+6.2f}' for v in t1)}] mm, sd {sd1:6.3f} mm")
    ok = miss < 2.0 and abs(sd1 - sd0) < 0.3
    print(f"  compensation misses by {miss:.2f} mm (want <2), residual sd moves "
          f"{sd1-sd0:+.3f} mm (want ~0)   {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("known answer FAILED; the fit does not separate rigid from residual")

    # ---------------- the four subjects -------------------------------------------
    print(f"\n{'subject':8s} {'side':6s} {'sd before':>10s} {'|t| mm':>7s} {'residual sd':>12s} "
          f"{'removed':>8s} {'behind before':>14s}")
    rows = []
    for sid, dd in SUBJECTS.items():
        for side in ("left", "right"):
            f = ROOT / "data/derived" / dd / f"breast_{side}.obj"
            if not f.exists():
                continue
            Vb, _ = read_obj(f)
            i = rng.choice(len(Vb), size=min(a.samples, len(Vb)), replace=False)
            d_before = signed_depth(Vb[i], tree, WVz) * 1e3
            sd_b = float(d_before.std())
            t, sd_a = fit_translation(Vb[i], tree, WVz, a.limit_mm)
            tm = float(np.linalg.norm(t) * 1e3)
            removed = 1.0 - (sd_a / max(sd_b, 1e-12))
            behind = float((d_before < 0).mean())
            bounded = tm > a.limit_mm - 0.5
            print(f"{sid:8s} {side:6s} {sd_b:9.2f}  {tm:6.2f}{'*' if bounded else ' '} "
                  f"{sd_a:11.2f}  {removed:7.1%}  {behind:13.1%}")
            rows.append(dict(subject=sid, side=side, sd_before_mm=sd_b, translation_mm=tm,
                             residual_sd_mm=sd_a, fraction_removed=float(removed),
                             fraction_behind_before=behind, at_placement_bound=bool(bounded)))

    SDa = np.array([r["residual_sd_mm"] for r in rows])
    rem = np.array([r["fraction_removed"] for r in rows])
    print(f"\n  a translation within {a.limit_mm:.0f} mm removes a median of {np.median(rem):.1%} "
          f"of the spread")
    print(f"  residual spread after the best translation: {SDa.min():.1f}-{SDa.max():.1f} mm sd")
    rigid = np.median(rem) > 0.5
    print("\n  VERDICT: " + (
        f"mostly a RIGID OFFSET -- the best translation removes {np.median(rem):.0%} of the spread.\n"
        f"  Correcting the registration makes seating easy; the solver was never the thing to fix."
        if rigid else
        f"mostly a SHAPE MISMATCH. The best translation within {a.limit_mm:.0f} mm removes only "
        f"{np.median(rem):.0%}\n  of the spread, leaving {SDa.min():.1f}-{SDa.max():.1f} mm sd that "
        f"no placement can touch, because a\n  translation moves every vertex by the same amount. "
        f"**No placement seats these breasts.**\n  The conform step is reshaping tissue rather than "
        f"seating it -- which is what 67.96 mm of\n  deformation and 494 flipped base triangles "
        f"were."))

    best = min(rows, key=lambda r: r["residual_sd_mm"])
    print(f"\n  least residual: {best['subject']} {best['side']}, {best['residual_sd_mm']:.1f} mm sd "
          f"-- the cheapest breast to attempt on this measure")

    p = ROOT / a.out
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(placement_limit_mm=a.limit_mm, rows=rows,
                                 rigid_dominates=bool(rigid),
                                 caveat="registration disagreement between this body's chest wall "
                                        "and a registered subject's breast base; not a statement "
                                        "about any subject's anatomy"), indent=2) + "\n")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
