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

THE DECOMPOSITION. For the base nodes, take the signed distance to the chest wall along each
node's own outward normal. Find the translation `t` minimising the sum of squared *penetrations*
that remain after moving by `t` -- penetration only, since lifting a node clear of the wall costs
nothing and pulling it further in does. Report the penetration before, the penetration after the
best translation, and the fraction removed. **The fraction removed is the answer.**

A BOUND ON THE TRANSLATION, fixed before any number: `docs/BODY_PARAMETERS.md` already records
that a breast needing more than **25 mm** of translation is a REGISTRATION failure for that
subject, not a seating. The same 25 mm bound is used here, and a subject whose best translation
sits on the boundary is reported as bounded rather than solved -- the earlier rigid search ran to
the corner of its box at 41.5 mm magnitude and that is exactly the failure mode to name, not to
absorb.

KNOWN ANSWER, and it breaks the symmetry it tests: displace a breast by a KNOWN translation and
the decomposition must recover it, returning the original penetration and removing essentially
all of the added part. A decomposition that cannot recover a translation it was handed cannot be
trusted to say a translation is absent.

THE CAVEAT THAT TRAVELS: these breasts are registered from four subjects onto a male-derived
body. What is measured here is the disagreement between that body's chest wall and a registered
subject's breast base -- a statement about the registration, not about any subject's anatomy.
"""
import argparse, gzip, json, sys
from pathlib import Path
import numpy as np

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


def penetration(Vb, Nb, WV, t=np.zeros(3), chunk=3000):
    """Depth each base node sits BEHIND the wall, after translating by t. 0 if clear."""
    P = Vb + t
    best_d = np.full(len(P), np.inf); best_i = np.zeros(len(P), int)
    for s in range(0, len(WV), chunk):
        d = np.linalg.norm(P[:, None, :] - WV[None, s:s + chunk, :], axis=2)
        j = d.argmin(1); v = d.min(1)
        m = v < best_d
        best_d[m] = v[m]; best_i[m] = j[m] + s
    # signed along the node's own outward normal: negative means behind the wall
    sgn = np.einsum('ij,ij->i', Nb, P - WV[best_i])
    return np.maximum(-sgn, 0.0)


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

    def best_translation(Vb, Nb):
        """Coordinate search for the translation minimising residual penetration, within the bound."""
        t = np.zeros(3)
        cost = lambda tt: float((penetration(Vb, Nb, WV, tt) ** 2).sum())
        c = cost(t)
        step = a.limit_mm / 1000.0 * 8
        while step > a.limit_mm / 1000.0 * 0.05:
            improved = False
            for k in range(3):
                for s in (+step, -step):
                    tt = t.copy(); tt[k] += s / 1000.0
                    if np.linalg.norm(tt) * 1000.0 > a.limit_mm:
                        continue
                    cc = cost(tt)
                    if cc < c - 1e-12:
                        t, c, improved = tt, cc, True
            if not improved:
                step *= 0.5
        return t, c

    # ---------------- known answer ------------------------------------------------
    sid0 = "s1159"
    Vb0, Fb0 = read_obj(ROOT / "data/derived" / SUBJECTS[sid0] / "breast_left.obj")
    Nb0 = vertex_normals(Vb0, Fb0)
    idx = rng.choice(len(Vb0), size=min(a.samples, len(Vb0)), replace=False)
    Vs, Ns = Vb0[idx], Nb0[idx]
    p0 = penetration(Vs, Ns, WV)
    known = np.array([0.0, 0.0, -0.012])        # 12 mm deeper, straight back
    pk = penetration(Vs + known, Ns, WV)
    tk, _ = best_translation(Vs + known, Ns)
    rec = np.linalg.norm((tk - (-known))) * 1000.0
    print("KNOWN ANSWER: displace the breast 12.0 mm further into the wall; the decomposition")
    print("must recover that translation and return the original penetration.")
    print(f"  penetration median before {np.median(p0)*1e3:6.2f} mm -> displaced "
          f"{np.median(pk)*1e3:6.2f} mm")
    print(f"  recovered translation misses the truth by {rec:.2f} mm   "
          f"{'PASS' if rec < 4.0 else 'FAIL -- the decomposition cannot recover a known offset'}")
    if rec >= 4.0:
        sys.exit("known answer FAILED; no decomposition below is interpretable")

    # ---------------- the four subjects -------------------------------------------
    print(f"\n{'subject':8s} {'side':6s} {'penetration mm':>22s}  {'after best t':>20s}  "
          f"{'removed':>8s}  {'|t| mm':>7s}")
    rows = []
    for sid, d in SUBJECTS.items():
        for side in ("left", "right"):
            f = ROOT / "data/derived" / d / f"breast_{side}.obj"
            if not f.exists():
                continue
            Vb, Fb = read_obj(f)
            Nb = vertex_normals(Vb, Fb)
            i = rng.choice(len(Vb), size=min(a.samples, len(Vb)), replace=False)
            before = penetration(Vb[i], Nb[i], WV)
            t, _ = best_translation(Vb[i], Nb[i])
            after = penetration(Vb[i] + t, Nb[i], WV)
            tm = float(np.linalg.norm(t) * 1000.0)
            removed = 1.0 - after.sum() / max(before.sum(), 1e-30)
            bounded = tm > a.limit_mm - 0.5
            print(f"{sid:8s} {side:6s} med {np.median(before)*1e3:6.2f} max {before.max()*1e3:6.2f}"
                  f"   med {np.median(after)*1e3:6.2f} max {after.max()*1e3:6.2f}"
                  f"   {removed:7.1%}  {tm:6.2f}{'*' if bounded else ' '}")
            rows.append(dict(subject=sid, side=side,
                             penetration_median_mm=float(np.median(before) * 1e3),
                             penetration_max_mm=float(before.max() * 1e3),
                             residual_median_mm=float(np.median(after) * 1e3),
                             residual_max_mm=float(after.max() * 1e3),
                             fraction_removed=float(removed), translation_mm=tm,
                             at_placement_bound=bool(bounded)))
    print(f"  * translation sits on the {a.limit_mm:.0f} mm bound: a REGISTRATION failure for that "
          f"breast, reported as bounded rather than solved")

    frac = np.array([r["fraction_removed"] for r in rows])
    print(f"\n  rigid part removes a median of {np.median(frac):.1%} of the penetration "
          f"(range {frac.min():.1%}-{frac.max():.1%})")
    rigid = np.median(frac) > 0.5
    print("\n  VERDICT: " + (
        f"the penetration is mostly a PLACEMENT error -- a translation within {a.limit_mm:.0f} mm\n"
        f"  removes {np.median(frac):.0%} of it. Correct the registration and the seating problem\n"
        f"  largely goes away; the solver was never the thing to fix." if rigid else
        f"the penetration is mostly a SHAPE MISMATCH -- the best translation within "
        f"{a.limit_mm:.0f} mm\n  removes only {np.median(frac):.0%} of it. **No placement seats "
        f"these breasts.** The conform\n  step is reshaping the tissue rather than seating it, "
        f"which is what the 67.96 mm of\n  deformation and the 494 flipped triangles were. A "
        f"different subject or a different chest\n  wall is the honest next move, not a better "
        f"solver."))

    p = ROOT / a.out
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(placement_limit_mm=a.limit_mm, known_answer_miss_mm=rec,
                                 rows=rows, rigid_dominates=bool(rigid),
                                 caveat="registration disagreement between this body's chest wall "
                                        "and a registered subject's breast base; not a statement "
                                        "about any subject's anatomy"), indent=2) + "\n")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
