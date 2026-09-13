#!/usr/bin/env python3
"""Is the ~19 mm the ANTERIOR CHEST -- the place no bone centroid looks?

Pre-registered in `docs/BODY_PARAMETERS.md`, 2026-09-13, and committed before this ran.

`register_female_torso.py` fits bone CENTROIDS. A rib's centroid is dominated by its long
lateral and posterior arc, not by its anterior end, so a fit can match all 34-39 centroids to
11-16 mm and still leave the two anterior chest walls far apart -- almost none of the fitted
quantity lives there. `register.log`'s gate d1 is consistent with it: ribs contained at
0.97-0.999 while the sternum, the one entirely anterior bone, is contained at 0.09-0.767 with
95-98% of the escaping vertices anterior.

THE MEASUREMENT. Map each subject's bones by that registration's own transform. Over a grid of
(lateral, superoinferior) cells covering the breast footprint, compare the most-anterior bone
surface of this body against the mapped subject's -- and, in the SAME cells over the SAME bone
union, the most-posterior one.

THE GATE, two halves:
  (i)  magnitude    -- the anterior difference is >= 50% of that subject's measured anterior
                       correction in >= 3 of 4 subjects.
  (ii) localisation -- the anterior difference is >= 2x the posterior one in >= 3 of 4.
  CONFIRMED only if both hold; REFUTED otherwise, with the fractions stated.

Half (ii) is the half that can fail, and that is the point. A global AP offset in the fit would
satisfy (i) and say nothing about the anterior chest -- it moves front and back equally, and the
anisotropy entry already showed the fit is displaced somewhat. Only (ii) separates a chest wall
that is a different SHAPE from a fit that is bodily DISPLACED.

KNOWN ANSWER 1, symmetry-breaking on the axis under test. Displace this body's bones +10.0 mm
anteriorly. BOTH differences must increase by 10.0 +/- 0.5, because a translation moves
everything equally. A measurement where only the anterior one moves is not a signed distance.

KNOWN ANSWER 2, the null. This body's bones against themselves through the identity map, with
the identical sampling, must return 0.0 +/- 0.2 mm in both regions. Anything else is the
sampler's own bias and every number below it is that bias.

KNOWN ANSWER 3, the axis convention, reused from `measure_registration_anisotropy.py` rather
than re-assumed: T1 -> T12 is superoinferior, left -> right clavicle is lateral, the remaining
axis is anterior. Five AP/ML selection failures on this line say not to assume it.

CAVEAT: four CT subjects registered onto a male-derived body. A statement about the
registration and this body's thorax, not about any subject's anatomy.
"""
import argparse, importlib.util, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("mra", ROOT / "scripts/measure_registration_anisotropy.py")
MRA = importlib.util.module_from_spec(spec); spec.loader.exec_module(MRA)

SRC_ROOT = ROOT / "data/derived/female-torso-totalsegmentator-v1"
REG = MRA.REG
MEASURED = MRA.MEASURED_ANTERIOR_MM
ORD = MRA.ORD
RIBS = range(2, 8)                 # the ribs a breast sits on
VERTS = range(2, 8)
CELL_M = 0.005                     # 5 mm grid over the footprint
RADIUS_M = 0.008                   # bone vertices within 8 mm of a cell centre
MIN_PER_CELL = 4                   # both sides must have this many, or the cell is dropped
FRACTION_MIN, LOCALISATION_MIN, N_MIN = 0.50, 2.0, 3


def bone_sets():
    """(source label -> this body's entity names) for the ribs, sternum and vertebrae a breast
    sits on. Same naming the registration uses, so the two sides are the same bones."""
    pairs = {"sternum": ["manubrium", "body of sternum", "xiphoid process"]}
    for side in ("left", "right"):
        for i in RIBS:
            pairs[f"rib_{side}_{i}"] = [f"{side} {ORD[i-1]} rib"]
    for i in VERTS:
        pairs[f"vertebrae_T{i}"] = [f"{ORD[i-1]} thoracic vertebra"]
    return pairs


def body_vertices(pairs, by_name):
    out = []
    for names in pairs.values():
        for nm in names:
            cand = by_name.get(nm, [])
            if len(cand) != 1:
                raise SystemExit(f"expected one entity named {nm!r}, found {len(cand)}")
            out.append(MRA.entity_mesh(cand[0])[0])
    return np.vstack(out)


def source_vertices(sid, pairs):
    out = []
    for lab in pairs:
        p = SRC_ROOT / sid / "meshes" / f"{lab}.obj"
        if not p.exists():                     # a label the scan cut; the fit dropped it too
            continue
        out.append(MRA.read_obj(p)[0])
    return np.vstack(out)


def envelope(V, cells, ax_ml, ax_si, ax_ap):
    """Per cell, the most-anterior and most-posterior bone coordinate, or NaN if too few.

    A cell is a column through the thorax: over the SAME bone union, its maximum along the
    anterior axis is the front of the chest and its minimum is the back. Taking both from one
    population is what makes the localisation half of the gate a like-for-like comparison.
    """
    key = np.stack([np.round(V[:, ax_ml] / CELL_M), np.round(V[:, ax_si] / CELL_M)], 1)
    ant = np.full(len(cells), np.nan)
    post = np.full(len(cells), np.nan)
    for i, (cx, cy) in enumerate(cells):
        m = (np.abs(V[:, ax_ml] - cx * CELL_M) <= RADIUS_M) & (np.abs(V[:, ax_si] - cy * CELL_M) <= RADIUS_M)
        if m.sum() < MIN_PER_CELL:
            continue
        z = V[m, ax_ap]
        ant[i], post[i] = z.max(), z.min()
    return ant, post


def footprint_cells(sid, ax_ml, ax_si):
    """The (lateral, superoinferior) cells the two registered breasts cover."""
    cs = set()
    for side in ("left", "right"):
        p = ROOT / "data/derived" / REG[sid] / f"breast_{side}.obj"
        V = MRA.read_obj(p)[0]
        for cx, cy in np.unique(np.stack([np.round(V[:, ax_ml] / CELL_M),
                                          np.round(V[:, ax_si] / CELL_M)], 1), axis=0):
            cs.add((int(cx), int(cy)))
    return sorted(cs)


def compare(Vb, Vs, cells, ax):
    ax_ml, ax_si, ax_ap = ax
    a_b, p_b = envelope(Vb, cells, ax_ml, ax_si, ax_ap)
    a_s, p_s = envelope(Vs, cells, ax_ml, ax_si, ax_ap)
    ok = ~(np.isnan(a_b) | np.isnan(a_s))
    if ok.sum() < 20:
        raise SystemExit(f"only {ok.sum()} usable cells; the footprint sampling is not comparable")
    return (float(np.median(a_b[ok] - a_s[ok]) * 1e3),
            float(np.median(p_b[ok] - p_s[ok]) * 1e3), int(ok.sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/anterior_chest_offset.json")
    a = ap.parse_args()

    print("Is the ~19 mm the ANTERIOR CHEST -- the place no bone centroid looks?")
    print(f"pre-registered gate: anterior difference >= {FRACTION_MIN:.0%} of the measured "
          f"correction in >= {N_MIN} of 4, AND >= {LOCALISATION_MIN:.0f}x the posterior "
          f"difference in >= {N_MIN} of 4\n")

    pairs = bone_sets()
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    by_name = {}
    for e in ents:
        if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-"):
            by_name.setdefault(e["name"], []).append(e)

    # ---- known answer 3: the axis convention, from the registration's own correspondences
    src0, dst0, labels0, named0 = MRA.build_correspondences("s0790")
    ax_ml, ax_si, ax_ap, si, ml = MRA.identify_axes(named0)
    print(f"KNOWN ANSWER 3  AXES: SI = {ax_si} (T1->T12 {np.round(si*1e3,1)} mm), "
          f"ML = {ax_ml} (L->R clavicle {np.round(ml*1e3,1)} mm), AP = {ax_ap}   PASS")
    ax = (ax_ml, ax_si, ax_ap)

    Vb = body_vertices(pairs, by_name)
    cells = footprint_cells("s0790", ax_ml, ax_si)
    print(f"\nbone union: {len(Vb):,} vertices over {len(pairs)} labels "
          f"(sternum, ribs 2-7, T2-T7); footprint {len(cells)} cells of {CELL_M*1e3:.0f} mm")

    # ---- known answer 2: the null
    n_a, n_p, n_cells = compare(Vb, Vb.copy(), cells, ax)
    ok2 = abs(n_a) < 0.2 and abs(n_p) < 0.2
    print(f"\nKNOWN ANSWER 2  NULL (this body against itself, {n_cells} cells): "
          f"anterior {n_a:+.4f} mm, posterior {n_p:+.4f} mm (both |.| < 0.2) -> "
          f"{'PASS' if ok2 else 'FAIL'}")

    # ---- known answer 2b: the null above CANNOT FAIL -- it compares an array with its own
    # copy, so it tests determinism and nothing else. The bias it should have tested is real:
    # a per-cell MAXIMUM rises with the number of samples in the cell, so the denser of two
    # meshes reads as more anterior. This body carries 65,478 vertices against the subjects'
    # 88,611-97,696, i.e. the SUBJECTS are denser, so any such bias pushes `body - subject`
    # negative and works AGAINST the hypothesis. Measured on a thinning ladder against self.
    rng = np.random.default_rng(0)
    ladder = []
    for frac in (0.5, 0.25, 0.1):
        keep = rng.random(len(Vb)) < frac
        d_a, d_p, _ = compare(Vb, Vb[keep], cells, ax)
        ladder.append(dict(fraction=frac, n=int(keep.sum()), anterior_mm=d_a, posterior_mm=d_p))
        print(f"KNOWN ANSWER 2b DENSITY, dense vs {frac:.0%}-thinned SELF ({keep.sum():,} vtx): "
              f"anterior {d_a:+.3f} mm, posterior {d_p:+.3f} mm")
    print("     this body is the SPARSER mesh (0.67-0.74x the subjects'), so the bias is "
          "conservative here")

    # ---- known answer 1: a +10 mm anterior translation must move BOTH by +10
    Vt = Vb.copy(); Vt[:, ax_ap] += 0.010
    t_a, t_p, _ = compare(Vt, Vb, cells, ax)
    ok1 = abs(t_a - 10.0) < 0.5 and abs(t_p - 10.0) < 0.5
    print(f"KNOWN ANSWER 1  +10.0 mm ANTERIOR TRANSLATION: anterior {t_a:+.3f} mm, "
          f"posterior {t_p:+.3f} mm (both 10.0 +/- 0.5) -> {'PASS' if ok1 else 'FAIL'}")
    if not (ok1 and ok2):
        raise SystemExit("a known answer FAILED; this measurement is void")

    rows = []
    for sid in sorted(REG):
        M = np.array(json.loads((ROOT / "data/derived" / REG[sid] / "manifest.json").read_text())["transform"])
        Vs = source_vertices(sid, pairs) @ M[:3, :3].T + M[:3, 3]
        cells_s = footprint_cells(sid, ax_ml, ax_si)
        d_ant, d_post, n = compare(Vb, Vs, cells_s, ax)
        rows.append(dict(subject=sid, n_cells=n, anterior_mm=d_ant, posterior_mm=d_post,
                         measured_anterior_mm=MEASURED[sid],
                         fraction=d_ant / MEASURED[sid],
                         localisation=abs(d_ant) / max(abs(d_post), 1e-9)))

    print(f"\n{'subject':8s} {'cells':>6s} {'anterior':>10s} {'posterior':>10s} {'measured':>9s} "
          f"{'fraction':>9s} {'ant/post':>9s}")
    for r in rows:
        print(f"{r['subject']:8s} {r['n_cells']:6d} {r['anterior_mm']:+10.2f} "
              f"{r['posterior_mm']:+10.2f} {r['measured_anterior_mm']:+9.2f} "
              f"{r['fraction']:9.1%} {r['localisation']:9.2f}")

    n_frac = sum(r["fraction"] >= FRACTION_MIN for r in rows)
    n_loc = sum(r["localisation"] >= LOCALISATION_MIN for r in rows)
    confirmed = n_frac >= N_MIN and n_loc >= N_MIN
    print(f"\n  (i)  anterior difference reaches {FRACTION_MIN:.0%} of the measured correction "
          f"in {n_frac} of 4 (need >= {N_MIN})")
    print(f"  (ii) anterior is >= {LOCALISATION_MIN:.0f}x posterior in {n_loc} of 4 "
          f"(need >= {N_MIN})")
    print("\n  VERDICT: " + (
        "CONFIRMED -- the ~19 mm is THIS BODY'S OWN anterior thoracic shape, not a registration\n"
        "  defect. No registration fixes it: the honest remedies are a different chest wall or an\n"
        "  explicit anterior-chest correction, declared as authored."
        if confirmed else
        "REFUTED on the gate as written. The fractions above stand as they are and the gate is\n"
        "  not moved; the next candidate is named rather than this one rescued."))

    p = ROOT / a.out
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(
        gate=dict(fraction_min=FRACTION_MIN, localisation_min=LOCALISATION_MIN, n_min=N_MIN),
        confirmed=bool(confirmed), n_fraction_met=int(n_frac), n_localisation_met=int(n_loc),
        known_answers=dict(null_anterior_mm=n_a, null_posterior_mm=n_p,
                           translation_anterior_mm=t_a, translation_posterior_mm=t_p,
                           density_ladder=ladder,
                           null_cannot_fail="the self-vs-copy null tests determinism only; the "
                                             "density ladder is the null that can fail"),
        rows=rows,
        caveat="four CT subjects registered onto a male-derived body; a statement about the "
               "registration and this body's thorax, not about any subject's anatomy"), indent=2) + "\n")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
