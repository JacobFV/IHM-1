#!/usr/bin/env python3
"""Is recorded human motion admissible under THIS body's declared joint ranges?

IBM-1 `docs/LOG.md`, 2026-09-11, fixed the gates before any of this was written. The claim
under test is the one the mocap corpus was fetched on -- *"recorded human motion is admissible
by construction"* -- and it is load-bearing, because this body has no admissible motion of its
own: **0 of 67** stored crawl trajectories respect `crawl.declared_ranges()`, and crawl-best's
left knee is outside its limit for 95.6% of 1,600 frames.

The route taken here is better than the ASF/AMC retargeting the pre-registration anticipated,
because `data/derived/pose-corpus/` already holds 68 recorded motions -- walking, running,
jumping, crouching -- **already expressed in OpenSim coordinate space**, with sha256
provenance to the `.mot`/`.sto` they came from. No retargeting step means no retargeting
confound.

IT DOES NOT MEAN NO CONFOUND. Two were found in the data before any figure was computed, and
either one alone would have manufactured a total violation out of nothing:

  SIGN.  `knee_angle_*` is negative in ~100% of frames in every gait2392-derived file, while
         `engineering_stance_v1` declares `knee_angle_* in [0, 2.443]`. gait2392 and the
         `walker_knee` used here disagree on which direction of knee rotation is positive. A
         pure convention flip reads as a 100% violation.
  UNITS. `runningModel_Kinematics_q` carries `knee_angle_r` down to **-114.0**: degrees,
         unconverted, while its siblings hold radians -- and every file's `in_degrees` flag
         describes the SOURCE `.mot`, not the JSON. The corpus is not uniformly converted.

So this keeps the pre-registration's G1 in a new guise, and G1 still gates everything:

  G1  Each motion is checked against the declared ranges of **the model it was generated
      from**, found in the same source tree as its `.mot`. OpenSim produced that motion with
      that model, so a correct reader MUST agree with it: >= 99% of (frame, coordinate) pairs
      inside. A sign flip, a degrees-for-radians read or a mis-parsed file breaks the
      agreement. **A motion that fails G1 is not scored on G2** -- its violations are mine.
  G2  Motions that pass G1 are then checked against `crawl.declared_ranges()`, this body's
      own. A motion is admissible if >= 99% of its (frame, coordinate) pairs lie inside; the
      corpus is admissible if >= 90% of motions are.
  G3  The comparison that gives G2 meaning: the stored crawl trajectories score 0 of 67 on
      exactly this rule.

WHY THE CONVENTION IS NOT INFERRED FROM THE ANSWER. The tempting shortcut is to pick, per
file, whichever sign and unit make the motion fit `declared_ranges` -- which would fit the
instrument to produce admissibility and guarantee a pass. The source model's own declaration
is used instead: it is independent of this body's ranges and was written by whoever built that
model. Where the source and this body disagree, that disagreement is REPORTED as a fact about
two models rather than resolved in either direction.

G1 IS STRUCTURALLY BLIND TO A CONVENTION MISMATCH, and the first run of this script proved it.
G1 passed 42 of 48 motions and G2 then failed at 4.8% against a 90% bar -- but the failure was
almost entirely `knee_angle_*` (40/42 and 32/42 motions, median 93% and 82% of frames outside),
and reading the two declarations side by side settles why:

    gait2392 family        knee_angle_r in [-120.0, +10.0] deg
    engineering_stance_v1  knee_angle_r in [   0.0, +140.0] deg

Mirror images. The two families disagree about which direction of knee rotation is positive,
and the motion sits inside its own model's range in ITS convention. G1 cannot see this,
because it compares each motion against a model that shares the motion's convention, so the
flip cancels on both sides. A gate that compares like with like is blind to a difference
between the two likes -- and reporting "declared_ranges is too narrow" off that first run
would have been wrong in a way the pre-registration's fork did not anticipate.

THE CONVENTION MAP, derived from the two models' DECLARATIONS and never from the motion. A
coordinate is treated as sign-flipped when the source's declared range is approximately the
negation of this body's, i.e. `lo_src ~= -hi_tgt` and `hi_src ~= -lo_tgt` within
`--flip-tol` of the span. That is a statement about two XML files, fixed before any frame is
read, and it cannot be steered by whether the answer comes out well. Symmetric ranges are
exempt because negating one changes nothing. Every flip is printed.

Per CLAUDE.md: an instrument may change after a FAIL; a threshold may not. G1, G2 and the
corpus bar are exactly as pre-registered. What changed is that the comparison now maps
conventions before comparing, which the first run showed it must.

THE FORK, from the pre-registration, decided before the numbers: G1 pass + G2 fail does not
mean the motion is rejected. Real humans produced these trajectories, so a surviving failure
indicts `declared_ranges` -- ranges that are inherited OpenSim gait2392/Moco defaults rather
than a measurement, and that `docs/NATIVE_JOINT_LIMITS.md` already records disagreeing with the
model's OWN passive stops. G2's threshold does not move in that event.

The first run added a third branch the pre-registration did not have: G1 pass + G2 fail **for a
reason that is neither the ranges nor the motion**. A convention mismatch survives G1 and looks
exactly like a range that is too narrow. It is now mapped out above, so what remains after the
map is the real disagreement -- and any coordinate that still fails is failing on its merits.
"""
import argparse, importlib.util, json, math, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
POSE = ROOT / "data/derived/pose-corpus"

_spec = importlib.util.spec_from_file_location("crawl", ROOT / "scripts/crawl.py")
_crawl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_crawl)
declared_ranges = _crawl.declared_ranges

# a rotational coordinate in radians cannot reach 2*pi on any of these joints; the widest
# range this body declares is pelvis_* at +/- pi.  this separates degrees from radians
# WITHOUT consulting the ranges under test.
RAD_CEILING = 2 * math.pi


def model_ranges(path):
    """The same parse crawl.declared_ranges uses, against an arbitrary .osim.

    Deliberately NOT filtered by the sentinel/blacklist rules: for G1 we want every range the
    source model actually declares, including ones this body drops.
    """
    text = Path(path).read_text(errors="replace")
    out = {}
    for m in re.finditer(r'<Coordinate name="([^"]+)">(.*?)</Coordinate>', text, re.S):
        r = re.search(r"<range>([^<]*)</range>", m.group(2))
        if not r:
            continue
        try:
            lo, hi = (float(v) for v in r.group(1).split())
        except ValueError:
            continue
        if hi - lo > 18.0:            # the "no range" sentinel, same bar crawl.py uses
            continue
        out[m.group(1)] = (lo, hi)
    return out


def find_source_model(source_field):
    """The .osim in the same source tree as this motion's .mot, nearest first.

    OpenSim distributions keep a pipeline's model beside or just above its results, so walking
    up from the motion file and taking the closest .osim finds the model that produced it.
    Returns None when the tree holds none, and such a motion is reported, never guessed at.
    """
    p = ROOT / source_field
    for up in list(p.parents)[:5]:
        if not up.is_dir():
            continue
        here = sorted(up.glob("*.osim"))
        if here:
            return here, up
    return None, None


def convention_map(src, tgt, tol=0.10):
    """Coordinates whose SOURCE declaration is the negation of this body's.

    Decided from the two XML declarations alone -- no motion data is consulted, so this cannot
    be steered by whether the result comes out well. A symmetric range negates to itself and is
    excluded: there is nothing to detect and nothing to correct.

    THE FIRST VERSION OF THIS DETECTED NOTHING, and the reason is worth keeping. It asked
    whether `lo_src ~= -hi_tgt` and `hi_src ~= -lo_tgt` -- that the two ranges be negations of
    EQUAL MAGNITUDE. Mirrored models need not agree on how far the joint goes: the source
    allows 120 deg of knee flexion and this body 140, so `-2.094` was compared against `-2.443`
    and missed by 0.349 rad against a 0.244 bar. The ranges are mirrored; they are not equal.

    What identifies a flip is the SIDE OF ZERO the joint's travel is on, not how far it goes.
    For an asymmetric range, take the sign of whichever bound is larger in magnitude -- the
    direction the joint actually moves in. If the source and this body disagree on that sign,
    the conventions are opposite:

        knee_angle_r  source [-120.0, +10.0] deg -> dominant bound -120.0, sign -
                      here   [   0.0, +140.0] deg -> dominant bound +140.0, sign +   FLIP

    Still two XML files and no motion data, so it cannot be steered by the answer.
    """
    flips = {}
    for name, (lo_t, hi_t) in tgt.items():
        if name not in src:
            continue
        lo_s, hi_s = src[name]
        span = max(hi_t - lo_t, 1e-9)
        if abs(lo_t + hi_t) < tol * span:        # symmetric: a flip is undetectable AND harmless
            continue
        if abs(lo_s + hi_s) < tol * max(hi_s - lo_s, 1e-9):
            continue                             # source symmetric: same reasoning
        s_src = -1 if abs(lo_s) > abs(hi_s) else 1
        s_tgt = -1 if abs(lo_t) > abs(hi_t) else 1
        if s_src != s_tgt:
            flips[name] = ((lo_s, hi_s), (lo_t, hi_t))
    return flips


def series(frames, coords):
    """(n_frames, n_coords) array of coordinate values, NaN where a frame omits one."""
    a = np.full((len(frames), len(coords)), np.nan)
    for i, f in enumerate(frames):
        c = f["coordinates"]
        for j, name in enumerate(coords):
            if name in c:
                a[i, j] = c[name]
    return a


def outside(a, coords, ranges):
    """Fraction of non-NaN (frame, coordinate) pairs outside range, and the per-coordinate detail.

    Only coordinates the range table actually declares are scored; a coordinate with no
    declared range cannot be violated and must not be counted as passing either.
    """
    tot = bad = 0
    detail = {}
    for j, name in enumerate(coords):
        if name not in ranges:
            continue
        lo, hi = ranges[name]
        v = a[:, j]
        m = ~np.isnan(v)
        if not m.any():
            continue
        o = (v[m] < lo) | (v[m] > hi)
        tot += int(m.sum()); bad += int(o.sum())
        if o.any():
            exc = np.maximum(lo - v[m], v[m] - hi).max()
            detail[name] = {"frac_out": float(o.mean()), "worst_rad": float(exc),
                            "range": [lo, hi], "span": [float(np.nanmin(v)), float(np.nanmax(v))]}
    return (bad / tot if tot else float("nan")), tot, detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g1", type=float, default=0.99, help="PRE-REGISTERED; do not change")
    ap.add_argument("--g2", type=float, default=0.99, help="PRE-REGISTERED; do not change")
    ap.add_argument("--corpus", type=float, default=0.90, help="PRE-REGISTERED; do not change")
    ap.add_argument("--flip-tol", type=float, default=0.10,
                    help="how close a source range must be to the negation of this body's to "
                         "count as a sign-convention flip, as a fraction of the span")
    ap.add_argument("--out", default="out/recorded_motion_admissibility.json")
    a = ap.parse_args()

    target = declared_ranges()
    print(f"this body declares ranges for {len(target)} coordinates "
          f"(engineering_stance_v1/model.osim, via crawl.declared_ranges)")

    idx = json.loads((POSE / "index.json").read_text())
    motions = idx["motions"]
    print(f"pose corpus: {len(motions)} motions\n")

    rows, no_model, deg_files = [], [], []
    flip_counts = defaultdict(int)
    for m in motions:
        f = POSE / f"{m['id']}.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        frames = d.get("frames") or []
        if not frames:
            continue
        coords = sorted({k for fr in frames[:50] for k in fr["coordinates"]})
        rot = [c for c in coords if not c.startswith("pelvis_t")]
        A = series(frames, coords)
        ri = [coords.index(c) for c in rot]

        # UNITS, decided without consulting either range table
        peak = np.nanmax(np.abs(A[:, ri])) if ri else 0.0
        in_deg = bool(peak > RAD_CEILING)
        if in_deg:
            A[:, ri] = np.deg2rad(A[:, ri])
            deg_files.append((m["id"], float(peak)))

        cand, tree = find_source_model(m["source"])
        if not cand:
            no_model.append(m["id"])
            continue
        # the source model is the one whose declared coordinates best cover this motion
        best, best_n = None, -1
        for c in cand:
            r = model_ranges(c)
            n = len(set(r) & set(coords))
            if n > best_n:
                best, best_n, best_r = c, n, r
        # G1 is scored in the motion's OWN convention -- untouched.
        g1_out, g1_n, g1_detail = outside(A, coords, best_r)
        # G2 maps conventions first, from the two models' declarations.
        flips = convention_map(best_r, target, a.flip_tol)
        B = A.copy()
        for name in flips:
            if name in coords:
                B[:, coords.index(name)] *= -1.0
        g2_out, g2_n, g2_detail = outside(B, coords, target)
        for name in flips:
            flip_counts[name] += 1
        rows.append({"id": m["id"], "n_frames": len(frames), "in_degrees_json": in_deg,
                     "source_model": str(best.relative_to(ROOT)), "source_coords_matched": best_n,
                     "sign_flipped": sorted(flips),
                     "g1_frac_outside_own_model": g1_out, "g1_pairs": g1_n, "g1_detail": g1_detail,
                     "g2_frac_outside_this_body": g2_out, "g2_pairs": g2_n, "g2_detail": g2_detail})

    if flip_counts:
        print("SIGN CONVENTION, read from the two models' declarations before any frame:")
        for k, n in sorted(flip_counts.items(), key=lambda kv: -kv[1]):
            print(f"  {k:20s} flipped in {n:3d} motion(s) -- the source declares the negation of "
                  f"this body's range")
        print()
    if deg_files:
        print(f"UNITS: {len(deg_files)} motion(s) held DEGREES in the JSON despite the corpus's "
              f"`in_degrees` flag describing the source file. Converted:")
        for i, p in deg_files[:6]:
            print(f"  {i[:46]:46s} peak |angle| {p:8.1f}  -> {math.radians(p):.3f} rad")
    if no_model:
        print(f"\nNO SOURCE MODEL found beside the .mot for {len(no_model)} motion(s); "
              f"not scored rather than guessed at:\n  {', '.join(no_model[:8])}"
              + (" ..." if len(no_model) > 8 else ""))

    print(f"\n{'='*78}\nG1 -- each motion against ITS OWN source model's declared ranges")
    print(f"     (OpenSim made the motion with that model, so a correct reader must agree)")
    print(f"     bar: >= {a.g1:.0%} of pairs inside\n")
    g1_pass = [r for r in rows if 1 - r["g1_frac_outside_own_model"] >= a.g1]
    g1_fail = [r for r in rows if 1 - r["g1_frac_outside_own_model"] < a.g1]
    print(f"  PASS {len(g1_pass)}/{len(rows)}   FAIL {len(g1_fail)}")
    for r in sorted(g1_fail, key=lambda r: -r["g1_frac_outside_own_model"])[:10]:
        worst = sorted(r["g1_detail"].items(), key=lambda kv: -kv[1]["frac_out"])[:2]
        print(f"    FAIL {r['id'][:40]:40s} {r['g1_frac_outside_own_model']:6.1%} outside  "
              f"[{r['source_model'].split('/')[-1]}]  worst: "
              + ", ".join(f"{k} {v['frac_out']:.0%}" for k, v in worst))
    if not g1_pass:
        print("\n  G1 PASSES FOR NOTHING. The reader disagrees with every source model, so the")
        print("  fault is here and no G2 number below is interpretable. Stopping.")
        sys.exit(1)

    print(f"\n{'='*78}\nG2 -- the {len(g1_pass)} G1-passing motions against THIS BODY's declared ranges")
    print(f"     bar: a motion is admissible if >= {a.g2:.0%} of pairs inside; "
          f"the corpus if >= {a.corpus:.0%} of motions are\n")
    adm = [r for r in g1_pass if 1 - r["g2_frac_outside_this_body"] >= a.g2]
    frac = len(adm) / len(g1_pass)
    print(f"  admissible: {len(adm)}/{len(g1_pass)} = {frac:.1%}   "
          f"{'PASS' if frac >= a.corpus else 'FAIL'} against the {a.corpus:.0%} bar")

    agg = defaultdict(list)
    for r in g1_pass:
        for k, v in r["g2_detail"].items():
            agg[k].append(v["frac_out"])
    if agg:
        print(f"\n  which coordinates fall outside, over the {len(g1_pass)} scored motions:")
        print(f"    {'coordinate':20s} {'motions':>8s}  {'median % of frames out':>22s}  {'declared range (deg)':>22s}")
        for k, v in sorted(agg.items(), key=lambda kv: (-len(kv[1]), -np.median(kv[1]))):
            lo, hi = target[k]
            print(f"    {k:20s} {len(v):>4d}/{len(g1_pass):<3d}  {np.median(v):>21.1%}  "
                  f"{math.degrees(lo):>10.1f} {math.degrees(hi):>10.1f}")

    print(f"\n  G3, for comparison: the stored crawl trajectories score 0 of 67 on this rule.")
    print(f"\n  VERDICT: " + (
        f"recorded human motion IS admissible under this body's declared ranges "
        f"({frac:.1%} >= {a.corpus:.0%})." if frac >= a.corpus else
        f"recorded human motion is NOT admissible under this body's declared ranges "
        f"({frac:.1%} < {a.corpus:.0%}).\n"
        f"  Real humans produced these trajectories and the reader agrees with each motion's own\n"
        f"  source model, so per the pre-registration this indicts `declared_ranges`, not the\n"
        f"  motion. Those ranges are inherited OpenSim gait2392/Moco defaults, not a measurement,\n"
        f"  and docs/NATIVE_JOINT_LIMITS.md already records them disagreeing with the model's own\n"
        f"  passive stops. The G2 threshold is NOT moved."))

    p = ROOT / a.out
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"gates": {"g1": a.g1, "g2": a.g2, "corpus": a.corpus},
                             "n_motions": len(rows), "g1_pass": len(g1_pass),
                             "g2_admissible": len(adm), "g2_fraction": frac,
                             "degrees_files": deg_files, "no_source_model": no_model,
                             "motions": rows}, indent=2) + "\n")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
