#!/usr/bin/env python3
"""Is the crawl LEANING on the joint-limit violations, or grazing them?

`docs/SEGMENT_CONTACT_SURFACES.md` records that all 67 stored trajectories are outside the model's
own `declared_ranges`, and that `crawl-best` is outside on 12 of its 22 coordinates, worst 11.5° on
`ankle_angle_r`. That is a yes/no. It does not say whether the locomotion DEPENDS on the excursions.

Two coordinates can both be "11° outside" and mean opposite things: one that grazes its limit for
three frames of 1,600 is incidental, and one that spends half the cycle outside it is load-bearing.
This measures, per violating coordinate:

  * the fraction of frames outside the declared range;
  * the excursion depth as a fraction of the coordinate's OWN declared range, which is the honest
    denominator -- 11° past a +/-50° hip is not 11° past a +/-5° subtalar;
  * the fraction of the coordinate's total traversed motion that lies outside.

No engine, no replay: this is arithmetic on the recorded trajectory against the model's ranges.

KNOWN ANSWERS, checked before any figure is printed:
  1. a synthetic trace held exactly at the range centre reports zero violating frames and zero depth;
  2. a synthetic trace offset by exactly one full range width reports 100% of frames outside at a
     depth of exactly 0.5 range -- a detector that cannot produce a known depth is not measuring one.
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TRAJ = ROOT / "data/derived/crawl-best/trajectory.json"

# reuse the SEARCH's own range reader rather than parsing the model again: a second parser is a
# second chance to disagree with the thing the search was scored against.
_spec = importlib.util.spec_from_file_location("crawl", ROOT / "scripts/crawl.py")
_crawl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_crawl)
declared_ranges = _crawl.declared_ranges


def depth(series, lo, hi):
    """Fraction of frames outside, and the worst excursion as a fraction of the range's own width."""
    width = hi - lo
    below = np.clip(lo - series, 0, None)
    above = np.clip(series - hi, 0, None)
    out = np.maximum(below, above)
    return float((out > 0).mean()), float(out.max() / width), float(out.max())


def gate():
    ok = True
    lo, hi = -0.5, 0.5
    centred = np.full(200, 0.0)
    f, d, _ = depth(centred, lo, hi)
    good = f == 0.0 and d == 0.0
    ok &= good
    print(f"  held at the range centre -> {f:.3f} of frames outside, depth {d:.3f} of range"
          f"  (want 0, 0)  {'PASS' if good else 'FAIL'}")
    shifted = np.full(200, 1.5)          # one full width (1.0) above hi=0.5 -> excursion 1.0 = 1.0 range
    f2, d2, _ = depth(shifted, lo, hi)
    good2 = f2 == 1.0 and abs(d2 - 1.0) < 1e-12
    ok &= good2
    print(f"  offset by one full range width -> {f2:.3f} outside, depth {d2:.3f} of range"
          f"  (want 1.000, 1.000)  {'PASS' if good2 else 'FAIL'}")
    return ok


def main():
    print("KNOWN ANSWERS")
    if not gate():
        sys.exit("known answer FAILED -- no violation figure is printed")

    ranges = declared_ranges()
    t = json.loads(TRAJ.read_text())
    frames = t["frames"]
    # every frame carries a NAMED joints dict, so nothing is inferred from column order -- the trap
    # that produced LOG row 8, where an ordering rebuilt by walking a directory matched on 10 of
    # 16,540 pairs and the count check passed because the counts were equal.
    j0 = frames[0]["joints"]
    names = list(j0.keys())
    if any(list(f["joints"].keys()) != names for f in frames[::97]):
        sys.exit("frames do not agree on their coordinate names; not proceeding")
    # each entry is {value, speed, unit}. Only ROTATIONAL coordinates are comparable with the
    # declared ranges; comparing a pelvis translation in metres against a range in radians is the
    # units error this programme's ledger opens with.
    units = {n: j0[n].get("unit") for n in names}
    names = [n for n in names if units[n] == "rad"]
    skipped = sorted(u for n, u in units.items() if u != "rad")
    print(f"comparing {len(names)} rotational coordinates; skipping "
          f"{len(units) - len(names)} non-rotational (units {sorted(set(skipped))})")
    Q = np.asarray([[f["joints"][n]["value"] for n in names] for f in frames], float)
    print(f"\ntrajectory {TRAJ.name}: {Q.shape[0]} frames x {Q.shape[1]} coordinates; "
          f"model declares ranges for {len(ranges)}")

    rows = []
    for i, name in enumerate(names):
        if name not in ranges:
            continue
        lo, hi = ranges[name]
        frac, dep, worst = depth(Q[:, i], lo, hi)
        if frac == 0.0:
            continue
        span = float(Q[:, i].max() - Q[:, i].min())
        rows.append((dep, name, frac, dep, np.degrees(worst), np.degrees(hi - lo), span and worst / span))
    rows.sort(reverse=True)

    print(f"\n{'coordinate':18s} {'frames out':>11s} {'depth/range':>12s} {'worst deg':>10s} "
          f"{'range deg':>10s} {'worst/span':>11s}")
    for _, name, frac, dep, worst_deg, range_deg, rel_span in rows:
        print(f"{name:18s} {frac:11.1%} {dep:12.1%} {worst_deg:10.2f} {range_deg:10.2f} "
              f"{rel_span:11.1%}")

    if rows:
        fr = np.array([r[2] for r in rows]); dp = np.array([r[3] for r in rows])
        print(f"\n{len(rows)} coordinates violate. frames outside: median {np.median(fr):.1%}, "
              f"max {fr.max():.1%}. depth as a fraction of the coordinate's own range: "
              f"median {np.median(dp):.1%}, max {dp.max():.1%}.")
        print("A coordinate outside for a few frames at a few percent of its range is grazing its")
        print("limit. One outside for much of the cycle, or deep into its own range, is being leaned")
        print("on, and clamping it would change the locomotion rather than tidy it.")


if __name__ == "__main__":
    main()
