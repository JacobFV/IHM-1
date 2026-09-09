"""a library of forced pose trajectories, and what the real muscles feel in each.

`collect_forced_gait_corpus.py` did this for one motion.  the strategy needs
many: every reasonable pose trajectory becomes 3D-animation-style forced motion
on the body, and the body-state trajectory it produces is pretraining data.
reinforcement learning takes over afterwards and has to learn to move with the
real muscles -- but it starts from a periphery that has already seen what these
motions feel like, instead of from nothing.

**prescribed kinematics need no contact, and that is the point.**  the body
currently has contact on two feet and a fall plane and nothing else, so prone
motion cannot be SIMULATED (docs/DISCONNECTS.md item 9).  it can still be
IMPOSED: a pose trajectory is read at each frame and the muscle geometry
evaluated there.  so a crawl corpus is collectable today on a body that cannot
yet physically crawl, and the contact work only blocks the reinforcement phase.

what is read per frame, per muscle, from the analytic path polynomials:

  path length          musculotendon path, metres
  excursion            that length in units of the muscle's own optimal fibre
                       length -- the spindle-relevant quantity, and NOT a
                       normalised fibre length (tendon is a quarter of the path;
                       see the sibling script for the measurement)
  rate                 its time derivative
  moment arms          per coordinate, so the torque a unit activation produces
                       is recoverable

sources are OpenSim's own .mot/.sto motion files.  ground-reaction files are
excluded by name -- a GRF file has force columns, not coordinates, and reading
one as a pose would silently produce nonsense.  every corpus records the file it
came from and its hash.
"""
from __future__ import annotations

import argparse, hashlib, json, math, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import FittedMomentArms

PATHSET = ("data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/"
           "example3DWalking/subject_walk_scaled_FunctionBasedPathSet.xml")
#: files whose columns are forces or markers, never coordinates.  reading one as
#: a pose is silent nonsense, so exclude by name rather than hoping the column
#: check catches it.
EXCLUDE = ("_grf", "_forces", "_marker", "_GRF", "ReactionLoads", "_states",
           # these carry coordinate-NAMED columns holding something that is not
           # an angle: dudt is acceleration, _u is velocity, power is watts.  the
           # excursion gate below caught them at 1e22 optimal fibre lengths, but
           # a name filter is cheaper than parsing them first.
           "_dudt", "_Actuation", "_power", "_Kinematics_u", "_speeds",
           "_Velocity", "_Acceleration", "_controls", "_activation")


def read_motion(path: Path):
    """OpenSim storage/motion: header to `endheader`, then names, then rows."""
    lines = path.read_text(errors="replace").splitlines()
    try:
        end = next(i for i, l in enumerate(lines) if l.strip().lower() == "endheader")
    except StopIteration:
        return None
    meta = dict(l.split("=", 1) for l in lines[:end] if "=" in l)
    names = lines[end + 1].split()
    rows = []
    for l in lines[end + 2:]:
        if not l.strip():
            continue
        try:
            v = [float(x) for x in l.split()]
        except ValueError:
            return None
        if len(v) == len(names):
            rows.append(v)
    if len(rows) < 2 or not names or names[0].lower() != "time":
        return None
    data = np.asarray(rows)
    if not np.isfinite(data).all():
        return None
    # degrees vs radians is declared in the header and getting it wrong is a
    # 57x error that still looks like a plausible pose.
    deg = str(meta.get("inDegrees", "")).strip().lower() == "yes"
    return names, data, deg, meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--search", default="data/raw/anatomy/opensim-models")
    ap.add_argument("--extra-search", default="data/raw/mechanics")
    ap.add_argument("--min-coords", type=int, default=6,
                    help="motions matching fewer model coordinates are skipped")
    ap.add_argument("--max-frames", type=int, default=400)
    ap.add_argument("--mass", type=float, default=77.6122029)
    ap.add_argument("--out", default="data/derived/pose-corpus")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    out = ROOT / a.out
    out.mkdir(parents=True, exist_ok=True)

    stream = NativeMechanicalStream(ROOT, out / f"plant-{int(time.time())}",
                                    environment="upright", target_mass_kg=a.mass)
    s0 = stream.snapshot()
    muscles, model_coords = s0["muscles"], set(s0["coordinates"])
    geom = FittedMomentArms(ROOT / PATHSET)
    fitted = [m for m in muscles if m in geom.paths]
    stream.close()
    print(f"model: {len(model_coords)} coordinates, {len(fitted)} muscles with "
          f"a fitted path", flush=True)

    files = []
    for base in (a.search, a.extra_search):
        d = ROOT / base
        if d.exists():
            files += [p for p in d.rglob("*.mot")] + [p for p in d.rglob("*.sto")]
    files = sorted({p for p in files if not any(x in p.name for x in EXCLUDE)})
    print(f"candidate motion files: {len(files)}", flush=True)

    index, kept, skipped = [], 0, {}
    for p in files:
        if a.limit and kept >= a.limit:
            break
        r = read_motion(p)
        if r is None:
            skipped["unparseable"] = skipped.get("unparseable", 0) + 1
            continue
        names, data, deg, meta = r
        cols = {n: i for i, n in enumerate(names) if n in model_coords}
        if len(cols) < a.min_coords:
            skipped["too few coordinates"] = skipped.get("too few coordinates", 0) + 1
            continue

        t = data[:, 0]
        step = max(1, len(t) // a.max_frames)
        frames, prev = [], None
        for r_i in range(0, len(t), step):
            q = {}
            for n, i in cols.items():
                v = float(data[r_i, i])
                q[n] = math.radians(v) if (deg and "_t" not in n[-3:]) else v
            length, arms = {}, {}
            for m in fitted:
                try:
                    L, A = geom.length_and_moment_arms(m, q)
                except KeyError:
                    continue
                length[m], arms[m] = L, A
            if not length:
                break
            exc = {m: length[m] / float(muscles[m]["optimal_fiber_length_m"])
                   for m in length}
            dt = (float(t[r_i]) - frames[-1]["time_s"]) if frames else 0.0
            rate = ({m: (exc[m] - prev[m]) / dt for m in exc if m in prev}
                    if prev is not None and dt > 0 else {m: 0.0 for m in exc})
            prev = exc
            frames.append({"time_s": float(t[r_i]), "coordinates": q,
                           "path_length_m": length,
                           "path_over_optimal_fiber": exc,
                           "rate_per_s": rate, "moment_arms": arms})
        if len(frames) < 2:
            skipped["no muscle path resolved"] = skipped.get("no muscle path resolved", 0) + 1
            continue

        # the excursion each muscle actually undergoes.  a motion in which every
        # muscle is static teaches nothing, so this is the corpus's own quality
        # measure and it is recorded per motion rather than assumed.
        # A CASE WHOSE ANSWER IS KNOWN.  a musculotendon path is a few times its
        # own optimal fibre length -- never a thousandth of it and never a
        # million times it.  files holding accelerations or powers under
        # coordinate names parse cleanly and produce excursions of 1e22, which is
        # how they were caught.  gate on physiology, not just on filenames.
        allv = [v for f in frames for v in f["path_over_optimal_fiber"].values()]
        lo_all, hi_all = min(allv), max(allv)
        if not (0.2 < lo_all and hi_all < 20.0):
            skipped["implausible path length (not a pose file)"] = skipped.get(
                "implausible path length (not a pose file)", 0) + 1
            continue

        rng = {m: (min(f["path_over_optimal_fiber"][m] for f in frames),
                   max(f["path_over_optimal_fiber"][m] for f in frames))
               for m in frames[0]["path_over_optimal_fiber"]}
        swing = {m: hi - lo for m, (lo, hi) in rng.items()}
        total = float(sum(swing.values()))
        name = p.stem
        rec = {
            "id": name, "source": str(p.relative_to(ROOT)),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "in_degrees": deg, "n_frames": len(frames),
            "duration_s": float(t[-1] - t[0]),
            "coordinates_matched": sorted(cols), "n_coordinates": len(cols),
            "muscles": len(frames[0]["path_over_optimal_fiber"]),
            "total_excursion_optimal_fiber": round(total, 3),
            "largest_excursions": {m: round(swing[m], 4) for m in
                                   sorted(swing, key=swing.get, reverse=True)[:5]},
        }
        (out / f"{name}.json").write_text(json.dumps(
            {**rec, "prescribed": "forced pose trajectory; no contact, no "
                                  "dynamics. this is pretraining data, not "
                                  "locomotion.", "frames": frames}))
        index.append(rec)
        kept += 1
        print(f"  {name:44s} {len(frames):4d} fr  {len(cols):2d} coords  "
              f"excursion {total:7.2f}", flush=True)

    index.sort(key=lambda r: -r["total_excursion_optimal_fiber"])
    (out / "index.json").write_text(json.dumps(
        {"n_motions": len(index), "skipped": skipped,
         "note": "prescribed pose trajectories. muscle geometry is analytic from "
                 "the fitted path polynomials; path length is MUSCULOTENDON and "
                 "the excursion, not the absolute value, is the spindle-relevant "
                 "quantity.", "motions": index}, indent=1))
    print(f"\n  kept {len(index)} motions, skipped {skipped}", flush=True)
    print(f"  richest by total muscle excursion:", flush=True)
    for r in index[:8]:
        print(f"    {r['id']:44s} {r['total_excursion_optimal_fiber']:8.2f}  "
              f"{r['n_frames']:4d} fr  {r['duration_s']:6.2f}s", flush=True)
    print(f"  wrote {out}/index.json", flush=True)


if __name__ == "__main__":
    main()
