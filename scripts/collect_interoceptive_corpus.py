#!/usr/bin/env python3
"""Collect a visceral afferent corpus from recorded native systemic trajectories.

The twin of `collect_body_corpus.py`, for the other half of the body.  That one
makes movement one more paired corpus in IBM-1's supervision soup --
(muscle state, motor command).  This one makes the viscera one:
(visceral afferent rates, the physiological scalars they are evidence about).

**The physics has already run.**  `ihm/assembly/systemic.py` recorded seven
BioGears trajectories under five protocols -- rest, hydration, meal, exercise,
meal+exercise -- at 5 s and 30 s cadence, each 180 native quantities wide.  They
are on disk.  Re-running them would cost hours of a machine that is oversubscribed
and would produce the same numbers, so this reads them and transduces.  Nothing
here integrates physiology and nothing here invents a trajectory; every rate is
`ihm.assembly.interoception` applied to a frame BioGears wrote.

**Why several protocols and not one long rest run.**  A corpus recorded at
equilibrium has almost no variance and cloning it is cloning a constant -- the
lesson `collect_body_corpus.py` records for the motor side, where an unperturbed
stance corpus made every ablation arm lose to predicting the mean.  The visceral
equivalent of a perturbation is a meal and a bout of exercise, and the protocol
set is chosen so that both the substrate and the demand move: metabolic rate
spans 82-232 W, arterial lactate 1.5-86 mg/dL, gastric contents 0-1098 mL.

**The split is defined HERE, not by the consumer.**  `sensation` is a projection
whose basis must be fitted without seeing test data, so the builder fits it and
writes it down alongside the split it was fitted on.  A consumer that chose its
own split would silently fit the basis on data it then evaluated against.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.interoception import (  # noqa: E402
    CHANNELS, SCALAR_NAMES, channel_table, derived_scalars, source_keys,
    visceral_afferent_rates)

#: the recorded runs this reads.  named explicitly rather than globbed: a glob
#: would silently pick up a future run recorded under a different engine variant
#: or protocol and mix it into the corpus without anything saying so.
SOURCES = (
    ("exertion_v3", "rest"), ("exertion_v3", "hydration"),
    ("exertion_v3", "meal"), ("exertion_v3", "exercise"),
    ("exertion_v3", "meal_exercise"),
    ("six_hour_v4", "meal"), ("six_hour_v4", "hydration"),
)


def read_sequence(run: str, protocol: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """One recorded trajectory -> (afferent rates, scalars, times, dt)."""
    path = ROOT / "data/derived/systemic" / run / protocol / "frames.jsonl"
    if not path.exists():
        raise SystemExit(f"no recorded trajectory at {path}")
    keys = None
    rates, scal, times = [], [], []
    for line in path.open():
        frame = json.loads(line)
        v = frame["values"]
        r = visceral_afferent_rates(v)          # raises on any absent/null source
        if keys is None:
            keys = sorted(r)
        elif sorted(r) != keys:
            raise ValueError(f"{run}/{protocol}: channel set changed mid-run")
        d = derived_scalars(v, r)
        rates.append([r[k] for k in keys])
        scal.append([d[s] for s in SCALAR_NAMES])
        times.append(frame["time_s"])
    t = np.asarray(times, np.float64)
    dt = float(np.median(np.diff(t)))
    if not np.allclose(np.diff(t), dt, rtol=1e-6, atol=1e-6):
        raise ValueError(f"{run}/{protocol}: non-uniform sample clock")
    return (np.asarray(rates, np.float32), np.asarray(scal, np.float32), t, dt)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/derived/intero-corpus")
    ap.add_argument("--also", default="",
                    help="a second directory to write the same arrays into, so "
                         "the consuming repo does not have to reach across a "
                         "filesystem at training time")
    ap.add_argument("--train-frac", type=float, default=0.8)
    ap.add_argument("--guard-frac", type=float, default=0.02,
                    help="fraction of each sequence discarded at the split "
                         "boundary.  visceral state at 5 s cadence is heavily "
                         "autocorrelated and a test sample one step after a "
                         "training sample is not independent of it")
    ap.add_argument("--sensation-dim", type=int, default=3)
    a = ap.parse_args()

    out = ROOT / a.out
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)

    keys = [f"{c.trunk}/{c.name}" for c in CHANNELS]
    keys_sorted = sorted(keys)

    X, Y, T, SEQ, meta_seq = [], [], [], [], []
    for i, (run, protocol) in enumerate(SOURCES):
        x, y, t, dt = read_sequence(run, protocol)
        X.append(x); Y.append(y); T.append(t)
        SEQ.append(np.full(len(x), i, np.int32))
        meta_seq.append({"index": i, "run": run, "protocol": protocol,
                         "n": int(len(x)), "dt_s": dt,
                         "duration_s": float(t[-1] - t[0])})
        print(f"  {run}/{protocol:14s} {len(x):5d} frames  dt {dt:5.1f} s  "
              f"{t[-1]-t[0]:7.0f} s", flush=True)
    X = np.concatenate(X); Y = np.concatenate(Y)
    T = np.concatenate(T); SEQ = np.concatenate(SEQ)

    # contiguous split WITHIN each sequence, with a guard band.  splitting
    # across sequences instead would put whole protocols on one side, which is a
    # different and harder question -- worth asking, but not the same one, and
    # mixing the two is how a result stops meaning anything.
    is_train = np.zeros(len(X), bool)
    is_test = np.zeros(len(X), bool)
    for s in meta_seq:
        m = SEQ == s["index"]
        n = int(m.sum())
        cut = int(n * a.train_frac)
        guard = max(1, int(n * a.guard_frac))
        idx = np.nonzero(m)[0]
        is_train[idx[:cut - guard]] = True
        is_test[idx[cut + guard:]] = True

    # a channel that never moves anywhere in the corpus is reported, not dropped.
    dead = [k for j, k in enumerate(keys_sorted) if X[:, j].std() == 0.0]

    # `sensation`: a low-dimensional linear projection of the afferent vector,
    # basis fitted on the TRAIN split only.  standardised first, because the
    # channels differ in max rate by 3x and an unstandardised PCA would report
    # the loudest channel rather than the dominant pattern.
    mu = X[is_train].mean(0)
    sd = X[is_train].std(0)
    sd_safe = np.where(sd > 0, sd, 1.0)
    Z = (X - mu) / sd_safe
    u, s, vt = np.linalg.svd(Z[is_train] - Z[is_train].mean(0), full_matrices=False)
    basis = vt[:a.sensation_dim].astype(np.float32)          # (dim, n_channels)
    var_explained = (s[:a.sensation_dim] ** 2 / (s ** 2).sum()).tolist()
    S = (Z @ basis.T).astype(np.float32)

    np.save(out / "afferent.npy", X)
    np.save(out / "scalars.npy", Y)
    np.save(out / "sensation.npy", S)
    np.save(out / "time.npy", T.astype(np.float32))
    np.save(out / "sequence.npy", SEQ)
    np.save(out / "is_train.npy", is_train)
    np.save(out / "is_test.npy", is_test)

    meta = {
        "schema": "ihm.intero-corpus.v1",
        "channels": keys_sorted,
        "n_channels": len(keys_sorted),
        "channel_table": channel_table(),
        "source_keys": list(source_keys()),
        "scalars": list(SCALAR_NAMES),
        "sensation_dim": a.sensation_dim,
        "sensation_basis": basis.tolist(),
        "sensation_standardize_mean": mu.tolist(),
        "sensation_standardize_sd": sd.tolist(),
        "sensation_variance_explained": var_explained,
        "sensation_basis_fitted_on": "train split only",
        "sequences": meta_seq,
        "n": int(len(X)),
        "n_train": int(is_train.sum()),
        "n_test": int(is_test.sum()),
        "train_frac": a.train_frac, "guard_frac": a.guard_frac,
        "silent_channels": dead,
        "engine": "BioGears via ihm/assembly/systemic.py, variant "
                  "whole_body_integrity_depletion; frames recorded 2026-09-05",
        "physics_ran_at_collection": True,
        "scalar_semantics":
            "endurance_h and discomfort are NAMED PROJECTIONS OF MEASURED "
            "PHYSIOLOGICAL STATE, not claims about experience.  endurance_h is "
            "hours of carbohydrate substrate at the current metabolic rate; "
            "discomfort is a fibre-class-weighted sum of the afferent rates in "
            "this same file, normalised by its own saturated ceiling.  "
            "sensation is a linear projection of the afferent vector.  none of "
            "the three is evidence about anything felt.",
        "limitations": [
            "Receptor laws (threshold, saturation, max rate) are authored "
            "anatomical priors, not measured from these organs.",
            "Native gut is a lumped stomach/small-intestine nutrient model; no "
            "chewing, peristalsis, colon or stool process.",
            "The vagal and splanchnic routes carry no native transduction; this "
            "corpus is the first thing to put a rate on them.",
            "Gastric volume is derived from stomach water plus macronutrient "
            "mass at nominal densities; BioGears exposes no gastric volume.",
        ],
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")

    print(f"\n{len(X)} frames, {len(keys_sorted)} channels, "
          f"{len(SCALAR_NAMES)} scalars, sensation dim {a.sensation_dim}")
    print(f"  train {is_train.sum()}  test {is_test.sum()}  "
          f"(guard band discards {len(X) - is_train.sum() - is_test.sum()})")
    print(f"  afferent variance (per channel, corpus)  "
          f"min {X.var(0).min():.4f}  max {X.var(0).max():.4f}")
    print(f"  sensation variance explained: "
          f"{', '.join(f'{v:.3f}' for v in var_explained)}")
    if dead:
        print(f"  SILENT across the whole corpus: {', '.join(dead)}")
        print("    a high-threshold channel the protocols never drove into "
              "range.  reported, not dropped -- dropping it would hide that "
              "the corpus does not exercise it.")
    for s_ in SCALAR_NAMES:
        j = list(SCALAR_NAMES).index(s_)
        print(f"  {s_:14s} min {Y[:, j].min():9.4f}  max {Y[:, j].max():9.4f}  "
              f"sd {Y[:, j].std():9.4f}")
    print(f"wrote {out}")

    if a.also:
        dst = Path(a.also).expanduser()
        shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(out, dst)
        print(f"copied to {dst}")


if __name__ == "__main__":
    main()
