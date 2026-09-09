#!/usr/bin/env python3
"""The moment the derived ligaments actually make about each declared coordinate.

The brief's gate was: with ligaments carrying load, the joint stops should become
less necessary.  Measured on a 2 s prone drop they do not -- the worst excursion
past the model's declared ranges goes 17.4 deg without ligaments to 32.8 deg with
them, and adding them to the stops takes 5.3 deg to 10.6 deg.  This is the
measurement that says why, coordinate by coordinate, and it is arithmetic rather
than a rollout so it is exact and free.

A line element under tension `F` makes a generalized force `-F dl/dq` on
coordinate `q`.  Summing that over the 117 derived elements gives the passive
moment the ligament set makes about each coordinate, as a function of that
coordinate, with everything else held at the binding's reference pose.  Put
beside the `CoordinateLimitForce` joint stop's own moment at the same excursion,
it says for each joint whether a ligament restraint exists at all and where it
starts.

Two things this separates that the rollout mixed:

* **A coordinate the ligaments restrain.** The moment rises steeply outside some
  range; where that range sits relative to the model's DECLARED range is the
  whole question, because the ligament's slack length encodes the ANATOMY's pose
  and the declared range is a modelling convention of the Rajagopal source.
  Where they disagree, a ligament pulls the joint out of its declared range and
  a stop pushes it back in.
* **A coordinate the ligaments cannot restrain.** If no derived element has a
  moment arm about that axis, the moment is flat at zero however far the
  coordinate goes, and no amount of ligament stiffness will hold it.

The forward kinematics is this repository's own `OsimModel.forward`, gated to
3.3e-16 m against the engine's own `Blankevoort1991Ligament` path lengths in
`scripts/measure_tissue_mechanics.py`.  The force law is Blankevoort's, written
out here and checked against the engine at the same states.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BUNDLE = ROOT / "data/models/engineering_stance_v1"
TISSUE = ROOT / "data/derived/tissue-force-elements-v1"
BINDING = ROOT / "data/derived/anatomy-segment-binding/binding.json"
STEP_RAD = 1e-5
MARGIN_RAD = math.radians(30.0)
SAMPLES = 61


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def spring_force(length, slack, stiffness, transition):
    """Blankevoort 1991: zero below slack, quadratic toe, then linear."""
    strain = (length - slack) / slack
    if strain <= 0:
        return 0.0
    if strain <= transition:
        return 0.5 * stiffness * strain * strain / transition
    return stiffness * (strain - transition / 2.0)


def lengths(model, elements, pose):
    transforms = model.forward(pose)
    out = np.empty(len(elements))
    for i, row in enumerate(elements):
        a = transforms[row["body1"]][:3, :3] @ row["_p1"] + transforms[row["body1"]][:3, 3]
        b = transforms[row["body2"]][:3, :3] @ row["_p2"] + transforms[row["body2"]][:3, 3]
        out[i] = np.linalg.norm(a - b)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True)
    parser.add_argument("--samples", type=int, default=SAMPLES)
    args = parser.parse_args()
    started = time.time()

    crawl = load_module("crawl", ROOT / "scripts/crawl.py")
    render = load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")
    model = render.OsimModel(BUNDLE / "model.osim")
    elements = json.loads((TISSUE / "ligaments.json").read_text())["elements"]
    for row in elements:
        row["_p1"] = np.asarray(row["point1_m"], float)
        row["_p2"] = np.asarray(row["point2_m"], float)
    forces = np.empty(len(elements))
    slack = np.array([r["slack_length_m"] for r in elements])
    stiffness = np.array([r["linear_stiffness_n"] for r in elements])
    transition = np.array([r["transition_strain"] for r in elements])

    reference = dict(json.loads(BINDING.read_text())["reference_pose_rad"])
    ranges = crawl.declared_ranges()
    stop = crawl.JOINT_STOP

    rows = []
    for name, (low, high) in sorted(ranges.items()):
        values = np.linspace(low - MARGIN_RAD, high + MARGIN_RAD, args.samples)
        samples = []
        for value in values:
            pose = dict(reference)
            pose[name] = float(value)
            base = lengths(model, elements, pose)
            pose[name] = float(value) + STEP_RAD
            plus = lengths(model, elements, pose)
            derivative = (plus - base) / STEP_RAD
            for i in range(len(elements)):
                forces[i] = spring_force(base[i], slack[i], stiffness[i], transition[i])
            moment = float(-(forces * derivative).sum())
            past = max(low - value, value - high, 0.0)
            samples.append(dict(value_rad=float(value), ligament_moment_nm=moment,
                                past_declared_range_rad=float(past),
                                joint_stop_moment_nm=float(
                                    -math.copysign(stop["stiffness_nm_per_rad"] * past,
                                                   value - (high if value > high else low)))
                                if past > 0 else 0.0,
                                elements_loaded=int((forces > 1e-9).sum())))
        # Does any element move at all with this coordinate?
        moving = int(np.count_nonzero(
            np.abs(lengths(model, elements, {**reference, name: float(values[0])})
                   - lengths(model, elements, {**reference, name: float(values[-1])})) > 1e-6))
        outside = [s for s in samples if s["past_declared_range_rad"] > 0]
        at_limit = max((abs(s["ligament_moment_nm"]) for s in outside), default=0.0)
        rows.append(dict(coordinate=name, low_rad=low, high_rad=high,
                         reference_value_rad=reference.get(name),
                         reference_inside_declared_range=(
                             None if name not in reference else bool(low <= reference[name] <= high)),
                         elements_whose_length_changes=moving,
                         maximum_ligament_moment_outside_range_nm=at_limit,
                         joint_stop_moment_at_same_excursion_nm=float(
                             stop["stiffness_nm_per_rad"] * max(
                                 (s["past_declared_range_rad"] for s in outside), default=0.0)),
                         samples=samples))
        print("  %-20s elements moving %3d   max |ligament moment| outside range %8.2f N.m   "
              "stop at same excursion %6.2f N.m%s"
              % (name, moving, at_limit,
                 rows[-1]["joint_stop_moment_at_same_excursion_nm"],
                 "" if rows[-1]["reference_inside_declared_range"] in (True, None)
                 else "   REFERENCE POSE IS OUTSIDE THE DECLARED RANGE"), flush=True)

    silent = [r["coordinate"] for r in rows if r["elements_whose_length_changes"] == 0]
    weak = [r["coordinate"] for r in rows
            if r["elements_whose_length_changes"] > 0
            and r["maximum_ligament_moment_outside_range_nm"]
            < r["joint_stop_moment_at_same_excursion_nm"]]
    disagree = [r["coordinate"] for r in rows if r["reference_inside_declared_range"] is False]

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(
        schema="ihm.ligament-moment-measurement.v1",
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        wall_s=time.time() - started,
        elements=len(elements), joint_stop=stop,
        held_pose="the binding's reference pose, with one coordinate swept",
        finite_difference_step_rad=STEP_RAD, margin_rad=MARGIN_RAD,
        coordinates_with_no_ligament_spanning_them=silent,
        coordinates_where_the_stop_is_stronger_than_the_ligaments=weak,
        coordinates_whose_reference_pose_is_outside_the_declared_range=disagree,
        coordinates=rows), indent=2) + "\n")
    print()
    print("coordinates no derived ligament spans:            %d  %s" % (len(silent), silent))
    print("coordinates where the 30 N.m/rad stop is stronger: %d  %s" % (len(weak), weak))
    print("coordinates whose reference pose is out of range:  %d  %s" % (len(disagree), disagree))
    print(str(out))


if __name__ == "__main__":
    main()
