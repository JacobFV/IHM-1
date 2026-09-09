#!/usr/bin/env python3
"""What the tissue force elements do to the plant, and what they cost.

The plant carried force elements for none of the body's 625 tissue structures.
`scripts/build_tissue_force_elements.py` derives two-ended attachments for the
ones that span two of the scaffold's rigid bodies; this measures what happens
when they are switched on.

Two experiments, because the two questions are different.

**Statics** -- the same stance pose, the same mass, the same excitation, with
only the force set differing.  Adding internal forces must not invent or lose
external force, so:

* **Standing weight.** 761.38 N = 77.6122029 x 9.81, to the last digit.
  Ligaments are internal, so this number is not allowed to move.
* **Momentum balance.** `m(a_com - g) - contact - external`.  An internal force
  cannot appear in it.  Reported against the same arm without ligaments, because
  the residual scales with the forces in the plant and an absolute threshold
  would be comparing against the wrong thing.
* **Cost.** Median and worst wall seconds per advance.  `advance` is an
  error-controlled Simbody integration whose cost is not bounded by dt: a
  previous plant went 0.112 s -> 24.3 s when it entered a bad region.  If tissue
  makes the plant unaffordable, that is the finding.

**Range excursion** -- the question the brief asks.  Nothing in this plant
enforces the coordinate ranges the model declares; the standing joint stops are
a `CoordinateLimitForce` at a stated engineering constant, and real ligaments
and capsules are what those stops are standing in for.  So: drop the body prone,
drive nothing, and measure how far past its own declared ranges it goes -- with
ligaments, without, with stops, with both.  If ligaments do nothing, the
excursion is unchanged and that is the result.

Gates that are not the arms restated
------------------------------------
* **Reference-pose strain.** Every ligament's slack length IS its length at the
  binding's reference pose, so at that pose its strain must be exactly zero.
  Checked in closed form against this script's own forward kinematics.
* **Cross-implementation length.** The engine's reported ligament length against
  the length this script computes from the same coordinates through the model's
  own joint definitions.  Two independent implementations of the same geometry.
* **Ultimate strain.** Quapp and Weiss report human MCL failing at 17.1% strain.
  Any element past that at the run's pose is a derivation that has put an
  attachment where the scaffold cannot carry it, and is listed by name.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import re
import shutil
import statistics
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402

BUNDLE = ROOT / "data/models/engineering_stance_v1"
REGISTRATION = "data/models/engineering_stance_v1/registration.json"
TISSUE = "data/derived/tissue-force-elements-v1"
TARGET_MASS_KG = 77.6122029
WEIGHT_N = TARGET_MASS_KG * 9.81
# Quapp and Weiss 1998, the same publication the modulus comes from: human MCL
# ultimate strain 17.1 +/- 1.5%.
ULTIMATE_STRAIN = 0.171
PRONE_POSE = {"pelvis_tilt": -1.5708, "pelvis_ty": 0.25}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


crawl = load_module("crawl", ROOT / "scripts/crawl.py")
render = load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")

# The force sets under test.  `stops` is what the plant does today; `ligament`
# and `joint_capsule` are the derived tissue classes.
ARMS = {
    "bare": dict(stops=False, classes=None),
    "stops": dict(stops=True, classes=None),
    "ligaments": dict(stops=False, classes=["ligament"]),
    "ligaments_capsules": dict(stops=False, classes=["ligament", "joint_capsule"]),
    "stops_ligaments": dict(stops=True, classes=["ligament"]),
    "stops_ligaments_capsules": dict(stops=True, classes=["ligament", "joint_capsule"]),
}


def open_stream(out, pose, arm):
    spec = ARMS[arm]
    return NativeMechanicalStream(
        ROOT, out, environment="upright", target_mass_kg=TARGET_MASS_KG,
        initial_pose=pose, augmented_registration=REGISTRATION,
        coordinate_limits=crawl.joint_stops() if spec["stops"] else None,
        tissue_ligaments=None if spec["classes"] is None else TISSUE,
        tissue_ligament_classes=spec["classes"])


def ligament_summary(state, spec_by_element):
    rows = state.get("tissue_ligaments")
    if not rows:
        return None
    strains = np.array([r["strain"] for r in rows])
    forces = np.array([r["total_force_n"] for r in rows])
    over = [dict(element=r["name"], structure=spec_by_element[r["name"]]["name"],
                 tissue_class=r["tissue_class"], body1=r["body1"], body2=r["body2"],
                 strain=r["strain"], force_n=r["total_force_n"])
            for r in rows if r["strain"] > ULTIMATE_STRAIN]
    over.sort(key=lambda r: -r["strain"])
    return dict(elements=len(rows),
                loaded=int((forces > 1e-9).sum()),
                total_tension_n=float(forces.sum()),
                median_strain=float(np.median(strains)),
                p90_strain=float(np.percentile(strains, 90)),
                maximum_strain=float(strains.max()),
                maximum_force_n=float(forces.max()),
                over_ultimate_strain=len(over),
                ultimate_strain=ULTIMATE_STRAIN,
                over_ultimate=over[:20],
                by_class=_by_class(rows))


def _by_class(rows):
    out = {}
    for row in rows:
        entry = out.setdefault(row["tissue_class"], dict(elements=0, loaded=0, tension_n=0.0))
        entry["elements"] += 1
        if row["total_force_n"] > 1e-9:
            entry["loaded"] += 1
            entry["tension_n"] += row["total_force_n"]
    return out


def statics_row(state, spec_by_element):
    total = np.zeros(3)
    for contact in state["contacts"]:
        total += np.asarray(contact["force_n"])
    return dict(vertical_contact_force_n=float(total[1]),
                weight_n=WEIGHT_N,
                weight_error_n=float(total[1] - WEIGHT_N),
                momentum_balance_residual_norm_n=float(
                    np.linalg.norm(state["momentum_balance_residual_n"])),
                pelvis_ty_m=state["coordinates"]["pelvis_ty"]["value"],
                worst_range_excursion_rad=crawl.worst_range_excursion(state)[0],
                worst_range_coordinate=crawl.worst_range_excursion(state)[1],
                tissue=ligament_summary(state, spec_by_element))


def run_statics(arm, pose, steps, dt, work, spec_by_element):
    out = work / ("statics-" + arm)
    shutil.rmtree(out, ignore_errors=True)
    stream = open_stream(out, pose, arm)
    try:
        initial = stream.snapshot()
        costs = []
        state = initial
        for _ in range(steps):
            mark = time.monotonic()
            state = stream.advance(dt)
            costs.append(time.monotonic() - mark)
        return dict(arm=arm, steps=steps, dt_s=dt,
                    initial=statics_row(initial, spec_by_element),
                    final=statics_row(state, spec_by_element),
                    cost=dict(advances=len(costs), median_s=statistics.median(costs),
                              worst_s=max(costs), best_s=min(costs),
                              total_wall_s=sum(costs)))
    finally:
        stream.close()


def run_excursion(arm, seconds, dt, work, spec_by_element, excitation=0.02):
    """Drop the body prone with nothing driving it, and see where it ends up.

    This is the configuration the withdrawn 973 mm crawl came from: the source
    model's `PassiveAnkleDamping` is literally `-0.1*qdot`, so an unloaded foot
    plantarflexes to 2.52 rad against a declared +/-0.873 and nothing objects.
    """
    out = work / ("excursion-" + arm)
    shutil.rmtree(out, ignore_errors=True)
    stream = open_stream(out, PRONE_POSE, arm)
    try:
        state = stream.snapshot()
        muscles = sorted(state["muscles"])
        targets = {k: 0.0 for k in
                   [s + "_" + t for s in crawl.ARM_PORT for t in "rl"] + list(crawl.LUMBAR_PORT)}
        worst, culprit, track, costs = 0.0, None, [], []
        for i in range(int(round(seconds / dt))):
            mark = time.monotonic()
            state = stream.advance(dt, actuation={m: excitation for m in muscles},
                                   coordinate_actuation=crawl.limb_pd(state, targets))
            costs.append(time.monotonic() - mark)
            excursion, name = crawl.worst_range_excursion(state)
            if excursion > worst:
                worst, culprit = excursion, name
            if i % 20 == 0:
                track.append(dict(time_s=state["time_s"],
                                  worst_excursion_rad=excursion, coordinate=name,
                                  pelvis_ty_m=state["coordinates"]["pelvis_ty"]["value"],
                                  kinetic_energy_j=state["kinetic_energy_j"]))
        return dict(arm=arm, seconds=seconds, dt_s=dt,
                    worst_range_excursion_rad=worst,
                    worst_range_excursion_deg=math.degrees(worst),
                    worst_range_coordinate=culprit,
                    final=statics_row(state, spec_by_element),
                    coordinates_outside_range=_outside(state),
                    track=track,
                    cost=dict(advances=len(costs), median_s=statistics.median(costs),
                              worst_s=max(costs), total_wall_s=sum(costs)))
    finally:
        stream.close()


def _outside(state):
    ranges = crawl.declared_ranges()
    out = []
    for name, (low, high) in sorted(ranges.items()):
        value = state["coordinates"][name]["value"]
        past = max(low - value, value - high, 0.0)
        if past > 0:
            out.append(dict(coordinate=name, value_rad=value, low_rad=low, high_rad=high,
                            past_rad=past, past_deg=math.degrees(past)))
    out.sort(key=lambda r: -r["past_rad"])
    return out


def cross_implementation_gate(state, elements):
    """The engine's ligament length against this script's own kinematics.

    Two independent implementations of the same geometry: the engine builds a
    GeometryPath between two station points on two OpenSim Bodies; this walks the
    model XML's own joint definitions.  They must agree to machine precision.
    """
    model = render.OsimModel(BUNDLE / "model.osim")
    transforms = model.forward({k: v["value"] for k, v in state["coordinates"].items()})
    spec = {r["element"]: r for r in elements}
    errors = []
    for row in state.get("tissue_ligaments", []):
        record = spec[row["name"]]
        a = transforms[record["body1"]][:3, :3] @ np.asarray(record["point1_m"]) \
            + transforms[record["body1"]][:3, 3]
        b = transforms[record["body2"]][:3, :3] @ np.asarray(record["point2_m"]) \
            + transforms[record["body2"]][:3, 3]
        errors.append(abs(float(np.linalg.norm(a - b)) - row["length_m"]))
    return dict(elements=len(errors),
                maximum_absolute_error_m=max(errors) if errors else None,
                mean_absolute_error_m=float(np.mean(errors)) if errors else None)


def reference_pose_gate(elements, binding):
    """At the binding's reference pose every element must be exactly slack."""
    model = render.OsimModel(BUNDLE / "model.osim")
    transforms = model.forward(binding["reference_pose_rad"])
    errors = []
    for record in elements:
        a = transforms[record["body1"]][:3, :3] @ np.asarray(record["point1_m"]) \
            + transforms[record["body1"]][:3, 3]
        b = transforms[record["body2"]][:3, :3] @ np.asarray(record["point2_m"]) \
            + transforms[record["body2"]][:3, 3]
        errors.append(abs(float(np.linalg.norm(a - b)) - record["slack_length_m"]))
    return dict(elements=len(errors), maximum_absolute_error_m=max(errors),
                basis="slack length is BY CONSTRUCTION the reference-pose separation, so this "
                      "gate tests the frame conventions, not the physics")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--steps", type=int, default=25)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--excursion-seconds", type=float, default=2.0)
    parser.add_argument("--arm", action="append", default=[])
    parser.add_argument("--skip-excursion", action="store_true")
    parser.add_argument("--work", default="data/derived/tissue-mechanics")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    elements = json.loads((ROOT / TISSUE / "ligaments.json").read_text())["elements"]
    spec_by_element = {r["element"]: r for r in elements}
    binding = json.loads((ROOT / "data/derived/anatomy-segment-binding/binding.json").read_text())
    pose = json.loads((BUNDLE / "initial_pose.json").read_text())
    work = ROOT / args.work
    work.mkdir(parents=True, exist_ok=True)
    work = work / ("run-" + str(int(time.time())))
    work.mkdir()

    selected = args.arm or list(ARMS)
    statics, excursion = [], []
    cross = None
    for arm in selected:
        print("== statics " + arm, flush=True)
        try:
            row = run_statics(arm, pose, args.steps, args.dt, work, spec_by_element)
        except Exception as error:
            row = dict(arm=arm, failed=repr(error))
        statics.append(row)
        print(json.dumps({k: v for k, v in row.items() if k in ("cost", "failed")}, indent=2),
              flush=True)
    if not args.skip_excursion:
        for arm in selected:
            print("== excursion " + arm, flush=True)
            try:
                row = run_excursion(arm, args.excursion_seconds, args.dt, work, spec_by_element)
            except Exception as error:
                row = dict(arm=arm, failed=repr(error))
            excursion.append(row)
            print(json.dumps({k: v for k, v in row.items()
                              if k in ("worst_range_excursion_deg", "worst_range_coordinate",
                                       "cost", "failed")}, indent=2), flush=True)

    # The cross-implementation gate needs one state that carries ligaments.
    gate_out = work / "gate"
    shutil.rmtree(gate_out, ignore_errors=True)
    stream = open_stream(gate_out, pose, "ligaments_capsules")
    try:
        cross = cross_implementation_gate(stream.snapshot(), elements)
    finally:
        stream.close()

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(
        schema="ihm.tissue-mechanics-measurement.v1",
        pose="data/models/engineering_stance_v1/initial_pose.json",
        prone_pose=PRONE_POSE, registration=REGISTRATION, tissue_bundle=TISSUE,
        target_mass_kg=TARGET_MASS_KG, weight_n=WEIGHT_N,
        steps=args.steps, dt_s=args.dt, excursion_seconds=args.excursion_seconds,
        work=str(work.relative_to(ROOT)),
        joint_stop=crawl.JOINT_STOP,
        gates=dict(cross_implementation_length=cross,
                   reference_pose_slack=reference_pose_gate(elements, binding)),
        statics=statics, excursion=excursion), indent=2) + "\n")
    print(str(out))


if __name__ == "__main__":
    main()
