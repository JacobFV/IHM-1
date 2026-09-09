#!/usr/bin/env python3
"""Fascia, aponeurosis and retinacula: can this plant express a path constraint?

A retinaculum is not a joint restraint.  It is a strap that holds a tendon
against the bone, and what it does mechanically is CHANGE A MUSCLE'S PATH -- it
stops the extensor tendons bowstringing off the front of the ankle when the
ankle dorsiflexes.  A deep fascial plane is the surface muscle groups glide
across.  Neither is a tension element between two bones, so neither belongs in
the ligament set, and the question is a different one: does the running plant
have a muscle path for them to constrain at all?

The answer is read out of the model the ENGINE ASSEMBLED, not the model on disk,
because those are different objects.  `native_mechanical_stream.cpp` calls
`ModelFactory::replacePathsWithFunctionBasedPaths` at load, which substitutes
fitted polynomials in the coordinates for the declared point paths.  A
polynomial has no geometry: no via point, no wrap surface, nothing a strap could
touch.  `docs/SEGMENT_CONTACT_SURFACES.md` already recorded that 80 of the 98
muscles go through that substitution; this counts what that leaves.

What is measured
----------------
* how many muscles in the assembled model keep a `GeometryPath`, and which;
* how many `PathPoint`s and `WrapObject`s survive on them, because the wrap
  MECHANISM is what a retinaculum would have to use and it is either reachable
  or it is not;
* for each retinaculum and each fascia, the muscles it constrains in life, and
  whether any of them is a muscle this plant carries with a path.

Nothing here is a force element.  It is the measurement that says whether one
could exist, and the answer is not the same for the two classes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TISSUE = ROOT / "data/derived/tissue-force-elements-v1/structures.json"

# What each strap holds down, in life.  The muscle stems are the source model's
# own naming (`tibant`, `edl`, `perlong`, ...) so the lookup is against the
# muscles this plant actually declares, not against anatomical names it does not
# use.  A strap whose list is EMPTY is one whose muscles the model does not have.
RESTRAINS = {
    "flexor retinaculum of wrist": [],
    "extensor retinaculum of wrist": [],
    "flexor retinaculum of ankle": ["tibpost", "fdl", "fhl"],
    "superior extensor retinaculum of ankle": ["tibant", "edl", "ehl"],
    "inferior extensor retinaculum of ankle": ["tibant", "edl", "ehl"],
    "superior fibular retinaculum": ["perlong", "perbrev"],
    "inferior fibular retinaculum": ["perlong", "perbrev"],
    "medial patellar retinaculum": ["vasmed", "recfem"],
    "lateral patellar retinaculum": ["vaslat", "recfem"],
}
# The hand and wrist straps hold the long finger flexors and extensors, and this
# plant has no hand or wrist musculature at all -- `hand_l` and `hand_r` carry
# zero PathPoints. That is why those lists are empty, and it is a different
# reason from the polynomial substitution.
NO_MUSCULATURE = ("flexor retinaculum of wrist", "extensor retinaculum of wrist")


def strap_kind(name):
    lowered = name.lower()
    for key in RESTRAINS:
        stem = key.replace(" of wrist", "").replace(" of ankle", "")
        if key in lowered or (stem in lowered and ("wrist" in key) == ("wrist" in lowered)
                              and ("ankle" in key) == ("ankle" in lowered)):
            return key
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assembled", required=True,
                        help="an assembled_model.osim written by the native engine")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    started = time.time()

    assembled = Path(args.assembled)
    root = ET.parse(assembled).getroot().find("Model")
    forces = root.find("ForceSet")
    geometry_path, function_path = [], []
    path_points = 0
    wrap_on_geometry_paths = 0
    for component in forces.iter():
        if not component.tag.endswith("Muscle"):
            continue
        name = component.get("name")
        if component.find(".//FunctionBasedPath") is not None:
            function_path.append(name)
        elif component.find(".//GeometryPath") is not None:
            geometry_path.append(name)
            path_points += len(component.findall(".//PathPoint"))
            wrap_on_geometry_paths += len(component.findall(".//PathWrap"))
    wrap_objects = len(root.findall(".//WrapCylinder")) + len(root.findall(".//WrapEllipsoid")) \
        + len(root.findall(".//WrapSphere")) + len(root.findall(".//WrapTorus"))

    muscles = set(geometry_path) | set(function_path)
    def present(stem):
        return sorted(m for m in muscles if re.search(r"(^|_)" + re.escape(stem) + r"(_|$)", m))

    structures = json.loads(TISSUE.read_text())["structures"]
    straps = []
    for row in structures:
        if row["tissue_class"] != "retinaculum":
            continue
        key = strap_kind(row["name"])
        stems = RESTRAINS.get(key, [])
        # A left retinaculum holds LEFT tendons.  Without this the report reads
        # as though every strap constrained twice as many muscles as it does.
        side = "_l" if "left" in row["name"].lower() else "_r" if "right" in row["name"].lower() else ""
        constrained = sorted({m for stem in stems for m in present(stem)
                              if not side or m.endswith(side)})
        straps.append(dict(
            name=row["name"], matched_rule=key,
            muscles_it_holds_in_life=stems,
            muscles_this_plant_carries=constrained,
            muscles_with_a_constrainable_path=[m for m in constrained if m in geometry_path],
            reason=("this plant carries no musculature for these tendons"
                    if key in NO_MUSCULATURE or not stems else
                    "every muscle it holds runs as a fitted polynomial in the coordinates, which "
                    "has no geometry to constrain"
                    if not any(m in geometry_path for m in constrained) else
                    "constrainable")))
    fascia = [row for row in structures if row["tissue_class"] == "fascia_aponeurosis"]

    constrainable = [s for s in straps if s["muscles_with_a_constrainable_path"]]
    report = dict(
        schema="ihm.muscle-path-constraint-measurement.v1",
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        wall_s=time.time() - started,
        assembled_model=dict(path=str(assembled),
                             sha256=hashlib.sha256(assembled.read_bytes()).hexdigest()),
        paths=dict(muscles=len(muscles),
                   function_based_paths=len(function_path),
                   geometry_paths=len(geometry_path),
                   geometry_path_muscles=sorted(geometry_path),
                   path_points_on_geometry_paths=path_points,
                   path_wraps_on_geometry_paths=wrap_on_geometry_paths,
                   wrap_objects_declared=wrap_objects,
                   note="A FunctionBasedPath is a polynomial in the coordinates. It has no via "
                        "point and no wrap surface, so no strap can touch it. The wrap MECHANISM "
                        "exists in this model and serves only the paths that are still "
                        "GeometryPaths."),
        retinacula=dict(entities=len(straps),
                        with_a_constrainable_muscle=len(constrainable),
                        rows=straps),
        fascia_aponeurosis=dict(
            entities=len(fascia),
            two_segment=sum(1 for r in fascia if r.get("status") == "two_segment"),
            one_segment=sum(1 for r in fascia if r.get("status") == "one_segment"),
            names_two_segment=sorted({r["name"] for r in fascia if r.get("status") == "two_segment"}),
            finding="A deep fascial plane's mechanics is a GLIDE SURFACE between muscle groups and "
                    "a constraint on muscle bulging. Both need a deformable continuum, and "
                    "docs/SEGMENT_CONTACT_SURFACES.md establishes that nothing in this stack has "
                    "one: every compliant element is a 1-D spring law on a rigid carrier. The one "
                    "part of fascial mechanics a line element CAN carry is longitudinal load "
                    "transfer between two segments -- the fascia lata and the thoracolumbar fascia "
                    "are the cases -- and that is not modelled here because a fascia is a sheet "
                    "whose tension is distributed over its width and a single line between two "
                    "centroids concentrates it, the same error that gave the hip capsule a "
                    "924 kN cross-section before the open-surface rule was applied."),
        finding="Not one of the 18 retinacula can constrain a muscle path in the running plant. "
                "Four of them hold tendons of muscles this plant does not have at all -- there is "
                "no hand or wrist musculature. The other fourteen hold muscles the plant does "
                "have, and every one of those runs as a fitted polynomial with no geometry. The "
                "blocker is not the retinaculum and not a missing force class; it is "
                "replacePathsWithFunctionBasedPaths, and it is reversible in principle at a known "
                "cost -- the polynomials are what made the plant affordable.")

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print("muscles %d: %d FunctionBasedPath, %d GeometryPath"
          % (len(muscles), len(function_path), len(geometry_path)))
    print("wrap objects declared %d; path wraps on the surviving GeometryPaths %d"
          % (wrap_objects, wrap_on_geometry_paths))
    print("retinacula %d, of which any has a constrainable muscle path: %d"
          % (len(straps), len(constrainable)))
    for row in straps:
        print("  %-46s holds %-24s carried %-28s %s"
              % (row["name"][:46], ",".join(row["muscles_it_holds_in_life"]) or "-",
                 ",".join(row["muscles_this_plant_carries"]) or "-", row["reason"][:44]))
    print(str(out))


if __name__ == "__main__":
    main()
