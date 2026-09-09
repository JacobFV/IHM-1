#!/usr/bin/env python3
"""Turn the body's TISSUE SURFACES into force elements, or say why they cannot be.

625 of the 4,000 bound entities are tissue structures -- ligaments, capsules,
menisci, discs, cartilage, bursae, fascia, retinacula, adipose -- and every one
of them has real mesh geometry.  The plant carries force elements for NONE of
them.  They exist as geometry and not as mechanics.

The recorded blocker was structural: `data/derived/anatomy-segment-binding/`
assigns each whole entity to exactly ONE segment by nearest-bone-group vertex
vote, and a ligament by definition spans TWO bones, so the binding artifact as
it stands cannot supply a two-ended attachment.  It cannot, but the vote it is
built from can: the SAME per-vertex nearest-bone-group distances that produce
the single winner also produce the runner-up's vertices, and the centroid of
each side's own vertices is an attachment site on each of two bodies.

That is what this script derives.  It is a CONSTRUCTION from mesh geometry and
a published modulus, not an anatomical measurement of insertion sites, and the
gates below are what stop it being taken for one.

What comes out
--------------
`ligaments.json`    every two-segment structure, with its two attachment points
                    in the two segments' own frames, its slack length, its
                    derived cross-section and its linear stiffness.
`ligaments.txt`     the same rows in the native engine's spec format, consumed
                    by `scripts/native_mechanical_stream.cpp` as
                    `Blankevoort1991Ligament` force elements.
`inventory.json`    every tissue class, how many entities it holds, how many are
                    two-segment on THIS scaffold, and -- the number that
                    matters -- how many are one-segment because the scaffold has
                    no joint where the structure acts.

How an attachment is derived
----------------------------
1. The 22 segments' anatomical bone groups are named from the atlas' own labels
   by `scripts/bind_anatomy_to_segments.py`'s own rule, imported, not restated.
2. Every vertex of the structure's surface takes the segment whose bone group is
   nearest.  The two segments with the most votes are its two ends.
3. Each end's attachment point is the centroid of the quarter of that end's
   vertices FURTHEST from the other end's centroid -- the tip, not the middle.
   A centroid-of-half attachment measures 20.9 mm on an ACL whose published
   length is 32; the tip quantile measures 33.4.
4. The points are carried into each segment's own frame through the binding's
   OWN recorded similarity and reference pose.  No new registration is fitted,
   so this inherits the binding's 24.7 mm RMS residual and adds nothing to it.
5. The slack length is the straight-line separation at that reference pose, so
   every ligament is exactly slack in the pose its geometry was registered in
   and its strain anywhere else is a real kinematic strain.
6. The linear stiffness is E*A: `E` is the body's OWN declared ligament
   along-fibre modulus (332.2 MPa, Quapp and Weiss 1998, human MCL, carried in
   `data/derived/tissue-material-candidate-v1/materials.json`), and `A` is the
   structure's own mesh volume divided by its own derived length.

Gates -- against cases whose answer is known
--------------------------------------------
* **Named pairs.** A ligament everybody can name has two bones everybody knows.
  ACL and PCL must come out femur+tibia, the collaterals femur+tibia, the
  iliofemoral/pubofemoral/ischiofemoral pelvis+femur, the talofibular
  tibia+talus, the annular ligament radius+ulna.  And the negative control:
  sacrotuberous, sacrospinous and inguinal must come out ONE segment, because
  the sacrum and the ischium are both inside the scaffold's single `pelvis`
  body -- a two-segment answer there would be the derivation inventing a joint.
* **Published lengths.** Six ligaments with a published straight-line length.
  Reported as a ratio, not asserted as agreement.
* **Published stiffness.** The derived ACL linear stiffness against the value
  the Blankevoort knee literature uses. The ratio is REPORTED; it is not tuned.
* **Convexity control on the vote.** A tooth is inside one bone, so it must
  return a second-segment share of zero. Run and recorded.

Usage
-----
    python scripts/build_tissue_force_elements.py --out data/derived/tissue-force-elements-v1
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ANATOMY = ROOT / "data/derived/canonical/anatomy.json"
BINDING = ROOT / "data/derived/anatomy-segment-binding/binding.json"
MODEL = ROOT / "data/models/engineering_stance_v1/model.osim"
MATERIALS = ROOT / "data/derived/tissue-material-candidate-v1/materials.json"
OUT = ROOT / "data/derived/tissue-force-elements-v1"

# First match wins.  The order is the one `scripts/audit_joint_substrate.py`
# already uses, so the two lanes bucket the same name the same way: the knee
# soft tissue named after its meniscus leaves the plain ligament bucket first.
CLASSES = (
    ("bursa_synovium", r"bursa|synovi"),
    ("retinaculum", r"retinacul"),
    ("meniscus_disc", r"meniscus|meniscotibial|meniscopatellar|intervertebral disc"
                      r"|articular disc|nucleus pulposus|interpubic disc|symphysis"),
    ("joint_capsule", r"articular capsule|joint capsule|zona orbicularis|frenula capsulae"),
    ("fascia_aponeurosis", r"fascia|aponeuros"),
    ("cartilage", r"cartilage"),
    ("tendon", r"tendon"),
    ("ligament", r"ligament"),
    ("adipose", r"fat pad|adipose"),
    ("skin", r"\bskin\b"),
)
# Classes whose mechanics is a tension element between two bones.  Everything
# else is assessed, not instantiated.
TENSILE = ("ligament", "joint_capsule")

# The end quantile, and the floor under how many vertices an end must have for
# its centroid to mean anything.
TIP_QUANTILE = 25.0
MIN_END_VERTICES = 4
# A second segment claiming less than this share of the vertices is not an end.
# The count is flat from 1% to 10% (112 / 112 / 106 / 103 ligaments) and only
# starts falling at 20%, so the operating point is not a cliff edge; the whole
# curve is written into the report.
MIN_SECOND_SHARE = 0.05
SHARE_SWEEP = (0.01, 0.02, 0.05, 0.10, 0.20, 0.30)

# Blankevoort's own transition strain, and the damping convention the knee
# literature uses with it: c = 0.003 * k, in N.s/strain.
TRANSITION_STRAIN = 0.06
NORMALIZED_DAMPING_S = 0.003

# Cases with an answer anatomy fixes and this script never sees.  `expect` is
# the pair a correct derivation must name; a one-element set means the structure
# must NOT come out two-segment on this scaffold.
NAMED_PAIRS = [
    ("anterior cruciate ligament", {"femur", "tibia"}),
    ("posterior cruciate ligament", {"femur", "tibia"}),
    ("superficial part of tibial collateral ligament", {"femur", "tibia"}),
    ("deep part of tibial collateral ligament", {"femur", "tibia"}),
    ("fibular collateral ligament", {"femur", "tibia"}),
    ("descending part of iliofemoral ligament", {"pelvis", "femur"}),
    ("pubofemoral ligament", {"pelvis", "femur"}),
    ("ischiofemoral ligament", {"pelvis", "femur"}),
    ("anterior talofibular ligament", {"tibia", "talus"}),
    ("posterior talofibular ligament", {"tibia", "talus"}),
    ("calcaneofibular ligament", {"tibia", "calcn"}),
    ("annular ligament of radius", {"radius", "ulna"}),
    # Negative controls: both bones are inside ONE scaffold segment.
    ("sacrotuberous ligament", {"pelvis"}),
    ("sacrospinous ligament", {"pelvis"}),
    ("inguinal ligament", {"pelvis"}),
]

# Straight-line ligament lengths from the anatomical literature, in mm.  These
# are population means for a different specimen; they are a magnitude check on
# the derivation, not a target.
PUBLISHED_LENGTH_MM = {
    "anterior cruciate ligament": (32.0, "Odensten and Gillquist 1985: ACL mean length 31 mm; "
                                         "Duthon 2006 gives 32 mm"),
    "posterior cruciate ligament": (38.0, "Girgis 1975 / Amis 2006: PCL mean length 38 mm"),
    "superficial part of tibial collateral ligament": (
        95.0, "LaPrade 2007: superficial MCL femoral attachment to distal tibial attachment "
              "about 95 mm"),
    "fibular collateral ligament": (60.0, "LaPrade 2003: FCL mean length 59.9 mm"),
    "anterior talofibular ligament": (20.0, "Golano 2010: ATFL 15-20 mm"),
    "calcaneofibular ligament": (25.0, "Golano 2010: CFL 20-25 mm"),
}
# Published Blankevoort-form linear stiffness, N per unit strain.
PUBLISHED_STIFFNESS_N = {
    "anterior cruciate ligament": (5000.0, "Blankevoort and Huiskes 1991 / Smith 2016: ACL "
                                           "modelled as two bundles of 2500 N each"),
    "posterior cruciate ligament": (9000.0, "Smith 2016: PCL bundles 3500 + 5500 N"),
    "fibular collateral ligament": (2000.0, "Smith 2016: LCL 2000 N"),
}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git_sha():
    sha = os.environ.get("IHM_GIT_SHA") or os.environ.get("IBM_GIT_SHA")
    if sha:
        return sha
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "git-unknown"


def tissue_class(name):
    lowered = name.lower()
    for label, pattern in CLASSES:
        if re.search(pattern, lowered):
            return label
    return None


def side_stripped(segment):
    return segment[:-2] if segment.endswith(("_l", "_r")) else segment


def mesh_volume(vertices, indices):
    """Signed volume and area of the surface, after a lossless weld."""
    import trimesh
    mesh = trimesh.Trimesh(vertices=vertices, faces=indices, process=True)
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.update_faces(mesh.unique_faces())
    mesh.remove_unreferenced_vertices()
    return float(abs(mesh.volume)), float(mesh.area), bool(mesh.is_watertight)


def tissue_volume(name, role, signed_m3, area_m2, watertight):
    """What the surface holds as TISSUE, by the repository's own rule.

    An OPEN surface's signed integral is not a volume.  `tissue_materials` owns
    that argument already -- 84 open connective sheets were carrying 29.6 L of
    the muscle they wrap -- and it is the reason this function exists rather
    than `abs(mesh.volume)`.  It matters here more than it did there, because a
    cross-section is a volume divided by a length: the left hip capsule is an
    OPEN sleeve whose signed integral is the femoral head and neck inside it,
    and taking that as capsule tissue gives a 2,781 mm2 cross-section and a
    924 kN ligament.

    The rule is `tissue_materials.material_volume`, imported, with one addition:
    a surface that is open and is NOT on that module's membrane token list still
    cannot have its signed integral believed, so the same open-sheet bound is
    applied to it and the rule is recorded as `open_surface_not_a_membrane`.
    """
    from ihm.assembly.tissue_materials import material_volume, THICKNESS
    volume, record = material_volume(name, role, signed_m3, area_m2, watertight)
    if not watertight and record["rule"] == "solid":
        sheet = area_m2 * THICKNESS["membrane"][0]
        volume = min(signed_m3, sheet) if signed_m3 > 0 else sheet
        record = dict(rule="open_surface_not_a_membrane",
                      thickness_m=THICKNESS["membrane"][0],
                      thickness_tier=THICKNESS["membrane"][1],
                      open_surface_signed_integral_m3=signed_m3,
                      basis="the surface is open, so its signed integral is not the tissue it "
                            "holds; the same open-sheet bound tissue_materials applies to "
                            "membranes is applied here")
    return volume, record


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    started = time.time()

    bind = load_module("bind_anatomy", ROOT / "scripts/bind_anatomy_to_segments.py")
    render = load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")

    anatomy = json.loads(ANATOMY.read_text())
    entities = anatomy["entities"]
    geometry_path = {e["id"]: ROOT / e["reference_geometry"]["path"]
                     for e in entities if e.get("reference_geometry")}
    bind.GEOMETRY_PATH.update(geometry_path)

    # ---------------------------------------------------------------- segments
    groups = {}
    for entity in entities:
        if entity["role"] != "rigid_bone":
            continue
        segment = bind.named_segment(entity["name"])
        if segment is None:
            raise ValueError("unnamed bone: " + entity["name"])
        groups.setdefault(segment, []).append(entity["id"])
    segments = sorted(groups)
    clouds = {s: np.concatenate([bind.read_geometry(i) for i in groups[s]]) for s in segments}
    trees = [cKDTree(clouds[s]) for s in segments]
    print(f"segments {len(segments)}; bone vertices "
          f"{sum(len(v) for v in clouds.values())}", flush=True)

    # ------------------------------------------------ the binding's own frames
    binding = json.loads(BINDING.read_text())
    similarity = np.asarray(binding["similarity_atlas_from_opensim_ground"], float)
    inverse = np.linalg.inv(similarity)
    scale = float(binding["registration"]["scale"])
    reference_pose = binding["reference_pose_rad"]
    model = render.OsimModel(MODEL)
    rest = model.forward(reference_pose)
    parent = {j["child"]: j["parent"] for j in model.joints}

    def to_local(point_atlas, segment):
        world = inverse[:3, :3] @ point_atlas + inverse[:3, 3]
        transform = rest[segment]
        return np.linalg.inv(transform[:3, :3]) @ (world - transform[:3, 3])

    def to_ground(point_local, segment):
        transform = rest[segment]
        return transform[:3, :3] @ point_local + transform[:3, 3]

    def joints_between(a, b):
        """How many scaffold joints separate two segments (0 = same body)."""
        def chain(x):
            out = [x]
            while x in parent:
                x = parent[x]
                out.append(x)
            return out
        ca, cb = chain(a), chain(b)
        common = next((x for x in ca if x in cb), None)
        if common is None:
            return None
        return ca.index(common) + cb.index(common)

    # ------------------------------------------------------------- the modulus
    materials = json.loads(MATERIALS.read_text())["materials"]["ligament"]
    modulus = materials["linear"]["young_modulus_along_fibre"]
    if modulus["unit"] != "Pa":
        raise ValueError("declared ligament modulus is not in Pa")
    young_pa = float(modulus["value"])

    # ------------------------------------------------------------- the derivation
    rows = []
    counts = {}
    print("voting...", flush=True)
    for n, entity in enumerate(entities):
        label = tissue_class(entity["name"])
        if label is None:
            continue
        counts[label] = counts.get(label, 0) + 1
        if not entity.get("reference_geometry"):
            rows.append(dict(id=entity["id"], name=entity["name"], tissue_class=label,
                             role=entity["role"], status="no_geometry"))
            continue
        payload = json.loads(gzip.decompress(geometry_path[entity["id"]].read_bytes()))
        vertices = np.asarray(payload["positions"], float).reshape(-1, 3)
        indices = np.asarray(payload["indices"], int).reshape(-1, 3)
        distances = np.stack([tree.query(vertices)[0] for tree in trees], axis=1)
        nearest = np.argmin(distances, axis=1)
        votes = np.bincount(nearest, minlength=len(segments))
        order = np.argsort(-votes)
        first, second = int(order[0]), int(order[1])
        share_first = float(votes[first]) / len(vertices)
        share_second = float(votes[second]) / len(vertices)
        signed, area, watertight = mesh_volume(vertices, indices)
        volume, volume_rule = tissue_volume(entity["name"], entity["role"], signed, area, watertight)
        row = dict(id=entity["id"], name=entity["name"], tissue_class=label,
                   role=entity["role"], vertices=int(len(vertices)),
                   segment=segments[first], segment_share=share_first,
                   runner_up=segments[second], runner_up_share=share_second,
                   atlas_signed_volume_m3=signed, atlas_tissue_volume_m3=volume,
                   atlas_surface_area_m2=area, watertight=watertight,
                   volume_rule=volume_rule)
        ends = [k for k in (first, second)
                if votes[k] >= MIN_END_VERTICES and float(votes[k]) / len(vertices) >= MIN_SECOND_SHARE]
        if len(ends) < 2:
            row["status"] = "one_segment"
            row["reason"] = ("the whole surface votes for one scaffold segment: the two bones it "
                             "joins are the same rigid body here, so no relative motion exists "
                             "for it to resist")
            rows.append(row)
            continue
        a, b = ends
        va, vb = vertices[nearest == a], vertices[nearest == b]
        centre_a, centre_b = va.mean(0), vb.mean(0)

        def tip(side, other_centre):
            radial = np.linalg.norm(side - other_centre, axis=1)
            threshold = np.percentile(radial, 100.0 - TIP_QUANTILE)
            chosen = side[radial >= threshold]
            return (chosen if len(chosen) else side).mean(0)

        point_a = to_local(tip(va, centre_b), segments[a])
        point_b = to_local(tip(vb, centre_a), segments[b])
        length = float(np.linalg.norm(to_ground(point_a, segments[a])
                                      - to_ground(point_b, segments[b])))
        crossed = joints_between(segments[a], segments[b])
        # Atlas metric -> OpenSim metric.  The similarity carries a 0.963 scale,
        # so a volume read in the atlas is scale^3 of the same volume on the
        # body the engine integrates.
        volume_model = volume / scale ** 3
        section = volume_model / length if length > 0 else float("nan")
        stiffness = young_pa * section
        row.update(status="two_segment", body1=segments[a], body2=segments[b],
                   point1_m=[float(v) for v in point_a], point2_m=[float(v) for v in point_b],
                   slack_length_m=length, model_volume_m3=volume_model,
                   cross_section_m2=section, linear_stiffness_n=stiffness,
                   damping_n_s_per_strain=NORMALIZED_DAMPING_S * stiffness,
                   transition_strain=TRANSITION_STRAIN,
                   scaffold_joints_crossed=crossed,
                   end_vertices=[int(votes[a]), int(votes[b])])
        rows.append(row)
        if n % 500 == 0:
            print(f"  {n}/{len(entities)}", flush=True)

    # ------------------------------------------------------------------- gates
    by_name = {r["name"]: r for r in rows}
    named = []
    for stem, expect in NAMED_PAIRS:
        for side in ("left", "right"):
            key = side + " " + stem
            row = by_name.get(key)
            if row is None:
                continue
            if row.get("status") == "two_segment":
                got = {side_stripped(row["body1"]), side_stripped(row["body2"])}
            else:
                got = {side_stripped(row["segment"])}
            named.append(dict(name=key, expected=sorted(expect), derived=sorted(got),
                              status=row.get("status"), passed=got == expect,
                              runner_up_share=row.get("runner_up_share")))
    lengths = []
    for stem, (published, source) in PUBLISHED_LENGTH_MM.items():
        for side in ("left", "right"):
            row = by_name.get(side + " " + stem)
            if row is None or row.get("status") != "two_segment":
                continue
            lengths.append(dict(name=side + " " + stem, derived_mm=row["slack_length_m"] * 1000,
                                published_mm=published, ratio=row["slack_length_m"] * 1000 / published,
                                source=source))
    stiffnesses = []
    for stem, (published, source) in PUBLISHED_STIFFNESS_N.items():
        for side in ("left", "right"):
            row = by_name.get(side + " " + stem)
            if row is None or row.get("status") != "two_segment":
                continue
            stiffnesses.append(dict(name=side + " " + stem, derived_n=row["linear_stiffness_n"],
                                    published_n=published,
                                    ratio=row["linear_stiffness_n"] / published, source=source))

    # Control: a structure wholly inside one bone must give a second-segment
    # share of zero.  Teeth are the cleanest case in the atlas.
    control = []
    for entity in entities:
        if "molar tooth" not in entity["name"] or not entity.get("reference_geometry"):
            continue
        payload = json.loads(gzip.decompress(geometry_path[entity["id"]].read_bytes()))
        vertices = np.asarray(payload["positions"], float).reshape(-1, 3)
        distances = np.stack([tree.query(vertices)[0] for tree in trees], axis=1)
        votes = np.bincount(np.argmin(distances, axis=1), minlength=len(segments))
        order = np.argsort(-votes)
        control.append(dict(name=entity["name"], segment=segments[int(order[0])],
                            runner_up_share=float(votes[int(order[1])]) / len(vertices)))
        if len(control) >= 6:
            break

    # ---------------------------------------------------------------- inventory
    inventory = {}
    for label, _ in CLASSES:
        members = [r for r in rows if r["tissue_class"] == label]
        two = [r for r in members if r.get("status") == "two_segment"]
        inventory[label] = dict(
            entities=len(members),
            two_segment=len(two),
            one_segment=sum(1 for r in members if r.get("status") == "one_segment"),
            no_geometry=sum(1 for r in members if r.get("status") == "no_geometry"),
            instantiated_as_force_element=len(two) if label in TENSILE else 0,
            one_segment_by_segment=_tally(r["segment"] for r in members
                                          if r.get("status") == "one_segment"),
            two_segment_pairs=_tally("+".join(sorted((r["body1"], r["body2"]))) for r in two),
            joints_crossed=_tally(str(r.get("scaffold_joints_crossed")) for r in two))
    sweep = {}
    for threshold in SHARE_SWEEP:
        sweep[str(threshold)] = sum(
            1 for r in rows if r["tissue_class"] == "ligament"
            and r.get("runner_up_share", 0.0) >= threshold)

    # ----------------------------------------------------------------- outputs
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    elements = [r for r in rows if r.get("status") == "two_segment" and r["tissue_class"] in TENSILE]
    elements.sort(key=lambda r: r["id"])
    for index, row in enumerate(elements):
        row["element"] = "lig_%03d_%s" % (index, re.sub(r"[^A-Za-z0-9]", "", row["id"])[-12:])
    spec = ["IHM_TISSUE_LIGAMENTS_V1 %d" % len(elements)]
    for row in elements:
        spec.append(" ".join(map(str, [
            row["element"], row["body1"], *row["point1_m"], row["body2"], *row["point2_m"],
            row["linear_stiffness_n"], row["slack_length_m"], row["transition_strain"],
            row["damping_n_s_per_strain"]])))
    (out / "ligaments.txt").write_text("\n".join(spec) + "\n")
    (out / "ligaments.json").write_text(json.dumps(dict(
        schema="ihm.tissue-force-elements.v1", elements=elements), indent=2) + "\n")
    (out / "structures.json").write_text(json.dumps(dict(
        schema="ihm.tissue-force-elements.v1", structures=rows), indent=2) + "\n")

    manifest = dict(
        schema="ihm.tissue-force-elements.v1",
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        wall_s=time.time() - started,
        provenance=dict(git_sha=git_sha(), script="scripts/build_tissue_force_elements.py",
                        script_sha256=sha256(__file__),
                        anatomy={"path": str(ANATOMY.relative_to(ROOT)), "sha256": sha256(ANATOMY)},
                        binding={"path": str(BINDING.relative_to(ROOT)), "sha256": sha256(BINDING)},
                        model={"path": str(MODEL.relative_to(ROOT)), "sha256": sha256(MODEL)},
                        materials={"path": str(MATERIALS.relative_to(ROOT)),
                                   "sha256": sha256(MATERIALS)}),
        registration=dict(
            basis="the binding's OWN recorded similarity and reference pose; no new fit",
            scale=scale, reference_pose_rad=reference_pose,
            inherited_residual_note="24.7 mm RMS on bone-group centroids, from "
                                    "data/derived/anatomy-segment-binding/report.json. This "
                                    "derivation inherits it and adds no fit of its own."),
        material=dict(young_modulus_along_fibre_pa=young_pa,
                      source=modulus["source"], tier=modulus["tier"], note=modulus["note"],
                      stiffness_rule="linear_stiffness (N per unit strain) = E * A, with A the "
                                     "structure's own welded mesh volume divided by its own "
                                     "derived slack length, both in the model's metric",
                      transition_strain=TRANSITION_STRAIN,
                      damping_rule="c = %g * k, N.s per unit strain, the convention the "
                                   "Blankevoort knee literature uses" % NORMALIZED_DAMPING_S),
        derivation=dict(tip_quantile_percent=TIP_QUANTILE,
                        minimum_end_vertices=MIN_END_VERTICES,
                        minimum_second_segment_share=MIN_SECOND_SHARE,
                        second_segment_share_sweep_ligaments=sweep),
        inventory=inventory,
        totals=dict(tissue_entities=sum(counts.values()),
                    classified=counts,
                    force_elements_written=len(elements)),
        gates=dict(named_pairs=named,
                   named_pairs_passed=sum(1 for g in named if g["passed"]),
                   named_pairs_total=len(named),
                   published_length=lengths,
                   published_stiffness=stiffnesses,
                   inside_one_bone_control=control),
        scope="A CONSTRUCTION: attachment sites are centroids of a vertex-vote partition of the "
              "structure's own surface, not measured insertion footprints; the stiffness is a "
              "published cadaver modulus times a derived cross-section, not a measured structural "
              "stiffness for this specimen. What it is NOT: an anatomical attachment, a "
              "subject-specific ligament property, or a claim that a structure absent from this "
              "list has no mechanics.")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print()
    print("tissue entities classified: %d" % sum(counts.values()))
    for label, _ in CLASSES:
        row = inventory[label]
        print("  %-20s %4d entities  %4d two-segment  %4d one-segment  %4d written"
              % (label, row["entities"], row["two_segment"], row["one_segment"],
                 row["instantiated_as_force_element"]))
    print("named-pair gate: %d/%d" % (manifest["gates"]["named_pairs_passed"],
                                      manifest["gates"]["named_pairs_total"]))
    for gate in named:
        if not gate["passed"]:
            print("  FAILED %s expected %s got %s" % (gate["name"], gate["expected"],
                                                      gate["derived"]))
    for row in lengths:
        print("  length %-52s %5.1f mm vs %5.1f published (%.2f)"
              % (row["name"][:52], row["derived_mm"], row["published_mm"], row["ratio"]))
    for row in stiffnesses:
        print("  stiffness %-49s %7.0f N vs %7.0f published (%.2f)"
              % (row["name"][:49], row["derived_n"], row["published_n"], row["ratio"]))
    print("force elements written: %d -> %s" % (len(elements), out / "ligaments.txt"))


def _tally(values):
    out = {}
    for value in values:
        out[value] = out.get(value, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


if __name__ == "__main__":
    main()
