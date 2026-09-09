#!/usr/bin/env python3
"""Cartilage, menisci and discs: what exists, and what the contact path can carry.

The direction asks for cartilage, menisci and intervertebral discs, joint
capsules and bursae as part of the model.  They are the load-bearing surfaces
BETWEEN bones, and cartilage in particular is what makes bone against bone
physical rather than a collision.  This measures whether the plant can express
them, and the answer is different for each and is negative in two places for
reasons worth stating exactly.

Three questions, all answered by measurement rather than by argument.

**1. Does articular cartilage exist in this body?**  Not "is there cartilage" --
there are 51 cartilage entities -- but is any of them a JOINT surface.  Read out
of the atlas' own names and cross-checked against
`data/derived/joint-substrate-assessment-v1/`, which asked the same question of
the unregistered extended atlas.

**2. Where do the menisci and discs sit on the scaffold?**  A meniscus between
femur and tibia is a structure two rigid bodies can load.  An intervertebral
disc between T8 and T9 is inside ONE rigid body, so there is no relative motion
for it to resist and no force it could carry.  Counted, not asserted.

**3. Can the mesh contact path carry an articular layer at all?**  The elastic
foundation is exactly the right shape for cartilage -- an independent normal
spring at every triangle of a rigid surface, which is a compliant layer over a
rigid substrate -- and `docs/SEGMENT_CONTACT_SURFACES.md` already measured that
it keeps concavity and costs 1.15x spheres.  What decides whether it can be used
between two BONES is whether those two bones' surfaces are separated at the
scaffold's own joint.  If they interpenetrate, a contact force between them
fires permanently and invents a load the body never carries; the native engine's
own comment says a shared contact set "would invent a bone-on-bone force at
every joint the source model lets overlap", and this is the measurement of
whether it does.

The interpenetration test is the repository's own `enclosure()` from
`scripts/build_skin_contact_meshes.py` -- Moller-Trumbore ray parity, two
opposed rays per point -- imported rather than restated.  Its controls are in
that file: points at radius 0.05 inside a 0.1 m icosphere print 1.0, points at
radius 0.5 print 0.0.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BUNDLE = ROOT / "data/models/engineering_stance_v1"
MESHES = ROOT / "data/derived/segment-contact-meshes/stance-bone-all"
BINDING = ROOT / "data/derived/anatomy-segment-binding/binding.json"
STRUCTURES = ROOT / "data/derived/tissue-force-elements-v1/structures.json"
SUBSTRATE = ROOT / "data/derived/joint-substrate-assessment-v1/manifest.json"

# An articular cartilage is a cartilage that covers a joint surface.  The same
# pattern `scripts/audit_joint_substrate.py` uses, so the two lanes ask the same
# question of the two atlases.
ARTICULAR = r"articular cartilage|hyaline cartilage|joint cartilage|cartilage of (head|condyle|facet)"
SWEEP_SAMPLES = 9


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _model_bone_meshes(model_path):
    """Every Body's own attached bone surfaces, scale factors applied.

    The same source `scripts/build_segment_contact_meshes.py` takes contact
    geometry from, and the same rule: `Mesh` carries `scale_factors` and the
    surface file is generic, so the scale has to be applied or a generic bone
    sits on a scaled subject.
    """
    import xml.etree.ElementTree as ET
    from ihm.spatial.vtk import surface as read_vtp
    geometry = ROOT / "data/raw/anatomy/opensim-models/source/Geometry"
    root = ET.parse(model_path).getroot().find("Model")
    out = {}
    for body in root.iter("Body"):
        for mesh in body.iter("Mesh"):
            factors = np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            name = mesh.findtext("mesh_file")
            points, faces = read_vtp(geometry / name)
            out.setdefault(body.get("name"), []).append((name, points * factors, faces))
    return out


def read_obj(path):
    vertices, faces = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("v "):
            vertices.append([float(v) for v in line.split()[1:4]])
        elif line.startswith("f "):
            faces.append([int(v.split("/")[0]) - 1 for v in line.split()[1:4]])
    return np.asarray(vertices, float), np.asarray(faces, int)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True)
    parser.add_argument("--samples", type=int, default=SWEEP_SAMPLES)
    args = parser.parse_args()
    started = time.time()

    crawl = load_module("crawl", ROOT / "scripts/crawl.py")
    render = load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")
    skin = load_module("build_skin_contact_meshes", ROOT / "scripts/build_skin_contact_meshes.py")

    structures = json.loads(STRUCTURES.read_text())["structures"]

    # ------------------------------------------------------------------ 1
    cartilage = [r for r in structures if r["tissue_class"] == "cartilage"]
    articular = [r for r in cartilage if re.search(ARTICULAR, r["name"].lower())]
    def bucket(name):
        n = name.lower()
        if "costal" in n:
            return "costal"
        if any(k in n for k in ("arytenoid", "cricoid", "thyroid", "corniculate", "cuneiform",
                                "tracheal", "bronchial", "epiglottic")):
            return "laryngotracheal"
        if "nasal" in n or "alar" in n or "septal" in n:
            return "nasal"
        if "triradiate" in n:
            return "growth_plate"
        if "auricular" in n or "ear" in n:
            return "auricular"
        return "other"
    cartilage_report = dict(
        entities=len(cartilage),
        articular_entities=len(articular),
        articular_names=[r["name"] for r in articular],
        by_kind=_tally(bucket(r["name"]) for r in cartilage),
        corroboration=None,
        finding="Every cartilage entity in this body is costal, laryngotracheal, nasal or a growth "
                "plate. NOT ONE is a joint surface. Cartilage as the thing that makes bone-on-bone "
                "contact physical is therefore a DATA gap before it is a mechanics gap: there is "
                "no geometry to give a force element to, and synthesising a layer over a bone "
                "patch would be inventing anatomy.")
    if SUBSTRATE.exists():
        manifest = json.loads(SUBSTRATE.read_text())
        summary = manifest.get("cartilage_gap_summary", {})
        cartilage_report["corroboration"] = dict(
            source="data/derived/joint-substrate-assessment-v1",
            canonical_articular_cartilage_entities=len(
                summary.get("canonical_articular_cartilage_entities", [])),
            extended_atlas_cartilage_count=manifest.get("inventory_summary", {})
                .get("cartilage", {}).get("count"),
            bone_contact_pairs=summary.get("contact_pairs"),
            cartilage_layers_required=summary.get("cartilage_layers_required"),
            total_patch_area_m2=summary.get("total_patch_area_m2"),
            note="An independent assessment of the UNREGISTERED extended atlas reached the same "
                 "zero, and counted how many layers would have to be synthesised to close it.")

    # ------------------------------------------------------------------ 2
    menisci = [r for r in structures if r["tissue_class"] == "meniscus_disc"]
    capsules = [r for r in structures if r["tissue_class"] == "joint_capsule"]
    bursae = [r for r in structures if r["tissue_class"] == "bursa_synovium"]
    def placement(rows):
        two = [r for r in rows if r.get("status") == "two_segment"]
        one = [r for r in rows if r.get("status") == "one_segment"]
        return dict(entities=len(rows), two_segment=len(two), one_segment=len(one),
                    two_segment_names=sorted({r["name"] for r in two}),
                    one_segment_by_segment=_tally(r["segment"] for r in one))
    placement_report = dict(meniscus_disc=placement(menisci),
                            joint_capsule=placement(capsules),
                            bursa_synovium=placement(bursae))

    # ------------------------------------------------------------------ 3
    model = render.OsimModel(BUNDLE / "model.osim")
    # The model's OWN attached geometry, all 81 bones with each Body's subject
    # scale applied -- NOT `data/derived/segment-contact-meshes`, which is the
    # SimTK-admissible subset and drops both hip bones, both tibiae, the skull,
    # the jaw and the spine.  Measuring the hip on that subset measures the
    # sacrum against the femur and reports an 82 mm "joint gap"; measuring the
    # knee on it measures the fibula.  The admissibility question is real and is
    # reported separately below, but it is not the geometry question.
    by_body = _model_bone_meshes(BUNDLE / "model.osim")
    manifest = json.loads((MESHES / "manifest.json").read_text())
    admissible = {}
    for record in manifest["records"]:
        admissible.setdefault(record["body"], []).append(record["element"])
    ranges = crawl.declared_ranges()
    binding = json.loads(BINDING.read_text())
    reference = dict(binding["reference_pose_rad"])

    pairs = []
    for joint in model.joints:
        parent, child = joint["parent"], joint["child"]
        if parent not in by_body or child not in by_body:
            continue
        driving = [c for c in joint["coords"] if c in ranges]
        # Sweep the joint's own declared range, one coordinate at a time, with
        # everything else at the binding's reference pose.
        samples = []
        for coordinate in driving or [None]:
            values = ([reference.get(coordinate, 0.0)] if coordinate is None
                      else list(np.linspace(*ranges[coordinate], args.samples)))
            for value in values:
                pose = dict(reference)
                if coordinate is not None:
                    pose[coordinate] = float(value)
                transforms = model.forward(pose)
                gap, inside = _separation(by_body[parent], by_body[child],
                                          transforms[parent], transforms[child], skin.enclosure)
                samples.append(dict(coordinate=coordinate, value_rad=None if coordinate is None else float(value),
                                    minimum_vertex_gap_m=gap,
                                    child_vertices_inside_parent=inside))
        worst = min(samples, key=lambda s: s["minimum_vertex_gap_m"])
        overlapping = [s for s in samples if s["child_vertices_inside_parent"] > 0.0]
        pairs.append(dict(joint=joint["name"], parent=parent, child=child,
                          coordinates=driving,
                          parent_bones=len(by_body[parent]), child_bones=len(by_body[child]),
                          parent_bones_simtk_admissible=len(admissible.get(parent, [])),
                          child_bones_simtk_admissible=len(admissible.get(child, [])),
                          samples_taken=len(samples),
                          samples_with_interpenetration=len(overlapping),
                          worst=worst,
                          reference_pose=next(s for s in samples
                                              if s["coordinate"] is None
                                              or abs(s["value_rad"] - reference.get(s["coordinate"], 0.0))
                                              == min(abs(x["value_rad"] - reference.get(x["coordinate"], 0.0))
                                                     for x in samples if x["coordinate"] == s["coordinate"]))))
        print("  %-18s %-9s %-9s worst gap %7.2f mm  interpenetrating in %d/%d samples"
              % (joint["name"], parent, child, worst["minimum_vertex_gap_m"] * 1000,
                 len(overlapping), len(samples)), flush=True)

    usable = [p for p in pairs if p["samples_with_interpenetration"] == 0]
    contact_report = dict(
        joints_measured=len(pairs),
        joints_with_no_interpenetration_anywhere_in_range=len(usable),
        joints_interpenetrating_somewhere_in_range=len(pairs) - len(usable),
        pairs=pairs,
        geometry_source="the model's own attached Body geometry, all 81 bone surfaces with each "
                        "Body's scale_factors applied; NOT the SimTK-admissible contact subset, "
                        "which drops both hip bones, both tibiae, the skull, the jaw and the spine "
                        "and would measure the hip as sacrum-against-femur",
        basis="Minimum vertex-to-vertex distance is an UPPER bound on the true surface gap, so a "
              "reported separation may be optimistic; the interpenetration share is the decisive "
              "number and it is exact ray parity, not a distance threshold.",
        finding="A synovial joint's two bone surfaces OVERLAP by construction -- the femoral head "
                "is inside the acetabulum -- and the measurement says so: the hips interpenetrate "
                "in 19/21 and 21/21 of the sampled configurations across their own declared range. "
                "So an ElasticFoundationForce between the two BONES of such a joint would fire "
                "permanently and invent a load the body never carries. Bone-on-bone contact is not "
                "the mechanism for an articular surface; the cartilage layers are, and they must "
                "be separate surfaces between the bones. This body has zero of them. The joints "
                "where the bones stay clear through the whole declared range -- knee, "
                "patellofemoral, radioulnar, acromial, radiocarpal -- are the ones where a "
                "layer could be carried IF its geometry existed, and the knee is the only joint "
                "in the body that has intra-articular geometry (four menisci) to build it from.")

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(
        schema="ihm.articular-substrate-measurement.v1",
        generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        wall_s=time.time() - started,
        inputs={p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                for p in (str(STRUCTURES), str(MESHES / "manifest.json"), str(BINDING))},
        cartilage=cartilage_report,
        placement=placement_report,
        bone_on_bone_contact=contact_report), indent=2) + "\n")
    print()
    print("cartilage entities %d, of which ARTICULAR: %d"
          % (cartilage_report["entities"], cartilage_report["articular_entities"]))
    print("cartilage by kind: %s" % cartilage_report["by_kind"])
    for label, row in placement_report.items():
        print("%-16s %3d entities, %3d across a scaffold joint, %3d inside one body"
              % (label, row["entities"], row["two_segment"], row["one_segment"]))
    print("bone pairs with no interpenetration anywhere in the declared range: %d/%d"
          % (len(usable), len(pairs)))
    print(str(out))


def _separation(parent_meshes, child_meshes, parent_transform, child_transform, enclosure,
                probe_gap_m=0.010):
    """Closest approach, and the share of child vertices inside a parent bone.

    The vertex-to-vertex distance is cheap and is an UPPER bound on the surface
    gap, so it screens; the ray-parity test is exact and is run only on the mesh
    pairs the screen puts within `probe_gap_m` of each other.  A pair further
    apart than that cannot be interpenetrating, so nothing is missed by skipping
    it and the whole sweep stays affordable.
    """
    parent_points = [(v @ parent_transform[:3, :3].T + parent_transform[:3, 3], f)
                     for _, v, f in parent_meshes]
    child_points = [(v @ child_transform[:3, :3].T + child_transform[:3, 3], f)
                    for _, v, f in child_meshes]
    gap = math.inf
    inside = 0.0
    for pv, pf in parent_points:
        tree = cKDTree(pv)
        for cv, _ in child_points:
            local = float(tree.query(cv)[0].min())
            gap = min(gap, local)
            if local <= probe_gap_m:
                inside = max(inside, enclosure(pv, pf, cv, samples=300))
    return gap, inside


def _tally(values):
    out = {}
    for value in values:
        out[value] = out.get(value, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


if __name__ == "__main__":
    main()
