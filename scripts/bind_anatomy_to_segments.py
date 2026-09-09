#!/usr/bin/env python3
"""Bind every anatomical entity to the OpenSim body segment it rides on.

The simulated body carries 4,000 BodyParts3D-derived surfaces.  The controller
drives 22 OpenSim bodies.  Nothing connected the two, so brain-driven motion
could not reach the anatomy: every embodied video was the 22-capsule skeleton.
This script builds the missing map and writes it as its own artifact.

What it does
------------
1.  Names the 22 segments' *anatomical* bone groups from the atlas' own labels
    ("right femur" -> femur_r, every rib and vertebra -> torso, ...).  334 of
    the 4,000 entities are bones and every one of them is named.  Nothing
    geometric decides these; they are the reference shapes everything else is
    measured against.

2.  Registers the OpenSim model to the atlas.  Both the reference *pose* (33
    coordinates) and one global similarity (scale, rotation, translation) are
    fitted by least squares against the 22 named bone groups, because the
    OpenSim specimen is not the BodyParts3D specimen and its default pose is
    not the atlas pose.  The residual is reported, not hidden.

3.  Assigns the remaining 3,666 entities by vertex vote: each entity's surface
    is sampled, every sample takes the segment whose anatomical bone group is
    nearest, and the entity takes the majority.  The vote share is kept as
    `coherence` -- an entity that straddles a joint has a low one, and a rigid
    binding will tear it.

Gates (printed and stored; see --report)
----------------------------------------
* held-out bones: 241 named bones sit in a segment that has other bones, so
  each can be dropped from its own segment's reference shape and re-assigned
  purely geometrically.  That is a real test of the assignment rule.
* muscles: 80 native muscles have a canonical entity and a declared pair of
  attachment bodies (catalog.json).  A binding is counted correct when it
  names one of them.  These entities are never bound by name.
* an independent second opinion: the same vote run against the *OpenSim* bone
  meshes pushed through the registration instead of the atlas' own bones.
  Agreement between two different reference geometries is corroboration;
  disagreement localises the doubt.
* hand-checked cases with a known answer are listed in KNOWN_CASES.

Usage
-----
    python scripts/bind_anatomy_to_segments.py
    python scripts/bind_anatomy_to_segments.py --out data/derived/my-binding
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.spatial.vtk import surface as read_vtp  # noqa: E402

ANATOMY = ROOT / "data/derived/canonical/anatomy.json"
GEOMETRY = ROOT / "data/derived/canonical/geometry"
MODEL = ROOT / "data/models/engineering_stance_v1/model.osim"
CATALOG = ROOT / "data/models/engineering_stance_v1/catalog.json"
MECHANICS = ROOT / "data/derived/canonical/mechanics.json"
OSIM_GEOMETRY = ROOT / "data/raw/anatomy/opensim-models/source/Geometry"
OUT = ROOT / "data/derived/anatomy-segment-binding"

MAX_VERTEX_SAMPLES = 512


# --------------------------------------------------------------------------
# 1. The 22 segments, named from the atlas' own labels
# --------------------------------------------------------------------------
CARPALS = ("scaphoid", "lunate", "triquetral", "triquetrum", "pisiform",
           "trapezium", "trapezoid", "capitate", "hamate")
CRANIAL = ("frontal bone", "parietal bone", "occipital bone", "temporal bone",
           "sphenoid", "ethmoid", "vomer", "maxilla", "palatine", "nasal bone",
           "lacrimal", "zygomatic", "inferior nasal concha", "mandible",
           "hyoid", "incus", "malleus", "stapes", "sinus of")
DENTAL = ("tooth", "incisor", "canine", "molar", "premolar")


def named_segment(name: str) -> str | None:
    """Segment of a *bone*, from its atlas label alone.  None if not a rule.

    The rules follow the OpenSim model's own geometry: `torso` is the HAT body
    (hat_skull, hat_jaw, hat_spine, hat_ribs_scap), so skull, jaw, teeth,
    every vertebra, every rib, sternum, scapula and clavicle ride it -- the
    model has no neck or shoulder-girdle joint.  `calcn` carries r_foot.vtp,
    which is calcaneus + tarsals + metatarsals; `toes` carries r_bofoot.vtp,
    the phalanges only.
    """
    n = name.lower()
    side = "r" if "right" in n else ("l" if "left" in n else None)

    if "hip bone" in n or n in ("sacrum", "coccyx"):
        return "pelvis"
    if side is None:
        if ("vertebra" in n or n in ("atlas", "axis", "atlas c1", "axis c2")
                or "sternum" in n or "manubrium" in n or "xiphoid" in n):
            return "torso"
        if any(k in n for k in CRANIAL) or any(k in n for k in DENTAL):
            return "torso"
        return None
    if n.endswith("femur"):
        return "femur_" + side
    if n.endswith("patella"):
        return "patella_" + side
    if n.endswith("tibia") or n.endswith("fibula"):
        return "tibia_" + side
    if n.endswith("talus"):
        return "talus_" + side
    if ("calcaneus" in n or "cuboid" in n or "navicular" in n
            or "cuneiform" in n or "sesamoid" in n or "metatarsal" in n):
        return "calcn_" + side
    if "toe" in n or "finger of foot" in n:
        return "toes_" + side
    if n.endswith("humerus"):
        return "humerus_" + side
    if n.endswith("ulna"):
        return "ulna_" + side
    if n.endswith("radius"):
        return "radius_" + side
    if (any(c in n for c in CARPALS) or "metacarpal" in n
            or "finger of hand" in n or "thumb" in n or "finger" in n):
        return "hand_" + side
    if "rib" in n or "scapula" in n or "clavicle" in n or "vertebra" in n:
        return "torso"
    if any(k in n for k in CRANIAL) or any(k in n for k in DENTAL):
        return "torso"
    return None


# Entities with an answer that anatomy fixes and this script never sees.
# `accept` is the set a correct binding may name (a muscle spanning a joint
# has two defensible answers; a rigid binding must pick one of them).
KNOWN_CASES = [
    ("right femur", {"femur_r"}),
    ("left femur", {"femur_l"}),
    ("right calcaneus", {"calcn_r"}),
    ("left twelfth rib", {"torso"}),
    ("stomach", {"torso"}),
    ("spleen", {"torso"}),
    ("left kidney", {"torso"}),
    ("superior lobe of right lung", {"torso"}),
    ("midbrain", {"torso"}),
    ("urinary bladder", {"pelvis", "torso"}),
    ("right vastus lateralis", {"femur_r"}),
    ("left vastus medialis", {"femur_l"}),
    ("right soleus", {"tibia_r"}),
    ("left tibialis anterior", {"tibia_l"}),
    ("right gluteus maximus", {"pelvis", "femur_r"}),
    ("long head of right biceps brachii", {"humerus_r"}),
    ("acromial part of left deltoid", {"humerus_l", "torso"}),
    ("right calcaneal tendon", {"tibia_r", "calcn_r"}),
    ("right femoral artery", {"femur_r"}),
    ("right brachial artery", {"humerus_r"}),
    ("descending thoracic aorta", {"torso"}),
]


# --------------------------------------------------------------------------
# 2. OpenSim registration
# --------------------------------------------------------------------------
def load_osim():
    spec = importlib.util.spec_from_file_location(
        "_render_body_3d", ROOT / "scripts/render_body_3d.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.OsimModel(MODEL)


def osim_body_meshes():
    """Body-local point clouds of every OpenSim body, scale factors applied."""
    root = ET.parse(MODEL).getroot().find("Model")
    out = {}
    for body in root.iter("Body"):
        parts = []
        for mesh in body.iter("Mesh"):
            frame = (mesh.findtext("socket_frame") or "..").strip()
            if frame != "..":
                raise ValueError(f"mesh on an offset frame: {body.get('name')}")
            factors = np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            points, _ = read_vtp(OSIM_GEOMETRY / mesh.findtext("mesh_file"))
            parts.append(points * factors)
        if parts:
            out[body.get("name")] = np.concatenate(parts)
    return out


def principal_frame(points):
    mean = points.mean(0)
    values, vectors = np.linalg.eigh(np.cov((points - mean).T))
    return mean, vectors[:, ::-1], values[::-1]


def umeyama(a, b):
    ca, cb = a.mean(0), b.mean(0)
    x, y = a - ca, b - cb
    u, s, vt = np.linalg.svd(x.T @ y)
    d = np.diag([1.0, 1.0, np.sign(np.linalg.det(vt.T @ u.T))])
    r = vt.T @ d @ u.T
    scale = float((s @ np.diag(d)).sum() / (x ** 2).sum())
    return scale, r, cb - scale * r @ ca


def register(model, local_clouds, anatomy_clouds, segments):
    """Fit the reference pose and one similarity osim-ground -> atlas frame.

    The atlas is a different specimen in a different pose, so a pure rigid map
    cannot exist.  Both the pose and the similarity are free; the residual is
    the honest statement of what is left over.
    """
    root = ET.parse(MODEL).getroot().find("Model")
    default = {c.get("name"): float(c.findtext("default_value") or 0.0)
               for c in root.iter("Coordinate")}

    local = {s: principal_frame(local_clouds[s]) for s in segments}
    atlas = {s: principal_frame(anatomy_clouds[s]) for s in segments}
    # Only principal axes that are actually separated carry orientation.
    axes = {s: [i for i in range(3)
                if (local[s][2][i] / local[s][2][i + 1] if i < 2
                    else local[s][2][1] / local[s][2][2]) > 1.3]
            for s in segments}

    rest = model.forward(default)
    a = np.array([rest[s][:3, :3] @ local[s][0] + rest[s][:3, 3] for s in segments])
    b = np.array([atlas[s][0] for s in segments])
    scale0, rot0, t0 = umeyama(a, b)
    signs = {s: np.where(np.sum(((rot0 @ rest[s][:3, :3]) @ local[s][1])
                                * atlas[s][1], axis=0) >= 0, 1.0, -1.0)
             for s in segments}

    # *_beta coordinates are driven by a coupler constraint, and the pelvis
    # translations are already carried by the similarity's translation.
    free = [k for k in default
            if not k.endswith("_beta")
            and k not in ("pelvis_tx", "pelvis_ty", "pelvis_tz")]
    arm, weight = 0.05, 0.5

    def unpack(x):
        q = dict(default)
        q.update({k: x[7 + i] for i, k in enumerate(free)})
        for k in ("pelvis_tx", "pelvis_ty", "pelvis_tz"):
            q[k] = 0.0
        return float(np.exp(x[0])), Rotation.from_rotvec(x[1:4]).as_matrix(), x[4:7], q

    def residual(x):
        scale, rot, t, q = unpack(x)
        transforms = model.forward(q)
        rows = []
        for s in segments:
            tb = transforms[s]
            point = scale * rot @ (tb[:3, :3] @ local[s][0] + tb[:3, 3]) + t
            rows.append(point - atlas[s][0])
            oriented = rot @ (tb[:3, :3] @ local[s][1])
            for i in axes[s]:
                rows.append(weight * arm
                            * (oriented[:, i] - atlas[s][1][:, i] * signs[s][i]))
        return np.concatenate(rows)

    x0 = np.concatenate([[np.log(scale0)], Rotation.from_matrix(rot0).as_rotvec(),
                         t0, [default[k] for k in free]])
    lo = np.full(x0.shape, -np.inf)
    hi = np.full(x0.shape, np.inf)
    lo[7:] = np.array([default[k] for k in free]) - 1.0
    hi[7:] = np.array([default[k] for k in free]) + 1.0
    solution = least_squares(residual, x0, bounds=(lo, hi), xtol=1e-13, ftol=1e-13)
    scale, rot, t, q = unpack(solution.x)

    transforms = model.forward(q)
    per_segment = {}
    for s in segments:
        tb = transforms[s]
        point = scale * rot @ (tb[:3, :3] @ local[s][0] + tb[:3, 3]) + t
        per_segment[s] = float(np.linalg.norm(point - atlas[s][0]))
    values = np.array(list(per_segment.values()))
    fit = {
        "method": "joint least squares over 33 model coordinates and one "
                  "similarity (scale, rotation, translation); residual is "
                  "22 bone-group centroids plus separated principal axes",
        "scale": scale,
        "rotation": rot.tolist(),
        "translation_m": t.tolist(),
        "reference_pose_rad": q,
        "reference_pose_delta_from_model_default_deg": {
            k: float((q[k] - default[k]) * 180 / np.pi) for k in free
            if abs(q[k] - default[k]) > 1e-4},
        "segment_centroid_residual_m": per_segment,
        "segment_centroid_residual_rms_m": float(np.sqrt((values ** 2).mean())),
        "segment_centroid_residual_max_m": float(values.max()),
        "interpretation": "Cross-specimen geometric correspondence residual. "
                          "The OpenSim model is a scaled Rajagopal subject; "
                          "the atlas is the BodyParts3D adult male. This is "
                          "not measurement error and not subject variation.",
    }
    return scale, rot, t, q, fit


# --------------------------------------------------------------------------
# 3. Assignment
# --------------------------------------------------------------------------
GEOMETRY_PATH = {}


def read_geometry(entity_id):
    """Vertices of an entity's reference surface, in the atlas frame.

    The path comes from the atlas record, never from the id: a promoted
    Z-Anatomy surface does not live at `<id>.json.gz`.
    """
    path = GEOMETRY_PATH.get(entity_id)
    if path is None:
        path = GEOMETRY / f"{entity_id}.json.gz"
    payload = json.loads(gzip.decompress(Path(path).read_bytes()))
    return np.asarray(payload["positions"], float).reshape(-1, 3)


def sample(points, limit=MAX_VERTEX_SAMPLES):
    if len(points) <= limit:
        return points
    step = len(points) // limit
    return points[::step][:limit]


def vote(points, trees, segments):
    """Nearest segment per sampled vertex; the entity takes the majority."""
    distances = np.stack([tree.query(points)[0] for tree in trees], axis=1)
    nearest = np.argmin(distances, axis=1)
    counts = np.bincount(nearest, minlength=len(segments))
    winner = int(np.argmax(counts))
    return (segments[winner],
            float(counts[winner]) / len(points),
            float(np.mean(distances[:, winner])),
            float(np.min(distances[:, winner])),
            segments[int(np.argsort(-counts)[1])] if (counts > 0).sum() > 1 else None)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git_sha():
    sha = os.environ.get("IBM_GIT_SHA")
    if sha:
        return sha
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "git-unknown"


def segment_path(model, a, b):
    """Segments on the kinematic chain between two bodies, endpoints included.

    A two-joint muscle's belly does not lie on either bone it attaches to; it
    lies on the bone in between.  `attachment_bodies` names the endpoints, so
    the endpoints alone are the wrong acceptance set for such a muscle.
    """
    parent = {j["child"]: j["parent"] for j in model.joints}

    def chain(x):
        out = [x]
        while x in parent:
            x = parent[x]
            out.append(x)
        return out

    ca, cb = chain(a), chain(b)
    common = next((x for x in ca if x in cb), None)
    if common is None:
        return set(ca) | set(cb)
    return set(ca[:ca.index(common) + 1]) | set(cb[:cb.index(common) + 1])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--samples", type=int, default=MAX_VERTEX_SAMPLES)
    args = parser.parse_args()
    started = time.time()

    anatomy = json.loads(ANATOMY.read_text())
    entities = anatomy["entities"]
    GEOMETRY_PATH.update({e["id"]: ROOT / e["reference_geometry"]["path"]
                          for e in entities if e.get("reference_geometry")})
    print(f"atlas: {len(entities)} entities, frame {anatomy['frame']['id']}")

    # -- named bone groups ------------------------------------------------
    groups = {}
    for entity in entities:
        if entity["role"] != "rigid_bone":
            continue
        segment = named_segment(entity["name"])
        if segment is None:
            raise ValueError("unnamed bone: " + entity["name"])
        groups.setdefault(segment, []).append(entity["id"])
    segments = sorted(groups)
    print(f"segments: {len(segments)}; named bones: "
          f"{sum(len(v) for v in groups.values())}")

    model = load_osim()
    missing = set(b for b in model.bodies if b != "ground") - set(segments)
    if missing:
        raise ValueError(f"model bodies with no named bone group: {sorted(missing)}")

    bone_points = {i: read_geometry(i) for ids in groups.values() for i in ids}
    anatomy_clouds = {s: np.concatenate([bone_points[i] for i in groups[s]])
                      for s in segments}

    # -- registration -----------------------------------------------------
    local_clouds = osim_body_meshes()
    scale, rot, translation, pose, fit = register(model, local_clouds,
                                                  anatomy_clouds, segments)
    print(f"registration: scale {scale:.4f}, centroid residual "
          f"RMS {fit['segment_centroid_residual_rms_m']*1000:.1f} mm, "
          f"max {fit['segment_centroid_residual_max_m']*1000:.1f} mm")

    rest = model.forward(pose)
    similarity = np.eye(4)
    similarity[:3, :3] = scale * rot
    similarity[:3, 3] = translation
    osim_clouds = {}
    for s in segments:
        world = (rest[s][:3, :3] @ local_clouds[s].T).T + rest[s][:3, 3]
        osim_clouds[s] = (similarity[:3, :3] @ world.T).T + similarity[:3, 3]

    # Joint centres, in the atlas frame, and how far each lands from the
    # anatomical joint it is supposed to be.  A rigid binding rotates the
    # anatomy about *these* points, so this is the error that matters.
    joint_rows = []
    for joint in model.joints:
        parent, child = joint["parent"], joint["child"]
        if parent not in segments or child not in segments:
            continue
        centre = rest[parent] @ joint["Xpf"]
        centre = similarity[:3, :3] @ centre[:3, 3] + similarity[:3, 3]
        tree_p = cKDTree(anatomy_clouds[parent])
        tree_c = cKDTree(anatomy_clouds[child])
        pairs = tree_p.query_ball_tree(tree_c, r=0.02)
        contact = [anatomy_clouds[parent][i] for i, hits in enumerate(pairs) if hits]
        if contact:
            proxy = np.mean(contact, axis=0)
        else:
            distance, index = tree_c.query(anatomy_clouds[parent])
            k = int(np.argmin(distance))
            proxy = (anatomy_clouds[parent][k] + anatomy_clouds[child][index[k]]) / 2
        joint_rows.append({
            "joint": joint["name"], "parent": parent, "child": child,
            "centre_atlas_m": centre.tolist(),
            "anatomical_proxy_m": proxy.tolist(),
            "offset_m": float(np.linalg.norm(centre - proxy)),
        })
    offsets = np.array([r["offset_m"] for r in joint_rows])
    print(f"joint centres: median offset from the anatomical joint "
          f"{np.median(offsets)*1000:.1f} mm, max {offsets.max()*1000:.1f} mm")

    # -- assignment -------------------------------------------------------
    trees = [cKDTree(anatomy_clouds[s]) for s in segments]
    osim_trees = [cKDTree(osim_clouds[s]) for s in segments]
    owner = {i: s for s in segments for i in groups[s]}

    rows = {}
    agree = disagree = 0
    print("assigning...", flush=True)
    for n, entity in enumerate(entities):
        ident = entity["id"]
        points = bone_points.get(ident)
        if points is None:
            points = read_geometry(ident)
        points = sample(points, args.samples)
        segment, coherence, mean_d, min_d, runner = vote(points, trees, segments)
        second, _, _, _, _ = vote(points, osim_trees, segments)
        if second == segment:
            agree += 1
        else:
            disagree += 1
        rows[ident] = {
            "segment": owner.get(ident, segment),
            "basis": "named_bone" if ident in owner else "nearest_bone_group_vertex_vote",
            "geometric_segment": segment,
            "coherence": coherence,
            "mean_distance_m": mean_d,
            "min_distance_m": min_d,
            "runner_up": runner,
            "opensim_reference_segment": second,
            "name": entity["name"],
            "system": entity["system"],
            "role": entity["role"],
            "sampled_vertices": int(len(points)),
        }
        if n % 500 == 0:
            print(f"  {n}/{len(entities)}  {time.time()-started:.0f}s", flush=True)

    # -- gate: held-out bones --------------------------------------------
    heldout = {"tested": 0, "correct": 0, "failures": []}
    for segment in segments:
        ids = groups[segment]
        if len(ids) < 2:
            continue
        for ident in ids:
            rest_points = np.concatenate([bone_points[j] for j in ids if j != ident])
            local_trees = [cKDTree(rest_points) if s == segment else trees[k]
                           for k, s in enumerate(segments)]
            got, _, _, _, _ = vote(sample(bone_points[ident], args.samples),
                                   local_trees, segments)
            heldout["tested"] += 1
            if got == segment:
                heldout["correct"] += 1
            else:
                share = len(bone_points[ident]) / len(anatomy_clouds[segment])
                heldout["failures"].append(
                    {"id": ident, "name": rows[ident]["name"],
                     "truth": segment, "assigned": got,
                     "vertex_share_of_its_segment": float(share),
                     "note": "the held-out bone is this segment's dominant "
                             "shape; removing it removes the segment"
                             if share > 0.3 else "genuine confusion"})
    heldout["misassignment_rate"] = (
        1 - heldout["correct"] / heldout["tested"]) if heldout["tested"] else None
    print(f"held-out bones: {heldout['correct']}/{heldout['tested']} correct, "
          f"misassignment {heldout['misassignment_rate']*100:.1f}%")

    # -- gate: muscles with declared attachment bodies --------------------
    catalog = {m["id"]: m for m in json.loads(CATALOG.read_text())}
    native = json.loads(MECHANICS.read_text())["native_muscles"]
    muscle = {"tested": 0, "endpoint_correct": 0, "path_correct": 0,
              "criterion": "The muscle's bound segment must be one of the "
                           "bodies the OpenSim path attaches to (endpoint), "
                           "or one of the segments on the kinematic chain "
                           "between them (path). A two-joint muscle's belly "
                           "lies on an intervening bone, so `path` is the "
                           "criterion a rigid binding should be held to and "
                           "`endpoint` is reported for comparison.",
              "endpoint_failures": [], "path_failures": []}
    for record in native:
        ident = record["canonical_entity_id"]
        entry = catalog.get(record["source_name"])
        if entry is None or ident not in rows:
            continue
        bodies = list(entry["attachment_bodies"])
        accept = set(bodies)
        path = set()
        for i, a in enumerate(bodies):
            for b in bodies[i + 1:]:
                path |= segment_path(model, a, b)
        path |= accept
        got = rows[ident]["segment"]
        muscle["tested"] += 1
        muscle["endpoint_correct"] += int(got in accept)
        muscle["path_correct"] += int(got in path)
        row = {"id": ident, "name": rows[ident]["name"],
               "source_name": record["source_name"],
               "attachment_bodies": sorted(accept), "path": sorted(path),
               "assigned": got, "coherence": rows[ident]["coherence"]}
        if got not in accept:
            muscle["endpoint_failures"].append(row)
        if got not in path:
            muscle["path_failures"].append(row)
    muscle["endpoint_misassignment_rate"] = (
        1 - muscle["endpoint_correct"] / muscle["tested"]) if muscle["tested"] else None
    muscle["misassignment_rate"] = (
        1 - muscle["path_correct"] / muscle["tested"]) if muscle["tested"] else None
    print(f"muscles vs declared attachments: endpoint "
          f"{muscle['endpoint_correct']}/{muscle['tested']} "
          f"({muscle['endpoint_misassignment_rate']*100:.1f}% miss), "
          f"kinematic path {muscle['path_correct']}/{muscle['tested']} "
          f"({muscle['misassignment_rate']*100:.1f}% miss)")
    for f in muscle["path_failures"]:
        print(f"    OFF-PATH {f['source_name']}: {f['name']} -> "
              f"{f['assigned']}, attaches {f['attachment_bodies']}")

    # -- gate: hand-checked cases -----------------------------------------
    by_name = {}
    for ident, row in rows.items():
        by_name.setdefault(row["name"].lower(), []).append(ident)
    hand = {"tested": 0, "correct": 0, "cases": []}
    for name, accept in KNOWN_CASES:
        ids = by_name.get(name, [])
        if not ids:
            hand["cases"].append({"name": name, "status": "no such entity"})
            continue
        ident = ids[0]
        ok = rows[ident]["segment"] in accept
        hand["tested"] += 1
        hand["correct"] += int(ok)
        hand["cases"].append({"name": name, "id": ident, "accept": sorted(accept),
                              "assigned": rows[ident]["segment"],
                              "basis": rows[ident]["basis"],
                              "coherence": rows[ident]["coherence"], "ok": ok})
    print(f"hand-checked cases: {hand['correct']}/{hand['tested']} correct")
    for case in hand["cases"]:
        if not case.get("ok", True):
            print(f"    MISS {case['name']}: {case['assigned']} "
                  f"not in {case['accept']}")

    # -- distributions and honesty flags ----------------------------------
    geometric = [r for r in rows.values() if r["basis"] != "named_bone"]
    distance = np.array([r["min_distance_m"] for r in geometric])
    coherence = np.array([r["coherence"] for r in geometric])
    percentiles = [50, 75, 90, 95, 99, 100]
    by_system = {}
    for row in rows.values():
        s = by_system.setdefault(row["system"], {"n": 0, "torn": 0, "far": 0})
        s["n"] += 1
        s["torn"] += int(row["coherence"] < 0.6)
        s["far"] += int(row["min_distance_m"] > 0.05)

    torn = sorted((r for r in rows.values() if r["coherence"] < 0.6),
                  key=lambda r: r["coherence"])
    report = {
        "schema": "ihm.anatomy-segment-binding-report.v1",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "entities_total": len(entities),
        "entities_assigned": len(rows),
        "assignment_fraction": len(rows) / len(entities),
        "bound_by_name": sum(1 for r in rows.values() if r["basis"] == "named_bone"),
        "bound_geometrically": len(geometric),
        "distance_to_bound_segment_m": {
            f"p{p}": float(np.percentile(distance, p)) for p in percentiles},
        "distance_mean_m": float(distance.mean()),
        "coherence_low_percentiles": {
            f"p{p}": float(np.percentile(coherence, p))
            for p in (1, 5, 10, 25, 50)},
        "coherence_definition": "share of an entity's sampled vertices whose "
                                "own nearest segment is the one it was bound "
                                "to; 1.0 means the whole surface agrees",
        "coherence_below_0.6": int((coherence < 0.6).sum()),
        "opensim_reference_agreement": agree / (agree + disagree),
        "gate_heldout_bones": heldout,
        "gate_muscles": muscle,
        "gate_hand_checked": hand,
        "joint_centres": joint_rows,
        "joint_centre_offset_median_m": float(np.median(offsets)),
        "joint_centre_offset_max_m": float(offsets.max()),
        "by_system": by_system,
        "most_incoherent": [
            {k: t[k] for k in ("name", "system", "segment", "coherence", "runner_up")}
            for t in torn[:40]],
        "limitations": [
            "One rigid segment per entity. An entity whose surface spans a "
            "joint (a two-joint muscle, a long vessel or nerve, the skin) "
            "cannot ride one segment without tearing; `coherence` measures "
            "exactly how badly, and the renderer must not pretend otherwise.",
            "The skin, its three layers and the lymphatic network graph are "
            "excluded from rigid binding and listed in `excluded`: for the "
            "skin a continuous linear-blend attachment over the same 22 "
            "segments already exists at "
            "data/derived/canonical/continuous_surface_binding.json.gz.",
            "The OpenSim model has no neck, shoulder-girdle or spine joint. "
            "Skull, jaw, teeth, every vertebra, every rib, both scapulae and "
            "both clavicles ride `torso` as one rigid body, because that is "
            "what the model actually is. The head cannot nod.",
            "The registration is cross-specimen. Its centroid residual and "
            "the offset of every joint centre from the anatomical joint are "
            "reported above; neither is measurement error.",
        ],
    }

    # Entities that no single segment can honestly carry.
    excluded = {}
    for ident, row in rows.items():
        if row["role"] in ("skin", "skin_layer", "lymphatic_network"):
            excluded[ident] = {
                "name": row["name"], "role": row["role"],
                "reason": ("whole-body surface or graph; no single rigid "
                           "segment carries it"),
                "alternative": ("data/derived/canonical/"
                                "continuous_surface_binding.json.gz")
                if row["role"] in ("skin", "skin_layer") else None,
                "nearest_segment": row["segment"], "coherence": row["coherence"]}
    for ident in excluded:
        rows[ident]["segment"] = None
        rows[ident]["basis"] = "excluded"
    report["excluded"] = excluded
    report["entities_assigned"] = sum(1 for r in rows.values() if r["segment"])
    report["assignment_fraction"] = report["entities_assigned"] / len(entities)
    print(f"excluded {len(excluded)} whole-body entities: "
          f"{[v['name'] for v in excluded.values()]}")
    print(f"assigned {report['entities_assigned']}/{len(entities)} "
          f"({report['assignment_fraction']*100:.2f}%)")
    print(f"coherence < 0.6 (rigid binding will tear): "
          f"{report['coherence_below_0.6']} entities")

    counts = {}
    for row in rows.values():
        counts[row["segment"]] = counts.get(row["segment"], 0) + 1
    report["entities_per_segment"] = {k: v for k, v in
                                      sorted(counts.items(), key=lambda kv: -kv[1])}

    binding = {
        "schema": "ihm.anatomy-segment-binding.v1",
        "model_id": anatomy["model_id"],
        "frame": anatomy["frame"]["id"],
        "segments": segments,
        "segment_named_bones": groups,
        "registration": fit,
        "reference_pose_rad": pose,
        "similarity_atlas_from_opensim_ground": similarity.tolist(),
        "entities": rows,
        "centroids_m": {e["id"]: e["centroid_m"] for e in entities},
        "provenance": {
            "git_sha": git_sha(),
            "script": "scripts/bind_anatomy_to_segments.py",
            "script_sha256": digest(__file__),
            "anatomy": {"path": str(ANATOMY.relative_to(ROOT)),
                        "sha256": digest(ANATOMY)},
            "model": {"path": str(MODEL.relative_to(ROOT)), "sha256": digest(MODEL)},
            "catalog": {"path": str(CATALOG.relative_to(ROOT)),
                        "sha256": digest(CATALOG)},
            "opensim_geometry_dir": str(OSIM_GEOMETRY.relative_to(ROOT)),
            "vertex_samples_per_entity": args.samples,
        },
        "scope": "A geometric prior for which rigid segment carries each "
                 "anatomical surface. It is not an anatomical attachment, an "
                 "insertion site, or a claim that the segment boundary is a "
                 "tissue boundary.",
    }

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "binding.json").write_text(json.dumps(binding) + "\n")
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {args.out/'binding.json'} and {args.out/'report.json'} "
          f"in {time.time()-started:.0f}s")


if __name__ == "__main__":
    main()
