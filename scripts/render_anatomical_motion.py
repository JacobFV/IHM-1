#!/usr/bin/env python3
"""Play a joint trajectory on the real anatomy instead of the 22-capsule stick.

The controller produces 33 joint coordinates per frame.  `binding.json`
(scripts/bind_anatomy_to_segments.py) says which of the 22 OpenSim segments
each of the 3,995 anatomical surfaces rides.  This script closes the loop:
forward kinematics from `model.osim`, one rigid transform per segment, applied
to the BodyParts3D surfaces themselves.

Two outputs, both optional:

  --out FILE.mp4 / .png   a render of the real surfaces moving
  --emit FILE.json        an app-format trajectory --
                          `{centroids_m, frames:[{time_s, entities:{id:
                          {translation_m, rotation_matrix}}}]}` -- which is
                          what ihm/app/projection.py emits and what
                          app/src/state.js `bodyTransform` consumes, so the
                          viewer at 127.0.0.1:8765 can play it.

Geometry
--------
Let A be the similarity that carries OpenSim ground into the atlas frame and
T_ref[s] the segment's transform at the registered reference pose.  An entity
bound to s moves by

    M[s](t) = A . T_cur[s](t) . T_ref[s]^-1 . A^-1

which is exactly the identity at the reference pose, so the anatomy is never
deformed to make the model fit.  The scale in A cancels, so M is rigid and the
viewer's determinant check passes.  The render works in OpenSim ground metres
(`A^-1 . M[s]`), because that is the frame whose y = 0 is the floor the
contact model actually pushed against.

The skin is not bound rigidly.  It is the one surface the repository already
carries a continuous attachment for
(`data/derived/canonical/continuous_surface_binding.json.gz`, a
graph-regularised linear blend over the same 22 segments), and this script
uses those weights, so the integument stretches across a joint rather than
tearing at it.

What is rendered is stated, never implied
-----------------------------------------
`--select` and `--max-surfaces` choose a subset; the HUD and the printed
summary always name the number of surfaces drawn, out of how many are bound,
and what fraction of the bound surface area that is.

Usage
-----
    python scripts/render_anatomical_motion.py --select skin --out out.mp4
    python scripts/render_anatomical_motion.py --select bone,muscle --still 60
    python scripts/render_anatomical_motion.py --emit traj.json --select all
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BINDING = ROOT / "data/derived/anatomy-segment-binding/binding.json"
ANATOMY = ROOT / "data/derived/canonical/anatomy.json"
MODEL = ROOT / "data/models/engineering_stance_v1/model.osim"
MECHANICS = ROOT / "data/derived/canonical/mechanics.json"
SKIN_BINDING = ROOT / "data/derived/canonical/continuous_surface_binding.json.gz"
DEFAULT_TRAJECTORY = ROOT / "data/derived/gait-best/trajectory.json"

# Surface colour by anatomical system.  These are display conventions and
# carry no measurement; they exist so a vein does not read as an artery.
SYSTEM_COLOUR = {
    "skeletal": (0.86, 0.83, 0.75),
    "muscular": (0.62, 0.20, 0.19),
    "arterial": (0.72, 0.16, 0.16),
    "venous": (0.24, 0.32, 0.58),
    "nervous": (0.90, 0.87, 0.55),
    "connective": (0.74, 0.72, 0.62),
    "digestive": (0.72, 0.52, 0.36),
    "respiratory": (0.70, 0.45, 0.48),
    "cardiac": (0.68, 0.22, 0.24),
    "lymphatic": (0.52, 0.68, 0.58),
    "integumentary": (0.80, 0.63, 0.52),
    "urinary": (0.60, 0.56, 0.35),
    "reproductive": (0.66, 0.50, 0.52),
    "endocrine": (0.62, 0.55, 0.66),
    "sensory": (0.55, 0.62, 0.70),
}
SKIN_COLOUR = (0.80, 0.63, 0.52)
MUSCLE_COLD = np.array([0.44, 0.16, 0.16])
MUSCLE_HOT = np.array([1.00, 0.42, 0.24])

SELECTORS = {
    "bone": lambda r: r["role"] in ("rigid_bone", "cartilage"),
    "muscle": lambda r: r["role"] in ("muscle", "tendon"),
    "vessel": lambda r: r["role"] == "vascular",
    "nerve": lambda r: r["role"] == "nerve",
    "organ": lambda r: r["role"] in ("soft_organ", "fluid_cavity"),
    "connective": lambda r: r["role"] in ("connective_tissue", "ligament"),
    "all": lambda r: True,
}


def load_render_module():
    spec = importlib.util.spec_from_file_location(
        "_render_body_3d", ROOT / "scripts/render_body_3d.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_surface(path, orient=True):
    """Vertices and triangles of one atlas surface, wound outward.

    A handful of source surfaces are wound the other way; the rasteriser culls
    on winding and the shader lights on the normal, so an inverted one renders
    as a flat dark silhouette.  Orientation is decided by whether the face
    normals point away from the surface's own centre, which works for the open
    shells in the atlas as well as the closed ones.
    """
    payload = json.loads(gzip.decompress(Path(path).read_bytes()))
    v = np.asarray(payload["positions"], np.float64).reshape(-1, 3)
    f = np.asarray(payload["indices"], np.int64).reshape(-1, 3)
    if orient and len(f):
        a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
        normal = np.cross(b - a, c - a)
        outward = ((a + b + c) / 3.0) - v.mean(0)
        if float((normal * outward).sum()) < 0:
            f = f[:, ::-1]
    return v, f


def load_joint_frames(path):
    """(times, [{coordinate: value}]) from a gait trajectory or a stance trace."""
    payload = json.loads(Path(path).read_text())
    frames = payload["frames"] if isinstance(payload, dict) else payload
    times, coords, extra = [], [], []
    for n, frame in enumerate(frames):
        joints = frame.get("joints") or frame["coordinates"]
        coords.append({k: float(v["value"] if isinstance(v, dict) else v)
                       for k, v in joints.items()})
        times.append(float(frame.get("time_s", n * 0.01)))
        extra.append(frame)
    return np.asarray(times), coords, extra, payload


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    """Path as written in a receipt: repo-relative when it is inside the repo."""
    path = Path(path).resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def git_sha():
    sha = os.environ.get("IBM_GIT_SHA")
    if sha:
        return sha
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "git-unknown"


# --------------------------------------------------------------------------
# Segment motion
# --------------------------------------------------------------------------
class SegmentMotion:
    """M[s](t) = A . T_cur[s](t) . T_ref[s]^-1 . A^-1, in the atlas frame."""

    def __init__(self, binding, model):
        self.model = model
        self.segments = binding["segments"]
        self.A = np.asarray(binding["similarity_atlas_from_opensim_ground"], float)
        self.Ainv = np.linalg.inv(self.A)
        rest = model.forward(binding["reference_pose_rad"])
        self.rest_inverse = {s: np.linalg.inv(rest[s]) for s in self.segments}

    def atlas(self, q):
        current = self.model.forward(q)
        return {s: self.A @ current[s] @ self.rest_inverse[s] @ self.Ainv
                for s in self.segments}

    def render(self, q):
        """Atlas coordinates -> OpenSim ground metres, where y = 0 is the floor."""
        current = self.model.forward(q)
        return {s: current[s] @ self.rest_inverse[s] @ self.Ainv
                for s in self.segments}


def apply(transform, points):
    return points @ transform[:3, :3].T + transform[:3, 3]


# --------------------------------------------------------------------------
# What to draw
# --------------------------------------------------------------------------
def select_entities(rows, anatomy_by_id, spec, max_surfaces):
    tests = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if token == "skin":
            continue                       # handled by the continuous binding
        if token.startswith("system:"):
            name = token.split(":", 1)[1]
            tests.append(lambda r, n=name: r["system"] == n)
        elif token.startswith("role:"):
            name = token.split(":", 1)[1]
            tests.append(lambda r, n=name: r["role"] == n)
        elif token in SELECTORS:
            tests.append(SELECTORS[token])
        else:
            raise SystemExit(f"unknown selector {token!r}; known: "
                             f"{sorted(SELECTORS)} plus system:<x>, role:<x>, skin")
    if not tests:
        return [], 0.0, 0.0
    chosen = [i for i, r in rows.items()
              if r["segment"] and any(t(r) for t in tests)]
    area = {i: (anatomy_by_id[i].get("surface_area_m2") or 0.0) for i in chosen}
    total_area = sum(area.values())
    chosen.sort(key=lambda i: -area[i])
    if max_surfaces and len(chosen) > max_surfaces:
        chosen = chosen[:max_surfaces]
    kept_area = sum(area[i] for i in chosen)
    return chosen, kept_area, total_area


def muscle_excitation_map():
    """canonical entity id -> native muscle name, for the 80 driven muscles."""
    native = json.loads(MECHANICS.read_text())["native_muscles"]
    return {r["canonical_entity_id"]: r["source_name"] for r in native}



def write_motion_report(args, binding, rows, motion, model, coords, times,
                        index, by_id, geometry_path):
    """How far the anatomy moves, and what a rigid binding costs at each joint.

    Two quantities, both measured, neither a loss:

    * displacement -- the distance each entity's reference centroid travels
      across the clip.  The comparison that matters is the body's own
      canonical trajectory, which moves 6.35 mm at most over 30 s: it is
      perfusing and breathing, not moving.

    * joint opening -- at the reference pose the closest points across each
      joint are a fixed distance apart.  Bind the two sides to different rigid
      segments and that distance changes.  The maximum change is exactly the
      gap or interpenetration the rigid binding introduces, in metres, and it
      is a property of the binding and the motion, not of the anatomy.
    """
    from scipy.spatial import cKDTree

    ids = [i for i, r in rows.items() if r["segment"]]
    centroids = np.array([binding["centroids_m"][i] for i in ids])
    segment_of = np.array([motion.segments.index(rows[i]["segment"]) for i in ids])
    travel = np.zeros(len(ids))
    for k in index:
        matrices = np.stack([motion.atlas(coords[k])[s] for s in motion.segments])
        moved = (np.einsum("vij,vj->vi", matrices[segment_of][:, :3, :3], centroids)
                 + matrices[segment_of][:, :3, 3])
        travel = np.maximum(travel, np.linalg.norm(moved - centroids, axis=1))

    groups = binding["segment_named_bones"]
    clouds = {}
    for segment, bones in groups.items():
        points = [json.loads(gzip.decompress(Path(geometry_path[b]).read_bytes()))
                  ["positions"] for b in bones]
        clouds[segment] = np.concatenate(
            [np.asarray(x, float).reshape(-1, 3) for x in points])

    joints = []
    for joint in model.joints:
        parent, child = joint["parent"], joint["child"]
        if parent not in clouds or child not in clouds:
            continue
        a = clouds[parent]
        b = clouds[child]
        distance, index_b = cKDTree(b).query(a)
        pick = np.argsort(distance)[:64]
        pa, pb = a[pick], b[index_b[pick]]
        rest = np.linalg.norm(pa - pb, axis=1)
        worst = 0.0
        for k in index:
            transforms = motion.atlas(coords[k])
            qa = apply(transforms[parent], pa)
            qb = apply(transforms[child], pb)
            worst = max(worst, float(np.abs(np.linalg.norm(qa - qb, axis=1)
                                            - rest).max()))
        joints.append({"joint": joint["name"], "parent": parent, "child": child,
                       "rest_gap_m": float(rest.mean()),
                       "max_opening_m": worst})

    openings = np.array([j["max_opening_m"] for j in joints])
    report = {
        "schema": "ihm.anatomical-motion-report.v1",
        "trajectory": {"path": relative(args.trajectory),
                       "sha256": digest(args.trajectory),
                       "frames": len(index),
                       "duration_s": float(times[index[-1]] - times[index[0]])},
        "binding": {"path": relative(args.binding),
                    "sha256": digest(args.binding)},
        "entities": len(ids),
        "centroid_travel_m": {
            "max": float(travel.max()),
            "mean": float(travel.mean()),
            "median": float(np.median(travel)),
            "p90": float(np.percentile(travel, 90)),
            "over_10mm": int((travel > 0.010).sum()),
            "over_100mm": int((travel > 0.100).sum()),
        },
        "baseline_canonical_trajectory_max_travel_m": 0.00635,
        "baseline_note": "The body's own canonical run moves 6.35 mm at most "
                         "across 30 s (docs/DISCONNECTS.md item 1). It is "
                         "perfusing and breathing, not moving. That is the "
                         "number this render has to beat to mean anything.",
        "farthest_travelling": [
            {"id": ids[i], "name": rows[ids[i]]["name"],
             "segment": rows[ids[i]]["segment"], "travel_m": float(travel[i])}
            for i in np.argsort(-travel)[:10]],
        "joint_opening_m": {
            "max": float(openings.max()), "median": float(np.median(openings)),
            "definition": "largest change, over the clip, in the distance "
                          "between the 64 closest rest-pose surface point "
                          "pairs across the joint; a rigid binding cannot "
                          "keep this at zero",
            "per_joint": sorted(joints, key=lambda j: -j["max_opening_m"]),
        },
        "git_sha": git_sha(),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"motion: the anatomy travels up to "
          f"{travel.max()*1000:.0f} mm (median {np.median(travel)*1000:.0f} mm); "
          f"the canonical run's maximum was 6.35 mm over 30 s")
    print(f"joint opening: max {openings.max()*1000:.0f} mm at "
          f"{max(joints, key=lambda j: j['max_opening_m'])['joint']}, "
          f"median {np.median(openings)*1000:.0f} mm")
    print(f"wrote {args.report}")


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trajectory", type=Path, default=DEFAULT_TRAJECTORY)
    ap.add_argument("--binding", type=Path, default=BINDING)
    ap.add_argument("--select", default="bone,muscle",
                    help="comma list of bone, muscle, vessel, nerve, organ, "
                         "connective, skin, all, system:<x>, role:<x>")
    ap.add_argument("--max-surfaces", type=int, default=0,
                    help="keep only the N largest selected surfaces (0 = all)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--emit", type=Path, default=None,
                    help="write an app-format trajectory for the viewer")
    ap.add_argument("--emit-precision", type=int, default=6)
    ap.add_argument("--report", type=Path, default=None,
                    help="write how far the anatomy actually moved and how "
                         "badly the rigid binding tears at each joint")
    ap.add_argument("--still", type=int, default=None)
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=900)
    ap.add_argument("--fps", type=int, default=25)
    ap.add_argument("--ssaa", type=int, default=2)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--spin", type=float, default=0.0,
                    help="degrees of camera orbit across the clip")
    ap.add_argument("--title", default="IHM-1 - the real anatomy, driven by the controller")
    args = ap.parse_args()
    started = time.time()

    if args.out is None and args.emit is None and args.report is None:
        raise SystemExit("nothing to do: pass --out, --emit and/or --report")

    rb = load_render_module()
    binding = json.loads(args.binding.read_text())
    rows = binding["entities"]
    centroids = binding["centroids_m"]
    anatomy = json.loads(ANATOMY.read_text())
    by_id = {e["id"]: e for e in anatomy["entities"]}
    geometry_path = {e["id"]: ROOT / e["reference_geometry"]["path"]
                     for e in anatomy["entities"] if e.get("reference_geometry")}

    model = rb.OsimModel(MODEL)
    motion = SegmentMotion(binding, model)
    times, coords, raw_frames, payload = load_joint_frames(args.trajectory)
    # A trajectory that names only the coordinates it drives leaves the rest
    # undefined.  Held at zero they would drop the pelvis to the ground
    # origin and hang the body a metre below the floor, so the unspecified
    # ones are held at the registered reference pose instead.
    reference = binding["reference_pose_rad"]
    absent = sorted(set(reference) - set(coords[0]))
    coords = [{**reference, **frame} for frame in coords]
    if absent:
        print(f"trajectory does not drive {len(absent)} coordinates "
              f"({', '.join(absent[:6])}{'...' if len(absent) > 6 else ''}); "
              f"held at the reference pose")
    index = list(range(0, len(times), max(1, args.stride)))
    if args.max_frames:
        index = index[:args.max_frames]
    bound = sum(1 for r in rows.values() if r["segment"])
    print(f"binding: {bound}/{len(rows)} entities bound, "
          f"registration RMS "
          f"{binding['registration']['segment_centroid_residual_rms_m']*1000:.1f} mm")
    print(f"trajectory: {relative(args.trajectory)}  "
          f"{len(times)} frames, {times[-1]-times[0]:.2f} s")

    # ---------------------------------------------------------------- report
    if args.report:
        write_motion_report(args, binding, rows, motion, model, coords, times,
                            index, by_id, geometry_path)

    # ---------------------------------------------------------------- emit
    if args.emit:
        digits = args.emit_precision
        out_frames = []
        for k in index:
            transforms = motion.atlas(coords[k])
            entities = {}
            for ident, row in rows.items():
                segment = row["segment"]
                if segment is None:
                    continue
                m = transforms[segment]
                c = np.asarray(centroids[ident], float)
                entities[ident] = {
                    "translation_m": [round(float(v), digits)
                                      for v in (m[:3, :3] @ c + m[:3, 3] - c)],
                    "rotation_matrix": [[round(float(v), 9) for v in r]
                                        for r in m[:3, :3]],
                }
            out_frames.append({"time_s": float(times[k]), "entities": entities})
        emitted = {
            "schema": "ihm.body-trajectory.v1",
            "model_id": binding["model_id"],
            "frame": binding["frame"],
            "centroids_m": {i: centroids[i] for i in rows if rows[i]["segment"]},
            "frames": out_frames,
            "projection": {
                "kind": "segment_rigid_binding",
                "basis": "Forward kinematics of a controller joint trajectory "
                         "through the OpenSim model, applied to the anatomical "
                         "surfaces through data/derived/anatomy-segment-binding. "
                         "Every entity bound to one segment moves exactly "
                         "rigidly with it; deformation_gradient is identity "
                         "everywhere because none is computed.",
                "entities_moved": len(out_frames[0]["entities"]) if out_frames else 0,
                "entities_omitted": len(rows) - bound,
                "omitted_reason": "whole-body surfaces with no single rigid "
                                  "segment; see binding report `excluded`",
                "skin_note": "The skin entities are NOT in this trajectory. A "
                             "rigid segment cannot carry them; the continuous "
                             "linear-blend attachment is the correct vehicle "
                             "and this schema has no field for it.",
                "source_trajectory": {
                    "path": relative(args.trajectory),
                    "sha256": digest(args.trajectory),
                    "schema": payload.get("schema") if isinstance(payload, dict) else None,
                },
                "binding": {"path": relative(args.binding),
                            "sha256": digest(args.binding)},
                "git_sha": git_sha(),
            },
        }
        args.emit.parent.mkdir(parents=True, exist_ok=True)
        args.emit.write_text(json.dumps(emitted) + "\n")
        size = args.emit.stat().st_size / 1e6
        print(f"wrote {args.emit}  {len(out_frames)} frames, "
              f"{len(emitted['centroids_m'])} entities, {size:.1f} MB")
        if args.out is None:
            return

    if args.out is None:
        return

    # ---------------------------------------------------------------- render
    want_skin = "skin" in [t.strip() for t in args.select.split(",")]
    chosen, kept_area, total_area = select_entities(rows, by_id, args.select,
                                                    args.max_surfaces)
    print(f"selected {len(chosen)} surfaces "
          f"({kept_area:.3f} of {total_area:.3f} m2 selected area)")

    excite = muscle_excitation_map()
    meshes = []            # (segment, V_atlas, F, base colour, entity id)
    offset = 0
    verts, faces, colours, seg_index, entity_of = [], [], [], [], []
    for ident in chosen:
        v, f = read_surface(geometry_path[ident])
        row = rows[ident]
        colour = np.array(SYSTEM_COLOUR.get(row["system"], (0.6, 0.6, 0.6)), np.float32)
        verts.append(v)
        faces.append(f + offset)
        colours.append(np.repeat(colour[None, :], len(v), 0))
        seg_index.append(np.full(len(v), motion.segments.index(row["segment"]), np.int32))
        entity_of.append((ident, offset, offset + len(v)))
        offset += len(v)
    if not verts:
        if not want_skin:
            raise SystemExit("nothing selected to render")
        V0 = np.zeros((0, 3))
        F = np.zeros((0, 3), np.int32)
        VC0 = np.zeros((0, 3), np.float32)
        SEG = np.zeros(0, np.int32)
    else:
        V0 = np.concatenate(verts).astype(np.float64)
        F = np.concatenate(faces).astype(np.int32)
        VC0 = np.concatenate(colours).astype(np.float32)
        SEG = np.concatenate(seg_index)
    print(f"rigid surfaces: {len(V0)} vertices, {len(F)} triangles")

    skin = None
    if want_skin:
        payload_skin = json.loads(gzip.decompress(SKIN_BINDING.read_bytes()))
        order = [s["id"] for s in payload_skin["segments"]]
        weights = np.asarray(payload_skin["weights"], np.float64)
        positions = np.asarray(payload_skin["reference_positions_m"], np.float64)
        skin_id = payload_skin["surface_entity_ids"][0]
        _, skin_faces = read_surface(geometry_path[skin_id])
        skin = {"positions": positions, "faces": skin_faces.astype(np.int32),
                "weights": weights,
                "columns": [motion.segments.index(s) for s in order],
                "id": skin_id}
        print(f"skin: {len(positions)} vertices, {len(skin_faces)} triangles, "
              f"continuous linear blend over {weights.shape[1]} segments")
        # Topology never changes, so the draw list is built once.
        F = np.concatenate([F, skin["faces"] + len(V0)]).astype(np.int32)
        skin["colour"] = np.repeat(np.array(SKIN_COLOUR, np.float32)[None, :],
                                   len(positions), 0)

    ss = max(1, args.ssaa)
    width, height = args.width, args.height
    raster_w, raster_h = width * ss, height * ss
    aspect = width / height
    frames_to_do = [args.still] if args.still is not None else index

    writer = None
    if args.still is None:
        import imageio_ffmpeg
        cmd = [imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error",
               "-f", "rawvideo", "-pix_fmt", "rgb24",
               "-s", f"{width}x{height}", "-r", str(args.fps), "-i", "-",
               "-an", "-c:v", "libx264", "-preset", "slow", "-crf", "18",
               "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(args.out)]
        args.out.parent.mkdir(parents=True, exist_ok=True)
        writer = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    excitation_hi = 0.15
    values = [v for fr in raw_frames[::5] for v in fr.get("motor_excitations", {}).values()]
    if values:
        excitation_hi = max(float(np.percentile(values, 99.0)), 1e-3)

    for n, k in enumerate(frames_to_do):
        u = n / max(len(frames_to_do) - 1, 1)
        transforms = motion.render(coords[k])
        matrices = np.stack([transforms[s] for s in motion.segments])

        V = np.einsum("vij,vj->vi", matrices[SEG][:, :3, :3], V0) + matrices[SEG][:, :3, 3]
        VC = VC0.copy()
        ex = raw_frames[k].get("motor_excitations", {})
        if ex:
            for ident, lo, hi in entity_of:
                name = excite.get(ident)
                if name is None or name not in ex:
                    continue
                t = min(max(float(ex[name]) / excitation_hi, 0.0), 1.0)
                VC[lo:hi] = (MUSCLE_COLD + (MUSCLE_HOT - MUSCLE_COLD) * t).astype(np.float32)

        if skin is not None:
            blended = np.zeros_like(skin["positions"])
            for column, segment_index in enumerate(skin["columns"]):
                w = skin["weights"][:, column]
                active = w > 1e-8
                if not active.any():
                    continue
                m = matrices[segment_index]
                blended[active] += (w[active, None]
                                    * (skin["positions"][active] @ m[:3, :3].T + m[:3, 3]))
            V = np.concatenate([V, blended])
            VC = np.concatenate([VC, skin["colour"]])
        F_draw = F

        V = V.astype(np.float32)
        VN = rb.vertex_normals(V.astype(np.float64), F_draw).astype(np.float32)

        if n == 0:
            # One framing for the whole clip, taken from the first frame, so
            # that what moves on screen is the body and not the camera.
            lo, hi = V.min(0), V.max(0)
            centre = np.array([float((lo[0] + hi[0]) / 2), 0.0,
                               float((lo[2] + hi[2]) / 2)])
            target = np.array([centre[0], float((lo[1] + hi[1]) / 2), centre[2]])
            # 27 deg vertical fov: the body fills ~72% of the frame height,
            # which keeps it clear of the title and caption bars.
            dist = float(max(hi[1] - lo[1], 1.0)) * 2.90
        az = math.radians(-24.0 + args.spin * u)
        el = math.radians(8.0)
        eye = target + dist * np.array([math.cos(el) * math.cos(az), math.sin(el),
                                        math.cos(el) * math.sin(az)])
        eye[1] = max(eye[1], 0.35)
        cam = rb.Camera(eye, target, fov_deg=27.0, aspect=aspect)

        extent = (centre[0] - 3.0, centre[0] + 3.0, centre[2] - 3.0, centre[2] + 3.0)
        shadow = rb.shadow_map(V.astype(np.float64), F_draw, rb.KEY_DIR, extent, res=320)
        ground, depth = rb.render_ground(cam, width, height, shadow, extent,
                                         np.array([target[0], 1.0, target[2]]),
                                         [(centre[0], centre[2], (0.55, 0.62, 0.72))])
        # Two-sided shading.  A few atlas shells are wound the other way and
        # `read_surface` cannot always tell; a normal pointing away from the
        # eye would light them as flat silhouettes.  Flipping toward the eye
        # is a no-op on every correctly wound closed surface.
        towards = cam.eye.astype(np.float32)[None, :] - V
        VN = np.where(((VN * towards).sum(1) < 0)[:, None], -VN, VN)
        colour, coverage = rb.shade_body(cam, V, F_draw, VN, VC, raster_w, raster_h,
                                         depth.reshape(-1), ss)
        image = ground * (1.0 - coverage[..., None]) + colour
        image = rb.tonemap(rb.vignette(image))
        image = overlay(rb, image, width, height, args.title, binding, rows,
                        len(chosen), kept_area, total_area, skin, float(times[k]),
                        raw_frames[k], (n + 1) / len(frames_to_do), args.trajectory,
                        absent)

        buffer = (np.clip(image, 0, 1) * 255).astype(np.uint8)
        if writer is None:
            from PIL import Image
            path = args.out if args.out.suffix == ".png" else args.out.with_suffix(".png")
            path.parent.mkdir(parents=True, exist_ok=True)
            Image.fromarray(buffer).save(path)
            print(f"wrote still {path}  ({time.time()-started:.0f}s)")
            return
        writer.stdin.write(buffer.tobytes())
        if n % 5 == 0 or n == len(frames_to_do) - 1:
            elapsed = time.time() - started
            rate = (n + 1) / max(elapsed, 1e-6)
            print(f"  frame {n+1:4d}/{len(frames_to_do)}  {rate:4.2f} fps  "
                  f"eta {(len(frames_to_do)-n-1)/max(rate,1e-6)/60:4.1f} min", flush=True)

    writer.stdin.close()
    if writer.wait() != 0:
        print("ffmpeg failed", file=sys.stderr)
        sys.exit(1)
    print(f"wrote {args.out}  {width}x{height}  {len(frames_to_do)} frames  "
          f"{args.out.stat().st_size/1e6:.1f} MB  ({time.time()-started:.0f}s)")


def overlay(rb, image, width, height, title, binding, rows, drawn, kept_area,
            total_area, skin, t, frame, progress, trajectory, absent=()):
    """Caption that says exactly what is on screen, including what is not."""
    from PIL import Image, ImageDraw
    pil = Image.fromarray((np.clip(image, 0, 1) * 255).astype(np.uint8))
    draw = ImageDraw.Draw(pil, "RGBA")
    big = rb._font(int(height * 0.026), mono=False, bold=True)
    small = rb._font(int(height * 0.0155))
    tiny = rb._font(int(height * 0.0135))

    draw.rectangle([0, 0, width, int(height * 0.075)], fill=(0, 0, 0, 130))
    draw.rectangle([0, int(height * 0.865), width, height], fill=(0, 0, 0, 150))
    pad = int(width * 0.018)
    draw.text((pad, int(height * 0.018)), title, font=big, fill=(238, 242, 248))

    bound = sum(1 for r in rows.values() if r["segment"])
    fit = binding["registration"]
    lines = [
        f"t = {t:6.2f} s   " + ("FALLEN   " if frame.get("fallen") else "")
        + "   ".join(
            f"{name} {value:.2f}" if isinstance(value, (int, float))
            else f"{name} {value}"
            for name, value in (("phase", frame.get("phase")),
                                ("swing", frame.get("swing_side")),
                                )
            if value is not None and not isinstance(value, bool)),
        (f"{drawn} anatomical surfaces drawn of {bound} bound "
         f"({kept_area/max(total_area,1e-9)*100:.0f}% of the selected surface "
         f"area)" if drawn else "skin surface only; the other "
         f"{bound} bound surfaces are inside it and not drawn")
        + ("  +  skin, continuous linear blend" if skin and drawn else ""),
        f"22 OpenSim segments  ·  registration residual "
        f"{fit['segment_centroid_residual_rms_m']*1000:.0f} mm RMS  ·  "
        f"cross-specimen, not measurement error",
        (f"{trajectory.name}: colour on the 80 driven muscles is their own "
         f"excitation; the other muscle surfaces are not driven" if drawn
         else f"{trajectory.name}: skin position is a linear blend of the 22 "
              f"segment transforms, solved once on the rest topology"),
    ]
    if any(c.startswith("pelvis_t") for c in absent):
        lines.append("this trajectory carries no pelvis translation, so the "
                     "gait is played in place")
    y = int(height * 0.874)
    for n, line in enumerate(lines):
        draw.text((pad, y + n * int(height * 0.0215)), line,
                  font=small if n < 2 else tiny,
                  fill=(226, 232, 240) if n < 2 else (150, 162, 178))

    bar = int(height * 0.004)
    draw.rectangle([0, height - bar, int(width * progress), height],
                   fill=(120, 190, 230, 220))
    return np.asarray(pil).astype(np.float32) / 255.0


if __name__ == "__main__":
    main()
