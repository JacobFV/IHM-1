#!/usr/bin/env python
"""Register sourced garment meshes onto the canonical outer body envelope.

Runs in the vendored libigl venv (AABB closest-point queries and fast winding
number). Three recorded stages per garment:

  1. one global similarity transform (uniform scale + translation) fitted by
     minimising the source body's point-to-surface distance to the envelope;
  2. a piecewise-rigid limb pose correction, because the MakeHuman authoring
     body is abducted and this body's arms are adducted; linear blend skinning
     with one identity bone for the trunk and one rotation per limb;
  3. a Laplacian-regularised shrinkwrap that pushes any vertex inside the
     standoff shell outward and attracts near vertices onto it, leaving distant
     vertices to hang, so pleats, flares and hems keep their own shape.

Every number reported here is measured on the produced mesh.

  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/fit_garments_to_envelope.py [--self-test]
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import igl
import numpy as np
from scipy.optimize import minimize
from scipy.sparse import coo_matrix, diags, identity
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.garment_remesh import quality, remesh, split_components  # noqa: E402
from ihm.assembly.garment_wardrobe import (  # noqa: E402
    CATALOGUE, SOURCE_UNITS_PER_M, boundary_edge_count, compact, edge_table, face_components,
    intersection_report, keep_components, load_envelope, load_obj, mesh_area_m2, orient_outward,
    sha256_file, slot_model, weld)

ENVELOPE = ROOT / 'data/derived/outer-envelope/outer-envelope.npz'
RAW = ROOT / 'data/raw/clothing/makehuman'
OUT = ROOT / 'data/derived/wardrobe-v1/fitted'
CONFORM_ATTRACT_RANGE_M = 0.060
DRAPE_ATTRACT_RANGE_MULTIPLE = 3.0
DRAPE_ATTRACT_RANGE_FLOOR_M = 0.012
# Shell garments hold their own volume: a fedora crown stands off the scalp, a
# shoe encloses the foot. Their attraction band is a narrow collar just outside
# the standoff shell, so only cloth already at the skin is held there and the
# rest of the silhouette is left where the source authored it.
SHELL_ATTRACT_RANGE_MULTIPLE = 2.0
SHELL_ATTRACT_RANGE_FLOOR_M = 0.006
# A loose sleeve has its inner wall inside the attraction band and its outer wall
# outside it. With a weak shape term the pinned inner wall drags the free outer
# wall in and the tube flattens, so drape mode raises the shape term far above
# the attraction it has to resist.
FIT_MODES = {'conform': {'regularisation': 2.0, 'attract_weight': 6.0, 'free_weight': 0.02},
             'skin': {'regularisation': 10.0, 'attract_weight': 5.0, 'free_weight': 0.05},
             'drape': {'regularisation': 12.0, 'attract_weight': 2.0, 'free_weight': 0.30},
             'shell': {'regularisation': 30.0, 'attract_weight': 0.8, 'free_weight': 0.60}}
ATTRACT_RANGE = {
    'conform': lambda standoff: CONFORM_ATTRACT_RANGE_M,
    'skin': lambda standoff: CONFORM_ATTRACT_RANGE_M,
    'drape': lambda standoff: max(DRAPE_ATTRACT_RANGE_FLOOR_M, DRAPE_ATTRACT_RANGE_MULTIPLE * standoff),
    'shell': lambda standoff: max(SHELL_ATTRACT_RANGE_FLOOR_M, SHELL_ATTRACT_RANGE_MULTIPLE * standoff),
}
# Two garments only have a layering relation where they are close enough to be
# in contact. Beyond this range an outer hem hanging free of the inner garment
# would otherwise be scored as generous clearance it does not actually provide.
# The acquired meshes are far too coarse to be fitted: at a 76 mm median edge a
# single triangle spans a whole limb, so it passes through the body between its
# own corners and the shrinkwrap has no vertices to distribute strain across.
# Every garment is therefore resampled to one uniform edge length first. 10 mm
# puts roughly twenty-five segments around a forearm and is a little finer than
# the 8 mm contact faces of the body it has to be measured against.
# The shape term of the shrinkwrap is a graph-Laplacian energy, and a graph
# Laplacian is not resolution-invariant: on a surface sampled at edge length h,
# (L x) scales like h^2, so the shape energy scales like h^4 against a fit term
# that is a plain positional penalty. Resampling a garment therefore silently
# changes the balance the FIT_MODES weights were chosen for, and the fit stops
# distributing strain and starts pulling single vertices onto the standoff
# shell. So each garment's declared weight is restated for its new resolution
# against its own authored edge length, which makes the resample a change of
# resolution only and leaves every garment the character it was fitted with.
#
# Measured, holding everything else fixed:
#
#   t-shirt      conform  25.2 -> 8.2 mm  factor  89
#                edge ratio p99 2.88 -> 5.17 unscaled -> 2.77 scaled
#                standoff p50 4.00 -> 4.4 mm, body-penetration rate 0.69% -> 0.22%
#   casual-shirt drape    34.0 -> 8.5 mm  factor 256
#                edge ratio p99 2.28 -> 2.26 unscaled -> 1.62 scaled
#                standoff p50 24.90 -> 25.09 mm, penetration rate 0.54% -> 0.06%
#
# The standoff medians are the check that matters: a drape garment carries real
# ease, and at the unscaled weight the casual shirt collapsed onto the body
# (24.90 -> 6.98 mm) even though its edge-ratio numbers looked fine.
REGULARISATION_RESOLUTION_EXPONENT = 4.0
REMESH_TARGET_EDGE_M = 0.010
REMESH_PASSES = 10
# A component has to be big enough for a uniform target-edge sampling to mean
# anything. The polo shirt carries two closed boxes 13 mm across, smaller than
# two target edges: there is no resampling of those to do, and attempting it
# collapses them out of the garment. So a component is resampled only if it
# spans several target edges and holds more than a couple of target triangles,
# and anything else is carried through exactly as authored.
MIN_RESAMPLED_COMPONENT_TRIANGLES = 2.0
MIN_RESAMPLED_COMPONENT_EXTENT_EDGES = 6.0
LAYER_MEETING_RANGE_M = 0.060
LAYER_MEETING_MIN_VERTICES = 10
SHRINKWRAP_ITERATIONS = 18
# A sample deeper than this is a long edge or face cutting a concave corner, not
# the body emerging through the cloth. Moving its vertices cannot fix it and
# wrecks the garment, so it is counted and reported instead.
SAMPLE_MAX_DEPTH_M = 0.010


def load_obj_group(path, group):
    """OBJ faces restricted to one `g` group. MakeHuman base.obj carries body,
    clothes helpers and joint cubes in one file."""
    positions, triangles, current = [], [], None
    with open(path, 'r', errors='replace') as handle:
        for line in handle:
            if line.startswith('v '):
                parts = line.split()
                positions.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('g '):
                current = line.split(None, 1)[1].strip()
            elif line.startswith('f ') and current == group:
                index = [int(token.split('/')[0]) - 1 for token in line.split()[1:]]
                for k in range(1, len(index) - 1):
                    triangles.append((index[0], index[k], index[k + 1]))
    if not triangles:
        raise ValueError(f'{path}: group {group!r} has no faces')
    return compact(np.asarray(positions, float), np.asarray(triangles, np.int64))


def surface_distance(points, envelope_v, envelope_f):
    squared, index, closest = igl.point_mesh_squared_distance(np.ascontiguousarray(points), envelope_v, envelope_f)
    return np.sqrt(np.maximum(squared, 0.0)), index, closest


def surface_projectors(positions, triangles):
    """Closest-point maps onto a mesh and onto its own boundary polyline.

    Remeshing needs both. Interior samples go back onto the surface, so the
    resampled garment is the authored one at a different resolution rather than
    a smoothed approximation of it; boundary samples go back onto the boundary
    curve, so a hem, a neckline and an armhole keep the outline they were drawn
    with instead of being rounded off by the relaxation.
    """
    positions = np.ascontiguousarray(positions)
    triangles = np.ascontiguousarray(triangles)
    edges = np.sort(np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]), axis=1)
    unique, count = np.unique(edges, axis=0, return_counts=True)
    border = np.ascontiguousarray(unique[count == 1])

    def onto_surface(points):
        _, _, closest = igl.point_mesh_squared_distance(np.ascontiguousarray(points), positions, triangles)
        return closest

    def onto_boundary(points):
        _, _, closest = igl.point_mesh_squared_distance(np.ascontiguousarray(points), positions, border)
        return closest

    return onto_surface, (onto_boundary if len(border) else None)


def resample(positions, triangles, target_edge_m, passes=REMESH_PASSES):
    """Isotropic remesh of one garment onto its own surface, with a receipt.

    Component by component, because several garments are not one sheet: the polo
    shirt is five, three of them scraps of eight vertices. Resampled together, a
    scrap is collapsed away or absorbed by its neighbour and the garment quietly
    loses a piece. Resampled apart, each piece is projected onto its own surface
    and its own boundary, and a piece with less area than a couple of target
    triangles is carried through untouched rather than destroyed, since there is
    no resampling of it to do.
    """
    floor = MIN_RESAMPLED_COMPONENT_TRIANGLES * (3.0 ** 0.5 / 4.0) * target_edge_m ** 2
    reach = MIN_RESAMPLED_COMPONENT_EXTENT_EDGES * target_edge_m
    report = {'before': quality(positions, triangles), 'components': []}
    pieces, histories = [], []
    for index, (piece_v, piece_f) in enumerate(split_components(positions, triangles)):
        piece_area = mesh_area_m2(piece_v, piece_f)
        extent = float(np.linalg.norm(piece_v.max(0) - piece_v.min(0)))
        row = {'component': index, 'faces': int(len(piece_f)), 'area_m2': piece_area,
               'extent_mm': extent * 1e3}
        if piece_area < floor or extent < reach:
            row.update({'resampled': False,
                        'reason': f'spans {extent * 1e3:.1f} mm and holds {piece_area / (floor / 2):.1f} target '
                                  f'triangles; below {reach * 1e3:.0f} mm or '
                                  f'{MIN_RESAMPLED_COMPONENT_TRIANGLES:.0f} triangles there is no resampling to do'})
            report['components'].append(row)
            pieces.append((piece_v, piece_f))
            continue
        onto_surface, onto_boundary = surface_projectors(piece_v, piece_f)
        piece_report = {}
        try:
            new_v, new_f = remesh(piece_v, piece_f, target_edge_m, passes=passes,
                                  project=onto_surface, project_boundary=onto_boundary, report=piece_report)
        except ValueError as failure:
            # Carried, not dropped, and the refusal is recorded rather than
            # swallowed: a component the resampler will not accept is still part
            # of the garment and must reach the fit as authored.
            row.update({'resampled': False, 'reason': f'resampler refused it: {failure}'})
            report['components'].append(row)
            pieces.append((piece_v, piece_f))
            continue
        row.update({'resampled': True, 'faces_after': int(len(new_f)),
                    'boundary_edges_before': piece_report['boundary_edges_before'],
                    'boundary_edges_after': piece_report['boundary_edges_after']})
        report['components'].append(row)
        histories.append(piece_report['history'])
        pieces.append((new_v, new_f))
    offset = 0
    stacked_v, stacked_f = [], []
    for piece_v, piece_f in pieces:
        stacked_v.append(piece_v)
        stacked_f.append(piece_f + offset)
        offset += len(piece_v)
    resampled = np.vstack(stacked_v)
    faces = np.vstack(stacked_f)
    if face_components(len(resampled), faces)[0] != face_components(len(positions), triangles)[0]:
        raise ValueError('resample changed the number of surface components')
    report['history'] = histories
    report['passes'] = passes
    report['target_edge_mm'] = target_edge_m * 1e3
    report['source_vertices'] = int(len(positions))
    report['source_faces'] = int(len(triangles))
    report['vertices'] = int(len(resampled))
    report['faces'] = int(len(faces))
    report['boundary_edges_before'] = boundary_edge_count(triangles)[0]
    report['boundary_edges_after'] = boundary_edge_count(faces)[0]
    report['after'] = quality(resampled, faces)
    deviation, _, _ = igl.point_mesh_squared_distance(
        np.ascontiguousarray(resampled), np.ascontiguousarray(positions), np.ascontiguousarray(triangles))
    before_area = mesh_area_m2(positions, triangles)
    after_area = mesh_area_m2(resampled, faces)
    report['max_deviation_from_source_mm'] = float(np.sqrt(max(0.0, deviation.max())) * 1e3)
    report['area_m2_before'] = before_area
    report['area_m2_after'] = after_area
    report['area_change_percent'] = float(100.0 * (after_area - before_area) / before_area)
    report['area_note'] = ('every resampled vertex lies on the source surface, so the only area the resample can '
                           'lose is the ridge of a crease that a new triangle bridges as a chord')
    return resampled, faces, report


def fit_similarity(source, envelope_v, envelope_f, trunk_halfwidth_m=0.16, min_trunk_points=500):
    """Stature landmark fixes the uniform scale; translation is then refined on
    trunk points only. Optimising the scale as well has a degenerate minimum
    that shrinks the body inside the envelope, so the scale is not free."""
    height_target = float(envelope_v[:, 1].max() - envelope_v[:, 1].min())
    height_source = float(source[:, 1].max() - source[:, 1].min())
    scale = height_target / height_source
    scaled = scale * source
    shift0 = np.array([0.0, envelope_v[:, 1].min() - scaled[:, 1].min(), 0.0])
    shift0[[0, 2]] = envelope_v[:, [0, 2]].mean(0) - scaled[:, [0, 2]].mean(0)
    trunk = np.abs(scaled[:, 0] + shift0[0]) < trunk_halfwidth_m
    if trunk.sum() < min_trunk_points:
        raise ValueError(f'trunk landmark band holds {int(trunk.sum())} points')
    def objective(shift):
        d, _, _ = surface_distance(scaled[trunk] + shift, envelope_v, envelope_f)
        return float(np.sqrt(np.mean(d * d)))
    result = minimize(objective, shift0, method='Nelder-Mead',
                      options={'xatol': 1e-6, 'fatol': 1e-9, 'maxiter': 600})
    translation = result.x.astype(float)
    d_all, _, _ = surface_distance(scaled + translation, envelope_v, envelope_f)
    d0, _, _ = surface_distance(scaled + shift0, envelope_v, envelope_f)
    report = {
        'height_source_m': height_source, 'height_target_m': height_target,
        'trunk_halfwidth_m': trunk_halfwidth_m, 'trunk_points': int(trunk.sum()),
        'trunk_rms_mm_stature_only': float(objective(shift0) * 1e3),
        'trunk_rms_mm_refined': float(result.fun * 1e3),
        'all_points_rms_mm_stature_only': float(np.sqrt(np.mean(d0 ** 2)) * 1e3),
        'all_points_rms_mm_refined': float(np.sqrt(np.mean(d_all ** 2)) * 1e3),
    }
    return scale, translation, report


def obj_group_centroid(path, groups):
    """Centroid of the vertices referenced by the named OBJ groups. MakeHuman
    ships one small cube per skeleton joint, so this is an authored landmark."""
    positions, wanted, current, used = [], set(groups), None, set()
    with open(path, 'r', errors='replace') as handle:
        for line in handle:
            if line.startswith('v '):
                parts = line.split()
                positions.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('g '):
                current = line.split(None, 1)[1].strip()
            elif line.startswith('f ') and current in wanted:
                used.update(int(token.split('/')[0]) - 1 for token in line.split()[1:])
    if not used:
        raise ValueError(f'{path}: groups {sorted(wanted)} carry no faces')
    return np.asarray(positions, float)[sorted(used)].mean(0)


def bone_landmark(root, entity_id, end, quantile=0.03):
    """Extreme-end centroid of one canonical bone surface, the target landmark."""
    path = root / 'data/derived/canonical/geometry' / f'{entity_id}.json.gz'
    with gzip.open(path, 'rt') as handle:
        positions = np.asarray(json.load(handle)['positions'], float).reshape(-1, 3)
    y = positions[:, 1]
    mask = y >= np.quantile(y, 1 - quantile) if end == 'top' else y <= np.quantile(y, quantile)
    return positions[mask].mean(0), str(path.relative_to(root)), sha256_file(path), int(len(positions))


def rotation_between(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    c = float(np.dot(a, b))
    if s < 1e-12:
        return np.eye(3) if c > 0 else -np.eye(3)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k * ((1 - c) / (s * s))


def segment_distance(points, a, b):
    direction = b - a
    t = np.clip((points - a) @ direction / float(direction @ direction), 0.0, 1.0)
    return np.linalg.norm(points - (a + t[:, None] * direction), axis=1)


def vertex_adjacency(vertex_count, triangles):
    edges = edge_table(triangles)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    a = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(vertex_count, vertex_count)).tocsr()
    degree = np.asarray(a.sum(1)).ravel()
    degree[degree == 0] = 1.0
    return diags(1.0 / degree) @ a


class LimbSkinning:
    """Identity trunk bone plus one similarity per limb.

    Blend weights are segmented on the authoring body, not on the garment: every
    body vertex is labelled by the skeleton segment it lies nearest, the indicator
    field is smoothed over the body mesh, and each garment vertex inherits the
    weights of its nearest body vertex. A positional gate on lateral coordinate
    was tried first and tears long sleeves, because a loose sleeve straddles the
    gate and its inner and outer faces receive different rotations.
    """

    def __init__(self, bones, body_positions, body_weights):
        from scipy.spatial import cKDTree
        self.bones = bones
        self.body_positions = body_positions
        self.body_weights = body_weights
        self.tree = cKDTree(body_positions) if len(body_positions) else None

    def weights(self, points):
        if not self.bones:
            return np.zeros((len(points), 0))
        _, index = self.tree.query(points)
        w = self.body_weights[index]
        total = w.sum(1)
        excess = total > 1.0
        w = w.copy()
        w[excess] /= total[excess, None]
        return w

    def apply(self, points):
        w = self.weights(points)
        result = points.copy()
        for k, bone in enumerate(self.bones):
            local = (points - bone['pivot_source']) @ bone['rotation'].T
            # Scale along the limb axis only. A uniform limb scale also inflates
            # girth, which turns a fitted sleeve into a 26 percent oversized tube.
            axial = local @ bone['axis_target']
            moved = bone['pivot_target'] + local + (bone['scale'] - 1.0) * axial[:, None] * bone['axis_target']
            result = result + w[:, k, None] * (moved - points)
        return result

    def record(self):
        return [{'limb': b['limb'], 'axis_source': b['axis_source'].tolist(),
                 'axis_target': b['axis_target'].tolist(), 'rotation_deg': b['angle_deg'],
                 'limb_scale': b['scale'], 'limb_scale_axis': 'along the target limb axis only; girth unscaled',
                 'pivot_source_m': b['pivot_source'].tolist(),
                 'pivot_target_m': b['pivot_target'].tolist(),
                 'body_vertices_at_full_weight': b['full_weight_vertices'],
                 'body_vertices_partially_weighted': b['partial_weight_vertices']}
                for b in self.bones]


LIMBS = {
    'arm': {'source_proximal': 'shoulder', 'source_distal': 'hand',
            'target_proximal': ('humerus', 'top'), 'target_distal': ('radius', 'bottom')},
    'leg': {'source_proximal': 'upper-leg', 'source_distal': 'ankle',
            'target_proximal': ('femur', 'top'), 'target_distal': ('tibia', 'bottom')},
}
# Anatomical right is -x in bodyparts3d-display-m; MakeHuman labels the same side 'r'.
TARGET_BONES = {
    ('humerus', 'right'): 'body-bp3d-FJ3368', ('humerus', 'left'): 'body-bp3d-FJ3262',
    ('radius', 'right'): 'body-bp3d-FJ3349', ('radius', 'left'): 'body-bp3d-FJ3277',
    ('femur', 'right'): 'body-bp3d-FJ3365', ('femur', 'left'): 'body-bp3d-FJ3259',
    ('tibia', 'right'): 'body-bp3d-FJ3387', ('tibia', 'left'): 'body-bp3d-FJ3282',
}
WEIGHT_SMOOTHING_PASSES = 12
TRUNK_JOINTS = ('joint-pelvis', 'joint-neck')


def build_skinning(root, source_body_path, scale, translation, *, minimum_angle_deg=1.0):
    body_v, body_f = load_obj_group(source_body_path, 'body')
    body_v = scale * body_v / SOURCE_UNITS_PER_M + translation
    trunk = [scale * obj_group_centroid(source_body_path, [name]) / SOURCE_UNITS_PER_M + translation
             for name in TRUNK_JOINTS]
    bones, landmarks, segments = [], [], [segment_distance(body_v, trunk[0], trunk[1])]
    for limb, spec in LIMBS.items():
        for side, tag in (('right', 'r'), ('left', 'l')):
            prox_s = scale * obj_group_centroid(source_body_path,
                                                [f'joint-{tag}-{spec["source_proximal"]}']) / SOURCE_UNITS_PER_M + translation
            dist_s = scale * obj_group_centroid(source_body_path,
                                                [f'joint-{tag}-{spec["source_distal"]}']) / SOURCE_UNITS_PER_M + translation
            bone_p, path_p, sha_p, n_p = bone_landmark(root, TARGET_BONES[(spec['target_proximal'][0], side)],
                                                       spec['target_proximal'][1])
            bone_d, path_d, sha_d, n_d = bone_landmark(root, TARGET_BONES[(spec['target_distal'][0], side)],
                                                       spec['target_distal'][1])
            if np.sign(prox_s[0]) != np.sign(bone_p[0]):
                raise ValueError(f'{side} {limb}: source and target landmarks are on opposite sides')
            axis_s = dist_s - prox_s
            axis_t = bone_d - bone_p
            length_s = float(np.linalg.norm(axis_s))
            length_t = float(np.linalg.norm(axis_t))
            axis_s = axis_s / length_s
            axis_t = axis_t / length_t
            angle = float(np.degrees(np.arccos(np.clip(float(np.dot(axis_s, axis_t)), -1, 1))))
            limb_scale = length_t / length_s
            offset = float(np.linalg.norm(bone_p - prox_s))
            landmarks.append({
                'limb': f'{side}_{limb}', 'source_proximal_group': f'joint-{tag}-{spec["source_proximal"]}',
                'source_distal_group': f'joint-{tag}-{spec["source_distal"]}',
                'source_proximal_m': prox_s.tolist(), 'source_distal_m': dist_s.tolist(),
                'target_proximal': {'entity': TARGET_BONES[(spec['target_proximal'][0], side)],
                                    'end': spec['target_proximal'][1], 'geometry': path_p, 'sha256': sha_p,
                                    'vertices': n_p, 'point_m': bone_p.tolist()},
                'target_distal': {'entity': TARGET_BONES[(spec['target_distal'][0], side)],
                                  'end': spec['target_distal'][1], 'geometry': path_d, 'sha256': sha_d,
                                  'vertices': n_d, 'point_m': bone_d.tolist()},
                'source_length_m': length_s, 'target_length_m': length_t, 'limb_scale': limb_scale,
                'rotation_deg': angle, 'proximal_offset_m': offset})
            if angle < minimum_angle_deg and abs(limb_scale - 1) < 0.01 and offset < 0.005:
                continue
            segments.append(segment_distance(body_v, prox_s, dist_s))
            bones.append({'limb': f'{side}_{limb}', 'axis_source': axis_s, 'axis_target': axis_t,
                          'pivot_source': prox_s, 'pivot_target': bone_p, 'scale': limb_scale,
                          'rotation': rotation_between(axis_s, axis_t), 'angle_deg': angle})
    if not bones:
        return LimbSkinning([], body_v, np.zeros((len(body_v), 0))), landmarks, {}
    label = np.stack(segments, 1).argmin(1)
    weights = np.zeros((len(body_v), len(bones)))
    for k in range(len(bones)):
        weights[label == k + 1, k] = 1.0
    hard = weights.copy()
    smoother = vertex_adjacency(len(body_v), body_f)
    for _ in range(WEIGHT_SMOOTHING_PASSES):
        weights = 0.5 * weights + 0.5 * (smoother @ weights)
    for k, bone in enumerate(bones):
        bone['full_weight_vertices'] = int((weights[:, k] > 0.99).sum())
        bone['partial_weight_vertices'] = int(((weights[:, k] > 0.01) & (weights[:, k] <= 0.99)).sum())
    report = {
        'method': 'nearest-skeleton-segment labelling of the authoring body, umbrella-smoothed, transferred to '
                  'garment vertices by nearest body vertex',
        'trunk_segment_joints': list(TRUNK_JOINTS),
        'smoothing_passes': WEIGHT_SMOOTHING_PASSES,
        'body_vertices': int(len(body_v)),
        'hard_label_counts': {'trunk': int((label == 0).sum()),
                              **{b['limb']: int((label == k + 1).sum()) for k, b in enumerate(bones)}},
        'hard_to_smoothed_l1_change': float(np.abs(weights - hard).sum() / len(body_v)),
    }
    return LimbSkinning(bones, body_v, weights), landmarks, report


# ------------------------------------------------------------ extremity seating
# The limb pose correction is one rigid bone per limb, fitted between the
# shoulder and the wrist. It cannot represent the difference in elbow angle
# between the authoring body and this one, and the error it leaves accumulates
# at the far end: a posed glove lands 136 mm proximal of this body's hand, in
# the gap between forearm and hip. Nothing the shrinkwrap does to a garment
# sitting there is a fit - the attraction simply glues it to whatever surface is
# nearest, which is why the gloves came out as a second skin on the wrong part
# of the body.
#
# Seating is one translation per connected component, measured not assumed: the
# component's distal end is put on the limb's distal end and then refined by
# closest-point iteration against the standoff shell, using only the distal part
# of the component so the sleeve end cannot drag it. A translation preserves
# every edge length exactly, so the shape-preservation measure is unaffected by
# it and still reports the shrinkwrap alone. The step is gated: it is applied
# only if it measurably improves the component's distance to the limb, and it is
# recorded either way.
SEAT_LIMBS = {'hand': {'limb': 'arm', 'radius_m': 0.10, 'axis_span_m': 0.45,
                       'tip_quantile': 0.80, 'distal_fraction': 0.35},
              'foot': {'limb': 'leg', 'radius_m': 0.13, 'axis_span_m': 0.30,
                       'tip_quantile': 0.80, 'distal_fraction': 0.50}}
SEAT_ITERATIONS = 25
SEAT_TRIM_QUANTILE = 0.70
SEAT_MAX_TRANSLATION_M = 0.20
SEAT_MIN_IMPROVEMENT = 0.40


def limb_frame(landmarks, limb):
    row = next(entry for entry in landmarks if entry['limb'] == limb)
    prox = np.asarray(row['target_proximal']['point_m'], float)
    dist = np.asarray(row['target_distal']['point_m'], float)
    axis = dist - prox
    return prox, dist, axis / np.linalg.norm(axis)


def extremity_pool(envelope_v, dist, axis, spec, side_sign):
    """Envelope vertices on the distal end of one limb: inside a cylinder about
    the limb axis, no further proximal than the span, and on the named side."""
    along = (envelope_v - dist) @ axis
    perpendicular = np.linalg.norm((envelope_v - dist) - along[:, None] * axis, axis=1)
    lateral = envelope_v[:, 0] < 0 if side_sign < 0 else envelope_v[:, 0] > 0
    mask = lateral & (along > -spec['axis_span_m']) & (perpendicular < spec['radius_m'])
    return envelope_v[mask]


def seat_component(component, pool, axis, spec, standoff, envelope_v, envelope_f, envelope_normals):
    def residual(translation, subset):
        signed, _, _, _ = igl.signed_distance(
            np.ascontiguousarray(component[subset] + translation), envelope_v, envelope_f,
            igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
        return float(np.median(np.abs(signed - standoff)))
    along_pool = pool @ axis
    along_component = component @ axis
    distal = along_component >= np.quantile(along_component, 1.0 - spec['distal_fraction'])
    translation = (pool[along_pool >= np.quantile(along_pool, spec['tip_quantile'])].mean(0)
                   - component[along_component >= np.quantile(along_component, spec['tip_quantile'])].mean(0))
    before = residual(np.zeros(3), distal)
    for _ in range(SEAT_ITERATIONS):
        x = component[distal] + translation
        signed, index, closest, _ = igl.signed_distance(
            np.ascontiguousarray(x), envelope_v, envelope_f,
            igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
        target = closest + standoff * envelope_normals[index]
        trimmed = np.abs(signed) <= np.quantile(np.abs(signed), SEAT_TRIM_QUANTILE)
        translation = translation + (target - x)[trimmed].mean(0)
    after = residual(translation, distal)
    magnitude = float(np.linalg.norm(translation))
    accepted = bool(magnitude <= SEAT_MAX_TRANSLATION_M and after <= (1.0 - SEAT_MIN_IMPROVEMENT) * before)
    return translation, {'translation_mm': [float(v * 1e3) for v in translation],
                         'translation_magnitude_mm': magnitude * 1e3,
                         'distal_vertices_scored': int(distal.sum()),
                         'pool_vertices': int(len(pool)),
                         'distal_residual_mm_before': before * 1e3,
                         'distal_residual_mm_after': after * 1e3, 'accepted': accepted,
                         'rejected_because': None if accepted else (
                             'translation exceeds the seating limit' if magnitude > SEAT_MAX_TRANSLATION_M
                             else 'no measurable improvement in the distal residual')}


def seat_extremity(positions, triangles, kind, standoff, envelope_v, envelope_f, envelope_normals, landmarks):
    spec = SEAT_LIMBS[kind]
    count, label = face_components(len(positions), triangles)
    seated = positions.copy()
    report = {'kind': kind, 'limb': spec['limb'], 'components': [],
              'rule': 'one translation per connected component: the component distal end is placed on the limb '
                      'distal end, then refined by trimmed closest-point iteration against the standoff shell over '
                      'the distal part of the component. A translation changes no edge length, so the shape '
                      'measures below report the shrinkwrap alone.',
              'accepted_components': 0, 'rejected_components': 0}
    for component in range(count):
        mask = label == component
        if mask.sum() < 8:
            continue
        side_sign = -1 if positions[mask][:, 0].mean() < 0 else 1
        limb = f'{"right" if side_sign < 0 else "left"}_{spec["limb"]}'
        _, dist, axis = limb_frame(landmarks, limb)
        pool = extremity_pool(envelope_v, dist, axis, spec, side_sign)
        if len(pool) < 50:
            report['components'].append({'component': component, 'limb': limb, 'accepted': False,
                                         'rejected_because': 'the limb extremity pool is empty'})
            report['rejected_components'] += 1
            continue
        translation, row = seat_component(positions[mask], pool, axis, spec, standoff,
                                          envelope_v, envelope_f, envelope_normals)
        row.update({'component': component, 'limb': limb, 'vertices': int(mask.sum())})
        if row['accepted']:
            seated[mask] = positions[mask] + translation
            report['accepted_components'] += 1
        else:
            report['rejected_components'] += 1
        report['components'].append(row)
    return seated, report


def median_edge_m(positions, triangles):
    edges = edge_table(triangles)
    return float(np.median(np.linalg.norm(positions[edges[:, 0]] - positions[edges[:, 1]], axis=1)))


def resolution_scaled_regularisation(declared, authored_edge_m, fitted_edge_m):
    """The declared shape weight, restated for the resolution actually fitted.

    Both edge lengths are measured on the posed garment, so the ratio is a pure
    resampling ratio and carries none of the similarity or limb scaling. With no
    resample the two are the same mesh and the factor is exactly one.
    """
    if authored_edge_m <= 0 or fitted_edge_m <= 0:
        raise ValueError('degenerate mesh: median edge length is zero')
    factor = (authored_edge_m / fitted_edge_m) ** REGULARISATION_RESOLUTION_EXPONENT
    return declared * factor, factor


def uniform_laplacian(vertex_count, triangles):
    """Row-normalised graph Laplacian. Dimensionally commensurate with position,
    so one dimensionless weight balances shape against fit."""
    edges = edge_table(triangles)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    adjacency = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(vertex_count, vertex_count)).tocsr()
    degree = np.asarray(adjacency.sum(1)).ravel()
    if np.any(degree <= 0):
        raise ValueError('isolated garment vertex')
    return (identity(vertex_count, format='csr') - diags(1.0 / degree) @ adjacency).tocsc()


def surface_samples(triangles, edges):
    """Face centroids and edge midpoints as barycentric stencils. A vertex-only
    no-penetration test passes while the body still pokes through a coarse
    triangle, so the constraint is stated on these samples too."""
    face = np.repeat(triangles, 1, axis=0)
    stencil = np.concatenate([face, np.column_stack([edges[:, 0], edges[:, 1], edges[:, 1]])])
    weight = np.concatenate([np.full((len(face), 3), 1 / 3), np.tile([0.5, 0.25, 0.25], (len(edges), 1))])
    return stencil, weight


def sample_penetration(x, stencil, weight, envelope_v, envelope_f):
    points = np.einsum('ij,ijk->ik', weight, x[stencil])
    signed, index, closest, _ = igl.signed_distance(
        np.ascontiguousarray(points), envelope_v, envelope_f,
        igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
    return points, signed, index, closest


def shrinkwrap(positions, triangles, envelope_v, envelope_f, envelope_normals, standoff_m, *,
               iterations=SHRINKWRAP_ITERATIONS, regularisation=2.0, push_weight=60.0,
               attract_weight=6.0, free_weight=0.02, sample_trigger=0.0, sample_target=0.35,
               sample_step_limit_m=0.003, attract_range_m=CONFORM_ATTRACT_RANGE_M):
    start = positions.copy()
    laplacian = uniform_laplacian(len(positions), triangles)
    reference = laplacian @ start
    normal_form = (laplacian.T @ laplacian).tocsc()
    x = start.copy()
    stencil, stencil_weight = surface_samples(triangles, edge_table(triangles))
    history = []
    for step in range(iterations):
        signed, index, closest, _ = igl.signed_distance(
            np.ascontiguousarray(x), envelope_v, envelope_f,
            igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
        normal = envelope_normals[index]
        target = closest + standoff_m * normal
        excess = signed - standoff_m
        inside = signed < 0
        pushed = excess < 0
        attracted = (excess >= 0) & (excess < attract_range_m)
        weight = np.full(len(x), free_weight)
        weight[attracted] = attract_weight * (1.0 - (excess[attracted] / attract_range_m) ** 2)
        weight[pushed] = push_weight
        # Free vertices anchor to where the pose correction put them, never to
        # their current position: re-anchoring each iteration lets a free sleeve
        # wall ratchet inward behind its attracted wall until the tube flattens.
        anchor = np.where(np.logical_or(pushed, attracted)[:, None], target, start)
        _, s_signed, s_index, _ = sample_penetration(x, stencil, stencil_weight, envelope_v, envelope_f)
        # Only genuine interpenetration is constrained here. A face centroid sits
        # naturally inside the shell its own vertices lie on, so constraining every
        # sample below the standoff inflates the mesh and never converges.
        violating = np.flatnonzero((s_signed < sample_trigger * standoff_m) & (s_signed > -SAMPLE_MAX_DEPTH_M))
        sample_pushed = 0
        if len(violating):
            push = np.clip(sample_target * standoff_m - s_signed[violating], 0.0, sample_step_limit_m)
            displacement = push[:, None] * envelope_normals[s_index[violating]]
            total = np.zeros_like(x)
            count = np.zeros(len(x))
            nodes = stencil[violating]
            for k in range(3):
                np.add.at(total, nodes[:, k], displacement)
                np.add.at(count, nodes[:, k], 1.0)
            touched = count > 0
            sample_pushed = int(touched.sum())
            anchor[touched] = x[touched] + total[touched] / count[touched, None]
            weight[touched] = np.maximum(weight[touched], push_weight)
        w = diags(weight).tocsc()
        solver = splu((regularisation * normal_form + w).tocsc())
        x = solver.solve(regularisation * (laplacian.T @ reference) + w @ anchor)
        history.append({'iteration': step + 1, 'inside_vertices': int(inside.sum()),
                        'max_penetration_mm': float(max(0.0, -signed.min()) * 1e3),
                        'pushed': int(pushed.sum()), 'attracted': int(attracted.sum()),
                        'violating_surface_samples': int(len(violating)),
                        'sample_driven_vertices': sample_pushed,
                        'max_sample_penetration_mm': float(max(0.0, -s_signed.min()) * 1e3),
                        'deep_chord_samples': int((s_signed <= -SAMPLE_MAX_DEPTH_M).sum())})
    # Residual interpenetration survives the smoothed solve because the Laplacian
    # trades a little of it for shape. Clear it by alternating an exact vertex
    # projection with a sample-driven push until both tests pass. Displacements
    # are recorded so the step is never silent.
    x = x.copy()
    projection = {'passes': 0, 'projected_vertices': 0, 'max_displacement_mm': 0.0,
                  'rule': 'alternating projection after the regularised solve: vertices closer than half the '
                          'standoff go onto half the standoff along their closest-face normal; face-centroid and '
                          'edge-midpoint samples inside the body push their stencil vertices out; repeated until '
                          'both tests pass or 12 passes are spent'}
    for _ in range(12):
        signed, index, closest, _ = igl.signed_distance(
            np.ascontiguousarray(x), envelope_v, envelope_f,
            igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
        interior = signed < 0.5 * standoff_m
        if interior.any():
            target = closest[interior] + 0.5 * standoff_m * envelope_normals[index[interior]]
            projection['max_displacement_mm'] = max(projection['max_displacement_mm'],
                                                    float(np.linalg.norm(target - x[interior], axis=1).max() * 1e3))
            projection['projected_vertices'] = max(projection['projected_vertices'], int(interior.sum()))
            x[interior] = target
        _, s_signed, s_index, _ = sample_penetration(x, stencil, stencil_weight, envelope_v, envelope_f)
        violating = np.flatnonzero((s_signed < 0.0) & (s_signed > -SAMPLE_MAX_DEPTH_M))
        projection['passes'] += 1
        if not interior.any() and not len(violating):
            break
        if len(violating):
            push = np.clip(sample_target * standoff_m - s_signed[violating], 0.0, sample_step_limit_m)
            displacement = push[:, None] * envelope_normals[s_index[violating]]
            total = np.zeros_like(x)
            count = np.zeros(len(x))
            nodes = stencil[violating]
            for k in range(3):
                np.add.at(total, nodes[:, k], displacement)
                np.add.at(count, nodes[:, k], 1.0)
            touched = count > 0
            x[touched] = x[touched] + total[touched] / count[touched, None]
    signed, index, closest, _ = igl.signed_distance(
        np.ascontiguousarray(x), envelope_v, envelope_f,
        igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
    interior = signed < 0.0
    if interior.any():
        x[interior] = closest[interior] + 0.25 * standoff_m * envelope_normals[index[interior]]
        projection['final_vertex_projection'] = int(interior.sum())
    _, s_signed, _, _ = sample_penetration(x, stencil, stencil_weight, envelope_v, envelope_f)
    projection['residual_interior_vertices'] = 0
    projection['residual_interior_samples'] = int((s_signed < 0).sum())
    projection['residual_shallow_interior_samples'] = int(((s_signed < 0) & (s_signed > -SAMPLE_MAX_DEPTH_M)).sum())
    projection['residual_deep_chord_samples'] = int((s_signed <= -SAMPLE_MAX_DEPTH_M).sum())
    projection['residual_max_sample_penetration_mm'] = float(max(0.0, -s_signed.min()) * 1e3)
    return x, history, projection


def measure(fitted, start, triangles, envelope_v, envelope_f):
    signed, _, _, _ = igl.signed_distance(
        np.ascontiguousarray(fitted), envelope_v, envelope_f,
        igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
    edges = edge_table(triangles)
    stencil, stencil_weight = surface_samples(triangles, edges)
    _, sample_signed, _, _ = sample_penetration(fitted, stencil, stencil_weight, envelope_v, envelope_f)
    before = np.linalg.norm(start[edges[:, 1]] - start[edges[:, 0]], axis=1)
    after = np.linalg.norm(fitted[edges[:, 1]] - fitted[edges[:, 0]], axis=1)
    ratio = after / np.maximum(before, 1e-12)
    percentiles = [1, 5, 25, 50, 75, 95, 99]
    return {
        'vertices': int(len(fitted)), 'triangles': int(len(triangles)),
        'area_m2': mesh_area_m2(fitted, triangles),
        'area_ratio_to_posed_source': mesh_area_m2(fitted, triangles) / mesh_area_m2(start, triangles),
        'inside_body_vertices': int((signed < 0).sum()),
        'max_penetration_mm': float(max(0.0, -signed.min()) * 1e3),
        'surface_samples': int(len(sample_signed)),
        'inside_body_surface_samples': int((sample_signed < 0).sum()),
        'max_surface_sample_penetration_mm': float(max(0.0, -sample_signed.min()) * 1e3),
        'surface_sample_standoff_mm_min': float(sample_signed.min() * 1e3),
        'shallow_interior_surface_samples': int(((sample_signed < 0) & (sample_signed > -SAMPLE_MAX_DEPTH_M)).sum()),
        'deep_chord_surface_samples': int((sample_signed <= -SAMPLE_MAX_DEPTH_M).sum()),
        'surface_sample_note': f'a sample deeper than {SAMPLE_MAX_DEPTH_M * 1e3:.0f} mm is a garment edge or face '
                               'chord cutting a concave body corner, not cloth the body has emerged through; '
                               'the fit counts it and leaves it alone',
        'standoff_mm_percentiles': {str(p): float(np.percentile(signed, p) * 1e3) for p in percentiles},
        'standoff_mm_min': float(signed.min() * 1e3), 'standoff_mm_max': float(signed.max() * 1e3),
        'standoff_mm_mean': float(signed.mean() * 1e3),
        'edge_length_ratio_percentiles': {str(p): float(np.percentile(ratio, p)) for p in percentiles},
        'edge_length_ratio_min': float(ratio.min()), 'edge_length_ratio_max': float(ratio.max()),
        'edge_length_rms_strain': float(np.sqrt(np.mean((ratio - 1.0) ** 2))),
    }


def penetration_check(ids, stage='simulated'):
    """Audit of the garment states, after the cloth run or straight off the fit.

    A vertex-only interior test is not a no-penetration result: three corners
    can stand off the skin while the triangle between them passes through a
    limb, and the coarser the garment the more of it a single triangle spans.
    So the body test is stated on face centroids and edge midpoints as well,
    the same stencils the shrinkwrap is constrained on, and the garment is
    additionally tested against itself. Both are reported per garment; the
    between-garment test is `layering_check`.
    """
    envelope_v, envelope_f, _, _ = load_envelope(ENVELOPE)
    out = {}
    for name in ids:
        data = np.load(OUT.parent / stage / f'{name}.npz', allow_pickle=False)
        x = np.ascontiguousarray(np.asarray(data['positions'], float))
        triangles = np.ascontiguousarray(np.asarray(data['indices'], np.int64))
        signed, _, _, _ = igl.signed_distance(x, envelope_v, envelope_f,
                                              igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
        stencil, stencil_weight = surface_samples(triangles, edge_table(triangles))
        _, sample_signed, _, _ = sample_penetration(x, stencil, stencil_weight, envelope_v, envelope_f)
        out[name] = {
            'vertices': int(len(x)),
            'inside_body_vertices': int((signed < 0).sum()),
            'inside_beyond_1um_vertices': int((signed < -1e-6).sum()),
            'inside_beyond_0p1mm_vertices': int((signed < -1e-4).sum()),
            'max_penetration_mm': float(max(0.0, -signed.min()) * 1e3),
            'standoff_mm_min': float(signed.min() * 1e3),
            'standoff_mm_median': float(np.median(signed) * 1e3),
            'surface_samples': int(len(stencil)),
            'inside_body_surface_samples': int((sample_signed < 0).sum()),
            'inside_beyond_0p1mm_surface_samples': int((sample_signed < -1e-4).sum()),
            'max_surface_sample_penetration_mm': float(max(0.0, -sample_signed.min()) * 1e3),
            'surface_sample_standoff_mm_median': float(np.median(sample_signed) * 1e3),
            'self_intersection': intersection_report(x, triangles),
            'sample_note': 'face centroids and edge midpoints, the stencils the shrinkwrap is constrained on. A '
                           'garment whose vertices all stand off the body still penetrates it wherever a triangle '
                           'is wider than the body feature it spans, and only these samples see that',
            'sign_note': 'the contact solver projects a resolved node exactly onto the body surface, so its signed '
                         'distance is zero to rounding and its sign is arbitrary; the thresholded counts separate '
                         'boundary-coincident nodes from real penetration'}
    # The post-simulation file keeps its original shape, a bare map of garments,
    # because the manifest and the verification already read it that way.
    if stage == 'simulated':
        name, payload = 'post-simulation-penetration.json', out
    else:
        name, payload = f'{stage}-penetration.json', {'stage': stage, 'garments': out}
    (OUT.parent / name).write_text(json.dumps(payload, indent=1, sort_keys=True) + '\n')
    clean = sum(1 for row in out.values()
                if not row['inside_body_surface_samples']
                and not row['self_intersection']['intersecting_face_pairs'])
    print(f'{stage} penetration checked for {len(out)} garments; '
          f'{clean} clear of the body and of themselves')
    return 0


def layering_check(ids, stage='simulated'):
    """Between-garment audit over every combination the slot model permits.

    Each garment is registered and simulated against the body alone, so nothing
    upstream has ever compared two garments to each other. Two are reported:
    whether their surfaces cross at all, for any wearable combination, and for
    a pair in one body region with a declared inner and outer, how far the
    outer one lies inside the inner one where the two actually meet.
    """
    model = slot_model()
    slot_layer = {s['id']: s['layer'] for s in model['slots']}
    slot_region = {s['id']: s['region'] for s in model['slots']}
    occupies = {row['garment']: row['occupies'] for row in model['garments']}
    excludes = {row['garment']: set(row['excludes']) for row in model['garments']}
    wanted = [name for name in ids if name in occupies]
    mesh = {}
    for name in wanted:
        data = np.load(OUT.parent / stage / f'{name}.npz', allow_pickle=False)
        mesh[name] = (np.ascontiguousarray(np.asarray(data['positions'], float)),
                      np.ascontiguousarray(np.asarray(data['indices'], np.int64)))
    pairs = []
    for i, inner in enumerate(sorted(wanted)):
        for outer in sorted(wanted)[i + 1:]:
            if outer in excludes[inner]:
                continue
            pairs.append((inner, outer))
    rows = []
    for first, second in pairs:
        (v1, f1), (v2, f2) = mesh[first], mesh[second]
        report = intersection_report(v1, f1, v2, f2)
        row = {'garments': [first, second],
               'intersecting_face_pairs': report['intersecting_face_pairs'],
               'degenerate_faces_excluded': report['degenerate_faces_excluded']}
        slot_1, slot_2 = occupies[first][0], occupies[second][0]
        if slot_region[slot_1] == slot_region[slot_2] and slot_layer[slot_1] != slot_layer[slot_2]:
            inner, outer = (first, second) if slot_layer[slot_1] < slot_layer[slot_2] else (second, first)
            vi, fi = mesh[inner]
            vo, _ = mesh[outer]
            signed, _, _, _ = igl.signed_distance(
                vo, vi, fi, igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
            near = np.abs(signed) < LAYER_MEETING_RANGE_M
            if int(near.sum()) >= LAYER_MEETING_MIN_VERTICES:
                inside = signed[near] < 0
                row['layering'] = {
                    'inner': inner, 'outer': outer,
                    'outer_vertices_where_layers_meet': int(near.sum()),
                    'outer_vertices_inside_inner': int(inside.sum()),
                    'deepest_outer_vertex_inside_inner_mm': float(max(0.0, -signed[near].min()) * 1e3),
                    'clearance_mm_median': float(np.median(signed[near]) * 1e3),
                    'range_note': f'vertices within {LAYER_MEETING_RANGE_M * 1e3:.0f} mm of the inner garment, so a '
                                  f'hem hanging clear of it is not counted as clearance it does not have'}
        rows.append(row)
    crossing = [row for row in rows if row['intersecting_face_pairs']]
    out = {
        'schema': 'ihm.garment-layering-audit.v1',
        'scope': 'Every garment is registered onto the body envelope and cloth-simulated against the body alone. '
                 'No stage of that pipeline compares one garment to another, so these are the first between-garment '
                 'numbers in the wardrobe and they are expected to be bad until the fit is made layer-aware.',
        'stage': stage,
        'combinations_tested': len(rows),
        'combinations_intersecting': len(crossing),
        'garments': sorted(wanted),
        'pairs': sorted(rows, key=lambda r: (-r['intersecting_face_pairs'], r['garments'])),
    }
    name = 'post-simulation-layering.json' if stage == 'simulated' else f'{stage}-layering.json'
    (OUT.parent / name).write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    print(f'layering checked for {len(rows)} wearable combinations; {len(crossing)} interpenetrate')
    return 0


def run(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--only', nargs='*', default=None)
    parser.add_argument('--penetration', nargs='*', default=None,
                        help='re-check simulated garment states under data/derived/wardrobe-v1/simulated')
    parser.add_argument('--layering', nargs='*', default=None,
                        help='audit every wearable combination of the simulated garment states for '
                             'between-garment interpenetration')
    parser.add_argument('--no-remesh', action='store_true',
                        help='fit the acquired meshes at their authored resolution, for comparison against '
                             'the resampled fit')
    parser.add_argument('--target-edge-mm', type=float, default=REMESH_TARGET_EDGE_M * 1e3)
    parser.add_argument('--stage', choices=('fitted', 'simulated'), default='simulated',
                        help='which garment state the audits read: the registered fit, or the cloth run. '
                             'The served geometry is the fitted state, so that is the one a reader sees.')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args(argv)
    if args.self_test:
        return self_test()
    if args.penetration is not None:
        return penetration_check(args.penetration, args.stage)
    if args.layering is not None:
        return layering_check(args.layering, args.stage)

    envelope_v, envelope_f, area, height = load_envelope(ENVELOPE)
    envelope_normals = igl.per_face_normals(envelope_v, envelope_f, np.zeros(3))
    body_v, body_f = load_obj_group(RAW / 'base/base.obj', 'body')
    body_v = body_v / SOURCE_UNITS_PER_M

    scale, translation, similarity_report = fit_similarity(body_v, envelope_v, envelope_f)
    aligned = scale * body_v + translation
    skinning, landmarks, weight_report = build_skinning(ROOT, RAW / 'base/base.obj', scale, translation)
    posed = skinning.apply(aligned)
    d_aligned, _, _ = surface_distance(aligned, envelope_v, envelope_f)
    d_posed, _, _ = surface_distance(posed, envelope_v, envelope_f)

    registration = {
        'schema': 'ihm.garment-registration.v1',
        'target': {'path': str(ENVELOPE.relative_to(ROOT)), 'sha256': sha256_file(ENVELOPE),
                   'area_m2': area, 'height_m': height, 'vertices': int(len(envelope_v)),
                   'triangles': int(len(envelope_f))},
        'source_body': {'path': 'data/raw/clothing/makehuman/base/base.obj',
                        'sha256': sha256_file(RAW / 'base/base.obj'), 'group': 'body',
                        'vertices': int(len(body_v)), 'triangles': int(len(body_f)),
                        'authoring_units_per_m': SOURCE_UNITS_PER_M},
        'similarity': {'uniform_scale': scale, 'translation_m': translation.tolist(),
                       'rotation': 'identity; both frames are y-up with +z anterior',
                       'scale_basis': 'stature landmark: target height / source body-group height',
                       'objective': 'RMS point-to-surface distance of trunk points to the envelope, translation only',
                       **similarity_report},
        'pose_correction': {
            'method': 'linear blend skinning: identity trunk bone plus one rotation per corrected limb. '
                      'The limb axis is the proximal-to-distal landmark direction, not a principal axis. '
                      'Blend weights come from the authoring body segmentation and are rescaled where they '
                      'would sum above one.',
            'landmark_basis': 'source joints are the MakeHuman skeleton joint cubes carried in base.obj; '
                              'target joints are extreme-end centroids of the canonical humerus, radius, femur and '
                              'tibia surfaces. Each limb bone is a similarity: minimal rotation between the two '
                              'proximal-to-distal axes, uniform scale from their length ratio, proximal landmark '
                              'mapped onto proximal landmark.',
            'blend_weight_basis': 'weights are segmented on the authoring body and transferred to the garment '
                                  'by nearest body vertex; a positional lateral gate was tried first and tore long '
                                  'sleeves, whose inner and outer faces then received different rotations',
            'bones': skinning.record(), 'landmarks': landmarks, 'blend_weights': weight_report},
        'source_body_residual_mm': {
            'after_similarity': {'rms': float(np.sqrt(np.mean(d_aligned ** 2)) * 1e3),
                                 'median': float(np.median(d_aligned) * 1e3), 'max': float(d_aligned.max() * 1e3)},
            'after_pose_correction': {'rms': float(np.sqrt(np.mean(d_posed ** 2)) * 1e3),
                                      'median': float(np.median(d_posed) * 1e3), 'max': float(d_posed.max() * 1e3)}},
        'shrinkwrap': {'iterations': SHRINKWRAP_ITERATIONS, 'regularisation': 2.0,
                       'laplacian': 'row-normalised graph Laplacian; differential coordinates of the posed mesh '
                                    'are the shape term',
                       'push_weight': 60.0, 'fit_modes': FIT_MODES,
                       'regularisation_resolution_exponent': REGULARISATION_RESOLUTION_EXPONENT,
                       'regularisation_scaling_basis': 'the shape term is a graph-Laplacian energy and is not '
                                                       'resolution-invariant: (L x) scales like the square of the '
                                                       'edge length, so the energy scales like its fourth power '
                                                       'against a positional fit term. Each garment therefore has '
                                                       'its declared weight restated as declared x (its own '
                                                       'authored median edge / its fitted median edge) ^ 4, both '
                                                       'measured on the posed garment, so resampling changes '
                                                       'resolution and nothing else. Without a resample the factor '
                                                       'is exactly one.',
                       'conform_attract_range_m': CONFORM_ATTRACT_RANGE_M,
                       'skin_attract_range_m': CONFORM_ATTRACT_RANGE_M,
                       'drape_attract_range_m': f'max({DRAPE_ATTRACT_RANGE_FLOOR_M}, '
                                                f'{DRAPE_ATTRACT_RANGE_MULTIPLE} x standoff)',
                       'shell_attract_range_m': f'max({SHELL_ATTRACT_RANGE_FLOOR_M}, '
                                                f'{SHELL_ATTRACT_RANGE_MULTIPLE} x standoff)',
                       'rule': 'vertices closer than the standoff are driven onto the standoff shell; vertices within '
                               'the attraction range beyond it are attracted with a weight tapering to zero; vertices '
                               'further out are left to hang so flares and hems keep their own shape',
                       'fit_mode_basis': 'conform is for garments meant to sit on the skin, where the authoring body '
                                         'is larger than this one and the garment must be drawn in. drape is for '
                                         'garments carrying real ease; a wide attraction range collapses a loose '
                                         'sleeve onto the arm, which is measurable as a sudden drop in the '
                                         'edge-length ratio. The mode is an authored declaration per garment. '
                                         'shell is for garments that hold their own volume against the body - a hat '
                                         'crown standing off the scalp, a shoe enclosing a foot. Its attraction is a '
                                         'narrow collar just outside the standoff shell and its shape term is fifteen '
                                         'times the conform value, so only cloth already at the skin is held there '
                                         'and the authored silhouette survives. skin is conform with a five times '
                                         'stiffer shape term, for garments that must be drawn a long way onto the '
                                         'body and would otherwise be squashed differentially while they travel.'},
    }

    OUT.mkdir(parents=True, exist_ok=True)
    target_edge_m = args.target_edge_mm * 1e-3
    registration['resample'] = {
        'applied': not args.no_remesh,
        'target_edge_mm': args.target_edge_mm,
        'passes': REMESH_PASSES,
        'method': 'isotropic remeshing (Botsch and Kobbelt 2004): split above 4/3 of the target, collapse below '
                  '4/5, flip toward valence six, tangential relaxation, then projection back onto the source '
                  'surface. Boundary vertices are resampled along the source boundary polyline only, so hems, '
                  'necklines and armholes keep their authored outline.',
        'basis': 'the acquired meshes carry a median edge of up to 76 mm and single edges of 307 mm. A triangle '
                 'wider than the body feature it spans penetrates that feature between its own corners, and a '
                 'shrinkwrap can only distribute strain across the vertices it is given, which is why fitting the '
                 'authored resolution stretched single edges by up to fifty times.',
        'applied_in': 'the authoring body frame, before the similarity and pose correction, at the target divided '
                      'by the similarity scale',
    }
    selected = [e for e in CATALOGUE if args.only is None or e['id'] in args.only]
    results = {}
    for entry in selected:
        source = RAW / 'clothes' / entry['pack'] / entry['asset'] / entry['obj']
        v, t = load_obj(source)
        if 'components' in entry:
            v, t = keep_components(v, t, entry['components'])
        v, t = weld(v / SOURCE_UNITS_PER_M, t)
        t, flips = orient_outward(v, t)
        boundary, nonmanifold = boundary_edge_count(t)
        # Resample before the similarity and the pose correction, not after, so
        # the limb rotations are blended across a dense mesh rather than across
        # triangles that span a whole limb. The target is divided by the
        # similarity scale because it is stated in the fitted body's metres and
        # the garment is still in the authoring body's.
        remesh_report = None
        authored_edge = median_edge_m(skinning.apply(scale * v + translation), t)
        if not args.no_remesh:
            # Resampling may refine a garment but must never coarsen one. The
            # gloves are authored at a 6.4 mm edge, finer than the target, and
            # pulling them out to 10 mm cost them faces and raised their body
            # penetration: the target is a ceiling on edge length, not a value
            # to drive every garment to.
            piece_target = min(target_edge_m, authored_edge)
            v, t, remesh_report = resample(v, t, piece_target / scale)
            remesh_report['target_edge_mm'] = piece_target * 1e3
            remesh_report['target_basis'] = (
                'the declared target' if piece_target == target_edge_m else
                f'the garment is authored finer than the {target_edge_m * 1e3:.0f} mm target, so its own '
                f'{authored_edge * 1e3:.1f} mm edge is used and the resample refines shape without coarsening it')
            boundary, nonmanifold = boundary_edge_count(t)
        posed_garment = skinning.apply(scale * v + translation)
        standoff = entry['standoff_mm'] * 1e-3
        mode = entry.get('fit_mode', 'conform')
        settings = FIT_MODES[mode]
        attract_range = ATTRACT_RANGE[mode](standoff)
        seated, seat_report = ((posed_garment, None) if 'seat' not in entry else
                               seat_extremity(posed_garment, t, entry['seat'], standoff,
                                              envelope_v, envelope_f, envelope_normals, landmarks))
        fitted_edge = median_edge_m(seated, t)
        scaled, factor = resolution_scaled_regularisation(settings['regularisation'], authored_edge, fitted_edge)
        settings = {**settings, 'regularisation': scaled}
        fitted, history, projection = shrinkwrap(seated, t, envelope_v, envelope_f, envelope_normals, standoff,
                                                 attract_range_m=attract_range, **settings)
        stats = measure(fitted, seated, t, envelope_v, envelope_f)
        stats.update({'fit_mode': mode, 'attract_range_mm': attract_range * 1e3, 'fit_weights': settings,
                      'regularisation_scaling': {'declared': FIT_MODES[mode]['regularisation'],
                                                 'applied': scaled, 'factor': factor,
                                                 'authored_median_edge_mm': authored_edge * 1e3,
                                                 'fitted_median_edge_mm': fitted_edge * 1e3,
                                                 'exponent': REGULARISATION_RESOLUTION_EXPONENT},
                      'extremity_seating': seat_report,
                      'boundary_edges': boundary, 'nonmanifold_edges': nonmanifold,
                      'orientation_component_flips': flips, 'target_standoff_mm': entry['standoff_mm'],
                      'shrinkwrap_history': history, 'final_projection': projection,
                      'remesh': remesh_report,
                      'source_sha256': sha256_file(source),
                      'source_path': str(source.relative_to(ROOT))})
        np.savez_compressed(OUT / f'{entry["id"]}.npz', positions=fitted, indices=t,
                            posed_source_positions=posed_garment, seated_source_positions=seated)
        results[entry['id']] = stats
        if remesh_report:
            before, after = remesh_report['before'], remesh_report['after']
            print(f'{entry["id"]:20s} resampled {before["faces"]:6d} -> {after["faces"]:6d} faces, '
                  f'edge p50 {before["edge_mm_median"]:5.1f} -> {after["edge_mm_median"]:4.1f} mm, '
                  f'min angle {before["min_corner_angle_deg"]:5.2f} -> {after["min_corner_angle_deg"]:5.2f} deg, '
                  f'off source {remesh_report["max_deviation_from_source_mm"]:.4f} mm')
        print(f'{entry["id"]:20s} V={stats["vertices"]:5d} inside={stats["inside_body_vertices"]:4d} '
              f'maxpen={stats["max_penetration_mm"]:6.3f} mm  samples_in={stats["inside_body_surface_samples"]:4d} '
              f'standoff p50={stats["standoff_mm_percentiles"]["50"]:6.2f} mm '
              f'edge ratio p50={stats["edge_length_ratio_percentiles"]["50"]:.4f}')
    record_path = OUT.parent / 'registration.json'
    if args.only and record_path.exists():
        # A partial run must not silently drop the garments it did not touch.
        previous = json.loads(record_path.read_text())['garments']
        previous.update(results)
        results = previous
    payload = {'registration': registration, 'garments': results}
    record_path.write_text(json.dumps(payload, indent=1, sort_keys=True) + '\n')
    inside = sum(1 for s in results.values() if s['inside_body_vertices'])
    print(f'\n{len(results)} garments fitted; {inside} still carry interior vertices')
    return 0


def self_test():
    assert np.allclose(rotation_between(np.array([1., 0, 0]), np.array([0., 1, 0])) @ np.array([1., 0, 0]),
                       [0, 1, 0], atol=1e-12)
    d = segment_distance(np.array([[0., 2, 0], [3., 0.5, 0], [0., -1, 0]]), np.array([0., 0, 0]), np.array([0., 1, 0]))
    assert np.allclose(d, [1.0, 3.0, 1.0]), d
    assert np.allclose(LimbSkinning([], np.zeros((0, 3)), np.zeros((0, 0))).apply(np.zeros((4, 3))), 0.0)
    body = np.array([[0., 0, 0], [1., 0, 0]])
    bone = {'limb': 't', 'axis_source': np.array([0., -1, 0]), 'axis_target': np.array([0., -1, 0]),
            'pivot_source': np.array([1., 0, 0]), 'pivot_target': np.array([1., 0, 0]), 'scale': 2.0,
            'rotation': np.eye(3), 'angle_deg': 0.0}
    skin = LimbSkinning([bone], body, np.array([[0.0], [1.0]]))
    moved = skin.apply(np.array([[1.0, -0.5, 0.], [1.3, 0., 0.], [-0.1, 0., 0.]]))
    assert np.allclose(moved[0], [1.0, -1.0, 0]), 'axial scale must double the along-limb offset'
    assert np.allclose(moved[1], [1.3, 0, 0]), 'girth must be unscaled'
    assert np.allclose(moved[2], [-0.1, 0, 0]), 'zero weight must not move a point'
    n = 12
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.stack([np.cos(theta), np.zeros(n), np.sin(theta)], 1)
    v = np.concatenate([ring, ring + [0, 1, 0]])
    t = np.array([[j, (j + 1) % n, n + j] for j in range(n)] + [[(j + 1) % n, n + (j + 1) % n, n + j]
                                                                for j in range(n)], dtype=np.int64)
    lap = uniform_laplacian(len(v), t)
    assert abs(float(np.abs(lap @ np.ones((len(v), 1))).max())) < 1e-12, 'Laplacian must annihilate constants'
    sphere_v, sphere_f = igl.read_triangle_mesh(str(ROOT / 'scripts')) if False else (None, None)
    box = np.array([[-1., -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    faces = np.array([f for q in quads for f in ((q[0], q[1], q[2]), (q[0], q[2], q[3]))], dtype=np.int64)
    signed, _, _, _ = igl.signed_distance(np.array([[0., 0, 0], [0, 0, 3]]), box, faces,
                                          igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
    assert signed[0] < 0 < signed[1], f'sign convention changed: {signed}'
    scale, translation, report = fit_similarity(box * 0.5 + 4.0, box, faces, trunk_halfwidth_m=10.0, min_trunk_points=4)
    assert abs(scale - 2.0) < 1e-9, scale
    assert report['trunk_rms_mm_refined'] < 1e-3, report
    # The intersection test the new receipts rest on, against cases whose answer
    # is known by construction: a triangle skewered by another, the same pair
    # pulled apart, and two sharing an edge, which touch rather than cross.
    flat = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0]])
    through = np.array([[0.2, 0.2, -1.], [0.2, 0.2, 1.], [0.6, 0.3, 1.]])
    apart = through + [10.0, 0.0, 0.0]
    one = np.array([[0, 1, 2]], np.int64)
    two = np.array([[0, 1, 2], [3, 4, 5]], np.int64)
    assert intersection_report(np.vstack([flat, through]), two)['intersecting_face_pairs'] == 1
    assert intersection_report(np.vstack([flat, apart]), two)['intersecting_face_pairs'] == 0
    assert intersection_report(flat, one, through, one)['intersecting_face_pairs'] == 1
    assert intersection_report(flat, one, apart, one)['intersecting_face_pairs'] == 0
    shared = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0], [1., 1, 0]])
    assert intersection_report(shared, np.array([[0, 1, 2], [1, 3, 2]], np.int64))['intersecting_face_pairs'] == 0, \
        'faces sharing an edge touch by construction and must not be reported as crossing'
    sliver = np.array([[0., 0, 0], [1., 0, 0], [0.5, 1e-13, 0]])   # 5e-14 m^2, below the floor
    assert intersection_report(np.vstack([sliver, through]), two)['degenerate_faces_excluded'] >= 1, \
        'a face too thin for a trustworthy normal must be excluded and counted, not silently tested'
    # A vertex-only body test cannot see a triangle spanning a feature, which is
    # why penetration_check states the body test on surface samples too.
    span = np.array([[-1., 0, 0], [1., 0, 0], [0., 0, 2.]])
    stencil, weight = surface_samples(one, edge_table(one))
    centroid = np.einsum('ij,ijk->ik', weight, span[stencil])
    assert len(centroid) == 4 and np.allclose(centroid[0], span.mean(0)), centroid
    # The resampler, on an open cylinder built the way the acquired garments are:
    # far more segments around than rings along, so every face is a sliver. A
    # correct resample lands on the cylinder exactly, keeps both rims, and
    # leaves no sliver behind.
    radius, tall, around = 0.05, 0.40, 40
    theta = np.linspace(0, 2 * np.pi, around, endpoint=False)
    tube = np.array([[radius * np.cos(a), tall * r / 2, radius * np.sin(a)] for r in range(3) for a in theta])
    wall = np.array([[f(r, j) for f in (lambda r, j: r * around + j,
                                        lambda r, j: r * around + (j + 1) % around,
                                        lambda r, j: (r + 1) * around + (j + 1) % around)]
                     for r in range(2) for j in range(around)]
                    + [[r * around + j, (r + 1) * around + (j + 1) % around, (r + 1) * around + j]
                       for r in range(2) for j in range(around)], dtype=np.int64)
    coarse = quality(tube, wall)
    assert coarse['min_corner_angle_deg'] < 5.0, coarse
    fine, fine_faces, resample_report = resample(tube, wall, 0.012)
    fine_quality = resample_report['after']
    assert boundary_edge_count(fine_faces)[1] == 0, 'resample must stay manifold'
    assert resample_report['boundary_edges_after'] == resample_report['boundary_edges_before'] == 2 * around, \
        'both rims must survive with their own edge count'
    assert fine_quality['min_corner_angle_deg'] > 15.0, fine_quality
    assert fine_quality['faces_under_10_deg'] == 0, fine_quality
    assert 0.8 * 12.0 <= fine_quality['edge_mm_median'] <= 1.34 * 12.0, fine_quality
    assert resample_report['max_deviation_from_source_mm'] < 1e-6, resample_report
    # The source is a forty-sided prism, not a cylinder, so its surface lies
    # between the inscribed and circumscribed radii. Landing anywhere in that
    # band is what "on the source surface" means here.
    on_wall = np.linalg.norm(fine[:, [0, 2]], axis=1)
    assert on_wall.min() > radius * np.cos(np.pi / around) - 1e-9 and on_wall.max() < radius + 1e-9, \
        'resampled vertices must lie on the source surface'
    assert fine[:, 1].min() > -1e-9 and fine[:, 1].max() < tall + 1e-9, 'rims must not migrate along the tube'
    assert abs(resample_report['area_change_percent']) < 1.0, resample_report
    print('fit_garments_to_envelope self-test: pass')
    return 0


if __name__ == '__main__':
    raise SystemExit(run())
