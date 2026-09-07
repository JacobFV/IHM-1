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
from ihm.assembly.garment_wardrobe import (  # noqa: E402
    CATALOGUE, SOURCE_UNITS_PER_M, boundary_edge_count, compact, edge_table, keep_components,
    load_envelope, load_obj, mesh_area_m2, orient_outward, sha256_file, weld)

ENVELOPE = ROOT / 'data/derived/outer-envelope/outer-envelope.npz'
RAW = ROOT / 'data/raw/clothing/makehuman'
OUT = ROOT / 'data/derived/wardrobe-v1/fitted'
CONFORM_ATTRACT_RANGE_M = 0.060
DRAPE_ATTRACT_RANGE_MULTIPLE = 3.0
DRAPE_ATTRACT_RANGE_FLOOR_M = 0.012
# A loose sleeve has its inner wall inside the attraction band and its outer wall
# outside it. With a weak shape term the pinned inner wall drags the free outer
# wall in and the tube flattens, so drape mode raises the shape term far above
# the attraction it has to resist.
FIT_MODES = {'conform': {'regularisation': 2.0, 'attract_weight': 6.0, 'free_weight': 0.02},
             'drape': {'regularisation': 12.0, 'attract_weight': 2.0, 'free_weight': 0.30}}
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


def penetration_check(ids):
    """Interior-vertex check on the post-simulation garment states."""
    envelope_v, envelope_f, _, _ = load_envelope(ENVELOPE)
    out = {}
    for name in ids:
        data = np.load(OUT.parent / 'simulated' / f'{name}.npz', allow_pickle=False)
        x = np.ascontiguousarray(np.asarray(data['positions'], float))
        signed, _, _, _ = igl.signed_distance(x, envelope_v, envelope_f,
                                              igl.SignedDistanceType.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
        out[name] = {
            'vertices': int(len(x)),
            'inside_body_vertices': int((signed < 0).sum()),
            'inside_beyond_1um_vertices': int((signed < -1e-6).sum()),
            'inside_beyond_0p1mm_vertices': int((signed < -1e-4).sum()),
            'max_penetration_mm': float(max(0.0, -signed.min()) * 1e3),
            'standoff_mm_min': float(signed.min() * 1e3),
            'standoff_mm_median': float(np.median(signed) * 1e3),
            'sign_note': 'the contact solver projects a resolved node exactly onto the body surface, so its signed '
                         'distance is zero to rounding and its sign is arbitrary; the thresholded counts separate '
                         'boundary-coincident nodes from real penetration'}
    (OUT.parent / 'post-simulation-penetration.json').write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    print(f'post-simulation penetration checked for {len(out)} garments')
    return 0


def run(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--only', nargs='*', default=None)
    parser.add_argument('--penetration', nargs='*', default=None,
                        help='re-check simulated garment states under data/derived/wardrobe-v1/simulated')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args(argv)
    if args.self_test:
        return self_test()
    if args.penetration is not None:
        return penetration_check(args.penetration)

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
                       'conform_attract_range_m': CONFORM_ATTRACT_RANGE_M,
                       'drape_attract_range_m': f'max({DRAPE_ATTRACT_RANGE_FLOOR_M}, '
                                                f'{DRAPE_ATTRACT_RANGE_MULTIPLE} x standoff)',
                       'rule': 'vertices closer than the standoff are driven onto the standoff shell; vertices within '
                               'the attraction range beyond it are attracted with a weight tapering to zero; vertices '
                               'further out are left to hang so flares and hems keep their own shape',
                       'fit_mode_basis': 'conform is for garments meant to sit on the skin, where the authoring body '
                                         'is larger than this one and the garment must be drawn in. drape is for '
                                         'garments carrying real ease; a wide attraction range collapses a loose '
                                         'sleeve onto the arm, which is measurable as a sudden drop in the '
                                         'edge-length ratio. The mode is an authored declaration per garment.'},
    }

    OUT.mkdir(parents=True, exist_ok=True)
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
        posed_garment = skinning.apply(scale * v + translation)
        standoff = entry['standoff_mm'] * 1e-3
        mode = entry.get('fit_mode', 'conform')
        settings = FIT_MODES[mode]
        attract_range = (CONFORM_ATTRACT_RANGE_M if mode == 'conform'
                         else max(DRAPE_ATTRACT_RANGE_FLOOR_M, DRAPE_ATTRACT_RANGE_MULTIPLE * standoff))
        fitted, history, projection = shrinkwrap(posed_garment, t, envelope_v, envelope_f, envelope_normals, standoff,
                                                 attract_range_m=attract_range, **settings)
        stats = measure(fitted, posed_garment, t, envelope_v, envelope_f)
        stats.update({'fit_mode': mode, 'attract_range_mm': attract_range * 1e3, 'fit_weights': settings,
                      'boundary_edges': boundary, 'nonmanifold_edges': nonmanifold,
                      'orientation_component_flips': flips, 'target_standoff_mm': entry['standoff_mm'],
                      'shrinkwrap_history': history, 'final_projection': projection,
                      'source_sha256': sha256_file(source),
                      'source_path': str(source.relative_to(ROOT))})
        np.savez_compressed(OUT / f'{entry["id"]}.npz', positions=fitted, indices=t,
                            posed_source_positions=posed_garment)
        results[entry['id']] = stats
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
    print('fit_garments_to_envelope self-test: pass')
    return 0


if __name__ == '__main__':
    raise SystemExit(run())
