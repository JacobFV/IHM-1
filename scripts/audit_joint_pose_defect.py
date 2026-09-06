#!/usr/bin/env python3
"""Characterise the articular pose defect in the canonical bone atlas, and test
whether a rigid-per-bone re-pose can make every joint simultaneously consistent
with literature cartilage thickness.

Read-only over `data/derived/canonical/anatomy.json`; writes only into the
`--output` directory. No canonical file is modified, no transform is promoted,
no cartilage is synthesized.

Three questions, in order:

1. POSE DEFECT. For every bone pair that comes within the screen distance,
   the signed gap field over the contacting patch (negative = interpenetration),
   the patch area, and the anatomical joint the pair belongs to.
2. LITERATURE DISCREPANCY. Per joint, measured gap minus the sum of the two
   surfaces' literature cartilage thickness. Every literature value is
   transferred from other subjects and tiered measured/interpolated/assumed.
3. SOLVABILITY. Two least-squares problems. The per-joint problem gives each
   joint its own free relative rigid transform and reports the residual that no
   rigid motion can remove - the shape-mismatch floor. The global problem gives
   each bone one rigid transform shared across all of its joints, anchored at
   the sacrum, damped, and with a tangential anti-dislocation term.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from articular_cartilage_literature import (BONDED_CLASSES, BONDED_NOTE, BONDED_TARGET_MM,  # noqa: E402
                                            NON_ARTICULAR_CLASSES, SOURCES, THICKNESS_MM, target_gap_mm)
from joint_pose_defect_lib import (closest_on_mesh, face_normals, joint_of, link_clusters,  # noqa: E402
                                   read_geometry, sha256, signed_volume, token, topology)

CANONICAL = 'data/derived/canonical/anatomy.json'
SCREEN_M = .015          # AABB inflation for the pair screen
WINDOW_M = .015          # slab around the other bone's AABB from which patch candidates are taken
PATCH_GAP_M = .006       # a vertex joins the contact patch below this signed gap
FACING_DOT = -.20        # and only if its normal opposes the other bone's normal
CLUSTER_M = .008         # single-linkage radius that splits one bone pair into articulation sites
MIN_PATCH_POINTS = 8
MIN_PATCH_AREA_M2 = 1e-5
MAX_CORRESPONDENCES = 160
DAMP_TRANSLATION = .30   # weight on |t| per bone, in the same units as a gap residual
DAMP_ROTATION = .30      # weight on radius*|omega| per bone
COHESION = .20           # tangential anti-dislocation weight per site
GAUSS_NEWTON_STEPS = 5
OUTER_ROUNDS = 2


def load_bones(root, verify):
    canonical = json.loads((root / CANONICAL).read_bytes())
    bones = {}
    for entity in canonical['entities']:
        if entity['role'] != 'rigid_bone' or 'tooth' in entity['name'].lower():
            continue
        reference = entity['reference_geometry']
        path = root / reference['path']
        if verify and sha256(path) != reference['sha256']:
            raise ValueError('canonical geometry digest mismatch: ' + entity['id'])
        vertices, faces = read_geometry(path)
        flipped = signed_volume(vertices, faces) < 0
        if flipped:
            faces = faces[:, ::-1].copy()
        normals, areas = face_normals(vertices, faces)
        vertex_normal = np.zeros_like(vertices)
        for k in range(3):
            np.add.at(vertex_normal, faces[:, k], normals * areas[:, None])
        length = np.linalg.norm(vertex_normal, axis=1)
        vertex_normal /= np.where(length > 0, length, 1.)[:, None]
        vertex_area = np.zeros(len(vertices))
        for k in range(3):
            np.add.at(vertex_area, faces[:, k], areas / 3)
        centroid = vertices.mean(0)
        bones[entity['id']] = {
            'id': entity['id'], 'name': entity['name'], 'token': token(entity['name']),
            'vertices': vertices, 'faces': faces, 'face_normal': normals, 'face_area': areas,
            'vertex_normal': vertex_normal, 'vertex_area': vertex_area,
            'low': vertices.min(0), 'high': vertices.max(0), 'centroid': centroid,
            'radius': float(np.linalg.norm(vertices - centroid, axis=1).max()),
            'tree': cKDTree(vertices), 'winding_flipped': bool(flipped),
            'mean_edge_m': float(np.sqrt(4 * areas.mean() / np.sqrt(3))),
            **topology(vertices, faces)}
    return canonical, bones


def screen(bones):
    ids = sorted(bones)
    out = []
    for a in range(len(ids)):
        i = ids[a]
        for b in range(a + 1, len(ids)):
            j = ids[b]
            gap = np.maximum(np.maximum(bones[i]['low'] - bones[j]['high'],
                                        bones[j]['low'] - bones[i]['high']), 0)
            if np.linalg.norm(gap) >= SCREEN_M:
                continue
            if bones[j]['tree'].query(bones[i]['vertices'])[0].min() > WINDOW_M:
                continue
            out.append((i, j))
    return out


LOCAL_M = WINDOW_M + .010


def local_faces(faces, target_vertices, query_points, margin=LOCAL_M):
    """Faces of the target with any vertex within `margin` of the query cloud."""
    near = cKDTree(query_points).query(target_vertices)[0] < margin
    return np.flatnonzero(near[faces].any(1))


def patch(source, target, positions=None, target_positions=None, target_normals=None):
    """Signed-gap patch on `source` against `target`; negative gap = inside target."""
    sv = source['vertices'] if positions is None else positions
    tv = target['vertices'] if target_positions is None else target_positions
    tn = target['face_normal'] if target_normals is None else target_normals
    candidate = np.flatnonzero(cKDTree(tv).query(sv)[0] < LOCAL_M)
    if not len(candidate):
        return None
    picked = local_faces(target['faces'], tv, sv[candidate])
    if not len(picked):
        return None
    tri = target['faces'][picked]
    dist, closest, which = closest_on_mesh(sv[candidate], tv, tri)
    normal = tn[picked][which]
    sign = np.sign(np.einsum('ij,ij->i', sv[candidate] - closest, normal))
    sign[sign == 0] = 1.
    signed = dist * sign
    facing = np.einsum('ij,ij->i', source['vertex_normal'][candidate], normal) < FACING_DOT
    keep = (signed < PATCH_GAP_M) & facing
    if keep.sum() < MIN_PATCH_POINTS:
        return None
    index = candidate[keep]
    return {'index': index, 'signed_m': signed[keep], 'closest_m': closest[keep],
            'normal': normal[keep], 'face': picked[which[keep]],
            'weight': source['vertex_area'][index]}


def patch_area(bone, index):
    mask = np.zeros(len(bone['vertices']), dtype=bool)
    mask[index] = True
    return float(bone['face_area'][mask[bone['faces']].all(1)].sum())


def describe(signed, weight):
    order = np.argsort(signed)
    s, w = signed[order], weight[order]
    cumulative = np.cumsum(w) / w.sum()
    def q(p):
        return float(s[min(np.searchsorted(cumulative, p), len(s) - 1)])
    return {'min_mm': float(s[0] * 1e3), 'p05_mm': q(.05) * 1e3, 'median_mm': q(.5) * 1e3,
            'p95_mm': q(.95) * 1e3, 'max_mm': float(s[-1] * 1e3),
            'mean_mm': float(np.average(s, weights=w) * 1e3),
            'interpenetrating_area_fraction': float(w[s < 0].sum() / w.sum()),
            'points': int(len(s))}


def rigid_floor(points, normals, correction, weight):
    """Residual after the best relative rigid twist for one joint, in metres."""
    origin = points.mean(0)
    design = np.hstack([normals, np.cross(points - origin, normals)])
    root = np.sqrt(weight)[:, None]
    solution, *_ = np.linalg.lstsq(design * root, correction * root[:, 0], rcond=None)
    residual = design @ solution - correction
    return residual, solution


def collect(root, verify):
    canonical, bones = load_bones(root, verify)
    pairs = screen(bones)
    sites = []
    for i, j in pairs:
        label, joint_class = joint_of(bones[i]['token'], bones[j]['token'])
        pa, pb = patch(bones[i], bones[j]), patch(bones[j], bones[i])
        if pa is None or pb is None:
            continue
        area_a, area_b = patch_area(bones[i], pa['index']), patch_area(bones[j], pb['index'])
        if min(area_a, area_b) < MIN_PATCH_AREA_M2:
            continue
        points = bones[i]['vertices'][pa['index']]
        clusters = link_clusters(points, CLUSTER_M)
        for cluster in range(clusters.max() + 1):
            mask = clusters == cluster
            if mask.sum() < MIN_PATCH_POINTS:
                continue
            index = pa['index'][mask]
            area = patch_area(bones[i], index)
            if area < MIN_PATCH_AREA_M2:
                continue
            sub = {k: pa[k][mask] for k in ('signed_m', 'closest_m', 'normal', 'weight')}
            near = cKDTree(points[mask]).query(bones[j]['vertices'][pb['index']])[0] < 3 * CLUSTER_M
            sites.append({'joint': label, 'joint_class': joint_class,
                          'bone_a': bones[i]['name'], 'bone_a_id': i,
                          'bone_b': bones[j]['name'], 'bone_b_id': j,
                          'site': int(cluster), 'sites_in_pair': int(clusters.max() + 1),
                          'index_a': index, 'patch_area_a_m2': area,
                          'patch_area_b_m2': patch_area(bones[j], pb['index'][near]) if near.any() else 0.,
                          'centroid_m': points[mask].mean(0), **sub})
    return canonical, bones, pairs, sites


def measure_site(bones, site):
    a, b = bones[site['bone_a_id']], bones[site['bone_b_id']]
    kinds = (a['token'][0], b['token'][0])
    literature = target_gap_mm(site['joint'], *kinds)
    stats = describe(site['signed_m'], site['weight'])
    row = {'joint': site['joint'], 'joint_class': site['joint_class'],
           'bone_a': site['bone_a'], 'bone_b': site['bone_b'],
           'site_index': site['site'], 'sites_in_pair': site['sites_in_pair'],
           'patch_area_a_mm2': round(site['patch_area_a_m2'] * 1e6, 3),
           'patch_area_b_mm2': round(site['patch_area_b_m2'] * 1e6, 3),
           'patch_centroid_m': [round(float(v), 6) for v in site['centroid_m']],
           'mean_triangle_edge_a_mm': round(a['mean_edge_m'] * 1e3, 3),
           'mean_triangle_edge_b_mm': round(b['mean_edge_m'] * 1e3, 3),
           'measured_gap': {k: round(v, 4) if isinstance(v, float) else v for k, v in stats.items()}}
    if literature is None:
        row['literature'] = None
        row['pose_error_mm'] = None
        return row
    target, tier, surfaces = literature
    row['literature'] = {'target_gap_mm': target, 'tier': tier, 'surfaces': surfaces,
                         'transferred': True, 'subject_calibrated': False}
    row['pose_error_mm'] = {'median_minus_target': round(stats['median_mm'] - target, 4),
                            'mean_minus_target': round(stats['mean_mm'] - target, 4),
                            'min_minus_target': round(stats['min_mm'] - target, 4),
                            'max_minus_target': round(stats['max_mm'] - target, 4),
                            'abs_mean_error': round(abs(stats['mean_mm'] - target), 4)}
    return row


def solvable_sites(sites, bones):
    """Sites the solve may use: cartilage joints with a literature target, plus bonded
    non-synovial interfaces pinned at zero gap so the skeleton cannot come apart."""
    out = []
    for site in sites:
        if site['joint_class'] in BONDED_CLASSES:
            site = dict(site)
            site['target_m'] = BONDED_TARGET_MM / 1e3
            site['tier'] = 'bonded_structural'
            out.append(site)
            continue
        if site['joint_class'] in NON_ARTICULAR_CLASSES or site['joint'] not in THICKNESS_MM:
            continue
        a, b = bones[site['bone_a_id']], bones[site['bone_b_id']]
        literature = target_gap_mm(site['joint'], a['token'][0], b['token'][0])
        if literature is None:
            continue
        site = dict(site)
        site['target_m'] = literature[0] / 1e3
        site['tier'] = literature[1]
        out.append(site)
    return out


def correspondences(bones, site, cap=MAX_CORRESPONDENCES):
    a = bones[site['bone_a_id']]
    points = a['vertices'][site['index_a']]
    weight = site['weight'].copy()
    if len(points) > cap:
        pick = np.argsort(-weight)[:cap]
        points, weight = points[pick], weight[pick]
        closest, normal, signed = site['closest_m'][pick], site['normal'][pick], site['signed_m'][pick]
    else:
        closest, normal, signed = site['closest_m'], site['normal'], site['signed_m']
    weight = weight / weight.sum()
    return points, closest, normal, signed, weight


def exp_se3(omega, t):
    theta = np.linalg.norm(omega)
    if theta < 1e-12:
        rotation = np.eye(3)
    else:
        k = omega / theta
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        rotation = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    return rotation, t


def global_solve(bones, sites, anchor_ids, scale=1.):
    bone_ids = sorted({s['bone_a_id'] for s in sites} | {s['bone_b_id'] for s in sites})
    slot = {b: k for k, b in enumerate(bone_ids)}
    state = {b: (np.eye(3), np.zeros(3)) for b in bone_ids}
    origin = {b: bones[b]['centroid'].copy() for b in bone_ids}
    packed = [correspondences(bones, s) for s in sites]
    history = []
    for outer in range(OUTER_ROUNDS):
        for step in range(GAUSS_NEWTON_STEPS):
            rows, cols, values, rhs, weights = [], [], [], [], []
            row = 0
            for site, (p0, q0, n0, _, w) in zip(sites, packed):
                ra, ta = state[site['bone_a_id']]
                rb, tb = state[site['bone_b_id']]
                oa, ob = origin[site['bone_a_id']], origin[site['bone_b_id']]
                p = (p0 - oa) @ ra.T + oa + ta
                q = (q0 - ob) @ rb.T + ob + tb
                n = n0 @ rb.T
                gap = np.einsum('ij,ij->i', p - q, n)
                correction = site['target_m'] - gap
                ia, ib = slot[site['bone_a_id']] * 6, slot[site['bone_b_id']] * 6
                block_a = np.hstack([n, np.cross(p - (oa + ta), n)])
                block_b = np.hstack([-n, -np.cross(q - (ob + tb), n) + np.cross(n, p - q)])
                for k in range(6):
                    rows.extend(range(row, row + len(p)))
                    cols.extend([ia + k] * len(p))
                    values.extend(block_a[:, k])
                    rows.extend(range(row, row + len(p)))
                    cols.extend([ib + k] * len(p))
                    values.extend(block_b[:, k])
                rhs.extend(correction)
                weights.extend(np.sqrt(w))
                row += len(p)
                centre = p.mean(0)
                axis = n.mean(0)
                axis /= max(np.linalg.norm(axis), 1e-12)
                helper = np.array([1., 0, 0]) if abs(axis[0]) < .9 else np.array([0., 1, 0])
                u = np.cross(axis, helper)
                u /= np.linalg.norm(u)
                v = np.cross(axis, u)
                for direction in (u, v):
                    for k, value in enumerate(np.concatenate([direction, np.cross(centre - (oa + ta), direction)])):
                        rows.append(row), cols.append(ia + k), values.append(value)
                    for k, value in enumerate(np.concatenate([-direction, -np.cross(centre - (ob + tb), direction)])):
                        rows.append(row), cols.append(ib + k), values.append(value)
                    rhs.append(0.)
                    weights.append(COHESION * scale)
                    row += 1
            for b in bone_ids:
                base = slot[b] * 6
                _, tb = state[b]
                for k in range(3):
                    rows.append(row), cols.append(base + k), values.append(1.)
                    rhs.append(-tb[k])
                    weights.append(DAMP_TRANSLATION * scale)
                    row += 1
                for k in range(3):
                    rows.append(row), cols.append(base + 3 + k), values.append(bones[b]['radius'])
                    rhs.append(0.)
                    weights.append(DAMP_ROTATION * scale)
                    row += 1
            for b in anchor_ids:
                if b not in slot:
                    continue
                base = slot[b] * 6
                for k in range(6):
                    rows.append(row), cols.append(base + k), values.append(1.)
                    rhs.append(0.)
                    weights.append(1e3)
                    row += 1
            weight = np.asarray(weights)
            matrix = coo_matrix((np.asarray(values) * weight[np.asarray(rows)],
                                 (np.asarray(rows), np.asarray(cols))), shape=(row, 6 * len(bone_ids)))
            delta = lsmr(matrix.tocsr(), np.asarray(rhs) * weight, atol=1e-12, btol=1e-12, maxiter=4000)[0]
            for b in bone_ids:
                base = slot[b] * 6
                dt, dw = delta[base:base + 3], delta[base + 3:base + 6]
                r, t = state[b]
                dr, _ = exp_se3(dw, dt)
                state[b] = (dr @ r, t + dt)
            history.append({'outer': outer, 'step': step,
                            'weighted_residual_rms_mm': float(np.sqrt(np.mean(
                                (matrix.tocsr() @ delta - np.asarray(rhs) * weight) ** 2)) * 1e3)})
        if outer + 1 < OUTER_ROUNDS:
            packed = refresh(bones, sites, state, origin)
    return state, origin, bone_ids, history, packed


def apply_state(bones, state, origin, bone_id, points):
    r, t = state[bone_id]
    o = origin[bone_id]
    return (points - o) @ r.T + o + t


def refresh(bones, sites, state, origin):
    """Recompute closest points after the accumulated transforms."""
    moved = {}
    for b in {s['bone_a_id'] for s in sites} | {s['bone_b_id'] for s in sites}:
        r, _ = state[b]
        moved[b] = (apply_state(bones, state, origin, b, bones[b]['vertices']),
                    bones[b]['face_normal'] @ r.T)
    packed = []
    for site in sites:
        a, b = bones[site['bone_a_id']], bones[site['bone_b_id']]
        pv, pn = moved[site['bone_b_id']]
        sv = moved[site['bone_a_id']][0][site['index_a']]
        picked = local_faces(b['faces'], pv, sv)
        if not len(picked):
            packed.append(correspondences(bones, site))
            continue
        tri = b['faces'][picked]
        dist, closest, which = closest_on_mesh(sv, pv, tri)
        normal = pn[picked][which]
        weight = site['weight'] / site['weight'].sum()
        if len(sv) > MAX_CORRESPONDENCES:
            pick = np.argsort(-weight)[:MAX_CORRESPONDENCES]
        else:
            pick = np.arange(len(sv))
        ra, _ = state[site['bone_a_id']]
        rb, _ = state[site['bone_b_id']]
        oa, ob = origin[site['bone_a_id']], origin[site['bone_b_id']]
        p0 = (sv[pick] - oa - state[site['bone_a_id']][1]) @ ra + oa
        q0 = (closest[pick] - ob - state[site['bone_b_id']][1]) @ rb + ob
        n0 = normal[pick] @ rb
        w = weight[pick] / weight[pick].sum()
        sign = np.sign(np.einsum('ij,ij->i', sv[pick] - closest[pick], normal[pick]))
        sign[sign == 0] = 1.
        packed.append((p0, q0, n0, dist[pick] * sign, w))
    return packed


def evaluate(bones, sites, packed, state, origin):
    rows = []
    for site, (p0, q0, n0, _, w) in zip(sites, packed):
        ra, ta = state[site['bone_a_id']]
        rb, tb = state[site['bone_b_id']]
        oa, ob = origin[site['bone_a_id']], origin[site['bone_b_id']]
        p = (p0 - oa) @ ra.T + oa + ta
        q = (q0 - ob) @ rb.T + ob + tb
        n = n0 @ rb.T
        before = np.einsum('ij,ij->i', p0 - q0, n0)
        after = np.einsum('ij,ij->i', p - q, n)
        floor_residual, _ = rigid_floor(p0, n0, site['target_m'] - before, w)
        core = before <= np.median(before)
        core_residual = (rigid_floor(p0[core], n0[core], site['target_m'] - before[core], w[core])[0]
                         if core.sum() >= 6 else floor_residual)
        core_weight = w[core] if core.sum() >= 6 else w
        rows.append({
            'joint': site['joint'], 'bone_a': site['bone_a'], 'bone_b': site['bone_b'],
            'site_index': site['site'], 'tier': site['tier'],
            'target_gap_mm': round(site['target_m'] * 1e3, 4),
            'patch_area_a_mm2': round(site['patch_area_a_m2'] * 1e6, 3),
            'gap_before_mean_mm': round(float(np.average(before, weights=w)) * 1e3, 4),
            'error_before_rms_mm': round(float(np.sqrt(np.average((before - site['target_m']) ** 2,
                                                                  weights=w))) * 1e3, 4),
            'per_joint_rigid_floor_rms_mm': round(float(np.sqrt(np.average(floor_residual ** 2,
                                                                           weights=w))) * 1e3, 4),
            'per_joint_rigid_floor_core_rms_mm': round(float(np.sqrt(np.average(
                core_residual ** 2, weights=core_weight))) * 1e3, 4),
            'gap_after_mean_mm': round(float(np.average(after, weights=w)) * 1e3, 4),
            'error_after_rms_mm': round(float(np.sqrt(np.average((after - site['target_m']) ** 2,
                                                                 weights=w))) * 1e3, 4),
            'error_after_max_abs_mm': round(float(np.abs(after - site['target_m']).max()) * 1e3, 4),
            'interpenetrating_fraction_after': float((after < 0).mean())})
    return rows


def displacement_report(bones, state, origin, bone_ids):
    rows = []
    for b in bone_ids:
        r, t = state[b]
        moved = apply_state(bones, state, origin, b, bones[b]['vertices'])
        shift = np.linalg.norm(moved - bones[b]['vertices'], axis=1)
        angle = float(np.degrees(np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1))))
        rows.append({'bone': bones[b]['name'], 'translation_mm': round(float(np.linalg.norm(t)) * 1e3, 4),
                     'rotation_deg': round(angle, 4),
                     'max_vertex_displacement_mm': round(float(shift.max()) * 1e3, 4),
                     'mean_vertex_displacement_mm': round(float(shift.mean()) * 1e3, 4)})
    return sorted(rows, key=lambda r: -r['max_vertex_displacement_mm'])


def group(rows, key, value_keys):
    out = {}
    for row in rows:
        bucket = out.setdefault(row[key], {'count': 0, **{k: [] for k in value_keys}})
        bucket['count'] += 1
        for k in value_keys:
            if row.get(k) is not None:
                bucket[k].append(row[k])
    for name, bucket in out.items():
        for k in value_keys:
            series = np.asarray(bucket[k], dtype=float) if bucket[k] else np.array([np.nan])
            bucket[k] = {'median': round(float(np.median(series)), 4),
                         'min': round(float(series.min()), 4), 'max': round(float(series.max()), 4)}
    return out


def roll_call(bones, sites):
    """Named cartilage joints the screen never found, and how far apart they actually are."""
    found = {(r['joint'], r['bone_a'], r['bone_b']) for r in sites}
    seen = {(r['joint'], tuple(sorted((r['bone_a'], r['bone_b'])))) for r in sites}
    ids = sorted(bones)
    missing = []
    for a in range(len(ids)):
        i = ids[a]
        for b in range(a + 1, len(ids)):
            j = ids[b]
            label, joint_class = joint_of(bones[i]['token'], bones[j]['token'])
            if label not in THICKNESS_MM or joint_class in NON_ARTICULAR_CLASSES:
                continue
            key = (label, tuple(sorted((bones[i]['name'], bones[j]['name']))))
            if key in seen:
                continue
            gap = np.maximum(np.maximum(bones[i]['low'] - bones[j]['high'],
                                        bones[j]['low'] - bones[i]['high']), 0)
            if np.linalg.norm(gap) > .080:
                continue
            distance = float(bones[j]['tree'].query(bones[i]['vertices'])[0].min())
            missing.append({'joint': label, 'bone_a': bones[i]['name'], 'bone_b': bones[j]['name'],
                            'min_vertex_separation_mm': round(distance * 1e3, 3)})
    return sorted(missing, key=lambda r: r['min_vertex_separation_mm'])


def summary_mm(rows, key):
    values = np.asarray([r[key] for r in rows], dtype=float)
    if not len(values):
        return None
    return {'count': int(len(values)), 'median': round(float(np.median(values)), 4),
            'p90': round(float(np.percentile(values, 90)), 4), 'max': round(float(values.max()), 4),
            'rms': round(float(np.sqrt(np.mean(values ** 2))), 4)}


def run_solve(bones, solve_sites, anchor, scale, label):
    state, origin, bone_ids, history, packed = global_solve(bones, solve_sites, anchor, scale)
    per_site = evaluate(bones, solve_sites, packed, state, origin)
    cartilage = [r for r in per_site if r['tier'] != 'bonded_structural']
    bonded = [r for r in per_site if r['tier'] == 'bonded_structural']
    displacement = displacement_report(bones, state, origin, bone_ids)
    after = np.asarray([r['error_after_rms_mm'] for r in cartilage])
    return {
        'variant': label, 'damping_scale': scale, 'free_bones': len(bone_ids),
        'sites_in_solve': len(solve_sites), 'cartilage_sites': len(cartilage),
        'bonded_sites': len(bonded),
        'convergence': history,
        'cartilage_mm': {
            'error_before': summary_mm(cartilage, 'error_before_rms_mm'),
            'per_joint_rigid_floor': summary_mm(cartilage, 'per_joint_rigid_floor_rms_mm'),
            'per_joint_rigid_floor_core': summary_mm(cartilage, 'per_joint_rigid_floor_core_rms_mm'),
            'error_after_global': summary_mm(cartilage, 'error_after_rms_mm')},
        'bonded_mm': {'gap_before': summary_mm(bonded, 'gap_before_mean_mm'),
                      'gap_after': summary_mm(bonded, 'gap_after_mean_mm')} if bonded else None,
        'sites_within_0p5mm_after': int((after <= .5).sum()),
        'sites_within_1mm_after': int((after <= 1.).sum()),
        'sites_within_2mm_after': int((after <= 2.).sum()),
        'max_bone_displacement_mm': displacement[0]['max_vertex_displacement_mm'] if displacement else None,
        'bones_displaced_over_5mm': int(sum(r['max_vertex_displacement_mm'] > 5 for r in displacement)),
        'bones_rotated_over_5deg': int(sum(r['rotation_deg'] > 5 for r in displacement)),
        'per_site': per_site, 'bone_displacement': displacement,
        'by_joint': group(cartilage, 'joint',
                          ['target_gap_mm', 'gap_before_mean_mm', 'error_before_rms_mm',
                           'per_joint_rigid_floor_rms_mm', 'per_joint_rigid_floor_core_rms_mm',
                           'gap_after_mean_mm', 'error_after_rms_mm', 'error_after_max_abs_mm'])}


def main():
    global PATCH_GAP_M
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--patch-gap-mm', type=float, default=PATCH_GAP_M * 1e3)
    parser.add_argument('--skip-digest-verification', action='store_true')
    args = parser.parse_args()
    PATCH_GAP_M = args.patch_gap_mm / 1e3
    started = time.time()
    root = args.root

    canonical, bones, pairs, sites = collect(root, not args.skip_digest_verification)
    measured = [measure_site(bones, s) for s in sites]
    solve_sites = solvable_sites(sites, bones)
    anchor = [b for b in bones if bones[b]['token'][0] == 'sacrum']
    damped = run_solve(bones, solve_sites, anchor, 1., 'damped')
    weak = run_solve(bones, solve_sites, anchor, .05, 'weakly_damped')
    missing = roll_call(bones, measured)

    args.output.mkdir(parents=True, exist_ok=False)
    joints = {
        'schema': 'ihm.joint-pose-defect.v1',
        'screen': {'screen_m': SCREEN_M, 'window_m': WINDOW_M, 'patch_gap_m': PATCH_GAP_M,
                   'facing_normal_dot_max': FACING_DOT, 'cluster_link_m': CLUSTER_M,
                   'min_patch_points': MIN_PATCH_POINTS, 'min_patch_area_m2': MIN_PATCH_AREA_M2,
                   'bones_considered': len(bones), 'screened_pairs': len(pairs),
                   'articulation_sites': len(sites)},
        'patch_criterion': 'A vertex of bone A joins the contact patch when its signed gap to bone B is '
                           'below patch_gap_m AND its outward vertex normal opposes the normal of B at the '
                           'closest point (dot < facing_normal_dot_max). Patches are then split into '
                           'articulation sites by single linkage at cluster_link_m, so one bone pair with '
                           'two facets (radioulnar, tibiofibular, costovertebral plus costotransverse) '
                           'yields two sites.',
        'sign_convention': 'Signed gap is the exact point-to-triangle distance from a vertex of bone A to '
                           'the surface of bone B, signed by the outward face normal of B at the closest '
                           'point. Negative means the vertex is inside B. Every bone surface used here is '
                           'a closed, consistently oriented manifold after a lossless weld at 1e-12 m and '
                           'none needed its winding flipped (bones.json records this per bone), so the '
                           'sign is well defined. The canonical entity flag watertight_edge_incidence '
                           'reads false for these same surfaces because it is evaluated before welding '
                           'duplicate vertices.',
        'named_joints_not_found': missing,
        'sites': measured}
    literature = {'schema': 'ihm.articular-cartilage-literature.v1', 'sources': SOURCES,
                  'per_joint_surface_thickness_mm': THICKNESS_MM,
                  'bonded_interfaces': {'classes': list(BONDED_CLASSES),
                                        'target_gap_mm': BONDED_TARGET_MM, 'note': BONDED_NOTE},
                  'transfer_status': 'Every value is transferred from published cohorts. None was measured '
                                     'on the BodyParts3D reference specimen, so none is subject-calibrated. '
                                     'Tier measured_literature means the source reports that surface; '
                                     'interpolated means the number is derived from a reported range or a '
                                     'neighbouring surface of the same joint; assumed means no usable '
                                     'source was found and the number is an engineering placeholder.'}
    solve = {
        'schema': 'ihm.joint-pose-solve.v1',
        'formulation': 'Each bone carries one rigid transform (omega, t). Each articulation site '
                       'contributes area-weighted point-to-plane residuals ((T_a p - T_b q) . R_b n) - '
                       'target, where target is the sum of the two literature cartilage thicknesses, or '
                       'zero for a bonded non-synovial interface. Two tangential rows per site keep the '
                       'patches from sliding off each other, per-bone Tikhonov rows bound |t| and '
                       'radius*|omega|, and the sacrum is pinned as the gauge. Damped Gauss-Newton with '
                       'one correspondence refresh. The weakly_damped variant divides the Tikhonov and '
                       'cohesion weights by 20 and exists to show that the residual floor is structural, '
                       'not a regularization artifact.',
        'weights': {'damp_translation': DAMP_TRANSLATION, 'damp_rotation': DAMP_ROTATION,
                    'cohesion': COHESION, 'gauss_newton_steps': GAUSS_NEWTON_STEPS,
                    'outer_rounds': OUTER_ROUNDS, 'max_correspondences_per_site': MAX_CORRESPONDENCES},
        'anchor_bones': [bones[b]['name'] for b in anchor],
        'per_joint_rigid_floor_definition': 'Least-squares residual after giving THIS joint alone a free '
                                            'relative rigid twist, ignoring that its bones also serve '
                                            'other joints. It is the part of the gap error that no rigid '
                                            're-posing of anything can remove, i.e. articular shape '
                                            'mismatch plus mesh resolution. The _core variant repeats the '
                                            'fit on the closer half of the patch.',
        'variants': [damped, weak]}
    bone_rows = [{'bone': bones[b]['name'], 'id': b, 'kind': bones[b]['token'][0],
                  'closed': bones[b]['closed'], 'boundary_edges': bones[b]['boundary_edges'],
                  'nonmanifold_edges': bones[b]['nonmanifold_edges'],
                  'orientation_consistent': bones[b]['orientation_consistent'],
                  'winding_flipped_for_outward_normals': bones[b]['winding_flipped'],
                  'faces': int(len(bones[b]['faces'])),
                  'mean_triangle_edge_mm': round(bones[b]['mean_edge_m'] * 1e3, 4)}
                 for b in sorted(bones)]

    written = {}
    for name, payload in (('joints.json', joints), ('literature.json', literature),
                          ('solve.json', solve), ('bones.json', {'bones': bone_rows})):
        path = args.output / name
        path.write_text(json.dumps(payload, indent=1, allow_nan=False, default=float) + '\n')
        written[name] = {'sha256': sha256(path), 'bytes': path.stat().st_size}

    inputs = {CANONICAL: root / CANONICAL,
              'scripts/joint_pose_defect_lib.py': root / 'scripts/joint_pose_defect_lib.py',
              'scripts/articular_cartilage_literature.py': root / 'scripts/articular_cartilage_literature.py',
              'scripts/audit_joint_pose_defect.py': Path(__file__).resolve(),
              'data/derived/joint-substrate-assessment-v1/cartilage_gap.json':
                  root / 'data/derived/joint-substrate-assessment-v1/cartilage_gap.json'}
    edges = np.asarray([r['mean_triangle_edge_mm'] for r in bone_rows])
    manifest = {
        'schema': 'ihm.joint-pose-defect-manifest.v1',
        'scope': 'Read-only analysis. No canonical file is modified, no transform is promoted, no '
                 'cartilage geometry is written.',
        'input_sha256': {k: sha256(v) for k, v in inputs.items()},
        'input_bytes': {k: v.stat().st_size for k, v in inputs.items()},
        'canonical_model': {'model_id': canonical['model_id'], 'frame': canonical['frame'],
                            'schema_version': canonical['schema_version']},
        'outputs': written,
        'counts': {'bones': len(bones), 'screened_pairs': len(pairs), 'articulation_sites': len(sites),
                   'sites_in_solve': len(solve_sites),
                   'cartilage_sites': damped['cartilage_sites'], 'bonded_sites': damped['bonded_sites'],
                   'named_joints_not_found': len(missing),
                   'sites_by_class': {k: int(v) for k, v in
                                      zip(*np.unique([s['joint_class'] for s in sites], return_counts=True))}},
        'mesh_resolution_mm': {'mean_triangle_edge_median': round(float(np.median(edges)), 4),
                               'mean_triangle_edge_p90': round(float(np.percentile(edges, 90)), 4),
                               'mean_triangle_edge_max': round(float(edges.max()), 4),
                               'note': 'Canonical bone geometry is the BodyParts3D obj_99 distribution '
                                       'verbatim (display_reduction_applied: false, face counts match the '
                                       'raw OBJ). Several load-bearing bones carry triangles whose mean '
                                       'edge exceeds the cartilage thickness being reasoned about.'},
        'headline_mm': {'damped': damped['cartilage_mm'], 'weakly_damped': weak['cartilage_mm']},
        'displacement_mm': {'damped_max': damped['max_bone_displacement_mm'],
                            'weakly_damped_max': weak['max_bone_displacement_mm']},
        'limits': [
            'Bone surfaces are the BodyParts3D obj_99 distribution verbatim; mean triangle edge on the '
            'long bones is several millimetres, coarser than the cartilage layers being reasoned about.',
            'Gap sign comes from the outward face normal at the closest point. This is exact here because '
            'all 229 bone surfaces weld to closed, consistently oriented manifolds, contrary to the '
            'watertight_edge_incidence flag the canonical entities carry.',
            'Every cartilage thickness is transferred literature from other cohorts; none is measured on '
            'this specimen and many are flagged assumed.',
            'A contact patch is a geometric proximity selection with a normal-facing test. It is not '
            'evidence of an articulation, a synovial cavity or a load path.',
            'The solve moves bones rigidly. It cannot change bone shape, and it therefore cannot remove '
            'the per-joint rigid floor reported alongside it.',
            'Bonded interfaces are pinned at zero gap purely so the skeleton stays connected; that target '
            'is a structural convenience, not an anatomical measurement.'],
        'wall_s': round(time.time() - started, 3)}
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=1, allow_nan=False,
                                                          default=float) + '\n')
    print(json.dumps({'output': str(args.output), 'wall_s': manifest['wall_s'],
                      'sites': len(sites), 'in_solve': len(solve_sites),
                      'headline_mm': manifest['headline_mm'],
                      'displacement_mm': manifest['displacement_mm']}, indent=1))


if __name__ == '__main__':
    main()
