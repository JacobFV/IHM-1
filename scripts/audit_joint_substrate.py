#!/usr/bin/env python3
"""Assess the joint substrate available for emergent (contact/constraint) joints.

Three bounded questions, no promotion into the canonical frame and no writes
outside the chosen output directory:

1. What joint-relevant surfaces does the unregistered extended atlas hold, and
   are they closed after lossless vertex welding?
2. How does the extended frame relate to `bodyparts3d-display-m`, and what is
   the residual of the recorded transform on bones it was never fitted to?
3. Does any articular cartilage exist in either atlas, and if not, which
   canonical bone-surface patches would have to carry synthesized layers?

Surface residuals are exact point-to-triangle distances on strided vertex
samples. They are geometric correspondence diagnostics between two atlases that
share ancestry; they are not measurement uncertainty and not validation.
"""
import argparse
import gzip
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.anatomy import LandmarkRegistration, normalized_name

CANONICAL = 'data/derived/canonical/anatomy.json'
SOURCE_INDEX = 'data/derived/anatomy/extended/source_index.json'
FRAGMENT = 'data/derived/anatomy/extended/manifest_fragment.json'
COVERAGE = 'data/derived/anatomy/extended/coverage.json'
ROTATION = np.array([[1., 0, 0], [0, 0, 1], [0, -1, 0]])
WELD_DECIMALS = 12
CARTILAGE_GAP_M = .004
AABB_SCREEN_M = .006
MAX_QUERY_POINTS = 800

# Joint-relevant classes, first match wins; the order separates knee soft tissue
# named after its meniscus from the plain ligament bucket.
CLASSES = (
    ('nucleus_pulposus', r'nucleus pulposus'),
    ('intervertebral_disc', r'intervertebral disc'),
    ('meniscus', r'meniscus|meniscotibial|meniscopatellar|transverse ligament of knee'),
    ('cruciate_ligament', r'cruciate'),
    ('labrum', r'labrum'),
    ('articular_capsule', r'articular capsule|frenula capsulae|zona orbicularis'),
    ('articular_disc', r'articular disc|interpubic disc'),
    ('symphysis', r'symphysis'),
    ('synovial_bursa', r'bursa'),
    ('tendon_sheath', r'tendon sheath|synovial sheath|sheath of'),
    ('retinaculum', r'retinacul'),
    ('fascia_aponeurosis', r'fascia|aponeuros'),
    ('interosseous_membrane', r'interosseous membrane|obturator membrane|intercostal membrane'
                              r'|quadrangular membrane|oblique cord'),
    ('cartilage', r'cartilage'),
    ('fat_pad', r'fat pad'),
    ('ligament', r'ligament'),
)
NOT_JOINT = r'tensor fasciae|node of |lymph'
ARTICULAR_CARTILAGE = r'articular cartilage|hyaline cartilage|joint cartilage|cartilage of (head|condyle|facet)'
ORDINALS = ('first second third fourth fifth sixth seventh eighth ninth tenth eleventh twelfth').split()
DIGITS = {'hand': dict(zip(ORDINALS, ['thumb', 'index finger', 'middle finger', 'ring finger', 'little finger'])),
          'foot': dict(zip(ORDINALS, ['big toe', 'second toe', 'third toe', 'fourth toe', 'little toe']))}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def joint_class(name, system):
    lowered = name.lower()
    if system == 'lymphatic' or re.search(NOT_JOINT, lowered):
        return None
    for label, pattern in CLASSES:
        if re.search(pattern, lowered):
            return label
    return None


def synonym(name):
    """Z-Anatomy bone names the canonical name normalizer cannot already match."""
    side = ''
    if name.endswith('.l'):
        side, name = 'left ', name[:-2]
    elif name.endswith('.r'):
        side, name = 'right ', name[:-2]
    name = name.strip()
    spine = re.fullmatch(r'Vertebra (C|T|L)(\d+)', name)
    if spine:
        region = {'C': 'cervical', 'T': 'thoracic', 'L': 'lumbar'}[spine.group(1)]
        return f'{ORDINALS[int(spine.group(2)) - 1]} {region} vertebra'
    carpal = re.fullmatch(r'(Capitate|Hamate|Lunate|Pisiform|Scaphoid|Trapezium|Trapezoid|Triquetrum) bone', name)
    if carpal:
        return side + {'Triquetrum': 'triquetral'}.get(carpal.group(1), carpal.group(1).lower())
    phalanx = re.fullmatch(r'(Distal|Middle|Proximal) phalanx of (\w+) finger of (hand|foot)', name)
    if phalanx:
        return f'{phalanx.group(1).lower()} phalanx of {side}{DIGITS[phalanx.group(3)][phalanx.group(2)]}'
    fixed = {'Atlas (C1)': 'atlas', 'Axis (C2)': 'axis', 'Manubrium of sternum': 'manubrium',
             'Ethmoid bone': 'ethmoid', 'Navicular bone': f'navicular bone of {side.strip()} foot',
             'Inferior nasal concha bone': side + 'inferior nasal concha'}
    return fixed.get(name)


def welded_topology(vertices, faces):
    """Lossless weld on exact rounded coordinates, then edge-incidence topology."""
    _, inverse = np.unique(np.round(vertices, WELD_DECIMALS), axis=0, return_inverse=True)
    welded = inverse[faces]
    kept = (welded[:, 0] != welded[:, 1]) & (welded[:, 1] != welded[:, 2]) & (welded[:, 0] != welded[:, 2])
    welded, original = welded[kept], faces[kept]
    directed = np.stack([welded[:, [0, 1]], welded[:, [1, 2]], welded[:, [2, 0]]], 1).reshape(-1, 2)
    _, directed_counts = np.unique(directed, axis=0, return_counts=True)
    undirected, counts = np.unique(np.sort(directed, axis=1), axis=0, return_counts=True)
    a, b, c = (vertices[original[:, i]] for i in range(3))
    cross = np.cross(b - a, c - a)
    return {'welded_vertices': int(inverse.max()) + 1, 'welded_triangles': int(len(welded)),
            'degenerate_triangles_dropped': int((~kept).sum()),
            'duplicate_vertices_merged': int(len(vertices) - inverse.max() - 1),
            'boundary_edges': int((counts == 1).sum()), 'nonmanifold_edges': int((counts > 2).sum()),
            'edge_manifold': bool((counts == 2).all()), 'orientation_consistent': bool((directed_counts == 1).all()),
            'euler_characteristic': int((inverse.max() + 1) - len(undirected) + len(welded)),
            'closed_oriented_manifold': bool((counts == 2).all() and (directed_counts == 1).all()),
            'surface_area_m2': float(np.linalg.norm(cross, axis=1).sum() / 2),
            'signed_volume_m3': float(np.einsum('ij,ij->', a, cross) / 6)}


def _block_distance(points, p0, p1, p2):
    """Exact unsigned point-to-triangle distance for one cache-sized block."""
    ab, ac, bc = p1 - p0, p2 - p0, p2 - p1
    d1 = np.einsum('kj,ikj->ik', ab, points[:, None, :] - p0[None])
    d2 = np.einsum('kj,ikj->ik', ac, points[:, None, :] - p0[None])
    d3 = np.einsum('kj,ikj->ik', ab, points[:, None, :] - p1[None])
    d4 = np.einsum('kj,ikj->ik', ac, points[:, None, :] - p1[None])
    d5 = np.einsum('kj,ikj->ik', ab, points[:, None, :] - p2[None])
    d6 = np.einsum('kj,ikj->ik', ac, points[:, None, :] - p2[None])
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    total = va + vb + vc
    safe = np.where(total != 0, total, 1.)
    closest = p0[None] + (vb / safe)[..., None] * ab[None] + (vc / safe)[..., None] * ac[None]
    for mask, corner in (((d1 <= 0) & (d2 <= 0), p0), ((d3 >= 0) & (d4 <= d3), p1),
                         ((d6 >= 0) & (d5 <= d6), p2)):
        closest = np.where(mask[..., None], np.broadcast_to(corner[None], closest.shape), closest)
    for mask, origin, edge, num, den in (((vc <= 0) & (d1 >= 0) & (d3 <= 0), p0, ab, d1, d1 - d3),
                                         ((vb <= 0) & (d2 >= 0) & (d6 <= 0), p0, ac, d2, d2 - d6),
                                         ((va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0), p1, bc, d4 - d3,
                                          (d4 - d3) + (d5 - d6))):
        t = np.where(den != 0, num / np.where(den != 0, den, 1.), 0.)
        closest = np.where(mask[..., None], origin[None] + t[..., None] * edge[None], closest)
    return np.linalg.norm(points[:, None, :] - closest, axis=2).min(1)


def triangle_distance(points, vertices, faces, point_block=256, face_block=1024):
    """Exact unsigned point-to-triangle distance, blocked in both operands.

    A k-d tree on face centroids bounds each point block, so only faces that can
    still improve the running minimum enter the quadratic inner product.
    """
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    centroids = (a + b + c) / 3
    radii = np.maximum.reduce([np.linalg.norm(x - centroids, axis=1) for x in (a, b, c)])
    tree = cKDTree(centroids)
    best = np.empty(len(points))
    for start in range(0, len(points), point_block):
        block = points[start:start + point_block]
        centre = block.mean(0)
        seed = tree.query(block)[0].max() + radii.max() + np.linalg.norm(block - centre, axis=1).max()
        candidates = np.asarray(tree.query_ball_point(centre, seed), dtype=np.int64)
        if not len(candidates):
            candidates = np.arange(len(faces))
        result = np.full(len(block), np.inf)
        for inner in range(0, len(candidates), face_block):
            picked = candidates[inner:inner + face_block]
            result = np.minimum(result, _block_distance(block, a[picked], b[picked], c[picked]))
        best[start:start + point_block] = result
    return best


def stride(points, limit=MAX_QUERY_POINTS):
    return points[::max(1, len(points) // limit)]


def symmetric_residual(source_vertices, source_faces, target_vertices, target_faces):
    both = np.concatenate([triangle_distance(stride(source_vertices), target_vertices, target_faces),
                           triangle_distance(stride(target_vertices), source_vertices, source_faces)])
    return {'surface_rms_m': float(np.sqrt(np.mean(both ** 2))), 'surface_max_m': float(both.max()),
            'surface_mean_m': float(both.mean()), 'surface_p95_m': float(np.percentile(both, 95)),
            'sample_points': int(len(both))}


def read_geometry(path):
    payload = json.loads(gzip.open(path).read())
    return (np.asarray(payload['positions'], dtype=float).reshape(-1, 3),
            np.asarray(payload['indices'], dtype=np.int64).reshape(-1, 3))


def read_source(path):
    with np.load(path) as mesh:
        return mesh['vertices'].astype(float), mesh['faces'].astype(np.int64)


def inventory(root, index, verify_digests):
    rows, unverified = [], []
    for mesh in index['meshes']:
        label = joint_class(mesh['name'], mesh['system'])
        if label is None:
            continue
        path = root / mesh['source_geometry_path']
        if verify_digests and sha256(path) != mesh['source_geometry_sha256']:
            raise ValueError('Extended source geometry digest mismatch: ' + mesh['id'])
        vertices, faces = read_source(path)
        low, high = np.asarray(mesh['bounds'])
        rows.append({'id': mesh['id'], 'name': mesh['name'], 'joint_class': label, 'system': mesh['system'],
                     'in_joints_collection': '3: Joints' in mesh['source_collections'],
                     'source_vertices': mesh['source_vertices'], 'source_triangles': mesh['source_triangles'],
                     'source_geometry_sha256': mesh['source_geometry_sha256'],
                     'blender_world_bounds_m': [low.tolist(), high.tolist()],
                     'extent_m': (high - low).tolist(), **welded_topology(vertices, faces)})
        if re.search(ARTICULAR_CARTILAGE, mesh['name'], re.I):
            unverified.append(mesh['name'])
    return rows, unverified


def registration_audit(root, canonical, index, recorded):
    entities = {e['id']: e for e in canonical['entities']}
    by_name = {}
    for entity in canonical['entities']:
        if entity['provenance'].get('source_ids', [''])[0].startswith('bp3d'):
            by_name.setdefault(normalized_name(entity['name']), entity)
    landmarks = recorded['landmarks']
    source = np.array([l['source_rotated_m'] for l in landmarks])
    target = np.array([l['target_m'] for l in landmarks])
    fit = LandmarkRegistration(source, target)
    report = fit.report()
    folds = np.arange(len(source)) % 5
    fold_fits = [LandmarkRegistration(source[folds != f], target[folds != f]) for f in range(5)]
    held = np.concatenate([np.linalg.norm(fold_fits[f].transform(source[folds == f]) - target[folds == f], axis=1)
                           for f in range(5)])
    reproduced = {k: bool(abs(report[k] - recorded[k]) <= 1e-12) for k in
                  ('affine_rms_m', 'fit_rms_m', 'fit_max_m', 'affine_determinant')}
    reproduced['held_out_rms_m'] = bool(abs(float(np.sqrt(np.mean(held ** 2))) - recorded['held_out_rms_m']) <= 1e-12)

    by_id = {mesh['id']: mesh for mesh in index['meshes']}
    in_sample, out_of_sample = [], []
    for position, landmark in enumerate(landmarks):
        mesh = by_id[landmark['source_id']]
        vertices, faces = read_source(root / mesh['source_geometry_path'])
        vertices = vertices @ ROTATION.T
        canonical_vertices, canonical_faces = read_geometry(root / entities[landmark['target_id']]['reference_geometry']['path'])
        for tag, transform, bucket in (('in_sample', fit, in_sample), ('fold_held_out', fold_fits[folds[position]], out_of_sample)):
            moved = transform.transform(vertices)
            bucket.append({'z_name': mesh['name'], 'canonical_name': landmark['name'], 'basis': tag,
                           'centroid_residual_m': float(np.linalg.norm(
                               transform.transform(source[position:position + 1])[0] - target[position])),
                           **symmetric_residual(moved, faces, canonical_vertices, canonical_faces)})

    unfitted = []
    for mesh in index['meshes']:
        if mesh['system'] != 'skeletal' or normalized_name(mesh['name']) in by_name:
            continue
        name = synonym(mesh['name'])
        if name is None or name not in by_name:
            continue
        entity = by_name[name]
        vertices, faces = read_source(root / mesh['source_geometry_path'])
        vertices = vertices @ ROTATION.T
        canonical_vertices, canonical_faces = read_geometry(root / entity['reference_geometry']['path'])
        centre = np.mean(np.asarray(mesh['bounds']), axis=0) @ ROTATION.T
        unfitted.append({'z_name': mesh['name'], 'canonical_name': name, 'canonical_id': entity['id'],
                         'basis': 'never_used_in_fit',
                         'centroid_residual_m': float(np.linalg.norm(
                             fit.transform(centre[None])[0] - np.asarray(entity['centroid_m']))),
                         **symmetric_residual(fit.transform(vertices), faces, canonical_vertices, canonical_faces)})
    return fit, {'recorded': {k: v for k, v in recorded.items() if k not in ('landmarks', 'held_out_residuals_m')},
            'reproduced_exactly': reproduced, 'landmark_bones': in_sample,
            'fold_held_out_bones': out_of_sample, 'never_fitted_bones': unfitted,
            'never_fitted_note': 'Vertebrae, carpals, phalanges, navicular, manubrium, ethmoid and inferior '
                                 'nasal conchae are absent from the recorded 99-landmark set only because the '
                                 'two atlases name them differently. They are a true out-of-sample test of the '
                                 'recorded transform at the sites where spinal, carpal and digital joint '
                                 'structures would have to land.'}


ATTACHMENT_PROBES = (
    ('Anterior cruciate ligament.l', ('left femur', 'left tibia')),
    ('Posterior cruciate ligament.l', ('left femur', 'left tibia')),
    ('Medial meniscus.l', ('left femur', 'left tibia')),
    ('Lateral meniscus.l', ('left femur', 'left tibia')),
    ('Articular capsule of knee joint.l', ('left femur', 'left tibia', 'left patella')),
    ('Articular capsule of hip joint.l', ('left hip bone', 'left femur')),
    ('Acetabular labrum.l', ('left hip bone', 'left femur')),
    ('Ligament of head of femur.l', ('left hip bone', 'left femur')),
    ('Glenoid labrum.l', ('left scapula', 'left humerus')),
    ('Articular capsule of glenohumeral joint.l', ('left scapula', 'left humerus')),
    ('Articular capsule of elbow joint.l', ('left humerus', 'left ulna', 'left radius')),
    ('Articular capsule of radiocarpal joint.l', ('left radius', 'left lunate', 'left scaphoid')),
    ('Scapholunate interosseous ligament.l', ('left scaphoid', 'left lunate')),
    ('Anterior talofibular ligament.l', ('left fibula', 'left talus')),
    ('Calcaneofibular ligament.l', ('left fibula', 'left calcaneus')),
    ('Articular disc of temporomandibular joint.l', ('mandible', 'left temporal bone')),
    ('Intervertebral disc L4-L5', ('fourth lumbar vertebra', 'fifth lumbar vertebra')),
    ('Nucleus pulposus L4-L5', ('fourth lumbar vertebra', 'fifth lumbar vertebra')),
    ('Intervertebral disc C5-C6', ('fifth cervical vertebra', 'sixth cervical vertebra')),
    ('Anterior sacro-iliac ligament.l', ('sacrum', 'left hip bone')),
    ('Pubic symphysis', ('left hip bone', 'right hip bone')),
)


def attachment_landing(root, canonical, index, fit):
    """Do registered joint surfaces reach the canonical bones they must attach to?"""
    by_name, by_z = {}, {mesh['name']: mesh for mesh in index['meshes']}
    for entity in canonical['entities']:
        by_name.setdefault(entity['name'], entity)
    rows = []
    for z_name, bones in ATTACHMENT_PROBES:
        vertices, _ = read_source(root / by_z[z_name]['source_geometry_path'])
        moved = stride(fit.transform(vertices @ ROTATION.T), 1500)
        row = {'structure': z_name, 'bones': []}
        for bone in bones:
            bone_vertices, bone_faces = read_geometry(root / by_name[bone]['reference_geometry']['path'])
            distances = triangle_distance(moved, bone_vertices, bone_faces)
            row['bones'].append({'bone': bone, 'min_gap_m': float(distances.min()),
                                 'median_gap_m': float(np.median(distances)),
                                 'fraction_within_2mm': float((distances < .002).mean())})
        rows.append(row)
    return {'basis': 'Structure surface sampled at up to 1500 strided vertices; distance to the canonical bone '
                     'surface after the recorded transform. Median gap is span across the joint, not attachment '
                     'error; the attachment evidence is the minimum gap and the near-surface fraction.',
            'probes': rows,
            'max_min_gap_m': float(max(b['min_gap_m'] for r in rows for b in r['bones']))}


def summarize(rows, key):
    values = np.array([r[key] for r in rows])
    return {'count': int(len(values)), 'rms_m': float(np.sqrt(np.mean(values ** 2))),
            'median_m': float(np.median(values)), 'max_m': float(values.max()),
            'p95_m': float(np.percentile(values, 95))}


def cartilage_gap(root, canonical):
    bones = [e for e in canonical['entities'] if e['role'] == 'rigid_bone']
    meshes, low, high = {}, {}, {}
    for bone in bones:
        meshes[bone['id']] = read_geometry(root / bone['reference_geometry']['path'])
        low[bone['id']] = np.asarray(bone['bounds_m']['min'])
        high[bone['id']] = np.asarray(bone['bounds_m']['max'])
    trees = {i: cKDTree(meshes[i][0]) for i in meshes}
    names = {b['id']: b['name'] for b in bones}
    identifiers = [b['id'] for b in bones]
    pairs = []
    for a in range(len(identifiers)):
        for b in range(a + 1, len(identifiers)):
            i, j = identifiers[a], identifiers[b]
            gap = np.maximum(np.maximum(low[i] - high[j], low[j] - high[i]), 0)
            if np.linalg.norm(gap) < AABB_SCREEN_M:
                pairs.append((i, j))
    def near(mesh_id, other_id, margin):
        """Vertices and faces of one bone inside the other's inflated bounds."""
        vertices, faces = meshes[mesh_id]
        inside = ((vertices >= low[other_id] - margin) & (vertices <= high[other_id] + margin)).all(1)
        return np.flatnonzero(inside), faces[inside[faces].any(1)]

    contacts = []
    for i, j in pairs:
        if trees[j].query(meshes[i][0])[0].min() > 2 * CARTILAGE_GAP_M + AABB_SCREEN_M:
            continue
        window = 4 * CARTILAGE_GAP_M
        patches, gaps = [], []
        for me, other in ((i, j), (j, i)):
            vertices, faces = meshes[me]
            selected_vertices, _ = near(me, other, window)
            _, other_faces = near(other, me, window)
            if not len(selected_vertices) or not len(other_faces):
                gaps.append(np.inf)
                patches.append({'bone': names[me], 'bone_id': me, 'patch_faces': 0, 'patch_area_m2': 0.})
                continue
            distances = triangle_distance(vertices[selected_vertices], meshes[other][0], other_faces)
            gaps.append(float(distances.min()))
            close = np.zeros(len(vertices), dtype=bool)
            close[selected_vertices[distances < CARTILAGE_GAP_M]] = True
            patch = faces[close[faces].all(1)]
            area = 0. if not len(patch) else float(np.linalg.norm(np.cross(
                vertices[patch[:, 1]] - vertices[patch[:, 0]],
                vertices[patch[:, 2]] - vertices[patch[:, 0]]), axis=1).sum() / 2)
            patches.append({'bone': names[me], 'bone_id': me, 'patch_faces': int(len(patch)),
                            'patch_area_m2': area})
        if min(gaps) > CARTILAGE_GAP_M:
            continue
        contacts.append({'bones': [names[i], names[j]], 'min_surface_gap_m': float(min(gaps)),
                         'interpenetrating': bool(min(gaps) == 0.), 'patches': patches})
    return {'bone_count': len(bones), 'aabb_screened_pairs': len(pairs), 'contact_pairs': len(contacts),
            'contact_threshold_m': CARTILAGE_GAP_M, 'pairs': contacts,
            'cartilage_layers_required': 2 * len(contacts),
            'total_patch_area_m2': float(sum(p['patch_area_m2'] for c in contacts for p in c['patches']))}


def audit(root, verify_digests=True):
    started = time.time()
    inputs = {p: root / p for p in (CANONICAL, SOURCE_INDEX, FRAGMENT, COVERAGE)}
    digests = {p: sha256(f) for p, f in inputs.items()}
    canonical = json.loads(inputs[CANONICAL].read_bytes())
    index = json.loads(inputs[SOURCE_INDEX].read_bytes())
    fragment = json.loads(inputs[FRAGMENT].read_bytes())
    coverage = json.loads(inputs[COVERAGE].read_bytes())

    rows, suspicious = inventory(root, index, verify_digests)
    classes = {}
    for row in rows:
        entry = classes.setdefault(row['joint_class'], {'count': 0, 'triangles': 0, 'closed_oriented_manifold': 0,
                                                        'with_boundary': 0, 'nonmanifold': 0})
        entry['count'] += 1
        entry['triangles'] += row['source_triangles']
        entry['closed_oriented_manifold'] += row['closed_oriented_manifold']
        entry['with_boundary'] += row['boundary_edges'] > 0
        entry['nonmanifold'] += row['nonmanifold_edges'] > 0

    promoted = {e['provenance']['source_ids'][0] for e in canonical['entities']
                if e['evidence_kind'] == 'registered_geometry'}
    fit, registration = registration_audit(root, canonical, index, canonical['registrations']['z_anatomy'])
    registration['joint_structure_landing'] = attachment_landing(root, canonical, index, fit)
    gap = cartilage_gap(root, canonical)

    canonical_cartilage = [e['name'] for e in canonical['entities'] if e['role'] == 'cartilage']
    canonical_discs = [e['name'] for e in canonical['entities'] if 'intervertebral disk' in e['name']]
    articular = [n for n in canonical_cartilage if re.search(ARTICULAR_CARTILAGE, n, re.I)]

    return {
        'schema': 'ihm.joint-substrate-assessment.v1',
        'scope': 'Read-only assessment. No structure is promoted, no canonical file is modified, '
                 'no transform is committed.',
        'input_sha256': digests,
        'input_bytes': {p: f.stat().st_size for p, f in inputs.items()},
        'extended_atlas': {
            'model': fragment['models'][0]['id'], 'frame': fragment['models'][0]['frame'],
            'source_frame': fragment['models'][0]['source_frame'],
            'source_units': fragment['models'][0]['source_units'],
            'structure_count': coverage['structure_count'], 'source_triangles': coverage['source_triangles'],
            'independent_subject_count': coverage['independent_subject_count'],
            'source_revision': fragment['models'][0]['source_revision'],
            'display_transform': fragment['models'][0]['display_transform'],
            'annotation_only_articular_features': sorted({re.sub(r'\.[a-z]$', '', e['name'])
                for e in index['excluded'] if e['type'] == 'MESH' and re.search(
                    r'articular (facet|surface|circumference|cavity)|facet for|glenoid cavity|acetabul', e['name'], re.I)}),
        },
        'inventory': {'classes': classes, 'total_structures': len(rows),
                      'total_triangles': int(sum(r['source_triangles'] for r in rows)),
                      'already_promoted_into_canonical': sorted(promoted & {r['id'] for r in rows}),
                      'structures': rows},
        'registration': {
            **registration,
            'landmark_summary': {
                'in_sample_centroid': summarize(registration['landmark_bones'], 'centroid_residual_m'),
                'in_sample_surface_rms': summarize(registration['landmark_bones'], 'surface_rms_m'),
                'fold_held_out_centroid': summarize(registration['fold_held_out_bones'], 'centroid_residual_m'),
                'fold_held_out_surface_rms': summarize(registration['fold_held_out_bones'], 'surface_rms_m'),
                'never_fitted_centroid': summarize(registration['never_fitted_bones'], 'centroid_residual_m'),
                'never_fitted_surface_rms': summarize(registration['never_fitted_bones'], 'surface_rms_m'),
                'never_fitted_surface_max': summarize(registration['never_fitted_bones'], 'surface_max_m')},
            'comparison_bar': {'rejected_cervical_recipe_c1_centroid_residual_m': .036229,
                               'rejected_cervical_recipe_rms_m': .014896,
                               'inherited_global_rigid_fit_rms_m': .062356,
                               'source': 'docs/research/CERVICAL_REGISTRATION_RECIPE.md and '
                                         'docs/research/ANATOMICAL_REGISTRATION_RESIDUALS.md'}},
        'cartilage_gap': {
            **gap,
            'canonical_cartilage_entities': canonical_cartilage,
            'canonical_articular_cartilage_entities': articular,
            'canonical_intervertebral_disk_entities': canonical_discs,
            'canonical_intervertebral_disk_role': sorted({e['role'] for e in canonical['entities']
                                                          if 'intervertebral disk' in e['name']}),
            'extended_articular_cartilage_names': suspicious,
            'synthesis_operation': 'No articular cartilage surface exists in either atlas. Each contact patch '
                                   'listed here would need an offset shell raised off the bone surface along its '
                                   'outward normal, capped at the patch boundary and thickness-assigned from '
                                   'literature; both the patch delimitation and every thickness value would be '
                                   'invented, not sourced.'},
        'limits': [
            'Z-Anatomy derives in part from BodyParts3D; agreement between the two atlases is shared ancestry, '
            'not independent validation, and no second human is added.',
            'Surface residuals are exact point-to-triangle distances on strided vertex samples; they bound '
            'correspondence between two authored atlases, not anatomical accuracy.',
            'Bounding-box centers are engineering proxies, not homologous measured fiducials.',
            'Bone contact proximity is a geometric screen on open, largely non-watertight source surfaces; '
            'it does not establish an articulation, a synovial cavity or a load path.',
            'No mass, material property, joint centre or degree of freedom is created by this audit.'],
        'wall_s': round(time.time() - started, 3)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, required=True, help='Fresh output directory under data/derived')
    parser.add_argument('--skip-digest-verification', action='store_true')
    args = parser.parse_args()
    result = audit(args.root, verify_digests=not args.skip_digest_verification)
    args.output.mkdir(parents=True, exist_ok=False)
    written = {}
    for name, payload in (('inventory.json', result['inventory']),
                          ('registration.json', result['registration']),
                          ('cartilage_gap.json', result['cartilage_gap'])):
        path = args.output / name
        path.write_text(json.dumps(payload, indent=1, allow_nan=False) + '\n')
        written[name] = {'sha256': sha256(path), 'bytes': path.stat().st_size}
    manifest = {k: v for k, v in result.items() if k not in ('inventory', 'registration', 'cartilage_gap')}
    manifest['outputs'] = written
    manifest['inventory_summary'] = result['inventory']['classes']
    manifest['registration_summary'] = result['registration']['landmark_summary']
    manifest['cartilage_gap_summary'] = {k: v for k, v in result['cartilage_gap'].items() if k != 'pairs'}
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=1, allow_nan=False) + '\n')
    print(json.dumps({'output': str(args.output), 'wall_s': result['wall_s'],
                      'structures': result['inventory']['total_structures'],
                      'contact_pairs': result['cartilage_gap']['contact_pairs']}, indent=1))


if __name__ == '__main__':
    main()
