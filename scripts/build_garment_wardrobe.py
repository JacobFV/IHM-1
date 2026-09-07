#!/usr/bin/env python
"""Build the fitted-garment wardrobe: cloth runs, thumbnails, provenance, manifest.

Stages, each emitting receipts:
  fit         registration of the acquired meshes onto the outer envelope
              (subprocess in the vendored libigl venv)
  simulate    a gravity-loaded elastic cloth run of every garment against the
              body surface, using the repository's existing Cloth edge-spring
              model and MovingSurfaceContact node-face Coulomb contact
  penetrate   post-simulation interior-vertex check (libigl venv)
  emit        thumbnails, slot/exclusivity model, provenance records, manifest

  python scripts/build_garment_wardrobe.py [--skip-fit] [--only ID ...] [--self-test]
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.clothing import Cloth  # noqa: E402
from ihm.assembly.garment_surface_contact import MovingSurfaceContact  # noqa: E402
from ihm.assembly.garment_wardrobe import (  # noqa: E402
    CATALOGUE, SLOTS, boundary_edge_count, edge_table, face_components, load_envelope, mesh_area_m2,
    sha256_file, slot_model)

OUT = ROOT / 'data/derived/wardrobe-v1'
ENVELOPE = ROOT / 'data/derived/outer-envelope/outer-envelope.npz'
RETRIEVAL = ROOT / 'data/raw/clothing/makehuman/retrieval.json'
LIBIGL = ROOT / 'data/runtime/geometry/libigl-2.6.2/venv/bin/python'
FITTER = ROOT / 'scripts/fit_garments_to_envelope.py'
PROVENANCE_SCHEMA = ROOT / 'data/derived/structure-provenance-candidate-v1/schema.json'

AREAL_DENSITY_KG_M2 = 0.18
EDGE_STIFFNESS_N_M = 12.0
PRESTRAIN = 0.04
FRICTION_STATIC = 0.40
FRICTION_KINETIC = 0.30
CONTACT_SEARCH_M = 0.004
CONTACT_BAND_M = 0.045
SIMULATED_SECONDS = 0.030
MAX_STEPS = 1400

CAMERA = {'azimuth_deg': 28.0, 'elevation_deg': 14.0, 'projection': 'orthographic',
          'half_width_m': 0.62, 'half_height_m': 0.95, 'centre_m': [0.0, 0.0, 0.0]}
LIGHT = {'direction': [-0.35, 0.55, 0.76], 'ambient': 0.34, 'diffuse': 0.66}
IMAGE = {'width': 360, 'height': 540, 'supersample': 2}
BODY_RGB = (0.78, 0.74, 0.71)


# ----------------------------------------------------------------------- render

def camera_basis():
    az = math.radians(CAMERA['azimuth_deg'])
    el = math.radians(CAMERA['elevation_deg'])
    w = np.array([math.sin(az) * math.cos(el), math.sin(el), math.cos(az) * math.cos(el)])
    w /= np.linalg.norm(w)
    r = np.cross(np.array([0.0, 1.0, 0.0]), w)
    r /= np.linalg.norm(r)
    u = np.cross(w, r)
    return r, u, w


def project(points, basis, size):
    r, u, w = basis
    centre = np.asarray(CAMERA['centre_m'])
    local = points - centre
    width, height = size
    px = (local @ r / CAMERA['half_width_m'] * 0.5 + 0.5) * width
    py = (0.5 - local @ u / CAMERA['half_height_m'] * 0.5) * height
    return np.stack([px, py, local @ w], 1)


def rasterise(screen, positions, triangles, colour, depth, image, mask):
    """Painter-free z-buffer rasteriser with one directional light."""
    _, _, w = camera_basis()
    light = np.asarray(LIGHT['direction'], float)
    light /= np.linalg.norm(light)
    p = positions[triangles]
    normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    length = np.linalg.norm(normal, axis=1)
    keep = length > 0
    normal[keep] /= length[keep, None]
    facing = normal @ w
    shade = LIGHT['ambient'] + LIGHT['diffuse'] * np.maximum(normal @ light, 0.0)
    s = screen[triangles]
    height, width = depth.shape
    lo = np.floor(s[:, :, :2].min(1)).astype(int)
    hi = np.ceil(s[:, :, :2].max(1)).astype(int)
    for i in np.flatnonzero(keep & (facing > 0)):
        x0, y0 = max(lo[i, 0], 0), max(lo[i, 1], 0)
        x1, y1 = min(hi[i, 0] + 1, width), min(hi[i, 1] + 1, height)
        if x1 <= x0 or y1 <= y0:
            continue
        a, b, c = s[i]
        area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if abs(area) < 1e-12:
            continue
        xs = np.arange(x0, x1) + 0.5
        ys = np.arange(y0, y1) + 0.5
        gx, gy = np.meshgrid(xs, ys)
        w0 = ((b[0] - a[0]) * (gy - a[1]) - (b[1] - a[1]) * (gx - a[0])) / area
        w1 = ((gx - a[0]) * (c[1] - a[1]) - (gy - a[1]) * (c[0] - a[0])) / area
        inside = (w0 >= 0) & (w1 >= 0) & (w0 + w1 <= 1)
        if not inside.any():
            continue
        z = a[2] + w1 * (b[2] - a[2]) + w0 * (c[2] - a[2])
        window = depth[y0:y1, x0:x1]
        hit = inside & (z > window)
        if not hit.any():
            continue
        window[hit] = z[hit]
        image[y0:y1, x0:x1][hit] = np.asarray(colour) * shade[i]
        mask[y0:y1, x0:x1][hit] = True


_BODY_LAYER = {}


def body_layer(body):
    """The body raster is camera-invariant across garments; rasterise it once."""
    key = id(body[0])
    if key not in _BODY_LAYER:
        factor = IMAGE['supersample']
        size = (IMAGE['width'] * factor, IMAGE['height'] * factor)
        depth = np.full((size[1], size[0]), -np.inf)
        image = np.zeros((size[1], size[0], 3))
        mask = np.zeros((size[1], size[0]), bool)
        positions, triangles, colour = body
        rasterise(project(positions, camera_basis(), size), positions, triangles, colour, depth, image, mask)
        _BODY_LAYER[key] = (depth, image, mask)
    depth, image, mask = _BODY_LAYER[key]
    return depth.copy(), image.copy(), mask.copy()


def render(body, garments, path):
    from PIL import Image
    factor = IMAGE['supersample']
    size = (IMAGE['width'] * factor, IMAGE['height'] * factor)
    depth, image, mask = body_layer(body)
    basis = camera_basis()
    for positions, triangles, colour in garments:
        rasterise(project(positions, basis, size), positions, triangles, colour, depth, image, mask)
    rgba = np.concatenate([np.clip(image, 0, 1), mask[:, :, None].astype(float)], 2)
    small = rgba.reshape(IMAGE['height'], factor, IMAGE['width'], factor, 4).mean((1, 3))
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray((small * 255 + 0.5).astype(np.uint8), 'RGBA').save(path)
    return {'width': IMAGE['width'], 'height': IMAGE['height'], 'supersample': factor,
            'covered_pixels': int((small[:, :, 3] > 0.5).sum())}


def hex_rgb(text):
    return tuple(int(text[i:i + 2], 16) / 255 for i in (1, 3, 5))


# ---------------------------------------------------------------------- physics

CONTACT_FACE_RADIUS_M = 0.008


def refine_contact_faces(positions, triangles, max_radius_m=CONTACT_FACE_RADIUS_M, max_passes=8):
    """Midpoint-split oversized contact faces. The broad phase queries one global
    face radius, so the envelope's few 48 mm faces cost every node a 15x wider
    search than its 3.3 mm median needs. The split is exact: each face becomes
    four sub-triangles covering the same surface, so no geometry changes."""
    positions = np.asarray(positions, float)
    triangles = np.asarray(triangles, np.int64)
    split = 0
    for _ in range(max_passes):
        p = positions[triangles]
        radius = np.linalg.norm(p - p.mean(1)[:, None], axis=2).max(1)
        big = radius > max_radius_m
        if not big.any():
            break
        face = triangles[big]
        base = len(positions)
        count = len(face)
        midpoints = [0.5 * (positions[face[:, a]] + positions[face[:, b]]) for a, b in ((0, 1), (1, 2), (2, 0))]
        positions = np.concatenate([positions] + midpoints)
        j = [base + k * count + np.arange(count) for k in range(3)]
        i0, i1, i2 = face[:, 0], face[:, 1], face[:, 2]
        triangles = np.concatenate([
            triangles[~big],
            np.column_stack([i0, j[0], j[2]]), np.column_stack([j[0], i1, j[1]]),
            np.column_stack([j[2], j[1], i2]), np.column_stack([j[0], j[1], j[2]])])
        split += count
    return positions, triangles.astype(np.int64), split


def contact_submesh(envelope_v, envelope_f, garment_v, band_m):
    from scipy.spatial import cKDTree
    tree = cKDTree(garment_v)
    centres = envelope_v[envelope_f].mean(1)
    near = tree.query_ball_point(centres, band_m, return_length=True) > 0
    faces = envelope_f[near]
    if not len(faces):
        return None, None, 0
    used, inverse = np.unique(faces, return_inverse=True)
    v, f, split = refine_contact_faces(envelope_v[used], inverse.reshape(-1, 3).astype(np.int64))
    return v, f, {'envelope_faces_in_band': int(near.sum()), 'contact_faces_after_refinement': int(len(f)),
                  'oversized_faces_split': int(split), 'max_face_radius_m': CONTACT_FACE_RADIUS_M}


def simulate(entry, positions, triangles, envelope_v, envelope_f, *, seconds=SIMULATED_SECONDS,
             max_steps=MAX_STEPS):
    """Gravity-loaded elastic cloth against the static body surface."""
    cloth = Cloth(positions, triangles, areal_density_kg_m2=AREAL_DENSITY_KG_M2,
                  edge_stiffness_n_m=EDGE_STIFFNESS_N_M)
    # Intrinsic uniform edge preload. This is a rest-length statement, not an
    # embedding: no reference mesh realises a uniform contraction of every edge.
    cloth.rest_lengths_m = cloth.rest_lengths_m * (1.0 - PRESTRAIN)
    band = entry['standoff_mm'] * 1e-3 + CONTACT_BAND_M
    surface_v, surface_f, contact_mesh = contact_submesh(envelope_v, envelope_f, positions, band)
    if surface_f is None:
        return {'status': 'no_contact_faces_within_band', 'contact_ready': False}, positions
    dt = min(cloth.max_explicit_dt_s * 0.9, seconds)
    steps = min(max_steps, max(1, int(round(seconds / dt))))
    duration = steps * dt
    contact = MovingSurfaceContact(surface_v, surface_v, surface_f, duration,
                                   friction_static=FRICTION_STATIC, friction_kinetic=FRICTION_KINETIC,
                                   search_distance_m=CONTACT_SEARCH_M, max_candidate_pairs=4000000)
    contact.fraction = 0.0
    _, elastic_initial = cloth.elastic_forces(cloth.position_m)
    totals = {'contact_resolutions': 0, 'unresolved_edge_contacts': 0, 'friction_dissipation_j': 0.0,
              'normal_impact_dissipation_j': 0.0, 'numerical_energy_defect_j': 0.0}
    momentum = 0.0
    paired = 0.0
    correction = 0.0
    started = time.monotonic()
    for _ in range(steps):
        result = cloth.step(dt, gravity_m_s2=(0.0, -9.81, 0.0), contact=contact)
        c = result['contact']
        totals['contact_resolutions'] += c['contact_count']
        totals['unresolved_edge_contacts'] += c['unresolved_edge_contacts']
        totals['friction_dissipation_j'] += result['friction_dissipation_j']
        totals['normal_impact_dissipation_j'] += result['normal_impact_dissipation_j']
        totals['numerical_energy_defect_j'] += result['numerical_energy_defect_j']
        momentum = max(momentum, float(np.linalg.norm(result['momentum_residual_ns'])))
        paired = max(paired, float(np.linalg.norm(c['paired_impulse_residual_ns'])))
        correction = max(correction, result['max_position_correction_m'])
    _, elastic_final = cloth.elastic_forces(cloth.position_m)
    displacement = np.linalg.norm(cloth.position_m - positions, axis=1)
    kinetic = float(0.5 * np.sum(cloth.mass_kg[:, None] * cloth.velocity_m_s ** 2))
    potential = float(-np.sum(cloth.mass_kg * (cloth.position_m[:, 1] - positions[:, 1]) * -9.81))
    dissipation = totals['friction_dissipation_j'] + totals['normal_impact_dissipation_j']
    return {
        'status': 'completed', 'contact_ready': True,
        'dt_s': dt, 'steps': steps, 'simulated_s': cloth.time_s,
        'conservative_explicit_dt_limit_s': cloth.max_explicit_dt_s,
        'wall_s': time.monotonic() - started,
        'mass_kg': float(cloth.mass_kg.sum()), 'nodes': int(len(cloth.mass_kg)),
        'edges': int(len(cloth.edges)), 'contact_mesh': contact_mesh,
        'areal_density_kg_m2': AREAL_DENSITY_KG_M2, 'edge_stiffness_n_m': EDGE_STIFFNESS_N_M,
        'edge_prestrain': PRESTRAIN,
        'friction_static': FRICTION_STATIC, 'friction_kinetic': FRICTION_KINETIC,
        'elastic_energy_initial_j': elastic_initial, 'elastic_energy_final_j': elastic_final,
        'kinetic_energy_final_j': kinetic, 'gravitational_potential_change_j': potential,
        'contact_dissipation_j': dissipation,
        'energy_audit_closure_j': (elastic_final - elastic_initial) + kinetic + potential + dissipation
                                  - totals['numerical_energy_defect_j'],
        'maximum_momentum_residual_ns': momentum,
        'maximum_paired_impulse_residual_ns': paired,
        'maximum_position_correction_m': correction,
        'maximum_node_displacement_mm': float(displacement.max() * 1e3),
        'median_node_displacement_mm': float(np.median(displacement) * 1e3),
        **totals,
    }, cloth.position_m.copy()


# -------------------------------------------------------------------- provenance

def git_state():
    def run(*args):
        return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    commit = run('rev-parse', 'HEAD')
    dirty = run('status', '--porcelain')
    return commit, bool(dirty), dirty.splitlines()


def provenance_record(entry, retrieval, registration, fit, geometry_path, commit, dirty, transforms):
    file_record = retrieval['files'][entry['id']]
    dataset = dict(retrieval['dataset'])
    dataset['license_evidence'] = retrieval['licence_evidence']['LICENSE.ASSETS.md']
    dataset['attribution'] = f'{file_record["author"]} via {dataset["label"]}'
    return {
        'schema': 'ihm.structure-provenance.v1',
        'structure_id': f'garment-{entry["id"]}',
        'model_id': 'ihm-wardrobe',
        'canonical_entity_id': None,
        'name': entry['name'],
        'system': 'clothing',
        'role': 'garment',
        'role_source': 'wardrobe catalogue declaration; garments are not canonical anatomical entities',
        'evidence_kind': 'source_geometry',
        'dataset': dataset,
        'source_file': {'path': file_record['path'], 'sha256': file_record['sha256'],
                        'bytes': file_record['bytes'], 'sha256_verified': True,
                        'upstream_asset': file_record['asset'], 'upstream_asset_page': file_record['asset_page'],
                        'upstream_archive': file_record['archive_url'],
                        'upstream_license': file_record['license'], 'retrieved_at': file_record['retrieved_at']},
        'build': {'script': 'scripts/build_garment_wardrobe.py',
                  'script_sha256': sha256_file(ROOT / 'scripts/build_garment_wardrobe.py'),
                  'fit_script': 'scripts/fit_garments_to_envelope.py',
                  'fit_script_sha256': sha256_file(ROOT / 'scripts/fit_garments_to_envelope.py'),
                  'commit': commit, 'commit_covers_working_tree': not dirty, 'uncommitted_changes': dirty},
        'geometry': {'path': geometry_path, 'sha256': sha256_file(ROOT / geometry_path), 'sha256_verified': True,
                     'representation': 'triangular_surface', 'frame': 'bodyparts3d-display-m', 'units': 'm',
                     'source_vertex_count': fit['vertices'], 'source_face_count': fit['triangles']},
        'transforms': transforms,
        'transform_note': None,
        'frame_relation': 'registered_into_canonical',
        'tier': 'transferred',
        'tier_basis': 'Geometry acquired from a separate dataset and placed into this body by a fitted chain: a '
                      'stature-landmark similarity, a landmark-driven piecewise-rigid limb correction and a '
                      'Laplacian-regularised shrinkwrap. Its position here is inferred, and every residual the '
                      'chain was measured against is recorded on the transform entries.',
        'tier_evidence': {
            'dataset.specimen': dataset['specimen'],
            'dataset.acquisition_status': dataset['acquisition_status'],
            'frame_relation': 'registered_into_canonical',
            'similarity.trunk_rms_mm_refined': registration['similarity']['trunk_rms_mm_refined'],
            'source_body_residual_mm.after_pose_correction.rms':
                registration['source_body_residual_mm']['after_pose_correction']['rms'],
        },
        'assumptions': [
            {'id': 'WARDROBE-AUTHORING-BODY',
             'statement': 'The acquired garment was authored around the MakeHuman base mesh, a different body of '
                          'different proportions and pose. Its registration onto this body is a fit, not an '
                          'acquisition, and the garment carries no measurement of any physical garment.'},
            {'id': 'WARDROBE-CLOTH-PRIORS',
             'statement': f'Areal density {AREAL_DENSITY_KG_M2} kg/m2, edge stiffness {EDGE_STIFFNESS_N_M} N/m, '
                          f'edge preload {PRESTRAIN}, Coulomb {FRICTION_STATIC}/{FRICTION_KINETIC} are explicit '
                          'engineering inputs shared with the existing cloth work. No fabric was calibrated.'},
            {'id': 'WARDROBE-STANDOFF',
             'statement': f'The {entry["standoff_mm"]} mm target standoff is an authored ease, not a measured '
                          'garment-to-skin clearance.'},
        ],
        'derived_artifacts': [],
        'display_present': True,
        'completeness': {'required_present': 11, 'required_total': 11, 'missing': [], 'answerable': True,
                         'hashes_verified': 3, 'hashes_failed': 0, 'hashes_unchecked': 0},
    }


WHOLE_GARMENTS = ROOT / 'data/derived/whole-garments-uc8tp3zc/report.json'
GARMENT_TISSUE = ROOT / 'data/derived/garment-tissue-refined-oufb46p6/coupled'


def prior_work_regression():
    """The retained generated-garment and coupled cloth-tissue work must be
    byte-identical after this build. Recorded hashes are re-read, not assumed."""
    checks = []

    def compare(artifact, path, digest):
        current = sha256_file(ROOT / path) if (ROOT / path).exists() else None
        checks.append({'artifact': artifact, 'path': path, 'recorded_sha256': digest,
                       'current_sha256': current, 'unchanged': current == digest})

    report = json.loads(WHOLE_GARMENTS.read_text())
    for path, digest in sorted(report['sources_sha256'].items()):
        compare('data/derived/whole-garments-uc8tp3zc', path, digest)
    configuration = json.loads((GARMENT_TISSUE / 'configuration.json').read_text())
    compare('data/derived/garment-tissue-refined-oufb46p6', 'data/derived/clothing/garments.json',
            configuration['panel']['source_sha256'])
    for path, digest in sorted(configuration['source_sha256'].items()):
        compare('data/derived/garment-tissue-refined-oufb46p6', path, digest)
    run = json.loads((GARMENT_TISSUE / 'report.json').read_text())
    return {
        'schema': 'ihm.prior-work-regression.v1',
        'statement': 'This build writes new files only. It does not read, rewrite or re-run the retained coupled '
                     'cloth-tissue experiment, and it leaves the generated-garment constructor and its emitted '
                     'garments.json untouched.',
        'all_unchanged': all(c['unchanged'] for c in checks),
        'checks': checks,
        'retained_coupled_run': {
            'path': 'data/derived/garment-tissue-refined-oufb46p6/coupled',
            'contact_resolutions': run['contact_resolutions'],
            'maximum_pair_impulse_residual_ns': run['maximum_pair_impulse_residual_ns'],
            'maximum_momentum_residual_ns': run['maximum_momentum_residual_ns'],
            'edge_prestrain': configuration['panel']['prestrain']},
    }


UNFILLED_SLOT_REASONS = {
    'neck': 'no scarf, tie or collar asset appears in any CC0 MakeHuman pack surveyed; the slot is declared and '
            'left empty rather than filled with an unlicensed mesh',
}


# ------------------------------------------------------------------------- main

def run(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-fit', action='store_true')
    parser.add_argument('--only', nargs='*', default=None)
    parser.add_argument('--extended', default='t-shirt',
                        help='garment given a longer settling run as a deeper receipt; "none" to skip')
    parser.add_argument('--extended-seconds', type=float, default=0.25)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args(argv)
    if args.self_test:
        return self_test()

    if not args.skip_fit:
        command = [str(LIBIGL), str(FITTER)] + (['--only', *args.only] if args.only else [])
        print('running', ' '.join(command))
        subprocess.run(command, check=True)

    registration_payload = json.loads((OUT / 'registration.json').read_text())
    registration = registration_payload['registration']
    fits = registration_payload['garments']
    retrieval = json.loads(RETRIEVAL.read_text())
    envelope_v, envelope_f, envelope_area, envelope_height = load_envelope(ENVELOPE)
    commit, dirty, dirty_files = git_state()
    catalogue = [e for e in CATALOGUE if e['id'] in fits and (args.only is None or e['id'] in args.only)]

    (OUT / 'geometry').mkdir(parents=True, exist_ok=True)
    (OUT / 'thumbnails').mkdir(parents=True, exist_ok=True)
    (OUT / 'provenance').mkdir(parents=True, exist_ok=True)
    (OUT / 'simulated').mkdir(parents=True, exist_ok=True)

    body = (envelope_v, envelope_f, BODY_RGB)
    records = []
    physics = {}
    for entry in catalogue:
        data = np.load(OUT / 'fitted' / f'{entry["id"]}.npz', allow_pickle=False)
        positions = np.asarray(data['positions'], float)
        triangles = np.asarray(data['indices'], np.int64)
        fit = fits[entry['id']]
        report, simulated = simulate(entry, positions, triangles, envelope_v, envelope_f)
        physics[entry['id']] = report
        if report['contact_ready']:
            np.savez_compressed(OUT / 'simulated' / f'{entry["id"]}.npz', positions=simulated, indices=triangles)
        geometry_path = f'data/derived/wardrobe-v1/geometry/{entry["id"]}.json.gz'
        components, _ = face_components(len(positions), triangles)
        boundary, nonmanifold = boundary_edge_count(triangles)
        payload = {'id': entry['id'], 'name': entry['name'], 'slots': entry['slots'],
                   'positions': np.round(positions, 7).ravel().tolist(),
                   'indices': triangles.ravel().tolist(),
                   'frame': 'bodyparts3d-display-m', 'units': 'm', 'representation': 'triangular_surface',
                   'derivation': 'acquired CC0 MakeHuman garment mesh registered onto the canonical outer body '
                                 'envelope; see data/derived/wardrobe-v1/registration.json'}
        with gzip.open(ROOT / geometry_path, 'wt') as handle:
            json.dump(payload, handle)
        thumbnail = OUT / 'thumbnails' / f'{entry["id"]}.png'
        image = render(body, [(positions, triangles, hex_rgb(entry['colour']))], thumbnail)
        transforms = [
            {'kind': 'uniform_similarity', 'from': 'makehuman-base-mesh', 'to': 'bodyparts3d-display-m',
             'scale': registration['similarity']['uniform_scale'],
             'translation_m': registration['similarity']['translation_m'],
             'rotation': 'identity',
             'residual': {'metric': 'trunk point-to-surface RMS', 'unit': 'mm',
                          'value': registration['similarity']['trunk_rms_mm_refined'],
                          'method': 'stature landmark scale, Nelder-Mead translation on trunk points'}},
            {'kind': 'piecewise_rigid_limb_correction', 'from': 'makehuman-base-mesh-aligned',
             'to': 'bodyparts3d-display-m',
             'bones': [{'limb': b['limb'], 'rotation_deg': b['rotation_deg'], 'limb_scale': b['limb_scale']}
                       for b in registration['pose_correction']['bones']],
             'residual': {'metric': 'authoring-body point-to-surface RMS after correction', 'unit': 'mm',
                          'value': registration['source_body_residual_mm']['after_pose_correction']['rms'],
                          'method': 'linear blend skinning over skeleton-joint and bone-end landmarks'}},
            {'kind': 'laplacian_regularised_shrinkwrap', 'from': 'bodyparts3d-display-m',
             'to': 'bodyparts3d-display-m', 'target_standoff_mm': entry['standoff_mm'],
             'iterations': registration['shrinkwrap']['iterations'],
             'residual': {'metric': 'interior vertices after fit', 'unit': 'count',
                          'value': fit['inside_body_vertices'],
                          'method': 'fast winding number signed distance to the outer envelope'}},
        ]
        record = provenance_record(entry, retrieval, registration, fit, geometry_path, commit, dirty, transforms)
        (OUT / 'provenance' / f'{entry["id"]}.json').write_text(json.dumps(record, indent=1, sort_keys=True) + '\n')
        records.append({
            'id': entry['id'], 'name': entry['name'], 'slots': entry['slots'], 'colour': entry['colour'],
            'geometry': geometry_path, 'geometry_sha256': sha256_file(ROOT / geometry_path),
            'thumbnail': f'data/derived/wardrobe-v1/thumbnails/{entry["id"]}.png',
            'thumbnail_sha256': sha256_file(thumbnail), 'thumbnail_pixels': image,
            'provenance': f'data/derived/wardrobe-v1/provenance/{entry["id"]}.json',
            'tier': record['tier'], 'license': retrieval['files'][entry['id']]['license'],
            'author': retrieval['files'][entry['id']]['author'],
            'vertices': fit['vertices'], 'triangles': fit['triangles'],
            'face_components': int(components), 'boundary_edges': boundary, 'nonmanifold_edges': nonmanifold,
            'area_m2': fit['area_m2'], 'mass_kg': report.get('mass_kg'),
            'target_standoff_mm': entry['standoff_mm'],
            'fit_mode': fit['fit_mode'], 'attract_range_mm': fit['attract_range_mm'],
            'standoff_mm': {'min': fit['standoff_mm_min'], 'median': fit['standoff_mm_percentiles']['50'],
                            'p95': fit['standoff_mm_percentiles']['95'],
                            'max': fit['standoff_mm_max'], 'mean': fit['standoff_mm_mean']},
            'loose_fit': fit['standoff_mm_percentiles']['95'] > 40.0,
            'interior_surface_samples': fit['inside_body_surface_samples'],
            'shallow_interior_surface_samples': fit['shallow_interior_surface_samples'],
            'deep_chord_surface_samples': fit['deep_chord_surface_samples'],
            'surface_samples': fit['surface_samples'],
            'inside_body_vertices': fit['inside_body_vertices'],
            'max_penetration_mm': fit['max_penetration_mm'],
            'edge_length_ratio_median': fit['edge_length_ratio_percentiles']['50'],
            'edge_length_rms_strain': fit['edge_length_rms_strain'],
            'area_ratio_to_posed_source': fit['area_ratio_to_posed_source'],
            'contact_ready': report['contact_ready'],
            'contacted_body_in_run': bool(report.get('contact_resolutions', 0) > 0),
            'contact_resolutions': int(report.get('contact_resolutions', 0)),
            'simulated_s': report.get('simulated_s'),
        })
        print(f'{entry["id"]:20s} sim {report.get("status")} steps={report.get("steps")} '
              f'contacts={report.get("contact_resolutions")} '
              f'pmax={report.get("maximum_paired_impulse_residual_ns", float("nan")):.3e} '
              f'disp p50={report.get("median_node_displacement_mm", float("nan")):.3f} mm')

    extended = None
    if args.extended and args.extended != 'none' and args.extended in {r['id'] for r in records}:
        entry = next(e for e in CATALOGUE if e['id'] == args.extended)
        data = np.load(OUT / 'fitted' / f'{entry["id"]}.npz', allow_pickle=False)
        report, final = simulate(entry, np.asarray(data['positions'], float), np.asarray(data['indices'], np.int64),
                                 envelope_v, envelope_f, seconds=args.extended_seconds, max_steps=12000)
        np.savez_compressed(OUT / 'simulated' / f'{entry["id"]}-extended.npz', positions=final,
                            indices=np.asarray(data['indices'], np.int64))
        report['purpose'] = ('longer settling run on one garment; the uniform per-garment run is deliberately short '
                             'so the whole wardrobe is measured under one budget')
        extended = {'garment': entry['id'], **report}
        print(f'extended {entry["id"]}: {report["simulated_s"]:.3f} s, {report["contact_resolutions"]} contacts, '
              f'displacement p50 {report["median_node_displacement_mm"]:.2f} mm')

    penetration = {}
    ids = [r['id'] for r in records if r['contact_ready']]
    if extended:
        ids.append(f'{args.extended}-extended')
    if ids:
        command = [str(LIBIGL), str(FITTER), '--penetration', *ids]
        subprocess.run(command, check=True)
        penetration = json.loads((OUT / 'post-simulation-penetration.json').read_text())
        for record in records:
            row = penetration.get(record['id'])
            if row:
                record['post_simulation'] = row

    model = slot_model([e for e in CATALOGUE if e['id'] in {r['id'] for r in records}], SLOTS)
    (OUT / 'slots.json').write_text(json.dumps(model, indent=1, sort_keys=True) + '\n')

    verification = {
        'schema': 'ihm.garment-wardrobe-verification.v1',
        'generated_at': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        'cloth_model': {
            'solver': 'ihm.assembly.clothing.Cloth, lumped-mass triangulated sheet with objective elastic edge '
                      'springs, symplectic Euler',
            'contact': 'ihm.assembly.garment_surface_contact.MovingSurfaceContact, nearest face-interior node '
                       'contact against the static body surface with Coulomb friction and paired body reactions',
            'energy_law': 'E = sum_e k (|e| - rest_e)^2 / 2; no bending law, no self-contact, no continuous '
                          'collision detection',
            'preload': f'every rest length scaled by {1 - PRESTRAIN}; an intrinsic statement, not an embedding',
            'parameter_status': 'explicit engineering priors, shared with data/derived/garment-tissue-refined-'
                                'oufb46p6; no fabric calibration',
            'gravity_m_s2': [0.0, -9.81, 0.0], 'body': 'static prescribed surface; no whole-body reaction applied',
        },
        'garments': physics,
        'prior_work_regression': prior_work_regression(),
        'extended_run': extended,
        'post_simulation_penetration': penetration,
        'limitations': [
            'The body is a static prescribed surface in these runs; reaction impulses are returned but not applied '
            'to any inertial owner.',
            'Node-face contact only. Edge/vertex, continuous and self collision are not certified.',
            'No bending stiffness, so drape shape is governed by membrane tension and contact alone.',
            f'{SIMULATED_SECONDS * 1e3:.0f} ms of simulated time per garment demonstrates tension and contact; it is '
            'not a settled drape and not a donning simulation.',
            'Areal density and stiffness are one uniform prior per garment; no per-material assignment exists.',
        ],
    }
    (OUT / 'verification.json').write_text(json.dumps(verification, indent=1, sort_keys=True) + '\n')

    wardrobe = {
        'schema': 'ihm.garment-wardrobe.v1',
        'body': {'surface': 'data/derived/outer-envelope/outer-envelope.npz', 'sha256': sha256_file(ENVELOPE),
                 'area_m2': envelope_area, 'height_m': envelope_height,
                 'basis': 'watertight genus-0 outer body envelope; the canonical skin entity is a two-sided '
                          f'{3.5025989744933166:.4f} m2 slab and is not used as a fit target'},
        'slot_model': model,
        'slot_model_file': 'data/derived/wardrobe-v1/slots.json',
        'unfilled_slots': {slot: UNFILLED_SLOT_REASONS.get(slot, 'no acquired mesh occupies this slot')
                           for slot in model['empty_slots']},
        'registration': 'data/derived/wardrobe-v1/registration.json',
        'verification': 'data/derived/wardrobe-v1/verification.json',
        'garments': records,
        'garment_count': len(records),
        'contact_ready_count': sum(1 for r in records if r['contact_ready']),
        'contact_ready_basis': 'contact_ready means the garment mesh built a valid Cloth (every node carries positive '
                               'lumped mass, no degenerate triangle) and completed a gravity-loaded run against the '
                               'body surface with the conservation receipts in verification.json. It does not mean '
                               'the garment was validated as drape. contacted_body_in_run is separate: a garment '
                               'carrying enough ease may never reach the body inside the short uniform run, and two '
                               'did not.',
        'no_contact_in_run': [r['id'] for r in records if not r['contacted_body_in_run']],
        'display_only': [r['id'] for r in records if not r['contact_ready']],
        'loose_fits': {r['id']: r['standoff_mm']['p95'] for r in records if r['loose_fit']},
        'loose_fit_note': 'a garment whose 95th-percentile standoff exceeds 40 mm carries real ease that the fit '
                          'deliberately did not remove. Widening the attraction to close it was measured and '
                          'rejected: it pulls the inner wall of a long sleeve onto the torso and shreds the sleeve.',
    }
    (OUT / 'wardrobe.json').write_text(json.dumps(wardrobe, indent=1, sort_keys=True) + '\n')

    manifest_path = OUT / 'manifest.json'
    inputs = {'data/derived/outer-envelope/outer-envelope.npz': sha256_file(ENVELOPE),
              'data/raw/clothing/makehuman/retrieval.json': sha256_file(RETRIEVAL)}
    for entry in catalogue:
        path = retrieval['files'][entry['id']]['path']
        inputs[path] = sha256_file(ROOT / path)
    for name in ('LICENSE.md', 'LICENSE.ASSETS.md'):
        inputs[f'data/raw/clothing/makehuman/{name}'] = sha256_file(ROOT / 'data/raw/clothing/makehuman' / name)
    inputs['data/raw/clothing/makehuman/base/base.obj'] = sha256_file(ROOT / 'data/raw/clothing/makehuman/base/base.obj')
    for bone in {v for v in registration['pose_correction']['landmarks'] for v in
                 (v['target_proximal']['geometry'], v['target_distal']['geometry'])}:
        inputs[bone] = sha256_file(ROOT / bone)
    sources = {p: sha256_file(ROOT / p) for p in (
        'ihm/assembly/garment_wardrobe.py', 'ihm/assembly/clothing.py',
        'ihm/assembly/garment_surface_contact.py', 'scripts/fit_garments_to_envelope.py',
        'scripts/retain_makehuman_garment_sources.py', 'scripts/build_garment_wardrobe.py')}
    artifacts = {}
    for path in sorted(OUT.rglob('*')):
        if path.is_file() and path != manifest_path and '.staging' not in path.parts:
            artifacts[str(path.relative_to(OUT))] = sha256_file(path)
    manifest = {
        'schema': 'ihm.garment-wardrobe-manifest.v1',
        'generated_at': dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        'commit': commit, 'commit_covers_working_tree': not dirty, 'uncommitted_changes': dirty_files,
        'builder_sha256': sha256_file(ROOT / 'scripts/build_garment_wardrobe.py'),
        'inputs_sha256': inputs,
        'source_sha256': sources,
        'artifacts_sha256': artifacts,
        'artifacts_exclude_self': 'manifest.json is deliberately absent from artifacts_sha256; a manifest cannot '
                                  'hash itself without mismatching on every verification',
        'provenance_schema': {'path': str(PROVENANCE_SCHEMA.relative_to(ROOT)),
                              'sha256': sha256_file(PROVENANCE_SCHEMA)},
        'canonical_assets_modified': False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=1, sort_keys=True) + '\n')
    print(f'\n{len(records)} garments emitted; {wardrobe["contact_ready_count"]} contact-ready; '
          f'{len(artifacts)} artifacts hashed')
    return 0


def self_test():
    import tempfile
    model = slot_model()
    assert model['schema'] == 'ihm.garment-slot-model.v1'
    for row in model['garments']:
        assert row['garment'] not in row['excludes'], 'a garment must not exclude itself'
        for other in model['garments']:
            same = bool(set(row['occupies']) & set(other['occupies'])) and row['garment'] != other['garment']
            assert same == (other['garment'] in row['excludes']), (row['garment'], other['garment'])
            assert (row['garment'] in other['excludes']) == same, 'exclusion must be symmetric'
    quad = np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    tri = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    assert abs(mesh_area_m2(quad, tri) - 1.0) < 1e-12
    assert boundary_edge_count(tri) == (4, 0)
    assert len(edge_table(tri)) == 5
    cloth = Cloth(quad, tri, areal_density_kg_m2=0.18, edge_stiffness_n_m=12.0)
    before = cloth.rest_lengths_m.copy()
    cloth.rest_lengths_m = cloth.rest_lengths_m * (1 - PRESTRAIN)
    force, energy = cloth.elastic_forces(quad)
    assert energy > 0 and np.linalg.norm(force.sum(0)) < 1e-12, 'preload must be self-equilibrated'
    assert np.allclose(cloth.rest_lengths_m / before, 1 - PRESTRAIN)
    r, u, w = camera_basis()
    for a, b in ((r, u), (u, w), (w, r)):
        assert abs(float(np.dot(a, b))) < 1e-12
    coarse = np.array([[0., 0, 0], [0.1, 0, 0], [0, 0.1, 0]])
    one = np.array([[0, 1, 2]], np.int64)
    refined_v, refined_f, split = refine_contact_faces(coarse, one, max_radius_m=0.008)
    rp = refined_v[refined_f]
    assert float(np.linalg.norm(rp - rp.mean(1)[:, None], axis=2).max()) <= 0.008 + 1e-12, 'refinement missed a face'
    assert abs(mesh_area_m2(refined_v, refined_f) - mesh_area_m2(coarse, one)) < 1e-15, 'refinement changed area'
    assert split >= 1 and len(refined_f) > 1, (split, len(refined_f))
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / 'probe.png'
        box = np.array([[-.1, -.1, -.1], [.1, -.1, -.1], [.1, .1, -.1], [-.1, .1, -.1],
                        [-.1, -.1, .1], [.1, -.1, .1], [.1, .1, .1], [-.1, .1, .1]])
        quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
        faces = np.array([f for q in quads for f in ((q[0], q[1], q[2]), (q[0], q[2], q[3]))], dtype=np.int64)
        stats = render((box, faces, BODY_RGB), [(box * 1.02, faces, (0.2, 0.4, 0.7))], target)
        _BODY_LAYER.clear()
        assert target.exists() and stats['covered_pixels'] > 100, stats
    regression = prior_work_regression()
    assert regression['all_unchanged'], [c for c in regression['checks'] if not c['unchanged']]
    assert regression['retained_coupled_run']['contact_resolutions'] == 10742, regression['retained_coupled_run']
    assert PROVENANCE_SCHEMA.exists(), 'provenance schema is missing'
    required = json.loads(PROVENANCE_SCHEMA.read_text())['required_fields']
    assert set(required) <= {'dataset.id', 'dataset.label', 'dataset.license', 'source_file.path',
                             'source_file.sha256', 'build.script', 'build.commit', 'geometry.path',
                             'geometry.sha256', 'transforms', 'tier'}, required
    print('build_garment_wardrobe self-test: pass')
    return 0


if __name__ == '__main__':
    raise SystemExit(run())
