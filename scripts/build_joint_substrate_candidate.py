"""Promotable candidate: 630 registered, repaired, TetGen-proved joint structures from the extended atlas.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_joint_substrate_candidate.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_joint_substrate_candidate.py \
  --output data/derived/joint-substrate-candidate-v1 --stage all

Nothing under data/derived/canonical is read for writing and no canonical file is touched. The result is
a candidate for review, not a promotion.

Stages, each resumable and each leaving its own receipt under stages/:
  registration  reproduce registrations.z_anatomy from its own recorded landmarks before applying it
  build         rotate + affine + thin-plate-spline into bodyparts3d-display-m, then the repair order
                imported verbatim from scripts/build_muscle_tet_ready_surfaces.py
  tetgen        forked stdout-capturing harness from scripts/verify_muscle_tet_ready_surfaces.py, run on
                the welded control and on the repaired surface, gated on accounted volume
  placement     signed distance to the repaired canonical bone union, plus name-derived expected bones
                and a duplicate screen against entities canonical already holds
  finalize      entities.jsonl, summary.json, manifest.json

No step fills a hole. An open sleeve is an open sleeve: it is reported as such and left open.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
import re
import shutil
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import igl
import igl.copyleft.cgal as cgal
from igl.copyleft import tetgen

from ihm.assembly.anatomy import LandmarkRegistration, ROTATION
from audit_joint_substrate import joint_class
from build_muscle_tet_ready_surfaces import (sha, write, normalize, repair, diagnose, signed_volume,
                                             weld_exact, drop_repeated_index, drop_zero_area)
from verify_muscle_tet_ready_surfaces import tetrahedralize, welded_control, VOLUME_GATE

CANONICAL = 'data/derived/canonical/anatomy.json'
SOURCE_INDEX = 'data/derived/anatomy/extended/source_index.json'
BONE_BUILD = 'data/derived/entity-tet-ready-v1'
SCHEMA = 'ihm.joint-substrate-candidate.v1'
TETGEN_FLAGS = 'pYq1.414'
PLACEMENT_SAMPLES = 3000
EMPTY_SPACE_M = .010
DEEP_PENETRATION_M = .002
DUPLICATE_SURFACE_M = .003
DUPLICATE_CENTROID_M = .020

# Name fragment -> generic canonical bone names. Side comes from the .l/.r suffix; an unsided
# structure gets both sides. Inferred association, not a curated attachment table.
BONE_TOKENS = (
    (r'femor|femur', ('femur',)),
    (r'tibio|tibial|tibia', ('tibia',)),
    (r'fibular|fibula|peroneal', ('fibula',)),
    (r'patell', ('patella',)),
    (r'humer', ('humerus',)),
    (r'radiocarpal', ('radius', 'scaphoid', 'lunate')),
    (r'radio|radial|radius', ('radius',)),
    (r'ulnar|ulna', ('ulna',)),
    (r'scapul|glenoid|glenohumeral|coracoid|coraco|acromio', ('scapula',)),
    (r'clavic', ('clavicle',)),
    (r'sternoclavicular', ('manubrium', 'clavicle')),
    (r'sterno|sternal|sternum|xiphoid', ('manubrium', 'body of sternum', 'xiphoid process')),
    (r'ilio|iliac|ischio|ischial|ischium|pubo|pubic|acetabul|obturator|inguinal|hip joint|coxal',
     ('hip bone',)),
    (r'sacro|sacral|sacrum', ('sacrum',)),
    (r'calcane', ('calcaneus',)),
    (r'talo|talar|talus', ('talus',)),
    (r'navicul', ('navicular bone of left foot', 'navicular bone of right foot')),
    (r'cuboid', ('cuboid bone',)),
    (r'cuneonavicular|cuneiform bone|intercuneiform|cuneometatarsal',
     ('medial cuneiform bone', 'intermediate cuneiform bone', 'lateral cuneiform bone')),
    (r'metatars', tuple('%s metatarsal bone' % n for n in
                        ('first', 'second', 'third', 'fourth', 'fifth'))),
    (r'metacarp', tuple('%s metacarpal bone' % n for n in
                        ('first', 'second', 'third', 'fourth', 'fifth'))),
    (r'scapho', ('scaphoid',)),
    (r'lunate', ('lunate',)),
    (r'triquetr', ('triquetral',)),
    (r'pisiform|pisi', ('pisiform',)),
    (r'trapezium', ('trapezium',)),
    (r'trapezoid', ('trapezoid',)),
    (r'capitate', ('capitate',)),
    (r'hamate|hamulus', ('hamate',)),
    (r'temporomandibular', ('mandible', 'temporal bone')),
    (r'mandib', ('mandible',)),
    (r'maxill', ('maxilla',)),
    (r'stylo|temporal bone|petro|tympan', ('temporal bone',)),
    (r'zygomat', ('zygomatic bone',)),
    (r'sphenoid|pterygo', ('sphenoid bone',)),
    (r'occipit|nuchal|foramen magnum', ('occipital bone',)),
    (r'frontal bone|supraorbital', ('frontal bone',)),
    (r'parietal', ('parietal bone',)),
    (r'nasal bone|nasomaxillary', ('nasal bone',)),
    (r'palatin', ('palatine bone',)),
    (r'lacrimal', ('lacrimal bone',)),
    (r'hyoid|hyo-|thyrohyoid|stylohyoid', ('hyoid bone',)),
    (r'atlanto|atlas', ('atlas',)),
    (r'\baxis\b|odontoid|dens', ('axis',)),
    (r'knee', ('femur', 'tibia', 'patella')),
    (r'elbow', ('humerus', 'ulna', 'radius')),
    (r'shoulder', ('scapula', 'humerus')),
    (r'wrist|carpal tunnel', ('radius', 'ulna', 'scaphoid', 'lunate', 'triquetral')),
    (r'ankle|talocrural', ('tibia', 'fibula', 'talus')),
    (r'costotransverse|costovertebral|costal|intercostal|\brib\b', ()),
)
CARPALS = ('scaphoid', 'lunate', 'triquetral', 'pisiform', 'trapezium', 'trapezoid', 'capitate', 'hamate')
METACARPALS = tuple('%s metacarpal bone' % n for n in ('first', 'second', 'third', 'fourth', 'fifth'))
METATARSALS = tuple('%s metatarsal bone' % n for n in ('first', 'second', 'third', 'fourth', 'fifth'))
CUNEIFORMS = ('medial cuneiform bone', 'intermediate cuneiform bone', 'lateral cuneiform bone')
NAVICULARS = ('navicular bone of left foot', 'navicular bone of right foot')

# Z-Anatomy source collection -> generic canonical bone names. The collection names the joint the
# structure was authored under, which is stronger evidence than wording alone.
JOINT_COLLECTIONS = {
    'Knee joint': ('femur', 'tibia', 'patella'),
    'Hip joint': ('hip bone', 'femur'),
    'Iliofemoral ligament': ('hip bone', 'femur'),
    'Sacro-iliac joint': ('sacrum', 'hip bone'),
    'Sacrotuberous ligament': ('sacrum', 'hip bone'),
    'Pubic symphysis': ('hip bone',),
    'Obturator membrane': ('hip bone',),
    'Elbow joint': ('humerus', 'ulna', 'radius'),
    'Glenohumeral joint': ('scapula', 'humerus'),
    'Glenohumeral ligaments': ('scapula', 'humerus'),
    'Acromioclavicular joint': ('scapula', 'clavicle'),
    'Coracoclavicular ligament': ('scapula', 'clavicle'),
    'Sternoclavicular joint': ('manubrium', 'clavicle'),
    'Distal radio-ulnar joint': ('radius', 'ulna'),
    'Radio-ulnar syndesmoses': ('radius', 'ulna'),
    'Radiocarpal joint': ('radius', 'ulna') + CARPALS[:3],
    'Palmar radiocarpal ligament': ('radius',) + CARPALS[:3],
    'Palmar ulnocarpal ligament': ('ulna', 'lunate', 'triquetral'),
    'Intercarpal joints': CARPALS,
    'Palmar intercarpal ligaments': CARPALS,
    'Dorsal intercarpal ligaments': CARPALS,
    'Interosseus intercarpal ligaments': CARPALS,
    'Pisiform joint': ('pisiform', 'triquetral'),
    'Carpometacarpal joints': CARPALS + METACARPALS,
    'Intermetacarpal joints': METACARPALS,
    'Metacarpophalangeal joints': METACARPALS,
    'Superior tibiofibular joint': ('tibia', 'fibula'),
    'Tibiofibular syndesmosis': ('tibia', 'fibula'),
    'Tibial collateral ligament': ('tibia', 'femur'),
    'Subtalar joint': ('talus', 'calcaneus'),
    'Talocalcaneonavicular joint': ('talus', 'calcaneus') + NAVICULARS,
    'Transverse tarsal joint': ('talus', 'calcaneus', 'cuboid bone') + NAVICULARS,
    'Calcaneocuboid joint': ('calcaneus', 'cuboid bone'),
    'Cuboidonavicular joint': ('cuboid bone',) + NAVICULARS,
    'Cuneonavicular joint': CUNEIFORMS + NAVICULARS,
    'Cuneocuboid joint': CUNEIFORMS + ('cuboid bone',),
    'Intercuneiform joints': CUNEIFORMS,
    'Bifurcate ligament': ('calcaneus', 'cuboid bone') + NAVICULARS,
    'Tarsometatarsal joints': CUNEIFORMS + ('cuboid bone',) + METATARSALS,
    'Intermetatarsal joints': METATARSALS,
    'Metatarsophalangeal joints': METATARSALS,
    'Temporomandibular joint': ('mandible', 'temporal bone'),
    'Thyrohyoid membrane': ('hyoid bone',),
}
SPINE_LEVEL = re.compile(r'\b(C|T|L)(\d{1,2})\b')
# A bursa or sheath is named for the muscle it serves, not for the bone it lies on, so wording
# tokens are muscle names there and only the authoring collection may speak.
EPONYMOUS_MUSCLE = re.compile(r'bursa|sheath|\bmuscle\b|tendinous|tendon of')
SPINE_REGION = {'C': 'cervical', 'T': 'thoracic', 'L': 'lumbar'}
ORDINALS = 'first second third fourth fifth sixth seventh eighth ninth tenth eleventh twelfth'.split()


def read_source(path):
    with np.load(path) as mesh:
        return mesh['vertices'].astype(float), np.ascontiguousarray(mesh['faces'].astype(np.int64))


def read_geometry(path):
    payload = json.loads(gzip.decompress(Path(path).read_bytes()))
    return (np.ascontiguousarray(np.asarray(payload['positions'], float).reshape(-1, 3)),
            np.ascontiguousarray(np.asarray(payload['indices'], np.int64).reshape(-1, 3)))


def side_of(name):
    if name.endswith('.l'):
        return 'left'
    if name.endswith('.r'):
        return 'right'
    return None


def expected_bones(name, collections, bone_names):
    """Attachment expectation inferred from structure wording and authoring collection.

    Inferred, not a curated attachment table: absence of an expectation is not evidence of
    misplacement, and a wrong token match would produce a wrong expectation.
    """
    lowered = name.lower()
    side = side_of(name)
    generic, matched = [], []
    for collection in collections or ():
        if collection in JOINT_COLLECTIONS:
            matched.append('collection:' + collection)
            generic.extend(JOINT_COLLECTIONS[collection])
    if not EPONYMOUS_MUSCLE.search(lowered):
        for pattern, bones in BONE_TOKENS:
            if re.search(pattern, lowered):
                matched.append(pattern)
                generic.extend(bones)
    for letter, number in SPINE_LEVEL.findall(name):
        index = int(number)
        if letter == 'C' and index == 1:
            generic.append('atlas')
        elif letter == 'C' and index == 2:
            generic.append('axis')
        elif 1 <= index <= 12:
            generic.append('%s %s vertebra' % (ORDINALS[index - 1], SPINE_REGION[letter]))
            matched.append('spine_level')
    sides = [side] if side else ['left', 'right']
    resolved = []
    for base in dict.fromkeys(generic):
        if base in bone_names:
            # An already-sided generic name (navicular, sesamoid) must still match the structure's side.
            if side and ('left ' in base or 'right ' in base) and side + ' ' not in base:
                continue
            candidates = (base,)
        else:
            candidates = tuple(s + ' ' + base for s in sides)
        for candidate in candidates:
            if candidate in bone_names and candidate not in resolved:
                resolved.append(candidate)
    return resolved, sorted(set(matched))


def stride(points, limit):
    return points[::max(1, len(points) // limit)]


# ---------------------------------------------------------------- stage: registration

def fit_registration(recorded):
    source = np.array([l['source_rotated_m'] for l in recorded['landmarks']])
    target = np.array([l['target_m'] for l in recorded['landmarks']])
    fit = LandmarkRegistration(source, target, smoothing=recorded['smoothing_prior'])
    return fit, source, target


def stage_registration(out, canonical):
    recorded = canonical['registrations']['z_anatomy']
    fit, source, target = fit_registration(recorded)
    report = fit.report()
    folds = np.arange(len(source)) % 5
    fold_fits = [LandmarkRegistration(source[folds != f], target[folds != f],
                                      smoothing=recorded['smoothing_prior']) for f in range(5)]
    held = np.concatenate([np.linalg.norm(fold_fits[f].transform(source[folds == f]) - target[folds == f], axis=1)
                           for f in range(5)])
    reproduced = {k: bool(abs(report[k] - recorded[k]) <= 1e-12)
                  for k in ('affine_rms_m', 'fit_rms_m', 'fit_max_m', 'affine_determinant')}
    reproduced['held_out_rms_m'] = bool(abs(float(np.sqrt(np.mean(held ** 2))) - recorded['held_out_rms_m']) <= 1e-12)
    reproduced['held_out_max_m'] = bool(abs(float(held.max()) - recorded['held_out_max_m']) <= 1e-12)
    reproduced['affine_4x3'] = bool(np.abs(np.asarray(report['affine_4x3']) -
                                           np.asarray(recorded['affine_4x3'])).max() <= 1e-12)
    reproduced['landmark_count'] = report['landmark_count'] == recorded['landmark_count']
    payload = dict(schema=SCHEMA + '.registration',
                   recorded={k: v for k, v in recorded.items()
                             if k not in ('landmarks', 'held_out_residuals_m', 'affine_4x3')},
                   recomputed={k: report[k] for k in ('affine_rms_m', 'fit_rms_m', 'fit_max_m',
                                                      'affine_determinant', 'landmark_count')},
                   recomputed_held_out_rms_m=float(np.sqrt(np.mean(held ** 2))),
                   recomputed_held_out_max_m=float(held.max()),
                   max_absolute_affine_coefficient_difference=float(np.abs(
                       np.asarray(report['affine_4x3']) - np.asarray(recorded['affine_4x3'])).max()),
                   reproduced_exactly=reproduced,
                   all_reproduced=bool(all(reproduced.values())),
                   never_fitted_bar_from_audit=dict(
                       source='data/derived/joint-substrate-assessment-v1/registration.json',
                       never_fitted_centroid_rms_m=.006180582966048701,
                       never_fitted_surface_rms_m=.004124695382644262,
                       rejected_precedent_centroid_residual_m=.036229),
                   interpretation='Reproducing the recorded fit proves the transform is applied as recorded. '
                                  'It is not independent validation of the registration itself.')
    if not payload['all_reproduced']:
        raise ValueError('Recorded registration did not reproduce; refusing to transform')
    write(out / 'stages/registration.json', payload)
    return fit


# ---------------------------------------------------------------- stage: build

_WORKER = {}


def _init_worker(out, affine, source_landmarks, smoothing, target_landmarks):
    _WORKER['out'] = Path(out)
    _WORKER['fit'] = LandmarkRegistration(np.asarray(source_landmarks), np.asarray(target_landmarks),
                                          smoothing=smoothing)
    assert np.abs(np.asarray(_WORKER['fit'].affine) - np.asarray(affine)).max() <= 1e-12


def _build_one(job):
    out, fit = _WORKER['out'], _WORKER['fit']
    started = time.monotonic()
    path = ROOT / job['source_geometry_path']
    digest = sha(path)
    if digest != job['source_geometry_sha256']:
        raise ValueError('Extended source geometry digest mismatch: ' + job['id'])
    v, f = read_source(path)
    if not np.isfinite(v).all() or f.min() < 0 or f.max() >= len(v):
        raise ValueError('Invalid source geometry: ' + job['id'])
    source_bounds = [v.min(0).tolist(), v.max(0).tolist()]
    moved = np.ascontiguousarray(fit.transform(v @ ROTATION.T))
    jacobian = fit.jacobian_determinants(stride(np.ascontiguousarray(v @ ROTATION.T), 400))
    raw = diagnose(moved, f, intersections=False)
    before = diagnose(*weld_exact(moved, f))
    rv, rf, log = repair(moved, f)
    after = diagnose(rv, rf)
    drift = after['signed_volume_m3'] - before['signed_volume_m3']
    payload = dict(schema=SCHEMA + '.surface', structure_id=job['id'], name=job['name'],
                   joint_class=job['joint_class'], units='m', frame='bodyparts3d-display-m',
                   representation='triangular_surface', source_path=job['source_geometry_path'],
                   source_sha256=digest, registration='registrations.z_anatomy',
                   positions=[float(x) for x in rv.ravel()], indices=[int(i) for i in rf.ravel()])
    geometry = out / 'geometry' / (job['id'] + '.json.gz')
    geometry.write_bytes(gzip.compress(json.dumps(payload, allow_nan=False).encode(), mtime=0))
    return dict(structure_id=job['id'], name=job['name'], joint_class=job['joint_class'],
                system=job['system'], in_joints_collection=job['in_joints_collection'],
                source_collections=job['source_collections'],
                source_path=job['source_geometry_path'], source_sha256=digest,
                source_bounds_blender_m=source_bounds,
                transform=dict(rotation='ROTATION', registration='registrations.z_anatomy',
                               jacobian_determinant_min=float(jacobian.min()),
                               jacobian_determinant_max=float(jacobian.max())),
                raw=raw, before=before, after=after, operations=log,
                output_path=str(geometry.relative_to(out)), output_sha256=sha(geometry),
                signed_volume_drift_m3=drift,
                signed_volume_relative_drift=drift / before['signed_volume_m3'] if before['signed_volume_m3'] else None,
                hole_filling_required=after['boundary_edges'] > 0,
                open_sheet_or_sleeve=after['boundary_edges'] > 0,
                encloses_volume=bool(after['boundary_edges'] == 0 and after['nonmanifold_edges'] == 0
                                     and after['signed_volume_m3'] > 0),
                vertex_links_manifold=after['nonmanifold_vertices'] == 0,
                tet_ready=bool(after['closed'] and after['self_intersection_free']
                               and after['orientation_consistent'] and after['signed_volume_m3'] > 0),
                repair_seconds=time.monotonic() - started)


def stage_build(out, canonical, index, fit, workers, limit):
    jobs = []
    for mesh in index['meshes']:
        label = joint_class(mesh['name'], mesh['system'])
        if label is None:
            continue
        jobs.append(dict(id=mesh['id'], name=mesh['name'], joint_class=label, system=mesh['system'],
                         in_joints_collection='3: Joints' in mesh['source_collections'],
                         source_collections=mesh['source_collections'],
                         source_geometry_path=mesh['source_geometry_path'],
                         source_geometry_sha256=mesh['source_geometry_sha256'],
                         triangles=mesh['source_triangles']))
    jobs.sort(key=lambda j: (-j['triangles'], j['id']))
    if limit:
        jobs = jobs[:limit]
    (out / 'geometry').mkdir(parents=True, exist_ok=True)
    recorded = canonical['registrations']['z_anatomy']
    args = (str(out), fit.affine.tolist(), [l['source_rotated_m'] for l in recorded['landmarks']],
            recorded['smoothing_prior'], [l['target_m'] for l in recorded['landmarks']])
    records, started = [], time.monotonic()
    path = out / 'stages/build-records.jsonl'
    path.parent.mkdir(parents=True, exist_ok=True)
    with mp.get_context('fork').Pool(workers, initializer=_init_worker, initargs=args) as pool:
        with path.open('w') as handle:
            for record in pool.imap_unordered(_build_one, jobs, chunksize=1):
                records.append(record)
                handle.write(json.dumps(record, allow_nan=False) + '\n')
                handle.flush()
                if len(records) % 25 == 0:
                    print('repaired %d/%d  %.0fs' % (len(records), len(jobs), time.monotonic() - started), flush=True)
    print('repaired %d/%d  %.0fs' % (len(records), len(jobs), time.monotonic() - started), flush=True)
    return sorted(records, key=lambda r: r['structure_id'])


# ---------------------------------------------------------------- stage: tetgen

def gate(outcome, surface_volume):
    """The harness's volume gate, made binding. A zero exit with tets that account for no volume fails."""
    tet_volume = outcome.get('tet_volume_m3')
    error = ((tet_volume - surface_volume) / surface_volume
             if tet_volume is not None and surface_volume else None)
    outcome['tet_volume_vs_surface_relative_error'] = error
    outcome['surface_signed_volume_m3'] = surface_volume
    outcome['succeeded_exit_and_tets_only'] = bool(outcome.get('status') == 0 and outcome.get('tets', 0) > 0)
    outcome['succeeded'] = bool(outcome['succeeded_exit_and_tets_only'] and surface_volume > 0
                                and error is not None and abs(error) <= VOLUME_GATE)
    return outcome


def _tetgen_one(record):
    out = _WORKER['out']
    payload = json.loads(gzip.decompress((out / record['output_path']).read_bytes()))
    rv = np.ascontiguousarray(np.asarray(payload['positions'], float).reshape(-1, 3))
    rf = np.ascontiguousarray(np.asarray(payload['indices'], np.int64).reshape(-1, 3))
    v, f = read_source(ROOT / record['source_path'])
    moved = np.ascontiguousarray(_WORKER['fit'].transform(v @ ROTATION.T))
    cv, cf = welded_control(moved, f)
    control = gate(tetrahedralize(cv, cf, TETGEN_FLAGS), signed_volume(cv, cf))
    repaired = gate(tetrahedralize(rv, rf, TETGEN_FLAGS), signed_volume(rv, rf))
    return dict(structure_id=record['structure_id'], name=record['name'],
                joint_class=record['joint_class'],
                faces_welded_control=int(len(cf)), faces_repaired=int(len(rf)),
                self_intersecting_pairs_before=record['before']['self_intersecting_face_pairs'],
                boundary_edges_after=record['after']['boundary_edges'],
                nonmanifold_edges_after=record['after']['nonmanifold_edges'],
                nonmanifold_vertices_after=record['after']['nonmanifold_vertices'],
                face_components_after=record['after']['face_components'],
                signed_volume_before_m3=record['before']['signed_volume_m3'],
                signed_volume_after_m3=record['after']['signed_volume_m3'],
                signed_volume_relative_drift=record['signed_volume_relative_drift'],
                tet_volume_vs_surface_relative_error=repaired['tet_volume_vs_surface_relative_error'],
                welded_control=control, repaired=repaired)


def stage_tetgen(out, records, fit, canonical, workers):
    recorded = canonical['registrations']['z_anatomy']
    args = (str(out), fit.affine.tolist(), [l['source_rotated_m'] for l in recorded['landmarks']],
            recorded['smoothing_prior'], [l['target_m'] for l in recorded['landmarks']])
    trials, started = [], time.monotonic()
    path = out / 'stages/tetgen-trials.jsonl'
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records, key=lambda r: (-r['after']['faces'], r['structure_id']))
    with mp.get_context('fork').Pool(workers, initializer=_init_worker, initargs=args) as pool:
        with path.open('w') as handle:
            for trial in pool.imap_unordered(_tetgen_one, ordered, chunksize=1):
                trials.append(trial)
                handle.write(json.dumps(trial, allow_nan=False) + '\n')
                handle.flush()
                if len(trials) % 25 == 0:
                    print('tetgen %d/%d  %.0fs' % (len(trials), len(ordered), time.monotonic() - started), flush=True)
    print('tetgen %d/%d  %.0fs' % (len(trials), len(ordered), time.monotonic() - started), flush=True)
    return sorted(trials, key=lambda t: t['structure_id'])


# ---------------------------------------------------------------- stage: placement

NOT_TRUE_BONE = re.compile(r'intervertebral disk|tooth')


def bone_union(canonical, build, entities=None):
    """Concatenate the already-repaired canonical bone shells into one winding-number domain."""
    vertices, faces, owner, names, meshes, offset = [], [], [], [], [], 0
    bones = entities if entities is not None else [e for e in canonical['entities']
                                                   if e['role'] == 'rigid_bone']
    for position, bone in enumerate(sorted(bones, key=lambda e: e['id'])):
        v, f = read_geometry(build / 'geometry' / (bone['id'] + '.json.gz'))
        vertices.append(v)
        faces.append(f + offset)
        owner.append(np.full(len(f), position, np.int64))
        names.append((bone['id'], bone['name']))
        meshes.append((v, f))
        offset += len(v)
    return (np.ascontiguousarray(np.concatenate(vertices)),
            np.ascontiguousarray(np.concatenate(faces)),
            np.concatenate(owner), names, meshes)


def canonical_neighbours(canonical):
    keep = {'rigid_bone', 'cartilage', 'ligament', 'tendon', 'connective_tissue', 'soft_organ'}
    rows = []
    for entity in canonical['entities']:
        if entity['role'] not in keep or entity['reference_geometry'].get('representation') != 'triangular_surface':
            continue
        rows.append(dict(id=entity['id'], name=entity['name'], role=entity['role'],
                         centroid=np.asarray(entity['centroid_m'], float),
                         low=np.asarray(entity['bounds_m']['min'], float),
                         high=np.asarray(entity['bounds_m']['max'], float),
                         path=entity['reference_geometry']['path']))
    return rows


def stage_placement(out, records, canonical, build):
    rigid = [e for e in canonical['entities'] if e['role'] == 'rigid_bone']
    true_bone = [e for e in rigid if not NOT_TRUE_BONE.search(e['name'])]
    bone_v, bone_f, bone_owner, bone_names, bone_meshes = bone_union(canonical, build, rigid)
    skeletal_v, skeletal_f, skeletal_owner, skeletal_names, _ = bone_union(canonical, build, true_bone)
    name_to_index = {name: i for i, (_, name) in enumerate(bone_names)}
    neighbours = canonical_neighbours(canonical)
    neighbour_centroids = np.stack([n['centroid'] for n in neighbours])
    neighbour_cache = {}
    ordered = sorted(records, key=lambda r: r['structure_id'])
    surfaces, samples, slices, cursor = {}, [], {}, 0
    for record in ordered:
        payload = json.loads(gzip.decompress((out / record['output_path']).read_bytes()))
        v = np.ascontiguousarray(np.asarray(payload['positions'], float).reshape(-1, 3))
        f = np.ascontiguousarray(np.asarray(payload['indices'], np.int64).reshape(-1, 3))
        surfaces[record['structure_id']] = (v, f)
        if not len(v):
            continue
        picked = np.ascontiguousarray(stride(v, PLACEMENT_SAMPLES))
        samples.append(picked)
        slices[record['structure_id']] = (cursor, cursor + len(picked))
        cursor += len(picked)
    points = np.ascontiguousarray(np.concatenate(samples))
    started = time.monotonic()
    print('signed distance: %d points against %d bone faces' % (len(points), len(bone_f)), flush=True)
    all_signed, all_facet, _, _ = igl.signed_distance(points, bone_v, bone_f,
                                                      igl.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
    true_signed, true_facet, _, _ = igl.signed_distance(points, skeletal_v, skeletal_f,
                                                        igl.SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER)
    print('signed distance done in %.0fs' % (time.monotonic() - started), flush=True)
    rows = []
    path = out / 'stages/placement.jsonl'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        for record in ordered:
            if record['structure_id'] not in slices:
                continue
            v, faces = surfaces[record['structure_id']]
            lo, hi = slices[record['structure_id']]
            structure_points = points[lo:hi]
            signed, facet = all_signed[lo:hi], all_facet[lo:hi]
            owner = bone_owner[facet]
            gaps = np.maximum(signed, 0.)
            nearest = {int(i): float(gaps[owner == i].min()) for i in np.unique(owner)}
            ranked = sorted(nearest.items(), key=lambda kv: kv[1])[:4]
            skeletal_signed = true_signed[lo:hi]
            deepest = int(np.argmin(skeletal_signed))
            deepest_bone = skeletal_names[int(skeletal_owner[true_facet[lo:hi][deepest]])][1]
            names, patterns = expected_bones(record['name'], record['source_collections'], name_to_index)
            expected = []
            for bone in names:
                index = name_to_index[bone]
                bv, bf = bone_meshes[index]
                squared, _, _ = igl.point_mesh_squared_distance(structure_points, bv, bf)
                distances = np.sqrt(np.maximum(squared, 0.))
                expected.append(dict(bone=bone, bone_id=bone_names[index][0],
                                     min_gap_m=float(distances.min()),
                                     median_gap_m=float(np.median(distances)),
                                     fraction_within_2mm=float((distances < .002).mean())))
            centroid = v.mean(0)
            low, high = v.min(0), v.max(0)
            near = np.flatnonzero(np.linalg.norm(neighbour_centroids - centroid, axis=1) < DUPLICATE_CENTROID_M)
            duplicates = []
            for index in near:
                other = neighbours[index]
                if (low > other['high']).any() or (high < other['low']).any():
                    continue
                if index not in neighbour_cache:
                    neighbour_cache[int(index)] = read_geometry(ROOT / other['path'])
                ov, of = neighbour_cache[int(index)]
                forward = np.sqrt(np.maximum(igl.point_mesh_squared_distance(
                    np.ascontiguousarray(stride(v, 600)), ov, of)[0], 0.))
                backward = np.sqrt(np.maximum(igl.point_mesh_squared_distance(
                    np.ascontiguousarray(stride(ov, 600)), v, faces)[0], 0.))
                both = np.concatenate([forward, backward])
                if float(both.mean()) <= DUPLICATE_SURFACE_M:
                    duplicates.append(dict(canonical_id=other['id'], canonical_name=other['name'],
                                           canonical_role=other['role'],
                                           symmetric_mean_surface_m=float(both.mean()),
                                           symmetric_max_surface_m=float(both.max()),
                                           centroid_distance_m=float(np.linalg.norm(other['centroid'] - centroid))))
            row = dict(structure_id=record['structure_id'], name=record['name'],
                       joint_class=record['joint_class'],
                       sample_points=int(len(points)),
                       centroid_m=centroid.tolist(),
                       min_gap_to_any_bone_m=float(gaps.min()),
                       median_gap_to_nearest_bone_m=float(np.median(gaps)),
                       fraction_within_2mm_of_bone=float((gaps < .002).mean()),
                       inside_rigid_bone_role_fraction=float((signed < 0).mean()),
                       max_rigid_bone_role_penetration_m=float(max(0., -signed.min())),
                       inside_true_bone_fraction=float((skeletal_signed < 0).mean()),
                       max_true_bone_penetration_m=float(max(0., -skeletal_signed.min())),
                       deepest_true_bone=deepest_bone,
                       nearest_bones=[dict(bone=bone_names[i][1], bone_id=bone_names[i][0], min_gap_m=g)
                                      for i, g in ranked],
                       expected_bones=expected,
                       expected_bone_source='name_and_collection_tokens' if names else 'none_derivable',
                       expected_bone_patterns=patterns,
                       min_expected_bone_min_gap_m=min((e['min_gap_m'] for e in expected), default=None),
                       max_expected_bone_min_gap_m=max((e['min_gap_m'] for e in expected), default=None),
                       duplicates_canonical=duplicates)
            row['flag_no_bone_contact'] = bool(row['min_gap_to_any_bone_m'] > EMPTY_SPACE_M)
            row['flag_reaches_no_expected_bone'] = bool(
                expected and row['min_expected_bone_min_gap_m'] > EMPTY_SPACE_M)
            row['flag_some_expected_bone_far'] = bool(
                expected and row['max_expected_bone_min_gap_m'] > EMPTY_SPACE_M)
            row['flag_deep_inside_bone'] = bool(row['inside_true_bone_fraction'] > .5
                                                and row['max_true_bone_penetration_m'] > DEEP_PENETRATION_M)
            row['flag_duplicates_canonical'] = bool(duplicates)
            rows.append(row)
            handle.write(json.dumps(row, allow_nan=False) + '\n')
            handle.flush()
            if len(rows) % 50 == 0:
                print('placement %d/%d  %.0fs' % (len(rows), len(records), time.monotonic() - started), flush=True)
    print('placement %d/%d  %.0fs' % (len(rows), len(records), time.monotonic() - started), flush=True)
    return rows


# ---------------------------------------------------------------- stage: finalize

def class_table(records, trials, placement):
    by_trial = {t['structure_id']: t for t in trials}
    by_place = {p['structure_id']: p for p in placement}
    table = {}
    for record in records:
        row = table.setdefault(record['joint_class'], dict(
            count=0, source_triangles=0, repaired_faces=0,
            before_closed=0, before_self_intersection_free=0, before_orientation_consistent=0,
            before_nonmanifold_edges=0, before_boundary_edges=0, before_self_intersecting_pairs=0,
            after_closed=0, after_self_intersection_free=0, after_nonmanifold_edges=0,
            after_boundary_edges=0, after_nonmanifold_vertex_links=0,
            open_after_repair=0, encloses_volume=0, tet_ready=0,
            arrangement_applied=0, self_union_applied=0, patches_flipped_outward=0,
            control_tetrahedralized=0, repaired_tetrahedralized=0,
            repaired_zero_exit_but_no_volume=0, tets=0,
            flag_no_bone_contact=0, flag_reaches_no_expected_bone=0, flag_some_expected_bone_far=0,
            flag_deep_inside_bone=0, flag_duplicates_canonical=0, expected_bones_derived=0))
        row['count'] += 1
        row['source_triangles'] += record['raw']['faces']
        row['repaired_faces'] += record['after']['faces']
        for stage in ('before', 'after'):
            row[stage + '_closed'] = row.get(stage + '_closed', 0) + bool(record[stage]['closed'])
        row['before_self_intersection_free'] += bool(record['before']['self_intersection_free'])
        row['before_orientation_consistent'] += bool(record['before']['orientation_consistent'])
        row['before_nonmanifold_edges'] += record['before']['nonmanifold_edges']
        row['before_boundary_edges'] += record['before']['boundary_edges']
        row['before_self_intersecting_pairs'] += record['before']['self_intersecting_face_pairs']
        row['after_self_intersection_free'] += bool(record['after']['self_intersection_free'])
        row['after_nonmanifold_edges'] += record['after']['nonmanifold_edges']
        row['after_boundary_edges'] += record['after']['boundary_edges']
        row['after_nonmanifold_vertex_links'] += record['after']['nonmanifold_vertices'] > 0
        row['open_after_repair'] += bool(record['open_sheet_or_sleeve'])
        row['encloses_volume'] += bool(record['encloses_volume'])
        row['tet_ready'] += bool(record['tet_ready'])
        row['arrangement_applied'] += any(s['stage'] == 'resolve_self_intersections' for s in record['operations'])
        row['self_union_applied'] += any(s['stage'] == 'extract_outer_manifold' for s in record['operations'])
        row['patches_flipped_outward'] += sum(s.get('patches_flipped_outward', 0) for s in record['operations'])
        trial = by_trial.get(record['structure_id'])
        if trial:
            row['control_tetrahedralized'] += bool(trial['welded_control']['succeeded'])
            row['repaired_tetrahedralized'] += bool(trial['repaired']['succeeded'])
            row['repaired_zero_exit_but_no_volume'] += bool(trial['repaired']['succeeded_exit_and_tets_only']
                                                            and not trial['repaired']['succeeded'])
            row['tets'] += trial['repaired'].get('tets', 0) if trial['repaired']['succeeded'] else 0
        place = by_place.get(record['structure_id'])
        if place:
            for flag in ('flag_no_bone_contact', 'flag_reaches_no_expected_bone',
                         'flag_some_expected_bone_far', 'flag_deep_inside_bone',
                         'flag_duplicates_canonical'):
                row[flag] += bool(place[flag])
            row['expected_bones_derived'] += bool(place['expected_bones'])
    return {k: table[k] for k in sorted(table)}


def rediagnose(out, records):
    """Recompute every stored diagnostic from the emitted file rather than trusting the build."""
    mismatches = []
    for record in records:
        payload = json.loads(gzip.decompress((out / record['output_path']).read_bytes()))
        v = np.ascontiguousarray(np.asarray(payload['positions'], float).reshape(-1, 3))
        f = np.ascontiguousarray(np.asarray(payload['indices'], np.int64).reshape(-1, 3))
        current = diagnose(v, f)
        differing = {k: (record['after'][k], current[k]) for k in current
                     if (abs(current[k] - record['after'][k]) > 1e-18 if k == 'signed_volume_m3'
                         else current[k] != record['after'][k])}
        if differing:
            mismatches.append(dict(structure_id=record['structure_id'], differing=differing))
    return mismatches


def stage_finalize(out, records, trials, placement, canonical_path, index_path, build, elapsed):
    mismatches = rediagnose(out, records)
    if mismatches:
        write(out / 'stages/diagnostic-mismatches.json', mismatches)
    by_trial = {t['structure_id']: t for t in trials}
    by_place = {p['structure_id']: p for p in placement}
    merged = []
    with (out / 'entities.jsonl').open('w') as handle:
        for record in sorted(records, key=lambda r: r['structure_id']):
            trial = by_trial.get(record['structure_id'])
            place = by_place.get(record['structure_id'])
            row = dict(record)
            row['tetgen'] = dict(flags=TETGEN_FLAGS, volume_gate=VOLUME_GATE,
                                 welded_control=trial['welded_control'] if trial else None,
                                 repaired=trial['repaired'] if trial else None) if trial else None
            row['placement'] = place
            row['promotable'] = bool(record['tet_ready'] and trial and trial['repaired']['succeeded']
                                     and place and not place['flag_no_bone_contact']
                                     and not place['flag_reaches_no_expected_bone']
                                     and not place['flag_deep_inside_bone']
                                     and not place['flag_duplicates_canonical'])
            merged.append(row)
            handle.write(json.dumps(row, allow_nan=False) + '\n')
    successes = [t for t in trials if t['repaired']['succeeded']]
    errors = [abs(t['tet_volume_vs_surface_relative_error']) for t in successes
              if t['tet_volume_vs_surface_relative_error'] is not None]
    drifts = [abs(r['signed_volume_relative_drift']) for r in records
              if r['signed_volume_relative_drift'] is not None]
    summary = dict(
        schema=SCHEMA + '.summary', structures=len(records),
        canonical_modified=False,
        stage_meaning=dict(
            raw='registered into the canonical frame, rendering topology, seam-duplicated vertices',
            before='after lossless exact-coordinate welding only; the defect baseline',
            after='after the full repair order imported from build_muscle_tet_ready_surfaces.py'),
        registration=json.loads((out / 'stages/registration.json').read_text())['reproduced_exactly'],
        recorded_diagnostics_reproduced_from_emitted_geometry=not mismatches,
        diagnostic_mismatches=len(mismatches),
        totals=dict(
            source_triangles=sum(r['raw']['faces'] for r in records),
            repaired_faces=sum(r['after']['faces'] for r in records),
            before_closed=sum(bool(r['before']['closed']) for r in records),
            after_closed=sum(bool(r['after']['closed']) for r in records),
            before_self_intersection_free=sum(bool(r['before']['self_intersection_free']) for r in records),
            after_self_intersection_free=sum(bool(r['after']['self_intersection_free']) for r in records),
            before_total_self_intersecting_pairs=sum(r['before']['self_intersecting_face_pairs'] for r in records),
            after_total_self_intersecting_pairs=sum(r['after']['self_intersecting_face_pairs'] for r in records),
            before_total_nonmanifold_edges=sum(r['before']['nonmanifold_edges'] for r in records),
            after_total_nonmanifold_edges=sum(r['after']['nonmanifold_edges'] for r in records),
            before_total_boundary_edges=sum(r['before']['boundary_edges'] for r in records),
            after_total_boundary_edges=sum(r['after']['boundary_edges'] for r in records),
            tet_ready=sum(bool(r['tet_ready']) for r in records),
            open_after_repair=sum(bool(r['open_sheet_or_sleeve']) for r in records),
            encloses_volume=sum(bool(r['encloses_volume']) for r in records),
            arrangement_applied=sum(any(s['stage'] == 'resolve_self_intersections' for s in r['operations'])
                                    for r in records),
            self_union_applied=sum(any(s['stage'] == 'extract_outer_manifold' for s in r['operations'])
                                   for r in records),
            self_union_failed=sum(any(s['stage'] == 'extract_outer_manifold' and 'failed' in s
                                      for s in r['operations']) for r in records)),
        tetgen=dict(
            flags=TETGEN_FLAGS, volume_gate=VOLUME_GATE, trials=len(trials),
            repaired_tetrahedralized=len(successes),
            welded_control_tetrahedralized=sum(t['welded_control']['succeeded'] for t in trials),
            repaired_zero_exit_but_volume_gate_failed=sum(
                t['repaired']['succeeded_exit_and_tets_only'] and not t['repaired']['succeeded'] for t in trials),
            control_zero_exit_but_volume_gate_failed=sum(
                t['welded_control']['succeeded_exit_and_tets_only'] and not t['welded_control']['succeeded']
                for t in trials),
            repaired_success_rate=len(successes) / len(trials) if trials else 0.,
            total_tets=sum(t['repaired'].get('tets', 0) for t in successes),
            median_tets=float(np.median([t['repaired']['tets'] for t in successes])) if successes else 0.,
            max_absolute_tet_vs_surface_volume_error=max(errors, default=0.),
            median_absolute_tet_vs_surface_volume_error=float(np.median(errors)) if errors else 0.,
            repaired_process_aborts=sum(t['repaired'].get('exception') == 'ProcessAborted' for t in trials),
            control_process_aborts=sum(t['welded_control'].get('exception') == 'ProcessAborted' for t in trials)),
        volume_drift=dict(
            max_absolute_relative=max(drifts, default=0.),
            median_absolute_relative=float(np.median(drifts)) if drifts else 0.,
            over_1pct=[r['structure_id'] for r in records if r['signed_volume_relative_drift'] is not None
                       and abs(r['signed_volume_relative_drift']) > .01]),
        placement=dict(
            structures=len(placement),
            with_name_derived_expected_bones=sum(bool(p['expected_bones']) for p in placement),
            flag_no_bone_contact=sum(p['flag_no_bone_contact'] for p in placement),
            flag_reaches_no_expected_bone=sum(p['flag_reaches_no_expected_bone'] for p in placement),
            flag_some_expected_bone_far=sum(p['flag_some_expected_bone_far'] for p in placement),
            flag_deep_inside_bone=sum(p['flag_deep_inside_bone'] for p in placement),
            flag_duplicates_canonical=sum(p['flag_duplicates_canonical'] for p in placement),
            empty_space_threshold_m=EMPTY_SPACE_M, deep_penetration_threshold_m=DEEP_PENETRATION_M,
            duplicate_symmetric_mean_threshold_m=DUPLICATE_SURFACE_M,
            median_min_gap_to_any_bone_m=float(np.median([p['min_gap_to_any_bone_m'] for p in placement]))
            if placement else None,
            median_expected_bone_min_gap_m=float(np.median(
                [e['min_gap_m'] for p in placement for e in p['expected_bones']]))
            if any(p['expected_bones'] for p in placement) else None,
            deep_inside_basis='canonical rigid_bone role minus intervertebral disks and teeth; the '
                              'rigid_bone role itself carries 23 intervertebral disks and 28 teeth',
            expected_bone_basis='structure wording plus Z-Anatomy authoring collection; over-inclusive '
                                'by construction, so the binding flag is reaching no expected bone at all'),
        promotable=sum(r['promotable'] for r in merged),
        per_class=class_table(records, trials, placement),
        hole_filling_required=[r['structure_id'] for r in records if r['hole_filling_required']],
        tetgen_failures=[dict(structure_id=t['structure_id'], name=t['name'], joint_class=t['joint_class'],
                              faces=t['faces_repaired'], boundary_edges=t['boundary_edges_after'],
                              nonmanifold_edges=t['nonmanifold_edges_after'],
                              status=t['repaired'].get('status'), tets=t['repaired'].get('tets'),
                              tet_volume_vs_surface_relative_error=t['tet_volume_vs_surface_relative_error'],
                              exception=t['repaired'].get('exception'), message=t['repaired'].get('message'),
                              tetgen_stdout=t['repaired']['tetgen_stdout'][-10:])
                         for t in sorted(trials, key=lambda t: t['structure_id'])
                         if not t['repaired']['succeeded']],
        elapsed_s=elapsed)
    write(out / 'summary.json', summary)

    sources = {r['source_path']: r['source_sha256'] for r in records}
    artifacts = {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
    write(out / 'manifest.json', dict(
        schema=SCHEMA, packages={n: importlib.metadata.version(n) for n in ('libigl', 'numpy', 'scipy')},
        python=sys.version,
        native_cgal_sha256=sha(Path(cgal.pyigl_copyleft_cgal.__file__)),
        native_tetgen_sha256=sha(Path(tetgen.pyigl_copyleft_tetgen.__file__)),
        anatomy_sha256=sha(canonical_path), source_index_sha256=sha(index_path),
        bone_reference_build=str(build.relative_to(ROOT)),
        bone_reference_manifest_sha256=sha(build / 'manifest.json'),
        muscle_builder_sha256=sha(ROOT / 'scripts/build_muscle_tet_ready_surfaces.py'),
        muscle_verifier_sha256=sha(ROOT / 'scripts/verify_muscle_tet_ready_surfaces.py'),
        joint_audit_sha256=sha(ROOT / 'scripts/audit_joint_substrate.py'),
        source_geometry_sha256=sources, artifacts_sha256=artifacts,
        registration='data/derived/canonical/anatomy.json registrations.z_anatomy, reproduced before use',
        method='rotate by the recorded source_rotation; apply the recorded affine plus regularized '
               'thin-plate-spline residual; exact weld; drop repeated-index and exactly-collinear faces; '
               'collapse coincident faces (opposite-winding pairs removed as zero-volume fins); bfs_orient '
               'patches flipped to positive signed volume; CGAL remesh_self_intersections(stitch_all); CGAL '
               'self-union outer shell only where the arrangement leaves nonmanifold edges; then TetGen '
               '%s in a forked process with native stdout captured, gated on accounted volume' % TETGEN_FLAGS,
        canonical_geometry_modified=False, promotion_performed=False,
        limitations=[
            'This is a candidate for review. No canonical file is written and no structure is promoted.',
            'No hole is filled. Structures still carrying boundary edges after repair are open sheets or '
            'open sleeves; an open sleeve is not a solid and no volumetric claim is made for it.',
            'Expected attachment bones are inferred from structure wording, not from a curated attachment '
            'table; structures with no name-derived expectation are judged only by nearest-bone distance.',
            'The duplicate screen is geometric and name-blind; it finds overlap with entities canonical '
            'already holds but does not decide which copy should survive.',
            'Reproducing the recorded registration proves correct application, not independent validation. '
            'The 6.18 mm never-fitted centroid RMS is inherited from the prior audit, not recomputed here.',
            'A TetGen success proves the surface is a valid closed PLC whose tets account for its volume. '
            'It establishes nothing about anatomical correctness, inter-structure disjointness, material '
            'ownership, or fitness for a specific materialization.']))
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ('tetgen_failures', 'hole_filling_required', 'per_class')}, indent=2))
    return summary


# ---------------------------------------------------------------- self test

def self_test():
    canonical = json.loads((ROOT / CANONICAL).read_text())
    fit, source, target = fit_registration(canonical['registrations']['z_anatomy'])
    report = fit.report()
    for key in ('affine_rms_m', 'fit_rms_m', 'fit_max_m', 'affine_determinant'):
        assert abs(report[key] - canonical['registrations']['z_anatomy'][key]) <= 1e-12, key

    v = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [0, 0, 1.]])
    f = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]], np.int64)
    good = gate(tetrahedralize(v, f, TETGEN_FLAGS), signed_volume(v, f))
    assert good['succeeded'] and good['succeeded_exit_and_tets_only']
    assert abs(good['tet_volume_vs_surface_relative_error']) < 1e-9
    starved = gate(dict(good, tet_volume_m3=1e-12), signed_volume(v, f))
    assert starved['succeeded_exit_and_tets_only'] and not starved['succeeded'], \
        'a zero-exit run that meshed no volume is not a success'
    sleeve = gate(dict(status=0, tets=5, tet_volume_m3=1e-6), 0.)
    assert not sleeve['succeeded'], 'an open sleeve encloses no volume and cannot pass the gate'
    assert not gate(tetrahedralize(v, f[:3], TETGEN_FLAGS), signed_volume(v, f[:3]))['succeeded']

    bones = {'left femur', 'right femur', 'left tibia', 'left patella', 'atlas', 'axis',
             'fourth lumbar vertebra', 'fifth lumbar vertebra', 'mandible', 'left temporal bone',
             'left fibula', 'left talus'}
    assert expected_bones('Anterior cruciate ligament.l', ['3: Joints', 'Knee joint'], bones)[0] == [
        'left femur', 'left tibia', 'left patella']
    assert expected_bones('Anterior talofibular ligament.l', ['3: Joints'], bones)[0] == [
        'left fibula', 'left talus']
    assert set(expected_bones('Intervertebral disc L4-L5', ['Intervertebral disc'], bones)[0]) == {
        'fourth lumbar vertebra', 'fifth lumbar vertebra'}
    assert expected_bones('Articular disc of temporomandibular joint.l',
                          ['Temporomandibular joint'], bones)[0] == ['mandible', 'left temporal bone']
    assert expected_bones('Some unnamed sheet', ['Fascia'], bones)[0] == []
    assert expected_bones('Subtendinous bursa of tibialis anterior.l', ['Bursae'], bones)[0] == [], \
        'a bursa is named for its muscle, not for a bone'
    assert expected_bones('Suprapatellar bursa.l', ['Knee joint'], bones)[0] == [
        'left femur', 'left tibia', 'left patella'], 'the authoring collection still speaks'

    log = []
    welded = normalize(v[f].reshape(-1, 3), np.arange(12, dtype=np.int64).reshape(-1, 3), log, 'probe')
    assert len(welded[0]) == 4 and log[0]['welded_vertices'] == 4
    assert diagnose(*welded)['closed']
    inverted = repair(v, f[:, [0, 2, 1]])
    assert inverted[2][0]['patches_flipped_outward'] == 1 and signed_volume(*inverted[:2]) > 0
    print('PASS registration reproduction, volume-gated TetGen predicate, attachment lexicon, '
          'imported repair order')


# ---------------------------------------------------------------- driver

def run(output, stage, workers, limit):
    out = Path(output).resolve()
    if not out.is_relative_to(ROOT / 'data/derived'):
        raise ValueError('Output must live under data/derived')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'stages').mkdir(exist_ok=True)
    (out / 'inputs').mkdir(exist_ok=True)
    canonical_path, index_path = ROOT / CANONICAL, ROOT / SOURCE_INDEX
    build = ROOT / BONE_BUILD
    canonical = json.loads(canonical_path.read_text())
    index = json.loads(index_path.read_text())
    if importlib.metadata.version('libigl') != '2.6.2':
        raise ValueError('This build pins libigl 2.6.2')
    for name in ('build_muscle_tet_ready_surfaces.py', 'verify_muscle_tet_ready_surfaces.py',
                 'audit_joint_substrate.py', Path(__file__).name):
        shutil.copyfile(ROOT / 'scripts' / name, out / 'inputs' / name)
    native = Path(cgal.pyigl_copyleft_cgal.__file__)
    shutil.copyfile(native, out / 'inputs' / native.name)

    started = time.monotonic()
    fit = (stage_registration(out, canonical)
           if stage in ('all', 'registration') or not (out / 'stages/registration.json').exists()
           else fit_registration(canonical['registrations']['z_anatomy'])[0])
    if stage == 'registration':
        print(json.dumps(json.loads((out / 'stages/registration.json').read_text())['reproduced_exactly'], indent=2))
        return

    records_path = out / 'stages/build-records.jsonl'
    if stage in ('all', 'build') or not records_path.exists():
        records = stage_build(out, canonical, index, fit, workers, limit)
    else:
        records = sorted((json.loads(l) for l in records_path.read_text().splitlines()),
                         key=lambda r: r['structure_id'])
    if stage == 'build':
        return

    trials_path = out / 'stages/tetgen-trials.jsonl'
    if stage in ('all', 'tetgen') or not trials_path.exists():
        trials = stage_tetgen(out, records, fit, canonical, workers)
    else:
        trials = sorted((json.loads(l) for l in trials_path.read_text().splitlines()),
                        key=lambda t: t['structure_id'])
    if stage == 'tetgen':
        return

    placement_path = out / 'stages/placement.jsonl'
    if stage in ('all', 'placement') or not placement_path.exists():
        placement = stage_placement(out, records, canonical, build)
    else:
        placement = sorted((json.loads(l) for l in placement_path.read_text().splitlines()),
                           key=lambda p: p['structure_id'])
    if stage == 'placement':
        return

    stage_finalize(out, records, trials, placement, canonical_path, index_path, build,
                   time.monotonic() - started)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--output', type=Path)
    p.add_argument('--stage', choices=('all', 'registration', 'build', 'tetgen', 'placement', 'finalize'),
                   default='all')
    p.add_argument('--workers', type=int, default=6)
    p.add_argument('--limit', type=int, default=0)
    a = p.parse_args()
    if a.self_test:
        self_test()
    if a.output:
        run(a.output, a.stage, a.workers, a.limit)
    if not a.self_test and not a.output:
        p.error('Select --self-test or --output')
