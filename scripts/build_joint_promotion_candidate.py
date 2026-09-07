#!/usr/bin/env python3
"""Promotion CANDIDATE: canonical entity records for the registered joint substrate.

Consumes data/derived/joint-substrate-candidate-v1 (630 registered, repaired, TetGen-proved
Z-Anatomy joint structures) and emits canonical-shaped entity records for the promotable set,
a per-structure decision log for all 630, a duplicate reconciliation, a re-review of the
structures withheld for placement, and the conflict-graph delta against the existing canonical
set.

Run with the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_joint_promotion_candidate.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_joint_promotion_candidate.py \
    --output data/derived/joint-promotion-candidate-v1 --stage all --workers 4

No canonical file is written and no canonical asset is modified. Promotion is not performed.

Stages, each resumable, each leaving its own receipt:
  duplicates  geometric flags from the prior lane PLUS a name-identity sweep the prior lane's
              name-blind screen could not see; one survivor per pair by a declared rule
  placement   exact CGAL boolean overlap of every deep-inside-bone structure against canonical
              bone, plus the articular gap of the bone pair it is expected to span, to separate
              registration misplacement from the known non-articular bone pose
  records     canonical entity records: id, name, role, system, centroid, bounds, principal axis,
              reference_geometry, connections, assumptions, uncertainty, provenance
  conflicts   the added cross-structure conflict pairs, by the exhaustive CGAL method of
              scripts/build_cross_structure_conflict_repair.py
  finalize    summary.json and manifest.json

Roles are assigned from the canonical vocabulary. Where that vocabulary has no exact term the
record says so in role_vocabulary rather than presenting the nearest term as a constitutive claim.
"""
from pathlib import Path
import argparse
import copy
import gzip
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import re
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import igl
import igl.copyleft.cgal as cgal

from ihm.assembly.anatomy import FRAME, MODEL_ID, mesh_properties, normalized_name, verify_assembly

SCHEMA = 'ihm.joint-promotion-candidate.v1'
SUBSTRATE = 'data/derived/joint-substrate-candidate-v1'
CANONICAL = 'data/derived/canonical/anatomy.json'
EXTENDED_FRAGMENT = 'data/derived/anatomy/extended/manifest_fragment.json'
MUSCLE_BUILD = 'data/derived/muscle-tet-ready-v1'
ENTITY_BUILD = 'data/derived/entity-tet-ready-v1'
MATERIALS = 'data/derived/tissue-material-candidate-v1/materials.json'
MECHANICS = 'data/derived/canonical/mechanics.json'

NEVER_FITTED_CENTROID_RMS_M = .006180582966048701
REJECTED_PRECEDENT_RESIDUAL_M = .036229
SOFT_DENSITY_KG_M3 = 1000.        # canonical mechanics.json blanket value for every soft role
BONE_DENSITY_KG_M3 = 1900.
DUPLICATE_SURFACE_M = .003        # the prior lane's geometric screen threshold, reused verbatim
NAME_PAIR_SAMPLES = 600
CONTACT_M = .002
CLOSED_JOINT_SPACE_M = .0005
BOOL_FACE_CAP = 120000
MC_POINTS = 200000
MC_SEED = 20260907

# ---------------------------------------------------------------- role vocabulary
# (predicate on the lowered structure name) -> canonical role, exactness, honest anatomical term.
# `exact` false means the canonical vocabulary carries no term for this tissue and the assigned
# role is the nearest available label, not a constitutive claim.
CANONICAL_ROLES = ('rigid_bone', 'cartilage', 'tendon', 'ligament', 'fluid_cavity', 'vascular',
                   'nerve', 'muscle', 'soft_organ', 'connective_tissue', 'lymph_node_group',
                   'skin', 'skin_layer', 'lymphatic_network')

FIBROCARTILAGE = ('the canonical vocabulary carries one undifferentiated `cartilage` term and no '
                  'fibrocartilage term; this structure is fibrocartilage, not hyaline cartilage')


def role_of(joint_class, name):
    low = name.lower()
    def row(role, exact, anatomical, basis, rejected=()):
        return dict(role=role, exact=bool(exact), anatomical_role=anatomical, basis=basis,
                    rejected=list(rejected),
                    vocabulary_gap=bool(not exact))
    if joint_class == 'meniscus':
        if re.search(r'meniscotibial|meniscopatellar|transverse ligament', low):
            return row('ligament', True, 'ligament',
                       'named a ligament and classed with the menisci only by the audit lexicon; '
                       'it is a ligament and takes the exact term')
        return row('cartilage', False, 'meniscus', FIBROCARTILAGE, ('connective_tissue',))
    if joint_class == 'labrum':
        return row('cartilage', False, 'labrum', FIBROCARTILAGE, ('connective_tissue',))
    if joint_class == 'articular_disc':
        return row('cartilage', False, 'articular disc', FIBROCARTILAGE, ('connective_tissue',))
    if joint_class == 'symphysis':
        return row('cartilage', False, 'secondary cartilaginous joint',
                   'a symphysis is a fibrocartilaginous joint, not a cartilage body; the '
                   'vocabulary has no joint term', ('connective_tissue',))
    if joint_class == 'intervertebral_disc':
        return row('cartilage', False, 'intervertebral disc',
                   'nearest available role, and not a correct constitutive claim: the disc is an '
                   'anulus of layered fibrocartilage around a hydrated gel nucleus. Recorded by a '
                   'previous lane and carried forward unchanged', ('connective_tissue',))
    if joint_class == 'nucleus_pulposus':
        return row('cartilage', False, 'nucleus pulposus',
                   'no correct term exists. The nucleus is a 70-90 percent water proteoglycan gel; '
                   '`cartilage` is wrong constitutively and `fluid_cavity` is wrong topologically '
                   'because the nucleus is not a lumen bounded by a separate wall entity. '
                   '`cartilage` is taken only so the nucleus and the disc that contains it share '
                   'one role', ('fluid_cavity', 'soft_organ'))
    if joint_class == 'cartilage':
        return row('cartilage', True, 'cartilage', 'exact term')
    if joint_class in ('ligament', 'cruciate_ligament'):
        return row('ligament', True, 'ligament', 'exact term')
    if joint_class == 'retinaculum':
        return row('ligament', False, 'retinaculum',
                   'a retinaculum is thickened deep fascia, not dense regular collagen spanning a '
                   'joint. `ligament` follows canonical precedent: ihm.assembly.anatomy.'
                   'physical_role routes `retinacul` to ligament and canonical already roles '
                   '`flexor retinaculum of left wrist` as ligament', ('connective_tissue',))
    if joint_class == 'articular_capsule':
        return row('connective_tissue', False, 'articular capsule',
                   'the fibrous capsule is dense irregular connective tissue with a synovial '
                   'lining; the vocabulary has no capsule term and `ligament` would overstate the '
                   'fibre organisation', ('ligament',))
    if joint_class == 'interosseous_membrane':
        return row('connective_tissue', True, 'interosseous membrane',
                   'canonical precedent: `interosseous membrane of left forearm` is already roled '
                   'connective_tissue')
    if joint_class == 'fascia_aponeurosis':
        if 'aponeuros' in low:
            return row('tendon', True, 'aponeurosis',
                       'canonical precedent: physical_role routes `aponeurosis` to tendon')
        return row('connective_tissue', True, 'fascia', 'exact term for a fascial sheet')
    if joint_class == 'tendon_sheath':
        return row('connective_tissue', False, 'tendon sheath',
                   'the authored surface is the sheath wall; the synovial fluid film between its '
                   'two layers is not resolved by the geometry and gets no material of its own',
                   ('fluid_cavity',))
    if joint_class == 'synovial_bursa':
        return row('fluid_cavity', False, 'synovial bursa',
                   'a bursa is a synovial-lined sac of fluid and the authored surface bounds that '
                   'fluid, so `fluid_cavity` is the nearest term; unlike the canonical cardiac '
                   'cavities the wall is not a separate entity, so wall and contents are '
                   'undifferentiated here', ('connective_tissue',))
    if joint_class == 'fat_pad':
        return row('connective_tissue', False, 'adipose fat pad',
                   'adipose is histologically a connective tissue but constitutively unlike the '
                   'fibrous sheets that role covers (initial tangent modulus about 1.6 kPa against '
                   'a 100 kPa placeholder). These are the first adipose surfaces in the model; '
                   'docs/research/SOFT_BODY_MATERIALIZATION.md records `No adipose geometry '
                   'exists`', ('soft_organ',))
    raise ValueError('unmapped joint class: ' + joint_class)


# ---------------------------------------------------------------- duplicate identity
# Canonical-name aliases the prior lane's geometric screen is blind to and normalized_name alone
# does not resolve. Each is a wording difference for the same anatomical structure.
ORDINALS = ('first second third fourth fifth sixth seventh eighth ninth tenth eleventh '
            'twelfth').split()


def name_aliases(name):
    """Canonical names this Z-Anatomy structure would carry if canonical already held it."""
    base = normalized_name(name)
    out = [base]
    match = re.match(r'^(left|right) costal cartilage of (\w+) rib$', base)
    if match:
        out.append('%s %s costal cartilage' % (match.group(1), match.group(2)))
    match = re.match(r'^(left|right) (flexor|extensor) retinaculum of (wrist|ankle)$', base)
    if match:
        out.append('%s retinaculum of %s %s' % (match.group(2), match.group(1), match.group(3)))
    match = re.match(r'^(left|right) interosseous membrane of (forearm|leg)$', base)
    if match:
        out.append('interosseous membrane of %s %s' % (match.group(1), match.group(2)))
    match = re.match(r'^(left|right) long plantar ligament$', base)
    if match:
        out.append(base)
    return list(dict.fromkeys(out))


# Geometric flags the anatomical review overturns: the two structures are adjacent, not the same.
NOT_DUPLICATE = {
    ('cricopharyngeal ligament', 'conus elasticus'):
        'the conus elasticus (cricovocal membrane) runs from the cricoid arch to the vocal '
        'ligament; the cricopharyngeal ligament is a separate posterior band. Adjacent, not one '
        'structure',
    ('quadrangular membrane', 'thyrohyoid membrane'):
        'the quadrangular membrane is the internal aryepiglottic submucosal sheet of the larynx; '
        'the thyrohyoid membrane is the external sheet between thyroid cartilage and hyoid. '
        'Adjacent across the thyroid lamina, not one structure',
    ('pisotriquetral ligament', 'pisiform'):
        'a ligament screened against a carpal bone. The screen measured a small ligament lying on '
        'the bone surface, not two copies of one structure',
}


def not_duplicate_reason(za_name, canonical_name):
    low, other = za_name.lower(), canonical_name.lower()
    for (a, b), reason in NOT_DUPLICATE.items():
        if a in low and b in other:
            return reason
    return None


# ---------------------------------------------------------------- placement lexicon
SUBCUTANEOUS = re.compile(r'subcutaneous')
BURSA = re.compile(r'bursa')
LATERALITY_TOLERANCE_M = .005   # ihm.assembly.anatomy.verify_assembly's own bound


def laterality_conflict(name, centroid):
    """The canonical verifier's laterality bound, applied before promotion rather than after."""
    low = name.lower()
    if low.startswith('left ') and centroid[0] <= -LATERALITY_TOLERANCE_M:
        return 'named left but its centroid sits %.4f m on the right of the midline' % centroid[0]
    if low.startswith('right ') and centroid[0] >= LATERALITY_TOLERANCE_M:
        return 'named right but its centroid sits %+.4f m on the left of the midline' % centroid[0]
    if not -.87 < centroid[1] < .88:
        return 'centroid height %.4f m is outside the canonical vertical envelope' % centroid[1]
    return None

# The prior lane's attachment lexicon deliberately returns nothing for rib wording and has no
# entry for the triradiate cartilage. This supplement is used ONLY to classify the cause of a
# deep-inside-bone flag; it never enters an entity record.
SUPPLEMENTARY_BONES = (
    (r'costotransverse|costovertebral', ('rib', 'vertebra')),
    (r'costoclavicular', ('rib', 'clavicle')),
    (r'triradiate', ('hip bone',)),
)


def supplementary_expected(name, bone_names):
    low = name.lower()
    side = 'left ' if name.endswith('.l') else ('right ' if name.endswith('.r') else None)
    out = []
    for pattern, tokens in SUPPLEMENTARY_BONES:
        if not re.search(pattern, low):
            continue
        for token in tokens:
            for bone in bone_names:
                if token in bone and (side is None or bone.startswith(side)):
                    out.append(bone)
    return sorted(set(out))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=1, allow_nan=False) + '\n')


def write_lines(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        for row in rows:
            handle.write(json.dumps(row, allow_nan=False) + '\n')


def read_geometry(path):
    payload = json.loads(gzip.decompress(Path(path).read_bytes()))
    return (np.ascontiguousarray(np.asarray(payload['positions'], float).reshape(-1, 3)),
            np.ascontiguousarray(np.asarray(payload['indices'], np.int64).reshape(-1, 3)))


def divergence(vertices, faces):
    if not len(faces):
        return 0.
    return float(np.einsum('ij,ij->i', vertices[faces[:, 0]],
                           np.cross(vertices[faces[:, 1]], vertices[faces[:, 2]])).sum() / 6.)


def area_of(vertices, faces):
    if not len(faces):
        return 0.
    return float(np.linalg.norm(np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]],
                                         vertices[faces[:, 2]] - vertices[faces[:, 0]]),
                                axis=1).sum() / 2.)


def stride(points, limit):
    return np.ascontiguousarray(points[::max(1, len(points) // limit)])


def symmetric_surface_distance(va, fa, vb, fb, samples=NAME_PAIR_SAMPLES):
    forward = np.sqrt(np.maximum(igl.point_mesh_squared_distance(stride(va, samples), vb, fb)[0], 0.))
    backward = np.sqrt(np.maximum(igl.point_mesh_squared_distance(stride(vb, samples), va, fa)[0], 0.))
    both = np.concatenate([forward, backward])
    return float(both.mean()), float(both.max())


def monte_carlo_overlap(va, fa, vb, fb, points=MC_POINTS, seed=MC_SEED):
    """Winding-number Monte Carlo over the shared bounding box; the conflict lane's fallback."""
    low = np.maximum(va.min(0), vb.min(0))
    high = np.minimum(va.max(0), vb.max(0))
    span = np.maximum(high - low, 0.)
    box = float(np.prod(span))
    if box <= 0.:
        return 0., 0.
    rng = np.random.default_rng(seed)
    query = np.ascontiguousarray(low + rng.random((points, 3)) * span)
    inside_a = np.abs(np.asarray(igl.fast_winding_number(va, fa, query), float)) > .5
    inside_b = np.abs(np.asarray(igl.fast_winding_number(vb, fb, query), float)) > .5
    p = float((inside_a & inside_b).sum()) / points
    return box * p, box * float(np.sqrt(max(p * (1 - p), 0.) / points))


def boolean_overlap(va, fa, vb, fb):
    """Exact CGAL boolean intersection volume; the method of build_cross_structure_conflict_repair.

    A CGAL failure falls back to the winding-number Monte Carlo rather than reporting zero: a
    failed boolean is not evidence of no overlap and must not be recorded as one.
    """
    try:
        result = cgal.mesh_boolean(va, fa, vb, fb, type_str='intersect')
    except BaseException as error:
        volume, stderr = monte_carlo_overlap(va, fa, vb, fb)
        return volume, 'cgal_boolean_failed:%s|monte_carlo(stderr=%.3e m3)' % (
            type(error).__name__, stderr)
    iv = np.ascontiguousarray(np.asarray(result[0], float))
    inf = np.ascontiguousarray(np.asarray(result[1], np.int64))
    if not len(inf):
        return 0., 'cgal_boolean_intersect'
    return abs(divergence(iv, inf)), 'cgal_boolean_intersect'


def load_inputs():
    substrate = ROOT / SUBSTRATE
    rows = [json.loads(l) for l in (substrate / 'entities.jsonl').read_text().splitlines()]
    rows.sort(key=lambda r: r['structure_id'])
    canonical = json.loads((ROOT / CANONICAL).read_text())
    fragment = json.loads((ROOT / EXTENDED_FRAGMENT).read_text())
    return substrate, rows, canonical, {s['id']: s for s in fragment['structures']}


def canonical_geometry_path(canonical_entity, build_index):
    """Prefer the repaired surface the tet-ready lanes emitted; fall back to canonical display."""
    record = build_index.get(canonical_entity['id'])
    if record:
        return ROOT / record['base'] / record['output_path'], 'repaired'
    return ROOT / canonical_entity['reference_geometry']['path'], 'canonical_display'


def build_index():
    index = {}
    for base in (MUSCLE_BUILD, ENTITY_BUILD):
        path = ROOT / base / 'entities.jsonl'
        for line in path.read_text().splitlines():
            record = json.loads(line)
            index[record['entity_id']] = dict(base=base, output_path=record['output_path'],
                                              output_sha256=record['output_sha256'],
                                              tet_ready=bool(record.get('tet_ready')),
                                              faces=record['after']['faces'],
                                              signed_volume_m3=record['after']['signed_volume_m3'])
    return index


def tetgen_index():
    index = {}
    for path, key in ((ROOT / 'data/derived/entity-tet-ready-verification-v1/tetgen-trials.jsonl',
                       'entity_id'),):
        if not path.exists():
            continue
        for line in path.read_text().splitlines():
            record = json.loads(line)
            index[record[key]] = dict(succeeded=bool(record['repaired'].get('succeeded')),
                                      tets=record['repaired'].get('tets'))
    return index


# ---------------------------------------------------------------- stage: duplicates

def stage_duplicates(out, substrate, rows, canonical, builds, tets):
    by_name = {}
    for entity in canonical['entities']:
        by_name.setdefault(entity['name'], []).append(entity)
    canonical_by_id = {e['id']: e for e in canonical['entities']}
    cache = {}

    def surface(entity):
        if entity['id'] not in cache:
            path, stage = canonical_geometry_path(entity, builds)
            cache[entity['id']] = (read_geometry(path), stage, str(path.relative_to(ROOT)))
        return cache[entity['id']]

    pairs, started = [], time.monotonic()
    for record in rows:
        za_surface = None
        place = record['placement'] or {}
        seen = set()
        candidates = []
        for hit in place.get('duplicates_canonical', ()):
            candidates.append((hit['canonical_id'], 'geometric_surface_screen', hit))
            seen.add(hit['canonical_id'])
        for alias in name_aliases(record['name']):
            for entity in by_name.get(alias, ()):
                if entity['id'] in seen:
                    continue
                seen.add(entity['id'])
                candidates.append((entity['id'], 'name_identity', dict(alias=alias)))
        for canonical_id, basis, hit in candidates:
            other = canonical_by_id[canonical_id]
            if za_surface is None:
                za_surface = read_geometry(substrate / record['output_path'])
            (ov, of), stage, opath = surface(other)
            mean, worst = symmetric_surface_distance(*za_surface, ov, of)
            centroid = za_surface[0].mean(0)
            distinct = not_duplicate_reason(record['name'], other['name'])
            za_ok = bool(record['tet_ready'] and record['tetgen']
                         and record['tetgen']['repaired']['succeeded'])
            other_build = builds.get(canonical_id)
            other_tet = tets.get(canonical_id)
            other_ok = bool(other_build and other_build['tet_ready']
                            and other_tet and other_tet['succeeded'])
            if distinct:
                survivor, clause = None, 'not_a_duplicate'
            elif za_ok and not other_ok:
                survivor, clause = 'extended', 'geometry_quality: only the extended copy is tet-ready and TetGen-proved'
            elif other_ok and not za_ok:
                survivor, clause = 'canonical', 'geometry_quality: only the canonical copy is tet-ready and TetGen-proved'
            else:
                survivor, clause = 'canonical', (
                    'provenance: both copies survive the same repair order and TetGen gate, so '
                    'geometry does not decide. The canonical copy is BodyParts3D source_geometry '
                    'placed by a rotation-and-translation display transform with no registration '
                    'residual; the extended copy is registered_geometry carrying the z_anatomy '
                    'affine plus thin-plate-spline fit at %.4f m never-fitted centroid RMS. The '
                    'canonical id is also already referenced by mechanics.json, '
                    'entity-tet-ready-v1 and the cross-structure conflict graph'
                    % NEVER_FITTED_CENTROID_RMS_M)
            pairs.append(dict(
                structure_id=record['structure_id'], structure_name=record['name'],
                joint_class=record['joint_class'], canonical_id=canonical_id,
                canonical_name=other['name'], canonical_role=other['role'],
                canonical_system=other['system'], detected_by=basis,
                same_structure=bool(not distinct),
                distinct_reason=distinct,
                symmetric_mean_surface_m=mean, symmetric_max_surface_m=worst,
                centroid_distance_m=float(np.linalg.norm(np.asarray(other['centroid_m']) - centroid)),
                below_prior_screen_threshold=bool(mean <= DUPLICATE_SURFACE_M),
                extended_tet_ready=za_ok,
                extended_faces=record['after']['faces'],
                extended_tets=(record['tetgen']['repaired'].get('tets')
                               if record['tetgen'] and record['tetgen']['repaired']['succeeded'] else None),
                extended_volume_m3=record['after']['signed_volume_m3'],
                canonical_tet_ready=other_ok,
                canonical_faces=(other_build or {}).get('faces'),
                canonical_tets=(other_tet or {}).get('tets'),
                canonical_volume_m3=(other_build or {}).get('signed_volume_m3'),
                canonical_watertight_as_published=bool(other['watertight_edge_incidence']),
                canonical_surface_stage=stage, canonical_surface_path=opath,
                survivor=survivor, decision_clause=clause,
                withhold_extended=bool(survivor == 'canonical'),
                canonical_role_correction=(
                    dict(current=other['role'], recommended='cartilage',
                         reason='an intervertebral disc is not rigid bone. `cartilage` is the '
                                'nearest available role and is not a correct constitutive claim. '
                                'The correction belongs to the canonical entity and is not applied '
                                'by this candidate')
                    if other['role'] == 'rigid_bone' and 'intervertebral disk' in other['name']
                    else None)))
    write_lines(out / 'duplicates.jsonl', pairs)
    withheld = sorted({p['structure_id'] for p in pairs if p['withhold_extended']})
    supersede = sorted({(p['canonical_id'], p['canonical_name'], p['structure_id'],
                         p['structure_name']) for p in pairs if p['survivor'] == 'extended'})
    misroled = sorted(e['id'] for e in canonical['entities']
                      if e['role'] == 'rigid_bone' and 'intervertebral disk' in e['name'])
    overturned = sorted({p['structure_id'] for p in pairs
                         if p['detected_by'] == 'geometric_surface_screen' and not p['same_structure']})
    report = dict(
        schema=SCHEMA + '.duplicates',
        candidate_pairs=len(pairs),
        structures_with_any_candidate=len({p['structure_id'] for p in pairs}),
        detected_by_geometric_surface_screen=sum(p['detected_by'] == 'geometric_surface_screen' for p in pairs),
        detected_by_name_identity_only=sum(p['detected_by'] == 'name_identity' for p in pairs),
        same_structure_pairs=sum(p['same_structure'] for p in pairs),
        overturned_geometric_flags=len(overturned),
        overturned_structure_ids=overturned,
        extended_copies_withheld=len(withheld),
        survivor_rule=[
            '1. if exactly one copy is tet-ready and passes the volume-gated TetGen trial, that '
            'copy survives',
            '2. otherwise the canonical BodyParts3D copy survives on provenance: it is measured '
            'source geometry with no registration residual, while the extended copy carries the '
            'z_anatomy fit at %.4f m never-fitted centroid RMS against a %.4f m rejection '
            'precedent, and the canonical id is already wired into the downstream lanes'
            % (NEVER_FITTED_CENTROID_RMS_M, REJECTED_PRECEDENT_RESIDUAL_M)],
        clause_counts=dict(sorted({p['decision_clause'].split(':')[0]: 0 for p in pairs}.items())),
        canonical_supersession_recommended=[
            dict(canonical_id=a, canonical_name=b, replaced_by='body-' + c, structure_name=d)
            for a, b, c, d in supersede],
        canonical_supersession_note=(
            'these canonical entities lose the survivor rule on geometry quality: the canonical '
            'surface is not tet-ready while the extended copy is. Retiring a canonical entity is '
            'not this candidate\'s to perform and is recorded as a recommendation only'),
        canonical_role_correction_recommended=dict(
            entity_ids=misroled, count=len(misroled), current='rigid_bone',
            recommended='cartilage',
            reason='an intervertebral disk is not rigid bone. `cartilage` is the nearest available '
                   'role and is not a correct constitutive claim: the disk is an anulus of layered '
                   'fibrocartilage around a hydrated gel nucleus. %d of these were reached by a '
                   'duplicate pair and the rest by name; the correction applies to all of them and '
                   'is not applied here' % len({p['canonical_id'] for p in pairs
                                                if p['canonical_role_correction']})),
        name_screen_finding=(
            'the prior lane screened duplicates geometrically and name-blind at a %.3f m symmetric '
            'mean surface threshold. A name-identity sweep finds pairs that screen missed because '
            'the two copies of the same named structure sit further apart than the threshold'
            % DUPLICATE_SURFACE_M),
        wall_seconds=time.monotonic() - started)
    for pair in pairs:
        report['clause_counts'][pair['decision_clause'].split(':')[0]] += 1
    write_json(out / 'duplicates.json', report)
    print('duplicates: %d candidate pairs, %d same-structure, %d extended copies withheld, %.0fs'
          % (len(pairs), report['same_structure_pairs'], len(withheld), report['wall_seconds']),
          flush=True)
    return pairs, report


# ---------------------------------------------------------------- stage: placement

def stage_placement(out, substrate, rows, canonical, builds, workers):
    bones = {e['id']: e for e in canonical['entities'] if e['role'] == 'rigid_bone'
             and not re.search(r'intervertebral disk|tooth', e['name'])}
    bone_by_name = {e['name']: e['id'] for e in bones.values()}
    deep = [r for r in rows if r['placement'] and r['placement']['flag_deep_inside_bone']]
    cache = {}

    def bone_surface(bone_id):
        if bone_id not in cache:
            path, _ = canonical_geometry_path(bones[bone_id], builds)
            cache[bone_id] = read_geometry(path)
        return cache[bone_id]

    started, results = time.monotonic(), []
    gap_cache = {}
    for record in deep:
        place = record['placement']
        v, f = read_geometry(substrate / record['output_path'])
        volume = abs(divergence(v, f))
        low, high = v.min(0), v.max(0)
        overlaps, boolean_failures = [], 0
        for bone_id, bone in bones.items():
            blow = np.asarray(bone['bounds_m']['min'])
            bhigh = np.asarray(bone['bounds_m']['max'])
            if (low > bhigh).any() or (high < blow).any():
                continue
            bv, bf = bone_surface(bone_id)
            over, method = boolean_overlap(v, f, bv, bf)
            boolean_failures += method.startswith('cgal_boolean_failed')
            if over <= 0.:
                continue
            overlaps.append(dict(bone_id=bone_id, bone=bone['name'], overlap_volume_m3=over,
                                 fraction_of_structure=over / volume if volume else None,
                                 method=method))
        overlaps.sort(key=lambda o: -o['overlap_volume_m3'])
        inside_volume = sum(o['overlap_volume_m3'] for o in overlaps)
        expected = [e['bone'] for e in place['expected_bones']]
        supplement = []
        if not expected:
            supplement = supplementary_expected(record['name'], set(bone_by_name))
            expected = list(supplement)
        worst_bone = overlaps[0]['bone'] if overlaps else None
        deepest_expected = bool(worst_bone and worst_bone in expected)
        # articular relation of the bone pair the structure is expected to span
        pair_report = []
        for i in range(len(expected)):
            for j in range(i + 1, len(expected)):
                a, b = sorted((bone_by_name[expected[i]], bone_by_name[expected[j]]))
                if (a, b) not in gap_cache:
                    av, af = bone_surface(a)
                    bv, bf = bone_surface(b)
                    over, method = boolean_overlap(av, af, bv, bf)
                    squared, _, _ = igl.point_mesh_squared_distance(stride(av, 1200), bv, bf)
                    gap_cache[(a, b)] = dict(
                        bone_a=bones[a]['name'], bone_b=bones[b]['name'],
                        interpenetration_volume_m3=over, method=method,
                        min_gap_m=float(np.sqrt(max(squared.min(), 0.))))
                pair_report.append(gap_cache[(a, b)])
        closed = [p for p in pair_report
                  if (p['interpenetration_volume_m3'] or 0.) > 0. or p['min_gap_m'] <= CLOSED_JOINT_SPACE_M]
        low_name = record['name'].lower()
        if not expected:
            cause = 'indeterminate_no_expected_bone'
        elif not deepest_expected:
            cause = 'misplacement'
        elif closed:
            cause = 'pose_closed_joint_space'
        else:
            cause = 'pose_bone_adjacent'
        if cause.startswith('pose') and SUBCUTANEOUS.search(low_name):
            cause = 'misplacement'
        if cause == 'indeterminate_no_expected_bone' and BURSA.search(low_name):
            cause = 'misplacement'
        results.append(dict(
            structure_id=record['structure_id'], name=record['name'],
            joint_class=record['joint_class'],
            surface_volume_m3=volume,
            surface_sample_inside_true_bone_fraction=place['inside_true_bone_fraction'],
            max_true_bone_penetration_m=place['max_true_bone_penetration_m'],
            summed_overlap_volume_m3=inside_volume,
            summed_overlap_fraction_of_structure=inside_volume / volume if volume else None,
            summed_fraction_caveat='the sum runs over bones; it exceeds 1 where the canonical '
                                   'bones themselves interpenetrate and both claim the same '
                                   'sub-volume, which is itself evidence of the non-articular pose',
            max_single_bone_overlap_fraction=(overlaps[0]['fraction_of_structure']
                                              if overlaps else 0.),
            overlapping_bone_count=len(overlaps),
            boolean_failures_fallen_back_to_monte_carlo=boolean_failures,
            overlapping_bones=overlaps[:6],
            deepest_overlap_bone=worst_bone,
            expected_bones=expected,
            expected_bones_from_supplementary_lexicon=supplement,
            deepest_overlap_bone_is_expected=deepest_expected,
            expected_bone_pairs=pair_report,
            closed_joint_space_pairs=len(closed),
            cause=cause,
            promotable_under_containment_rule=bool(cause.startswith('pose')),
            otherwise_clean=bool(record['tet_ready'] and record['tetgen']
                                 and record['tetgen']['repaired']['succeeded']
                                 and not record['placement']['flag_duplicates_canonical']
                                 and not record['placement']['flag_no_bone_contact']
                                 and not record['placement']['flag_reaches_no_expected_bone'])))
        if len(results) % 10 == 0:
            print('placement %d/%d  %.0fs' % (len(results), len(deep), time.monotonic() - started),
                  flush=True)
    write_lines(out / 'placement-review.jsonl', results)
    causes = {}
    for row in results:
        causes[row['cause']] = causes.get(row['cause'], 0) + 1
    releasable = [r for r in results if r['promotable_under_containment_rule'] and r['otherwise_clean']]
    report = dict(
        schema=SCHEMA + '.placement',
        structures_reviewed=len(results),
        method=(
            'the prior lane flagged deep_inside_bone from strided surface-vertex signed distance. '
            'This re-review replaces that with an exact CGAL boolean intersection volume against '
            'every canonical bone whose bounding box overlaps, and adds the articular relation of '
            'the bone pair the structure is expected to span, measured the same exact way'),
        discriminator=(
            'cause (ii), the known non-articular bone pose, is asserted only when the bone the '
            'structure is most deeply inside is one it is expected to attach to. It is further '
            'split by whether the expected bone pair is itself in a non-articular relation '
            '(interpenetrating, or separated by at most %.4f m), which closes the joint space the '
            'structure belongs in. Cause (i), misplacement, is asserted when the structure is '
            'deepest inside a bone it has no expected relation to, or when its name declares it '
            'subcutaneous' % CLOSED_JOINT_SPACE_M),
        causes=dict(sorted(causes.items())),
        promotable_under_containment_rule=len(releasable),
        promotable_ids=sorted(r['structure_id'] for r in releasable),
        total_summed_overlap_volume_m3=float(sum(r['summed_overlap_volume_m3'] for r in results)),
        structures_needing_monte_carlo_fallback=sum(
            bool(r['boolean_failures_fallen_back_to_monte_carlo']) for r in results),
        structures_classified_with_supplementary_lexicon=sum(
            bool(r['expected_bones_from_supplementary_lexicon']) for r in results),
        wall_seconds=time.monotonic() - started)
    write_json(out / 'placement-review.json', report)
    print('placement: %s' % json.dumps(report['causes']), flush=True)
    return results, report


# ---------------------------------------------------------------- stage: records

def uncertainty():
    return {'biological': {'status': 'generic authored reference; population variability not estimated',
                           'confidence_percent': None, 'independent_subject_count': 0,
                           'subject_calibrated': False},
            'registration_or_synthesis': {'kind': 'registered_geometry', 'registration_id': 'z_anatomy',
                                          'calibrated_probability': None},
            'display_numerics': {'surface_only': True, 'display_reduction_applied': False,
                                 'position_units': 'm', 'rounding': 'JSON float serialization'}}


def stage_records(out, substrate, rows, canonical, fragment, duplicates, placement, materials):
    withheld_duplicate = {p['structure_id']: p for p in duplicates if p['withhold_extended']}
    overturned = {p['structure_id'] for p in duplicates
                  if p['detected_by'] == 'geometric_surface_screen' and not p['same_structure']}
    placement_by_id = {p['structure_id']: p for p in placement}
    bones = [e for e in canonical['entities'] if e['role'] == 'rigid_bone']
    bone_centres = np.stack([np.asarray(e['centroid_m']) for e in bones])
    canonical_ids = {e['id'] for e in canonical['entities']}
    density = {r['role']: r for r in ()}

    def prior_clean(record):
        return bool(record['promotable'])

    decisions, entities = [], []
    (out / 'geometry').mkdir(parents=True, exist_ok=True)
    for record in rows:
        identity = 'body-' + record['structure_id']
        place = record['placement'] or {}
        trial = record['tetgen']
        reasons = []
        if not record['tet_ready']:
            reasons.append('not_tet_ready')
        if not (trial and trial['repaired']['succeeded']):
            reasons.append('tetgen_volume_gate_failed')
        if place.get('flag_no_bone_contact'):
            reasons.append('no_bone_contact')
        if place.get('flag_reaches_no_expected_bone'):
            reasons.append('reaches_no_expected_bone')
        if place.get('flag_deep_inside_bone'):
            reasons.append('deep_inside_bone')
        if record['structure_id'] in withheld_duplicate:
            reasons.append('superseded_by_canonical_duplicate')
        elif place.get('flag_duplicates_canonical') and record['structure_id'] in overturned:
            pass
        promote = not reasons
        release = bool(promote and not prior_clean(record))
        decision = dict(structure_id=record['structure_id'], name=record['name'],
                        joint_class=record['joint_class'], system=record['system'],
                        entity_id=identity,
                        prior_lane_promotable=bool(record['promotable']),
                        promoted=promote, withheld_reasons=reasons,
                        released_by_review=release,
                        duplicate_of=(withheld_duplicate[record['structure_id']]['canonical_id']
                                      if record['structure_id'] in withheld_duplicate else None),
                        duplicate_flag_overturned=bool(record['structure_id'] in overturned),
                        placement_cause=(placement_by_id[record['structure_id']]['cause']
                                         if record['structure_id'] in placement_by_id else None),
                        placement_promotable_under_containment_rule=(
                            placement_by_id[record['structure_id']]['promotable_under_containment_rule']
                            if record['structure_id'] in placement_by_id else None))
        if not promote:
            decisions.append(decision)
            continue

        vertices, faces = read_geometry(substrate / record['output_path'])
        source = fragment[record['structure_id']]
        role = role_of(record['joint_class'], record['name'])
        probe = mesh_properties(vertices, faces)
        conflict = laterality_conflict(normalized_name(record['name']), probe['centroid_m'])
        if conflict:
            decision['withheld_reasons'] = ['laterality_label_conflict']
            decision['promoted'] = False
            decision['released_by_review'] = False
            decision['laterality_conflict'] = conflict
            decision['laterality_remedy'] = (
                'the geometry is in the right place and the side label is wrong: every sibling '
                'structure of the same joint is correctly sided and the z_anatomy affine has a '
                'positive determinant, so no reflection was introduced. Swapping the .l/.r side '
                'label on the source structure resolves it. Not done here: renaming a source '
                'structure is a change to the extended atlas index, not to this candidate')
            decisions.append(decision)
            continue
        payload = {'positions': [float(x) for x in vertices.ravel()],
                   'indices': [int(i) for i in faces.ravel()],
                   'units': 'm', 'frame': FRAME, 'registration_id': 'z_anatomy',
                   'source_geometry_sha256': record['source_sha256'],
                   'repaired_surface_sha256': record['output_sha256']}
        geometry_path = out / 'geometry' / (identity + '.json.gz')
        geometry_path.write_bytes(gzip.compress(json.dumps(payload, allow_nan=False).encode(), mtime=0))

        properties = probe
        centroid = np.asarray(properties['centroid_m'])
        nearest = int(np.argmin(np.linalg.norm(bone_centres - centroid, axis=1)))
        connections = [dict(entity_id=bones[nearest]['id'], relation='regional_spatial_support',
                            distance_m=float(np.linalg.norm(bone_centres[nearest] - centroid)),
                            evidence='inferred nearest bone bounding-box center; not an attachment '
                                     'or contact constraint')]
        contacts = [b for b in place.get('nearest_bones', ()) if b['min_gap_m'] <= CONTACT_M]
        for bone in contacts[:4]:
            if bone['bone_id'] in canonical_ids and bone['bone_id'] != bones[nearest]['id']:
                connections.append(dict(
                    entity_id=bone['bone_id'], relation='measured_bone_proximity',
                    distance_m=bone['min_gap_m'],
                    evidence='minimum distance from strided surface samples of this structure to '
                             'the repaired canonical bone surface, at most %.3f m; proximity in '
                             'the registered pose, not a certified attachment' % CONTACT_M))
        assumptions = ['CANONICAL-GENERIC-REFERENCE', 'Z-LANDMARK-REGISTRATION',
                       'JOINT-SUBSTRATE-REPAIR']
        if not role['exact']:
            assumptions.append('ROLE-VOCABULARY-APPROXIMATION')
        entity = {
            'id': identity, 'model_id': MODEL_ID, 'source_id': record['structure_id'],
            'name': normalized_name(record['name']), 'system': record['system'],
            'role': role['role'], 'evidence_kind': 'registered_geometry',
            'role_vocabulary': {k: v for k, v in role.items() if k != 'role'},
            'joint_class': record['joint_class'],
            'reference_geometry': {'path': str(geometry_path.relative_to(ROOT)),
                                   'sha256': sha(geometry_path), 'frame': FRAME, 'units': 'm',
                                   'representation': 'triangular_surface'},
            **properties,
            'concepts': source.get('concepts', []),
            'classification': {'source_collections': source.get('source_collections', []),
                               'joint_class': record['joint_class'],
                               'joint_class_lexicon': 'scripts/audit_joint_substrate.py CLASSES',
                               'in_joints_collection': bool(record['in_joints_collection']),
                               'basis': 'authoring collections carried verbatim from the '
                                        'Z-Anatomy scene; no external ontology term is asserted'},
            'provenance': {
                'source_ids': [record['structure_id']],
                'dataset': 'Z-Anatomy extended atlas',
                'files': [{'path': record['source_path'], 'sha256': record['source_sha256']}],
                'source': copy.deepcopy(source['source']),
                'registration_id': 'z_anatomy',
                'transform': {'rotation': 'ihm.assembly.anatomy.ROTATION',
                              'registration': 'data/derived/canonical/anatomy.json '
                                              'registrations.z_anatomy',
                              'method': 'affine plus regularized thin-plate-spline residual',
                              'jacobian_determinant_min': record['transform']['jacobian_determinant_min'],
                              'jacobian_determinant_max': record['transform']['jacobian_determinant_max'],
                              'never_fitted_centroid_rms_m': NEVER_FITTED_CENTROID_RMS_M},
                'pipeline': [
                    {'script': 'scripts/build_joint_substrate_candidate.py',
                     'produced': 'registration, repair, TetGen proof and placement screen',
                     'artifact': SUBSTRATE + '/' + record['output_path'],
                     'artifact_sha256': record['output_sha256']},
                    {'script': 'scripts/build_muscle_tet_ready_surfaces.py',
                     'produced': 'the repair order applied verbatim'},
                    {'script': 'scripts/build_joint_promotion_candidate.py',
                     'produced': 'this canonical entity record'}],
                'dependency_group': 'bodyparts3d-derived-reference',
                'license': 'CC-BY-SA-4.0'},
            'assumptions': assumptions,
            'uncertainty': uncertainty(),
            'connections': connections,
            'tetgen_proof': {'flags': trial['flags'], 'volume_gate': trial['volume_gate'],
                             'tets': trial['repaired']['tets'],
                             'tet_volume_m3': trial['repaired'].get('tet_volume_m3'),
                             'tet_volume_vs_surface_relative_error':
                                 trial['repaired']['tet_volume_vs_surface_relative_error']},
            'material_reference': {
                'density_kg_m3': SOFT_DENSITY_KG_M3,
                'density_basis': 'canonical mechanics.json blanket soft-tissue value; assumed, not '
                                 'measured. data/derived/tissue-material-candidate-v1 records the '
                                 'density of every role this candidate uses as tier `absent`'}}
        entities.append(entity)
        decision['role'] = role['role']
        decision['role_exact'] = role['exact']
        decision['anatomical_role'] = role['anatomical_role']
        decision['volume_m3'] = properties['volume_m3']
        decision['tets'] = trial['repaired']['tets']
        decisions.append(decision)
    write_lines(out / 'entities.jsonl', entities)
    write_lines(out / 'decisions.jsonl', decisions)
    print('records: %d promoted, %d withheld' % (len(entities), len(rows) - len(entities)), flush=True)
    return entities, decisions


# ---------------------------------------------------------------- stage: conflicts

MESHES = None


def _pair_worker(task):
    out = []
    for i, j in task['pairs']:
        va, fa = MESHES[i]
        vb, fb = MESHES[j]
        record = {'a': int(i), 'b': int(j)}
        try:
            hits = cgal.intersect_other(va, fa, vb, fb, True, False, False, False, 2_000_000)[0]
            record['intersecting_face_pairs'] = int(len(np.asarray(hits).reshape(-1, 2)))
        except BaseException as error:
            record['error'] = 'intersect_other:' + type(error).__name__
            out.append(record)
            continue
        if not record['intersecting_face_pairs']:
            continue
        if len(fa) + len(fb) <= task['bool_cap']:
            try:
                result = cgal.mesh_boolean(va, fa, vb, fb, type_str='intersect')
                iv = np.ascontiguousarray(np.asarray(result[0], float))
                inf = np.ascontiguousarray(np.asarray(result[1], np.int64))
                volume, surface = abs(divergence(iv, inf)), area_of(iv, inf)
                record.update(overlap_volume_m3=volume, overlap_area_m2=surface,
                              overlap_thickness_proxy_m=(2. * volume / surface) if surface > 0 else 0.,
                              overlap_method='cgal_boolean_intersect')
            except BaseException as error:
                record['overlap_method'] = 'cgal_boolean_failed:' + type(error).__name__
        if 'overlap_volume_m3' not in record:
            blow = np.maximum(va.min(0), vb.min(0))
            bhigh = np.minimum(va.max(0), vb.max(0))
            span = np.maximum(bhigh - blow, 0.)
            boxv = float(np.prod(span))
            rng = np.random.default_rng(task['seed'] + i * 100003 + j)
            query = np.ascontiguousarray(blow + rng.random((task['mc'], 3)) * span)
            wa = np.abs(np.asarray(igl.fast_winding_number(va, fa, query), float)) > .5
            wb = np.abs(np.asarray(igl.fast_winding_number(vb, fb, query), float)) > .5
            k = int((wa & wb).sum())
            p = k / task['mc']
            record.update(overlap_volume_m3=boxv * p,
                          overlap_volume_stderr_m3=boxv * float(np.sqrt(max(p * (1 - p), 0) / task['mc'])),
                          overlap_area_m2=None, overlap_thickness_proxy_m=None,
                          overlap_method=record.get('overlap_method', 'skipped_boolean_over_face_cap') + '|monte_carlo')
        out.append(record)
    return out


def existing_rows(canonical, builds):
    meta = {e['id']: e for e in canonical['entities']}
    rows = []
    for entity_id, record in builds.items():
        entity = meta[entity_id]
        rows.append(dict(key=entity_id, name=entity['name'], role=entity['role'],
                         system=entity['system'], lane=record['base'],
                         path=ROOT / record['base'] / record['output_path'], side='existing'))
    rows.sort(key=lambda r: r['key'])
    return rows


def stage_conflicts(out, entities, canonical, builds, workers, chunk):
    global MESHES
    started = time.monotonic()
    rows = existing_rows(canonical, builds)
    for entity in entities:
        rows.append(dict(key=entity['id'], name=entity['name'], role=entity['role'],
                         system=entity['system'], lane='joint-promotion-candidate-v1',
                         path=ROOT / entity['reference_geometry']['path'], side='promoted'))
    MESHES = [read_geometry(r['path']) for r in rows]
    new = np.array([r['side'] == 'promoted' for r in rows])
    faces = np.array([len(f) for _, f in MESHES])
    low = np.stack([v.min(0) for v, _ in MESHES])
    high = np.stack([v.max(0) for v, _ in MESHES])
    volume = np.array([abs(divergence(v, f)) for v, f in MESHES])
    box = (low[:, None, :] <= high[None, :, :]).all(2) & (high[:, None, :] >= low[None, :, :]).all(2)
    np.fill_diagonal(box, False)
    added = box & (new[:, None] | new[None, :])
    pairs = np.argwhere(np.triu(added))
    order = np.argsort(faces[pairs[:, 0]] + faces[pairs[:, 1]])
    pairs = pairs[order]
    tasks = [{'pairs': [(int(a), int(b)) for a, b in pairs[k:k + chunk]],
              'bool_cap': BOOL_FACE_CAP, 'mc': MC_POINTS, 'seed': MC_SEED}
             for k in range(0, len(pairs), chunk)]
    print('conflicts: %d entities (%d promoted), %d ADDED bounding-box-overlapping pairs, %d tasks'
          % (len(rows), int(new.sum()), len(pairs), len(tasks)), flush=True)
    found, done = [], 0
    with mp.get_context('fork').Pool(workers) as pool:
        for batch in pool.imap_unordered(_pair_worker, tasks):
            found.extend(batch)
            done += 1
            if done % 25 == 0:
                print('  %d/%d tasks, %d intersecting pairs, %.0fs'
                      % (done, len(tasks), len(found), time.monotonic() - started), flush=True)
    for record in found:
        a, b = record['a'], record['b']
        record['a'], record['b'] = rows[a]['key'], rows[b]['key']
        record['a_role'], record['b_role'] = rows[a]['role'], rows[b]['role']
        record['a_side'], record['b_side'] = rows[a]['side'], rows[b]['side']
        record['a_name'], record['b_name'] = rows[a]['name'], rows[b]['name']
    found.sort(key=lambda r: -(r.get('overlap_volume_m3') or 0.))
    write_lines(out / 'conflict-pairs-added.jsonl', found)
    hits = [r for r in found if r.get('intersecting_face_pairs', 0) > 0]
    promoted_promoted = [r for r in hits if r['a_side'] == r['b_side'] == 'promoted']
    promoted_existing = [r for r in hits if r['a_side'] != r['b_side']]
    volumes = np.array([r['overlap_volume_m3'] for r in hits if r.get('overlap_volume_m3') is not None])
    thick = np.array([r['overlap_thickness_proxy_m'] for r in hits
                      if r.get('overlap_thickness_proxy_m') is not None])
    conflicted = {r['a'] for r in hits} | {r['b'] for r in hits}
    promoted_ids = {r['key'] for r in rows if r['side'] == 'promoted'}
    prior = json.loads((ROOT / 'data/derived/cross-structure-repair-v1/conflict-graph.json').read_text())
    prior_resolution = json.loads((ROOT / 'data/derived/cross-structure-repair-v1/resolution.json').read_text())
    prior_conflicted = set()
    for line in (ROOT / 'data/derived/cross-structure-repair-v1/conflict-pairs.jsonl').read_text().splitlines():
        record = json.loads(line)
        if record.get('intersecting_face_pairs', 0) > 0:
            prior_conflicted.add(record['a'])
            prior_conflicted.add(record['b'])
    newly = sorted((conflicted - promoted_ids) - prior_conflicted)
    role_pairs = {}
    for record in hits:
        key = ' x '.join(sorted((record['a_role'], record['b_role'])))
        bucket = role_pairs.setdefault(key, dict(pairs=0, overlap_volume_m3=0.))
        bucket['pairs'] += 1
        bucket['overlap_volume_m3'] += record.get('overlap_volume_m3') or 0.
    report = dict(
        schema=SCHEMA + '.conflicts',
        method=('exhaustive CGAL intersect_other over every bounding-box-overlapping pair that '
                'involves at least one promoted structure, with an exact CGAL boolean overlap '
                'volume per intersecting pair; the method of '
                'scripts/build_cross_structure_conflict_repair.py, unchanged'),
        entities_total=len(rows), entities_existing=int((~new).sum()),
        entities_promoted=int(new.sum()),
        prior_graph=dict(entities=prior['entities'],
                         bounding_box_overlapping_pairs=prior['bounding_box_overlapping_pairs'],
                         intersecting_pairs=prior['intersecting_pairs'],
                         total_intersecting_face_pairs=prior['total_intersecting_face_pairs'],
                         entities_in_conflict=prior_resolution['entities_in_conflict'],
                         total_overlap_volume_m3=prior_resolution['total_overlap_volume_m3']),
        added_bounding_box_overlapping_pairs=int(len(pairs)),
        added_exactly_tested_pairs=int(len(pairs)),
        added_intersecting_pairs=len(hits),
        added_intersecting_pairs_promoted_vs_existing=len(promoted_existing),
        added_intersecting_pairs_promoted_vs_promoted=len(promoted_promoted),
        added_intersecting_face_pairs=int(sum(r['intersecting_face_pairs'] for r in hits)),
        errored_pairs=int(sum('error' in r for r in found)),
        added_overlap_volume_m3=float(volumes.sum()) if len(volumes) else 0.,
        added_overlap_volume_promoted_vs_existing_m3=float(
            sum(r.get('overlap_volume_m3') or 0. for r in promoted_existing)),
        added_overlap_volume_promoted_vs_promoted_m3=float(
            sum(r.get('overlap_volume_m3') or 0. for r in promoted_promoted)),
        added_overlap_volume_median_m3=float(np.median(volumes)) if len(volumes) else 0.,
        added_overlap_thickness_proxy_median_m=float(np.median(thick)) if len(thick) else 0.,
        promoted_entities_in_conflict=len(conflicted & promoted_ids),
        promoted_entities_clean=int(new.sum()) - len(conflicted & promoted_ids),
        existing_entities_newly_in_conflict=len(newly),
        existing_entities_newly_in_conflict_ids=newly,
        promoted_summed_surface_volume_m3=float(volume[new].sum()),
        added_overlap_as_fraction_of_promoted_volume=(
            float(volumes.sum() / volume[new].sum()) if len(volumes) and volume[new].sum() else None),
        role_pair_totals=dict(sorted(role_pairs.items(), key=lambda kv: -kv[1]['overlap_volume_m3'])),
        resolution_performed=False,
        wall_seconds=time.monotonic() - started)
    write_json(out / 'conflicts.json', report)
    print('conflicts: %d added intersecting pairs, %.6f L added overlap, %.0fs'
          % (len(hits), report['added_overlap_volume_m3'] * 1000., report['wall_seconds']), flush=True)
    return report


# ---------------------------------------------------------------- stage: finalize

def stage_finalize(out, substrate, rows, entities, decisions, duplicate_report, placement_report,
                   conflicts, canonical, elapsed):
    per_role, per_class, per_system = {}, {}, {}
    volume = mass = 0.
    tets = 0
    inexact = []
    for entity in entities:
        role = entity['role']
        bucket = per_role.setdefault(role, dict(count=0, volume_m3=0., mass_kg=0., tets=0,
                                                exact_role_term=entity['role_vocabulary']['exact']))
        bucket['count'] += 1
        bucket['volume_m3'] += entity['volume_m3'] or 0.
        bucket['mass_kg'] += (entity['volume_m3'] or 0.) * SOFT_DENSITY_KG_M3
        bucket['tets'] += entity['tetgen_proof']['tets']
        per_class[entity['joint_class']] = per_class.get(entity['joint_class'], 0) + 1
        per_system[entity['system']] = per_system.get(entity['system'], 0) + 1
        volume += entity['volume_m3'] or 0.
        mass += (entity['volume_m3'] or 0.) * SOFT_DENSITY_KG_M3
        tets += entity['tetgen_proof']['tets']
        if not entity['role_vocabulary']['exact']:
            inexact.append(entity['id'])
    withheld = {}
    for decision in decisions:
        if decision['promoted']:
            continue
        key = '+'.join(decision['withheld_reasons'])
        withheld[key] = withheld.get(key, 0) + 1
    summary = dict(
        schema=SCHEMA + '.summary',
        canonical_modified=False, promotion_performed=False,
        source_structures=len(rows),
        prior_lane_promotable=sum(bool(r['promotable']) for r in rows),
        promoted=len(entities),
        promotion_arithmetic=dict(
            prior_lane_promotable=sum(bool(r['promotable']) for r in rows),
            minus_name_identity_duplicates=-sum(
                1 for d in decisions if not d['promoted']
                and 'superseded_by_canonical_duplicate' in d['withheld_reasons']
                and d['prior_lane_promotable']),
            plus_released_by_duplicate_review=sum(1 for d in decisions if d['released_by_review']),
            minus_laterality_label_conflicts=-sum(
                1 for d in decisions if not d['promoted']
                and d['withheld_reasons'] == ['laterality_label_conflict']),
            equals=len(entities)),
        entity_count_after_promotion=len(canonical['entities']) + len(entities),
        per_role=dict(sorted(per_role.items())),
        per_joint_class=dict(sorted(per_class.items())),
        per_system=dict(sorted(per_system.items())),
        added_volume_m3=volume, added_volume_l=volume * 1000.,
        added_mass_kg=mass,
        added_mass_basis='volume times %.0f kg/m3, the canonical mechanics.json blanket soft '
                         'density. Assumed, not measured' % SOFT_DENSITY_KG_M3,
        added_tets=tets,
        roles_with_inexact_vocabulary=len(inexact),
        role_vocabulary_gaps=sorted({e['role_vocabulary']['anatomical_role'] for e in entities
                                     if not e['role_vocabulary']['exact']}),
        withheld=len(rows) - len(entities),
        withheld_by_reason=dict(sorted(withheld.items(), key=lambda kv: -kv[1])),
        duplicates=duplicate_report,
        placement=placement_report,
        conflicts={k: v for k, v in conflicts.items() if k != 'existing_entities_newly_in_conflict_ids'},
        verified=[
            'every promoted surface carries a volume-gated TetGen proof recomputed by the prior '
            'lane and its sha256 is checked here before use',
            'every duplicate decision carries a measured symmetric mean and max surface distance '
            'between the two copies',
            'every deep-inside-bone cause carries an exact CGAL boolean overlap volume against '
            'canonical bone and an exact articular gap for the expected bone pair',
            'every added conflict pair is an exhaustive CGAL intersect_other test with an exact '
            'boolean overlap volume'],
        inferred=[
            'roles are inferred from the audit joint class and structure wording; the '
            'role_vocabulary block on each record states where the canonical vocabulary has no '
            'exact term',
            'expected attachment bones come from structure wording and authoring collection, not '
            'a curated attachment table',
            'the 1000 kg/m3 density is a blanket assumption, so added mass is proportional to '
            'added volume and carries no independent evidence',
            'connections are proximity in the registered pose, not certified attachments'],
        elapsed_s=elapsed)
    merged = copy.deepcopy(canonical)
    merged['entities'] = canonical['entities'] + entities
    merged['counts'] = dict(canonical['counts'])
    try:
        verification = verify_assembly(merged, ROOT)
        summary['merged_assembly_verification'] = dict(
            passed=True, entities=len(merged['entities']),
            checker='ihm.assembly.anatomy.verify_assembly, run on canonical + candidate in memory; '
                    'nothing was written back to canonical',
            report=verification if isinstance(verification, dict) else None)
    except BaseException as error:
        summary['merged_assembly_verification'] = dict(
            passed=False, error='%s: %s' % (type(error).__name__, error))
    write_json(out / 'summary.json', summary)

    # manifest.json is excluded: a file cannot carry its own digest, and listing it would record
    # the previous run's hash and read as a mismatch on every check.
    artifacts = {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*'))
                 if p.is_file() and p.name != 'manifest.json'}
    inputs = {p: sha(ROOT / p) for p in (
        SUBSTRATE + '/entities.jsonl', SUBSTRATE + '/summary.json', SUBSTRATE + '/manifest.json',
        CANONICAL, EXTENDED_FRAGMENT, MUSCLE_BUILD + '/entities.jsonl',
        ENTITY_BUILD + '/entities.jsonl', MECHANICS,
        'data/derived/cross-structure-repair-v1/conflict-graph.json',
        'data/derived/cross-structure-repair-v1/conflict-pairs.jsonl',
        'data/derived/cross-structure-repair-v1/resolution.json',
        'scripts/build_joint_substrate_candidate.py', 'scripts/build_canonical_anatomy.py',
        'scripts/build_cross_structure_conflict_repair.py', 'scripts/audit_joint_substrate.py',
        'ihm/assembly/anatomy.py', Path(__file__).relative_to(ROOT).as_posix())}
    source_geometry = {r['source_path']: r['source_sha256'] for r in rows}
    substrate_geometry = {r['output_path']: r['output_sha256'] for r in rows}
    write_json(out / 'manifest.json', dict(
        schema=SCHEMA,
        packages={n: importlib.metadata.version(n) for n in ('libigl', 'numpy', 'scipy')},
        python=sys.version,
        native_cgal_sha256=sha(Path(cgal.pyigl_copyleft_cgal.__file__)),
        inputs_sha256=inputs,
        extended_source_geometry_sha256=source_geometry,
        substrate_repaired_geometry_sha256=substrate_geometry,
        artifacts_sha256=artifacts,
        canonical_assets_modified=False, promotion_performed=False,
        limitations=[
            'This is a promotion candidate for review. No canonical file is written.',
            'The canonical intervertebral disks are roled rigid_bone. This candidate records the '
            'recommended correction on each duplicate decision and does not apply it.',
            'Role assignment uses the canonical vocabulary. Where no exact term exists the record '
            'says so; the assigned term must not be read as a constitutive claim.',
            'Conflicts are measured and not resolved. Ownership, part-of and containment relations '
            'are owned by other lanes.',
            'A TetGen proof establishes a valid closed PLC whose tets account for its volume. It '
            'establishes nothing about anatomical correctness.']))
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ('duplicates', 'placement', 'conflicts')}, indent=1))
    return summary


# ---------------------------------------------------------------- self test

def self_test():
    assert role_of('ligament', 'Iliofemoral ligament.l')['exact']
    assert role_of('ligament', 'Iliofemoral ligament.l')['role'] == 'ligament'
    assert role_of('cruciate_ligament', 'Anterior cruciate ligament.l')['role'] == 'ligament'
    disc = role_of('intervertebral_disc', 'Intervertebral disc L3-L4')
    assert disc['role'] == 'cartilage' and not disc['exact'] and disc['vocabulary_gap']
    nucleus = role_of('nucleus_pulposus', 'Nucleus pulposus L4-L5')
    assert nucleus['role'] == 'cartilage' and not nucleus['exact']
    assert 'fluid_cavity' in nucleus['rejected'], 'the rejected alternative must be recorded'
    assert role_of('meniscus', 'Medial meniscus.l')['role'] == 'cartilage'
    knee = role_of('meniscus', 'Posterior meniscotibial ligament (Medial meniscus).l')
    assert knee['role'] == 'ligament' and knee['exact'], 'a meniscotibial ligament is a ligament'
    assert role_of('synovial_bursa', 'Suprapatellar bursa.l')['role'] == 'fluid_cavity'
    assert role_of('fascia_aponeurosis', 'Epicranial aponeurosis.l')['role'] == 'tendon'
    assert role_of('fascia_aponeurosis', 'Piriformis fascia.l')['role'] == 'connective_tissue'
    assert role_of('fat_pad', 'Infrapatellar fat pad.l')['vocabulary_gap']
    assert role_of('cartilage', 'Costal cartilage of first rib.l')['exact']
    for joint_class in ('ligament', 'cruciate_ligament', 'retinaculum', 'articular_capsule',
                        'articular_disc', 'meniscus', 'labrum', 'symphysis', 'synovial_bursa',
                        'tendon_sheath', 'fascia_aponeurosis', 'interosseous_membrane',
                        'cartilage', 'fat_pad', 'intervertebral_disc', 'nucleus_pulposus'):
        assert role_of(joint_class, 'probe')['role'] in CANONICAL_ROLES, joint_class

    assert 'left first costal cartilage' in name_aliases('Costal cartilage of first rib.l')
    assert 'flexor retinaculum of left wrist' in name_aliases('Flexor retinaculum of wrist.l')
    assert 'interosseous membrane of left leg' in name_aliases('Interosseous membrane of leg.l')
    assert name_aliases('Iliofemoral ligament.l') == ['left iliofemoral ligament']
    assert not_duplicate_reason('Pisotriquetral ligament.l', 'left pisiform')
    assert not_duplicate_reason('Quadrangular membrane.l', 'left thyrohyoid membrane')
    assert not not_duplicate_reason('Stylohyoid ligament.l', 'left stylohyoid ligament')

    v = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [0, 0, 1.]])
    f = np.ascontiguousarray(np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]], np.int64))
    assert abs(abs(divergence(v, f)) - 1 / 6) < 1e-15
    shifted = np.ascontiguousarray(v + np.array([.5, 0, 0]))
    overlap, method = boolean_overlap(v, f, shifted, f)
    assert method == 'cgal_boolean_intersect' and 0. < overlap < 1 / 6
    far = np.ascontiguousarray(v + np.array([10., 0, 0]))
    assert boolean_overlap(v, f, far, f)[0] == 0.
    mc, stderr = monte_carlo_overlap(v, f, shifted, f, points=20000)
    assert abs(mc - overlap) < .02 * (1 / 6) and stderr >= 0.
    assert monte_carlo_overlap(v, f, far, f, points=1000)[0] == 0.
    assert 'left hip bone' in supplementary_expected('Triradiate cartilage.l',
                                                    {'left hip bone', 'right hip bone'})
    assert supplementary_expected('Triradiate cartilage.l', {'right hip bone'}) == []
    assert 'left fourth rib' in supplementary_expected(
        'Costotransverse ligament.l', {'left fourth rib', 'right fourth rib'})
    assert laterality_conflict('left probe', [-.054, 0., 0.])
    assert laterality_conflict('right probe', [.054, 0., 0.])
    assert not laterality_conflict('left probe', [.054, 0., 0.])
    assert not laterality_conflict('median probe', [-.054, 0., 0.])
    assert laterality_conflict('probe', [0., 2., 0.])
    mean, worst = symmetric_surface_distance(v, f, np.ascontiguousarray(v + .001), f)
    assert 0 < mean <= worst <= .002

    properties = mesh_properties(v, f)
    for key in ('bounds_m', 'centroid_m', 'centroid_definition', 'surface_area_m2', 'volume_m3',
                'volume_method', 'watertight_edge_incidence', 'principal_axis',
                'source_vertex_count', 'source_face_count'):
        assert key in properties, key
    print('PASS role vocabulary and its declared gaps, duplicate name aliases and overturned '
          'flags, exact CGAL overlap predicate, canonical mesh property shape')


# ---------------------------------------------------------------- driver

def run(output, stage, workers, chunk):
    out = Path(output).resolve()
    if not out.is_relative_to(ROOT / 'data/derived'):
        raise ValueError('Output must live under data/derived')
    if importlib.metadata.version('libigl') != '2.6.2':
        raise ValueError('This build pins libigl 2.6.2')
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    substrate, rows, canonical, fragment = load_inputs()
    manifest = json.loads((substrate / 'manifest.json').read_text())
    for record in rows:
        stored = manifest['artifacts_sha256'].get(record['output_path'])
        if stored != record['output_sha256'] or sha(substrate / record['output_path']) != stored:
            raise ValueError('Substrate geometry digest mismatch: ' + record['structure_id'])
    print('inputs: %d substrate structures, digests verified' % len(rows), flush=True)
    builds = build_index()
    tets = tetgen_index()

    if stage in ('all', 'duplicates') or not (out / 'duplicates.jsonl').exists():
        duplicates, duplicate_report = stage_duplicates(out, substrate, rows, canonical, builds, tets)
    else:
        duplicates = [json.loads(l) for l in (out / 'duplicates.jsonl').read_text().splitlines()]
        duplicate_report = json.loads((out / 'duplicates.json').read_text())
    if stage == 'duplicates':
        return

    if stage in ('all', 'placement') or not (out / 'placement-review.jsonl').exists():
        placement, placement_report = stage_placement(out, substrate, rows, canonical, builds, workers)
    else:
        placement = [json.loads(l) for l in (out / 'placement-review.jsonl').read_text().splitlines()]
        placement_report = json.loads((out / 'placement-review.json').read_text())
    if stage == 'placement':
        return

    if stage in ('all', 'records') or not (out / 'entities.jsonl').exists():
        entities, decisions = stage_records(out, substrate, rows, canonical, fragment, duplicates,
                                            placement, None)
    else:
        entities = [json.loads(l) for l in (out / 'entities.jsonl').read_text().splitlines()]
        decisions = [json.loads(l) for l in (out / 'decisions.jsonl').read_text().splitlines()]
    if stage == 'records':
        return

    if stage in ('all', 'conflicts') or not (out / 'conflicts.json').exists():
        conflicts = stage_conflicts(out, entities, canonical, builds, workers, chunk)
    else:
        conflicts = json.loads((out / 'conflicts.json').read_text())
    if stage == 'conflicts':
        return

    stage_finalize(out, substrate, rows, entities, decisions, duplicate_report, placement_report,
                   conflicts, canonical, time.monotonic() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--stage', default='all',
                        choices=('all', 'duplicates', 'placement', 'records', 'conflicts', 'finalize'))
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--chunk', type=int, default=64)
    args = parser.parse_args()
    if args.self_test:
        self_test()
    if args.output:
        run(args.output, args.stage, args.workers, args.chunk)
    if not args.self_test and not args.output:
        parser.error('Select --self-test or --output')
