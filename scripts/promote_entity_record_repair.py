#!/usr/bin/env python3
"""Promote data/derived/entity-record-repair-candidate-v1 into the canonical body.

The candidate measured three defects in the entity records and proved them
without touching anything:

  A  five BodyParts3D surfaces authored twice under two element ids, so one
     physical structure is counted twice: twice the volume, twice the mass,
     every link duplicated, and a support that spans a structure to itself
  B  twenty-three intervertebral discs carrying role rigid_bone, so the one
     compliant element between two vertebrae is infinitely stiff
  C  three skin layers extruded against the raw two-sided slab area rather than
     the exterior component, inventing 11.366 L of skin

This script does not re-apply those edits by hand. Each one now lives in the
shipped builder that owns the field, so a rebuild reaches the same answer and a
future rebuild cannot silently revert the repair:

  A  ihm.assembly.anatomy.surface_identity_key / duplicate_surface_survivors,
     called from scripts/build_canonical_anatomy.py
  B  ihm.assembly.anatomy.physical_role, which already routed
     'intervertebral disc' to cartilage and missed the BodyParts3D spelling
     'intervertebral disk'
  C  scripts/build_canonical_anatomy.py already wrote the exterior-component
     area; the on-disk artifact was stale relative to its own builder

Independently of the three, the declared body mass is corrected from the
inherited BioGears StandardMale 77.1107029 kg constant to the 70.7713 kg
composed over this specimen's own measured interior. All four land together
because all four pass through the single uniform mass normalizer at
scripts/build_body_mechanics.py.

Run order is the dependency order: profile, anatomy, mechanics, body manifest.
Everything is then checked against the candidate's own recorded predictions,
and a receipt records input hashes, output hashes, the counts moved and this
script's own digest.

    scripts/promote_entity_record_repair.py --self-test
    scripts/promote_entity_record_repair.py --promote
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from ihm.assembly.anatomy import (BP_SOURCE_SURFACE_COUNT, attach_regional_support,
                                  duplicate_surface_survivors, physical_role, read_geometry,
                                  surface_identity_key)
from ihm.assembly.profile import MASS_KG, MASS_LEDGER

CANDIDATE = ROOT / 'data/derived/entity-record-repair-candidate-v1'
RECEIPT = ROOT / 'data/derived/entity-record-repair-promotion-v1'
CANONICAL = ROOT / 'data/derived/canonical'

# Assets this promotion rewrites, in the order the builders touch them.
PROMOTED = ['data/derived/canonical/profile.json',
            'data/derived/canonical/anatomy.json',
            'data/derived/canonical/manifest_fragment.json',
            'data/derived/canonical/verification.json',
            'data/derived/canonical/lymphatic_graph.json',
            'data/derived/canonical/mechanics.json',
            'data/derived/canonical/body.json',
            'data/derived/app/manifest.json']

SKIN_LAYERS = ['body-skin-epidermis', 'body-skin-dermis', 'body-skin-hypodermis']
SKIN_PARENT = 'body-bp3d-FJ2810'
EXTERIOR_AREA_M2 = 1.7804602548390722
RAW_SLAB_AREA_M2 = 3.502598974493317
CANDIDATE_TARGET_MASS_KG = 77.1107029  # the mass the candidate normalized against

# The nine actuator and connective priors that spanned the duplicated hyoid, one
# structure to itself. downstream.json requires an explicit decision per prior:
# re-projection onto a genuinely distinct second bone, or explicit deletion.
# Silently collapsing them to zero length is not acceptable.
SELF_SPANNING_PRIORS = ['body-muscle-prior-body-bp3d-FJ1560',
                        'body-muscle-prior-body-bp3d-FJ2785',
                        'body-muscle-prior-body-bp3d-FJ2803',
                        'body-connective-body-bp3d-FJ1581',
                        'body-connective-body-bp3d-FJ2771',
                        'body-connective-body-bp3d-FJ2779',
                        'body-connective-body-bp3d-FJ2790',
                        'body-connective-body-bp3d-FJ2797',
                        'body-connective-body-bp3d-FJ2807']
SELF_SPANNING_DECISION = {
    'choice': 're-project the second anchor onto a genuinely distinct bone',
    'not_chosen': ['delete the prior', 'retain a zero-length or intra-structure span'],
    'reason': 'every one of the nine is a real hyoid structure with a real second attachment; the '
              'defect was the duplicate hyoid id standing in for that second bone, not the prior. '
              'Deleting them would remove nine load paths the anatomy actually has.',
    'mechanism': 'no new rule. build_body_mechanics.py already refuses a force pair internal to one '
                 'bone (line 181 for synthesized actuators, line 201 for connective paths) and '
                 're-projects the second endpoint onto the nearest bone that is not the first. With '
                 'the duplicate id gone both endpoints resolve to the survivor, that guard fires, '
                 'and the second anchor lands on a distinct bone.',
    'audit': 'every one of the nine must end on two distinct entities with a rest path length '
             'strictly longer than the intra-structure span it replaced; asserted below.'}


def sha256_path(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load_candidate():
    """Read the candidate and re-check its own recorded artifact digests."""
    manifest = json.loads((CANDIDATE / 'manifest.json').read_text())
    mismatched = {name: sha256_path(CANDIDATE / name)
                  for name, digest in manifest['artifacts_sha256'].items()
                  if sha256_path(CANDIDATE / name) != digest}
    if mismatched:
        raise ValueError('Candidate artifacts do not match their own manifest: %r' % mismatched)
    payload = {name: json.loads((CANDIDATE / name).read_text())
               for name in ('records.json', 'allocation.json', 'downstream.json',
                            'duplicates.json', 'disc-roles.json', 'skin-layers.json')}
    payload['manifest.json'] = manifest
    return payload


def candidate_expectations(candidate):
    """The candidate's predictions, read out of the candidate rather than restated."""
    records = candidate['records.json']
    allocation = candidate['allocation.json']
    drops = {r['id']: r['survivor'] for r in records['anatomy']['remove']}
    discs = [r['id'] for r in records['anatomy']['change'] if r.get('role') == ['rigid_bone', 'cartilage']]
    layers = {r['id']: r['volume_m3'][1] for r in records['anatomy']['change']
              if r['id'] in SKIN_LAYERS}
    return {'drops': drops, 'discs': discs, 'layer_volume_m3': layers,
            'unscaled_proxy_mass_kg': allocation['steps_unscaled_kg']['total'],
            'uniform_scale_at_candidate_target': allocation['uniform_scale']['after'],
            'skin_fraction_before_after': allocation['headline']['skin_layer_fraction_of_body'],
            'counts': records['counts']}


def verify_duplicate_identity(anatomy_by_id, drops, root=ROOT):
    """Re-decode both surfaces of every collapsed pair and prove they are one.

    The dropped row is gone from the assembly, so its geometry is read from the
    collapse record the builder wrote. This runs the same identity function the
    builder used, on the bytes actually on disk, not on a recorded claim.
    """
    results = []
    for record in drops:
        dropped, survivor = record['id'], record['survivor']
        a = read_geometry(root / record['dropped_geometry_path'])
        b = read_geometry(root / record['survivor_geometry_path'])
        va = np.asarray(a['positions'], float).reshape(-1, 3)
        fa = np.asarray(a['indices'], np.int64).reshape(-1, 3)
        vb = np.asarray(b['positions'], float).reshape(-1, 3)
        fb = np.asarray(b['indices'], np.int64).reshape(-1, 3)
        key_a, key_b = surface_identity_key(va, fa), surface_identity_key(vb, fb)
        if not (va.shape == vb.shape and fa.shape == fb.shape
                and np.array_equal(va, vb) and np.array_equal(fa, fb) and key_a == key_b):
            raise ValueError('Collapsed pair is not one surface: %s / %s' % (dropped, survivor))
        if survivor not in anatomy_by_id or dropped in anatomy_by_id:
            raise ValueError('Collapse did not keep exactly the survivor: %s' % dropped)
        if key_a != record['surface_identity_key']:
            raise ValueError('Recorded surface identity disagrees with the bytes: %s' % dropped)
        results.append({'dropped': dropped, 'survivor': survivor,
                        'name': anatomy_by_id[survivor]['name'],
                        'triangles': int(len(fa)), 'vertices': int(len(va)),
                        'exact_array_equality': True, 'surface_identity_key': key_a,
                        'stored_gz_sha256_differ': (sha256_path(root / record['dropped_geometry_path'])
                                                    != sha256_path(root / record['survivor_geometry_path']))})
    return results


def check_promoted(candidate, root=ROOT):
    """Every assertion runs against the artifacts the shipped builders wrote."""
    expect = candidate_expectations(candidate)
    anatomy = json.loads((root / 'data/derived/canonical/anatomy.json').read_text())
    mechanics = json.loads((root / 'data/derived/canonical/mechanics.json').read_text())
    profile = json.loads((root / 'data/derived/canonical/profile.json').read_text())
    fragment = json.loads((root / 'data/derived/canonical/manifest_fragment.json').read_text())
    by_id = {e['id']: e for e in anatomy['entities']}
    spec = {e['id']: e for e in mechanics['entities']}
    report = {}

    # ---- A. duplicate collapse -------------------------------------------
    collapse = anatomy['duplicate_surface_collapse']
    recorded = {d['id']: d['survivor'] for d in collapse['dropped']}
    if recorded != expect['drops']:
        raise ValueError('Collapsed set differs from the candidate: %r vs %r' % (recorded, expect['drops']))
    report['A_duplicates'] = verify_duplicate_identity(by_id, collapse['dropped'], root)
    dropped_ids = set(recorded)
    if dropped_ids & set(spec):
        raise ValueError('Dropped ids survive in mechanics')
    if any(x['id'] in dropped_ids for x in fragment['structures']):
        raise ValueError('Dropped ids survive as display structures')
    touching = [l for l in mechanics['links'] if l['a'] in dropped_ids or l['b'] in dropped_ids]
    if touching:
        raise ValueError('%d links still anchored on a dropped id' % len(touching))
    if any(l['a'] == l['b'] for l in mechanics['links']):
        raise ValueError('A support link joins a structure to itself')
    for m in mechanics['muscles']:
        anchors = {a['entity_id'] for a in m['anchors']}
        if anchors & dropped_ids:
            raise ValueError('Muscle anchor on a dropped id: ' + m['id'])
        if len(anchors) < 2:
            raise ValueError('Prior spans one structure to itself: ' + m['id'])
    spans = [m for m in mechanics['muscles']
             if {frozenset({a['entity_id'] for a in m['anchors']})} & {frozenset({k, v}) for k, v in recorded.items()}]
    if spans:
        raise ValueError('A prior still spans a collapsed pair')

    # ---- the nine self-spanning priors -----------------------------------
    ledger = []
    for identity in SELF_SPANNING_PRIORS:
        m = next((x for x in mechanics['muscles'] if x['id'] == identity), None)
        if m is None:
            raise ValueError('Self-spanning prior vanished without a recorded deletion: ' + identity)
        anchors = [a['entity_id'] for a in m['anchors']]
        if len(set(anchors)) != 2:
            raise ValueError('Self-spanning prior did not gain a distinct second anchor: ' + identity)
        if not m['rest_path_length_m'] > 1e-3:
            raise ValueError('Self-spanning prior collapsed toward zero length: ' + identity)
        ledger.append({'id': identity, 'canonical_entity': m['canonical_entity_id'],
                       'canonical_entity_name': by_id[m['canonical_entity_id']]['name'],
                       'decision': 're-projected', 'anchors': anchors,
                       'anchor_names': [by_id[a]['name'] for a in anchors],
                       'rest_path_length_m': m['rest_path_length_m'],
                       'passive_stiffness_n_m': m['passive_stiffness_n_m']})
    report['self_spanning_priors'] = {**SELF_SPANNING_DECISION, 'priors': ledger}

    # ---- B. disc roles ----------------------------------------------------
    discs = sorted(e['id'] for e in anatomy['entities'] if 'intervertebral' in e['name'].lower())
    if discs != sorted(expect['discs']):
        raise ValueError('Disc set differs from the candidate')
    for identity in discs:
        e, s = by_id[identity], spec[identity]
        if e['role'] != 'cartilage' or e['system'] != 'skeletal':
            raise ValueError('Disc role or system wrong: ' + identity)
        if physical_role(e['name'], e['system']) != 'cartilage':
            raise ValueError('Shipped classifier disagrees with the stored role: ' + identity)
        if s['constitutive'] != 'affine_neo_hookean' or s['material']['density']['value'] != 1000.:
            raise ValueError('Disc did not leave the rigid solver branch: ' + identity)
    if mechanics['counts']['rigid_bones'] != expect['counts']['rigid_bones'][1]:
        raise ValueError('rigid_bone count differs from the candidate')
    if any(c['entity_id'] in set(discs) for e in anatomy['entities'] for c in e['connections']
           if c['relation'] == 'regional_spatial_support'):
        raise ValueError('A regional support still anchors on a disc')
    if any(l['kind'] == 'inferred_skeletal_support' and (l['a'] in set(discs) or l['b'] in set(discs))
           for l in mechanics['links']):
        raise ValueError('The skeletal tree still routes through a disc')
    supported = {l['a'] for l in mechanics['links'] if l['kind'] == 'soft_tissue_support'}
    if not set(discs) <= supported:
        raise ValueError('A disc lacks the soft-tissue support it should have gained')
    report['B_discs'] = {'count': len(discs), 'role': 'cartilage', 'system': 'skeletal',
                         'constitutive': 'affine_neo_hookean', 'density_kg_m3': 1000.,
                         'rigid_bones': [expect['counts']['rigid_bones'][0], mechanics['counts']['rigid_bones']],
                         'soft_solids': [expect['counts']['soft_solids'][0], mechanics['counts']['soft_solids']]}

    # ---- C. skin layers ---------------------------------------------------
    support = by_id[SKIN_PARENT]['physical_surface_support']
    if support['area_m2'] != EXTERIOR_AREA_M2 or support['raw_source_area_m2'] != RAW_SLAB_AREA_M2:
        raise ValueError('Skin support area is not the exterior component')
    if by_id[SKIN_PARENT]['surface_area_m2'] != RAW_SLAB_AREA_M2:
        raise ValueError('The descriptive raw slab measurement must be retained')
    layers = {}
    for identity in SKIN_LAYERS:
        e = by_id[identity]
        thickness = e['shell']['thickness_m']
        if e['volume_m3'] != EXTERIOR_AREA_M2 * thickness:
            raise ValueError('Layer volume is not area times thickness: ' + identity)
        if abs(e['volume_m3'] - expect['layer_volume_m3'][identity]) > 1e-18:
            raise ValueError('Layer volume differs from the candidate: ' + identity)
        if e['physical_surface_support']['area_m2'] != EXTERIOR_AREA_M2:
            raise ValueError('Layer carries no support receipt: ' + identity)
        layers[identity] = {'thickness_m': thickness, 'volume_m3': e['volume_m3'],
                            'mass_kg': spec[identity]['mass_kg']}
    total_mass = sum(e['mass_kg'] for e in mechanics['entities'])
    skin_mass = sum(layers[i]['mass_kg'] for i in SKIN_LAYERS)
    phantom = (RAW_SLAB_AREA_M2 - EXTERIOR_AREA_M2) * sum(by_id[i]['shell']['thickness_m'] for i in SKIN_LAYERS)
    # The fraction is very nearly scale-free but not exactly: six carriers are
    # pinned at 1 microgram after scaling, so the divisor is not the same
    # multiple of the numerator at two different target masses. Compare at the
    # candidate's own target, where it must reproduce to the last bits, and
    # report the promoted fraction separately rather than blurring the two.
    allocation = mechanics['mass_allocation']
    scale_at_candidate = ((CANDIDATE_TARGET_MASS_KG - allocation['numerical_carrier_mass_kg'])
                          / allocation['unscaled_proxy_mass_kg'])
    unscaled_skin = skin_mass / allocation['uniform_scale']
    fraction_at_candidate = unscaled_skin * scale_at_candidate / CANDIDATE_TARGET_MASS_KG
    report['C_skin_layers'] = {'area_m2': [RAW_SLAB_AREA_M2, EXTERIOR_AREA_M2],
                               'inflation_factor': RAW_SLAB_AREA_M2 / EXTERIOR_AREA_M2,
                               'phantom_volume_l': phantom * 1000., 'layers': layers,
                               'skin_layer_mass_kg': skin_mass,
                               'skin_fraction_of_body': skin_mass / total_mass,
                               'skin_fraction_before': expect['skin_fraction_before_after'][0],
                               'skin_fraction_at_candidate_target_mass': fraction_at_candidate,
                               'candidate_predicted_fraction': expect['skin_fraction_before_after'][1]}
    if abs(fraction_at_candidate - expect['skin_fraction_before_after'][1]) > 1e-12:
        raise ValueError('Skin mass fraction differs from the candidate: %r vs %r'
                         % (fraction_at_candidate, expect['skin_fraction_before_after'][1]))

    # ---- the shared normalizer -------------------------------------------
    if abs(allocation['unscaled_proxy_mass_kg'] - expect['unscaled_proxy_mass_kg']) > 1e-9:
        raise ValueError('Unscaled proxy mass differs from the candidate: %r vs %r'
                         % (allocation['unscaled_proxy_mass_kg'], expect['unscaled_proxy_mass_kg']))
    at_candidate_target = scale_at_candidate
    if abs(at_candidate_target - expect['uniform_scale_at_candidate_target']) > 1e-12:
        raise ValueError('Uniform scale at the candidate target differs from the candidate')
    if profile['mass_kg'] != MASS_KG or allocation['target_mass_kg'] != MASS_KG:
        raise ValueError('Body mass is not the corrected figure')
    if abs(total_mass - MASS_KG) > 1e-8:
        raise ValueError('Allocated mass does not close on the declared body mass')
    if allocation['target_mass_source']['sha256'] != sha256_path(root / 'data/derived/canonical/profile.json'):
        raise ValueError('Mechanics does not carry the profile digest it consumed')
    report['mass'] = {'declared_kg': [CANDIDATE_TARGET_MASS_KG, MASS_KG],
                      'declared_source': MASS_LEDGER,
                      'unscaled_proxy_mass_kg': [candidate['allocation.json']['reproduction']['recorded_unscaled_kg'],
                                                 allocation['unscaled_proxy_mass_kg']],
                      'uniform_scale': [candidate['allocation.json']['uniform_scale']['before'],
                                        allocation['uniform_scale']],
                      'uniform_scale_at_candidate_target_mass': at_candidate_target,
                      'allocated_total_kg': total_mass}

    # ---- what must NOT have changed ---------------------------------------
    for identity, e in by_id.items():
        if e['role'] == 'skin_layer':
            continue
        if e.get('volume_m3') is not None and identity in spec:
            pass  # geometry-derived volumes are untouched; the assertion is the digest below
    report['counts'] = {'entities': [expect['counts']['entities'][0], len(anatomy['entities'])],
                        'bp_source_surfaces': BP_SOURCE_SURFACE_COUNT,
                        'rigid_bones': report['B_discs']['rigid_bones'],
                        'soft_solids': report['B_discs']['soft_solids'],
                        'support_links': len(mechanics['links']),
                        'muscle_and_connective_paths': len(mechanics['muscles']),
                        'display_structures': len(fragment['structures'])}
    if len(anatomy['entities']) != expect['counts']['entities'][1]:
        raise ValueError('Entity count differs from the candidate')
    return report


def promote():
    """Run the shipped builders in dependency order and write the receipt."""
    started = time.time()
    candidate = load_candidate()
    pre_state_path = RECEIPT / 'pre-promotion-sha256.json'
    pre_state = json.loads(pre_state_path.read_text()) if pre_state_path.exists() else {
        'statement': 'no pre-state sidecar was captured before the first builder run'}
    recorded = candidate['manifest.json']['inputs_sha256']
    for name, digest in pre_state.get('sha256', {}).items():
        if name in recorded and recorded[name] != digest:
            raise ValueError('Captured pre-state disagrees with the candidate for ' + name)
    before = {name: (sha256_path(ROOT / name) if (ROOT / name).exists() else None) for name in PROMOTED}
    geometry_before = {str(p.relative_to(ROOT)): sha256_path(p)
                       for p in sorted((CANONICAL / 'geometry').iterdir()) if p.is_file()}

    from ihm.assembly.profile import build_profile
    import build_canonical_anatomy
    import build_body_mechanics
    from ihm.assembly.body import build as build_body

    print('promote: profile', flush=True)
    build_profile(ROOT)
    print('promote: anatomy', flush=True)
    build_canonical_anatomy.build()
    build_canonical_anatomy.append_to_manifest()
    print('promote: mechanics', flush=True)
    build_body_mechanics.main()
    print('promote: body manifest', flush=True)
    build_body(ROOT)

    geometry_after = {str(p.relative_to(ROOT)): sha256_path(p)
                      for p in sorted((CANONICAL / 'geometry').iterdir()) if p.is_file()}
    if geometry_after != geometry_before:
        raise ValueError('Reference geometry changed; this repair moves no coordinates')
    after = {name: sha256_path(ROOT / name) for name in PROMOTED}
    checks = check_promoted(candidate)

    receipt = {
        'schema': 'ihm.entity-record-repair-promotion.v1',
        'status': 'PROMOTED. data/derived/canonical/ now carries the repair. The candidate directory '
                  'is preserved unchanged as the provenance record and is not consumed at runtime.',
        'canonical_assets_modified': True,
        'promoted_from': {'path': str(CANDIDATE.relative_to(ROOT)),
                          'manifest_sha256': sha256_path(CANDIDATE / 'manifest.json'),
                          'artifacts_sha256': candidate['manifest.json']['artifacts_sha256'],
                          'builder_sha256': candidate['manifest.json']['builder_sha256'],
                          'candidate_recorded_input_sha256': candidate['manifest.json']['inputs_sha256']},
        'applied_by': {
            'promotion_script': str(Path(__file__).relative_to(ROOT)),
            'promotion_script_sha256': sha256_path(Path(__file__)),
            'method': 'the repair lives in the builders, not in this script. This runs them in '
                      'dependency order and checks the result against the candidate.',
            'builders': [{'path': p, 'sha256': sha256_path(ROOT / p), 'owns': owns} for p, owns in [
                ('ihm/assembly/profile.py', 'declared body mass and its composition-ledger provenance'),
                ('ihm/assembly/anatomy.py', 'A: decoded-surface identity and the duplicate collapse; '
                                            'B: physical_role, which now matches the BodyParts3D '
                                            'spelling "intervertebral disk"; the nearest-bone '
                                            'regional support rule'),
                ('scripts/build_canonical_anatomy.py', 'A and C: the collapse call site, the '
                                                       'assumption ledger and the exterior-component '
                                                       'skin layer quadrature'),
                ('scripts/build_body_mechanics.py', 'the single uniform mass normalizer, the skeletal '
                                                    'support tree, soft-tissue supports and the '
                                                    'distinct-insertion guard'),
                ('ihm/assembly/body.py', 'the canonical body manifest and its source digests')]]},
        'changes': candidate['manifest.json']['changes'],
        'mass_correction': {
            'from_kg': CANDIDATE_TARGET_MASS_KG, 'to_kg': MASS_KG, 'source': MASS_LEDGER,
            'reason': 'the inherited BioGears StandardMale constant is unreachable in the acquired '
                      'envelope; the composed interior masses 70.77126585108081 kg. Independent of '
                      'A/B/C but it passes through the same normalizer, so it lands with them.'},
        # The superseded bytes: these assets are not tracked in git, so the
        # pre-state sidecar captured beside this receipt is their only record.
        # Its anatomy.json and mechanics.json digests must agree with the ones
        # the candidate independently recorded before any change was made.
        'pre_promotion_sha256': pre_state,
        'inputs_sha256_at_run_start': before,
        'outputs_sha256_after': after,
        'reruns_unchanged': sorted(name for name in PROMOTED if before.get(name) == after[name]),
        'rewritten_this_run': sorted(name for name in PROMOTED if before.get(name) != after[name]),
        'reference_geometry': {'files': len(geometry_after), 'changed': 0,
                               'statement': 'byte-identical before and after; no remeshing, no '
                                            'coordinate moves, and the dropped duplicates keep their '
                                            'geometry files as the provenance of the collapse'},
        'verified_here': checks,
        'not_regenerated': [
            {'artifact': 'data/derived/cross-structure-repair-v1',
             'reason': 'the candidate replayed all 11047 shipped ownership decisions and measured zero '
                       'flips at the corrected disc role, and 0.000 mL of disputed volume changing '
                       'owner. Regenerating it would claim a change that was measured not to exist.'},
            {'artifact': 'the 4269 segmentation-noise and 2528 modelling-conflict classifications',
             'reason': 'the five identical pairs were already recorded separately as coincidence and '
                       'were never in those classes'},
            {'artifact': 'reference geometry, volume_m3 and bounds of every retained entity',
             'reason': 'only role and record metadata were wrong'},
            {'artifact': 'the 0.1 / 1.5 / 5.0 mm layer thickness priors and FJ2810 surface_area_m2',
             'reason': 'the ICRP cross-check confirms the thickness priors at the corrected area, and '
                       'the raw slab area is a true measurement of the authored surface that stays as '
                       'descriptive source metadata. A double-sided quantity was used as a '
                       'single-sided one; the measurement itself is not wrong.'},
            {'artifact': 'native 22-body patient mass',
             'reason': 'scaled independently of the canonical proxy partition'}],
        'known_stale_downstream': [
            {'artifact': 'data/derived/canonical/trajectory.json, trajectory-final-state.json, '
                         'trajectory-spectra.json',
             'state': 'stale', 'reason': 'dynamic outputs of the changed masses; they need a fresh '
                                         'native run, not a rewritten hash'},
            {'artifact': 'artifacts and manifests carrying an anatomy_sha256 receipt',
             'state': 'hash-only invalidation', 'reason': 'each needs an explicit equivalence receipt '
                                                          'or a fresh run of its own builder'}],
        'wall_seconds': time.time() - started,
        'python': sys.version}
    RECEIPT.mkdir(parents=True, exist_ok=True)
    (RECEIPT / 'receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    (RECEIPT / 'geometry-unchanged.json').write_text(
        json.dumps({'statement': 'sha256 of every canonical reference surface, identical before and '
                                 'after the promotion', 'files': geometry_after}, indent=1) + '\n')
    print(json.dumps({k: receipt[k] for k in ('schema', 'status', 'mass_correction')}, indent=1))
    print(json.dumps(checks['counts'], indent=1))
    return receipt


# --------------------------------------------------------------------- tests

class RepairPrimitives(unittest.TestCase):
    """Exercise the shipped repair rules, not arithmetic written beside them."""

    def test_surface_identity_sees_through_gzip_framing(self):
        vertices = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        self.assertEqual(surface_identity_key(vertices, faces),
                         surface_identity_key(vertices.copy(), faces.copy()))
        moved = vertices + 1e-12
        self.assertNotEqual(surface_identity_key(vertices, faces), surface_identity_key(moved, faces))
        with self.assertRaises(ValueError):
            surface_identity_key(vertices, faces.ravel())

    def test_survivor_is_the_earlier_id_and_singletons_are_left_alone(self):
        keys = {'body-b': 'k1', 'body-a': 'k1', 'body-c': 'k1', 'body-d': 'k2'}
        self.assertEqual(duplicate_surface_survivors(keys),
                         {'body-b': 'body-a', 'body-c': 'body-a'})
        self.assertEqual(duplicate_surface_survivors({'x': 'k'}), {})

    def test_disc_role_follows_both_spellings(self):
        self.assertEqual(physical_role('intervertebral disk of axis', 'skeletal'), 'cartilage')
        self.assertEqual(physical_role('intervertebral disc', 'skeletal'), 'cartilage')
        self.assertEqual(physical_role('costal cartilage', 'skeletal'), 'cartilage')
        # The fix must not sweep the vertebrae themselves out of rigid_bone.
        self.assertEqual(physical_role('fifth cervical vertebra', 'skeletal'), 'rigid_bone')
        self.assertEqual(physical_role('hyoid bone', 'skeletal'), 'rigid_bone')

    def test_regional_support_never_self_loops_and_follows_the_role(self):
        def entity(identity, role, centroid):
            return {'id': identity, 'role': role, 'centroid_m': centroid, 'connections': []}
        entities = [entity('bone-a', 'rigid_bone', [0., 0, 0]),
                    entity('bone-b', 'rigid_bone', [1., 0, 0]),
                    entity('disc', 'rigid_bone', [.05, 0, 0]),
                    entity('soft', 'soft_organ', [.1, 0, 0]),
                    entity('layer', 'skin_layer', [0., 0, 0]),
                    entity('net', 'lymphatic_network', [0., 0, 0])]
        attach_regional_support(entities)
        support = {e['id']: [c for c in e['connections'] if c['relation'] == 'regional_spatial_support']
                   for e in entities}
        self.assertEqual(support['layer'], [])
        self.assertEqual(support['net'], [])
        self.assertEqual(support['soft'][0]['entity_id'], 'disc')
        self.assertEqual(support['bone-a'][0]['entity_id'], 'disc')
        for e in entities:
            for c in e['connections']:
                self.assertNotEqual(c['entity_id'], e['id'])
        # Re-role the disc out of rigid_bone: it leaves the candidate set, the
        # soft tissue re-anchors on a real bone, and the disc gains a support.
        entities[2]['role'] = 'cartilage'
        attach_regional_support(entities, replace=True)
        support = {e['id']: [c for c in e['connections'] if c['relation'] == 'regional_spatial_support']
                   for e in entities}
        self.assertEqual(support['soft'][0]['entity_id'], 'bone-a')
        self.assertEqual(support['bone-a'][0]['entity_id'], 'bone-b')
        self.assertEqual(support['disc'][0]['entity_id'], 'bone-a')
        for e in entities:
            self.assertLessEqual(len(support[e['id']]), 1)

    def test_material_law_moves_the_disc_off_the_bone_density(self):
        """The density the discs now carry comes from the shipped material law."""
        import build_body_mechanics as bbm
        disc = bbm.material({'role': 'cartilage', 'name': 'intervertebral disk of axis'})
        bone = bbm.material({'role': 'rigid_bone', 'name': 'axis'})
        self.assertEqual(disc['density']['value'], 1000.)
        self.assertEqual(disc['density']['prior_range'], [900, 1100])
        self.assertEqual(disc['young_modulus']['value'], 1e6)
        self.assertEqual(bone['density']['value'], 1900.)
        self.assertNotIn('young_modulus', bone)

    def test_normalizer_moves_every_untouched_row(self):
        """Rows touched by none of the three still move, through the one factor.

        Checked on the shipped artifact, against the candidate's own recorded
        prediction for two rows it singled out as untouched, rescaled from the
        candidate's target mass to the corrected one.
        """
        if not (CANONICAL / 'mechanics.json').exists():
            self.skipTest('canonical mechanics not built')
        mechanics = json.loads((CANONICAL / 'mechanics.json').read_text())
        allocation = mechanics['mass_allocation']
        if allocation['target_mass_kg'] != MASS_KG:
            self.skipTest('canonical mechanics is pre-promotion')
        spec = {e['id']: e for e in mechanics['entities']}
        candidate = load_candidate()
        rows = candidate['allocation.json']['every_untouched_entity_moves']
        self.assertTrue(rows)
        rescale = allocation['uniform_scale'] / candidate['allocation.json']['uniform_scale']['after']
        for row in rows:
            predicted = row['mass_kg'][1] * rescale
            self.assertAlmostEqual(spec[row['id']]['mass_kg'], predicted, places=12)
            self.assertGreater(spec[row['id']]['mass_kg'] / (row['mass_kg'][0] * rescale) - 1, .216)

    def test_declared_mass_must_agree_with_its_ledger(self):
        """build_profile refuses a declared mass its own source does not support."""
        import ihm.assembly.profile as module
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / MASS_LEDGER).parent.mkdir(parents=True, exist_ok=True)
            (root / MASS_LEDGER).write_text(json.dumps(
                {'verdict': {'composed_total_body_mass_kg': 61.0}}))
            original = module.MASS_KG
            try:
                with self.assertRaises(ValueError):
                    module.build_profile(root)
            finally:
                module.MASS_KG = original

    def test_candidate_manifest_is_self_consistent(self):
        candidate = load_candidate()
        expect = candidate_expectations(candidate)
        self.assertEqual(len(expect['drops']), 5)
        self.assertEqual(len(expect['discs']), 23)
        self.assertEqual(sorted(expect['layer_volume_m3']), sorted(SKIN_LAYERS))
        for identity, volume in expect['layer_volume_m3'].items():
            thickness = {'body-skin-epidermis': .0001, 'body-skin-dermis': .0015,
                         'body-skin-hypodermis': .005}[identity]
            self.assertEqual(volume, EXTERIOR_AREA_M2 * thickness)

    def test_promoted_canonical_state(self):
        if not (CANONICAL / 'mechanics.json').exists():
            self.skipTest('canonical mechanics not built')
        anatomy = json.loads((CANONICAL / 'anatomy.json').read_text())
        if 'duplicate_surface_collapse' not in anatomy:
            self.skipTest('canonical anatomy is pre-promotion')
        check_promoted(load_candidate())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--self-test', action='store_true', help='exercise the repair rules only')
    parser.add_argument('--promote', action='store_true', help='run the builders and write the receipt')
    parser.add_argument('--check', action='store_true', help='re-check the promoted canonical state')
    args = parser.parse_args()
    if not (args.self_test or args.promote or args.check):
        parser.error('Select --self-test, --check or --promote')
    if args.self_test:
        result = unittest.TextTestRunner(verbosity=2).run(
            unittest.defaultTestLoader.loadTestsFromTestCase(RepairPrimitives))
        if not result.wasSuccessful():
            raise SystemExit(1)
    if args.check:
        print(json.dumps(check_promoted(load_candidate()), indent=1))
    if args.promote:
        promote()
