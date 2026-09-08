#!/usr/bin/env python3
"""Promote the sourced tissue-material assignment into the canonical body.

The candidate at `data/derived/tissue-material-assignment-candidate-v1` measured
two defects in the canonical mass ledger and proved them without touching
anything:

  A  every non-bone entity carried a density of 1000 kg/m3 and every bone 1900.
     Neither is a tissue density -- 1000 is water and 1900 is marrow-free
     cortical bone, which cannot describe a surface that encloses the marrow
     cavity. The 4000 entities now draw from a sourced per-tissue table.

  B  a surface's signed integral was taken as its tissue volume whether or not
     the surface was closed. For an anatomical SHEET authored as one open
     surface -- 92 fascias, intermuscular septa, joint capsules and serosae --
     the signed integral is the volume of the structures the sheet WRAPS. The
     84 open connective-tissue sheets alone carried 29.6 L, a third of the whole
     unscaled proxy mass, and it was the muscle they wrap counted twice. A
     hollow viscus is one closed surface around wall AND lumen, so 1.92 L of gut,
     airway and urinary lumen was carried at tissue density.

Neither repair is applied by this script. Both live in the shipped module that
owns the field, `ihm/assembly/tissue_materials.py`, called from
`scripts/build_body_mechanics.py`, so a rebuild reaches the same answer and a
future rebuild cannot silently revert them. This runs the builders in dependency
order, checks the result against the candidate's own recorded predictions, and
writes a receipt with input hashes, output hashes and this script's own digest.

    scripts/promote_tissue_material_assignment.py --self-test
    scripts/promote_tissue_material_assignment.py --promote
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from ihm.assembly import tissue_materials as tm
from ihm.assembly.profile import MASS_KG, MASS_LEDGER

CANDIDATE = ROOT / 'data/derived/tissue-material-assignment-candidate-v1'
RECEIPT = ROOT / 'data/derived/tissue-material-assignment-promotion-v1'
CANONICAL = ROOT / 'data/derived/canonical'

PROMOTED = ['data/derived/canonical/mechanics.json',
            'data/derived/canonical/body.json']

# What the assignment supersedes. Recorded so the defect is not lost when the
# artifact that carried it is rewritten.
SUPERSEDED = {
    'density_non_bone_kg_m3': 1000.,
    'density_rigid_bone_kg_m3': 1900.,
    'unscaled_proxy_mass_kg': 93.47748928171038,
    'uniform_scale': 0.7570923710490256,
}


def sha256_path(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load_candidate():
    manifest = json.loads((CANDIDATE / 'manifest.json').read_text())
    mismatched = {name: sha256_path(CANDIDATE / name)
                  for name, digest in manifest['artifacts_sha256'].items()
                  if sha256_path(CANDIDATE / name) != digest}
    if mismatched:
        raise ValueError('Candidate artifacts disagree with their manifest: %r' % mismatched)
    if manifest['rule_owner_sha256'] != sha256_path(ROOT / 'ihm/assembly/tissue_materials.py'):
        raise ValueError('The shipped rule changed since the candidate was measured; rebuild the '
                         'candidate before promoting so the prediction is against this rule')
    return {'manifest.json': manifest,
            'assignment.json': json.loads((CANDIDATE / 'assignment.json').read_text()),
            'colocation.json': json.loads((CANDIDATE / 'colocation.json').read_text())}


def check_promoted(candidate, root=ROOT):
    """Every claim the candidate made, checked against the promoted assets."""
    mech = json.loads((root / 'data/derived/canonical/mechanics.json').read_text())
    profile = json.loads((root / 'data/derived/canonical/profile.json').read_text())
    allocation = mech['mass_allocation']
    predicted = candidate['assignment.json']['allocation']['assigned']
    checks = {}

    checks['proxy_mass_matches_the_candidate_prediction'] = abs(
        allocation['unscaled_proxy_mass_kg'] - predicted['unscaled_proxy_mass_kg']) < 1e-6
    checks['uniform_scale_is_above_one'] = allocation['uniform_scale'] > 1.0
    checks['uniform_scale_left_the_superseded_value'] = abs(
        allocation['uniform_scale'] - SUPERSEDED['uniform_scale']) > 0.1
    checks['entity_masses_still_sum_to_the_declared_mass'] = abs(
        sum(e['mass_kg'] for e in mech['entities']) - float(profile['mass_kg'])) < 1e-8
    checks['declared_mass_still_agrees_with_its_ledger'] = abs(
        json.loads((root / MASS_LEDGER).read_text())['verdict']['composed_total_body_mass_kg']
        - MASS_KG) < 5e-5

    # No entity may still carry water or marrow-free cortical bone as its tissue
    # density, and every density must name a source and a tier.
    unsourced, superseded_rows = [], []
    for e in mech['entities']:
        row = (e.get('material') or {}).get('density')
        if not row:
            unsourced.append(e['id'])
            continue
        if not row.get('sources') or not row.get('tier') or not row.get('tissue_class'):
            unsourced.append(e['id'])
        if row['value'] in (SUPERSEDED['density_non_bone_kg_m3'],
                            SUPERSEDED['density_rigid_bone_kg_m3']):
            superseded_rows.append(e['id'])
    checks['every_density_names_a_source_a_tier_and_a_tissue_class'] = not unsourced
    checks['no_entity_kept_a_superseded_density'] = not superseded_rows

    # The sheets actually lost their enclosures, and nothing gained volume.
    grew = [e['id'] for e in mech['entities']
            if (e.get('volume_representation') or {}).get('open_surface_signed_integral_m3')
            and e['volume_m3'] > (e['volume_representation']['open_surface_signed_integral_m3']
                                  + 1e-15)]
    checks['no_representation_rule_increased_a_volume'] = not grew
    sheets = [e for e in mech['entities']
              if (e.get('volume_representation') or {}).get('rule') == 'open_membrane_sheet']
    checks['open_sheets_are_reclassified'] = len(sheets) >= 80
    checks['fascia_lata_is_a_sheet_not_a_thigh'] = all(
        e['volume_m3'] < 0.0005 for e in mech['entities'] if 'fascia lata' in e['name'].lower())

    # The one adipose material in the model, which is what lets the fat budget
    # locate the subcutaneous depot at all.
    hypodermis = next((e for e in mech['entities'] if e['id'] == 'body-skin-hypodermis'), None)
    checks['hypodermis_carries_an_adipose_density'] = bool(
        hypodermis and hypodermis['material']['density']['value'] == tm.DENSITY['adipose'][0])

    detail = {'entities_without_a_sourced_density': unsourced[:20],
              'entities_still_at_a_superseded_density': superseded_rows[:20],
              'entities_whose_volume_grew': grew[:20],
              'open_membrane_sheets': len(sheets)}
    return checks, detail


def promote():
    started = time.time()
    candidate = load_candidate()
    before = {name: (sha256_path(ROOT / name) if (ROOT / name).exists() else None)
              for name in PROMOTED}
    geometry_before = {str(p.relative_to(ROOT)): sha256_path(p)
                       for p in sorted((CANONICAL / 'geometry').iterdir()) if p.is_file()}

    import build_body_mechanics
    from ihm.assembly.body import build as build_body

    print('promote: mechanics', flush=True)
    build_body_mechanics.main()
    print('promote: body manifest', flush=True)
    build_body(ROOT)

    geometry_after = {str(p.relative_to(ROOT)): sha256_path(p)
                      for p in sorted((CANONICAL / 'geometry').iterdir()) if p.is_file()}
    if geometry_after != geometry_before:
        raise ValueError('Reference geometry changed; this assignment moves no coordinates')
    after = {name: sha256_path(ROOT / name) for name in PROMOTED}
    checks, detail = check_promoted(candidate)
    mech = json.loads((CANONICAL / 'mechanics.json').read_text())

    receipt = {
        'schema': 'ihm.tissue-material-assignment-promotion.v1',
        'status': ('PROMOTED. data/derived/canonical/ now carries the sourced tissue densities and '
                   'the representation rule. The candidate directory is preserved unchanged as the '
                   'provenance record and is not consumed at runtime.'),
        'canonical_assets_modified': True,
        'promoted_from': {
            'path': str(CANDIDATE.relative_to(ROOT)),
            'manifest_sha256': sha256_path(CANDIDATE / 'manifest.json'),
            'artifacts_sha256': candidate['manifest.json']['artifacts_sha256'],
            'candidate_recorded_input_sha256': candidate['manifest.json']['inputs_sha256']},
        'applied_by': {
            'promotion_script': str(Path(__file__).relative_to(ROOT)),
            'promotion_script_sha256': sha256_path(Path(__file__)),
            'method': 'the assignment lives in the builders, not in this script. This runs them in '
                      'dependency order and checks the result against the candidate.',
            'builders': [{'path': p, 'sha256': sha256_path(ROOT / p), 'owns': owns} for p, owns in [
                ('ihm/assembly/tissue_materials.py',
                 'the sourced density table, the membrane and hollow-viscus representation rules, '
                 'and the tier on every one of them'),
                ('scripts/build_body_mechanics.py',
                 'the single uniform mass normalizer; calls the module for every entity volume '
                 'and every entity density'),
            ]]},
        'supersedes': SUPERSEDED,
        'mass_allocation_after': {k: mech['mass_allocation'][k] for k in
                                  ('unscaled_proxy_mass_kg', 'uniform_scale', 'target_mass_kg')},
        'material_assignment': mech['mass_allocation']['material_assignment'],
        'decomposition_of_the_superseded_excess': {
            k: v for k, v in candidate['assignment.json']['decomposition'].items()
            if not isinstance(v, dict) or k == 'volume_double_count_by_colocation'},
        'fat_budget': candidate['assignment.json']['fat_budget'],
        'colocation_resolution': {
            'considered': candidate['colocation.json']['considered'],
            'resolved_as_duplicate': candidate['colocation.json']['resolved_as_duplicate'],
            'unresolved': candidate['colocation.json']['unresolved'],
            'action_taken': candidate['colocation.json']['action_taken']},
        'checks': checks,
        'check_detail': detail,
        'inputs_sha256_at_run_start': before,
        'outputs_sha256_after': after,
        'reference_geometry_unchanged': True,
        'invalidated': [
            {'artifact': 'data/derived/interstitial-matrix-v1/composition.json',
             'state': 'superseded by data/derived/interstitial-matrix-composition-v2/'
                      'composition.json, rebuilt against the promoted mechanics in the same '
                      'session with the composition stage only. The v1 directory is left intact '
                      'because its geometry, lumen and provenance stages are unaffected and cost '
                      '25 minutes to reproduce; only the composition stage reads mechanics.json.'},
            {'artifact': 'data/derived/tissue-material-candidate-v1',
             'state': 'rebuilt in the same session. Its `audit` block censuses the canonical '
                      'material fields and had to be re-run to stop describing the superseded '
                      '1000/1900 pair; two of the conflicts it raised (soft_organ.lung.'
                      'density_effective and rigid_bone.density_effective) are resolved by this '
                      'assignment and its self-test now asserts that they are. Its `geometry_mass` '
                      'block is byte-identical before and after, because that half weighs the '
                      'material-domains voxel partition by SYSTEM_DENSITY and never reads a '
                      'canonical density -- which is why the profile mass gate at '
                      'interstitial-composition-prior-v1/ledger.json did not move, and it was '
                      're-run to confirm it still composes 70.77126585 kg.'},
            {'artifact': 'artifacts carrying a mechanics_sha256 receipt',
             'state': 'hash-only invalidation'},
        ],
        'not_established': candidate['assignment.json']['not_established'],
        'wall_seconds': time.time() - started,
    }
    RECEIPT.mkdir(parents=True, exist_ok=True)
    (RECEIPT / 'receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    if not all(checks.values()):
        raise ValueError('Promotion checks failed: %r' % {k: v for k, v in checks.items() if not v})
    print(json.dumps({k: receipt[k] for k in ('schema', 'status', 'mass_allocation_after', 'checks')},
                     indent=1))
    return receipt


class Tests(unittest.TestCase):
    def test_the_rule_only_ever_reduces_a_volume(self):
        for name, role, prior, area, wt in (
                ('left fascia lata', 'connective_tissue', 0.004086, 0.16932, False),
                ('stomach', 'soft_organ', 0.0005665, 0.0383, False),
                ('left biceps brachii', 'muscle', 0.0004, 0.02, False),
                ('pleura', 'soft_organ', 0.0002147, 0.8909, True)):
            value, _ = tm.material_volume(name, role, prior, area, wt)
            self.assertLessEqual(value, prior + 1e-15, name)

    def test_no_density_row_is_the_superseded_water_or_cortical_value(self):
        for key, (value, source, tier, note) in tm.DENSITY.items():
            self.assertNotEqual(value, SUPERSEDED['density_non_bone_kg_m3'], key)
            self.assertNotEqual(value, SUPERSEDED['density_rigid_bone_kg_m3'], key)

    def test_the_normalizer_still_lands_on_the_declared_mass(self):
        mech = json.loads((CANONICAL / 'mechanics.json').read_text())
        profile = json.loads((CANONICAL / 'profile.json').read_text())
        self.assertAlmostEqual(sum(e['mass_kg'] for e in mech['entities']),
                               float(profile['mass_kg']), places=8)

    def test_the_profile_mass_gate_is_untouched_by_this_lane(self):
        ledger = json.loads((ROOT / MASS_LEDGER).read_text())
        self.assertLess(abs(ledger['verdict']['composed_total_body_mass_kg'] - MASS_KG), 5e-5)

    def test_promoted_canonical_state(self):
        checks, detail = check_promoted(load_candidate())
        self.assertTrue(all(checks.values()), {k: v for k, v in checks.items() if not v})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--promote', action='store_true', help='run the builders and write the receipt')
    args = parser.parse_args()
    if args.promote:
        promote()
        return
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.TestLoader().loadTestsFromTestCase(Tests))
    raise SystemExit(0 if result.wasSuccessful() else 1)


if __name__ == '__main__':
    main()
