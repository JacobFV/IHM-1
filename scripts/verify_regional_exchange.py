#!/usr/bin/env python3
"""Small, read-only checks against retained native Python adapter snapshots."""
from copy import deepcopy
import json
import hashlib
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.regional_exchange import PREFIX, REGIONS, observe_regional_skin

FIXTURE = ROOT / 'data/derived/audits/regional-python-vvdedt24'
HASHES = {
    'initial': 'cb37d4e4f5460d48e2938a374ab5a92148890559f7070482926732b2903e0b4d',
    'loaded': 'fea8f66846d3a12f1d7cea48cecd1a6ff935f3220248c6528887b5d4efccd922',
    'released': 'a80f276241564280e431cde0436ba71a8794683e08a17486740f734f38e80337',
}


def snapshot(name='initial'):
    raw = (FIXTURE / (name + '.json')).read_bytes()
    if hashlib.sha256(raw).hexdigest() != HASHES[name]:
        raise ValueError('Changed retained native fixture: ' + name)
    return json.loads(raw)


class Checks(unittest.TestCase):
    def test_actual_native_lifecycle_and_ownership(self):
        for name, steps in [('initial', 0), ('loaded', 1), ('released', 2)]:
            source = snapshot(name)
            before = deepcopy(source)
            result = observe_regional_skin(source, native_identity={'fixture': name},
                                           source_hashes={'receipt': 'test'})
            self.assertEqual(source, before)
            self.assertTrue(result['available'])
            self.assertEqual(result['completed_native_steps'], steps)
            self.assertEqual(len(result['native_compartments']), 3)
            self.assertEqual(len(result['fluid_transfers']), 12)
            self.assertEqual(len(result['native_path_observations']), 27)
            self.assertEqual(result['native_identity'], {'fixture': name})
            parent = result['aggregate_views']['Skin.extracellular']
            self.assertFalse(parent['accounting_owner'])
            self.assertFalse(parent['independent_store'])
            self.assertNotIn('Skin.extracellular', result['native_compartments'])
            self.assertAlmostEqual(math.fsum(c['volume_ml'] for c in result['native_compartments'].values()), parent['volume_ml'])
            for sub, mass in parent['mass_g'].items():
                self.assertAlmostEqual(math.fsum(c['mass_g'][sub] for c in result['native_compartments'].values()), mass)
            self.assertTrue(all(c['status'] in ('passed', 'unobserved') for c in result['ownership_checks']))
            json.dumps(result, allow_nan=False)

    def test_units_sources_unknown_concentrations_and_no_allocated_current_state(self):
        source = snapshot('loaded')
        result = observe_regional_skin(source)
        child = result['native_compartments']['Skin.region_a.extracellular']
        key = PREFIX + 'region_a.species.Albumin.mass_ug'
        self.assertEqual(child['mass_g']['Albumin'], source['values'][key] * 1e-6)
        self.assertEqual(child['quantity_sources']['mass_g']['Albumin']['source_key'], key)
        self.assertTrue(all(c is None for c in child['concentration_g_per_l'].values()))
        self.assertTrue(all(c is None for c in child['ionic_molarity_mmol_per_l'].values()))
        self.assertFalse(child['independent_chemical_law'])
        self.assertEqual(result['solute_fluxes'], [])
        self.assertNotEqual(child['volume_ml'], .2 * result['aggregate_views']['Skin.extracellular']['volume_ml'])
        del source['values'][key]
        unknown = observe_regional_skin(source)
        self.assertIsNone(unknown['native_compartments']['Skin.region_a.extracellular']['mass_g']['Albumin'])
        check = next(c for c in unknown['ownership_checks'] if c['check'] == 'children_minus_aggregate_mass.Albumin')
        self.assertEqual(check['status'], 'unobserved')
        source['values'] = {k: v for k, v in source['values'].items()
                            if not (k.startswith(PREFIX) and '.species.' in k)}
        missing = observe_regional_skin(source)
        self.assertIsNone(missing['native_compartments']['Skin.region_a.extracellular']['mass_g']['Albumin'])
        self.assertIsNone(missing['aggregate_views']['Skin.extracellular']['mass_g']['Albumin'])

    def test_actual_signed_incidence_and_no_serial_double_count(self):
        source = snapshot('released')
        result = observe_regional_skin(source)
        for region in REGIONS:
            def q(path):
                return source['values'][PREFIX + region + '.path.' + path + '.flow_ml_per_s']
            expected = q('SkinE2ToSkinE3') - q('SkinE3ToSkinI') - q('SkinE3ToSkinL1') - q('SkinSweating')
            self.assertAlmostEqual(result['internal_volume_rate_ml_per_s']['Skin.' + region + '.extracellular'], expected)
        self.assertAlmostEqual(math.fsum(result['internal_volume_rate_ml_per_s'].values()) + result['external_volume_rate_ml_per_s']['external.sweat'], 0.)
        self.assertFalse(any('ToGround' in t['native_path'] or 'Valve' in t['native_path'] for t in result['fluid_transfers']))
        self.assertTrue(all(t['native_path'].startswith('IHM_') for t in result['fluid_transfers']))

    def test_missing_installation_and_invalid_partial_or_extra_owner(self):
        self.assertFalse(observe_regional_skin({'values': {}})['available'])
        for suffix, value in [('region_a.volume_ml', None), ('region_a.volume_ml', -1.),
                              ('region_b.fraction', .4), ('fourth.volume_ml', 1.),
                              ('region_a.path.SkinE2ToSkinE3.flow_ml_per_s', None),
                              ('completed_native_steps', .5), ('region_a.pressure_mmhg', True)]:
            source = snapshot()
            source['values'][PREFIX + suffix] = value
            with self.subTest(suffix=suffix, value=value), self.assertRaises(ValueError):
                observe_regional_skin(source)

    def test_corrupted_volume_species_and_native_residual_fail_closed(self):
        for suffix in ['region_a.volume_ml', 'region_b.species.Albumin.mass_ug',
                       'aggregate.volume_ownership_residual_ml',
                       'aggregate.species.Sodium.ownership_residual_ug',
                       'region_a.fluid_step_residual_ml']:
            source = snapshot('loaded')
            source['values'][PREFIX + suffix] += 1.
            with self.subTest(suffix=suffix), self.assertRaises(ValueError):
                observe_regional_skin(source)
        source = snapshot()
        source['values']['tissue.Skin.extracellular.volume_ml'] += 1.
        with self.assertRaises(ValueError):
            observe_regional_skin(source)


if __name__ == '__main__':
    unittest.main()
