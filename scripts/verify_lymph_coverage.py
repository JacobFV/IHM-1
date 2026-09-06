#!/usr/bin/env python3
"""Bounded lymph coverage checks using the existing tissue snapshot fixture."""
from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.lymph_coverage import NativeRegionalLymphCoverage
from verify_body_exchange import snapshot


def anatomy():
    return {'entities': [
        {'id': 'skin-a', 'name': 'skin a', 'role': 'skin', 'surface_area_m2': 1.,
         'reference_geometry': {'path': 'a.json', 'sha256': 'fixture-a'}},
        {'id': 'skin-b', 'name': 'skin b', 'role': 'skin', 'surface_area_m2': 3.,
         'reference_geometry': {'path': 'b.json', 'sha256': 'fixture-b'}},
    ]}


class Checks(unittest.TestCase):
    def model(self, **kwargs):
        return NativeRegionalLymphCoverage(anatomy(), organs=('Skin',),
                                           native_identity={'run': 'fixture'}, **kwargs)

    def test_conservation_and_no_second_lymph_store(self):
        source = snapshot()
        before = deepcopy(source)
        result = self.model().observe(source)
        rows = result['regional_views']
        self.assertEqual([r['fraction'] for r in rows], [.25, .75])
        for field, total in [('volume_ml', 100.), ('albumin_mass_g', .1)]:
            self.assertEqual(math.fsum(r['quantities'][field]['value'] for r in rows), total)
        self.assertTrue(all(not r['accounting_owner'] and not r['independent_store'] for r in rows))
        self.assertEqual([r['native_owner'] for r in result['native_owners']], ['Skin.extracellular', 'Lymph'])
        self.assertTrue(all(r['native_owner'] != 'Lymph' for r in rows))
        self.assertEqual(source, before)
        self.assertFalse(result['whole_body_mass_closure_claimed'])

    def test_identity_units_and_unknowns(self):
        source = snapshot()
        del source['values']['tissue.Skin.extracellular.Albumin.mass_g']
        source['values']['tissue.Lymph.pressure_mmhg'] = float('nan')
        result = self.model(source_hashes={'fixture': 'abc'}).observe(source)
        self.assertEqual(result['native_identity'], {'run': 'fixture'})
        self.assertEqual(result['source_hashes'], {'fixture': 'abc'})
        for row in result['regional_views']:
            q = row['quantities']['albumin_mass_g']
            self.assertIsNone(q['value'])
            self.assertEqual(q['unit'], 'g')
            self.assertEqual(q['source_key'], 'tissue.Skin.extracellular.Albumin.mass_g')
            self.assertIsNone(row['regional_protein_flux_g_per_s'])
        self.assertIn('tissue.Lymph.pressure_mmhg', result['unobserved_source_keys'])
        json.dumps(result, allow_nan=False)

    def test_unlocalized_and_absent_topology(self):
        result = NativeRegionalLymphCoverage({'entities': []}, organs=('Skin', 'Fat')).observe(snapshot())
        self.assertEqual(set(result['unlocalized_native_owners']), {'Lymph', 'Skin.extracellular', 'Fat.extracellular'})
        self.assertEqual(result['topology']['status'], 'absent_source_graph')
        self.assertTrue(all(r['classification'] == 'unlocalized_native_view' for r in result['regional_views']))
        fat = next(r for r in result['regional_views'] if r['native_owner'] == 'Fat.extracellular')
        self.assertIsNone(fat['quantities']['volume_ml']['value'])

    def test_source_graph_is_not_native_topology(self):
        result = self.model(topology={'model_id': 'fixture', 'nodes': [1, 2], 'edges': [1], 'directed': False}).observe(snapshot())
        self.assertEqual(result['topology']['nodes'], 2)
        self.assertEqual(result['topology']['native_regional_mapping'], 'absent')
        self.assertFalse(result['topology']['allocated_lymph_fluid_or_protein'])
        self.assertEqual(len(result['pathways']), 4)
        self.assertIsNone(result['pathways'][1]['flow']['value'])
        self.assertAlmostEqual(result['pathways'][0]['flow']['value'], .2)

    def test_invalid_owner_duplicate_support_and_quantity(self):
        for organs in [('Skin', 'Skin'), ('Invented',)]:
            with self.assertRaises(ValueError):
                NativeRegionalLymphCoverage(anatomy(), organs=organs)
        source = anatomy()
        source['entities'].append(deepcopy(source['entities'][0]))
        with self.assertRaises(ValueError):
            NativeRegionalLymphCoverage(source)
        for invalid in [-1., float('inf'), True, '1']:
            source = snapshot()
            source['values']['tissue.Lymph.volume_ml'] = invalid
            with self.assertRaises(ValueError):
                self.model().observe(source)

    def test_workspace_geometry_receipt_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ['data/derived/canonical', 'scripts', 'ihm/assembly']:
                (root / name).mkdir(parents=True)
            for name in ['scripts/native_tissue_ports.h', 'ihm/assembly/body_microstructure.py', 'ihm/assembly/lymph_coverage.py']:
                (root / name).write_text('fixture')
            (root / 'data/derived/canonical/anatomy.json').write_text(json.dumps(anatomy()))
            (root / 'a.json').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'Changed anatomical'):
                NativeRegionalLymphCoverage.from_workspace(root, {'run': 'fixture'}, organs=('Skin',))


if __name__ == '__main__':
    unittest.main()
