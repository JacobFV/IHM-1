#!/usr/bin/env python3
"""Bounded source-only inventory tests, including retained canonical skin bytes."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.lymph_registration import SkinTerritoryInventory
from audit_lymph_registration import audit


class Checks(unittest.TestCase):
    def fixture(self, directory):
        p = Path(directory) / 'skin.json'
        p.write_text(json.dumps({'positions': [0,0,0, 1,0,0, 0,1,0, 1,1,0],
                                 'indices': [0,1,2, 1,3,2]}))
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        row = {'id': 'candidate', 'triangle_ids': [0], 'geometry_sha256': sha,
               'annotation_receipt': {'kind': 'synthetic_test_only'}}
        return p, sha, row

    def test_actual_canonical_inventory_and_byte_receipts(self):
        result = audit()
        self.assertEqual(result['skin']['inventory']['triangle_count'], 203382)
        self.assertEqual(result['skin']['inventory']['regions'][0]['area_fraction'], 1.)
        self.assertEqual(result['graph']['named_nodes'], 0)
        self.assertEqual(result['graph']['validated_drainage_associations'], 0)
        self.assertEqual(len(result['lymph_node_groups']), 158)
        self.assertFalse(result['native_engineering_regions_rebound'])
        template = json.loads((ROOT / 'data/research/lymph_registration/partition_template.json').read_text())
        self.assertFalse(template['activation_allowed'])
        self.assertEqual(len(template['territories']), 8)
        self.assertTrue(all(t['triangle_ids'] is None and t['native_owner'] is None
                            and t['first_destination']['accepted_entity_id'] is None
                            for t in template['territories']))
        for r in json.loads((ROOT / 'data/research/lymph_registration/retrieval.json').read_text()):
            if 'path' in r:
                raw = (ROOT / r['path']).read_bytes()
                self.assertEqual(len(raw), r['bytes'])
                self.assertEqual(hashlib.sha256(raw).hexdigest(), r['sha256'])

    def test_exclusive_inventory_residual_and_force_conservation(self):
        with tempfile.TemporaryDirectory() as directory:
            p, sha, row = self.fixture(directory)
            model = SkinTerritoryInventory(p, sha, [row])
            self.assertEqual(sum(r['triangle_count'] for r in model.report()['regions']), 2)
            self.assertEqual(sum(r['area_fraction'] for r in model.report()['regions']), 1.)
            result = model.reduce_contacts([
                {'triangle_id': 0, 'area_m2': .1, 'normal_pressure_pa': 100., 'force_n': [0,10,0]},
                {'triangle_id': 1, 'area_m2': .2, 'normal_pressure_pa': 50., 'force_n': [3,10,0]}])
            self.assertEqual(result['regions']['candidate']['whole_territory_mean_pressure_pa'], 20.)
            self.assertEqual(sum(r['normal_force_n'] for r in result['regions'].values()), 20.)
            self.assertEqual(sum(r['force_n'][0] for r in result['regions'].values()), 3.)
            self.assertEqual(result['native_commands'], [])

    def test_duplicate_masks_invalid_indices_and_native_rebinding_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            p, sha, row = self.fixture(directory)
            for change in [{'triangle_ids': [0,0]}, {'triangle_ids': [True]},
                           {'triangle_ids': [2]}, {'geometry_sha256': 'changed'},
                           {'native_owner': 'Skin.region_a.extracellular'}, {'annotation_receipt': None}]:
                bad = {**row, **change}
                with self.subTest(change=change), self.assertRaises(ValueError):
                    SkinTerritoryInventory(p, sha, [bad])
            with self.assertRaises(ValueError):
                SkinTerritoryInventory(p, sha, [row, {**row, 'id': 'overlap'}])
            with self.assertRaises(ValueError):
                SkinTerritoryInventory(p, 'changed')

    def test_contact_overlap_and_excess_area_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            p, sha, row = self.fixture(directory)
            model = SkinTerritoryInventory(p, sha, [row])
            contact = {'triangle_id': 0, 'area_m2': .1, 'normal_pressure_pa': 100., 'force_n': [0,10,0]}
            for contacts in [[contact, contact], [{**contact, 'area_m2': .6}],
                             [{**contact, 'normal_pressure_pa': -1}],
                             [{**contact, 'force_n': [float('nan'),0,0]}]]:
                with self.assertRaises(ValueError):
                    model.reduce_contacts(contacts)


if __name__ == '__main__':
    unittest.main()
