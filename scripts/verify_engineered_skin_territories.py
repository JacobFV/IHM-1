#!/usr/bin/env python3
"""Bounded actual-source materialization and contact-accounting tests."""
import json
from pathlib import Path
import unittest

from materialize_engineered_skin_territories import ROOT, generate, three_region_proposal, reduce_contacts
from ihm.assembly.lymph_registration import SkinTerritoryInventory


class Checks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = generate()

    def test_reproducible_retained_inventory(self):
        retained = json.loads((ROOT/'data/research/engineered_skin_territories/materialization.json').read_text())
        self.assertEqual(self.result, retained)
        self.assertEqual(generate(), retained)
        ids=[i for r in retained['regions'] for i in r['triangle_ids']]
        self.assertEqual(len(ids),len(set(ids)))
        self.assertEqual(len(ids),7448)
        self.assertEqual(sum(r['triangle_count'] for r in retained['inventory']['regions']),203382)
        self.assertAlmostEqual(sum(r['area_fraction'] for r in retained['inventory']['regions']),1.)
        self.assertFalse(retained['inventory']['anatomical_drainage_registration_validated'])
        self.assertIsNone(retained['native_configuration'])
        diagnostic=retained['surface_diagnostic']
        self.assertEqual(len(diagnostic['face_component_ids']),203382)
        self.assertEqual(diagnostic['components'][0]['boundary_edges'],827)
        self.assertEqual(diagnostic['inconsistent_interior_edge_winding'],0)
        self.assertTrue(all(diagnostic['face_component_ids'][i]==0 for i in ids))
        self.assertEqual(len(retained['contact_eligible_triangle_ids']),109183)

    def test_explicit_three_union_proposal_is_conservative_and_unbound(self):
        left=[r['id'] for r in self.result['regions'] if r['id'].startswith('left_')]
        right=[r['id'] for r in self.result['regions'] if r['id'].startswith('right_')]
        proposal=three_region_proposal(self.result,left,right)
        self.assertEqual(sum(proposal['initial_area_weights']),1.)
        self.assertTrue(all(x>0 for x in proposal['initial_area_weights']))
        self.assertNotEqual(proposal['initial_area_weights'],[.2,.3,.5])
        self.assertIsNone(proposal['installed_native_binding'])
        for a,b in [(left,left),([],right),(['unknown'],right),(['unresolved'],right)]:
            with self.assertRaises(ValueError):three_region_proposal(self.result,a,b)

    def test_source_contact_maps_to_exact_material_owner(self):
        source=self.result['source_receipts'][0]
        model=SkinTerritoryInventory(ROOT/source['path'],source['sha256'],self.result['regions'])
        region=self.result['regions'][0]; triangle=region['triangle_ids'][0]
        area=float(model.area[triangle])*.5; force=1000*area
        contact=reduce_contacts(self.result,[{'triangle_id':triangle,'area_m2':area,
                                       'normal_pressure_pa':1000.,'force_n':[0.,force,0.]}])
        row=contact['regions'][region['id']]
        self.assertEqual(row['normal_force_n'],force)
        self.assertEqual(row['force_n'],[0.,force,0.])
        self.assertEqual(sum(r['normal_force_n'] for r in contact['regions'].values()),force)
        self.assertEqual(contact['native_commands'],[])
        inner=next(i for i,c in enumerate(self.result['surface_diagnostic']['face_component_ids']) if c==1)
        with self.assertRaises(ValueError):
            reduce_contacts(self.result,[{'triangle_id':inner,'area_m2':1e-8,'normal_pressure_pa':1.,'force_n':[0.,1e-8,0.]}])

    def test_parameters_are_bounded_and_change_mask_honestly(self):
        for parameters in [{'end_quantile':True,'boundary_band_m':.003,'maximum_radius_m':.13},
                           {'end_quantile':.05,'boundary_band_m':-.1,'maximum_radius_m':.13}]:
            with self.assertRaises(ValueError):generate(parameters=parameters)
        wider=generate(parameters={'end_quantile':.05,'boundary_band_m':.006,'maximum_radius_m':.13})
        self.assertGreater(wider['inventory']['regions'][0]['triangle_count'],
                           self.result['inventory']['regions'][0]['triangle_count'])


if __name__=='__main__':unittest.main()
