#!/usr/bin/env python3
"""Bounded actual-source materialization and contact-accounting tests."""
import hashlib
import json
from pathlib import Path
import unittest

from materialize_engineered_skin_territories import ROOT, generate, three_region_proposal, reduce_contacts
from ihm.assembly.lymph_registration import SkinTerritoryInventory


class Checks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = generate()

    # data/research/engineered_skin_territories/materialization.json is a frozen source
    # identity, not a regenerable derived artifact. Its sha256 is pinned by
    # ihm/assembly/regional_skin_configuration.py (MATERIALIZATION_SHA) and
    # ihm/assembly/skin_microvascular_patch.py (TERRITORY_SHA), the configuration it
    # produces is compiled into data/research/configured_regional_skin/prepared_v1/
    # native_configured_regional_skin.h as configuration_sha256, and the built
    # whole_body_integrity_regional_skin_graph_v2 library pins that header. It therefore
    # records anatomy_sha256 -- the hash of the canonical anatomy CONTAINER it was cut
    # from -- and rebuilding that container must not reopen the frozen artifact.
    #
    # The invariant this asserts is consequently not "the retained bytes are byte-identical
    # to a fresh generate()" but the two claims that actually govern the artifact:
    #   1. everything generate() derives is still reproduced exactly, and
    #   2. every byte the artifact consumes is still pinned by the LIVE anatomy and still
    #      hashes to the receipt it froze.
    # (2) is a stronger check than the old whole-document equality, which only compared the
    # container hash and never re-read the geometry the container points at.
    ANATOMY_PATH = 'data/derived/canonical/anatomy.json'

    def test_reproducible_retained_inventory(self):
        retained = json.loads((ROOT/'data/research/engineered_skin_territories/materialization.json').read_text())
        derived = lambda document: {k:v for k,v in document.items() if k!='anatomy_sha256'}
        self.assertEqual(derived(self.result), derived(retained))
        self.assertEqual(derived(generate()), derived(retained))
        anatomy_raw = (ROOT/self.ANATOMY_PATH).read_bytes()
        anatomy = json.loads(anatomy_raw)
        self.assertEqual(self.result['anatomy_sha256'], hashlib.sha256(anatomy_raw).hexdigest())
        entities = {e['id']:e for e in anatomy['entities']}
        self.assertEqual(len(retained['source_receipts']), 5)
        for receipt in retained['source_receipts']:
            reference = entities[receipt['entity_id']]['reference_geometry']
            self.assertEqual((reference['path'],reference['sha256'],reference['units'],reference['frame']),
                             (receipt['path'],receipt['sha256'],'m',anatomy['frame']['id']))
            raw = (ROOT/receipt['path']).read_bytes()
            self.assertEqual(len(raw), receipt['bytes'])
            self.assertEqual(hashlib.sha256(raw).hexdigest(), receipt['sha256'])
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
