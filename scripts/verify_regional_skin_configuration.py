#!/usr/bin/env python3
"""Lightweight configured ownership, law-identity and actual-face mapping tests."""
from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.assembly.regional_skin_configuration import build_configuration,validate_configuration,map_contacts
from prepare_configured_regional_skin import patch_header,PARENT,stage

BASE=ROOT/'data/research/configured_regional_skin'
MATERIAL=ROOT/'data/research/engineered_skin_territories/materialization.json'


class Checks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.selection=json.loads((BASE/'example_selection.json').read_text())
        cls.config=build_configuration(MATERIAL,**cls.selection)

    def test_frozen_configuration_and_full_inventory_normalization(self):
        retained=json.loads((BASE/'prepared_v1/configuration.json').read_text())
        self.assertEqual(self.config,retained)
        w=validate_configuration(retained)
        self.assertEqual(sum(w),1.)
        self.assertEqual(set(retained['face_to_native_index']),{-1,0,1,2})
        self.assertEqual(sum(i>=0 for i in retained['face_to_native_index']),109183)
        self.assertTrue(retained['full_native_inventory_preserved'])
        for total in [0.,100.,1255.07579750261,24299070.4128279]:
            a,b=total*w[0],total*w[1];self.assertAlmostEqual(math.fsum([a,b,total-a-b]),total)
        # Algebraic baseline compliance, conductance and flow-source preservation.
        self.assertAlmostEqual(sum(2.*x for x in w),2.)
        self.assertAlmostEqual(sum(1/(100./x) for x in w),1/100.)
        self.assertAlmostEqual(sum(.03*x for x in w),.03)

    def test_source_patch_changes_only_identity_constants_and_receipt(self):
        patched=patch_header((PARENT/'native_regional_skin.h').read_text(),self.config)
        self.assertEqual(patched,(BASE/'prepared_v1/native_configured_regional_skin.h').read_text())
        self.assertIn('engineered_left_lower_leg',patched)
        self.assertNotIn('"region_a"',patched)
        for scope in ['ResistanceBaseline','NextResistance','ComplianceBaseline','NextCompliance','FlowSourceBaseline','NextFlowSource']:
            self.assertIn(scope,patched)
        with tempfile.TemporaryDirectory() as directory:
            receipt=stage(self.config,Path(directory)/'new')
            self.assertFalse(receipt['native_compiled'])
            with self.assertRaises(FileExistsError):stage(self.config,Path(directory)/'new')

    def test_reject_invalid_selection_prior_and_tampered_receipt(self):
        for field,value in [('selected_unions',[['unknown'],self.selection['selected_unions'][1]]),
                            ('region_names',['region_a','region_b','residual']),
                            ('region_names',['same','same','other'])]:
            selection={**self.selection,field:value}
            with self.assertRaises(ValueError):build_configuration(MATERIAL,**selection)
        for value in [0.,-1.,True,float('nan')]:
            selection=deepcopy(self.selection);selection['priors'][0]['thickness_m']=value
            with self.assertRaises(ValueError):build_configuration(MATERIAL,**selection)
        changed=deepcopy(self.config);changed['regions'][0]['initial_fraction']+=.01
        with self.assertRaises(ValueError):validate_configuration(changed)

    def test_named_contact_mapping_and_inner_or_oversized_rejection(self):
        face=next(i for i,r in enumerate(self.config['face_to_native_index']) if r==0)
        contact={'triangle_id':face,'area_m2':1e-8,'normal_pressure_pa':100.,'force_n':[0.,1e-6,0.]}
        result=map_contacts(self.config,[contact],root=ROOT)
        self.assertEqual(result['regions']['engineered_left_lower_leg']['normal_force_n'],1e-6)
        self.assertEqual(result['native_commands'],[])
        inner=next(i for i,r in enumerate(self.config['face_to_native_index']) if r==-1)
        for contacts in [[{**contact,'triangle_id':inner}],[{**contact,'area_m2':1.}],[contact,contact]]:
            with self.assertRaises(ValueError):map_contacts(self.config,contacts,root=ROOT)


if __name__=='__main__':unittest.main()
