#!/usr/bin/env python3
"""Validator tests with labeled synthetic records; no native engine execution."""
from copy import deepcopy
import json
import unittest
from verify_configured_skin_engine import prepare,assess,ROOT


def synthetic_frames(receipt):
    names=receipt['region_names'];fractions=receipt['fractions']
    subs=json.loads((ROOT/'data/research/configured_skin_species/prepared_v2/preparation.json').read_text())['species_names']
    frames={}
    for mode in ['baseline','zero','load']:
        frames[mode]=[]
        for tick in range(13):
            v={'temperature_c':37.,'acceptance.liquid_owners.external_boundary_count':1,
               'acceptance.boundary.Ambient.volume_is_positive_infinity':1,'acceptance.liquid_owners.count':1 if mode=='baseline' else 3,
               'acceptance.liquid_owners.volume_ml':100.,'acceptance.skin.volume_ownership_residual_ml':0.}
            for sub in subs:
                v['acceptance.skin.species.'+sub+'.mass_ug']=1.
                v['acceptance.liquid_owners.species.'+sub+'.mass_ug']=2.
            if mode!='baseline':
                v['acceptance.original_path.SkinL1ToSkinL2.resistance_mmhg_s_per_ml']=100.
                for i,name in enumerate(names):
                    p='acceptance.region.'+name+'.';loaded=mode=='load' and tick in (6,7) and i==0
                    v[p+'initial_fraction']=fractions[i]
                    v[p+'path.SkinL1ToSkinL2.resistance_mmhg_s_per_ml']=100./fractions[i]
                    v[p+'external_pressure_pa']=133.322387415 if loaded else 0.
                    v[p+'pressure_mmhg']=1. if loaded else 0.
                    v[p+'lymph_flow_ml_per_s']=.2 if loaded else .1
                    for sub in subs:v[p+sub+'.mass_ug']=fractions[i]
            frames[mode].append({'tick':tick,'elapsed_s':tick*.02,'configuration_sha256':receipt['configuration_sha256'],'values':v})
    return frames


class Checks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.receipt=prepare()

    def test_accepts_explicit_valid_synthetic_records(self):
        result=assess(synthetic_frames(self.receipt),self.receipt)
        self.assertTrue(result['passed']);self.assertEqual(len(result['species_names']),28)
        self.assertFalse(result['whole_body_species_closure_claimed'])

    def test_rejects_duplicate_owner_missing_species_and_stale_law(self):
        name=self.receipt['region_names'][0]
        for key,value in [('acceptance.liquid_owners.count',4),
                          ('acceptance.region.'+name+'.Albumin.mass_ug',None),
                          ('acceptance.region.'+name+'.path.SkinL1ToSkinL2.resistance_mmhg_s_per_ml',100.),
                          ('temperature_c',None)]:
            frames=synthetic_frames(self.receipt);frames['zero'][3]['values'][key]=value
            with self.subTest(key=key):self.assertFalse(assess(frames,self.receipt)['passed'])

    def test_rejects_aggregate_mass_parity_and_configuration_drift(self):
        frames=synthetic_frames(self.receipt)
        frames['zero'][4]['values']['acceptance.liquid_owners.species.Albumin.mass_ug']+=.1
        self.assertFalse(assess(frames,self.receipt)['passed'])
        frames=synthetic_frames(self.receipt);frames['load'][7]['configuration_sha256']='changed'
        self.assertFalse(assess(frames,self.receipt)['passed'])


if __name__=='__main__':unittest.main()
