"""Small ownership, donor and energy contract regressions; synthetic priors only."""
from pathlib import Path
import sys, unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.skin_epithelial_inventory import bind_inventory, run_passive, exact_patch_support, native_skin_owner


def fixture():
    owners={}
    for key,origin in [('film','external_surface_film'),('ic','native_skin_intracellular'),('ec','native_skin_extracellular')]:
        owners[key]=dict(origin=origin,accounting_owner=True,volume_m3=1e-6,
            mass_g=dict(Sodium=.0023,Potassium=.0039,Chloride=.00355,Albumin=.01))
    return dict(owners=owners,owner_ids=['film','ic','ec'],epoch='receipt-sha:step0',area_m2=1e-4,
        support_identity='synthetic-test-faces',epithelial_depth_m=1e-5,water_fraction=.8,
        cellular_fraction=.5,film_depth_m=1e-5,molar_mass_g_mol=[23,39,35.5],
        capacitance_f=[1e-8,1e-8],initial_potential_v=[-.01,-.02,0],temperature_k=310.,
        prior_source='synthetic engineering test only',countercharge_owner='unresolved-background-charge',
        initial_field_energy_owner='explicit-initial-electrical-state',heat_owner='explicit-thermal-reservoir')

class Checks(unittest.TestCase):
    def test_exact_held_exterior_and_native_units(self):
        import json
        root=Path(__file__).resolve().parents[1]
        evidence=(root/'data/research/engineered_skin_territories/materialization.json').read_bytes()
        e=json.loads(evidence);ref=dict(e['source_receipts'][0],units='m',frame='bodyparts3d-display-m')
        geometry=(root/ref['path']).read_bytes()
        r=exact_patch_support(ref,geometry,evidence,e['contact_eligible_triangle_ids'][:1])
        self.assertGreater(r['area_m2'],0)
        self.assertLess(r['area_m2'],r['source_support']['area_m2'])
        with self.assertRaises(ValueError):exact_patch_support(ref,geometry,evidence,[-1])
        n=native_skin_owner(dict(native_owner='Skin.intracellular',accounting_owner=True,volume_ml=1,mass_g={}),compartment='intracellular')
        self.assertEqual(n['volume_m3'],1e-6)

    def test_partition_all_species_and_no_mutation(self):
        args=fixture(); b=bind_inventory(**args)
        for key in args['owner_ids']:
            for species,mass in args['owners'][key]['mass_g'].items():
                self.assertAlmostEqual(b['represented'][key]['mass_g'][species]+b['complement'][key]['mass_g'][species],mass)
        self.assertEqual(args['owners']['ic']['mass_g']['Albumin'],.01)
        self.assertFalse(b['native_commit_available'])
    def test_nonowner_and_unknown_species_rejected(self):
        args=fixture();args['owners']['ec']['accounting_owner']=False
        with self.assertRaises(ValueError):bind_inventory(**args)
        args=fixture();args['owners']['ic']['mass_g']['Albumin']=None
        with self.assertRaises(ValueError):bind_inventory(**args)
    def test_missing_external_owner_and_overallocation(self):
        args=fixture();args['owners']['film']['origin']='native_skin_extracellular'
        with self.assertRaises(ValueError):bind_inventory(**args)
        args=fixture();args['epithelial_depth_m']=1
        with self.assertRaises(ValueError):bind_inventory(**args)
    def test_charge_voltage_energy_and_passive_reassembly(self):
        b=bind_inventory(**fixture());r=run_passive(b,np.ones((3,3))*1e-12,.001,3)
        np.testing.assert_allclose(r['species_delta_total_mol'],0,atol=1e-22)
        self.assertAlmostEqual(r['tep_v'],.01,places=4)
        self.assertNotEqual(r['tep_v'],r['basal_membrane_v'])
        self.assertGreaterEqual(r['heat_credit_j'],0)
        self.assertEqual(r['external_active_work_j'],0)
        for k in b['owner_ids']:
            self.assertEqual(r['proposed_owner_mass_g'][k]['Albumin'],fixture()['owners'][k]['mass_g']['Albumin'])
        np.testing.assert_allclose(b['fixed_countercharge_c']+96485.33212*(b['patch'].initial_moles@np.array([1,1,-1])),b['patch'].capacitance_matrix@b['patch'].initial_potential_v,atol=1e-16)
    def test_nonfinite_prior_rejected(self):
        args=fixture();args['water_fraction']=float('nan')
        with self.assertRaises(ValueError):bind_inventory(**args)

if __name__=='__main__':unittest.main()
