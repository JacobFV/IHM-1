"""Small finite-reservoir conservation/energy and observation-scope checks."""
from pathlib import Path
import sys, unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.epithelial_electrodiffusion import EpithelialPatch, simulate, prior_ensemble, fit_voltage_observations, condition_tep_magnitude, fit_tape_strip_response

def patch(**kw):
    args=dict(volumes_m3=[1e-12,1e-15,1e-12], concentrations_mol_m3=[[140,4,110],[15,140,30],[140,4,110]],
        capacitance_f=[1e-11,1e-11], conductance_s=[[1e-11,1e-12,1e-11],[1e-12,1e-10,1e-11],[1e-12,1e-12,1e-12]],
        initial_potential_v=[-.03,-.024,0],temperature_k=310.15,
        provenance={'preparation':'synthetic_numeric_test','parameters':'explicit_unidentified_priors'})
    args.update(kw);return EpithelialPatch(**args)

class Checks(unittest.TestCase):
    def test_retained_human_extraction_and_bytes(self):
        import json,hashlib
        from build_epithelial_electrodiffusion import extract,DIRECTORY,ROOT
        rows=extract()
        self.assertEqual([r['tape_strips'] for r in rows],[0,1,4,7])
        self.assertAlmostEqual(rows[0]['reported_tep_mV'],-10.122183)
        receipt=json.loads((DIRECTORY/'source-receipt.json').read_text())
        self.assertEqual(receipt['participant_count'],1)
        for source in receipt['raw_files']:
            self.assertEqual(hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest(),source['sha256'])
    def test_conditioning_preserves_observable_scope(self):
        p=patch()
        row=dict(value_v=.01,quantity='tep_magnitude',preparation='synthetic_numeric_test',source_id='held')
        c=condition_tep_magnitude(p,row,basal_positive=True)
        self.assertAlmostEqual(simulate(c,.01,2)['tep_v'][0],.01)
        self.assertEqual(c.initial_potential_v[1],p.initial_potential_v[1])
        with self.assertRaises(ValueError):condition_tep_magnitude(p,dict(row,quantity='cell_vm'),basal_positive=True)
    def test_strip_fit_only_identifies_phenotype(self):
        rows=[dict(tape_strips=x,value_v=.01-.0002*x) for x in [0,1,4,7]]
        r=fit_tape_strip_response(rows)
        self.assertAlmostEqual(r['slope_v_per_strip'],-.0002)
        self.assertLess(abs(r['holdout_error_v']),1e-15)
        self.assertFalse(r['rc_parameters_identified'])
    def test_closed_species_and_charge(self):
        r=simulate(patch(),.2,21)
        self.assertLess(r['audit']['max_species_residual_mol'],1e-24)
        self.assertLess(r['audit']['max_charge_residual_c'],1e-22)
        self.assertGreater(np.min(r['moles']),0)
        self.assertFalse(r['ownership']['native_blood_mutated'])
    def test_passive_free_energy_dissipation(self):
        r=simulate(patch(),.2,21)
        self.assertTrue(np.all(np.diff(r['free_energy_change_j'])<=1e-20))
        self.assertGreater(r['dissipated_energy_j'][-1],0)
        self.assertLess(r['audit']['relative_energy_residual'],1e-5)
    def test_active_work_closes_energy(self):
        r=simulate(patch(active_flux_mol_s=[[0,0,0],[3e-18,-2e-18,0],[0,0,0]]),.2,21)
        self.assertGreater(abs(r['active_work_j'][-1]),1e-18)
        self.assertLess(r['audit']['relative_energy_residual'],1e-5)
    def test_membrane_and_tep_distinct_and_coupled(self):
        r=simulate(patch(),.2,21)
        self.assertAlmostEqual(r['tep_v'][0],.03)
        self.assertAlmostEqual(r['basal_membrane_v'][0],-.024)
        self.assertAlmostEqual(r['apical_membrane_v'][0],.006)
        self.assertNotAlmostEqual(r['tep_v'][-1],r['tep_v'][0],places=7)
    def test_equilibrium_and_no_flux(self):
        p=patch(concentrations_mol_m3=[[100,100,100]]*3,initial_potential_v=[0,0,0])
        r=simulate(p,.1,3)
        np.testing.assert_allclose(r['potential_v'],0,atol=1e-15)
        self.assertEqual(r['audit']['relative_energy_residual'],0)
    def test_missing_prior_and_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):patch(capacitance_f=[-1,1])
        with self.assertRaises(ValueError):patch(provenance={})
        with self.assertRaises(ValueError):simulate(patch(),float('nan'))
        with self.assertRaises(ValueError):patch(conductance_s=[[1,2,3]]*2)
    def test_ensemble_is_prior_spread_not_posterior(self):
        r=prior_ensemble([patch(),patch(capacitance_f=[2e-11,2e-11])],.1,4)
        self.assertEqual(r['uncertainty_kind'],'explicit_prior_scenario_envelope_not_posterior')
        self.assertEqual(r['member_count'],2)
        self.assertTrue(np.any(np.array(r['tep_max_v'])>np.array(r['tep_min_v'])))
    def test_observation_fit_cannot_identify_channels(self):
        rows=[dict(value_v=.03,quantity='tep',preparation='human_intact',source_id='x'),dict(value_v=.02,quantity='tep',preparation='human_intact',source_id='x')]
        r=fit_voltage_observations(rows)
        self.assertEqual(r['identified_parameters'],['voltage_location_v'])
        self.assertAlmostEqual(r['voltage_location_v'],.025)
        with self.assertRaises(ValueError):fit_voltage_observations([rows[0],dict(rows[1],quantity='cell_vm')])
        with self.assertRaises(ValueError):fit_voltage_observations([rows[0],dict(rows[1],preparation='cultured_keratinocyte')])
if __name__=='__main__':unittest.main()
