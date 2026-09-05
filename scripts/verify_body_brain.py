"""Exercise imported brain law, canonical registration, causal response and timestep convergence."""
from pathlib import Path
import json
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


class BodyBrainVerification(unittest.TestCase):
    def test_executable_model_exists(self):
        self.assertTrue((ROOT / 'ihm/assembly/brain.py').exists(), 'Executable IBM brain adapter is missing')

    def test_physiology_response_and_timestep_convergence(self):
        from ihm.assembly.brain import BodyBrain
        data = json.loads((ROOT / 'data/derived/canonical/brain.json').read_text())
        baseline = BodyBrain.from_dict(data)
        hypoxia = BodyBrain.from_dict(data)
        coarse = BodyBrain.from_dict(data, max_step_s=0.001)
        fine = BodyBrain.from_dict(data, max_step_s=0.0005)
        normal = {'mean_arterial_pressure_mmHg': 90., 'oxygen_saturation': .98, 'core_temperature_C': 37.}
        low = dict(normal, mean_arterial_pressure_mmHg=45., oxygen_saturation=.65)
        b = baseline.step(1., normal)
        h = hypoxia.step(1., low)
        c = coarse.step(.15, low)
        f = fine.step(.15, low)
        self.assertGreater(np.mean(b['regional_state']['activity_hz']), np.mean(h['regional_state']['activity_hz']))
        self.assertGreater(h['autonomic_commands']['sympathetic_fraction'], b['autonomic_commands']['sympathetic_fraction'])
        np.testing.assert_allclose(c['regional_state']['activity_hz'], f['regional_state']['activity_hz'], atol=2e-5, rtol=2e-5)
        self.assertFalse(h['autonomic_commands']['applied_to_body'])
        self.assertGreater(np.mean(b['regional_state']['activity_hz']), 1.)
        self.assertLess(np.max(b['regional_state']['activity_hz']), 100.)
        self.assertEqual(h['input_provenance']['mean_arterial_pressure_mmHg'], 'caller_supplied')
        self.assertAlmostEqual(h['time_s'], 1.)
        with self.assertRaises(ValueError):
            baseline.step(-1., normal)
        with self.assertRaises(ValueError):
            baseline.step(1., dict(normal, oxygen_saturation=98.))
        with self.assertRaises(ValueError):
            baseline.step(1., dict(normal, core_temperature_C=float('nan')))

    def test_source_preservation_registration_and_topology(self):
        from ihm.assembly.brain import BodyBrain, verify_sources
        data = json.loads((ROOT / 'data/derived/canonical/brain.json').read_text())
        verified = verify_sources(data)
        self.assertGreaterEqual(len(verified), 8)
        import copy
        corrupt = copy.deepcopy(data)
        corrupt['sources'][0]['sha256'] = '0' * 64
        with self.assertRaises(ValueError):
            BodyBrain.from_dict(corrupt)
        cortex = [n for n in data['nodes'] if n['kind'] == 'cortical_population']
        self.assertEqual(len(cortex), 68)
        self.assertEqual(len(set(n['id'] for n in data['nodes'])), len(data['nodes']))
        self.assertTrue(all(n['body_entity_id'].startswith('body-bp3d-') for n in data['nodes']))
        self.assertGreater(len(data['edges']), len(cortex))
        self.assertTrue(all(e['evidence_kind'] == 'geometric_prior' for e in data['edges']))
        self.assertGreater(data['registration']['landmark_rms_m'], 0.)
        self.assertLess(data['registration']['landmark_rms_m'], .06)
        positions = np.asarray([n['position_m'] for n in cortex])
        self.assertTrue(np.all(positions[:,1] > .6))
        left = [n['position_m'][0] for n in cortex if n['hemisphere']=='left']
        right = [n['position_m'][0] for n in cortex if n['hemisphere']=='right']
        self.assertGreater(np.mean(left), np.mean(right))
        model = BodyBrain.from_dict(data)
        independent = BodyBrain.from_dict(data)
        original = independent.state.copy()
        output = model.step(.03, {})
        np.testing.assert_array_equal(independent.state, original)
        self.assertEqual(len(output['regional_state']['activity_hz']), len(data['nodes']))
        self.assertTrue(np.isfinite(output['regional_state']['potential_mV']).all())
        self.assertEqual(output['input_provenance']['oxygen_saturation'], 'assumed_baseline')
        # Imported IBM law's initial membrane derivative for V=-65, adaptation=0, drive=2.
        derivative = model.source_rate_law({'neural.exc.potential': np.array([-65.]),
            'neural.exc.adaptation': np.array([0.]), 'neural.exc.ampa': np.array([2.]),
            'neural.exc.nmda': np.array([0.]), 'neural.exc.activity': np.array([0.])}, {})
        self.assertAlmostEqual(derivative['neural.exc.potential'][0], 133.33333333333334)
        self.assertAlmostEqual(derivative['neural.exc.activity'][0], 1517.163600424871)


if __name__ == '__main__':
    unittest.main()
