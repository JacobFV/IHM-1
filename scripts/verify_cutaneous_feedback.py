"""Bounded real pinned-IBM cutaneous fixtures; no native jobs."""
from copy import deepcopy
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def sites():
    return [{'id': 'left-forearm-fixture', 'position_m': [0., 0., 0.],
             'normal': [0., 0., 1.], 'contact_area_m2': .001,
             'stiffness_pa_per_m': 1e7, 'sensory_region': 'brain-rh-postcentral',
             'reference_temperature_C': 33., 'support_basis': 'explicit synthetic fixture; not anatomical registration'}]


class CutaneousTests(unittest.TestCase):
    def test_module_available(self):
        self.assertTrue((ROOT / 'ihm/assembly/cutaneous_feedback.py').is_file(),
                        'Executable cutaneous bridge missing')

    def make(self, **kwargs):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        return CutaneousFeedback(ROOT, sites=sites(), recruitment_hz_per_response=.1, **kwargs)

    def sample(self, model, force=1., temperature=33.):
        return {'time_s': model.time_s,
                'contacts': [{'id': 'left-forearm-fixture', 'force_n': [0., 0., -force]}],
                'skin_temperature_C': temperature}

    def test_physical_units_and_source_operator(self):
        model = self.make()
        result = model.step(.01, self.sample(model))
        sensed = result['sites'][0]
        self.assertEqual(sensed['normal_force_n'], 1.)
        self.assertEqual(sensed['pressure_pa'], 1000.)
        self.assertAlmostEqual(sensed['indentation_um'], 100.)
        self.assertGreater(sensed['rapid_response'], 0.)
        self.assertGreater(result['sensory_inputs_hz']['brain-rh-postcentral'], 0.)
        self.assertLess(model.audit['receptors']['rapid']['donor_relative_error'], 1e-9)
        self.assertEqual(model.audit['receptors']['thermal']['implementation'], 'thermoreceptor_static_dynamic')

    def test_delay_block_release_and_real_neural_effect(self):
        from ihm.assembly.brain import BodyBrain
        import json
        model = self.make(delay_s=.02)
        first = model.step(.01, self.sample(model))
        self.assertEqual(first['sensory_inputs_hz'], {})
        blocked = model.step(.01, self.sample(model), sensory_blocks=['left-forearm-fixture'])
        self.assertEqual(blocked['sensory_inputs_hz'], {})
        for _ in range(3):
            blocked = model.step(.01, self.sample(model, force=0.), sensory_blocks=['left-forearm-fixture'])
        self.assertEqual(blocked['sensory_inputs_hz'], {})
        restored = model.step(.01, self.sample(model))
        self.assertEqual(restored['sensory_inputs_hz'], {})
        for _ in range(3):
            restored = model.step(.01, self.sample(model))
        self.assertGreater(restored['sensory_inputs_hz']['brain-rh-postcentral'], 0.)
        for _ in range(4):
            released = model.step(.01, self.sample(model, force=0.))
        self.assertLess(released['sites'][0]['rapid_response'], 0.)
        brain_data = json.loads((ROOT / 'data/derived/canonical/brain.json').read_text())
        a, b = BodyBrain(brain_data, root=ROOT), BodyBrain(brain_data, root=ROOT)
        x = a.step(.01, sensory_inputs_hz=restored['sensory_inputs_hz'])
        y = b.step(.01)
        self.assertNotEqual(x['regional_state']['activity_hz'], y['regional_state']['activity_hz'])

    def test_thermal_transient_and_failure_atomicity(self):
        model = self.make()
        warm = model.step(.01, self.sample(model, force=0., temperature=35.))
        self.assertGreater(warm['sites'][0]['thermal_response'], 0.)
        self.assertEqual(warm['sites'][0]['temperature_change_C'], 2.)
        before = model.checkpoint()
        bad = self.sample(model)
        bad['contacts'][0]['id'] = 'unknown'
        with self.assertRaises(ValueError):
            model.step(.01, bad)
        self.assertEqual(before, model.checkpoint())
        a = model.step(.01, self.sample(model))
        model.restore(before)
        b = model.step(.01, self.sample(model))
        self.assertEqual(a, b)
        corrupted = deepcopy(before)
        corrupted['sites']['left-forearm-fixture']['rapid']['identity'] = 'wrong'
        with self.assertRaises(ValueError):
            model.restore(corrupted)
        self.assertEqual(b['time_s'], model.time_s)

    def native_site(self):
        native = sites()
        native[0].pop('stiffness_pa_per_m')
        native[0].update(mechanical_input='native_indentation',
            material_identity={'manifest_sha256':'a'*64,'quadrature_index':3,'triangle_index':7},
            indentation_basis='native confined-layer compression fixture',
            area_basis='projected reference quadrature area')
        return native

    def test_native_indentation_bypasses_stiffness_and_preserves_identity(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        native = self.native_site()
        model = CutaneousFeedback(ROOT, sites=native, recruitment_hz_per_response=.1)
        sample = self.sample(model)
        sample['contacts'][0].update(indentation_m=.00025,
            material_identity=native[0]['material_identity'],
            indentation_basis=native[0]['indentation_basis'])
        result = model.step(.01, sample)
        row = result['sites'][0]
        self.assertEqual(row['indentation_um'], 250.)
        self.assertEqual(row['pressure_pa'], 1000.)
        self.assertEqual(row['material_identity'], native[0]['material_identity'])
        self.assertEqual(row['mechanical_status'], 'native_modeled_indentation')
        self.assertGreater(row['rapid_response'], 0.)
        released = model.step(.01, {'time_s':model.time_s, 'contacts':[], 'skin_temperature_C':33.})
        self.assertEqual(released['sites'][0]['indentation_um'], 0.)

    def test_native_moving_site_uses_current_normal_and_keeps_reference_prior(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        native = self.native_site()
        model = CutaneousFeedback(ROOT, sites=native, recruitment_hz_per_response=.1)
        sample = self.sample(model)
        sample['contacts'][0].update(indentation_m=.00025,
            material_identity=native[0]['material_identity'],
            indentation_basis=native[0]['indentation_basis'],
            point_m=[1., 2., 3.], normal=[1., 0., 0.], force_n=[-2., 0., 0.])
        result = model.step(.01, sample)
        row = result['sites'][0]
        self.assertEqual(row['normal_force_n'], 2.)
        self.assertEqual(row['pressure_pa'], 2000.)
        self.assertEqual(row['position_m'], [1., 2., 3.])
        self.assertEqual(row['outward_normal'], [1., 0., 0.])
        self.assertEqual(row['reference_position_m'], [0., 0., 0.])
        self.assertEqual(row['position_basis'], 'current native material point; source receptor support remains reference prior')
        self.assertEqual(model.audit['receptors']['rapid']['sites_m'], [[0., 0., 0.]])
        self.assertEqual(model.sites[native[0]['id']]['position_m'], [0., 0., 0.])

    def test_native_moving_site_rejects_partial_or_invalid_geometry(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        native = self.native_site()
        model = CutaneousFeedback(ROOT, sites=native, recruitment_hz_per_response=.1)
        for geometry in ({'point_m':[1., 0., 0.]}, {'normal':[1., 0., 0.]},
                         {'point_m':[1., 0., 0.], 'normal':[2., 0., 0.]},
                         {'point_m':[float('nan'), 0., 0.], 'normal':[1., 0., 0.]}):
            sample = self.sample(model)
            sample['contacts'][0].update(indentation_m=.0001,
                material_identity=native[0]['material_identity'],
                indentation_basis=native[0]['indentation_basis'], **geometry)
            before = model.checkpoint()
            with self.assertRaises(ValueError):
                model.step(.01, sample)
            self.assertEqual(before, model.checkpoint())

    def test_native_indentation_does_not_require_pressure_or_stiffness(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        native = self.native_site()
        native[0]['contact_area_m2'] = None
        model = CutaneousFeedback(ROOT, sites=native, recruitment_hz_per_response=.1)
        sample = self.sample(model)
        sample['contacts'][0].update(indentation_m=.00025,
            material_identity=native[0]['material_identity'],
            indentation_basis=native[0]['indentation_basis'])
        result = model.step(.01, sample)
        self.assertIsNone(result['sites'][0]['pressure_pa'])
        self.assertEqual(result['sites'][0]['indentation_um'], 250.)
        self.assertGreater(result['sites'][0]['rapid_response'], 0.)

    def test_native_indentation_rejects_missing_or_mismatched_material_receipt(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        native = self.native_site()
        model = CutaneousFeedback(ROOT, sites=native, recruitment_hz_per_response=.1)
        for fields in ({}, {'indentation_m':-.1}, {'indentation_m':float('nan')},
                       {'indentation_m':.0001,'material_identity':{'triangle_index':8}},
                       {'indentation_m':.0001,'material_identity':native[0]['material_identity'],
                        'indentation_basis':'invented basis'}):
            before = model.checkpoint()
            sample = self.sample(model)
            sample['contacts'][0].update(fields)
            with self.assertRaises(ValueError):
                model.step(.01, sample)
            self.assertEqual(before, model.checkpoint())

    def test_missing_temperature_is_not_a_reference_temperature_sample(self):
        model = self.make()
        model.step(.01, self.sample(model, force=0., temperature=35.))
        missing = model.step(.01, self.sample(model, force=0., temperature=None))
        self.assertEqual(missing['sites'][0]['thermal_response'], 0.)
        self.assertIsNone(missing['sites'][0]['temperature_change_C'])

    def test_force_orientation_point_identity_and_nonfinite_inputs(self):
        model = self.make()
        outward = model.step(.01, self.sample(model, force=-1.))
        self.assertEqual(outward['sites'][0]['normal_force_n'], 0.)
        self.assertEqual(outward['sites'][0]['indentation_um'], 0.)
        before = model.checkpoint()
        for field, value in [('force_n', [0., 0., float('nan')]),
                             ('point_m', [1., 0., 0.])]:
            bad = self.sample(model)
            bad['contacts'][0][field] = value
            with self.assertRaises(ValueError):
                model.step(.01, bad)
            self.assertEqual(before, model.checkpoint())

    def test_unknown_area_preserves_force_without_pressure(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        unknown = sites()
        unknown[0]['contact_area_m2'] = None
        model = CutaneousFeedback(ROOT, sites=unknown, recruitment_hz_per_response=.1)
        result = model.step(.01, self.sample(model))
        self.assertEqual(result['sites'][0]['force_n'], [0., 0., -1.])
        self.assertIsNone(result['sites'][0]['pressure_pa'])
        self.assertEqual(result['sensory_inputs_hz'], {})

    def test_unknown_area_and_invalid_supports(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        for field, value in [('contact_area_m2', 0.), ('normal', [0., 0., 0.]),
                             ('sensory_region', 'invented'), ('stiffness_pa_per_m', 0.)]:
            bad = sites()
            bad[0][field] = value
            with self.assertRaises(ValueError):
                CutaneousFeedback(ROOT, sites=bad, recruitment_hz_per_response=.1)


if __name__ == '__main__':
    unittest.main()
