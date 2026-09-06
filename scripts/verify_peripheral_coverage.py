"""Lightweight peripheral identity/provenance checks; no native execution."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class PeripheralCoverageTests(unittest.TestCase):
    def test_module_available(self):
        self.assertTrue((ROOT / 'ihm/assembly/peripheral_coverage.py').is_file(),
                        'Peripheral coverage materializer is missing')

    def test_exact_registered_coverage_and_reflex_dependencies(self):
        from ihm.assembly.peripheral_coverage import build_peripheral_coverage
        report = build_peripheral_coverage(ROOT)
        self.assertEqual(report['summary']['effector_count'], 92)
        self.assertEqual(report['summary']['source_reflex_effector_count'], 8)
        self.assertEqual(report['summary']['anatomically_measured_connection_count'], 0)
        rows = {r['id']: r for r in report['effectors']}
        self.assertEqual(len(report['receptors']), 186)
        self.assertEqual(rows['arm26_BIClong_r']['spinal_reflex']['status'], 'unsupported')
        self.assertIn('soleus_r', rows['tibant_r']['spinal_reflex']['input_effectors'])
        self.assertEqual(rows['gasmed_r']['spinal_reflex']['input_effectors'], ['gasmed_r', 'gaslat_r'])
        self.assertFalse(report['biological_validation'])
        json.dumps(report, allow_nan=False)

    def test_reject_controller_identity_and_assignment_drift(self):
        from ihm.assembly.peripheral_coverage import build_peripheral_coverage
        from ihm.assembly.sensorimotor_catalog import whole_body_effector_catalog
        manifest = json.loads((ROOT / 'data/derived/mechanics/whole_body_arm26_v2/registration.json').read_text())
        catalog = whole_body_effector_catalog(ROOT, manifest)
        for mutate in (lambda c: c.pop(), lambda c: c[0].update(id='invented'),
                       lambda c: c[0].update(motor_region='brain-missing'),
                       lambda c: c[0].update(side='l'),
                       lambda c: c[0].update(attachment_bodies=['invented'])):
            bad = deepcopy(catalog)
            mutate(bad)
            with self.assertRaises(ValueError):
                build_peripheral_coverage(ROOT, controller_catalog=bad)

    def test_reject_stale_model_and_escape(self):
        from ihm.assembly.peripheral_coverage import build_peripheral_coverage
        manifest = json.loads((ROOT / 'data/derived/mechanics/whole_body_arm26_v2/registration.json').read_text())
        for field, value in [('model_sha256', '0' * 64), ('model_path', '../outside.osim')]:
            bad = deepcopy(manifest)
            bad[field] = value
            with self.assertRaises(ValueError):
                build_peripheral_coverage(ROOT, registration=bad)

    def test_materialized_roundtrip(self):
        from ihm.assembly.peripheral_coverage import write_peripheral_coverage
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'coverage.json'
            report = write_peripheral_coverage(ROOT, path)
            self.assertEqual(json.loads(path.read_text()), report)


if __name__ == '__main__':
    unittest.main()
