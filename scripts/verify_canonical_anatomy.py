#!/usr/bin/env python3
"""Verify canonical assembly identity, provenance, registration and connections."""
import argparse
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

class RegistrationTests(unittest.TestCase):
    def test_affine_recovery_and_inverse_handedness(self):
        from ihm.assembly.anatomy import LandmarkRegistration
        source = np.random.default_rng(21).normal(size=(20, 3))
        target = source * [1.02, .98, 1.05] + [.1, -.2, .3]
        fit = LandmarkRegistration(source, target)
        np.testing.assert_allclose(fit.transform(source), target, atol=1e-8)
        self.assertGreater(fit.report()['affine_determinant'], 0)

    def test_degenerate_landmarks_rejected(self):
        from ihm.assembly.anatomy import LandmarkRegistration
        with self.assertRaises(ValueError):
            LandmarkRegistration(np.zeros((5, 3)), np.zeros((5, 3)))

    def test_role_classification(self):
        from ihm.assembly.anatomy import physical_role
        self.assertEqual(physical_role('left femur', 'skeletal'), 'rigid_bone')
        self.assertEqual(physical_role('costal cartilage', 'skeletal'), 'cartilage')
        self.assertEqual(physical_role('Achilles tendon', 'muscular'), 'tendon')
        self.assertEqual(physical_role('skin', 'integumentary'), 'skin')
        self.assertEqual(physical_role('cavity of left ventricle', 'cardiac'), 'fluid_cavity')
        self.assertEqual(physical_role('Node of ligamentum arteriosum', 'lymphatic'), 'lymph_node_group')


def verify(path):
    from ihm.assembly.anatomy import verify_assembly
    return verify_assembly(json.loads(Path(path).read_text()), ROOT)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--unit-only', action='store_true')
    parser.add_argument('--path', default=str(ROOT/'data/derived/canonical/anatomy.json'))
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RegistrationTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    if not args.unit_only:
        report = verify(args.path)
        print(json.dumps(report, indent=2))
