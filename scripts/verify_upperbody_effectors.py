"""Small source/registration/controller tests; no native compilation."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
from scripts.build_upperbody_effectors import build_upperbody_model, ARM_SOURCE, TARGET_SOURCE
from ihm.assembly.sensorimotor_catalog import whole_body_effector_catalog
from ihm.assembly.sensorimotor import SensorimotorController
from scripts.verify_sensorimotor import observation
ROOT=Path(__file__).resolve().parents[1]

class UpperbodyTests(unittest.TestCase):
    def test_registered_source_preservation_and_controller(self):
        before={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (ARM_SOURCE,TARGET_SOURCE)}
        with tempfile.TemporaryDirectory(dir=ROOT/'data/derived/audits',prefix='arms-test-') as tmp:
            manifest=build_upperbody_model(ROOT,Path(tmp))
            catalog=whole_body_effector_catalog(ROOT,manifest)
            self.assertEqual(len(catalog),92)
            self.assertEqual(len({r['id'] for r in catalog}),92)
            model=ET.parse(ROOT/manifest['model_path']).getroot()
            target=ET.parse(ROOT/TARGET_SOURCE).getroot()
            self.assertEqual(ET.tostring(model.find('.//JointSet')),ET.tostring(target.find('.//JointSet')))
            masses=lambda t:{b.get('name'):b.findtext('mass') for b in t.findall('.//BodySet/objects/Body')}
            self.assertEqual(masses(model),masses(target))
            self.assertEqual(len(model.findall('.//Thelen2003Muscle')),12)
            self.assertEqual(len([w for w in model.findall('.//WrapObjectSet/objects/*') if w.get('name','').startswith('arm26_')]),8)
            for row in manifest['registrations']:
                self.assertLess(row['elbow_landmark_residual_m'],1e-12)
                self.assertAlmostEqual(abs(np.linalg.det(row['linear_map'])/(row['length_scale']**3)),1.)
            c=SensorimotorController.from_root(ROOT,muscle_catalog=catalog)
            for _ in range(5):
                obs=observation(c.time_s)
                for row in catalog:
                    obs['muscles'].setdefault(row['id'],dict(fiber_length_m=row['optimal_fiber_length_m'],optimal_fiber_length_m=row['optimal_fiber_length_m'],tendon_force_n=0.,max_isometric_force_n=row['max_isometric_force_n'],sensor_basis='source-shaped synthetic CE test'))
                out=c.step(.02,obs,descending={'arm26_BRA_r':1.})
            self.assertGreater(out['motor_excitations']['arm26_BRA_r'],0)
            self.assertEqual(out['motor_excitations']['arm26_BRA_l'],0)
            self.assertEqual(out['motor_excitations']['arm26_TRIlong_r'],0)
            path=ROOT/manifest['model_path'];path.write_text(path.read_text()+'\n')
            with self.assertRaises(ValueError):whole_body_effector_catalog(ROOT,manifest)
        self.assertEqual(before,{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in before})
if __name__=='__main__':unittest.main()
