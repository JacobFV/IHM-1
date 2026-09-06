"""Bounded tests for research-only hip restraint calibration."""
import hashlib,json,unittest
from pathlib import Path
import numpy as np
from ihm.assembly.hip_capsule import FixedPoseRestraint

ROOT=Path(__file__).resolve().parents[1]
class CapsuleTests(unittest.TestCase):
    def setUp(self):
        self.model=FixedPoseRestraint.from_file(ROOT/'data/research/hip_capsule/neutral_slice.json')
    def test_measured_endpoint_and_rad_units(self):
        m=self.model
        for sign in [-1,1]:
            q=m.center+sign*(m.half_slack+m.toe)
            e,t,k=m.evaluate(q,flexion_rad=0,adduction_rad=0)
            self.assertAlmostEqual(t,-sign*5)
            self.assertAlmostEqual(k,np.rad2deg(m.stiffness_per_degree))
            self.assertGreater(e,0)
    def test_energy_gradient_and_passivity(self):
        m=self.model
        for sign in [-1,1]:
            q=m.center+sign*(m.half_slack+.6*m.toe)
            e,t,k=m.evaluate(q,flexion_rad=0,adduction_rad=0)
            h=1e-7
            ep=m.evaluate(q+h,flexion_rad=0,adduction_rad=0)[0]
            em=m.evaluate(q-h,flexion_rad=0,adduction_rad=0)[0]
            self.assertAlmostEqual(t,-(ep-em)/(2*h),places=7)
            self.assertLess(t*(q-m.center),0)
            self.assertGreaterEqual(k,0)
    def test_slack_and_toe_continuity(self):
        m=self.model
        self.assertEqual(m.evaluate(m.center,flexion_rad=0,adduction_rad=0),(0.,0.,0.))
        for sign in [-1,1]:
            q=m.center+sign*m.half_slack
            self.assertEqual(m.evaluate(q,flexion_rad=0,adduction_rad=0),(0.,0.,0.))
            self.assertLess(abs(m.evaluate(q+sign*1e-8,flexion_rad=0,adduction_rad=0)[1]),1e-10)
    def test_reject_unmeasured_pose_and_extrapolation(self):
        m=self.model
        for f,a,q in [(1e-4,0,m.center),(0,1e-4,m.center),(0,0,float('nan')),(0,0,np.deg2rad(40)),(0,0,np.deg2rad(-40))]:
            with self.assertRaises(ValueError):m.evaluate(q,flexion_rad=f,adduction_rad=a)
        with self.assertRaises(ValueError):m.evaluate(m.center,flexion_rad=0,adduction_rad=0,compression_N=0)

    def test_receipts_and_no_native_permission(self):
        folder=ROOT/'data/research/hip_capsule'
        cards=json.loads((folder/'sources.json').read_bytes())
        self.assertFalse(cards['native_activation_allowed'])
        for source in cards['retained_files']:
            raw=(folder/source['file']).read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(),source['sha256'])
            self.assertEqual(len(raw),source['bytes'])

if __name__=='__main__':unittest.main()
