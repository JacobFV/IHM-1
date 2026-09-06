"""Short actual-geometry free trajectory refinement, no native or recoil."""
import unittest
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_free_dynamics import ThoracicFreeDynamics
from ihm.assembly.thoracic_trajectory import FreeThoraxTrajectory, rotation_increment

ROOT=Path(__file__).resolve().parents[1]

class TrajectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=ThoracicMechanism(ROOT/'data/research/thoracic_mechanism/v2/manifest.json')
        cls.flow=FreeThoraxTrajectory(ThoracicFreeDynamics(cls.model))
        q=np.zeros(26);q[0]=.002;q[24]=.0003;q[25]=.001
        u=np.linspace(-.01,.02,32);u[:6]=[.2,.1,-.15,.6,-.4,.5]
        for k in cls.model.locked:u[6+k]=0
        cls.initial={'q':q,'velocity':u,'translation_m':np.array([.3,-.2,.4]),'rotation':rotation_increment([.1,.2,-.1])}
        cls.start=cls.flow.invariants(cls.initial)
        cls.coarse=cls.flow.step(cls.initial,.05)
        cls.fine=cls.flow.step(cls.flow.step(cls.initial,.025),.025)
        cls.coarse_inv=cls.flow.invariants(cls.coarse);cls.fine_inv=cls.flow.invariants(cls.fine)

    def test_free_invariant_error_decreases_under_refinement(self):
        for key in ['kinetic_energy_J','world_linear_momentum_kg_m_per_s','world_angular_momentum_kg_m2_per_s']:
            coarse=np.linalg.norm(np.asarray(self.coarse_inv[key])-self.start[key])
            fine=np.linalg.norm(np.asarray(self.fine_inv[key])-self.start[key])
            self.assertLess(fine,.4*coarse+2e-13,(key,coarse,fine))
            self.assertLess(fine,2e-5)
        self.assertAlmostEqual(self.fine_inv['mass_kg'],20.235757397039237,places=10)

    def test_rotation_locks_and_input_ownership(self):
        for state in [self.coarse,self.fine]:
            r=state['rotation'];np.testing.assert_allclose(r.T@r,np.eye(3),atol=1e-14)
            self.assertAlmostEqual(np.linalg.det(r),1,places=14)
            for k in self.model.locked:self.assertEqual(state['q'][k],0);self.assertEqual(state['velocity'][6+k],0)
        np.testing.assert_array_equal(self.initial['translation_m'],[.3,-.2,.4])
        for dt in [0,-1,np.nan,np.inf]:
            with self.assertRaises(ValueError):self.flow.step(self.initial,dt)
        bad={**self.initial,'velocity':self.initial['velocity'].copy()};bad['velocity'][17]=.1
        with self.assertRaises(ValueError):self.flow.step(bad,.01)

if __name__=='__main__':unittest.main()
