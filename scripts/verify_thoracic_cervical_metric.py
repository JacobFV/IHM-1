"""Actual retained donor couplers, jets, coupled mass and mechanical power."""
import unittest,json
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_cervical_metric import CoupledThoracicMetric,CervicalKinematics,skew
from ihm.assembly.cervical_registration import reference_joint
ROOT=Path(__file__).resolve().parents[1]

class CervicalJetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.neck=CervicalKinematics(json.loads((ROOT/'data/research/cervical_registration/v2/recipe.json').read_bytes()))
    def test_affine_couplers_and_reference_frames(self):
        neck=self.neck;self.assertEqual(neck.n.shape,(24,6));self.assertEqual(np.linalg.matrix_rank(neck.n),6)
        for name,jet in neck.frames(np.zeros(6)).items():np.testing.assert_allclose(jet[0],neck.recipe['bodies'][name]['body_to_target_torso_reference'],atol=1e-12)
        q=np.linspace(-.01,.012,6);full=neck.coordinates(q)[1]
        for dep,(ind,a,b) in neck.constraints.items():self.assertAlmostEqual(full[dep],a*full[ind]+b,places=14)
        frames={'spine':np.array(neck.recipe['donor_spine_to_target_torso'])}
        for name,joint in neck.joints.items():frames[name]=frames[joint.findtext('parent_body')]@reference_joint(joint,full)
        for name,jet in neck.frames(q).items():np.testing.assert_allclose(jet[0],frames[name],atol=1e-12)
    def test_analytic_transform_derivatives(self):
        q=np.linspace(-.01,.012,6);s=np.linspace(.03,-.02,6);h=1e-4;frames=self.neck.frames(q,s)
        plus=self.neck.frames(q+h*s,s);minus=self.neck.frames(q-h*s,s)
        for name,(t,d,td,tdd) in frames.items():
            np.testing.assert_allclose(td,(plus[name][0]-minus[name][0])/(2*h),atol=2e-11)
            np.testing.assert_allclose(tdd,(plus[name][2]-minus[name][2])/(2*h),atol=2e-11)
        for index in range(6):
            delta=np.eye(6)[index]*h;plus=self.neck.frames(q+delta);minus=self.neck.frames(q-delta)
            for name in frames:np.testing.assert_allclose(frames[name][1][index],(plus[name][0]-minus[name][0])/(2*h),atol=2e-10)

class CoupledTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=CoupledThoracicMetric(ThoracicMechanism(ROOT/'data/research/thoracic_mechanism/v2/manifest.json'),ROOT/'data/research/thoracic_mechanism/native_composition_v1/plan.json')
        cls.q=np.zeros(32);cls.u=np.linspace(-.03,.04,38)
        for k in cls.model.thorax.locked:cls.u[6+k]=0
        cls.result=cls.model.evaluate(cls.q,cls.u,75)
    def test_reference_parent_mass_and_actual_pressure_power(self):
        r=self.result;mass=27.654676965260336;c=np.array([-.03,.32,0]);i=np.diag([1.520014507406313,.7788205903211117,1.4755841071015214]);a=np.c_[np.eye(3),-skew(c)];expected=mass*a.T@a;expected[3:,3:]+=i
        np.testing.assert_allclose(r['mass_matrix'][:6,:6],expected,atol=2e-12)
        self.assertAlmostEqual(r['mass_kg'],mass,places=11);self.assertEqual(len(r['active_eigenvalues']),36);self.assertGreater(r['active_eigenvalues'][0],0)
        self.assertGreater(np.linalg.norm(r['mass_matrix'][:6,32:]),.01)
        self.assertAlmostEqual(r['kinetic_energy_J'],r['direct_material_energy_J'],places=12)
        port=r['pressure_port'];self.assertAlmostEqual(float(self.u@port['generalized_force']),port['mechanical_power_W'],places=12)
        self.assertFalse(port['native_mapping_established'])
    def test_combined_energy_and_momentum_differential_balance(self):
        h=1e-5;u=self.u;q=self.q;model=self.model;r=self.result
        plus=model.evaluate(q+h*u[6:],u);minus=model.evaluate(q-h*u[6:],u);mdot=(plus['mass_matrix']-minus['mass_matrix'])/(2*h)
        m=r['mass_matrix'];a=r['velocity_derivative'];force=r['pressure_port']['generalized_force']
        self.assertLess(abs(u@m@a+.5*u@mdot@u-u@force),2e-9)
        p=m@u;pdot=mdot@u+m@a
        np.testing.assert_allclose(pdot[:3]+np.cross(u[3:6],p[:3]),0,atol=2e-9)
        np.testing.assert_allclose(pdot[3:6]+np.cross(u[3:6],p[3:6])+np.cross(u[:3],p[:3]),0,atol=2e-9)

if __name__=='__main__':unittest.main()
