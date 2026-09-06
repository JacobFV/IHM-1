"""Bounded operator proof with actual source metric and a declared native fixture."""
import unittest
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_cervical_metric import CoupledThoracicMetric
from ihm.assembly.thoracic_native_operator import rigid_parent_operator,compose_native_operator,native_frame_to_body_twist_map
ROOT=Path(__file__).resolve().parents[1]

class OperatorTests(unittest.TestCase):
    def test_nontrivial_parent_map_bias_replacement_and_power(self):
        metric=CoupledThoracicMetric(ThoracicMechanism(ROOT/'data/research/thoracic_mechanism/v2/manifest.json'),ROOT/'data/research/thoracic_mechanism/native_composition_v1/plan.json')
        j=np.eye(6);j[0,4]=.2;j[2,3]=-.1;beta=np.array([.002,-.003,.001,.004,.002,-.002]);native_u=np.linspace(-.03,.04,6);twist=j@native_u
        old,bold=rigid_parent_operator(metric.plan['current_native_torso'],twist)
        # Declared algebraic fixture: original actual torso plus positive native
        # remainder metric. This is not claimed to be extracted from native.
        remainder=np.diag([2,3,4,.1,.2,.3]);baseline=j.T@old@j+remainder;bias=j.T@(old@beta+bold)
        q=np.zeros(32);speed=np.linspace(-.01,.02,32)
        for k in metric.thorax.locked:speed[k]=0
        result=compose_native_operator(metric,baseline,bias,j,beta,native_u,q,speed,native_model_sha256=metric.plan['target_native_model_identity']['sha256'],pressure_pa=2.)
        c=result['composite_state_evaluation'];s=result['composite_velocity_pullback'];u=np.r_[native_u,speed];expected_energy=c['direct_material_energy_J']+.5*native_u@remainder@native_u
        self.assertAlmostEqual(float(.5*u@result['mass_matrix']@u),expected_energy,places=12)
        np.testing.assert_allclose(result['inertial_bias'],s.T@(c['inertial_bias']+c['mass_matrix']@np.r_[beta,np.zeros(32)]),atol=1e-12)
        self.assertAlmostEqual(result['pressure_power_W'],c['pressure_port']['mechanical_power_W'],places=13)
        self.assertFalse(result['native_acceleration_solved']);self.assertFalse(result['native_constraints_projected'])
        with self.assertRaises(ValueError):compose_native_operator(metric,baseline,bias,j,beta,native_u,q,speed,native_model_sha256='wrong')
    def test_world_native_spatial_to_body_quasivelocity(self):
        r=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        j=np.eye(6);u=np.array([.2,-.3,.4,1.,2.,3.]);a=np.array([.03,.02,.01,.1,.2,.3])
        body,beta=native_frame_to_body_twist_map(r,j,a,u);twist=body@u
        np.testing.assert_allclose(twist,np.r_[r.T@u[3:],r.T@u[:3]])
        np.testing.assert_allclose(r@(beta[:3]+np.cross(twist[3:],twist[:3])),a[3:])
        np.testing.assert_allclose(r@beta[3:],a[:3])

if __name__=='__main__':unittest.main()
