"""Bounded actual-source quasi-velocity dynamics and conservation checks."""
import unittest
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_free_dynamics import ThoracicFreeDynamics, pose_rates

ROOT=Path(__file__).resolve().parents[1]

class FreeDynamicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=ThoracicMechanism(ROOT/'data/research/thoracic_mechanism/v2/manifest.json')
        cls.dynamics=ThoracicFreeDynamics(cls.model)
        cls.q=np.zeros(26);cls.q[0]=.002;cls.q[24]=.0003;cls.q[25]=.001
        cls.u=np.linspace(-.03,.04,32)
        for k in cls.model.locked:cls.u[6+k]=0
        cls.free=cls.dynamics.evaluate(cls.q,cls.u)
        eps=1e-5
        plus=cls.model.kinetic(cls.q+eps*cls.u[6:])['mass_matrix']
        minus=cls.model.kinetic(cls.q-eps*cls.u[6:])['mass_matrix']
        cls.mdot=(plus-minus)/(2*eps)

    def test_free_energy_and_world_momentum_rates(self):
        r=self.free;u=self.u;m=r['mass_matrix'];a=r['velocity_derivative']
        self.assertLess(abs(u@m@a+.5*u@self.mdot@u),2e-10)
        momentum=m@u;rate=self.mdot@u+m@a
        np.testing.assert_allclose(rate[:3]+np.cross(u[3:6],momentum[:3]),0,atol=2e-9)
        np.testing.assert_allclose(rate[3:6]+np.cross(u[3:6],momentum[3:6])+np.cross(u[:3],momentum[:3]),0,atol=2e-9)
        for k in self.model.locked:self.assertEqual(a[6+k],0)
        self.assertAlmostEqual(r['mass_kg'],20.235757397039237,places=10)
        self.assertGreater(r['active_eigenvalues'][0],0)

    def test_pressure_and_external_force_energy_balance(self):
        pressure=self.model.cavity(self.q,pressure_pa=75)
        ident=next(iter(self.model.material))
        load=self.model.point_load(ident,0,self.q,[.4,-.2,.1])['generalized_force']
        force=load.copy();force[6:]+=pressure['pressure_generalized_force']
        r=self.dynamics.solve(self.free,force)
        energy_rate=self.u@r['mass_matrix']@r['velocity_derivative']+.5*self.u@self.mdot@self.u
        self.assertAlmostEqual(energy_rate,float(self.u@force),places=9)
        self.assertAlmostEqual(r['external_power_W'],float(self.u@force),places=12)
        self.assertAlmostEqual(float(self.u[6:]@pressure['pressure_generalized_force']),75*float(self.u[6:]@pressure['volume_jacobian_m3_per_coordinate']),places=12)

    def test_material_second_directional_derivative(self):
        eps=1e-4
        for ident,m in self.model.material.items():
            ids=np.unique(np.linspace(0,len(m['reference'])-1,4,dtype=int))
            x,j=self.model.material_state(ident,self.q,ids)
            plus,jplus=self.model.material_state(ident,self.q+eps*self.u[6:],ids)
            minus,jminus=self.model.material_state(ident,self.q-eps*self.u[6:],ids)
            h=self.dynamics.directional_curvature(j,self.u[6:])
            np.testing.assert_allclose(h,np.einsum('ijk,k->ij',(jplus-jminus)/(2*eps),self.u[6:]),atol=1e-11,rtol=1e-7)

    def test_locked_load_reaction_has_zero_power(self):
        force=np.zeros(32)
        for k in self.model.locked:force[6+k]=3.5
        loaded=self.dynamics.solve(self.free,force)
        np.testing.assert_allclose(loaded['velocity_derivative'],self.free['velocity_derivative'],atol=1e-15)
        for k in self.model.locked:self.assertEqual(loaded['constraint_generalized_reaction'][6+k],-3.5)
        self.assertEqual(loaded['external_power_W'],0)
        self.assertAlmostEqual(float(self.u@loaded['constraint_generalized_reaction']),0,places=15)

    def test_world_pose_rates_and_invalid_inputs(self):
        r=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        tdot,rdot=pose_rates(r,self.u)
        np.testing.assert_allclose(tdot,r@self.u[:3])
        np.testing.assert_allclose(rdot@np.array([.2,.3,.4]),r@np.cross(self.u[3:6],[.2,.3,.4]))
        np.testing.assert_allclose(r.T@rdot+rdot.T@r,0,atol=1e-15)
        with self.assertRaises(ValueError):pose_rates(r*2,self.u)
        bad=self.u.copy();bad[17]=.1
        with self.assertRaises(ValueError):self.dynamics.evaluate(self.q,bad)
        with self.assertRaises(ValueError):self.dynamics.solve(self.free,np.full(32,np.nan))

if __name__=='__main__':unittest.main()
