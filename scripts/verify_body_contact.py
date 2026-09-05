"""Regional FEM checks against analytic homogeneous deformation and constraints."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

class ContactTests(unittest.TestCase):
    def test_analytic_compression_and_independent_energy_gradient(self):
        from ihm.assembly.mechanics_backend import tetra_box,DeformableRegion
        x,t=tetra_box((.01,.01,.003),(2,2,2))
        mu,lam=2800.,25200.
        model=DeformableRegion(x,t,mu_pa=mu,lambda_pa=lam,density_kg_m3=1000.)
        stretch=.94;y=x*np.array([1.,1.,stretch]);energy,gradient=model.energy_gradient(y)
        density=mu/2*(stretch**2-1)-mu*np.log(stretch)+lam/2*np.log(stretch)**2
        self.assertAlmostEqual(energy,density*.01*.01*.003,places=14)
        top=np.flatnonzero(x[:,2]==.003)
        stress=mu*(stretch-1/stretch)+lam*np.log(stretch)/stretch
        self.assertAlmostEqual(gradient[top,2].sum(),stress*.01*.01,places=12)
        direction=np.random.default_rng(482).normal(size=x.shape);direction/=np.linalg.norm(direction)
        eps=1e-8
        derivative=(model.energy_gradient(y+eps*direction)[0]-model.energy_gradient(y-eps*direction)[0])/(2*eps)
        self.assertAlmostEqual(derivative,float(np.sum(gradient*direction)),places=8)

    def test_reference_volume_and_objectivity(self):
        from ihm.assembly.mechanics_backend import tetra_box,DeformableRegion
        x,t=tetra_box((.01,.01,.003),(2,2,2))
        model=DeformableRegion(x,t,mu_pa=20000.,lambda_pa=30000.,density_kg_m3=1000.)
        self.assertAlmostEqual(model.volumes.sum(),.01*.01*.003,places=14)
        rot=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        energy,gradient=model.energy_gradient(x@rot.T+[.1,.2,.3])
        self.assertLess(abs(energy),1e-14);self.assertLess(np.max(np.abs(gradient)),1e-10)

    def test_indentation_contact_and_release(self):
        from ihm.assembly.mechanics_backend import tetra_box,DeformableRegion
        x,t=tetra_box((.01,.01,.003),(4,4,2))
        model=DeformableRegion(x,t,mu_pa=20000.,lambda_pa=30000.,density_kg_m3=1000.)
        bottom=np.flatnonzero(x[:,2]==0)
        first=model.solve(fixed_nodes=bottom,indenter={'center_m':[.005,.005],'radius_m':.003,'height_m':.0029})
        self.assertGreater(first['indenter_reaction_n'],0)
        self.assertLess(first['maximum_penetration_m'],1e-10)
        self.assertGreater(first['minimum_jacobian'],0)
        self.assertLess(first['force_balance_residual_n'],1e-6)
        released=model.solve(fixed_nodes=bottom)
        self.assertLess(np.max(np.abs(np.array(released['positions_m'])-x)),1e-7)

    def test_invalid_solve_preserves_checkpoint(self):
        from ihm.assembly.mechanics_backend import tetra_box,DeformableRegion
        x,t=tetra_box((.01,.01,.003),(2,2,1))
        model=DeformableRegion(x,t,mu_pa=20000.,lambda_pa=30000.,density_kg_m3=1000.)
        before=model.checkpoint()
        with self.assertRaises(ValueError):model.solve(fixed_nodes=[0],indenter={'center_m':[0,0],'radius_m':-1,'height_m':0})
        np.testing.assert_array_equal(before,model.checkpoint())

if __name__=='__main__':unittest.main()
