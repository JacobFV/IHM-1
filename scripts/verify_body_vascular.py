from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.details import solve_network,microvascular_unit
class VascularTests(unittest.TestCase):
 def test_flow_power_and_heterogeneity(self):
  p,e,r,k=microvascular_unit(np.zeros(3),.004,64,10);s=solve_network(p,e,r,{0:4500.,1:1500.})
  self.assertLess(s['maximum_internal_residual_m3_s'],1e-20)
  self.assertGreater(len(np.unique(s['viscosity_pa_s'])),1)
  self.assertAlmostEqual(float(s['boundary_flow_m3_s']@s['pressure_pa'][s['boundary_nodes']])/s['dissipation_w'],1,places=10)
  self.assertTrue(np.all((s['pressure_pa']>=1500)&(s['pressure_pa']<=4500)))
 def test_invalid_geometry_fails(self):
  for e,r in [([[-1,1]],[1e-5]),([[0,1]],[-1e-5]),([[0,1]],[float('nan')])]:
   with self.assertRaises(ValueError):solve_network([[0,0,0],[1,0,0]],e,r,{0:10.,1:0.})
 def test_unanchored_component_fails(self):
  with self.assertRaises(ValueError):solve_network([[0,0,0],[1,0,0],[2,0,0]],[[0,1]],[.01],{0:10.})
 def test_single_tube_analytic(self):
  s=solve_network([[0,0,0],[.01,0,0]],[[0,1]],[5e-6],{0:100,1:0},.004)
  self.assertAlmostEqual(s['flow_m3_s'][0]/(100*np.pi*(5e-6)**4/(8*.004*.01)),1)
if __name__=='__main__':unittest.main()
