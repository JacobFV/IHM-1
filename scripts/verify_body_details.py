from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
class DetailTests(unittest.TestCase):
 def test_flow_conserves_at_branch(self):
  from ihm.assembly.details import solve_network
  p=np.array([[0,0,0],[1,0,0],[2,1,0],[2,-1,0]],float);e=np.array([[0,1],[1,2],[1,3]])
  state=solve_network(p,e,np.array([.01,.01,.01]),{0:100.,2:0.,3:0.})
  self.assertLess(state['maximum_internal_residual_m3_s'],1e-15)
  self.assertAlmostEqual(state['flow_m3_s'][0],sum(state['flow_m3_s'][1:]),places=15)
  self.assertGreater(state['dissipation_w'],0)
 def test_paired_tree_has_return_and_capillaries(self):
  from ihm.assembly.details import microvascular_unit,solve_network
  p,e,r,k=microvascular_unit(np.zeros(3),.01,32,seed=9)
  self.assertEqual(sum(k==2),32);self.assertTrue(np.all(r>0))
  self.assertLessEqual(np.max(np.linalg.norm(p,axis=1)),.01)
  s=solve_network(p,e,r,{0:100.,1:0.})
  self.assertLess(s['maximum_internal_residual_m3_s'],1e-17)
  self.assertTrue(np.all(s['flow_m3_s']>0))
 def test_skin_samples_keep_surface_attachment(self):
  from ihm.assembly.details import sample_hair
  vertices=np.array([[0,0,0],[.01,0,0],[0,.01,0]],float);faces=np.array([[0,1,2]])
  h=sample_hair(vertices,faces,np.array([100.]),seed=1)
  self.assertGreater(len(h['roots_m']),20)
  self.assertLess(len(h['roots_m']),80)
  self.assertTrue(np.allclose(h['roots_m'][:,2],0))
  self.assertTrue(np.all(h['roots_m'][:,:2].sum(1)<=.01))
  self.assertTrue(np.allclose(h['barycentric'].sum(1),1))
  self.assertTrue(np.all(h['radius_m']>0))
 def test_welded_surface_ignores_attribute_seams(self):
  from ihm.assembly.details import physical_mesh
  v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float);f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
  split=v[f].reshape(-1,3);m=physical_mesh(split,np.arange(12).reshape(-1,3))
  self.assertTrue(m.is_watertight);self.assertAlmostEqual(m.volume,1/6)
if __name__=='__main__':unittest.main()
