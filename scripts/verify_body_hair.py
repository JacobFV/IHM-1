from pathlib import Path
import sys, unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.details import sample_hair,physical_mesh,shaft_mesh
class HairTests(unittest.TestCase):
 def test_glabrous_and_unresolved_regions_are_excluded(self):
  from ihm.assembly.hair import regional_density
  points=np.array([[.32,-.02,.03],[.1,-.83,.04],[0,.66,.1],[0,-.05,.10],[.28,.22,.03],[.1,-.5,.0]])
  density,region=regional_density(points)
  self.assertTrue(np.all(density[:4]==0));self.assertEqual(density[4],40.3);self.assertEqual(density[5],15.6)
 def test_inner_skin_shell_is_not_hair_bearing(self):
  from ihm.assembly.hair import outward_surface_mask
  c=np.array([[.28,.22,.04],[.28,.22,.039]])
  self.assertEqual(outward_surface_mask(c,np.array([[0,0,1],[0,0,-1]]),np.array(["forearm","forearm"])).tolist(),[True,False])
 def test_identity_and_deformation_attachment(self):
  from ihm.assembly.hair import attached_roots,display_indices
  v=np.array([[0,0,0],[.01,0,0],[0,.01,0]],float);f=np.array([[0,1,2]])
  h=sample_hair(v,f,[100.],seed=4);before=h['ids'].copy()
  small=display_indices(h,3);large=display_indices(h,10)
  self.assertTrue(set(small)<=set(large));self.assertTrue(np.array_equal(before,h['ids']))
  self.assertTrue(np.allclose(attached_roots(v+[1,2,3],f,h),h['roots_m']+[1,2,3]))
 def test_cross_entity_and_nearby_vertices_not_welded(self):
  v=np.array([[0,0,0],[0,0,0],[1e-12,0,0]],float)
  m=physical_mesh(v,np.empty((0,3),int),['a','b','a']);self.assertEqual(len(m.vertices),3)
 def test_shaft_uses_physical_radius(self):
  g=shaft_mesh([[0,0,0]],[[0,0,1]],[.002],[1e-5]);v=np.array(g['positions']).reshape(-1,3)
  self.assertTrue(np.allclose(np.linalg.norm(v[:,:2],axis=1),1e-5))
if __name__=='__main__':unittest.main()
