"""Small exact-geometry/retained-surface vascular zoom fixtures."""
from pathlib import Path
import sys,unittest,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.muscle_microvascular_embedding import clip_segment_box, point_surface_certificate, materialize_patch
class Checks(unittest.TestCase):
 def test_exact_segment_box_parameters(self):
  self.assertEqual(clip_segment_box([-2,0,0],[2,0,0],[-1,-1,-1],[1,1,1]),(.25,.75))
  self.assertIsNone(clip_segment_box([-2,2,0],[2,2,0],[-1,-1,-1],[1,1,1]))
 def test_closed_tetrahedron_clearance(self):
  v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float);f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
  x=point_surface_certificate([.1,.1,.1],v,f)
  self.assertAlmostEqual(abs(x['winding_number']),1.);self.assertAlmostEqual(x['clearance_m'],.1)
  self.assertFalse(point_surface_certificate([2,2,2],v,f)['inside'])
  with self.assertRaises(ValueError):point_surface_certificate([.1,.1,.1],v,f[:-1])
 def scenario(self):return {'radius_mode':'human_quadriceps_1988_equal_area','diameter_state':'shrinkage_corrected','capillary_radius_cv':.2,'supply_radius_multiplier':2.,'return_radius_multiplier':2.4,'viscosity_pa_s':.003,'pressure_boundaries_pa':{0:4000.,1:1000.}}
 def patch(self,**kw):
  args=dict(entity_id='body-bp3d-FJ1442',position_m=[-.14,-.23,0.],resolution_m=50e-6,hydraulic_scenario=self.scenario(),seed=19);args.update(kw)
  return materialize_patch(ROOT,**args)
 def test_retained_registration_preserves_physics_and_zoom(self):
  g=self.patch();h=self.patch(resolution_m=25e-6)
  self.assertEqual(g['patch_id'],h['patch_id']);self.assertEqual(g['graph']['edges'],h['graph']['edges'])
  self.assertEqual(g['graph']['radius_m'],h['graph']['radius_m']);self.assertEqual(g['graph']['flow_solution'],h['graph']['flow_solution'])
  p=np.array(g['graph']['positions_m']);e=np.array(g['graph']['edges'])
  np.testing.assert_allclose(np.linalg.norm(p[e[:,1]]-p[e[:,0]],axis=1),g['graph']['edge_length_m'],rtol=1e-9)
  self.assertGreater(g['containment']['clearance_margin_m'],0.)
  self.assertFalse(g['independent_store']);self.assertIsNone(g['macro_boundary_correspondence'])
  self.assertEqual(len(g['zoom_edges']),len(e));self.assertTrue(all(x['radius_m']>0 for x in g['zoom_edges']))
 def test_zoom_clip_keeps_edge_id_and_pressure_interpolation(self):
  g=self.patch();p=np.array(g['graph']['positions_m']);center=p.mean(axis=0);half=40e-6
  h=self.patch(zoom_bounds_m=[(center-half).tolist(),(center+half).tolist()])
  self.assertEqual(g['patch_id'],h['patch_id']);self.assertEqual(g['graph'],h['graph'])
  self.assertGreater(len(h['zoom_edges']),0);self.assertTrue(any(e['clipped_view_only'] for e in h['zoom_edges']))
  for edge in h['zoom_edges']:
   a,b=g['graph']['edges'][edge['edge_index']];pressure=g['graph']['flow_solution']['pressure_pa'];t=edge['source_t']
   np.testing.assert_allclose(edge['endpoint_pressure_pa'],[pressure[a]+x*(pressure[b]-pressure[a]) for x in t])
 def test_outside_and_boundary_crossing_rejected(self):
  with self.assertRaises(ValueError):self.patch(position_m=[0,0,0])
  with self.assertRaises(ValueError):self.patch(extent_mm=100.)
if __name__=='__main__':unittest.main()
