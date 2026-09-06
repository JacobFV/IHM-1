"""Light kidney registration and tortuous-centerline clipping invariants."""
from pathlib import Path
import sys,unittest,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.kidney_microvascular_embedding import _polyline_views,materialize_registered_kidney
from ihm.app.microvascular_patch import _scenario
class Checks(unittest.TestCase):
 def test_reentry_clips_separate_runs_and_preserves_arclength(self):
  p=np.array([[-2,0,0],[0,0,0],[2,0,0],[0,.5,0],[-2,.5,0]],float)
  runs=_polyline_views(p,np.array([-1,-1,-1]),np.array([1,1,1]),.25)
  self.assertEqual(len(runs),2)
  for points,t in runs:
   self.assertTrue(np.all(np.array(points)>=np.array([-1,-1,-1])));self.assertTrue(np.all(np.array(points)<=np.array([1,1,1])))
   self.assertTrue(np.all(np.diff(t)>0))
 def test_curved_physics_rigid_registration_and_view_identity(self):
  args={'entity_id':'body-bp3d-FJ3147','position_m':[-.07,.23,-.012],'resolution_m':50e-6,'seed':19,'scenario':_scenario('human_kidney_control_v1')['hydraulic_scenario']}
  p=materialize_registered_kidney(ROOT,**args);args['resolution_m']=25e-6;q=materialize_registered_kidney(ROOT,**args)
  self.assertEqual(p['patch_id'],q['patch_id']);self.assertEqual(p['graph'],q['graph'])
  g=p['graph'];pressure=g['physical_node_pressure_pa']
  for polyline,expected in zip(g['edge_polylines_m'],g['edge_length_m']):
   self.assertAlmostEqual(float(np.linalg.norm(np.diff(polyline,axis=0),axis=1).sum()),expected,places=12)
  for row in p['zoom_edges']:
   a,b=g['edges'][row['edge_index']];t=row['source_t']
   np.testing.assert_allclose(row['endpoint_pressure_pa'],[pressure[a]+x*(pressure[b]-pressure[a]) for x in t])
  self.assertGreater(p['containment']['clearance_margin_m'],0.);self.assertFalse(p['containment']['cortical_location_verified'])
 def test_uncertified_position_rejected(self):
  with self.assertRaises(ValueError):materialize_registered_kidney(ROOT,entity_id='body-bp3d-FJ3145',position_m=[0,0,0],resolution_m=50e-6,seed=1,scenario={})
if __name__=='__main__':unittest.main()
