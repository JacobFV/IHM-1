"""Bounded read-only API validation and retained canonical fixtures."""
from pathlib import Path
import sys,unittest,json
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.app.microvascular_patch import MicrovascularPatchService,PatchRequestError
class Checks(unittest.TestCase):
 def service(self):return MicrovascularPatchService(ROOT)
 def test_selected_default_returns_actual_geometry(self):
  r=self.service().materialize({'entity_id':'body-bp3d-FJ1442'})
  self.assertEqual(r['selection']['position_m'],[-.14,-.23,0.]);self.assertIn('fixture',r['selection']['position_basis'])
  self.assertGreater(len(r['patch']['zoom_edges']),0);self.assertTrue(all(e['radius_m']>0 for e in r['patch']['zoom_edges']))
  self.assertFalse(r['patch']['independent_store']);self.assertIsNone(r['patch']['native_volume_allocated_ml'])
  self.assertLess(len(json.dumps(r).encode()),r['limits']['max_response_bytes'])
 def test_mirrored_supported_selection(self):
  r=self.service().materialize({'entity_id':'body-bp3d-FJ1442M'})
  self.assertEqual(r['selection']['name'],'left vastus lateralis');self.assertTrue(r['patch']['containment']['whole_patch_inside_authored_surface'])
 def test_kidney_returns_curved_geometry_and_boundary_ownership(self):
  r=self.service().materialize({'entity_id':'body-bp3d-FJ3147'})
  p=r['patch'];self.assertEqual(p['native_owner'],'RightKidney');self.assertEqual(len(p['graph']['edges']),122)
  self.assertEqual(set(p['graph']['boundary_nodes']),{'arterial','venous','bowman','interstitial'})
  self.assertTrue(p['containment']['whole_patch_inside_authored_surface']);self.assertFalse(p['containment']['cortical_location_verified'])
  self.assertTrue(any(len(x['centerline_samples_m'])>100 for x in p['zoom_edges']))
  self.assertLess(len(json.dumps(r).encode()),r['limits']['max_response_bytes'])
 def test_kidney_injury_and_left_selection(self):
  r=self.service().materialize({'entity_id':'body-bp3d-FJ3145','scenario_id':'human_kidney_injury_v1'})
  self.assertEqual(r['patch']['native_owner'],'LeftKidney');self.assertTrue(r['patch']['graph']['body_registered'])
  self.assertEqual(r['patch']['graph']['peritubular_radius_conditioning']['mode'],'human_injury_area_equivalent')
 def test_kidney_rejects_muscle_conditioning(self):
  with self.assertRaises(PatchRequestError):self.service().materialize({'entity_id':'body-bp3d-FJ3147','cohort':'young_men'})
 def test_bad_queries_rejected_before_materializer(self):
  for data in [{},{'entity_id':'other'},{'entity_id':'body-bp3d-FJ1442','native_volume_ml':4},{'entity_id':'body-bp3d-FJ1442','section_mm':[20,20]},{'entity_id':'body-bp3d-FJ1442','resolution_m':1e-12},{'entity_id':'body-bp3d-FJ1442','seed':True},{'entity_id':'body-bp3d-FJ1442','position_m':[float('nan'),0,0]}]:
   with self.subTest(data=data),patch('ihm.app.microvascular_patch.materialize_patch') as materializer:
    with self.assertRaises(PatchRequestError):self.service().materialize(data)
    materializer.assert_not_called()
 def test_source_and_containment_errors_structured(self):
  for error,code in [(KeyError('schema'),'evidence_schema'),(FileNotFoundError('missing source'),'evidence_unavailable'),(ValueError('Changed or non-SI canonical geometry'),'evidence_integrity'),(ValueError('Requested position is not certified inside canonical muscle'),'patch_not_supported')]:
   with patch('ihm.app.microvascular_patch.materialize_patch',side_effect=error):
    with self.assertRaises(PatchRequestError) as ctx:self.service().materialize({'entity_id':'body-bp3d-FJ1442'})
    self.assertEqual(ctx.exception.code,code)
 def test_busy_no_wait(self):
  import ihm.app.microvascular_patch as module
  with module._QUERY_LOCK:
   with self.assertRaises(PatchRequestError) as ctx:self.service().materialize({'entity_id':'body-bp3d-FJ1442'})
   self.assertEqual(ctx.exception.status,503);self.assertEqual(ctx.exception.code,'busy')
 def test_budget_postcheck(self):
  huge={'zoom_edges':[{}]*257,'graph':{'edges':[]}}
  with patch('ihm.app.microvascular_patch.materialize_patch',return_value=huge):
   with self.assertRaises(PatchRequestError) as ctx:self.service().materialize({'entity_id':'body-bp3d-FJ1442'})
   self.assertEqual(ctx.exception.status,413)
if __name__=='__main__':unittest.main()
