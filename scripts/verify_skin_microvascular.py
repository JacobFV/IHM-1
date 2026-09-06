#!/usr/bin/env python3
"""Bounded source-face, density and hydraulic regression checks."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.assembly.skin_microvascular_patch import materialize_skin_patch, default_scenario

class SkinTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.patch=materialize_skin_patch(ROOT,territory_id='left_lower_leg_anteromedial',scenario=default_scenario(),seed=3)
 def test_source_identity_and_projection(self):
  p=self.patch;reg=p['registration'];m=json.loads((ROOT/reg['territory_path']).read_text())
  self.assertIn(reg['source_triangle_id'],m['contact_eligible_triangle_ids'])
  self.assertIn(reg['source_triangle_id'],next(x for x in m['regions'] if x['id']==reg['territory_id'])['triangle_ids'])
  self.assertFalse(reg['physical_tissue_containment_verified'])
  b=np.array(p['node_surface_barycentric']);self.assertGreaterEqual(b.min(),0);np.testing.assert_allclose(b.sum(1),1)
  tri=np.array(reg['source_triangle_positions_m']);xyz=np.array(p['nodes_m']);normal=np.array(reg['outward_normal_prior'])
  np.testing.assert_allclose(xyz,b@tri-np.array(p['node_depth_m'])[:,None]*normal,atol=1e-14)
 def test_density_and_diameter(self):
  p=self.patch;c=p['constraints'];self.assertLessEqual(abs(c['loop_count']-c['target_loop_density_per_mm2']*c['patch_area_m2']*1e6),.5000001)
  caps=[e for e in p['zoom_edges'] if e['vessel_class']=='papillary_capillary_loop']
  self.assertEqual(len(caps),c['loop_count']);np.testing.assert_allclose([e['radius_m']*2e6 for e in caps],9.59)
 def test_flow_and_ownership(self):
  p=self.patch;h=p['hydraulics'];q=np.array(h['flow_m3_per_s']);press=np.array(h['pressure_pa']);edges=np.array(p['edges'])
  self.assertGreater(q.max(),0);self.assertLess(h['max_internal_residual_m3_per_s'],1e-19)
  self.assertLess(abs(sum(h['boundary_outflow_m3_per_s'].values())),1e-19)
  self.assertGreaterEqual(press.min(),1500-1e-7);self.assertLessEqual(press.max(),4500+1e-7)
  self.assertTrue(np.all(q*(press[edges[:,0]]-press[edges[:,1]])>=0))
  self.assertEqual(p['ownership']['additional_native_volume_m3'],0);self.assertIsNone(p['ownership']['native_boundary_binding']);self.assertFalse(p['lymphatics_materialized'])
 def test_reproducible_and_pressure_scaling(self):
  p=self.patch;other=materialize_skin_patch(ROOT,territory_id=p['registration']['territory_id'],scenario=default_scenario(),seed=3)
  self.assertEqual(p,other)
  s=default_scenario();s['arterial_pressure_pa']=7500
  high=materialize_skin_patch(ROOT,territory_id=p['registration']['territory_id'],scenario=s,seed=3)
  np.testing.assert_allclose(high['hydraulics']['flow_m3_per_s'],np.array(p['hydraulics']['flow_m3_per_s'])*2,rtol=1e-8,atol=1e-24)
 def test_all_territories_and_maximum_patch(self):
  m=json.loads((ROOT/self.patch['registration']['territory_path']).read_text())
  for region in m['regions']:
   p=materialize_skin_patch(ROOT,territory_id=region['id'],scenario=default_scenario(),patch_area_mm2=8)
   self.assertEqual(p['constraints']['loop_count'],56)
   self.assertEqual(len(p['edges']),172)
   self.assertLess(p['hydraulics']['max_internal_residual_m3_per_s'],1e-18)
 def test_finite_radius_stays_in_nominal_dermis(self):
  s=default_scenario();s['loop_apex_depth_m']=s['epidermal_thickness_m']+1e-7
  with self.assertRaises(ValueError):materialize_skin_patch(ROOT,territory_id='left_lower_leg_anteromedial',scenario=s)
 def test_invalid_requests(self):
  for kwargs in [dict(territory_id='unknown'),dict(territory_id='left_lower_leg_anteromedial',source_triangle_id=-1),dict(territory_id='left_lower_leg_anteromedial',patch_area_mm2=100),dict(territory_id='left_lower_leg_anteromedial',seed=True)]:
   with self.assertRaises(ValueError):materialize_skin_patch(ROOT,scenario=default_scenario(),**kwargs)
  s=default_scenario();s['loop_apex_depth_m']=0
  with self.assertRaises(ValueError):materialize_skin_patch(ROOT,territory_id='left_lower_leg_anteromedial',scenario=s)
if __name__=='__main__':unittest.main()
