"""Bounded migration checks; no native jobs or full geometry scans."""
from copy import deepcopy
import json
import unittest
import verify_skin_layer_support
from prepare_skin_layer_migration import migrate_copy, normalized, encode, digest, ANATOMY

class Checks(unittest.TestCase):
 def setUp(self):
  s=verify_skin_layer_support.Checks();s.setUp();self.geometry=s.raw
  ids=['body-skin-epidermis','body-skin-dermis','body-skin-hypodermis'];depth=0;entities=[]
  skin={'id':'skin','role':'skin','reference_geometry':s.ref}
  for name,t in zip(ids,[.0001,.0015,.005]):
   entities.append({'id':name,'role':'skin_layer','reference_geometry':deepcopy(s.ref),'volume_m3':2*t,'shell':{'thickness_m':t,'depth_interval_m':[depth,depth+t],'geometry_materialized':False}});depth+=t
  self.a={'entities':[skin,*entities],'assumption_ledger':[{'id':'SKIN-LAYER-PRIOR','statement':'old'}]}
  # An arbitrary fixture constraint, not the canonical body mass: migrate_copy
  # normalizes to whatever mass_allocation.target_mass_kg the payload it is
  # handed declares. The canonical profile is 70.7713 kg.
  self.m={'entities':[],'links':[],'muscles':[{'retained':'path'}],'mass_allocation':{'target_mass_kg':77.1107029},'source_files':{}}
  for e in self.a['entities']:
   self.m['entities'].append({'id':e['id'],'role':e['role'],'reference_geometry':deepcopy(s.ref),'volume_m3':e.get('volume_m3',.000002),'mass_role':'numerical_boundary_carrier' if e['role']=='skin' else 'material_partition_proxy','material':{'density':{'value':1000}},'bounds_m':{'min':[0,0,0],'max':[1,1,1]}})
  self.m['links']=[{'a':ids[0],'b':ids[1],'stiffness_n_m':100.,'parameter_basis':{'stiffness':'EA/L','damping_ratio':.25,'young_modulus_pa':100.,'effective_area_m2':1.,'effective_length_m':1.}}]
  normalized(self.m);self.araw=encode(self.a);self.m['source_files'][ANATOMY]=digest(self.araw);self.e=s.evidence;self.e['anatomy_sha256']=digest(self.araw)
 def run_candidate(self):return migrate_copy(self.araw,encode(self.m),self.geometry,encode(self.e))
 def test_candidate(self):
  before=deepcopy(self.m);a,m,r=map(json.loads,self.run_candidate());self.assertEqual(before,self.m)
  self.assertAlmostEqual(sum(e['mass_kg'] for e in m['entities']),77.1107029)
  self.assertEqual(m['muscles'],before['muscles']);self.assertFalse(r['published']);self.assertFalse(r['native_22body_mass_modified'])
  self.assertEqual(m['source_files'][ANATOMY],r['outputs']['anatomy.json'])
  for old,new in zip(self.a['entities'],a['entities']):self.assertEqual(old['reference_geometry'],new['reference_geometry'])
  for e in a['entities'][1:]:self.assertEqual(e['volume_m3'],e['shell']['thickness_m'])
  expected=deepcopy(m);normalized(expected);self.assertEqual(m,expected)
 def test_stale_receipt(self):
  self.m['source_files'][ANATOMY]='0'*64
  with self.assertRaises(ValueError):self.run_candidate()
 def test_stale_damping(self):
  self.m['links'][0]['damping_ns_m']*=2
  with self.assertRaises(ValueError):self.run_candidate()
 def test_stale_mass(self):
  self.m['entities'][1]['mass_kg']+=1
  with self.assertRaises(ValueError):self.run_candidate()
 def test_evidence_epoch(self):
  self.e['anatomy_sha256']='0'*64
  with self.assertRaises(ValueError):self.run_candidate()

if __name__=='__main__':unittest.main()
