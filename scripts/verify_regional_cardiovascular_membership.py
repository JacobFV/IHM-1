import unittest,tempfile,json
from pathlib import Path
from prepare_regional_cardiovascular_experiment import prepare
from audit_regional_cardiovascular_membership import audit,HEADER
class Tests(unittest.TestCase):
 def test_source_topology_separates_vascular_from_tissue_laws(self):
  r=audit();self.assertEqual(len(r['cloned_paths']),9);self.assertEqual(r['cloned_resistance_path_count'],6);self.assertEqual(r['cloned_direct_vascular_members'],0)
  self.assertEqual(set(r['retained_vascular_tone_paths']),{'Aorta1ToSkin1','Skin1ToSkin2'})
 def test_detached_law_owner_survives_clone_and_updates_before_solve(self):
  s=HEADER.read_text();self.assertIn('circuit.RemovePath(*path)',s);self.assertNotIn('DeleteFluidPath',s)
  self.assertIn('for(const auto& [n,source]:original_paths)copy_laws(*source,*paths[i].at(n),i);',s)
  self.assertIn('source.Get##Slot(FlowResistanceUnit::mmHg_s_Per_mL)/fractions[i]',s)
 def test_preparation_preserves_regional_object_owners(self):
  with tempfile.TemporaryDirectory() as temporary:
   out=prepare(Path(temporary)/'prepared');r=json.loads((out/'manifest.json').read_text());self.assertFalse(r['compiled'])
   for recipe in r['recipes'].values():
    self.assertEqual(recipe['unchanged_objects'],358)
    self.assertEqual({Path(row['replacement']).name for row in recipe['replaced']},{'Circuit.cpp.o','Cardiovascular.cpp.o'})
   self.assertIn('bg.install_regions();',(out/'probe.cpp').read_text())
if __name__=='__main__':unittest.main()
