"""Light source/legacy topology proof before separately granted native roundtrip."""
import hashlib,json,tempfile,unittest
from pathlib import Path
from patch_cardiovascular_region_io import ROOT,IO,patch_io,topology,migrate,prepare
class Tests(unittest.TestCase):
 def test_existing_optional_field_read_write_and_stale_clear(self):
  s=patch_io(IO.read_text());self.assertIn('out.InvalidateCardiovascularRegion();',s);self.assertIn('out.CardiovascularRegion().reset();',s)
  for region in ('Muscle','Cerebral','Splanchnic','Extrasplanchnic','Myocardium'):
   self.assertIn('case CDM::enumResistancePathType::'+region,s);self.assertIn('case SEResistancePathType::'+region,s)
  with self.assertRaises(ValueError):patch_io(IO.read_text()+' ')
 def test_native_delete_recreate_is_not_silently_reclassified(self):
  assigned,deleted=topology();self.assertEqual(len(assigned),35);self.assertEqual(assigned['Aorta1ToMuscle1']['region'],'Muscle');self.assertEqual(len(deleted),6)
  for name in ('Aorta1ToBrain1','Brain1ToBrain2','Aorta1ToLeftKidney1','Aorta1ToRightKidney1','LeftKidney1ToLeftKidney2','RightKidney1ToRightKidney2'):self.assertNotIn(name,assigned)
 def test_exact_legacy_migration_changes_only_missing_structural_labels(self):
  r=json.loads((ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json').read_text());raw=Path(r['configuration']['state_path']).read_bytes();result,receipt=migrate(raw,r['state_sha256'])
  self.assertEqual(result.count(b'<CardiovascularRegion>'),35);self.assertTrue(receipt['all_nonmetadata_bytes_preserved'])
  with self.assertRaises(ValueError):migrate(raw+b' ',r['state_sha256'])
  with self.assertRaises(ValueError):migrate(result,hashlib.sha256(result).hexdigest())
 def test_isolated_prepare_retains_original_and_sha_receipts(self):
  with tempfile.TemporaryDirectory() as d:
   out=prepare(Path(d)/'out');r=json.loads((out/'manifest.json').read_text());self.assertEqual(hashlib.sha256((out/'region_metadata.xml').read_bytes()).hexdigest(),r['migrated_state_sha256'])
if __name__=='__main__':unittest.main()
