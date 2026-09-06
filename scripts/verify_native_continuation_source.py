#!/usr/bin/env python3
"""Lightweight composed ABI/owner/source and explicit initialization checks."""
import copy,hashlib,json,tempfile,unittest,xml.etree.ElementTree as ET
from pathlib import Path
import build_native_continuation_variant as build
import prepare_gi_serialization_repair as gi
import prepare_gi_inventory_presence as gi_presence
import prepare_tissue_burn_state_patch as burn
import patch_cardiovascular_region_io as circuit
from prepare_native_nervous_sleep_patch import seed_legacy
from verify_native_nervous_sleep_fixture import initializations,frozen_state
NS='{uri:/mil/tatrc/physiology/datamodel}'

def initialized_state(raw,provenance):
 if not isinstance(provenance,str) or not provenance.strip():raise ValueError('Explicit experiment provenance required')
 migrated,metadata=circuit.migrate(raw,hashlib.sha256(raw).hexdigest())
 sleep_amount,seeds=initializations(raw);sleep_seed={**seeds['fresh_rest'],'provenance':provenance+'; new rested-awake sleep history, native defaults and actual patient sleep amount'}
 sleep_state=seed_legacy(migrated,sleep_seed)
 burn_declaration={'schema':'ihm.explicit-fresh-tissue-burn.v1','declaration':'fresh_no_prior_burn_or_escharotomy','provenance':provenance+'; explicitly new zero burn/escharotomy history'}
 burn_state=burn.seed_fresh_legacy(sleep_state,burn_declaration)
 gi_declaration={'schema':'ihm.explicit-fresh-gi-inventory.v1','declaration':'fresh_no_prior_oral_drug_inventory','provenance':provenance+'; explicitly new empty oral drug inventory, not recovered from absence'}
 result=gi_presence.seed_fresh_legacy(burn_state,gi_declaration)
 receipt={'kind':'explicit_new_generic_initial_conditions_not_history_restoration','source_state_sha256':hashlib.sha256(raw).hexdigest(),'output_state_sha256':hashlib.sha256(result).hexdigest(),'sleep':sleep_seed,'burn':burn_declaration,'gi':gi_declaration,'cardiovascular_metadata':metadata,'intermediate_sha256':{'metadata':hashlib.sha256(migrated).hexdigest(),'sleep':hashlib.sha256(sleep_state).hexdigest(),'burn':hashlib.sha256(burn_state).hexdigest()}}
 return result,receipt

class CompositionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.out=build.DEFAULT;cls.r=json.loads((cls.out/'build-state.json').read_text())
 def test_complete_layout_dependency_union(self):
  build.guard(self.out,self.r);self.assertEqual(len(self.r['jobs']),38);self.assertEqual(sum(j['family']=='cdm' for j in self.r['jobs']),3)
  expected={str(p) for p in build.CORE.rglob('*.o.d') if 'engine/Systems/Gastrointestinal.h' in p.read_text()}
  self.assertEqual(len(expected),30);self.assertTrue(expected<=set(self.r['dependency_receipts']))
  self.assertTrue((self.out/'include/biogears/engine/Systems/Gastrointestinal.h').is_file())
 def test_shared_io_and_schema_preserve_all_owners(self):
  text=(self.out/'sources/src/io/biogears/BioGearsPhysiology.cpp').read_text();parent=(build.PARENT/'BioGearsPhysiology.cpp').read_text()
  for signature in ['void BiogearsPhysiology::UnMarshall(const CDM::BioGearsNervousSystemData&','void BiogearsPhysiology::Marshall(const Nervous&']:
   self.assertEqual(gi.method(text,signature),gi.method(parent,signature))
  for token in ['IHMSleepState','IHMBurnHistory','GITransitInventoryVersion','out.DrugTransitStates().clear()','Cannot save unknown GI drug inventory']:self.assertIn(token,text)
  schema=ET.parse(self.out/'xsd/biogears/BioGearsPhysiology.xsd');xs='{http://www.w3.org/2001/XMLSchema}'
  for name in ['IHMSleepState','IHMBurnHistory','GITransitInventoryVersion']:self.assertEqual(len(schema.findall('.//'+xs+'element[@name="'+name+'"]')),1)
 def test_inherited_corrections_not_replaced_with_raw_sources(self):
  for original,donor in self.r['source_donor_receipts'].items():
   relative=Path(original).relative_to(build.SOURCE/'projects/biogears/libBiogears')
   self.assertEqual(build.sha(self.out/'sources'/relative),donor['source_sha256'])
  for owner in ['Energy.cpp','Cardiovascular.cpp','Nervous.cpp','Diffusion.cpp','Environment.cpp','Renal.cpp','Saturation.cpp']:
   donor=next(d for key,d in self.r['source_donor_receipts'].items() if key.endswith('/'+owner));self.assertIn('/variants/',donor['source'])
 def test_explicit_fresh_history_and_metadata_seed(self):
  raw=frozen_state().read_bytes();result,receipt=initialized_state(raw,'Declared generic continuation fixture, source-only preparation')
  root=ET.fromstring(result)
  self.assertTrue(root.find('.//'+NS+'IHMSleepState').text.startswith('IHM_SLEEP_V1:0:'))
  self.assertEqual(root.find('.//'+NS+'GITransitInventoryVersion').text,'1');self.assertEqual(len(root.findall('.//'+NS+'CardiovascularRegion')),35)
  self.assertEqual(burn.decode(root.find('.//'+NS+'IHMBurnHistory').text),([0.]*7,0));self.assertEqual(receipt['sleep']['sleep_time_min'],480.)
  self.assertEqual(receipt['source_state_sha256'],build.sha(frozen_state()))
  with self.assertRaises(ValueError):initialized_state(raw,'')
 def test_changed_external_input_fails_before_native(self):
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'header.h';path.write_text('before');r=copy.deepcopy(self.r);r['external_build_inputs'][str(path)]=build.sha(path);path.write_text('after')
   with self.assertRaisesRegex(AssertionError,'External composition input changed'):build.guard(self.out,r)
if __name__=='__main__':unittest.main()
