#!/usr/bin/env python3
"""Light source/codec contract checks; no schema regeneration or native build."""
import sys,unittest,struct,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import prepare_tissue_burn_state_patch as burn
import prepare_native_nervous_sleep_patch as sleep

class Checks(unittest.TestCase):
    def test_fresh_legacy_seed_is_byte_local_and_explicit(self):
        raw=b'<m:Root xmlns:m="uri:/mil/tatrc/physiology/datamodel" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><m:System xsi:type="m:BioGearsTissueSystemData"><m:FatigueRunningAverage/></m:System></m:Root>'
        declaration=dict(schema='ihm.explicit-fresh-tissue-burn.v1',declaration='fresh_no_prior_burn_or_escharotomy',provenance='fixture new initialization')
        result=burn.seed_fresh_legacy(raw,declaration)
        inserted=('<m:IHMBurnHistory>'+burn.encode([0]*7,0)+'</m:IHMBurnHistory>').encode()
        self.assertEqual(result.replace(inserted,b''),raw)
        with self.assertRaises(ValueError):burn.seed_fresh_legacy(result,declaration)
        with self.assertRaises(ValueError):burn.seed_fresh_legacy(raw,{})
        with self.assertRaises(ValueError):burn.seed_fresh_legacy(raw.replace(b'<m:FatigueRunningAverage/>',b'<m:Action xsi:type="m:BurnWoundData"/>'),declaration)
    def test_exact_payload_and_finite_domain(self):
        values=[.1,.2,.3,.4,.5,7.,12345.678901234567]
        text=burn.encode(values,0b10101);v,flags=burn.decode(text)
        self.assertEqual(flags,21)
        self.assertEqual([struct.pack('>d',x) for x in values],[struct.pack('>d',x) for x in v])
        for bad in [float('nan'),float('inf'),-1]:
            with self.assertRaises(ValueError):burn.encode([bad]+values[1:],0)
        with self.assertRaises(ValueError):burn.encode([0]*5+[.5,1],0)
        with self.assertRaises(ValueError):burn.encode([1]+[0]*6,0)
        with self.assertRaises(ValueError):burn.encode([0]*7,32)
        for bad in [text+':',text.replace('V1','V2'),text.upper()]:
            with self.assertRaises(ValueError):burn.decode(bad)
    def test_fresh_ecf_reference_matches_original_capture_timing(self):
        original=burn.read_inputs()['tissue']
        self.assertIn('m_baselineECFluidVolume_mL = 0.0;',burn.function_body(original,'void Tissue::SetUp()'))
        process=burn.function_body(original,'void Tissue::CalculateCompartmentalBurn()')
        self.assertIn('if (m_baselineECFluidVolume_mL == 0.0)',process)
        self.assertIn('m_baselineECFluidVolume_mL = GetExtracellularFluidVolume().GetValue(VolumeUnit::mL);',process)
        changed=burn.function_body(burn.patch_tissue(original),'void Tissue::CalculateCompartmentalBurn()')
        self.assertIn('m_baselineECFluidVolume_mL = GetExtracellularFluidVolume().GetValue(VolumeUnit::mL);',changed)
    def test_setup_preserves_initialize_creates_invalidate_rejects(self):
        text=burn.patch_tissue(burn.read_inputs()['tissue'])
        setup=burn.function_body(text,'void Tissue::SetUp()')
        initial=burn.function_body(text,'void Tissue::Initialize()')
        invalid=burn.function_body(text,'void Tissue::Invalidate()')
        self.assertNotIn('m_baselineECFluidVolume_mL = 0.0;',setup)
        self.assertIn('ihm_burn::fresh()',initial)
        self.assertIn('quiet_NaN',invalid)
        self.assertIn('ihm_burn::validate',burn.function_body(text,'void Tissue::CalculateCompartmentalBurn()'))
    def test_sleep_schema_and_io_composition(self):
        raw=burn.read_inputs()
        left=burn.patch_schema(sleep.patch_schema(raw['schema']))
        right=sleep.patch_schema(burn.patch_schema(raw['schema']))
        self.assertEqual(left,right)
        self.assertEqual(left.count('name="IHMBurnHistory"'),1)
        self.assertEqual(left.count('name="IHMSleepState"'),1)
        for text in [burn.patch_io(sleep.patch_io(raw['io'])),sleep.patch_io(burn.patch_io(raw['io']))]:
            self.assertIn('out.IHMBurnHistory(',text);self.assertIn('out.IHMSleepState(',text)
            block=burn.function_body(text,'void BiogearsPhysiology::UnMarshall(const CDM::BioGearsTissueSystemData&')
            self.assertLess(block.index('ihm_burn::decode'),block.index('out.Invalidate'))
            self.assertLess(block.index('out.BioGearsSystem::LoadState'),block.index('out.m_compartmentSyndromeCount='))
        with self.assertRaises(ValueError):burn.patch_schema(left)
    def test_hash_guard_and_no_silent_legacy_restore(self):
        raw=burn.read_inputs();hashes={k:hashlib.sha256(v.encode()).hexdigest() for k,v in raw.items()}
        burn.prepare_texts(raw,hashes)
        with self.assertRaises(ValueError):burn.prepare_texts(dict(raw,tissue=raw['tissue']+'\n'),hashes)
        text=burn.patch_io(raw['io'])
        self.assertIn('Legacy Tissue burn history unavailable',text)
        self.assertEqual(burn.decode(burn.encode([0]*7,0)),([0.]*7,0))
if __name__=='__main__':unittest.main()
