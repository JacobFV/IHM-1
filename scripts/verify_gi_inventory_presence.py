#!/usr/bin/env python3
"""Light source-transform and legacy-declaration rejection tests, not native proof."""
import tempfile
from pathlib import Path
import prepare_gi_inventory_presence as p
import prepare_gi_serialization_repair as repair

def rejects(fn):
 try:fn()
 except ValueError:return
 raise AssertionError('Expected rejection')

def main():
 with tempfile.TemporaryDirectory() as d:
  repair.prepare(Path(d))
  io=p.patch_io((Path(d)/repair.FILES[0]).read_text())
  assert io.index('if (!in.GITransitInventoryVersion().present()')<io.index('out.m_IHMDrugTransitInventoryInitialized = true;')
  assert 'Cannot save unknown GI drug inventory' in io
  assert 'out.GITransitInventoryVersion(1u)' in io
  rejects(lambda:p.patch_io(io))
 raw=(repair.BASE/'src/engine/Systems/Gastrointestinal.cpp').read_text()
 gi=p.patch_gi(raw)
 assert gi.count('m_IHMDrugTransitInventoryInitialized = false;')==1
 assert gi.count('m_IHMDrugTransitInventoryInitialized = true;')==1
 assert 'Fresh GI initialization requires empty drug inventory' in gi
 header=p.patch_header((repair.BASE/'include/biogears/engine/Systems/Gastrointestinal.h').read_text())
 assert 'bool m_IHMDrugTransitInventoryInitialized = false;' in header
 rejects(lambda:p.patch_header(header))
 ns='uri:/mil/tatrc/physiology/datamodel'
 raw=('<BioGearsState xmlns="'+ns+'" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><System xsi:type="BioGearsGastrointestinalSystemData"><ChymeAbsorptionRate value="0" unit="mL/s"/></System></BioGearsState>').encode()
 declaration={'schema':'ihm.explicit-fresh-gi-inventory.v1','declaration':'fresh_no_prior_oral_drug_inventory','provenance':'synthetic source test'}
 out=p.seed_fresh_legacy(raw,declaration)
 assert out.replace(b'<GITransitInventoryVersion>1</GITransitInventoryVersion>',b'')==raw
 prefixed=raw.replace(b'<BioGearsState xmlns=',b'<d:BioGearsState xmlns:d=').replace(b'</BioGearsState>',b'</d:BioGearsState>').replace(b'<System ',b'<d:System ').replace(b'</System>',b'</d:System>').replace(b'<ChymeAbsorptionRate ',b'<d:ChymeAbsorptionRate ')
 prefixed_out=p.seed_fresh_legacy(prefixed,declaration)
 assert prefixed_out.replace(b'<d:GITransitInventoryVersion>1</d:GITransitInventoryVersion>',b'')==prefixed
 rejects(lambda:p.seed_fresh_legacy(out,declaration))
 rejects(lambda:p.seed_fresh_legacy(raw,{}))
 rejects(lambda:p.seed_fresh_legacy(raw,{**declaration,'provenance':' '}))
 for token in [b'<DrugTransitStates/>',b'<SubstanceOralDose/>']:
  rejects(lambda token=token:p.seed_fresh_legacy(raw.replace(b'</System>',token+b'</System>'),declaration))
 print('PASS source transforms; explicit fresh declaration; byte-local insertion; legacy/drug/duplicate rejection (no native execution)')
if __name__=='__main__':main()
