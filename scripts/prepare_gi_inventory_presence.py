#!/usr/bin/env python3
"""Composable GI owner/codec guards; schema/build owned by composition lead.

Required schema field: optional xs:unsignedInt GITransitInventoryVersion.
No raw source mutation. Explicit legacy initialization is not historical recovery.
"""
from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
from xml.parsers import expat
from prepare_gi_serialization_repair import method

MARKER='GITransitInventoryVersion'
MEMBER='m_IHMDrugTransitInventoryInitialized'

def once(text,old,new):
 if text.count(old)!=1:raise ValueError('Expected exactly one source anchor: '+old)
 return text.replace(old,new,1)

def edit(text,signature,transform):
 old=method(text,signature)
 return once(text,old,transform(old))

def patch_header(text):
 if MEMBER in text:raise ValueError('GI owner marker already present')
 return once(text,'  // Serializable member variables (Set in Initialize and in schema)',
             '  // Serializable member variables (Set in Initialize and in schema)\n  bool '+MEMBER+' = false;')

def patch_gi(text):
 if MEMBER in text:raise ValueError('GI owner guard already present')
 text=edit(text,'void Gastrointestinal::Invalidate()',lambda b:once(b,'  SEGastrointestinalSystem::Invalidate();',
       '  '+MEMBER+' = false;\n  SEGastrointestinalSystem::Invalidate();'))
 def initialize(b):
  b=once(b,'  BioGearsSystem::Initialize();',
    '  if (!m_DrugTransitStates.empty()) throw CommonDataModelException("Fresh GI initialization requires empty drug inventory");\n  BioGearsSystem::Initialize();')
  return b[:-1]+'  '+MEMBER+' = true; // Native fresh initialization completed; no prior oral inventory.\n}'
 return edit(text,'void Gastrointestinal::Initialize()',initialize)

def patch_io(text):
 if MARKER in text:raise ValueError('GI IO presence guard already present')
 # Must run AFTER accepted three-TU serialization repair, never standalone.
 writer=method(text,'void BiogearsPhysiology::Marshall(const Gastrointestinal& in, CDM::BioGearsGastrointestinalSystemData& out)')
 if 'out.DrugTransitStates().clear();' not in writer:raise ValueError('Accepted GI inventory writer repair required first')
 def load(b):
  b=once(b,'    out.Invalidate();',
    '    if (!in.GITransitInventoryVersion().present() || in.GITransitInventoryVersion().get() != 1u)\n      throw CommonDataModelException("GI drug inventory history missing or unsupported; explicit fresh declaration or recovered authoritative state required");\n    out.Invalidate();')
  return b[:-1]+'  out.'+MEMBER+' = true; // All declared records decoded successfully.\n  }'
 text=edit(text,'void BiogearsPhysiology::UnMarshall(const CDM::BioGearsGastrointestinalSystemData& in, const SESubstanceManager& substances, Gastrointestinal& out)',load)
 def save(b):
  b=once(b,'    io::Physiology::Marshall(in, out);',
    '    if (!in.'+MEMBER+') throw CommonDataModelException("Cannot save unknown GI drug inventory");\n    io::Physiology::Marshall(in, out);')
  return b[:-1]+'  out.GITransitInventoryVersion(1u);\n  }'
 return edit(text,'void BiogearsPhysiology::Marshall(const Gastrointestinal& in, CDM::BioGearsGastrointestinalSystemData& out)',save)

def seed_fresh_legacy(raw,declaration):
 """Add marker only after explicit fresh-no-prior-inventory declaration.

Reject any transit record or oral action, even an empty/zero-looking one. Neither
its absence nor this check proves history: the caller supplies that declaration.
Preserve all existing bytes, including scalar precision, via one local insertion.
 """
 if (set(declaration)!={'schema','declaration','provenance'} or
     declaration.get('schema')!='ihm.explicit-fresh-gi-inventory.v1' or
     declaration.get('declaration')!='fresh_no_prior_oral_drug_inventory' or
     not isinstance(declaration.get('provenance'),str) or not declaration['provenance'].strip()):
  raise ValueError('Explicit fresh GI inventory declaration with provenance required')
 root=ET.fromstring(raw);ns='uri:/mil/tatrc/physiology/datamodel';xsi='http://www.w3.org/2001/XMLSchema-instance'
 for node in root.iter():
  tokens=[node.tag.rsplit('}',1)[-1],*node.attrib.values()]
  if any('DrugTransitState' in t or 'SubstanceOralDose' in t for t in tokens):
   raise ValueError('Existing drug transit record or oral action prohibits fresh declaration')
 nodes=[n for n in root if n.tag=='{'+ns+'}System' and n.attrib.get('{'+xsi+'}type','').split(':')[-1]=='BioGearsGastrointestinalSystemData']
 if len(nodes)!=1 or nodes[0].find('{'+ns+'}'+MARKER) is not None:
  raise ValueError('Exactly one legacy GI owner required; no marker overwrite')
 parser=expat.ParserCreate(namespace_separator='|');depth=0;target=None;ends=[]
 def start(name,attrs):
  nonlocal depth,target
  depth+=1
  if depth==2 and name==ns+'|System' and attrs.get(xsi+'|type','').split(':')[-1]=='BioGearsGastrointestinalSystemData':target=depth
 def end(name):
  nonlocal depth,target
  if target==depth:ends.append(parser.CurrentByteIndex);target=None
  depth-=1
 parser.StartElementHandler=start;parser.EndElementHandler=end;parser.Parse(raw,True)
 if len(ends)!=1 or raw[ends[0]:ends[0]+2]!=b'</':raise ValueError('Explicit nonempty GI element required')
 at=ends[0];closing=raw[at+2:raw.index(b'>',at)].strip().decode();prefix=closing.rsplit(':',1)[0]+':' if ':' in closing else ''
 tag=prefix+MARKER;result=raw[:at]+('<'+tag+'>1</'+tag+'>').encode()+raw[at:]
 ET.fromstring(result);return result

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--seed-fresh-legacy',type=Path,required=True);p.add_argument('--declaration-json',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 receipt_path=a.output.with_suffix(a.output.suffix+'.initialization.json')
 if a.output.exists() or receipt_path.exists():raise ValueError('New output/receipt paths required')
 raw=a.seed_fresh_legacy.read_bytes();declaration=json.loads(a.declaration_json.read_text());result=seed_fresh_legacy(raw,declaration)
 receipt={'schema':'ihm.explicit-fresh-gi-inventory-receipt.v1','kind':'fresh_empty_inventory_not_historical_recovery','declaration':declaration,'source_state_sha256':hashlib.sha256(raw).hexdigest(),'output_state_sha256':hashlib.sha256(result).hexdigest(),'native_acceptance':'unverified'}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_bytes(result);receipt_path.write_text(json.dumps(receipt,indent=2)+'\n');print(a.output)
if __name__=='__main__':main()
