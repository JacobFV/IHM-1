#!/usr/bin/env python3
"""Prepare composable isolated Tissue source/schema patches; never build/load.

A missing legacy payload is rejected. Fresh zero history belongs to Initialize,
not LoadState. Explicit initialization of old states is a separate decision.
"""
from pathlib import Path
import argparse,json,hashlib,math,struct,difflib,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/raw/physiology/biogears'
DIRECTORY=ROOT/'data/research/tissue_burn_serialization'
PATHS={'tissue':'projects/biogears/libBiogears/src/engine/Systems/Tissue.cpp','io':'projects/biogears/libBiogears/src/io/biogears/BioGearsPhysiology.cpp','schema':'share/xsd/biogears/BioGearsPhysiology.xsd'}
REGIONS=('trunk','leftArm','rightArm','leftLeg','rightLeg')
NUMERIC=tuple(f'm_{r}DeltaResistance_mmHg_s_Per_mL' for r in REGIONS)+('m_compartmentSyndromeCount','m_baselineECFluidVolume_mL')
FLAGS=tuple(f'm_{r}Escharotomy' for r in REGIONS)

def read_inputs():return {k:(SOURCE/p).read_bytes().decode() for k,p in PATHS.items()}
def digest(text):return hashlib.sha256(text.encode()).hexdigest()
def once(text,old,new):
 if text.count(old)!=1:raise ValueError('Source anchor mismatch: '+old[:90])
 return text.replace(old,new,1)
def span(text,signature):
 start=text.index(signature);begin=text.index('{',start);depth=0
 for i in range(begin,len(text)):
  if text[i]=='{':depth+=1
  elif text[i]=='}':
   depth-=1
   if depth==0:return begin,i+1
 raise ValueError('Unclosed function')
def function_body(text,signature):
 a,b=span(text,signature);return text[a:b]
def edit_function(text,signature,edit):
 a,b=span(text,signature);return text[:a]+edit(text[a:b])+text[b:]
def snapshot(obj):
 vals=', '.join(f'{obj}.{n}' for n in NUMERIC[:5]);flags=' | '.join(f'({obj}.{name}?{1<<i}u:0u)' for i,name in enumerate(FLAGS))
 return 'ihm_burn::State{{'+vals+'}, '+obj+'.'+NUMERIC[5]+', '+obj+'.'+NUMERIC[6]+', '+flags+'}'
def restore(obj,var,indent='  '):
 lines=[f'{obj}.{name}={var}.resistance_mmHg_s_per_mL[{i}];' for i,name in enumerate(NUMERIC[:5])]
 lines += [f'{obj}.{NUMERIC[5]}={var}.syndrome_count;',f'{obj}.{NUMERIC[6]}={var}.baseline_ecf_mL;']
 lines += [f'{obj}.{name}=({var}.escharotomy_flags&{1<<i}u)!=0;' for i,name in enumerate(FLAGS)]
 return '\n'.join(indent+x for x in lines)
def patch_tissue(text):
 if 'native_tissue_burn_state.h' in text:raise ValueError('Tissue burn patch already present')
 text='#include "native_tissue_burn_state.h"\n'+text
 def setup(block):
  for name in NUMERIC:block=once(block,'  '+name+' = 0.0;','  // '+name+' is initialized/restored, not reset in SetUp.')
  for name in FLAGS:block=once(block,'  '+name+' = false;','  // '+name+' is preserved across SetUp.')
  return block
 text=edit_function(text,'void Tissue::SetUp()',setup)
 text=edit_function(text,'void Tissue::Initialize()',lambda block:once(block,'  BioGearsSystem::Initialize();','  BioGearsSystem::Initialize();\n  const auto burn=ihm_burn::fresh();\n'+restore('(*this)','burn')))
 invalid='\n'.join('  '+name+'=std::numeric_limits<double>::quiet_NaN();' for name in NUMERIC)+'\n'+'\n'.join('  '+name+'=false;' for name in FLAGS)
 text=edit_function(text,'void Tissue::Invalidate()',lambda block:once(block,'  SETissueSystem::Invalidate();','  SETissueSystem::Invalidate();\n'+invalid))
 # Validate even when no active BurnWound: uninitialized hidden state must not
 # continue simply because the action is absent on the first resumed tick.
 text=edit_function(text,'void Tissue::CalculateCompartmentalBurn()',lambda block:once(block,'  if (!m_data.GetActions().GetPatientActions().HasBurnWound()) {','  ihm_burn::validate('+snapshot('(*this)')+');\n  if (!m_data.GetActions().GetPatientActions().HasBurnWound()) {'))
 return text

def patch_io(text):
 if 'native_tissue_burn_state.h' in text:raise ValueError('Tissue IO patch already present')
 text='#include "native_tissue_burn_state.h"\n'+text
 def load(block):
  block=once(block,'    out.Invalidate();','    if(!in.IHMBurnHistory().present())throw std::runtime_error("Legacy Tissue burn history unavailable; explicit fresh initialization or recovered history required");\n    const auto burn=ihm_burn::decode(in.IHMBurnHistory().get());\n    out.Invalidate();')
  return once(block,'    out.BioGearsSystem::LoadState();','    out.BioGearsSystem::LoadState();\n'+restore('out','burn','    '))
 text=edit_function(text,'void BiogearsPhysiology::UnMarshall(const CDM::BioGearsTissueSystemData&',load)
 return edit_function(text,'void BiogearsPhysiology::Marshall(const Tissue&',lambda block:once(block,'    io::Physiology::Marshall(in, out);','    const auto burn='+snapshot('in')+';\n    ihm_burn::validate(burn);\n    io::Physiology::Marshall(in, out);\n    out.IHMBurnHistory(ihm_burn::encode(burn));'))

def patch_schema(text):
 root=ET.fromstring(text);ns={'xs':'http://www.w3.org/2001/XMLSchema'}
 target=root.findall("xs:complexType[@name='BioGearsTissueSystemData']",ns)
 if len(target)!=1 or target[0].find(".//xs:element[@name='IHMBurnHistory']",ns) is not None:raise ValueError('Tissue schema target absent or already patched')
 start=text.index('  <xs:complexType name="BioGearsTissueSystemData">');end=text.index('  </xs:complexType>',start)+len('  </xs:complexType>')
 block=once(text[start:end],'        </xs:sequence>','          <xs:element name="IHMBurnHistory" type="xs:string" minOccurs="0" maxOccurs="1"/>\n        </xs:sequence>')
 result=text[:start]+block+text[end:];ET.fromstring(result);return result

def encode(values,flags):
 values=list(values)
 if len(values)!=7 or any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<0 for x in values):raise ValueError('Invalid numeric Tissue history')
 if not isinstance(flags,int) or isinstance(flags,bool) or not 0<=flags<=31 or values[5]!=math.floor(values[5]) or values[5]>2**53:raise ValueError('Invalid flags/count')
 if (any(x>0 for x in values[:6]) or flags) and values[6]<=0:raise ValueError('Historical baseline required')
 return 'IHM_BURN_V1:'+format(flags,'02x')+''.join(':'+struct.pack('>d',float(x)).hex() for x in values)
def decode(text):
 try:
  p=text.split(':')
  if len(p)!=9 or p[0]!='IHM_BURN_V1' or len(p[1])!=2 or any(len(x)!=16 for x in p[2:]):raise ValueError('Version/length')
  values=[struct.unpack('>d',bytes.fromhex(x))[0] for x in p[2:]];flags=int(p[1],16)
  if encode(values,flags)!=text:raise ValueError('Noncanonical history')
  return values,flags
 except (TypeError,struct.error,OverflowError) as error:raise ValueError('Malformed burn history') from error

def seed_fresh_legacy(raw,declaration):
 """Byte-local XML insertion; explicit new history, never historical recovery."""
 from xml.parsers import expat
 required={'schema','declaration','provenance'}
 if set(declaration)!=required or declaration['schema']!='ihm.explicit-fresh-tissue-burn.v1' or declaration['declaration']!='fresh_no_prior_burn_or_escharotomy' or not isinstance(declaration['provenance'],str) or not declaration['provenance'].strip():raise ValueError('Explicit complete fresh no-burn declaration required')
 root=ET.fromstring(raw)
 for node in root.iter():
  tokens=[node.tag.rsplit('}',1)[-1],*node.attrib.values()]
  if any('BurnWound' in t or 'Escharotomy' in t or 'CompartmentSyndrome' in t for t in tokens):raise ValueError('Active/historical burn markers prohibit fresh no-burn initialization')
 ns='uri:/mil/tatrc/physiology/datamodel';xsi='http://www.w3.org/2001/XMLSchema-instance'
 nodes=[n for n in root if n.tag=='{'+ns+'}System' and n.attrib.get('{'+xsi+'}type','').split(':')[-1]=='BioGearsTissueSystemData']
 if len(nodes)!=1 or nodes[0].find('{'+ns+'}IHMBurnHistory') is not None:raise ValueError('Expected one legacy Tissue owner; refusing overwrite')
 parser=expat.ParserCreate(namespace_separator='|');depth=0;target=None;ends=[]
 def start(name,attrs):
  nonlocal depth,target
  depth+=1
  if depth==2 and name==ns+'|System' and attrs.get(xsi+'|type','').split(':')[-1]=='BioGearsTissueSystemData':target=depth
 def end(name):
  nonlocal depth,target
  if target==depth:ends.append(parser.CurrentByteIndex);target=None
  depth-=1
 parser.StartElementHandler=start;parser.EndElementHandler=end;parser.Parse(raw,True)
 if len(ends)!=1 or raw[ends[0]:ends[0]+2]!=b'</':raise ValueError('Explicit nonempty Tissue element required')
 at=ends[0];closing=raw[at+2:raw.index(b'>',at)].strip().decode();prefix=closing.rsplit(':',1)[0]+':' if ':' in closing else ''
 payload=encode([0]*7,0);tag=prefix+'IHMBurnHistory';insertion=('<'+tag+'>'+payload+'</'+tag+'>').encode()
 result=raw[:at]+insertion+raw[at:];ET.fromstring(result)
 return result

def prepare_texts(raw,expected_hashes):
 if set(raw)!=set(PATHS) or set(expected_hashes)!=set(PATHS):raise ValueError('All three source inputs/hashes required')
 for key,text in raw.items():
  if digest(text)!=expected_hashes[key]:raise ValueError('Input source hash mismatch: '+key)
 return {'tissue':patch_tissue(raw['tissue']),'io':patch_io(raw['io']),'schema':patch_schema(raw['schema'])}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--parent-dir',type=Path);p.add_argument('--expected-hashes',type=Path);p.add_argument('--seed-fresh-legacy',type=Path);p.add_argument('--declaration-json',type=Path);a=p.parse_args()
 if a.output.exists():raise ValueError('New output path required')
 if a.seed_fresh_legacy:
  if not a.declaration_json or a.parent_dir or a.expected_hashes:raise ValueError('Explicit declaration JSON and separate NEW state output required')
  declaration=json.loads(a.declaration_json.read_text());raw_bytes=a.seed_fresh_legacy.read_bytes();result=seed_fresh_legacy(raw_bytes,declaration)
  receipt_path=a.output.with_suffix(a.output.suffix+'.initialization.json')
  if receipt_path.exists():raise ValueError('Initialization receipt already exists')
  a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_bytes(result)
  receipt_path.write_text(json.dumps(dict(kind='explicit_new_zero_burn_history_not_recovery',declaration=declaration,source_state_sha256=hashlib.sha256(raw_bytes).hexdigest(),output_state_sha256=hashlib.sha256(result).hexdigest(),native_acceptance='unverified'),indent=2)+'\n');print(a.output);return
 names={'tissue':'Tissue.cpp','io':'BioGearsPhysiology.cpp','schema':'BioGearsPhysiology.xsd'}
 if a.parent_dir:
  if not a.expected_hashes:raise ValueError('Composed parent requires explicit hashes for all three source inputs')
  raw={k:(a.parent_dir/name).read_bytes().decode() for k,name in names.items()};expected=json.loads(a.expected_hashes.read_text())['source_sha256']
 else:raw=read_inputs();expected=json.loads((DIRECTORY/'source_hashes.json').read_text())['source_sha256']
 patched=prepare_texts(raw,expected);a.output.mkdir(parents=True)
 receipt={'schema':'ihm.tissue-burn-source-repair.v1','compiled':False,'schema_generated':False,'native_roundtrip_verified':False,'source_sha256':expected,'patched_sha256':{},
  'owner_contract':'All 12 Tissue histories exact; missing legacy payload fails; fresh Initialize starts zero; no inferred restoration',
  'build_constraint':'Schema change requires matched dependent CDM/core regeneration and rebuild; never mix objects generated against sleep-only or original schema'}
 for key,name in names.items():
  (a.output/name).write_text(patched[key]);receipt['patched_sha256'][key]=digest(patched[key]);(a.output/(name+'.patch')).write_text(''.join(difflib.unified_diff(raw[key].splitlines(True),patched[key].splitlines(True),fromfile='parent/'+name,tofile=name)))
 for name in ['native_tissue_burn_state.h','native_tissue_burn_state_fixture.cpp']:
  (a.output/name).write_bytes((ROOT/'scripts'/name).read_bytes())
 receipt['helper_sha256']=hashlib.sha256((a.output/'native_tissue_burn_state.h').read_bytes()).hexdigest()
 (a.output/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(a.output)
if __name__=='__main__':main()
