#!/usr/bin/env python3
"""Prepare isolated source/XSD patches only. Never compiles or edits held variants."""
from pathlib import Path
import argparse,difflib,hashlib,json,struct,subprocess,tempfile,math,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT/'data/raw/physiology/biogears';RUNTIME=ROOT/'data/runtime/physiology'
REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'
FIELDS=('attention_lapses','biological_debt','reaction_time_s','tired_time_hr','wake_time_min','sleep_time_min')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def once(text,old,new):
 if text.count(old)!=1:raise ValueError('Native source anchor mismatch: '+old[:80])
 return text.replace(old,new,1)
def function(text,start,end,edit):
 i=text.index(start);j=text.index(end,i);return text[:i]+edit(text[i:j])+text[j:]
def snapshot(obj):
 return f'ihm_sleep::State{{{obj}.m_AttentionLapses, {obj}.m_BiologicalDebt, {obj}.m_ReactionTime_s, {obj}.m_TiredTime_hr, {obj}.GetWakeTime(TimeUnit::min), {obj}.GetSleepTime(TimeUnit::min), {obj}.GetSleepState()==SESleepState::Sleeping?1u:0u}}'
def patch_nervous(text):
 text='#include "native_nervous_sleep_state.h"\n'+text
 text=once(text,'  SENervousSystem::Invalidate();','  SENervousSystem::Invalidate();\n  m_AttentionLapses=m_BiologicalDebt=m_ReactionTime_s=m_TiredTime_hr=std::numeric_limits<double>::quiet_NaN();')
 text=function(text,'void Nervous::SetUp()','void Nervous::AtSteadyState()',lambda s:once(s,'  m_BiologicalDebt = 0.0;','  // Stateful sleep debt is initialized or restored, never reset by SetUp.'))
 def sleep(s):
  s=once(s,'  const double rIntercept = 300.0;','  const double rIntercept = 0.3; // seconds')
  s=once(s,'  const double rSlope = 16.67;','  const double rSlope = 0.01667; // seconds per tired hour')
  s=once(s,'GetReactionTime().SetValue(m_ReactionTime_s / 1000.0, TimeUnit::s);','GetReactionTime().SetValue(m_ReactionTime_s, TimeUnit::s);')
  s=s.replace('//reaction time is computed in ms but we don\'t spport that time unit','// private and public reaction time are seconds')
  return once(s,'  SESleepState sleepState = GetSleepState();','  SESleepState sleepState = GetSleepState();\n  if(sleepState!=SESleepState::Awake&&sleepState!=SESleepState::Sleeping)throw std::runtime_error("Uninitialized sleep mode");\n  ihm_sleep::validate('+snapshot('(*this)')+');')
 return function(text,'void Nervous::CalculateSleepEffects()','void biogears::Nervous::UpdateSleepState()',sleep)
def patch_io(text):
 text='#include "native_nervous_sleep_state.h"\n'+text
 def load(s):
  s=once(s,'    out.Invalidate();','    if(!in.IHMSleepState().present())throw std::runtime_error("Legacy native sleep history is missing; explicit initialization required");\n    const auto sleep=ihm_sleep::decode(in.IHMSleepState().get());\n    out.Invalidate();')
  restore='''    out.BioGearsSystem::LoadState();
    ihm_sleep::require_public(out.GetAttentionLapses().GetValue(),sleep.attention_lapses);
    ihm_sleep::require_public(out.GetBiologicalDebt().GetValue(),sleep.biological_debt);
    ihm_sleep::require_public(out.GetReactionTime(TimeUnit::s),sleep.reaction_time_s);
    ihm_sleep::require_public(out.GetWakeTime(TimeUnit::min),sleep.wake_time_min);
    ihm_sleep::require_public(out.GetSleepTime(TimeUnit::min),sleep.sleep_time_min);
    if(out.GetSleepState()!=(sleep.sleep_state?SESleepState::Sleeping:SESleepState::Awake))throw std::runtime_error("Sleep mode conflicts with owner snapshot");
    out.m_AttentionLapses=sleep.attention_lapses;out.m_BiologicalDebt=sleep.biological_debt;
    out.m_ReactionTime_s=sleep.reaction_time_s;out.m_TiredTime_hr=sleep.tired_time_hr;
    out.GetAttentionLapses().SetValue(sleep.attention_lapses);out.GetBiologicalDebt().SetValue(sleep.biological_debt);
    out.GetReactionTime().SetValue(sleep.reaction_time_s,TimeUnit::s);
    out.GetWakeTime().SetValue(sleep.wake_time_min,TimeUnit::min);out.GetSleepTime().SetValue(sleep.sleep_time_min,TimeUnit::min);'''
  return once(s,'    out.BioGearsSystem::LoadState();',restore)
 text=function(text,'  void BiogearsPhysiology::UnMarshall(const CDM::BioGearsNervousSystemData&','  void BiogearsPhysiology::Marshall(const Nervous&',load)
 def save(s):
  checks='''    if(in.GetSleepState()!=SESleepState::Awake&&in.GetSleepState()!=SESleepState::Sleeping)throw std::runtime_error("Cannot save uninitialized sleep state");
    const auto sleep='''+snapshot('in')+''';
    ihm_sleep::validate(sleep);
    if(in.GetAttentionLapses()!=sleep.attention_lapses||in.GetBiologicalDebt()!=sleep.biological_debt||in.GetReactionTime(TimeUnit::s)!=sleep.reaction_time_s)throw std::runtime_error("Native sleep public/private owners diverged");
    io::Physiology::Marshall(in, out);
    out.IHMSleepState(ihm_sleep::encode(sleep));'''
  return once(s,'    io::Physiology::Marshall(in, out);',checks)
 return function(text,'  void BiogearsPhysiology::Marshall(const Nervous&','  // class RenalSystem',save)
def patch_schema(text):
 return function(text,'  <xs:complexType name="BioGearsNervousSystemData">','  <xs:complexType name="BioGearsRenalSystemData">',lambda s:once(s,'        </xs:sequence>','          <xs:element name="IHMSleepState" type="xs:string" minOccurs="0" maxOccurs="1"/>\n        </xs:sequence>'))
def validate_seed(seed):
 expected={'schema','provenance','sleep_state',*FIELDS}
 if set(seed)!=expected or seed['schema']!='ihm.explicit-native-sleep-initialization.v1':raise ValueError('Explicit complete sleep initialization required')
 if not isinstance(seed['provenance'],str) or not seed['provenance'].strip():raise ValueError('Sleep initialization provenance required')
 if seed['sleep_state'] not in ('Awake','Sleeping'):raise ValueError('Invalid sleep mode')
 for key in FIELDS:
  value=seed[key]
  if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:raise ValueError('Invalid explicit '+key)
 if seed['reaction_time_s']<=0 or seed['sleep_time_min']<=0:raise ValueError('Nonpositive reaction/sleep time')
 return 'IHM_SLEEP_V1:'+str(int(seed['sleep_state']=='Sleeping'))+''.join(':'+struct.pack('>d',float(seed[key])).hex() for key in FIELDS)
def seed_legacy(raw,seed):
 payload=validate_seed(seed);root=ET.fromstring(raw);ns='{uri:/mil/tatrc/physiology/datamodel}'
 nodes=[n for n in root if n.tag==ns+'System' and n.attrib.get('{http://www.w3.org/2001/XMLSchema-instance}type','').endswith('BioGearsNervousSystemData')]
 if len(nodes)!=1:raise ValueError('Expected exact native Nervous system')
 n=nodes[0]
 if n.find(ns+'IHMSleepState') is not None:raise ValueError('Refusing to reseed versioned native sleep state')
 for key,name,unit in [('attention_lapses','AttentionLapses',None),('biological_debt','BiologicalDebt',None),('reaction_time_s','ReactionTime','s'),('wake_time_min','WakeTime','min'),('sleep_time_min','SleepTime','min')]:
  field=n.find(ns+name)
  if field is None:raise ValueError('Missing legacy public sleep field '+name)
  field.set('value',repr(float(seed[key])))
  if unit:field.set('unit',unit)
 mode=n.find(ns+'SleepState')
 if mode is None:raise ValueError('Missing legacy SleepState')
 mode.text=seed['sleep_state'];ET.SubElement(n,ns+'IHMSleepState').text=payload
 # ElementTree does not track prefixes used only inside xsi:type values.
 declarations={}
 for _,(prefix,uri) in ET.iterparse(__import__('io').BytesIO(raw),events=['start-ns']):
  if prefix in declarations and declarations[prefix]!=uri:raise ValueError('Namespace rebinding requires a dedicated migration')
  declarations[prefix]=uri
  if not prefix.startswith('ns'):ET.register_namespace(prefix,uri)
 output=ET.tostring(root,encoding='utf-8',xml_declaration=True)
 required={value.split(':',1)[0] for node in root.iter() for key,value in node.attrib.items() if key=='{http://www.w3.org/2001/XMLSchema-instance}type' and ':' in value}
 from xml.sax.saxutils import quoteattr
 for prefix in sorted(required):
  if prefix not in declarations:raise ValueError('Unbound source type prefix')
  declaration=('xmlns:'+prefix+'=').encode()
  if declaration not in output:
   insertion=output.index(b'>',output.index(b'?>')+2)
   attribute=(' xmlns:'+prefix+'='+quoteattr(declarations[prefix])).encode()
   output=output[:insertion]+attribute+output[insertion:]
 ET.fromstring(output) # Refuse malformed output before any file is written.
 return output

def prepare():
 original={}
 for key,name in [('nervous','projects/biogears/libBiogears/src/engine/Systems/Nervous.cpp'),('io','projects/biogears/libBiogears/src/io/biogears/BioGearsPhysiology.cpp'),('schema','share/xsd/biogears/BioGearsPhysiology.xsd')]:
  p=SOURCE/name;raw=p.read_bytes();assert raw==subprocess.check_output(['git','-C',str(SOURCE),'show',REVISION+':'+name]);original[key]=(p,raw)
 donor=RUNTIME/'variants/whole_body_integrity_signed_muscle_v2';dm=json.loads((donor/'manifest.json').read_text());n=donor/'Nervous.cpp'
 assert sha(n)==dm['sources']['Nervous']['patched_source_sha256']
 parent=RUNTIME/'variants/whole_body_integrity_gi_absorption';pm=json.loads((parent/'manifest.json').read_text());obj=str(donor/'Nervous.cpp.o');assert pm['object_sha256'][obj]==sha(Path(obj))
 assert sha(parent/'libbiogears.so.8.0.0')==pm['library_sha256']
 original['nervous']=(n,n.read_bytes())
 patched={'nervous':patch_nervous(original['nervous'][1].decode()),'io':patch_io(original['io'][1].decode()),'schema':patch_schema(original['schema'][1].decode())}
 return original,patched,{'parent_variant':parent.name,'parent_manifest_sha256':sha(parent/'manifest.json'),'parent_library_sha256':sha(parent/'libbiogears.so.8.0.0'),'inherited_nervous_manifest_sha256':sha(donor/'manifest.json'),'source_revision':REVISION}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepare',action='store_true');p.add_argument('--seed-legacy',type=Path);p.add_argument('--seed-json',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
 if a.seed_legacy:
  if not a.seed_json or not a.output:raise ValueError('Explicit seed JSON and NEW output path required')
  if a.output.exists() or a.output.with_suffix(a.output.suffix+'.initialization.json').exists():raise ValueError('Refusing to overwrite state or initialization receipt')
  raw=a.seed_legacy.read_bytes();seed=json.loads(a.seed_json.read_text());new=seed_legacy(raw,seed)
  a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_bytes(new)
  a.output.with_suffix(a.output.suffix+'.initialization.json').write_text(json.dumps({'schema':'ihm.native-sleep-seed-receipt.v1','kind':'explicit_new_sleep_initialization_not_history_restore','source_state_sha256':hashlib.sha256(raw).hexdigest(),'output_state_sha256':sha(a.output),'seed':seed,'native_acceptance':'unverified_requires_corrected_schema_variant'},indent=2)+'\n');return
 original,patched,receipt=prepare()
 if not a.prepare:print(json.dumps({'source_ready':True,'compiled':False,**receipt}));return
 out=Path(tempfile.mkdtemp(prefix='native-sleep-state-source-',dir=ROOT/'data/derived/audits'));files={'nervous':'Nervous.cpp','io':'BioGearsPhysiology.cpp','schema':'BioGearsPhysiology.xsd'}
 receipt.update(compiled=False,scope='Source patch preparation only; requires isolated schema codegen and matching CDM+engine rebuild, not linkable with held CDM ABI',sources={})
 for key,name in files.items():
  path,raw=original[key];(out/name).write_text(patched[key]);diff=''.join(difflib.unified_diff(raw.decode().splitlines(True),patched[key].splitlines(True),fromfile=str(path),tofile=name));(out/(name+'.patch')).write_text(diff)
  receipt['sources'][key]={'inherited_source':str(path),'inherited_sha256':sha(path),'patched_sha256':sha(out/name)}
 for name in ['native_nervous_sleep_state.h','native_signed_muscle_port.h']:(out/name).write_bytes((ROOT/'scripts'/name).read_bytes())
 receipt['helper_sha256']=sha(out/'native_nervous_sleep_state.h');receipt['builder_sha256']=sha(Path(__file__))
 (out/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(out)
if __name__=='__main__':main()
