#!/usr/bin/env python3
"""Light source/seed tests; --compile-codec runs only a pure C++ codec fixture."""
from pathlib import Path
import argparse,json,subprocess,tempfile,unittest,struct,math,hashlib,resource,xml.etree.ElementTree as ET
from prepare_native_nervous_sleep_patch import prepare,validate_seed,seed_legacy,FIELDS
ROOT=Path(__file__).resolve().parents[1]
SEED={'schema':'ihm.explicit-native-sleep-initialization.v1','provenance':'unit test explicit new history','attention_lapses':3.,'biological_debt':.125,'reaction_time_s':.3,'tired_time_hr':1.25,'wake_time_min':90.,'sleep_time_min':480.,'sleep_state':'Awake'}
XML=b'''<BioGearsState xmlns="uri:/mil/tatrc/physiology/datamodel" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><Patient><Name>untouched</Name></Patient><System xsi:type="BioGearsNervousSystemData"><SleepTime value="480" unit="min"/><WakeTime value="0" unit="min"/><BiologicalDebt value="0"/><ReactionTime value="0.0003" unit="s"/><AttentionLapses value="3"/><SleepState>Awake</SleepState><Other value="17"/></System></BioGearsState>'''
CODE=r'''
#include "native_nervous_sleep_state.h"
#include <cassert>
#include <iostream>
using namespace ihm_sleep;
template<class F>void reject(F f){bool hit=false;try{f();}catch(const std::exception&){hit=true;}assert(hit);}
int main(){
 auto fresh_state=fresh(480);assert(fresh_state.reaction_time_s==.3&&fresh_state.tired_time_hr==0&&fresh_state.sleep_state==0);
 auto s=fresh_state;s.attention_lapses=std::nextafter(3.,4.);s.biological_debt=.125;s.tired_time_hr=1.25;s.wake_time_min=90;
 auto restored=decode(encode(s));auto a=s.values(),b=restored.values();for(unsigned i=0;i<6;++i)assert(std::bit_cast<std::uint64_t>(a[i])==std::bit_cast<std::uint64_t>(b[i]));
 assert(public_compatible(3.,s.attention_lapses));reject([&]{require_public(3.001,s.attention_lapses);});
 reject([&]{decode("IHM_SLEEP_V2:0");});reject([&]{decode(encode(s)+"0");});
 auto noncanonical=encode(s);noncanonical[noncanonical.size()-1]='G';reject([&]{decode(noncanonical);});
 for(double invalid:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){auto x=s;x.tired_time_hr=invalid;reject([&]{encode(x);});}
 reject([&]{fresh(0);});auto unknown=s;unknown.sleep_state=2;reject([&]{encode(unknown);});
 auto uninitialized=s;uninitialized.reaction_time_s=std::numeric_limits<double>::quiet_NaN();reject([&]{encode(uninitialized);});
 // Small deterministic state evolution models only scalar persistence, not native physiology.
 auto advance=[](State& x){x.reaction_time_s=.3+.01667*x.tired_time_hr;x.attention_lapses=4.2+.9*x.tired_time_hr;x.tired_time_hr+=.02/3600;x.wake_time_min+=.02/60;};
 auto direct=s,resumed=decode(encode(s));for(int i=0;i<10;++i){advance(direct);advance(resumed);resumed=decode(encode(resumed));assert(encode(direct)==encode(resumed));}
 std::cout<<"{\"passed\":true,\"native_engine_tested\":false,\"exact_payload\":\""<<encode(s)<<"\"}\n";
}
'''
class ContractTests(unittest.TestCase):
 def test_source_lifecycle_patch(self):
  original,patched,_=prepare();n=patched['nervous'];io=patched['io']
  self.assertIn('quiet_NaN()',n);self.assertNotIn('m_ReactionTime_s / 1000.0',n)
  self.assertIn('rIntercept = 0.3;',n);self.assertIn('rSlope = 0.01667;',n)
  block=io[io.index('  void BiogearsPhysiology::UnMarshall(const CDM::BioGearsNervousSystemData&'):io.index('  void BiogearsPhysiology::Marshall(const Nervous&')]
  self.assertLess(block.index('IHMSleepState().present()'),block.index('out.Invalidate();'))
  self.assertGreater(block.index('out.m_TiredTime_hr=sleep.tired_time_hr'),block.index('out.BioGearsSystem::LoadState();'))
  self.assertEqual(patched['schema'].count('name="IHMSleepState"'),1)
  self.assertIn('native_signed_muscle_port.h',n)
  sleep=n[n.index('void Nervous::CalculateSleepEffects()'):n.index('void biogears::Nervous::UpdateSleepState()')]
  self.assertEqual(sleep.count('ihm_sleep::validate('),2)
  self.assertLess(sleep.rindex('ihm_sleep::validate('),sleep.index('GetSleepTime().SetValue(sleepTime'))
 def test_explicit_seed_required(self):
  for key in SEED:
   bad=dict(SEED);bad.pop(key)
   with self.assertRaises(ValueError):validate_seed(bad)
  for key in FIELDS:
   for value in [True,-1,float('nan'),float('inf')]:
    with self.assertRaises(ValueError):validate_seed({**SEED,key:value})
  with self.assertRaises(ValueError):validate_seed({**SEED,'provenance':' '})
 def test_seed_exact_bytes_and_preserved_other_owners(self):
  payload=validate_seed(SEED);bits=payload.split(':')[2:]
  self.assertEqual([struct.unpack('>d',bytes.fromhex(x))[0] for x in bits],[SEED[k] for k in FIELDS])
  result=ET.fromstring(seed_legacy(XML,SEED));ns='{uri:/mil/tatrc/physiology/datamodel}'
  self.assertEqual(result.find(ns+'Patient/'+ns+'Name').text,'untouched');self.assertEqual(result.find('.//'+ns+'Other').attrib,{'value':'17'})
  self.assertEqual(result.find('.//'+ns+'IHMSleepState').text,payload)
  with self.assertRaises(ValueError):seed_legacy(ET.tostring(result),SEED)
 def test_namespace_lexical_type_preserved(self):
  raw=XML.replace(b'xmlns:xsi=',b'xmlns:cdm="uri:/mil/tatrc/physiology/datamodel" xmlns:xsi=').replace(b'xsi:type="BioGears',b'xsi:type="cdm:BioGears')
  output=seed_legacy(raw,SEED);ET.fromstring(output);self.assertIn(b'xmlns:cdm=',output);self.assertIn(b'cdm:BioGearsNervousSystemData',output)
 def test_sleeping_xml_uses_native_schema_asleep(self):
  from prepare_native_nervous_sleep_patch import SOURCE
  result=ET.fromstring(seed_legacy(XML,{**SEED,'sleep_state':'Sleeping'}));ns='{uri:/mil/tatrc/physiology/datamodel}'
  self.assertEqual(result.find('.//'+ns+'SleepState').text,'Asleep')
  self.assertTrue(result.find('.//'+ns+'IHMSleepState').text.startswith('IHM_SLEEP_V1:1:'))
  mapper=(SOURCE/'projects/biogears/libBiogears/src/io/cdm/Physiology.cpp').read_text()
  self.assertIn('case CDM::enumSleepState::Asleep:\n        out = SESleepState::Sleeping;',mapper)
 def test_native_default_precision_is_not_lossless(self):
  x=math.nextafter(3.,4.);self.assertNotEqual(float(format(x,'.15g')),x)

def compile_codec():
 out=Path(tempfile.mkdtemp(prefix='native-sleep-codec-',dir=ROOT/'data/derived/audits'));source=out/'fixture.cpp';source.write_text(CODE);binary=out/'fixture'
 command=['c++','-std=c++20','-Wall','-Wextra','-Werror','-I'+str(ROOT/'scripts'),str(source),'-o',str(binary)]
 def caps():resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2));resource.setrlimit(resource.RLIMIT_CPU,(30,30))
 subprocess.run(['nice','-n','10']+command,check=True,timeout=35,preexec_fn=caps)
 result=subprocess.run([str(binary)],check=True,text=True,capture_output=True,timeout=5,preexec_fn=caps);record=json.loads(result.stdout)
 self_digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 record.update(header_sha256=self_digest(ROOT/'scripts/native_nervous_sleep_state.h'),fixture_sha256=self_digest(source),executable_sha256=self_digest(binary),command=command,scope='Pure C++ codec/state persistence only; no schema codegen, native compilation, patient or physiology advance')
 (out/'verification.json').write_text(json.dumps(record,indent=2)+'\n');print(out)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compile-codec',action='store_true');a=p.parse_args()
 result=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(ContractTests))
 if not result.wasSuccessful():raise SystemExit(1)
 if a.compile_codec:compile_codec()
