#!/usr/bin/env python3
"""Actual native sleep owner XML/continuation fixture; --run needs heavy slot."""
from pathlib import Path
import argparse,json,math,os,resource,shutil,subprocess,tempfile,time,signal
import xml.etree.ElementTree as ET
from prepare_native_nervous_sleep_patch import ROOT,RUNTIME,SOURCE,sha,seed_legacy,validate_seed
from build_native_nervous_sleep_variant import DEFAULT,guard
NAMESPACE='{uri:/mil/tatrc/physiology/datamodel}'
FILES=['native_nervous_sleep_fixture.cpp','native_nervous_sleep_state.h']
def initializations(raw):
 root=ET.fromstring(raw);amount=root.find(NAMESPACE+'Patient/'+NAMESPACE+'SleepAmount')
 if amount is None or amount.get('unit') not in ('s','min','hr'):raise ValueError('Source patient SleepAmount with explicit native time unit required')
 sleep=float(amount.attrib['value'])*{'s':1/60,'min':1.,'hr':60.}[amount.attrib['unit']]
 if not math.isfinite(sleep) or sleep<=0:raise ValueError('Invalid source patient sleep amount')
 fresh={'schema':'ihm.explicit-native-sleep-initialization.v1','provenance':'New generic rested-awake sleep initial-condition experiment on retained nonsleep physiology. Native Initialize source defaults; sleep_time_min is the retained Patient SleepAmount. Missing legacy history is not recovered.','attention_lapses':3.,'biological_debt':0.,'reaction_time_s':.3,'tired_time_hr':0.,'wake_time_min':0.,'sleep_time_min':sleep,'sleep_state':'Awake'}
 history={**fresh,'provenance':'Synthetic declared awake tired history for native sleep-owner persistence acceptance; not patient history.','attention_lapses':math.nextafter(6.,7.),'biological_debt':math.nextafter(.2,1.),'reaction_time_s':math.nextafter(.4,1.),'tired_time_hr':math.nextafter(2.25,3.),'wake_time_min':4*sleep}
 sleeping={**history,'provenance':'Synthetic declared sleeping history for native sleep-owner persistence acceptance; not patient history.','sleep_state':'Sleeping'}
 for seed in (fresh,history,sleeping):validate_seed(seed);ET.fromstring(seed_legacy(raw,seed))
 return sleep,{'fresh_rest':fresh,'awake_history':history,'sleeping_history':sleeping}
def frozen_state():
 manifest=ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json'
 state=Path(json.loads(manifest.read_text())['configuration']['state_path'])
 assert sha(state)=='cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3'
 return state
def run_fixture():
 variant=DEFAULT;r=json.loads((variant/'build-state.json').read_text());assert r['stage']=='linked','Variant must be linked first'
 guard(variant,r);assert (ROOT/'scripts/native_nervous_sleep_state.h').read_bytes()==(variant/'native_nervous_sleep_state.h').read_bytes(),'Fixture codec differs from built owner'
 assert sha(variant/'libbiogears.so.8.0.0')==r['library_sha256'];assert sha(variant/'libbiogears_cdm.so.8.0.0')==r['cdm_library_sha256']
 source=frozen_state();raw=source.read_bytes();sleep,seeds=initializations(raw)
 out=Path(tempfile.mkdtemp(prefix='native-sleep-owner-',dir=ROOT/'data/derived/audits'));print(out,flush=True)
 for name in FILES:shutil.copyfile(ROOT/'scripts'/name,out/name)
 shutil.copyfile(source,out/'legacy.xml')
 report={'schema':'ihm.native-sleep-owner-acceptance.v1','variant':variant.name,'library_sha256':r['library_sha256'],'cdm_library_sha256':r['cdm_library_sha256'],'variant_manifest_sha256':sha(variant/'manifest.json'),'original_state_sha256':sha(source),'source_sha256':{name:sha(out/name) for name in FILES},'verifier_sha256':sha(Path(__file__)),'scope':'Native Nervous owner IO and CalculateSleepEffects with explicit clock. Exact full variant core object set links directly to expose hidden IO; paired CDM DSO. No core DSO invocation, whole-engine advance, initialization or stabilization.','cases':{}}
 (out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
 build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';binary=out/'fixture'
 command=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/FILES[0])]
 includes=build/'projects/biogears/libBiogears/CMakeFiles/libbiogears.dir/includes_CXX.rsp'
 # Same complete held include/link family as accepted native GI codec fixture; paired generated schema takes precedence.
 command+=['-I',str(variant/'generated'),'@'+str(includes)]
 report['include_response_sha256']=sha(includes)
 report['export_header_sha256']=sha(SOURCE/'projects/biogears/libCDM/include/biogears/cdm-exports.h')
 report['common_static_library_sha256']=sha(lib/'libbiogears_common_st.a')
 command +=['@'+str(variant/'objects.rsp'),str(lib/'libbiogears_common_st.a'),str(variant/'libbiogears_cdm.so.8.0.0'),'-ldl',str(RUNTIME/'sysroot/usr/lib/aarch64-linux-gnu/libxerces-c-3.2.so'),'-o',str(binary)]
 total_started=time.monotonic()
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(variant)+':'+str(lib)}
 def run(command,label,cwd,memory=1024**3,cpu=30):
  remaining=90-(time.monotonic()-total_started)
  if remaining<12:raise RuntimeError('Total fixture budget exhausted before next job')
  cpu=min(cpu,int(remaining)-10)
  def cap():resource.setrlimit(resource.RLIMIT_AS,(memory,memory));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu))
  start=time.monotonic()
  with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
   result=subprocess.Popen(['/usr/bin/time','-v','-o',str(out/(label+'.resources')),'nice','-n','10']+command,cwd=cwd,env=env,stdout=stdout,stderr=stderr,preexec_fn=cap,start_new_session=True)
   try:result.wait(timeout=cpu+10)
   except subprocess.TimeoutExpired:
    os.killpg(result.pid,signal.SIGKILL);result.wait();raise RuntimeError(f'{label} wall limit reached; process group reaped, retained {out}')
  if result.returncode:raise RuntimeError(f'{label} failed {result.returncode}; retained {out}')
  return {'wall_s':time.monotonic()-start,'command':command,'resources':(out/(label+'.resources')).read_text()}
 report['compile']=run(command,'compile',out)
 linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
 assert str(variant/'libbiogears_cdm.so.8.0.0') in linkage,'Wrong paired CDM library resolved'
 assert 'libbiogears.so.' not in linkage,'Owner fixture must use exact core objects, not another core DSO'
 report['core_object_inventory_sha256']=__import__('hashlib').sha256(json.dumps(r['object_sha256'],sort_keys=True).encode()).hexdigest()
 for name,seed in seeds.items():
  work=out/name;work.mkdir();new=work/'explicit-initialization.xml';new.write_bytes(seed_legacy(raw,seed))
  provenance={'kind':'explicit_new_sleep_initialization_not_history_restore','source_state_sha256':sha(source),'output_state_sha256':sha(new),'seed':seed}
  new.with_suffix('.xml.initialization.json').write_text(json.dumps(provenance,indent=2)+'\n')
  for entry in ['patients','substances','environments','nutrition','config','ecg','UCEDefs.conf','BioGearsConfiguration.xml']:(work/entry).symlink_to(build/'runtime'/entry)
  (work/'xsd').symlink_to(variant/'xsd')
  execution=run([str(binary),str(work),str(out/'legacy.xml'),str(new),str(work/'final.xml'),repr(sleep)],name,work,cpu=15)
  results=[json.loads(line[7:]) for line in (out/(name+'.stdout')).read_text().splitlines() if line.startswith('RESULT ')]
  assert len(results)==1 and results[0]['passed'] and results[0]['whole_engine_advances']==0 and results[0]['finite_extreme_overflow_rejected']
  report['cases'][name]={'initialization':provenance,'execution':execution,'result':results[0],'output_sha256':sha(work/'final.xml')}
 guard(variant,r);assert sha(source)==report['original_state_sha256']
 report['passed']=True;(out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':True,'receipt':str(out/'verification.json')}))
def resume_sleeping(retained):
 retained=retained.resolve();variant=DEFAULT;r=json.loads((variant/'build-state.json').read_text());guard(variant,r)
 prior=json.loads((retained/'preparation.json').read_text())
 assert prior['library_sha256']==r['library_sha256'] and prior['cdm_library_sha256']==r['cdm_library_sha256']
 assert prior['variant_manifest_sha256']==sha(variant/'manifest.json')
 for name in FILES:assert sha(retained/name)==prior['source_sha256'][name] and (retained/name).read_bytes()==(ROOT/'scripts'/name).read_bytes()
 binary=retained/'fixture';binary_digest=sha(binary)
 assert (retained/'compile.stderr').read_text()==''
 cases={}
 for name in ['fresh_rest','awake_history']:
  rows=[json.loads(line[7:]) for line in (retained/(name+'.stdout')).read_text().splitlines() if line.startswith('RESULT ')]
  assert len(rows)==1 and rows[0]['passed'] and rows[0]['native_xml_roundtrips']==13 and rows[0]['finite_extreme_overflow_rejected']
  cases[name]={'result':rows[0],'retained_directory':str(retained/name),'stdout_sha256':sha(retained/(name+'.stdout')),'final_xml_sha256':sha(retained/name/'final.xml')}
 source=frozen_state();sleep,seeds=initializations(source.read_bytes());seed=seeds['sleeping_history']
 out=Path(tempfile.mkdtemp(prefix='native-sleep-owner-completion-',dir=ROOT/'data/derived/audits'));print(out,flush=True)
 work=out/'sleeping_history';work.mkdir();new=work/'explicit-initialization.xml';new.write_bytes(seed_legacy(source.read_bytes(),seed))
 provenance={'kind':'explicit_new_sleep_initialization_not_history_restore','source_state_sha256':sha(source),'output_state_sha256':sha(new),'seed':seed,'xml_sleep_mode':'Asleep','native_enum':'SESleepState::Sleeping'}
 (work/'initialization.json').write_text(json.dumps(provenance,indent=2)+'\n')
 build=RUNTIME/'biogears-build'
 for entry in ['patients','substances','environments','nutrition','config','ecg','UCEDefs.conf','BioGearsConfiguration.xml']:(work/entry).symlink_to(build/'runtime'/entry)
 (work/'xsd').symlink_to(variant/'xsd')
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(variant)+':'+str(build/'outputs/Release/lib')}
 report={**prior,'resumed_from':str(retained),'retained_binary_sha256':binary_digest,'retained_compile_resources':(retained/'compile.resources').read_text(),'retained_failed_sleeping_stderr':(retained/'sleeping_history.stderr').read_text(),'resume_verifier_sha256':sha(Path(__file__)),'seeder_sha256':sha(ROOT/'scripts/prepare_native_nervous_sleep_patch.py'),'cases':cases,'initialization':provenance}
 command=[str(binary),str(work),str(retained/'legacy.xml'),str(new),str(work/'final.xml'),repr(sleep)]
 report['command']=command;(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
 def cap():resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));resource.setrlimit(resource.RLIMIT_CPU,(15,15));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
 with (out/'sleeping.stdout').open('w') as stdout,(out/'sleeping.stderr').open('w') as stderr:
  process=subprocess.Popen(['/usr/bin/time','-v','-o',str(out/'sleeping.resources'),'nice','-n','10']+command,cwd=work,env=env,stdout=stdout,stderr=stderr,preexec_fn=cap,start_new_session=True)
  try:process.wait(timeout=25)
  except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait();raise RuntimeError('Sleeping fixture group timed out and reaped')
 if process.returncode:raise RuntimeError(f'Sleeping fixture rejected {process.returncode}; retained {out}')
 rows=[json.loads(line[7:]) for line in (out/'sleeping.stdout').read_text().splitlines() if line.startswith('RESULT ')]
 assert len(rows)==1 and rows[0]['passed'] and rows[0]['native_xml_roundtrips']==13 and rows[0]['finite_extreme_overflow_rejected']
 assert rows[0]['sleep_owner_exact'].startswith('IHM_SLEEP_V1:1:')
 report['cases']['sleeping_history']={'result':rows[0],'directory':str(work),'final_xml_sha256':sha(work/'final.xml'),'resources':(out/'sleeping.resources').read_text()}
 guard(variant,r);assert sha(binary)==binary_digest and sha(source)==prior['original_state_sha256']
 report['passed']=True;(out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':True,'receipt':str(out/'verification.json')}))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');p.add_argument('--resume-sleeping',type=Path);a=p.parse_args()
 if a.resume_sleeping:resume_sleeping(a.resume_sleeping)
 elif a.run:run_fixture()
 else:
  sleep,seeds=initializations(frozen_state().read_bytes());assert seeds['fresh_rest']['sleep_time_min']==sleep and len(seeds)==3
  print(json.dumps({'source_preparation_passed':True,'patient_sleep_amount_min':sleep,'cases':list(seeds),'native_run':False}))
if __name__=='__main__':main()
