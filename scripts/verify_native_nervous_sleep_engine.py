#!/usr/bin/env python3
"""Public paired sleep-DSO full-engine respiratory parity; --run needs native slot."""
from pathlib import Path
import argparse,json,os,resource,shutil,subprocess,tempfile,time,signal,xml.etree.ElementTree as ET
from prepare_native_nervous_sleep_patch import ROOT,RUNTIME,SOURCE,sha,seed_legacy
from build_native_nervous_sleep_variant import DEFAULT,guard
from verify_native_nervous_sleep_fixture import frozen_state,initializations
from verify_native_respiratory_engine import FILES,ENGINE,ENGINE_PIN,step_source,compare_state,comparator_tests
PRIOR=ROOT/'data/derived/audits/native-respiratory-engine-1s_h5ouy'
WRAPPER='native_nervous_sleep_engine_probe.cpp'
def prepare():
 prior=json.loads((PRIOR/'preparation.json').read_text());pins=prior['source_sha256']
 for name,digest in pins.items():assert sha(ROOT/name)==digest,'Previously accepted observer source changed: '+name
 assert sha(ENGINE)==ENGINE_PIN;comparator_tests()
 state=frozen_state();sleep,seeds=initializations(state.read_bytes());seed=seeds['fresh_rest'];new=seed_legacy(state.read_bytes(),seed)
 r=json.loads((DEFAULT/'build-state.json').read_text());assert r['stage']=='linked';guard(DEFAULT,r)
 for name,key in [('libbiogears.so.8.0.0','library_sha256'),('libbiogears_cdm.so.8.0.0','cdm_library_sha256')]:assert sha(DEFAULT/name)==r[key]
 return r,state,seed,new,step_source(ENGINE.read_text()),pins

def execute():
 r,state,seed,new,method,pins=prepare();variant=DEFAULT;build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
 out=Path(tempfile.mkdtemp(prefix='native-sleep-engine-',dir=ROOT/'data/derived/audits'));print(out,flush=True)
 for name in [*FILES,WRAPPER]:shutil.copyfile(ROOT/'scripts'/name,out/name)
 (out/'native_respiratory_engine_step.inc').write_text(method);shutil.copyfile(state,out/'legacy.xml');(out/'fresh-rested.xml').write_bytes(new)
 initialization={'kind':'explicit_new_generic_fresh_rested_sleep_initial_condition_not_recovered_history','source_state_sha256':sha(state),'output_state_sha256':sha(out/'fresh-rested.xml'),'seed':seed,'seeder_sha256':sha(ROOT/'scripts/prepare_native_nervous_sleep_patch.py')}
 (out/'initialization.json').write_text(json.dumps(initialization,indent=2)+'\n')
 report={'schema':'ihm.native-sleep-public-dso-acceptance.v1','variant_manifest_sha256':sha(variant/'manifest.json'),'library_sha256':r['library_sha256'],'cdm_library_sha256':r['cdm_library_sha256'],'generated_header_receipts':r['generated_receipts'],'source_sha256':{**pins,WRAPPER:sha(out/WRAPPER)},'generated_step_sha256':sha(out/'native_respiratory_engine_step.inc'),'verifier_sha256':sha(Path(__file__)),'comparator_source_sha256':sha(ROOT/'scripts/verify_native_respiratory_engine.py'),'prior_failed_full_state_receipt':str(PRIOR/'verification.json'),'initialization':initialization,'scope':'Public paired core/CDM DSOs: explicit new generic sleep initialization; original frozen observer/control code, 4 cases x8 native steps. Whole-state observer parity, not proof of whole-engine save/reload continuation.','cases':{}}
 (out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
 includes=build/'projects/biogears/libBiogears/CMakeFiles/libbiogears.dir/includes_CXX.rsp'
 binary=out/'probe';command=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-I',str(variant/'generated'),'@'+str(includes),str(out/WRAPPER),str(variant/'libbiogears.so.8.0.0'),str(variant/'libbiogears_cdm.so.8.0.0'),'-o',str(binary)]
 report['include_response_sha256']=sha(includes)
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(variant)+':'+str(lib)};started=time.monotonic()
 def run(command,label,cwd,expected=0,cpu=20):
  remaining=90-(time.monotonic()-started)
  if remaining<12:raise RuntimeError('Total acceptance budget exhausted before next process')
  cpu=min(cpu,int(remaining)-10)
  def cap():resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
  with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
   proc=subprocess.Popen(['/usr/bin/time','-v','-o',str(out/(label+'.resources')),'nice','-n','10']+command,cwd=cwd,env=env,stdout=stdout,stderr=stderr,preexec_fn=cap,start_new_session=True)
   try:proc.wait(timeout=cpu+10)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();raise RuntimeError('Acceptance process group timed out and reaped')
  if proc.returncode!=expected:raise RuntimeError(f'{label} returned {proc.returncode}, expected {expected}; retained {out}')
  return {'command':command,'returncode':proc.returncode,'resources':(out/(label+'.resources')).read_text()}
 def workdir(name):
  work=out/name;work.mkdir()
  for entry in ['patients','substances','environments','nutrition','config','ecg','UCEDefs.conf','BioGearsConfiguration.xml']:(work/entry).symlink_to(build/'runtime'/entry)
  (work/'xsd').symlink_to(variant/'xsd');return work
 def actual_mappings(label):
  rows=[line[len('LOADED_LIBRARY '):].split() for line in (out/(label+'.stdout')).read_text().splitlines() if line.startswith('LOADED_LIBRARY ')]
  expected={str(variant/name) for name in ['libbiogears.so.8.0.0','libbiogears_cdm.so.8.0.0']}
  assert {row[-1] for row in rows}==expected,'Actual process mapped different DSOs'
  for row in rows:assert int(row[4])==Path(row[-1]).stat().st_ino,'Mapped DSO inode changed'
  return sorted(expected)
 report['compile']=run(command,'compile',out,cpu=30);report['binary_sha256']=sha(binary)
 work=workdir('legacy_rejection');record=run([str(binary),'baseline',str(out/'legacy.xml'),str(work/'never-accepted.xml')],'legacy_rejection',work,expected=1)
 assert 'RESULT ' not in (out/'legacy_rejection.stdout').read_text() and not (work/'never-accepted.xml').exists()
 diagnostic=(out/'legacy_rejection.stdout').read_text()+(out/'legacy_rejection.stderr').read_text()+(work/'respiratory-engine.log').read_text()
 assert 'Legacy native sleep history is missing' in diagnostic,'Legacy rejection had a different cause'
 report['legacy_rejection']={**record,'actual_dso_paths':actual_mappings('legacy_rejection'),'terminal_process_reaped':True,'rollback_claim':False}
 frames={}
 for mode in ['baseline','observer','pulse_control','pulse_observer']:
  work=workdir(mode);record=run([str(binary),mode,str(out/'fresh-rested.xml'),str(work/'final.xml')],mode,work)
  rows=[json.loads(line[7:]) for line in (out/(mode+'.stdout')).read_text().splitlines() if line.startswith('RESULT ')]
  assert len(rows)==9 and [row['tick'] for row in rows]==list(range(9))
  for i,row in enumerate(rows):assert abs(row['time_s']-rows[0]['time_s']-.02*i)<1e-8
  snapshot=ET.parse(work/'final.xml').find('.//{uri:/mil/tatrc/physiology/datamodel}IHMSleepState');assert snapshot is not None and snapshot.text.startswith('IHM_SLEEP_V1:0:')
  frames[mode]=rows;report['cases'][mode]={**record,'actual_dso_paths':actual_mappings(mode),'final_state_sha256':sha(work/'final.xml'),'exact_sleep_snapshot':snapshot.text,'native_steps':8}
  (out/'progress.json').write_text(json.dumps(report,indent=2)+'\n')
 comparisons={}
 for control,observed in [('baseline','observer'),('pulse_control','pulse_observer')]:
  assert all(a['values']==b['values'] and a['time_s']==b['time_s'] for a,b in zip(frames[control],frames[observed])),'Observer altered native outputs'
  comparisons[control+'__'+observed]=compare_state(out/control/'final.xml',out/observed/'final.xml')
 report['state_comparisons']=comparisons
 for mode in ['observer','pulse_observer']:
  cumulative=0.
  for i,row in enumerate(frames[mode][1:],1):
   work=row['work'];assert work['external_pressure_pa']==(10 if mode=='pulse_observer' and i in (3,4) else 0)
   cumulative+=work['source_work_j'];assert cumulative==work['cumulative_source_work_j']
   assert abs(work['source_kcl_m3'])<1e-10 and abs(work['source_work_j']-work['generated_work_j']-work['external_work_j'])<1e-10
   for side in ['left','right']:assert abs(work[side+'.constitutive_residual_m3'])<1e-10
 pulse=frames['pulse_observer'];base=frames['observer'];assert pulse[2]['values']==base[2]['values']
 delta_volume=pulse[4]['values']['lung_volume_ml']-base[4]['values']['lung_volume_ml'];delta_flow=pulse[3]['values']['airway_flow_l_per_s']-base[3]['values']['airway_flow_l_per_s']
 assert abs(delta_volume)>1e-6 and abs(delta_flow)>1e-6
 assert pulse[5]['work']['external_pressure_pa']==0 and pulse[5]['work']['external_work_j']==0
 report['pulse_lung_delta_ml']=delta_volume;report['pulse_airway_flow_delta_l_s']=delta_flow;report['frames']=frames
 report['passed']=all(c['exact_state_values_equal'] for c in comparisons.values());report['total_native_advances']=32
 guard(variant,r);assert sha(ROOT/'scripts/verify_native_respiratory_engine.py')==report['comparator_source_sha256']
 assert sha(state)==initialization['source_state_sha256'] and sha(binary)==report['binary_sha256']
 (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':report['passed'],'receipt':str(out/'verification.json')}))
 if not report['passed']:raise SystemExit(1)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');a=p.parse_args()
 if a.run:execute()
 else:prepare();print('Public DSO fixture source preparation passed; no native job')
if __name__=='__main__':main()
