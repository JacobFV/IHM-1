"""Coordinated isolated native probe across baseline, reader-only and prior libs."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.native_environment_evidence import freeze_native_environment,materialize_native_resources,RESOURCES
RUNTIME=ROOT/'data/runtime/physiology'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--reader',required=True);p.add_argument('--prior',required=True);p.add_argument('--run',action='store_true');p.add_argument('--probe-build',type=Path);p.add_argument('--baseline',default='whole_body_integrity_gi_absorption');p.add_argument('--state-overlay',type=Path);a=p.parse_args()
 if not a.run:raise SystemExit('Explicit coordinated compile/run slot required')
 out=Path(tempfile.mkdtemp(prefix='signed-cardiovascular-native-',dir=ROOT/'data/derived/audits'));started=time.monotonic();build=a.probe_build.resolve() if a.probe_build else out/'build';
 if not a.probe_build:build.mkdir()
 names=('native_signed_cardiovascular_probe.cpp','native_signed_vascular_prior.h','native_signed_muscle_port.h','native_coupled_engine.h','native_body_ports.h','native_tissue_ports.h','native_tissue_compression.h')
 frozen={ROOT/'scripts'/n:sha(ROOT/'scripts'/n) for n in names}
 for n in names:
  if a.probe_build:assert sha(build/n)==frozen[ROOT/'scripts'/n]
  else:shutil.copyfile(ROOT/'scripts'/n,build/n)
 frozen[Path(__file__).resolve()]=sha(__file__)
 donor=ROOT/'data/raw/physiology/biogears';enginebuild=RUNTIME/'biogears-build';libs=enginebuild/'outputs/Release/lib';sysroot=RUNTIME/'sysroot';parent=RUNTIME/'variants'/a.baseline
 command=['prlimit','--as=4294967296','--','nice','-n','10','c++','-std=c++20','-O1',str(build/names[0])]
 for path in (donor/'projects/biogears/libBiogears/include',donor/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',enginebuild/'projects/biogears/generated/Release'):command+=['-I',str(path)]
 binary=build/'probe';command+=['-L',str(parent),'-L',str(libs),f'-Wl,-rpath,{parent}:{libs}','-l:libbiogears.so.8.0.0','-lbiogears_cdm','-ldl','-o',str(binary)]
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 if not a.probe_build:
  dump(build/'command.json',command)
  with (build/'compile.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,env=env,timeout=120)
 frozen[binary]=sha(binary)
 reference=json.loads((ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json').read_text());state=Path(reference['configuration']['state_path']);assert sha(state)==reference['state_sha256']
 if a.state_overlay:
  overlay=a.state_overlay.resolve();migration=json.loads((overlay/'manifest.json').read_text());assert migration['legacy_state_sha256']==reference['state_sha256'];state=overlay/'region_metadata.xml';assert sha(state)==migration['migrated_state_sha256'];assert migration['all_nonmetadata_bytes_preserved'];frozen[overlay/'manifest.json']=sha(overlay/'manifest.json');frozen[state]=sha(state)
  assert json.loads((parent/'manifest.json').read_text())['schema']=='ihm.cardiovascular-region-io-variant.v1'
 shutil.copyfile(state,out/'initial.xml')
 data={}
 for label,variant in [('baseline',parent),('reader',RUNTIME/'variants'/a.reader),('prior',RUNTIME/'variants'/a.prior)]:
  selected=json.loads((variant/'manifest.json').read_text());assert sha(variant/'libbiogears.so.8.0.0')==selected['library_sha256'];assert sha(build/'native_signed_muscle_port.h')==selected['header_sha256'];frozen[variant/'manifest.json']=sha(variant/'manifest.json');frozen[variant/'libbiogears.so.8.0.0']=selected['library_sha256']
  if label!='baseline':
   assert sha(variant/'Cardiovascular.cpp')==selected['source_sha256'];frozen[variant/'Cardiovascular.cpp']=selected['source_sha256']
   if label=='prior':
    assert sha(variant/'native_signed_vascular_prior.h')==sha(build/'native_signed_vascular_prior.h');frozen[variant/'native_signed_vascular_prior.h']=sha(build/'native_signed_vascular_prior.h')
  runenv=dict(env,LD_LIBRARY_PATH=str(variant)+':'+str(libs)+':'+str(sysroot/'usr/lib/aarch64-linux-gnu'))
  linkage=subprocess.check_output(['ldd',str(binary)],env=runenv,text=True);assert str(variant/'libbiogears.so.8.0.0') in linkage;assert 'not found' not in linkage
  setup=out/label;setup.mkdir();(setup/'ldd.txt').write_text(linkage)
  dependencies={str(Path(x).resolve()):sha(x) for x in linkage.split() if x.startswith('/') and Path(x).is_file()}
  dump(setup/'manifest.json',{'schema':'ihm.native-session.v1','command':[str(binary)],'executable_sha256':sha(binary),'dependency_sha256':dependencies})
  for name in RESOURCES:(setup/name).symlink_to(enginebuild/'runtime'/name,target_is_directory=(enginebuild/'runtime'/name).is_dir())
  freeze_native_environment(ROOT,setup,before_start=True);resources=materialize_native_resources(ROOT,setup);data[label]={}
  for mode in ('disabled','zero','positive','negative'):
   work=setup/mode;work.mkdir();(work/'states').mkdir()
   for name in RESOURCES:(work/name).symlink_to(resources/name,target_is_directory=(resources/name).is_dir())
   run=subprocess.run(['prlimit','--as=4294967296','--','nice','-n','10',str(binary),str(out/'initial.xml'),mode],cwd=work,env=runenv,capture_output=True,text=True,timeout=20)
   (work/'stdout.log').write_text(run.stdout);(work/'stderr.log').write_text(run.stderr)
   if run.returncode:raise RuntimeError('Native probe failed: '+str(work/'stderr.log'))
   rows=[json.loads(line[4:]) for line in run.stdout.splitlines() if line.startswith('IHM\t')];assert len(rows)==3;data[label][mode]=rows
 dump(out/'frames.json',data)
 def compare(left,right):
  assert set(left)==set(right)
  errors=[]
  for key,x in left.items():
   y=right[key]
   if x is None or y is None:assert x==y,key
   else:assert abs(x-y)<=1e-12*max(1,abs(x),abs(y)),(key,x,y);errors.append(abs(x-y))
  return max(errors or [0])
 zero_error=0
 for label in ('reader','prior'):
  for mode in ('disabled','zero'):
   for base,row in zip(data['baseline'][mode],data[label][mode]):zero_error=max(zero_error,compare(base['values'],row['values']))
 for mode in ('disabled','zero','positive','negative'):
  for base,row in zip(data['baseline'][mode],data['reader'][mode]):compare(base['values'],row['values'])
 for label in ('reader','prior'):
  for mode in ('zero','positive','negative'):
   for row in data[label][mode]:assert row['reader_count']==1 and row['heat_count']==row['tissue_count']==1
 for mode,sign in [('positive',-1),('negative',1)]:
  rows=data['prior'][mode];first=rows[0];t=first['trace'];assert sign*(t['muscle_after']-t['muscle_before'])>0;assert sign*(t['other_after']-t['other_before'])>0
  for row in rows[1:]:
   t=row['trace'];assert t['delta_w']==0 and t['ratio']==1;assert t['muscle_before']==t['muscle_after'] and t['other_before']==t['other_after']
  assert first['values']['resistance.Aorta1ToMuscle1']!=data['reader'][mode][0]['values']['resistance.Aorta1ToMuscle1']
 assert all(sha(p)==digest for p,digest in frozen.items())
 report={'passed':True,'zero_mode_max_abs_difference':zero_error,'wall_s':time.monotonic()-started,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'source_sha256':{str(p.relative_to(ROOT)):v for p,v in frozen.items()},'scope':'Actual native cardiovascular reader access and isolated transferred exercise resistance prior; positive/negative causal resistance and exact local modifier release, no calibrated response or complete equilibrium claim'};dump(out/'verification.json',report);print(json.dumps({'passed':True,'output':str(out),**{k:report[k] for k in ('wall_s','zero_mode_max_abs_difference','maximum_child_rss_kib')}},indent=2))
if __name__=='__main__':main()
