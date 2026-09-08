"""Explicit isolated regional vascular experiment; consumes composed libraries only."""
from pathlib import Path
import argparse,hashlib,json,math,os,resource,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.native_environment_evidence import freeze_native_environment,materialize_native_resources,RESOURCES
from audit_regional_cardiovascular_membership import audit
RUNTIME=ROOT/'data/runtime/physiology'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def verify_frames(data):
 def compare(a,b):
  assert set(a)==set(b)
  for key,x in a.items():
   y=b[key]
   if x is None or y is None:assert x==y,key
   else:assert math.isfinite(x) and math.isfinite(y) and abs(x-y)<=1e-12*max(1,abs(x),abs(y)),(key,x,y)
 for label in ('reader','prior'):
  for mode in ('disabled','zero'):
   for a,b in zip(data['baseline'][mode],data[label][mode]):compare(a['values'],b['values'])
 for mode in ('disabled','zero','positive','negative'):
  for a,b in zip(data['baseline'][mode],data['reader'][mode]):compare(a['values'],b['values'])
 membership=audit()
 for label,modes in data.items():
  for mode,rows in modes.items():
   assert len(rows)==3
   for tick,row in enumerate(rows,1):
    values=row['values'];assert values['tissue.regional_skin.completed_native_steps']==tick
    if label!='baseline' and mode!='disabled':assert row['reader_count']==row['heat_count']==row['tissue_count']==1
    for name in membership['retained_vascular_tone_paths']:assert values['has_cardiovascular_region.'+name]==1 and values['has_resistance.'+name]==1
    for law in membership['cloned_paths']:
     for index,name in enumerate(law['regional_paths']):
      assert values['has_cardiovascular_region.'+name]==0
      assert values['has_resistance.'+name]==int(law['has_resistance_baseline'])
      assert math.isfinite(values['flow.'+name])
      if law['has_resistance_baseline']:
       expected=values['law_source.resistance.'+law['source_path']];actual=values['resistance.'+name]*(.2,.3,.5)[index]
       assert abs(actual-expected)<=1e-12*max(1,abs(expected)),name
    for key,value in values.items():
     if key.startswith('tissue.regional_skin.') and key.endswith(('ownership_residual_ug','paired_mass_residual_mg','fluid_step_residual_ml','volume_ownership_residual_ml')):assert value is not None and math.isfinite(value) and abs(value)<1e-7,(key,value)
 for mode,sign in (('positive',-1),('negative',1)):
  first=data['prior'][mode][0];comparison=data['reader'][mode][0]
  for key in ('Aorta1ToSkin1','Skin1ToSkin2'):
   assert sign*(first['values']['resistance.'+key]-comparison['values']['resistance.'+key])>0,key
   assert first['values']['flow.'+key]!=comparison['values']['flow.'+key],key
  for row in data['prior'][mode][1:]:
   trace=row['trace'];assert trace['delta_w']==0 and trace['ratio']==1
   assert trace['muscle_before']==trace['muscle_after'] and trace['other_before']==trace['other_after']
 return {'passed':True,'regional_paths_checked':27,'retained_vascular_skin_paths_checked':2,'native_sessions':12,'steps_per_session':3,'scope':'Regional causal response to explicit transferred native exercise prior; no calibration, configured-fraction acceptance, or complete physiological restart claim'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--preparation',required=True,type=Path);p.add_argument('--baseline',required=True,type=Path);p.add_argument('--reader',required=True,type=Path);p.add_argument('--prior',required=True,type=Path);p.add_argument('--state',required=True,type=Path);p.add_argument('--state-sha256',required=True);p.add_argument('--cdm-library',required=True,type=Path);p.add_argument('--generated-include',required=True,type=Path);p.add_argument('--run',action='store_true');a=p.parse_args()
 if not a.run:raise SystemExit('Explicit coordinated native slot required')
 started=time.monotonic();out=Path(tempfile.mkdtemp(prefix='regional-cardiovascular-native-',dir=ROOT/'data/derived/audits'));prepared=a.preparation.resolve();preparation=json.loads((prepared/'manifest.json').read_text())
 assert sha(a.state)==a.state_sha256
 generated=a.generated_include.resolve();generated_files=tuple(generated.rglob('*.hxx'));assert generated_files,'Generated schema headers required'
 frozen={a.state.resolve():a.state_sha256,a.cdm_library.resolve():sha(a.cdm_library),prepared/'manifest.json':sha(prepared/'manifest.json'),Path(__file__).resolve():sha(__file__)}
 frozen.update({path:sha(path) for path in generated_files})
 build=out/'build';build.mkdir()
 for name,digest in preparation['files'].items():
  assert sha(prepared/name)==digest
  if name.endswith(('.h','.inc','.cpp')):shutil.copyfile(prepared/name,build/name);frozen[prepared/name]=digest
 probe=build/'probe.cpp';text=probe.read_text();anchor='auto values=body_ports(bg);';assert text.count(anchor)==1
 text=text.replace(anchor,anchor+'for(const auto& [name,path]:bg.regional->original_paths)if(path->HasResistance())values["law_source.resistance."+name]=path->GetResistance(FlowResistanceUnit::mmHg_s_Per_mL);')
 probe.write_text(text);frozen[probe]=sha(probe)
 donor=ROOT/'data/raw/physiology/biogears';enginebuild=RUNTIME/'biogears-build';libs=enginebuild/'outputs/Release/lib';sysroot=RUNTIME/'sysroot';binary=build/'probe';variants={label:getattr(a,label).resolve() for label in ('baseline','reader','prior')};manifests={}
 for label,variant in variants.items():
  manifest=json.loads((variant/'manifest.json').read_text());assert sha(variant/'libbiogears.so.8.0.0')==manifest['library_sha256'];manifests[label]=manifest
  declared=manifest.get('paired_cdm_library_sha256',manifest.get('cdm_library_sha256'));assert declared==sha(a.cdm_library),'Explicit paired CDM identity required'
  for relative,digest in manifest.get('generated_receipts',{}).items():
   if relative.endswith('.hxx'):
    schema=generated/relative.split('generated/',1)[1];assert sha(schema)==digest,'Generated schema/CDM mismatch'
  frozen[variant/'manifest.json']=sha(variant/'manifest.json');frozen[variant/'libbiogears.so.8.0.0']=manifest['library_sha256']
 command=['prlimit','--as=4294967296','--','nice','-n','10','c++','-std=c++20','-O1',str(probe)]
 for path in (a.generated_include.resolve(),donor/'projects/biogears/libBiogears/include',donor/'projects/biogears-common/include',donor/'projects/biogears/libCDM/include',sysroot/'usr/include',sysroot/'usr/include/eigen3'):command+=['-I',str(path)]
 command += [str(variants['baseline']/'libbiogears.so.8.0.0'),str(a.cdm_library.resolve()),'-ldl','-o',str(binary)]
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1');dump(build/'command.json',command)
 with (build/'compile.log').open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
 shutil.copyfile(a.state,out/'initial.xml');data={}
 for label,variant in variants.items():
  runenv=dict(env,LD_LIBRARY_PATH=f'{variant}:{a.cdm_library.resolve().parent}:{libs}:{sysroot}/usr/lib/aarch64-linux-gnu')
  linkage=subprocess.check_output(['ldd',str(binary)],env=runenv,text=True);assert str(variant/'libbiogears.so.8.0.0') in linkage and str(a.cdm_library.resolve()) in linkage and 'not found' not in linkage
  setup=out/label;setup.mkdir();(setup/'ldd.txt').write_text(linkage)
  dependencies={str(Path(x).resolve()):sha(x) for x in linkage.split() if x.startswith('/') and Path(x).is_file()};dump(setup/'manifest.json',{'schema':'ihm.native-session.v1','command':[str(binary)],'executable_sha256':sha(binary),'dependency_sha256':dependencies})
  for name in RESOURCES:(setup/name).symlink_to(enginebuild/'runtime'/name,target_is_directory=(enginebuild/'runtime'/name).is_dir())
  freeze_native_environment(ROOT,setup,before_start=True);resources=materialize_native_resources(ROOT,setup);data[label]={}
  for mode in ('disabled','zero','positive','negative'):
   work=setup/mode;work.mkdir();(work/'states').mkdir()
   for name in RESOURCES:(work/name).symlink_to(resources/name,target_is_directory=(resources/name).is_dir())
   run=subprocess.run(['prlimit','--as=4294967296','--','nice','-n','10',str(binary),str(out/'initial.xml'),mode],cwd=work,env=runenv,capture_output=True,text=True,timeout=5);(work/'stdout.log').write_text(run.stdout);(work/'stderr.log').write_text(run.stderr);run.check_returncode();data[label][mode]=[json.loads(line[4:]) for line in run.stdout.splitlines() if line.startswith('IHM\t')]
 dump(out/'frames.json',data);report=verify_frames(data);assert all(sha(path)==digest for path,digest in frozen.items());report.update(wall_s=time.monotonic()-started,maximum_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,source_sha256={str(path):digest for path,digest in frozen.items()});dump(out/'verification.json',report);print(json.dumps({'output':str(out),**report},indent=2))
if __name__=='__main__':main()
