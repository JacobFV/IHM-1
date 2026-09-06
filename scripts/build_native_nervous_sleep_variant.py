#!/usr/bin/env python3
"""Explicit staged isolated native sleep variant: prepare/codegen/compile/link.

Each compiler/codegen/link invocation requires the coordinated heavy slot.
No stage edits original sources, libraries, object files or generated bindings.
"""
from pathlib import Path
import argparse,hashlib,json,os,resource,shlex,shutil,subprocess,time,signal
from prepare_native_nervous_sleep_patch import ROOT,RUNTIME,SOURCE,prepare,sha
BUILD=RUNTIME/'biogears-build';CORE=BUILD/'projects/biogears/libBiogears';CDM=BUILD/'projects/biogears/libCDM'
PARENT=RUNTIME/'variants/whole_body_integrity_gi_absorption'
DEFAULT=RUNTIME/'variants/whole_body_integrity_nervous_sleep_v1'
def record(out,r):
 temporary=out/'build-state.json.tmp';temporary.write_text(json.dumps(r,indent=2)+'\n');temporary.replace(out/'build-state.json')
def prepare_build(out):
 if out.exists():raise ValueError('New build directory required')
 original,patched,receipt=prepare();pm=json.loads((PARENT/'manifest.json').read_text())
 core_objects=shlex.split((PARENT/'objects.rsp').read_text());assert set(core_objects)==set(pm['object_sha256'])
 for name in core_objects:assert sha(CORE/name)==pm['object_sha256'][name]
 cdm_objects=shlex.split((CDM/'CMakeFiles/libbiogears_cdm.dir/objects1.rsp').read_text())
 inventory={'core':{str((CORE/name).resolve()):sha(CORE/name) for name in core_objects},'cdm':{str((CDM/name).resolve()):sha(CDM/name) for name in cdm_objects}}
 out.mkdir(parents=True);shutil.copytree(BUILD/'projects/biogears/generated/Release',out/'generated');shutil.copytree(SOURCE/'share/xsd',out/'xsd')
 (out/'xsd/biogears/BioGearsPhysiology.xsd').write_text(patched['schema'])
 (out/'Nervous.cpp').write_text(patched['nervous']);(out/'BioGearsPhysiology.cpp').write_text(patched['io'])
 for name in ['native_nervous_sleep_state.h','native_signed_muscle_port.h']:(out/name).write_bytes((ROOT/'scripts'/name).read_bytes())
 jobs=[];dependencies={}
 for family,cwd,target in [('cdm',CDM,'libbiogears_cdm'),('core',CORE,'libbiogears')]:
  objects=list(cwd.rglob('*.o'));deps=list(cwd.rglob('*.o.d'));assert len(objects)==len(deps)
  found=sorted(p for p in deps if 'BioGearsPhysiology.hxx' in p.read_text())
  dependencies[family]={str(p):sha(p) for p in found}
  for dep in found:
   obj=Path(str(dep)[:-2]);relative=obj.relative_to(cwd/f'CMakeFiles/{target}.dir')
   if family=='cdm':
    source=out/'generated'/Path(str(relative).removeprefix('__/generated/Release/').removesuffix('.o'))
   else:
    source=SOURCE/'projects/biogears/libBiogears'/str(relative).removesuffix('.o')
    if source.name=='BioGearsPhysiology.cpp':source=out/source.name
   assert source.is_file();jobs.append({'family':family,'old_object':str(obj.resolve()),'source':str(source),'object':str(out/'objects'/f'{len(jobs):02d}.o')})
 old=next(name for name in inventory['core'] if name.endswith('/whole_body_integrity_signed_muscle_v2/Nervous.cpp.o'))
 jobs.append({'family':'core','old_object':old,'source':str(out/'Nervous.cpp'),'object':str(out/'objects'/f'{len(jobs):02d}.o')})
 assert len(jobs)==13 and sum(j['family']=='cdm' for j in jobs)==3
 (out/'objects').mkdir()
 frozen={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}
 r={**receipt,'variant':out.name,'stage':'prepared','jobs':jobs,'dependency_receipts':dependencies,'inherited_objects':inventory,'original_cdm_library':str(BUILD/'outputs/Release/lib/libbiogears_cdm.so.8.0.0'),'original_cdm_library_sha256':sha(BUILD/'outputs/Release/lib/libbiogears_cdm.so.8.0.0'),'frozen_inputs':frozen,'compiled':{},'builder_sha256':sha(Path(__file__)),'patched_sources':{k:hashlib.sha256(v.encode()).hexdigest() for k,v in patched.items()}}
 bind_inputs(out,r);print(json.dumps({'prepared':str(out),'jobs':jobs}))
def bind_inputs(out,r):
 # Bind donor sources and build recipes even when resuming a prepared build.
 paths=set()
 nervous_dependency=CORE/'CMakeFiles/libbiogears.dir/src/engine/Systems/Nervous.cpp.o.d'
 r['dependency_receipts']['core'].setdefault(str(nervous_dependency),sha(nervous_dependency))
 for job in r['jobs']:
  source=Path(job['source'])
  if not source.is_relative_to(out):
   relative=source.relative_to(SOURCE)
   held=subprocess.check_output(['git','-C',str(SOURCE),'show','3f16a5fa1dade9c511b88d923606fa51cc35e95d:'+str(relative)])
   assert source.read_bytes()==held,'Donor source differs from held Git: '+str(relative)
   paths.add(source)
 for cwd,target in [(CDM,'libbiogears_cdm'),(CORE,'libbiogears')]:
  directory=cwd/f'CMakeFiles/{target}.dir'
  paths.update(directory/name for name in ['flags.make','includes_CXX.rsp','link.txt'])
 # Compiler dependency receipts enumerate transitive headers, including sysroot.
 for family,dependencies in r['dependency_receipts'].items():
  cwd=CDM if family=='cdm' else CORE
  for name,digest in dependencies.items():
   dep=Path(name);assert sha(dep)==digest
   body=dep.read_text().replace('\\\n',' ').split(':',1)[1]
   for token in shlex.split(body):
    path=(cwd/token).resolve()
    if path.is_file() and not path.is_relative_to(out):paths.add(path)
 r['external_build_inputs']={str(p):sha(p) for p in sorted(paths)}
 r['input_binding_stage']=r['stage']
 r['input_binding_note']='Live donor sources checked against held Git; build recipes and transitive donor headers pinned at this stage. Existing compiler commands remain retained; no claim of retroactive input capture.'
 record(out,r)
def guard(out,r):
 for name,digest in r.get('external_build_inputs',{}).items():assert sha(Path(name))==digest,'External build input changed: '+name
 assert sha(PARENT/'manifest.json')==r['parent_manifest_sha256'] and sha(PARENT/'libbiogears.so.8.0.0')==r['parent_library_sha256']
 assert sha(Path(r['original_cdm_library']))==r['original_cdm_library_sha256']
 for family,inventory in r['inherited_objects'].items():
  for path,digest in inventory.items():assert sha(Path(path))==digest,'Inherited object changed'
 for name,digest in r['frozen_inputs'].items():assert sha(out/name)==digest,'Prepared source changed: '+name
 for name,digest in r.get('generated_receipts',{}).items():assert sha(out/name)==digest,'Generated source changed'
 for job in r['jobs']:
  if str(job['index']) in r['compiled']:assert sha(Path(job['object']))==r['compiled'][str(job['index'])]['sha256']
def run(out,r,command,label,cwd,memory=2*1024**3,cpu=60):
 def cap():resource.setrlimit(resource.RLIMIT_AS,(memory,memory));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu))
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(out)+':'+str(RUNTIME/'sysroot/lib')+':'+str(BUILD/'outputs/Release/lib')}
 started=time.monotonic()
 with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
  result=subprocess.Popen(['/usr/bin/time','-v','-o',str(out/(label+'.resources')),'nice','-n','10']+command,cwd=cwd,env=env,stdout=stdout,stderr=stderr,preexec_fn=cap,start_new_session=True)
  try:result.wait(timeout=cpu+10)
  except subprocess.TimeoutExpired:
   os.killpg(result.pid,signal.SIGKILL);result.wait();raise RuntimeError(f'{label} wall limit reached; process group reaped, logs retained in {out}')
 if result.returncode:raise RuntimeError(f'{label} failed ({result.returncode}); logs retained in {out}')
 for name,digest in r.get('external_build_inputs',{}).items():assert sha(Path(name))==digest,'External build input changed during stage: '+name
 r.setdefault('commands',{})[label]={'command':command,'builder_sha256':sha(Path(__file__)),'wall_s':time.monotonic()-started,'resources':(out/(label+'.resources')).read_text()}
def codegen(out,r):
 if r['stage']!='prepared':raise ValueError('Codegen requires prepared stage')
 target=out/'generated/biogears/schema/biogears'
 command=[str(RUNTIME/'sysroot/usr/bin/xsdcxx'),'cxx-tree','--show-sloc','--output-dir',str(target),'--options-file',str(out/'xsd/BioGearsDataModel.cfg'),str(out/'xsd/biogears/BioGearsPhysiology.xsd')]
 run(out,r,command,'codegen',out/'xsd',512*1024**2,30)
 names=['generated/biogears/schema/biogears/BioGearsPhysiology.hxx','generated/biogears/schema/biogears/BioGearsPhysiology.cxx']
 r['generated_receipts']={name:sha(out/name) for name in names}
 for name in names:r['frozen_inputs'].pop(name)
 r['stage']='generated';record(out,r)
def compile_one(out,r,index,replace_unbound=False):
 if r['stage'] not in ('generated','compiling'):raise ValueError('Compile requires generated stage')
 if str(index) in r['compiled']:
  if not replace_unbound or 'input_binding_sha256' in r['compiled'][str(index)]:raise ValueError('Object already compiled with binding, or explicit --replace-unbound required')
  job=r['jobs'][index];archive=out/'superseded-unbound';archive.mkdir(exist_ok=True)
  saved=archive/Path(job['object']).name
  if saved.exists():raise ValueError('Unbound archive already exists')
  shutil.copyfile(job['object'],saved)
  label=f'compile-{index:02d}'
  for suffix in ['stdout','stderr','resources']:shutil.copyfile(out/(label+'.'+suffix),archive/(label+'.'+suffix))
  r.setdefault('superseded_unbound',{})[str(index)]={'object':str(saved),'sha256':sha(saved),'command_receipt':r['commands'][label]}
  del r['compiled'][str(index)];record(out,r)
 elif replace_unbound:raise ValueError('No unbound compiled object to replace')
 job=r['jobs'][index];family=job['family'];cwd=CDM if family=='cdm' else CORE;target='libbiogears_cdm' if family=='cdm' else 'libbiogears'
 flags=(cwd/f'CMakeFiles/{target}.dir/flags.make').read_text();defines=shlex.split(next(line.split('=',1)[1] for line in flags.splitlines() if line.startswith('CXX_DEFINES =')))
 options=shlex.split(next(line.split('=',1)[1] for line in flags.splitlines() if line.startswith('CXX_FLAGS =')))
 includes=shlex.split((cwd/f'CMakeFiles/{target}.dir/includes_CXX.rsp').read_text())
 command=['c++',*defines,*options,'-I'+str(out),'-I'+str(out/'generated'),*includes,'-c',job['source'],'-o',job['object']]
 run(out,r,command,f'compile-{index:02d}',cwd);r['compiled'][str(index)]={'sha256':sha(Path(job['object'])),'input_binding_sha256':hashlib.sha256(json.dumps(r['external_build_inputs'],sort_keys=True).encode()).hexdigest()};r['stage']='compiling';record(out,r)
def link(out,r):
 if len(r['compiled'])!=len(r['jobs']):raise ValueError('Every dependent object must be compiled before link')
 binding=hashlib.sha256(json.dumps(r['external_build_inputs'],sort_keys=True).encode()).hexdigest()
 if any(c.get('input_binding_sha256')!=binding for c in r['compiled'].values()):raise ValueError('All objects must be compiled under the complete current input binding before link')
 for family,cwd,target in [('cdm',CDM,'libbiogears_cdm'),('core',CORE,'libbiogears')]:
  replacements={j['old_object']:j['object'] for j in r['jobs'] if j['family']==family}
  objects=[replacements.get(name,name) for name in r['inherited_objects'][family]];rsp=out/(family+'-objects.rsp');rsp.write_text(' '.join(shlex.quote(x) for x in objects)+'\n')
  command=shlex.split((cwd/f'CMakeFiles/{target}.dir/link.txt').read_text());command[command.index('-o')+1]=str(out/(target+'.so.8.0.0'))
  command=[('@'+str(rsp)) if x.startswith('@') else (str(out/'libbiogears_cdm.so.8.0.0') if x.endswith('/libbiogears_cdm.so.8.0.0') else x) for x in command]
  run(out,r,command,'link-'+family,cwd)
 r['stage']='linked';r['library_sha256']=sha(out/'libbiogears.so.8.0.0');r['cdm_library_sha256']=sha(out/'libbiogears_cdm.so.8.0.0')
 for family,key in [('core','object_sha256'),('cdm','cdm_object_sha256')]:
  replacements={j['old_object']:j['object'] for j in r['jobs'] if j['family']==family}
  r[key]={replacements.get(name,name):sha(Path(replacements.get(name,name))) for name in r['inherited_objects'][family]}
 shutil.copyfile(out/'core-objects.rsp',out/'objects.rsp')
 r['acceptance']='unverified; do not register runtime variant';record(out,r)
 (out/'manifest.json').write_text(json.dumps(r,indent=2)+'\n')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','bind-inputs','codegen','compile','link']);p.add_argument('--directory',type=Path,default=DEFAULT);p.add_argument('--index',type=int);p.add_argument('--replace-unbound',action='store_true');a=p.parse_args();out=a.directory.resolve()
 if a.stage=='prepare':prepare_build(out);return
 r=json.loads((out/'build-state.json').read_text())
 for i,j in enumerate(r['jobs']):j['index']=i
 if a.stage=='bind-inputs':
  guard(out,r);bind_inputs(out,r);print(json.dumps({'external_build_inputs':len(r['external_build_inputs'])}));return
 if 'external_build_inputs' not in r:raise ValueError('Run explicit bind-inputs before resuming this older prepared build')
 guard(out,r)
 if a.stage=='codegen':codegen(out,r)
 elif a.stage=='compile':
  if a.index is None or not 0<=a.index<len(r['jobs']):raise ValueError('Explicit valid compile index required')
  compile_one(out,r,a.index,a.replace_unbound)
 else:link(out,r)
 print(json.dumps({'stage':r['stage'],'compiled_objects':len(r['compiled']),'directory':str(out)}))
if __name__=='__main__':main()
