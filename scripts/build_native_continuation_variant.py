#!/usr/bin/env python3
"""Isolated source-composed sleep/GI/Circuit/Tissue continuation variant.

Explicit prepare/codegen/compile/link stages; each native stage needs root slot.
"""
from pathlib import Path
import argparse,hashlib,json,shlex,shutil,subprocess
import build_native_nervous_sleep_variant as sleep_build
import prepare_gi_serialization_repair as gi
import prepare_tissue_burn_state_patch as burn
import patch_cardiovascular_region_io as circuit
import prepare_gi_inventory_presence as gi_presence
ROOT=sleep_build.ROOT;RUNTIME=sleep_build.RUNTIME;SOURCE=sleep_build.SOURCE;CORE=sleep_build.CORE;CDM=sleep_build.CDM;BUILD=sleep_build.BUILD
PARENT=sleep_build.DEFAULT;DEFAULT=RUNTIME/'variants/whole_body_integrity_continuation_v3'
PARENT_MANIFEST='245de131d8ff022c33b46c12cf62ff510fd1fbc7b8aee36a6b225a604f2ab94a'
GI_MANIFEST=ROOT/'data/derived/audits/gi-codec-full-tu-3juo1cj5/objects_manifest.json'
sha=sleep_build.sha;record=sleep_build.record

def parent():
 assert sha(PARENT/'manifest.json')==PARENT_MANIFEST
 p=json.loads((PARENT/'manifest.json').read_text());sleep_build.guard(PARENT,p)
 assert sha(PARENT/'libbiogears.so.8.0.0')==p['library_sha256'] and sha(PARENT/'libbiogears_cdm.so.8.0.0')==p['cdm_library_sha256']
 return p

def source_composition(out,p):
 gi_out=out/'gi-provenance';gi.prepare(gi_out)
 accepted=json.loads(GI_MANIFEST.read_text())['objects'];gi_sources={}
 for rec in accepted.values():
  relative=Path(rec['patched_source']).relative_to(GI_MANIFEST.parent)
  generated=gi_out/relative;assert sha(generated)==rec['patched_source_sha256'];assert sha(Path(rec['object']))==rec['sha256']
  gi_sources[str(relative)]=generated.read_text()
 raw_io=gi.BASE/gi.FILES[0];signature='void BiogearsPhysiology::Marshall(const Gastrointestinal& in, CDM::BioGearsGastrointestinalSystemData& out)'
 native_io=(PARENT/'BioGearsPhysiology.cpp').read_text();old=gi.method(raw_io.read_text(),signature);assert gi.method(native_io,signature)==old
 native_io=native_io.replace(old,gi.method(gi_sources[gi.FILES[0]],signature),1)
 composed_io=gi_presence.patch_io(burn.patch_io(native_io))
 # Both sleep and burn source transformations must survive the shared TU.
 assert 'IHMSleepState' in composed_io and 'IHMBurnHistory' in composed_io and 'out.DrugTransitStates().clear()' in composed_io
 tissue_object=next(name for name in p['object_sha256'] if Path(name).name=='Tissue.cpp.o')
 tissue_source=Path(str(tissue_object).removesuffix('.o'));tissue_manifest=json.loads((tissue_source.parent/'manifest.json').read_text())
 assert sha(tissue_source)==tissue_manifest['sources']['Tissue']['patched_source_sha256']
 gi_parent=RUNTIME/'variants/whole_body_integrity_gi_absorption/Gastrointestinal.cpp'
 gi_receipt=json.loads((gi_parent.parent/'manifest.json').read_text());assert sha(gi_parent)==gi_receipt['patched_source_sha256']
 raw={'src/engine/Systems/Gastrointestinal.cpp':gi_presence.patch_gi(gi_parent.read_text()),gi.FILES[0]:composed_io,gi.FILES[1]:gi_sources[gi.FILES[1]],gi.FILES[2]:gi_sources[gi.FILES[2]],'src/io/cdm/Circuit.cpp':circuit.patch_io(circuit.IO.read_text()),'src/engine/Systems/Tissue.cpp':burn.patch_tissue(tissue_source.read_text())}
 for name,text in raw.items():
  path=out/'sources'/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
 schema=burn.patch_schema((PARENT/'xsd/biogears/BioGearsPhysiology.xsd').read_text())
 start=schema.index('  <xs:complexType name="BioGearsGastrointestinalSystemData">');end=schema.index('  </xs:complexType>',start)
 block=schema[start:end];assert block.count('        </xs:sequence>')==1
 block=block.replace('        </xs:sequence>','          <xs:element name="GITransitInventoryVersion" type="xs:unsignedInt" minOccurs="0" maxOccurs="1"/>\n        </xs:sequence>')
 (out/'xsd/biogears/BioGearsPhysiology.xsd').write_text(schema[:start]+block+schema[end:])
 header=SOURCE/'projects/biogears/libBiogears/include/biogears/engine/Systems/Gastrointestinal.h'
 patched_header=out/'include/biogears/engine/Systems/Gastrointestinal.h';patched_header.parent.mkdir(parents=True);patched_header.write_text(gi_presence.patch_header(header.read_text()))
 return {'gi_accepted_manifest_sha256':sha(GI_MANIFEST),'tissue_parent_source':str(tissue_source),'tissue_parent_source_sha256':sha(tissue_source),'composed_source_sha256':{name:hashlib.sha256(text.encode()).hexdigest() for name,text in raw.items()}}

def prepare(out):
 if out.exists():raise ValueError('New composition output directory required')
 p=parent();out.mkdir(parents=True)
 shutil.copytree(PARENT/'generated',out/'generated');shutil.copytree(PARENT/'xsd',out/'xsd')
 receipt=source_composition(out,p)
 for name in ['native_nervous_sleep_state.h','native_signed_muscle_port.h']:shutil.copyfile(PARENT/name,out/name)
 shutil.copyfile(ROOT/'scripts/native_tissue_burn_state.h',out/'native_tissue_burn_state.h')
 inventories={'core':dict(p['object_sha256']),'cdm':dict(p['cdm_object_sha256'])}
 for inventory in inventories.values():
  for name,digest in inventory.items():assert sha(Path(name))==digest
 old_to_parent={j['old_object']:j['object'] for j in p['jobs']}
 def inherited(family,obj):
  name=str(obj.resolve());candidate=old_to_parent.get(name,name)
  if candidate in inventories[family]:return candidate
  candidates={name for name in inventories[family] if Path(name).name==obj.name and not Path(name).is_relative_to(CORE)}
  candidates.update(mapped for old,mapped in old_to_parent.items() if Path(old).name==obj.name and mapped in inventories[family])
  candidates=list(candidates)
  if len(candidates)!=1:raise ValueError('Ambiguous inherited object for '+str(obj))
  return candidates[0]
 jobs=[];dependencies={};donors={}
 by_object={j['object']:j for j in p['jobs']}
 def donor_source(obj,original):
  if obj in by_object:
   path=Path(by_object[obj]['source']);manifest=PARENT/'manifest.json'
  elif Path(obj).is_relative_to(CORE):
   path=original;manifest=None
   held=subprocess.check_output(['git','-C',str(SOURCE),'show','3f16a5fa1dade9c511b88d923606fa51cc35e95d:'+str(path.relative_to(SOURCE))]);assert path.read_bytes()==held
  else:path=Path(obj.removesuffix('.o'));manifest=path.parent/'manifest.json'
  digest=sha(path)
  if manifest:assert digest in manifest.read_text(),'Inherited CPP lacks accepted source receipt: '+str(path)
  donors[str(original)]={'object':obj,'object_sha256':inventories['core'][obj],'source':str(path),'source_sha256':digest,'manifest':str(manifest) if manifest else None,'manifest_sha256':sha(manifest) if manifest else None}
  return path
 def add(family,dep):
  cwd=CDM if family=='cdm' else CORE;target='libbiogears_cdm' if family=='cdm' else 'libbiogears';obj=Path(str(dep)[:-2]);relative=obj.relative_to(cwd/f'CMakeFiles/{target}.dir')
  if family=='cdm':source=out/'generated'/str(relative).removeprefix('__/generated/Release/').removesuffix('.o');original=source
  else:
   original=SOURCE/'projects/biogears/libBiogears'/str(relative).removesuffix('.o');source=out/'sources'/str(relative).removesuffix('.o')
   if not source.exists():
    donor=donor_source(inherited(family,obj),original);source.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(donor,source)
  assert source.is_file();old=inherited(family,obj)
  if any(j['old_object']==old for j in jobs):return
  index=len(jobs);jobs.append({'index':index,'family':family,'old_object':old,'source':str(source),'original_source':str(original),'object':str(out/'objects'/f'{index:02}.o')})
  dependencies[str(dep)]=sha(dep)
 for family,cwd in [('cdm',CDM),('core',CORE)]:
  for dep in sorted(p for p in cwd.rglob('*.o.d') if 'BioGearsPhysiology.hxx' in p.read_text() or (family=='core' and 'engine/Systems/Gastrointestinal.h' in p.read_text())):add(family,dep)
 for name in sorted(receipt['composed_source_sha256']):add('core',CORE/f'CMakeFiles/libbiogears.dir/{name}.o.d')
 assert sum(j['family']=='cdm' for j in jobs)==3
 assert len([d for d in dependencies if 'engine/Systems/Gastrointestinal.h' in Path(d).read_text()])==30
 (out/'objects').mkdir()
 external=set()
 for name,digest in dependencies.items():
  dep=Path(name);assert sha(dep)==digest;external.add(dep)
  body=dep.read_text().replace('\\\n',' ').split(':',1)[1];cwd=CDM if dep.is_relative_to(CDM) else CORE
  for token in shlex.split(body):
   path=(cwd/token).resolve()
   if path.is_file() and not path.is_relative_to(out):external.add(path)
 for cwd,target in [(CORE,'libbiogears'),(CDM,'libbiogears_cdm')]:
  folder=cwd/f'CMakeFiles/{target}.dir'
  for name in ['flags.make','includes_CXX.rsp','link.txt']:external.add(folder/name)
  # Guard the original archives/shared objects named by the held link recipe.
  for token in shlex.split((folder/'link.txt').read_text()):
   path=(cwd/token).resolve()
   if path.is_file() and ('.so' in path.name or path.suffix=='.a'):external.add(path)
 for name in ['prepare_gi_serialization_repair.py','prepare_tissue_burn_state_patch.py','patch_cardiovascular_region_io.py','prepare_gi_inventory_presence.py','native_tissue_burn_state.h']:external.add(ROOT/'scripts'/name)
 external.add(GI_MANIFEST)
 for donor in donors.values():
  external.add(Path(donor['source']))
  if donor['manifest']:external.add(Path(donor['manifest']))
 frozen={str(path.relative_to(out)):sha(path) for path in out.rglob('*') if path.is_file()}
 r={**receipt,'schema':'ihm.native-continuation-composition.v1','variant':out.name,'stage':'prepared','parent_variant':PARENT.name,'parent_manifest_sha256':PARENT_MANIFEST,'parent_library_sha256':p['library_sha256'],'parent_cdm_library_sha256':p['cdm_library_sha256'],'jobs':jobs,'dependency_receipts':dependencies,'source_donor_receipts':donors,'inherited_objects':inventories,'frozen_inputs':frozen,'external_build_inputs':{str(path):sha(path) for path in sorted(external)},'compiled':{},'builder_sha256':sha(Path(__file__)),'acceptance':'Unverified composed ABI; no default runtime promotion'}
 record(out,r);guard(out,r);print(json.dumps({'directory':str(out),'jobs':jobs,'external_inputs':len(external)}))

def guard(out,r):
 p=parent();assert p['library_sha256']==r['parent_library_sha256'] and p['cdm_library_sha256']==r['parent_cdm_library_sha256']
 for inventory in r['inherited_objects'].values():
  for name,digest in inventory.items():assert sha(Path(name))==digest,'Inherited object changed'
 for name,digest in r['external_build_inputs'].items():assert sha(Path(name))==digest,'External composition input changed: '+name
 for name,digest in r['frozen_inputs'].items():assert sha(out/name)==digest,'Frozen composition input changed: '+name
 for name,digest in r.get('generated_receipts',{}).items():assert sha(out/name)==digest
 for key,c in r['compiled'].items():assert sha(Path(r['jobs'][int(key)]['object']))==c['sha256']

def run(out,r,command,label,cwd,**caps):
 sleep_build.run(out,r,command,label,cwd,**caps)
 r['commands'][label]['composition_builder_sha256']=sha(Path(__file__))

def codegen(out,r):
 if r['stage']!='prepared':raise ValueError('Expected prepared stage')
 command=[str(RUNTIME/'sysroot/usr/bin/xsdcxx'),'cxx-tree','--show-sloc','--output-dir',str(out/'generated/biogears/schema/biogears'),'--options-file',str(out/'xsd/BioGearsDataModel.cfg'),str(out/'xsd/biogears/BioGearsPhysiology.xsd')]
 run(out,r,command,'codegen',out/'xsd',memory=512*1024**2,cpu=30)
 names=['generated/biogears/schema/biogears/BioGearsPhysiology.hxx','generated/biogears/schema/biogears/BioGearsPhysiology.cxx']
 for name in names:r['frozen_inputs'].pop(name)
 r['generated_receipts']={name:sha(out/name) for name in names};r['stage']='generated';record(out,r)

def compile_one(out,r,index):
 if r['stage'] not in ('generated','compiling') or str(index) in r['compiled']:raise ValueError('Uncompiled job after codegen required')
 job=r['jobs'][index];family=job['family'];cwd=CDM if family=='cdm' else CORE;target='libbiogears_cdm' if family=='cdm' else 'libbiogears';folder=cwd/f'CMakeFiles/{target}.dir'
 flags=(folder/'flags.make').read_text();values=lambda key:shlex.split(next(line.split('=',1)[1] for line in flags.splitlines() if line.startswith(key+' =')))
 command=['c++',*values('CXX_DEFINES'),*values('CXX_FLAGS'),'-I',str(out/'include'),'-I',str(out),'-I',str(out/'generated'),'@'+str(folder/'includes_CXX.rsp'),'-iquote',str(Path(job['original_source']).parent),'-c',job['source'],'-o',job['object']]
 run(out,r,command,f'compile-{index:02}',cwd);guard(out,r)
 r['compiled'][str(index)]={'sha256':sha(Path(job['object'])),'source_sha256':sha(Path(job['source'])),'external_inputs_sha256':hashlib.sha256(json.dumps(r['external_build_inputs'],sort_keys=True).encode()).hexdigest(),'generated_receipts':dict(r['generated_receipts'])};r['stage']='compiling';record(out,r)

def link(out,r):
 if len(r['compiled'])!=len(r['jobs']):raise ValueError('Complete matched compile set required')
 for family,cwd,target in [('cdm',CDM,'libbiogears_cdm'),('core',CORE,'libbiogears')]:
  replacements={j['old_object']:j['object'] for j in r['jobs'] if j['family']==family};objects=[replacements.get(name,name) for name in r['inherited_objects'][family]]
  rsp=out/(family+'-objects.rsp');rsp.write_text(' '.join(shlex.quote(name) for name in objects)+'\n')
  command=shlex.split((cwd/f'CMakeFiles/{target}.dir/link.txt').read_text());command[command.index('-o')+1]=str(out/(target+'.so.8.0.0'))
  command=[('@'+str(rsp)) if token.startswith('@') else (str(out/'libbiogears_cdm.so.8.0.0') if token.endswith('/libbiogears_cdm.so.8.0.0') else token) for token in command]
  run(out,r,command,'link-'+family,cwd,memory=4*1024**3,cpu=50);guard(out,r)
  r['cdm_object_sha256' if family=='cdm' else 'object_sha256']={name:sha(Path(name)) for name in objects}
 shutil.copyfile(out/'core-objects.rsp',out/'objects.rsp');r['library_sha256']=sha(out/'libbiogears.so.8.0.0');r['cdm_library_sha256']=sha(out/'libbiogears_cdm.so.8.0.0');r['stage']='linked';record(out,r)
 (out/'manifest.json').write_text(json.dumps(r,indent=2)+'\n')
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['prepare','codegen','compile','link']);parser.add_argument('--directory',type=Path,default=DEFAULT);parser.add_argument('--index',type=int);args=parser.parse_args();out=args.directory.resolve()
 if args.stage=='prepare':prepare(out);return
 r=json.loads((out/'build-state.json').read_text());guard(out,r)
 if args.stage=='codegen':codegen(out,r)
 elif args.stage=='compile':
  if args.index is None or not 0<=args.index<len(r['jobs']):raise ValueError('Explicit valid job index required')
  compile_one(out,r,args.index)
 else:link(out,r)
 print(json.dumps({'directory':str(out),'stage':r['stage'],'compiled':len(r['compiled'])}))
if __name__=='__main__':main()
