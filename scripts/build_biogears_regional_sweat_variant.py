#!/usr/bin/env python3
"""Fix direct native sweat writes for aggregate regional Skin; immutable child."""
import argparse,difflib,json,os,shlex,subprocess
from pathlib import Path
from verify_native_regional_skin import ROOT,RUNTIME,sha
from build_biogears_regional_skin_variant import limits
PARENT=RUNTIME/'variants/whole_body_integrity_regional_skin_gi_v1'
LIBRARY_PIN='4a836f61e6c123bfb5a5df806c325d4a217f1c707637a215b6f7ff84bc5ef949'
MANIFEST_PIN='ec5250b325ad9950f2b5324524153200c331f399c6ff6e172956f4b0e3b2c31a'
ENERGY_SOURCE_PIN='78bc4dd852cd3e3d0a5ec5f5e3289b08cc10067ab71622cb8a61298735399824'


def corrected_source(before):
    after='#define IHM_REGIONAL_SPECIES_IMPLEMENTATION\n#include "native_regional_species.h"\n'+before
    for index,(species,variable,valence) in enumerate([('Sodium','sodium',1),('Potassium','potassium',1),('Chloride','chloride',-1)]):
        old=f'  m_Skin{species}->GetMass().IncrementValue(-{variable}Lost_mg, MassUnit::mg);\n  Get{species}LostToSweat().IncrementValue({variable}Lost_mg, MassUnit::mg);'
        new=f'  ihm_regional::withdraw_skin_ion(*m_data.GetCompartments().GetLiquidCompartment(BGE::ExtravascularCompartment::SkinExtracellular), *m_Skin{species}, {variable}Lost_mg, Get{species}LostToSweat(), {index}, {valence});'
        if after.count(old)!=1:raise ValueError('Unexpected sweat donor/waste source pair: '+species)
        after=after.replace(old,new)
    return after


def verify_parent():
    if sha(PARENT/'manifest.json')!=MANIFEST_PIN or sha(PARENT/'libbiogears.so.8.0.0')!=LIBRARY_PIN:raise ValueError('Regional parent changed')
    manifest=json.loads((PARENT/'manifest.json').read_text());objects=shlex.split((PARENT/'objects.rsp').read_text());cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    if len(objects)!=len(set(objects)) or set(objects)!=set(manifest['object_sha256']):raise ValueError('Parent inventory mismatch')
    for obj in objects:
        if sha(cwd/obj)!=manifest['object_sha256'][obj]:raise ValueError('Parent object changed: '+obj)
    return manifest,objects,cwd


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--build-library',action='store_true');parser.add_argument('--output',type=Path,default=RUNTIME/'variants/whole_body_integrity_regional_skin_sweat_v1');args=parser.parse_args()
    parent,objects,cwd=verify_parent();matches=[o for o in objects if o.endswith('/Energy.cpp.o')]
    if len(matches)!=1:raise ValueError('Expected one Energy object')
    source=(cwd/matches[0]).with_suffix('')
    if sha(source)!=ENERGY_SOURCE_PIN:raise ValueError('Inherited signed Energy source changed')
    before=source.read_text();after=corrected_source(before);out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    patched=out/'Energy.cpp';patched.write_text(after)
    inputs={}
    for name in ['native_signed_muscle_port.h','native_regional_skin.h']:
        p=PARENT/name;(out/name).write_bytes(p.read_bytes());inputs[name]=sha(p)
    header=ROOT/'scripts/native_regional_species.h';(out/header.name).write_bytes(header.read_bytes());inputs[header.name]=sha(header)
    patch=out/'regional_sweat_ownership.patch';patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='parent/Energy.cpp',tofile='regional_sweat/Energy.cpp')))
    obj=out/'Energy.cpp.o';objects[objects.index(matches[0])]=str(obj);response=out/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O2','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    link=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());library=out/'libbiogears.so.8.0.0';link[link.index('-o')+1]=str(library);token='@CMakeFiles/libbiogears.dir/objects1.rsp'
    if link.count(token)!=1:raise ValueError('Unexpected link response')
    link=['@'+str(response) if c==token else c for c in link]
    preparation={'variant':out.name,'parent_variant':PARENT.name,'parent_manifest_sha256':MANIFEST_PIN,'parent_library_sha256':LIBRARY_PIN,
        'header_sha256':parent['header_sha256'],'regional_skin_header_sha256':parent['regional_skin_header_sha256'],'regional_species_header_sha256':sha(header),
        'source_sha256':ENERGY_SOURCE_PIN,'patched_source_sha256':sha(patched),'header_receipts':inputs,'builder_sha256':sha(Path(__file__)),
        'replaced_object':matches[0],'compile_command':command,'link_command':link,'build_cwd':str(cwd),'production_ready':False,
        'scope':'Native sweat Na/K/Cl donor updates over owning leaves paired with unchanged native cumulative sweat counters; all other native laws inherited'}
    (out/'preparation.json').write_text(json.dumps(preparation,indent=2)+'\n');print(out,flush=True)
    if not args.build_library:return
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    for label,cmd in [('compile',command),('link',link)]:
        with (out/(label+'.log')).open('w') as stream:subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources'))]+cmd,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,preexec_fn=limits,timeout=210,check=True)
    verify_parent()
    if sha(source)!=ENERGY_SOURCE_PIN:raise ValueError('Inherited source changed during build')
    manifest={**preparation,'library_sha256':sha(library),'object_sha256':{o:sha(cwd/o) for o in objects},'inherited_object_count':len(objects)-1,'resources':{name:(out/(name+'.resources')).read_text() for name in ['compile','link']}}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'built':True,'library_sha256':sha(library),'manifest_sha256':sha(out/'manifest.json')}))


if __name__=='__main__':main()
