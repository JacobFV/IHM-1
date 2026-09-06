#!/usr/bin/env python3
"""Prepare an immutable regional-Skin variant over the corrected GI lineage.

Default prepares source/commands only. --build-library requires a root compile
slot. The header is opt-in and not installed into existing native sessions.
"""
import argparse
import difflib
import json
import os
from pathlib import Path
import resource
import shlex
import subprocess

from verify_native_regional_skin import ROOT,RUNTIME,SOURCE,sha
from build_biogears_shared_donor_variant import corrected_source as shared_donor_source

PARENT=RUNTIME/'variants/whole_body_integrity_gi_absorption'
LIBRARY_PIN='9792d857c47a5907f571a03495fe9f4f1144f114afd72c0451869e1a7049588b'
MANIFEST_PIN='de1aa254b868b4b51e3fbd4f90323e09370a409957d83ac7552a09bb0f3e7125'

LOOKUPS = ('fluidFluxPathName', 'fluidPathName')


def regional_source(before):
    after=before
    for name in LOOKUPS:
        old='m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath('+name+')'
        new='m_data.GetCircuits().GetFluidPath('+name+')'
        if after.count(old)!=1:
            raise ValueError('Unexpected source-law lookup count: '+name)
        after=after.replace(old,new)
    return after


def verify_parent():
    if sha(PARENT/'manifest.json')!=MANIFEST_PIN or sha(PARENT/'libbiogears.so.8.0.0')!=LIBRARY_PIN:
        raise ValueError('Immutable signed parent changed')
    manifest=json.loads((PARENT/'manifest.json').read_text())
    objects=shlex.split((PARENT/'objects.rsp').read_text())
    if len(objects)!=len(set(objects)) or set(objects)!=set(manifest['object_sha256']):
        raise ValueError('Parent object inventory mismatch')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    for obj in objects:
        if sha(cwd/obj)!=manifest['object_sha256'][obj]:raise ValueError('Changed parent object: '+obj)
    return manifest,objects,cwd


def limits():
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(180,180))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=RUNTIME/'variants/whole_body_integrity_regional_skin_gi_v1')
    parser.add_argument('--build-library',action='store_true')
    args=parser.parse_args()
    parent,objects,cwd=verify_parent()
    donor=SOURCE/'projects/biogears/libBiogears/src/engine/Systems/Diffusion.cpp'
    before=shared_donor_source(donor.read_text())
    matches=[o for o in objects if o.endswith('/Diffusion.cpp.o')]
    if len(matches)!=1:raise ValueError('Expected one inherited Diffusion object')
    inherited_source=(cwd/matches[0]).with_suffix('')
    if inherited_source.read_text()!=before:raise ValueError('Inherited shared-donor source mismatch')
    signed_header=PARENT/'native_signed_muscle_port.h'
    if not signed_header.exists():
        signed_header=RUNTIME/'variants'/parent['signed_port_abi_parent']/'native_signed_muscle_port.h'
    if sha(signed_header)!=parent['header_sha256']:raise ValueError('Inherited signed header changed')
    after=regional_source(before)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    patched=out/'Diffusion.cpp';patched.write_text(after)
    header=out/'native_regional_skin.h';header.write_bytes((ROOT/'scripts/native_regional_skin.h').read_bytes())
    (out/signed_header.name).write_bytes(signed_header.read_bytes())
    patch=out/'regional_skin_source_lookup.patch'
    patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='parent/Diffusion.cpp',tofile='regional_skin/Diffusion.cpp')))
    obj=out/'Diffusion.cpp.o';replacement=list(objects);replacement[replacement.index(matches[0])]=str(obj)
    response=out/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in replacement))
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O2','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    link_command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library=out/'libbiogears.so.8.0.0';link_command[link_command.index('-o')+1]=str(library)
    token='@CMakeFiles/libbiogears.dir/objects1.rsp'
    if link_command.count(token)!=1:raise ValueError('Unexpected native link response')
    link_command=['@'+str(response) if c==token else c for c in link_command]
    preparation={'variant':out.name,'parent_variant':PARENT.name,'parent_manifest_sha256':MANIFEST_PIN,
        'parent_library_sha256':LIBRARY_PIN,'original_source_sha256':sha(donor),'inherited_source_sha256':sha(inherited_source),
        'patched_source_sha256':sha(patched),'regional_skin_header_sha256':sha(header),
        'header_sha256':parent['header_sha256'],'patch_sha256':sha(patch),
        'builder_sha256':sha(Path(__file__)),'replaced_object':matches[0],
        'engineering_regions':{'region_a':.2,'region_b':.3,'residual':.5},
        'compile_command':compile_command,'link_command':link_command,'build_cwd':str(cwd),
        'scope':'Opt-in native parallel Skin circuit and aggregate-law lookup compatibility; no production adapter activation, local chemical-transport validation or measured fractions',
        'source_laws':'Unchanged native algebra; detached aggregate flow caches remain visible through circuit manager',
        'production_ready':False}
    (out/'preparation.json').write_text(json.dumps(preparation,indent=2)+'\n')
    print(out,flush=True)
    if not args.build_library:return
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    def run(command,label):
        with (out/(label+'.log')).open('w') as stream:
            subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources'))]+command,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,preexec_fn=limits,timeout=210,check=True)
    run(compile_command,'compile');run(link_command,'link')
    verify_parent()
    if sha(donor)!=preparation['original_source_sha256']:raise ValueError('Held donor changed')
    manifest={**preparation,'library_sha256':sha(library),'object_sha256':{o:sha(cwd/o) for o in replacement},
        'inherited_object_count':len(replacement)-1,'replacement_object_sha256':sha(obj),
        'resources':{name:(out/(name+'.resources')).read_text() for name in ['compile','link']}}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'built':True,'library_sha256':manifest['library_sha256']}))


if __name__=='__main__':main()
