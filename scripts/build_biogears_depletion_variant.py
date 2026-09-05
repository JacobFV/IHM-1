#!/usr/bin/env python3
"""Empty an explicitly exhausted mixed-unit stomach-water pool exactly."""
from pathlib import Path
import difflib
import hashlib
import json
import shlex
import subprocess

BASE=Path(__file__).resolve().parents[1]
RUNTIME=BASE/'data/runtime/physiology'
REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build():
    source=BASE/'data/raw/physiology/biogears'
    revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    if revision!=REVISION:
        raise RuntimeError('Unexpected source revision')
    original=source/'projects/biogears/libBiogears/src/engine/Systems/Gastrointestinal.cpp'
    committed=subprocess.check_output(['git','-C',str(source),'show',f'{REVISION}:{original.relative_to(source)}'])
    if original.read_bytes()!=committed:raise RuntimeError('Donor GI source changed')
    ancestor=RUNTIME/'variants/whole_body_integrity_gi_water'
    ancestor_manifest=json.loads((ancestor/'manifest.json').read_text())
    before=(ancestor/'Gastrointestinal.cpp').read_text()
    if sha(ancestor/'Gastrointestinal.cpp')!=ancestor_manifest['patched_source_sha256']:raise RuntimeError('GI-water source changed')
    wrong='    m_StomachContents->GetWater().IncrementValue(-digestedAmount, VolumeUnit::mL);'
    right="""    if (digestedAmount > 0.0 && digestedAmount == waterContent_mL) {
      // The paired transfer exhausts this known pool. Cross-unit subtraction
      // can round below zero; retain the established exact empty state.
      m_StomachContents->GetWater().SetValue(0.0, VolumeUnit::mL);
    } else {
      m_StomachContents->GetWater().IncrementValue(-digestedAmount, VolumeUnit::mL);
    }"""
    if before.count(wrong)!=1:raise RuntimeError('Expected exact water debit once')
    after=before.replace(wrong,right,1)
    parent=RUNTIME/'variants/whole_body_integrity_energy'
    parent_manifest=json.loads((parent/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=parent_manifest['library_sha256']:raise RuntimeError('Energy parent changed')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(parent_manifest['object_sha256']):
        raise RuntimeError('Parent response list changed')
    for obj in objects:
        if sha(cwd/obj)!=parent_manifest['object_sha256'][obj]:
            raise RuntimeError(f'Parent object changed: {obj}')
    variant=RUNTIME/'variants/whole_body_integrity_depletion'
    variant.mkdir(parents=True,exist_ok=True)
    patched=variant/'Gastrointestinal.cpp';patched.write_text(after)
    patch=variant/'known_water_depletion.patch'
    patch_text=''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='whole_body_integrity_gi_water/Gastrointestinal.cpp',tofile='whole_body_integrity_depletion/Gastrointestinal.cpp'))
    patch.write_text(patch_text)
    obj=variant/'Gastrointestinal.cpp.o'
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    indices=[i for i,o in enumerate(objects) if o==str(ancestor/'Gastrointestinal.cpp.o')]
    if len(indices)!=1:
        raise RuntimeError('Expected one inherited GI-water object')
    original_obj=cwd/objects[indices[0]]
    objects[indices[0]]=str(obj)
    response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library=variant/'libbiogears.so.8.0.0'
    command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command]
    subprocess.run(command,cwd=cwd,check=True)
    manifest=dict(variant=variant.name,source_revision=revision,parent_variant=parent.name,parent_manifest_sha256=sha(parent/'manifest.json'),parent_library_sha256=parent_manifest['library_sha256'],original_source_sha256=sha(original),patched_source_sha256=sha(patched),patch_sha256=sha(patch),patch_text=patch_text,original_gi_object_sha256=sha(original_obj),library_sha256=sha(library),compile_command=compile_command,link_command=command,build_cwd=str(cwd),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],object_sha256={o:sha(cwd/o) for o in objects},scope='Known positive full-pool water transfer sets the established remainder to exact zero in mL. Partial transfer and arbitrary invalid input rejection unchanged. Inherits all energy/GI/renal/thermal objects.')
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if original.read_bytes()!=committed or sha(parent/'libbiogears.so.8.0.0')!=parent_manifest['library_sha256']:
        raise RuntimeError('Donor or parent changed during build')
    print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2))

if __name__=='__main__':
    build()
