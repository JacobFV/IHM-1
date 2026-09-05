#!/usr/bin/env python3
"""Correct native renal Tmax mass transfer in a new layer over whole_body_integrity."""
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
    original=source/'projects/biogears/libBiogears/src/engine/Systems/Renal.cpp'
    committed=subprocess.check_output(['git','-C',str(source),'show',f'{REVISION}:{original.relative_to(source)}'])
    if original.read_bytes()!=committed:
        raise RuntimeError('Donor renal source differs from retained revision')
    before=original.read_text()
    start=before.index('void Renal::CalculateReabsorptionTransport(SESubstance& sub)')
    end=before.index('void Renal::CalculateExcretion(SESubstance& sub)',start)
    branch=before[start:end]
    wrong='      reabsorptionRate_mg_Per_s = std::min(reabsorptionRate_mg_Per_s, transportMaximum_mg_Per_s);'
    right='''      // Apply the rate limit to the paired mass transfer before glucose bookkeeping.
      massToMove_mg = std::min(massToMove_mg, transportMaximum_mg_Per_s * m_dt);
      reabsorptionRate_mg_Per_s = massToMove_mg / m_dt;'''
    if branch.count(wrong)!=1:
        raise RuntimeError('Expected exact capped-rate/uncapped-mass defect once')
    after=before[:start]+branch.replace(wrong,right)+before[end:]
    parent=RUNTIME/'variants/whole_body_integrity'
    parent_manifest=json.loads((parent/'manifest.json').read_text())
    for filename,key in [('libbiogears.so.8.0.0','library_sha256'),('Gastrointestinal.cpp','patched_source_sha256')]:
        if sha(parent/filename)!=parent_manifest[key]:
            raise RuntimeError(f'Prior variant integrity failure: {filename}')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(parent_manifest['object_sha256']):
        raise RuntimeError('Parent response list changed')
    for obj in objects:
        if sha(cwd/obj)!=parent_manifest['object_sha256'][obj]:
            raise RuntimeError(f'Parent object changed: {obj}')
    variant=RUNTIME/'variants/whole_body_integrity_renal'
    variant.mkdir(parents=True,exist_ok=True)
    patched=variant/'Renal.cpp';patched.write_text(after)
    patch=variant/'renal_transport_mass.patch'
    patch_text=''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='original/Renal.cpp',tofile='whole_body_integrity_renal/Renal.cpp'))
    patch.write_text(patch_text)
    obj=variant/'Renal.cpp.o'
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    indices=[i for i,o in enumerate(objects) if o.endswith('/src/engine/Systems/Renal.cpp.o')]
    if len(indices)!=1:
        raise RuntimeError('Expected one engine renal object')
    original_obj=cwd/objects[indices[0]]
    objects[indices[0]]=str(obj)
    response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library=variant/'libbiogears.so.8.0.0'
    command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command]
    subprocess.run(command,cwd=cwd,check=True)
    manifest=dict(variant=variant.name,source_revision=revision,parent_variant=parent.name,parent_manifest_sha256=sha(parent/'manifest.json'),parent_library_sha256=parent_manifest['library_sha256'],original_source_sha256=sha(original),patched_source_sha256=sha(patched),patch_sha256=sha(patch),patch_text=patch_text,original_renal_object_sha256=sha(original_obj),library_sha256=sha(library),compile_command=compile_command,link_command=command,build_cwd=str(cwd),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],object_sha256={o:sha(cwd/o) for o in objects},scope='Cap paired tubules-to-peritubular mass by native per-kidney Tmax*dt before glucose bookkeeping; report the actual rate. Preserve source availability, permeability, backflow, units, compartment routing and infinite-Tmax behavior. No physiological parameter calibration.')
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if original.read_bytes()!=committed or sha(parent/'libbiogears.so.8.0.0')!=parent_manifest['library_sha256']:
        raise RuntimeError('Donor or parent changed during build')
    print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2))

if __name__=='__main__':
    build()
