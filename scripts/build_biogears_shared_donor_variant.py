#!/usr/bin/env python3
"""Replace only native Diffusion.cpp.o over the retained substrate/thermal lineage."""
import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import resource
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT/'data/runtime/physiology'
SOURCE = ROOT/'data/raw/physiology/biogears'
DONOR = SOURCE/'projects/biogears/libBiogears/src/engine/Systems/Diffusion.cpp'
SOURCE_PIN = '9f90e429a87ae1ba3b7071fc14525905e768d23a27a0f21d15a2581a49b7d48b'
REVISION = '3f16a5fa1dade9c511b88d923606fa51cc35e95d'
PARENT = RUNTIME/'variants/whole_body_integrity_substrate_availability'
PARENT_LIBRARY_PIN = '429c63730b3be08b8cf158b16c34817185424cadf7a240abdcd49b4b98d687f7'
VARIANT = RUNTIME/'variants/whole_body_integrity_shared_donor'
BUDGET_PATCH = '''  // Both reverse VE and forward EI withdraw from the same starting donor.
  // Allocate them together before any recipient credit or donor debit.
  if (massToMoveVE_ug < 0.0 && massToMoveEI_ug > 0.0) {
    const double availableExtracellular_ug = extracellularSubQ->GetMass(MassUnit::ug);
    const double combinedWithdrawal_ug = massToMoveEI_ug - massToMoveVE_ug;
    if (combinedWithdrawal_ug > availableExtracellular_ug) {
      const double scale = availableExtracellular_ug / combinedWithdrawal_ug;
      massToMoveVE_ug *= scale;
      massToMoveEI_ug *= scale;
      // Rounding must not exceed the shared budget. Adjust the transfer pair,
      // never a donor-only clamp after recipients have already been credited.
      while (massToMoveEI_ug - massToMoveVE_ug > availableExtracellular_ug) {
        if (massToMoveEI_ug >= -massToMoveVE_ug)
          massToMoveEI_ug = std::nextafter(massToMoveEI_ug, 0.0);
        else
          massToMoveVE_ug = std::nextafter(massToMoveVE_ug, 0.0);
      }
    }
  }

'''


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def corrected_source(before):
    if hashlib.sha256(before.encode()).hexdigest() != SOURCE_PIN:
        raise RuntimeError('Unexpected held diffusion source')
    method = before.index('void DiffusionCalculator::CalculateFacilitatedDiffusion(')
    anchor = before.index('  //For all distribute methods, we use VolumeWeighted', method)
    return before[:anchor]+BUDGET_PATCH+before[anchor:]


def verify_parent(cwd):
    manifest_path = PARENT/'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if sha(PARENT/'libbiogears.so.8.0.0') != PARENT_LIBRARY_PIN or manifest['library_sha256'] != PARENT_LIBRARY_PIN:
        raise RuntimeError('Retained substrate parent library changed')
    objects = shlex.split((PARENT/'objects.rsp').read_text())
    if len(objects) != len(set(objects)) or set(objects) != set(manifest['object_sha256']):
        raise RuntimeError('Parent object inventory differs from manifest')
    for obj in objects:
        if sha(cwd/obj) != manifest['object_sha256'][obj]:
            raise RuntimeError('Parent object changed: '+obj)
    return manifest, objects


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (180,180))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args()
    revision=subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip()
    if revision != REVISION or sha(DONOR) != SOURCE_PIN: raise RuntimeError('Donor identity changed')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    parent,objects=verify_parent(cwd)
    before=DONOR.read_text();after=corrected_source(before)
    # Never overwrite an existing completed or partial isolated build.
    VARIANT.mkdir(parents=True,exist_ok=False)
    patched=VARIANT/'Diffusion.cpp';patched.write_text(after)
    patch=VARIANT/'shared_extracellular_budget.patch'
    patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='donor/Diffusion.cpp',tofile=VARIANT.name+'/Diffusion.cpp')))
    preparation={'source_revision':revision,'original_source_sha256':sha(DONOR),'patched_source_sha256':sha(patched),
                 'parent_manifest_sha256':sha(PARENT/'manifest.json'),'parent_library_sha256':PARENT_LIBRARY_PIN,
                 'scope':'Shared extracellular donor allocation only; all other parent objects preserved.'}
    (VARIANT/'preparation.json').write_text(json.dumps(preparation,indent=2)+'\n')
    print(VARIANT,flush=True)
    if args.prepare_only:return
    obj=VARIANT/'Diffusion.cpp.o'
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    def run(command,name):
        with (VARIANT/(name+'.log')).open('w') as log:
            subprocess.run(['/usr/bin/time','-v','-o',str(VARIANT/(name+'.resources'))]+command,cwd=cwd,env=env,
                           stdout=log,stderr=subprocess.STDOUT,check=True,preexec_fn=limits,timeout=210)
    run(compile_command,'compile')
    original_obj='CMakeFiles/libbiogears.dir/src/engine/Systems/Diffusion.cpp.o'
    if objects.count(original_obj)!=1:raise RuntimeError('Expected one original Diffusion object')
    inherited=[o for o in objects if o!=original_obj]
    objects[objects.index(original_obj)]=str(obj)
    response=VARIANT/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library=VARIANT/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(library)
    response_token='@CMakeFiles/libbiogears.dir/objects1.rsp'
    if command.count(response_token)!=1:raise RuntimeError('Expected one inherited link response')
    command=['@'+str(response) if c==response_token else c for c in command]
    run(command,'link')
    # Validate parent again after compilation/link to catch concurrent mutation.
    verify_parent(cwd)
    if sha(DONOR)!=SOURCE_PIN:raise RuntimeError('Donor changed during build')
    manifest={**preparation,'variant':VARIANT.name,'parent_variant':PARENT.name,'patch_sha256':sha(patch),
              'library_sha256':sha(library),'compile_command':compile_command,'link_command':command,'build_cwd':str(cwd),
              'object_sha256':{o:sha(cwd/o) for o in objects},'inherited_object_count':len(inherited),
              'replaced_object':original_obj,'original_object_sha256':parent['object_sha256'][original_obj],
              'replacement_object_sha256':sha(obj),'builder_sha256':sha(Path(__file__)),
              'limits':{'address_space_bytes':4*1024**3,'cpu_seconds_per_child':180,'wall_seconds_per_child':210,'threads':1},
              'resources':{name:(VARIANT/(name+'.resources')).read_text() for name in ['compile','link']}}
    (VARIANT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ['variant','library_sha256','patched_source_sha256','inherited_object_count']},indent=2))


if __name__=='__main__':main()
