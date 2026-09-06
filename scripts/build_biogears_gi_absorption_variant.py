#!/usr/bin/env python3
"""Conservative GI absorption remainders over immutable signed-muscle-v2 lineage."""
import difflib
import json
import os
from pathlib import Path
import shlex
import subprocess
from build_biogears_shared_donor_variant import ROOT,RUNTIME,sha,limits

PARENT=RUNTIME/'variants/whole_body_integrity_signed_muscle_v2'
PARENT_PIN='5a02e2942e6ac5ff3cdce352246e5073de7871b500a572b4ab2d8f50eb87b695'
ORIGINAL=RUNTIME/'variants/whole_body_integrity_depletion/Gastrointestinal.cpp'
SOURCE_PIN='13af1a9350d9da3c4344b6cae187c7ecde44d11cb57358dc39ba7fe5cc1d333e'
VARIANT=RUNTIME/'variants/whole_body_integrity_gi_absorption'


def corrected_source(before):
    import hashlib
    if hashlib.sha256(before.encode()).hexdigest()!=SOURCE_PIN:raise ValueError('Held GI source changed')
    start=before.index('  // move nutrients, glucose',before.index('void Gastrointestinal::AbsorbNutrients()'))
    end=before.index('  // compute absorption rate as a function of volume',start)
    replacement='''  // Glucose and sodium share one bounded extent. Preserve the held 2:1
  // glucose/sodium mass ratio, including a finite tail or exact equality.
  const double glucoseSodiumExtent_g = std::max(0.0, std::min(sodiumAbsorbed_g,
    std::min(m_SmallIntestineChymeSodium->GetMass(MassUnit::g),
      m_SmallIntestineChymeGlucose->GetMass(MassUnit::g) / 2.0)));
  glucoseAbsorbed_g = 2.0 * glucoseSodiumExtent_g;
  if (glucoseSodiumExtent_g > 0.0) {
    m_SmallIntestineChymeSodium->GetMass().IncrementValue(-glucoseSodiumExtent_g, MassUnit::g);
    m_SmallIntestineChymeGlucose->GetMass().IncrementValue(-glucoseAbsorbed_g, MassUnit::g);
    m_SmallIntestineVascularSodium->GetMass().IncrementValue(glucoseSodiumExtent_g, MassUnit::g);
    m_smallIntestineVascularGlucose->GetMass().IncrementValue(glucoseAbsorbed_g, MassUnit::g);
  }

  // Preserve sequential glucose-before-AA priority: read the remaining live
  // sodium donor here, not the starting sodium already allocated above.
  aminoAcidAbsorbed_g = std::max(0.0, std::min(aminoAcidAbsorbed_g,
    std::min(m_SmallIntestineChymeSodium->GetMass(MassUnit::g),
      m_SmallIntestineChymeAminoAcids->GetMass(MassUnit::g))));
  if (aminoAcidAbsorbed_g > 0.0) {
    m_SmallIntestineChymeSodium->GetMass().IncrementValue(-aminoAcidAbsorbed_g, MassUnit::g);
    m_SmallIntestineChymeAminoAcids->GetMass().IncrementValue(-aminoAcidAbsorbed_g, MassUnit::g);
    m_SmallIntestineVascularSodium->GetMass().IncrementValue(aminoAcidAbsorbed_g, MassUnit::g);
    m_smallIntestineVascularAminoAcids->GetMass().IncrementValue(aminoAcidAbsorbed_g, MassUnit::g);
  }

  // Fat has one donor; consume at most its actual positive remainder.
  triacylglycerolAbsorbed_mg = std::max(0.0, std::min(triacylglycerolAbsorbed_mg,
    m_SmallIntestineChymeTriacylglycerol->GetMass(MassUnit::mg)));
  if (triacylglycerolAbsorbed_mg > 0.0) {
    m_SmallIntestineChymeTriacylglycerol->GetMass().IncrementValue(-triacylglycerolAbsorbed_mg, MassUnit::mg);
    m_smallintestineVAscularTriacylglycerol->GetMass().IncrementValue(triacylglycerolAbsorbed_mg, MassUnit::mg);
    m_smallintestineVAscularTriacylglycerol->Balance(BalanceLiquidBy::Mass);
    m_SmallIntestineChymeTriacylglycerol->Balance(BalanceLiquidBy::Mass);
  }

'''
    return before[:start]+replacement+before[end:]


def parent_inventory(cwd):
    manifest=json.loads((PARENT/'manifest.json').read_text())
    if sha(PARENT/'libbiogears.so.8.0.0')!=PARENT_PIN or manifest['library_sha256']!=PARENT_PIN:raise RuntimeError('Signed parent library changed')
    objects=shlex.split((PARENT/'objects.rsp').read_text())
    if len(objects)!=len(set(objects)) or set(objects)!=set(manifest['object_sha256']):raise RuntimeError('Parent inventory changed')
    for obj in objects:
        if sha(cwd/obj)!=manifest['object_sha256'][obj]:raise RuntimeError('Parent object changed: '+obj)
    return manifest,objects


def main():
    if sha(ORIGINAL)!=SOURCE_PIN:raise RuntimeError('Inherited GI source changed')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    parent,objects=parent_inventory(cwd);parent_manifest_sha=sha(PARENT/'manifest.json')
    before=ORIGINAL.read_text();after=corrected_source(before)
    VARIANT.mkdir(parents=True,exist_ok=False);print(VARIANT,flush=True)
    patched=VARIANT/'Gastrointestinal.cpp';patched.write_text(after)
    patch=VARIANT/'bounded_absorption.patch';patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='inherited/Gastrointestinal.cpp',tofile=VARIANT.name+'/Gastrointestinal.cpp')))
    obj=VARIANT/'Gastrointestinal.cpp.o'
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    def run(command,name):
        with (VARIANT/(name+'.log')).open('w') as log:
            subprocess.run(['/usr/bin/time','-v','-o',str(VARIANT/(name+'.resources'))]+command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,preexec_fn=limits,timeout=210)
    run(compile_command,'compile')
    original_obj=str(ORIGINAL)+'.o'
    if objects.count(original_obj)!=1:raise RuntimeError('Expected one inherited GI object')
    objects[objects.index(original_obj)]=str(obj)
    response=VARIANT/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());library=VARIANT/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(library)
    response_token='@CMakeFiles/libbiogears.dir/objects1.rsp'
    if command.count(response_token)!=1:raise RuntimeError('Unexpected link response structure')
    command=['@'+str(response) if arg==response_token else arg for arg in command];run(command,'link')
    parent_inventory(cwd)
    if sha(ORIGINAL)!=SOURCE_PIN or sha(PARENT/'manifest.json')!=parent_manifest_sha:raise RuntimeError('Parent/source changed during build')
    manifest={'variant':VARIANT.name,'source_revision':parent['source_revision'],'parent_variant':PARENT.name,
              'parent_manifest_sha256':parent_manifest_sha,'parent_library_sha256':PARENT_PIN,
              'header_sha256':parent['header_sha256'],'signed_port_abi_parent':PARENT.name,
              'original_source_sha256':SOURCE_PIN,'patched_source_sha256':sha(patched),'patch_sha256':sha(patch),
              'library_sha256':sha(library),'original_gi_object_sha256':parent['object_sha256'][original_obj],
              'object_sha256':{o:sha(cwd/o) for o in objects},'inherited_object_count':len(objects)-1,
              'compile_command':compile_command,'link_command':command,'builder_sha256':sha(Path(__file__)),
              'limits':{'address_space_bytes':4*1024**3,'cpu_seconds_per_child':180,'wall_seconds_per_child':210,'threads':1},
              'resources':{n:(VARIANT/(n+'.resources')).read_text() for n in ['compile','link']},
              'scope':'Bound glucose/AA paired and TAG absorption by live donor availability including equality. Held coefficients, sequential glucose-before-AA order and independent sodium/water policy unchanged. All signed-v2/shared-donor/thermal/substrate objects inherited.'}
    (VARIANT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ['variant','library_sha256','patched_source_sha256','inherited_object_count']},indent=2))


if __name__=='__main__':main()
