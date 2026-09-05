#!/usr/bin/env python3
"""Source-preserving native aqueous sodium integrity layer over calcium/renal fixes."""
from pathlib import Path
import difflib
import hashlib
import json
import shlex
import subprocess

BASE=Path(__file__).resolve().parents[1]
RUNTIME=BASE/'data/runtime/physiology'
REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def build():
    donor=BASE/'data/raw/physiology/biogears'
    revision=subprocess.check_output(['git','-C',str(donor),'rev-parse','HEAD'],text=True).strip()
    if revision!=REVISION:raise RuntimeError('Unexpected donor revision')
    original=donor/'projects/biogears/libBiogears/src/engine/Systems/Gastrointestinal.cpp'
    committed=subprocess.check_output(['git','-C',str(donor),'show',f'{REVISION}:{original.relative_to(donor)}'])
    if original.read_bytes()!=committed:raise RuntimeError('Donor GI bytes changed')
    parent=RUNTIME/'variants/whole_body_integrity_renal'
    ancestor=RUNTIME/'variants/whole_body_integrity'
    manifest=json.loads((parent/'manifest.json').read_text())
    ancestor_manifest=json.loads((ancestor/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=manifest['library_sha256']:raise RuntimeError('Parent library changed')
    if sha(ancestor/'manifest.json')!=manifest['parent_manifest_sha256']:raise RuntimeError('Parent calcium lineage changed')
    parent_source=ancestor/'Gastrointestinal.cpp'
    if sha(parent_source)!=ancestor_manifest['patched_source_sha256']:raise RuntimeError('Parent calcium GI source changed')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    objects=shlex.split((parent/'objects.rsp').read_text())
    if len(objects)!=len(set(objects)) or set(objects)!=set(manifest['object_sha256']):raise RuntimeError('Parent object list changed')
    for obj in objects:
        if sha(cwd/obj)!=manifest['object_sha256'][obj]:raise RuntimeError(f'Parent object changed: {obj}')
    before=parent_source.read_text()
    needle='void Gastrointestinal::DigestNutrient()\n{\n'
    if before.count(needle)!=1:raise RuntimeError('GI entry point not unique')
    preflight='''  // Unknown water is not a dry, measured-zero pool. Reject invalid inputs
  // before any nutrient/ion mutation; do not create carrier fluid or infer zero.
  if (!m_StomachContents->HasWater())
    throw CommonDataModelException("GI digestion requires known stomach water");
  const double waterContent_mL = m_StomachContents->GetWater(VolumeUnit::mL);
  const double waterRate_mL_Per_s = m_WaterDigestionRate.GetValue(VolumePerTimeUnit::mL_Per_s);
  if (!std::isfinite(waterContent_mL) || waterContent_mL < 0.0
      || !std::isfinite(waterRate_mL_Per_s) || waterRate_mL_Per_s < 0.0
      || !std::isfinite(m_dT_s) || m_dT_s <= 0.0)
    throw CommonDataModelException("GI digestion requires finite nonnegative water and rate, and positive timestep");
  if (m_StomachContents->HasSodium()) {
    const double sodium_g = m_StomachContents->GetSodium(MassUnit::g);
    if (!std::isfinite(sodium_g) || sodium_g < 0.0)
      throw CommonDataModelException("GI digestion requires nonnegative finite sodium or an explicitly unknown pool");
  }
  const double waterToDigest_mL = waterRate_mL_Per_s * m_dT_s;
  if (!std::isfinite(waterToDigest_mL))
    throw CommonDataModelException("GI requested water transfer is not finite");

'''
    after=before.replace(needle,needle+preflight)
    after=after.replace('#include <biogears/engine/Systems/Gastrointestinal.h>','#include <biogears/engine/Systems/Gastrointestinal.h>\n#include <cmath>',1)
    start=after.index('  // sodium\n',after.index(needle))
    end=after.index('  digestedAmount = waterToDigest_mL < waterContent_mL',start)
    replacement='''  // Sodium follows existing aqueous transfer only when carrier water exists.
  // Retain hydrated kinetics in paired grams using a bounded carrier fraction.
  if (m_StomachContents->HasSodium()) {
    totalCheck = m_StomachContents->GetSodium(MassUnit::g);
    double digestedNa_g = 0.0;
    if (waterContent_mL > 0.0 && totalCheck > 0.0 && waterToDigest_mL > 0.0) {
      const double carrierFraction = std::min(1.0, waterToDigest_mL / waterContent_mL);
      digestedNa_g = totalCheck * carrierFraction;
    }
    if (digestedNa_g > 0.0) {
      if (m_DecrementNutrients) {
        if (digestedNa_g >= totalCheck) {
          // This exact paired transfer establishes a known zero, not missing data.
          m_StomachContents->GetSodium().SetValue(0.0, MassUnit::g);
          Info("Stomach is out of Sodium");
        } else {
          m_StomachContents->GetSodium().IncrementValue(-digestedNa_g, MassUnit::g);
        }
      }
#ifdef logDigest
      m_ss << "Digested " << digestedNa_g << "(g) of Sodium";
      Info(m_ss);
#endif
      m_SmallIntestineChymeSodium->GetMass().IncrementValue(digestedNa_g, MassUnit::g);
      // Preserve source balancing after the paired water update.
    }
  }

'''
    after=after[:start]+replacement+after[end:]
    variant=RUNTIME/'variants/whole_body_integrity_gi_water';variant.mkdir(parents=True,exist_ok=True)
    patched=variant/'Gastrointestinal.cpp';patched.write_text(after)
    patch_text=''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='whole_body_integrity/Gastrointestinal.cpp',tofile='whole_body_integrity_gi_water/Gastrointestinal.cpp'))
    patch=variant/'dry_gi_sodium.patch';patch.write_text(patch_text)
    obj=variant/'Gastrointestinal.cpp.o'
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    indices=[i for i,o in enumerate(objects) if o==str(ancestor/'Gastrointestinal.cpp.o')]
    if len(indices)!=1:raise RuntimeError('Expected exactly one calcium-corrected GI object')
    replaced=objects[indices[0]];objects[indices[0]]=str(obj)
    response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library=variant/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command]
    subprocess.run(command,cwd=cwd,check=True)
    for name,digest in manifest['object_sha256'].items():
        if sha(cwd/name)!=digest:raise RuntimeError('Parent object mutated during isolated build')
    if original.read_bytes()!=committed or sha(parent/'libbiogears.so.8.0.0')!=manifest['library_sha256']:raise RuntimeError('Donor or parent library changed during build')
    result=dict(variant=variant.name,source_revision=revision,parent_variant=parent.name,parent_manifest_sha256=sha(parent/'manifest.json'),
                parent_library_sha256=manifest['library_sha256'],parent_object_sha256=manifest['object_sha256'],
                original_source_sha256=sha(original),parent_gi_source_sha256=sha(parent_source),replaced_parent_gi_object_sha256=sha(replaced),
                patched_source_sha256=sha(patched),patch_sha256=sha(patch),patch_text=patch_text,library_sha256=sha(library),
                compile_command=compile_command,link_command=command,build_cwd=str(cwd),
                compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],object_sha256={o:sha(cwd/o) for o in objects},
                scope='Zero carrier water causes zero sodium transfer; positive transfers capped and paired in g. Verified sodium depletion remains numeric zero, unknown sodium remains unknown. Invalid water/rate/sodium fails before GI mutation. Hydrated kinetics, native stabilization reservoir mode and all parent corrections retained. No new biological coefficients or secretion process.')
    (variant/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['parent_object_sha256','object_sha256','link_command']},indent=2))
    return result

if __name__=='__main__':build()
