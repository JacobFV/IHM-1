#!/usr/bin/env python3
"""Layer an isolated paired calcium transfer correction over retained physics fixes."""
from pathlib import Path
import difflib
import hashlib
import json
import shlex
import subprocess

BASE = Path(__file__).resolve().parents[1]
RUNTIME = BASE/'data/runtime/physiology'
REVISION = '3f16a5fa1dade9c511b88d923606fa51cc35e95d'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    source = BASE/'data/raw/physiology/biogears'
    revision = subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    if revision != REVISION:
        raise RuntimeError('Unexpected donor revision')
    original = source/'projects/biogears/libBiogears/src/engine/Systems/Gastrointestinal.cpp'
    relative = str(original.relative_to(source))
    committed = subprocess.check_output(['git','-C',str(source),'show',f'{REVISION}:{relative}'])
    if original.read_bytes() != committed:
        raise RuntimeError('Donor GI source differs from retained revision')
    before = original.read_text()
    start = before.index('  // calcium\n', before.index('void Gastrointestinal::DigestNutrient()'))
    end = before.index('  // sodium\n',start)
    replacement = '''  // calcium: a local paired transfer in mg; preserve source absorption-fraction scaling.
  if (m_StomachContents->HasCalcium()) {
    totalCheck = std::max(m_StomachContents->GetCalcium(MassUnit::mg), 0.0);
    digestedAmount = m_CalciumDigestionRate.GetValue(MassPerTimeUnit::mg_Per_s) * m_dT_s;
    digestedAmount *= m_data.GetConfiguration().GetCalciumAbsorptionFraction();
    digestedAmount = std::min(totalCheck, std::max(digestedAmount, 0.0));
    if (digestedAmount > 0.0) {
#ifdef logDigest
      m_ss << "Digested " << digestedAmount << "(mg) of Calcium";
      Info(m_ss);
#endif
      m_StomachContents->GetCalcium().IncrementValue(-digestedAmount, MassUnit::mg);
      m_SmallIntestineChymeCalcium->GetMass().IncrementValue(digestedAmount, MassUnit::mg);
      m_SmallIntestineChymeCalcium->Balance(BalanceLiquidBy::Mass);
    }
  }

'''
    after = before[:start]+replacement+before[end:]
    parent = RUNTIME/'variants/saturation_bounds_heatflux_thermal_units'
    parent_manifest = json.loads((parent/'manifest.json').read_text())
    for filename, key in [('libbiogears.so.8.0.0','library_sha256'),('Saturation.cpp','patched_source_sha256'),('Environment.cpp','environment_patched_sha256')]:
        if sha(parent/filename) != parent_manifest[key]:
            raise RuntimeError(f'Prior variant integrity failure: {filename}')
    variant = RUNTIME/'variants/whole_body_integrity'
    variant.mkdir(parents=True,exist_ok=True)
    patched = variant/'Gastrointestinal.cpp'
    patched.write_text(after)
    patch_text = ''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='original/Gastrointestinal.cpp',tofile='whole_body_integrity/Gastrointestinal.cpp'))
    patch = variant/'gi_calcium.patch'
    patch.write_text(patch_text)
    cwd = RUNTIME/'biogears-build/projects/biogears/libBiogears'
    obj = variant/'Gastrointestinal.cpp.o'
    compile_command = ['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    # Start from the actual parent object list, preserving BOTH prior replacements.
    objects = shlex.split((parent/'objects.rsp').read_text())
    for filename in ['Saturation.cpp.o','Environment.cpp.o']:
        matches = [o for o in objects if o == str(parent/filename) or o.endswith('/src/engine/Systems/'+filename)]
        if matches != [str(parent/filename)]:
            raise RuntimeError(f'Parent replacement object missing: {filename}')
    gi_indices = [i for i,o in enumerate(objects) if o.endswith('/src/engine/Systems/Gastrointestinal.cpp.o')]
    if len(gi_indices) != 1:
        raise RuntimeError('Expected exactly one original GI object')
    original_obj = cwd/objects[gi_indices[0]]
    objects[gi_indices[0]] = str(obj)
    response = variant/'objects.rsp'
    response.write_text(' '.join(shlex.quote(o) for o in objects))
    command = shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library = variant/'libbiogears.so.8.0.0'
    command[command.index('-o')+1] = str(library)
    command = ['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command]
    subprocess.run(command,cwd=cwd,check=True)
    manifest = dict(variant=variant.name,source_revision=revision,parent_variant=parent.name,parent_manifest_sha256=sha(parent/'manifest.json'),parent_library_sha256=parent_manifest['library_sha256'],original_source_sha256=sha(original),patched_source_sha256=sha(patched),patch_sha256=sha(patch),patch_text=patch_text,original_gi_object_sha256=sha(original_obj),library_sha256=sha(library),compile_command=compile_command,link_command=command,build_cwd=str(cwd),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],object_sha256={o:sha(cwd/o) for o in objects},scope='Paired stomach-to-chyme calcium transfer in mg, capped by nonnegative availability and including positive remainder. Absorption fraction and all non-calcium equations unchanged. Existing negative input is preserved without credit, not silently repaired.')
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if original.read_bytes() != committed:
        raise RuntimeError('Donor source changed during build')
    print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2))

if __name__=='__main__':
    build()
