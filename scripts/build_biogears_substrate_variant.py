#!/usr/bin/env python3
"""Read actual remaining glucose mass before anaerobic metabolism."""
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
    original=source/'projects/biogears/libBiogears/src/engine/Systems/Tissue.cpp'
    committed=subprocess.check_output(['git','-C',str(source),'show',f'{REVISION}:{original.relative_to(source)}'])
    if original.read_bytes()!=committed:raise RuntimeError('Donor Tissue source changed')
    before=original.read_text()
    wrong='    totalEnergyAsTissueIntracellularGlucose_kcal = (TissueGlucose->GetMolarity(AmountPerVolumeUnit::mol_Per_L) * TissueVolume_L) * anaerobic_ATP_Per_Glucose * energyPerMolATP_kcal;'
    right='    totalEnergyAsTissueIntracellularGlucose_kcal = (TissueGlucose->GetMass(MassUnit::g) / m_Glucose->GetMolarMass(MassPerAmountUnit::g_Per_mol)) * anaerobic_ATP_Per_Glucose * energyPerMolATP_kcal;'
    if before.count(wrong)!=1:raise RuntimeError('Expected exact anaerobic glucose availability once')
    after=before.replace(wrong,right,1)
    parent=RUNTIME/'variants/whole_body_integrity_evaporation_humidity'
    parent_manifest=json.loads((parent/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=parent_manifest['library_sha256']:raise RuntimeError('RH parent changed')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(parent_manifest['object_sha256']):
        raise RuntimeError('Parent response list changed')
    for obj in objects:
        if sha(cwd/obj)!=parent_manifest['object_sha256'][obj]:
            raise RuntimeError(f'Parent object changed: {obj}')
    variant=RUNTIME/'variants/whole_body_integrity_substrate_availability'
    variant.mkdir(parents=True,exist_ok=True)
    patched=variant/'Tissue.cpp';patched.write_text(after)
    patch=variant/'remaining_glucose_availability.patch'
    patch_text=''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='upstream/Tissue.cpp',tofile='whole_body_integrity_substrate_availability/Tissue.cpp'))
    patch.write_text(patch_text)
    obj=variant/'Tissue.cpp.o'
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    indices=[i for i,o in enumerate(objects) if o=='CMakeFiles/libbiogears.dir/src/engine/Systems/Tissue.cpp.o']
    if len(indices)!=1:
        raise RuntimeError('Expected one inherited native Tissue object')
    original_obj=cwd/objects[indices[0]]
    objects[indices[0]]=str(obj)
    response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library=variant/'libbiogears.so.8.0.0'
    command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command]
    subprocess.run(command,cwd=cwd,check=True)
    manifest=dict(variant=variant.name,source_revision=revision,parent_variant=parent.name,parent_manifest_sha256=sha(parent/'manifest.json'),parent_library_sha256=parent_manifest['library_sha256'],original_source_sha256=sha(original),patched_source_sha256=sha(patched),patch_sha256=sha(patch),patch_text=patch_text,original_tissue_object_sha256=sha(original_obj),library_sha256=sha(library),compile_command=compile_command,link_command=command,build_cwd=str(cwd),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],object_sha256={o:sha(cwd/o) for o in objects},scope='Anaerobic glucose availability uses actual remaining mass after aerobic consumption; no rates, stoichiometry, gains or clamps changed. All RH predecessor objects preserved.')
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if original.read_bytes()!=committed or sha(parent/'libbiogears.so.8.0.0')!=parent_manifest['library_sha256']:
        raise RuntimeError('Donor or parent changed during build')
    print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2))

if __name__=='__main__':
    build()
