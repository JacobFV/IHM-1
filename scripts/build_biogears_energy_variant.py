#!/usr/bin/env python3
"""Compose bounded exercise demand with thermoregulation in an isolated native layer."""
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
    original=source/'projects/biogears/libBiogears/src/engine/Systems/Energy.cpp'
    committed=subprocess.check_output(['git','-C',str(source),'show',f'{REVISION}:{original.relative_to(source)}'])
    if original.read_bytes()!=committed:
        raise RuntimeError('Donor energy source differs from retained revision')
    before=original.read_text()
    after=before
    old='  Exercise();\n}'
    new='  Exercise();\n  m_temperatureGroundToCorePath->GetNextHeatSource().SetValue(GetTotalMetabolicRate(PowerUnit::W), PowerUnit::W);\n}'
    if after.count(old)!=1:raise RuntimeError('PreProcess layout changed')
    after=after.replace(old,new,1)
    old='  m_temperatureGroundToCorePath->GetNextHeatSource().SetValue(GetTotalMetabolicRate(PowerUnit::W), PowerUnit::W);'
    # Remove the old pre-exercise heat assignment, retaining the new final one.
    index=after.index(old,after.index('void Energy::CalculateMetabolicHeatGeneration()'))
    after=after[:index]+after[index:].replace(old,'  // Final thermal source is assigned after exercise demand is composed.',1)
    marker='void Energy::CalculateMetabolicHeatGeneration()\n{'
    partition=marker+"""
  // Total requested power is a sum of non-exercise and exercise components.
  // Remove the previous exercise component before updating thermoregulation.
  const double previousExercise_W = GetExerciseEnergyDemand(PowerUnit::W);
  const double nonExercise_W = GetTotalMetabolicRate(PowerUnit::W) - previousExercise_W;
  if (!std::isfinite(previousExercise_W) || previousExercise_W < 0.0 || !std::isfinite(nonExercise_W) || nonExercise_W < 0.0)
    throw CommonDataModelException("Inconsistent exercise power partition; legacy exercise states require a clean baseline");
  GetTotalMetabolicRate().SetValue(nonExercise_W, PowerUnit::W);
"""
    if after.count(marker)!=1:raise RuntimeError('Thermal source layout changed')
    after=after.replace(marker,partition,1)
    start=after.index('void Energy::Exercise()');end=after.index('//--------------------------------------------------------------------------------------------------',start)
    branch=after[start:end]
    no_action='    return;\n  }\n  exercise->GetGenericExercise()'
    if branch.count(no_action)!=1:raise RuntimeError('No-action branch changed')
    branch=branch.replace(no_action,'    GetExerciseEnergyDemand().SetValue(0.0, PowerUnit::W);\n    return;\n  }\n  exercise->GetGenericExercise()',1)
    old="""  const double TotalMetabolicRateSetPoint_kcal_Per_day = basalMetabolicRate_kcal_Per_day + (workRateDesired_W * kcal_Per_day_Per_Watt);
  const double exerciseEnergyIncrement_kcal_Per_day = MetabolicRateGain * (TotalMetabolicRateSetPoint_kcal_Per_day - currentMetabolicRate_kcal_Per_day);

  GetExerciseEnergyDemand().IncrementValue(exerciseEnergyIncrement_kcal_Per_day, PowerUnit::kcal_Per_day);
  GetTotalMetabolicRate().IncrementValue(exerciseEnergyIncrement_kcal_Per_day, PowerUnit::kcal_Per_day);"""
    new="""  // Preserve the source relaxation gain, applied to its own demand component.
  const double previousExercise_W = GetExerciseEnergyDemand(PowerUnit::W);
  const double exerciseDemand_W = previousExercise_W + MetabolicRateGain * (workRateDesired_W - previousExercise_W);
  GetExerciseEnergyDemand().SetValue(exerciseDemand_W, PowerUnit::W);
  GetTotalMetabolicRate().IncrementValue(exerciseDemand_W, PowerUnit::W);"""
    if branch.count(old)!=1:raise RuntimeError('Exercise update changed')
    branch=branch.replace(old,new,1)
    after=after[:start]+branch+after[end:]
    parent=RUNTIME/'variants/whole_body_integrity_gi_water'
    parent_manifest=json.loads((parent/'manifest.json').read_text())
    for filename,key in [('libbiogears.so.8.0.0','library_sha256')]:
        if sha(parent/filename)!=parent_manifest[key]:
            raise RuntimeError(f'Prior variant integrity failure: {filename}')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(parent_manifest['object_sha256']):
        raise RuntimeError('Parent response list changed')
    for obj in objects:
        if sha(cwd/obj)!=parent_manifest['object_sha256'][obj]:
            raise RuntimeError(f'Parent object changed: {obj}')
    variant=RUNTIME/'variants/whole_body_integrity_energy'
    variant.mkdir(parents=True,exist_ok=True)
    patched=variant/'Energy.cpp';patched.write_text(after)
    patch=variant/'energy_demand_partition.patch'
    patch_text=''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='original/Energy.cpp',tofile='whole_body_integrity_energy/Energy.cpp'))
    patch.write_text(patch_text)
    obj=variant/'Energy.cpp.o'
    compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    indices=[i for i,o in enumerate(objects) if o.endswith('/src/engine/Systems/Energy.cpp.o')]
    if len(indices)!=1:
        raise RuntimeError('Expected one engine energy object')
    original_obj=cwd/objects[indices[0]]
    objects[indices[0]]=str(obj)
    response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text())
    library=variant/'libbiogears.so.8.0.0'
    command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command]
    subprocess.run(command,cwd=cwd,check=True)
    manifest=dict(variant=variant.name,source_revision=revision,parent_variant=parent.name,parent_manifest_sha256=sha(parent/'manifest.json'),parent_library_sha256=parent_manifest['library_sha256'],original_source_sha256=sha(original),patched_source_sha256=sha(patched),patch_sha256=sha(patch),patch_text=patch_text,original_energy_object_sha256=sha(original_obj),library_sha256=sha(library),compile_command=compile_command,link_command=command,build_cwd=str(cwd),compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],object_sha256={o:sha(cwd/o) for o in objects},scope='Independent source-gain exercise power relaxation, explicit total/non-exercise partition before thermal update, clear demand when inactive, final requested heat source after composition. No tissue stoichiometry or clinical gains changed. Reject incompatible historical power partition.')
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if original.read_bytes()!=committed or sha(parent/'libbiogears.so.8.0.0')!=parent_manifest['library_sha256']:
        raise RuntimeError('Donor or parent changed during build')
    print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2))

if __name__=='__main__':
    build()
