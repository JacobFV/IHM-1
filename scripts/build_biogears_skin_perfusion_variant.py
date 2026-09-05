"""Partition the cardiac-cycle mean total skin inflow once across thermal regions."""
from pathlib import Path
import difflib,json,shlex,subprocess
from build_biogears_thermal_boundary_variant import BASE,RUNTIME,REVISION,sha

def build():
    donor=BASE/'data/raw/physiology/biogears';energy=donor/'projects/biogears/libBiogears/src/engine/Systems/Energy.cpp'
    cardiovascular=donor/'projects/biogears/libBiogears/src/engine/Systems/Cardiovascular.cpp'
    for path in (energy,cardiovascular):
        if path.read_bytes()!=subprocess.check_output(['git','-C',str(donor),'show',f'{REVISION}:{path.relative_to(donor)}']):raise ValueError('Donor changed')
    producer=cardiovascular.read_text()
    for proof in ('m_pAortaToSkin = m_CirculatoryCircuit->GetPath(BGE::CardiovascularPath::Aorta1ToSkin1)',
                  'double SkinFlow_mL_Per_s = m_pAortaToSkin->GetNextFlow(VolumePerTimeUnit::mL_Per_s)',
                  'm_CardiacCycleSkinFlow_mL_Per_s.Sample(SkinFlow_mL_Per_s)',
                  'GetMeanSkinFlow().SetValue(m_CardiacCycleSkinFlow_mL_Per_s.Value(), VolumePerTimeUnit::mL_Per_s)'):
        if producer.count(proof)!=1:raise ValueError('Temporal total-skin-flow producer proof changed')
    ancestor=RUNTIME/'variants/whole_body_integrity_energy';am=json.loads((ancestor/'manifest.json').read_text());before=(ancestor/'Energy.cpp').read_text()
    if sha(ancestor/'Energy.cpp')!=am['patched_source_sha256']:raise ValueError('Inherited energy source changed')
    begin=before.index('void Energy::UpdateHeatResistance()');end=before.index('  int index = 0;',begin)
    block=before[begin:end]
    old='6.0 * skinBloodFlow_m3_Per_s'
    if block.count(old)!=6:raise ValueError('Expected six regional flow multipliers')
    corrected=block.replace(old,'skinBloodFlow_m3_Per_s').replace('// mean * 6 compts * compartmental fraction for each segments blood flow','// MeanSkinFlow is the cardiac-cycle temporal mean of the total Aorta1ToSkin1 inflow.\n  // Regional fractions sum to one; do not multiply a whole-skin flow by six.')
    after=before[:begin]+corrected+before[end:]
    parent=RUNTIME/'variants/whole_body_integrity_thermal_boundary_v2';pm=json.loads((parent/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256']:raise ValueError('Thermal parent library changed')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(pm['object_sha256']) or any(sha(cwd/o)!=pm['object_sha256'][o] for o in objects):raise ValueError('Parent object inventory changed')
    variant=RUNTIME/'variants/whole_body_integrity_skin_perfusion'
    if variant.exists():raise ValueError('Preserve existing skin-perfusion variant')
    variant.mkdir();patched=variant/'Energy.cpp';patched.write_text(after);patch=variant/'skin_perfusion_partition.patch'
    patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='energy_partition/Energy.cpp',tofile='skin_perfusion/Energy.cpp')))
    obj=variant/'Energy.cpp.o';compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    indices=[i for i,o in enumerate(objects) if o==str(ancestor/'Energy.cpp.o')]
    if len(indices)!=1:raise ValueError('Inherited energy object ambiguous')
    objects[indices[0]]=str(obj);response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());library=variant/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command];subprocess.run(command,cwd=cwd,check=True)
    manifest={'variant':variant.name,'source_revision':REVISION,'parent_variant':parent.name,'parent_library_sha256':pm['library_sha256'],
        'parent_manifest_sha256':sha(parent/'manifest.json'),'patched_source_sha256':sha(patched),'patch_sha256':sha(patch),
        'library_sha256':sha(library),'object_sha256':{o:sha(cwd/o) for o in objects},'compile_command':compile_command,'link_command':command,
        'producer_source_path':str(cardiovascular.relative_to(BASE)),'producer_source_sha256':sha(cardiovascular),'original_energy_sha256':sha(energy),
        'builder_sha256':sha(__file__),'scope':'Remove sixfold duplication of total skin perfusion when partitioning thermal regions. Regional proportions and empirical alpha=.42 retained; no clothing, thermoregulation gain, patient or environment fitting.'}
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2))
if __name__=='__main__':build()
