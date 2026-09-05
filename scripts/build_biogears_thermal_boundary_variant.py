"""Isolate area-consistent dry thermal boundaries over the retained depletion layer."""
from pathlib import Path
import argparse,difflib,hashlib,json,re,shlex,subprocess
BASE=Path(__file__).resolve().parents[1]
RUNTIME=BASE/'data/runtime/physiology'
REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def build(variant_id='whole_body_integrity_thermal_boundary'):
    if not re.fullmatch(r'whole_body_integrity_thermal_boundary(?:_v[2-9][0-9]*)?',variant_id):raise ValueError('Versioned thermal boundary identity required')
    donor=BASE/'data/raw/physiology/biogears';original=donor/'projects/biogears/libBiogears/src/engine/Systems/Environment.cpp'
    committed=subprocess.check_output(['git','-C',str(donor),'show',f'{REVISION}:{original.relative_to(donor)}'])
    if original.read_bytes()!=committed:raise ValueError('Original donor Environment.cpp changed')
    parent=RUNTIME/'variants/whole_body_integrity_depletion';pm=json.loads((parent/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256']:raise ValueError('Depletion parent library changed')
    ancestor=RUNTIME/'variants/saturation_bounds_heatflux_thermal_units'
    am=json.loads((ancestor/'manifest.json').read_text());before=(ancestor/'Environment.cpp').read_text()
    if sha(ancestor/'Environment.cpp')!=am['environment_patched_sha256']:raise ValueError('Inherited thermal source changed')
    after='#include "native_thermal_boundary.h"\n'+before
    old='  double skinToClothingResistance = std::max(dClothingResistance_rsi / dSurfaceArea_m2, m_data.GetConfiguration().GetDefaultClosedHeatResistance(HeatResistanceUnit::K_Per_W));'
    new='  // Use the source evaporation area partition, whose fractions sum to one.\n  const double skinAreaFractions[] = {0.36, 0.07, 0.092, 0.092, 0.193, 0.193};'
    if after.count(old)!=1:raise ValueError('Expected source clothing setup once')
    after=after.replace(old,new)
    for i in range(6):
        old=f'skinToClothingResistance * m_cloSegmentation[{i}]'
        new=f'ihm_thermal::segment_resistance(dClothingResistance_rsi, dSurfaceArea_m2, skinAreaFractions[{i}], m_data.GetConfiguration().GetDefaultClosedHeatResistance(HeatResistanceUnit::K_Per_W))'
        if after.count(old)!=1:raise ValueError('Expected exactly one segment resistor')
        after=after.replace(old,new)
    # Replace the complete film resistance blocks. Keep the published coefficient
    # calculations and area, while removing the duplicate empirical clo factors.
    for coefficient,target,indent in [('Radiative','Enclosure','    '),('Convective','Environment','  ')]:
        marker=indent+'double dResistance_K_Per_W = 0.0;'
        method=after.index('void Environment::Calculate'+('Radiation' if coefficient=='Radiative' else 'Convection')+'()')
        start=after.index(marker,method)
        end=after.index('HeatResistanceUnit::K_Per_W);',after.index(f'm_ClothingTo{target}Path->GetNextResistance().SetValue(',start))+len('HeatResistanceUnit::K_Per_W);')
        replacement=indent+f'm_ClothingTo{target}Path->GetNextResistance().SetValue(ihm_thermal::film_resistance(dSurfaceArea_m2, d{coefficient}HeatTransferCoefficient_WPerM2_K, m_data.GetConfiguration().GetDefaultClosedHeatResistance(HeatResistanceUnit::K_Per_W), m_data.GetConfiguration().GetDefaultOpenHeatResistance(HeatResistanceUnit::K_Per_W)), HeatResistanceUnit::K_Per_W);'
        after=after[:start]+replacement+after[end:]
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(pm['object_sha256']) or any(sha(cwd/o)!=pm['object_sha256'][o] for o in objects):raise ValueError('Parent object inputs changed')
    variant=RUNTIME/'variants'/variant_id
    if variant.exists():raise ValueError('Thermal boundary build already retained; preserve or choose a new variant version')
    variant.mkdir();patched=variant/'Environment.cpp';patched.write_text(after)
    header=variant/'native_thermal_boundary.h';header.write_bytes((BASE/'scripts/native_thermal_boundary.h').read_bytes())
    patch=variant/'thermal_boundary.patch';patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='thermal_units/Environment.cpp',tofile='thermal_boundary/Environment.cpp')))
    obj=variant/'Environment.cpp.o'
    compile_cmd=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_cmd,cwd=cwd,check=True)
    selected=[i for i,o in enumerate(objects) if o==str(ancestor/'Environment.cpp.o')]
    if len(selected)!=1:raise ValueError('Inherited environment object is ambiguous')
    objects[selected[0]]=str(obj);response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());library=variant/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command];subprocess.run(command,cwd=cwd,check=True)
    manifest={'variant':variant.name,'source_revision':REVISION,'parent_variant':parent.name,'parent_library_sha256':pm['library_sha256'],
        'parent_manifest_sha256':sha(parent/'manifest.json'),'original_source_sha256':sha(original),'patched_source_sha256':sha(patched),
        'boundary_header_sha256':sha(header),'patch_sha256':sha(patch),'library_sha256':sha(library),'object_sha256':{o:sha(cwd/o) for o in objects},
        'compile_command':compile_cmd,'link_command':command,'build_cwd':str(cwd),'builder_sha256':sha(__file__),
        'scope':'Uniform source whole-body clothing Icl distributed by local source area; independent physical convection/radiation films. No fitting to temperature, no bedding or supine posture added. Evaporation and source control heuristics unchanged.'}
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if original.read_bytes()!=committed or sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256']:raise ValueError('Parent or donor changed during build')
    print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2));return variant
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--variant-id',default='whole_body_integrity_thermal_boundary');build(p.parse_args().variant_id)
