"""Use ambient relative humidity once and actual local skin vapor pressure."""
from pathlib import Path
import difflib,json,shlex,subprocess
from build_biogears_thermal_boundary_variant import BASE,RUNTIME,REVISION,sha

def build():
    parent=RUNTIME/'variants/whole_body_integrity_sweat_evaporation';pm=json.loads((parent/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256'] or sha(parent/'Environment.cpp')!=pm['patched_source_sha256']:raise ValueError('Sweat-area parent changed')
    before=(parent/'Environment.cpp').read_text()
    # SupplementalValues stores saturation pressure in this field. Its separate
    # local dPressureOfWaterVapor_Pa applies RH for air density without changing
    # the field; applying RH at the evaporative consumer is therefore once only.
    for proof in ('m_dWaterVaporPressureInAmbientAir_Pa = Convert(dWaterVaporPressureInAmbientAir_mmHg, PressureUnit::mmHg, PressureUnit::Pa);',
                  'double dPressureOfWaterVapor_Pa = GetConditions().GetRelativeHumidity().GetValue() * m_dWaterVaporPressureInAmbientAir_Pa;'):
        if before.count(proof)!=1:raise ValueError('Ambient saturation/actual vapor producer proof changed')
    old='        double dMaxEvaporativePotential = (1.0 / 1000.0) * (m_dWaterVaporPressureAtSkin_Pa - m_dWaterVaporPressureInAmbientAir_Pa) / (dClothingResistance_m2_kPa_Per_W + 1.0 / (fCl * dEvaporativeHeatTransferCoefficient_W_Per_m2_kPa));'
    if before.count(old)!=1:raise ValueError('Expected single evaporative vapor-gradient consumer')
    new='''        // This path begins at the actual regional external skin thermal
        // node; do not reuse the torso saturation pressure for every region.
        const double regionalSkin_C = envSkinToGround->GetSourceNode().GetTemperature(TemperatureUnit::C);
        const double regionalSkinSaturation_Pa = Convert(GeneralMath::AntoineEquation(regionalSkin_C), PressureUnit::mmHg, PressureUnit::Pa);
        const double ambientActualVapor_Pa = GetConditions().GetRelativeHumidity().GetValue() * m_dWaterVaporPressureInAmbientAir_Pa;
        double dMaxEvaporativePotential = (1.0 / 1000.0) * (regionalSkinSaturation_Pa - ambientActualVapor_Pa) / (dClothingResistance_m2_kPa_Per_W + 1.0 / (fCl * dEvaporativeHeatTransferCoefficient_W_Per_m2_kPa));'''
    after=before.replace(old,new)
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(pm['object_sha256']) or any(sha(cwd/o)!=pm['object_sha256'][o] for o in objects):raise ValueError('Inherited objects changed')
    variant=RUNTIME/'variants/whole_body_integrity_evaporation_humidity'
    if variant.exists():raise ValueError('Preserve existing humidity variant')
    variant.mkdir();patched=variant/'Environment.cpp';patched.write_text(after)
    for name in pm['header_sha256']:
        if sha(parent/name)!=pm['header_sha256'][name]:raise ValueError('Inherited helper changed')
        (variant/name).write_bytes((parent/name).read_bytes())
    patch=variant/'evaporation_humidity.patch';patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='sweat_evaporation/Environment.cpp',tofile='evaporation_humidity/Environment.cpp')))
    obj=variant/'Environment.cpp.o';compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    selected=[i for i,o in enumerate(objects) if o==str(parent/'Environment.cpp.o')]
    if len(selected)!=1:raise ValueError('Inherited environment object ambiguous')
    objects[selected[0]]=str(obj);response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());library=variant/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command];subprocess.run(command,cwd=cwd,check=True)
    manifest=dict(variant=variant.name,source_revision=REVISION,parent_variant=parent.name,parent_library_sha256=pm['library_sha256'],parent_manifest_sha256=sha(parent/'manifest.json'),
        patched_source_sha256=sha(patched),patch_sha256=sha(patch),library_sha256=sha(library),object_sha256={o:sha(cwd/o) for o in objects},
        compile_command=compile_command,link_command=command,header_sha256=pm['header_sha256'],builder_sha256=sha(__file__),
        scope='Evaporative vapor gradient uses ambient RH times saturation pressure once, and each actual regional skin temperature. Dry boundary, sweat-area/capacity and skin-perfusion corrections inherited; no parameter fitting or garment water storage. Source regional garment permeability remains an assumption.')
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2));return variant

if __name__=='__main__':build()
