"""Conserve whole-body sweat latent power across native regional skin areas."""
from pathlib import Path
import difflib,json,shlex,subprocess
from build_biogears_thermal_boundary_variant import BASE,RUNTIME,REVISION,sha

def build():
    donor=BASE/'data/raw/physiology/biogears'
    producer=donor/'projects/biogears/libBiogears/src/engine/Systems/Energy.cpp'
    if producer.read_bytes()!=subprocess.check_output(['git','-C',str(donor),'show',f'{REVISION}:{producer.relative_to(donor)}']):raise ValueError('Original sweat producer changed')
    for proof in ('double massLost_kg = sweatRate_kg_Per_s * m_dT_s;',
                  'm_Patient->GetWeight().IncrementValue(-massLost_kg, MassUnit::kg);',
                  'GetSweatRate().SetValue(sweatRate_kg_Per_s, MassPerTimeUnit::kg_Per_s);',
                  'm_skinExtravascularToSweatingGroundPath->GetNextFlowSource().SetValue(sweatRate_mL_Per_s, VolumePerTimeUnit::mL_Per_s);'):
        if producer.read_text().count(proof)!=1:raise ValueError('Whole-body sweat producer identity changed')
    parent=RUNTIME/'variants/whole_body_integrity_skin_perfusion';pm=json.loads((parent/'manifest.json').read_text())
    ancestor=RUNTIME/'variants/whole_body_integrity_thermal_boundary_v2';am=json.loads((ancestor/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256'] or sha(ancestor/'Environment.cpp')!=am['patched_source_sha256']:raise ValueError('Inherited thermal variant changed')
    before=(ancestor/'Environment.cpp').read_text();after='#include "native_sweat_evaporation.h"\n'+before
    start=after.index('        double dSweatingControlMechanisms = dSweatRate_kgPers')
    end=after.index('        index += 1;',start)
    replacement='''        // SweatRate is whole-body kg/s, as shown by the single mass-loss
        // and SkinSweating fluid-path producer. Distribute its latent power
        // once by area and limit evaporation to the local vapor capacity.
        const auto evaporation = ihm_thermal::evaporation_budget(
            dSweatRate_kgPers * m_dHeatOfVaporizationOfWater_J_Per_kg,
            m_Patient->GetSkinSurfaceArea(AreaUnit::m2),
            segmentedSkinSurfaceAreaPercents[index],
            dMaxEvaporativePotential, skinWettednessDiffusion);
        // Unevaporated sweat already left the fluid circuit; it supplies no
        // latent cooling. Clothing liquid storage/runoff is not modeled here.
        envSkinToGround->GetNextHeatSource().SetValue(evaporation.total_w, PowerUnit::W);
'''
    after=after[:start]+replacement+after[end:]
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(pm['object_sha256']) or any(sha(cwd/o)!=pm['object_sha256'][o] for o in objects):raise ValueError('Inherited object inventory changed')
    variant=RUNTIME/'variants/whole_body_integrity_sweat_evaporation'
    if variant.exists():raise ValueError('Preserve existing evaporation variant')
    variant.mkdir();patched=variant/'Environment.cpp';patched.write_text(after)
    for name in ('native_thermal_boundary.h','native_sweat_evaporation.h'):
        source=ancestor/name if name=='native_thermal_boundary.h' else BASE/'scripts'/name
        (variant/name).write_bytes(source.read_bytes())
    patch=variant/'sweat_evaporation.patch';patch.write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='thermal_boundary_v2/Environment.cpp',tofile='sweat_evaporation/Environment.cpp')))
    obj=variant/'Environment.cpp.o';compile_command=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)]
    subprocess.run(compile_command,cwd=cwd,check=True)
    selected=[i for i,o in enumerate(objects) if o==str(ancestor/'Environment.cpp.o')]
    if len(selected)!=1:raise ValueError('Inherited environment object ambiguous')
    objects[selected[0]]=str(obj);response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());library=variant/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(library)
    command=['@'+str(response) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command];subprocess.run(command,cwd=cwd,check=True)
    manifest=dict(variant=variant.name,source_revision=REVISION,parent_variant=parent.name,parent_library_sha256=pm['library_sha256'],
        parent_manifest_sha256=sha(parent/'manifest.json'),patched_source_sha256=sha(patched),patch_sha256=sha(patch),library_sha256=sha(library),
        object_sha256={o:sha(cwd/o) for o in objects},compile_command=compile_command,link_command=command,
        producer_source_path=str(producer.relative_to(BASE)),producer_source_sha256=sha(producer),builder_sha256=sha(__file__),
        header_sha256={n:sha(variant/n) for n in ('native_thermal_boundary.h','native_sweat_evaporation.h')},
        scope='Whole-body sweat latent power distributed by area once, local evaporative-capacity and wettedness bounds; source perfusion/dry-boundary fixes inherited. Source humidity omission, regional vapor-pressure and garment-permeability assumptions remain unchanged. No sweat mass/controller change, wet garment storage or condensation model.')
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2));return variant

if __name__=='__main__':build()
