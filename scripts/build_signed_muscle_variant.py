#!/usr/bin/env python3
"""Layer narrow signed muscle hooks over an immutable, inventoried native variant."""
from pathlib import Path
import argparse,difflib,json,shlex,subprocess
from build_biogears_substrate_variant import BASE,RUNTIME,REVISION,sha

def once(text,old,new):
    if text.count(old)!=1:raise RuntimeError('Source anchor not unique: '+old[:100])
    return text.replace(old,new,1)

def patch_source(name,before):
    after='#include "native_signed_muscle_port.h"\n'+before
    if name=='Energy':
        after='#define IHM_SIGNED_MUSCLE_IMPLEMENTATION\n'+after
        # The inherited energy correction moved the final source write to
        # the end of PreProcess, after Exercise; bind to that exact write.
        source_write='  m_temperatureGroundToCorePath->GetNextHeatSource().SetValue(GetTotalMetabolicRate(PowerUnit::W), PowerUnit::W);'
        after=once(after,source_write,source_write+'''
  if (ihm_signed::active()) {
    if (m_data.GetActions().GetPatientActions().HasExercise())
      throw std::runtime_error("signed muscle port conflicts with generic Exercise");
    const double coupledHeat_W = ihm_signed::heat(m_temperatureGroundToCorePath->GetNextHeatSource().GetValue(PowerUnit::W));
    if (ihm_signed::active()->delta_h_W != 0)
      m_temperatureGroundToCorePath->GetNextHeatSource().SetValue(coupledHeat_W, PowerUnit::W);
  }''')
    elif name=='Tissue':
        anchor='  for (SETissueCompartment* tissue : m_ConsumptionProdutionTissues) {\n    vascular = m_TissueToVascular[tissue];'
        preflight='''  // Read-only muscle domain preflight before any tissue reaction in this call.
  double signedMuscleIncrement_kcal = 0;
  if (ihm_signed::active()) {
    auto* mv = m_TissueToVascular[m_MuscleTissue];
    auto& mi = m_data.GetCompartments().GetIntracellularFluid(*m_MuscleTissue);
    const double fraction = mv->HasInFlow() && totalFlowRate_mL_Per_min > 0
      ? mv->GetInFlow(VolumePerTimeUnit::mL_Per_min) / totalFlowRate_mL_Per_min : 0;
    const double mandatory = mandatoryMuscleAnaerobicFraction * nonbrainNeededEnergy_kcal * fraction;
    const double discretionary = nonbrainNeededEnergy_kcal * fraction - mandatory
      + 0.8 * exerciseEnergyRequested_kcal + 0.5 * otherEnergyDemandAboveBasal_kcal;
    double hormone = Hepatic::CalculateRelativeHormoneChange(GetLiverInsulinSetPoint().GetValue(AmountPerVolumeUnit::mmol_Per_L) * 1e9, GetLiverGlucagonSetPoint().GetValue(MassPerVolumeUnit::mg_Per_mL) * 1e9, mv->GetSubstanceQuantity(*m_Insulin), mv->GetSubstanceQuantity(*m_Glucagon), m_data);
    BLIM(hormone, -2, 0);
    const double requestedAA = (GeneralMath::LinearInterpolator(0, 2, 15, 110, -hormone) * time_s * fraction) / (24 * 3600 * m_AminoAcids->GetMolarMass(MassPerAmountUnit::g_Per_mol));
    const double availableAA = mi.GetSubstanceQuantity(*m_AminoAcids)->GetMolarity(AmountPerVolumeUnit::mol_Per_L) * mi.GetVolume(VolumeUnit::L);
    const double oxygenAA = mi.GetSubstanceQuantity(*m_O2)->GetMass(MassUnit::g) / m_O2->GetMolarMass(MassPerAmountUnit::g_Per_mol) / O2_Per_AA;
    const double pendingAA = std::min(std::min(requestedAA, availableAA), oxygenAA) * ATP_Per_AA * energyPerMolATP_kcal / AA_CellularEfficiency;
    signedMuscleIncrement_kcal = ihm_signed::budget(discretionary, pendingAA, mandatory, time_s, Convert(1.0, PowerUnit::W, PowerUnit::kcal_Per_s));
    if (signedMuscleIncrement_kcal != 0) totalEnergyRequested_kcal += signedMuscleIncrement_kcal;
  }

'''
        after=once(after,anchor,preflight+anchor)
        anchor='      double creatinineProductionRate_mg_Per_s = 2.0e-5;'
        after=once(after,anchor,'      if (signedMuscleIncrement_kcal != 0) tissueNeededEnergy_kcal += signedMuscleIncrement_kcal;\n\n'+anchor)
        anchor='    // Balance everything\n    TissueO2->Balance'
        after=once(after,anchor,'''    if (ihm_signed::active() && tissue == m_MuscleTissue)
      ihm_signed::active()->unmet_muscle_kcal = tissueNeededEnergy_kcal > 0 ? tissueNeededEnergy_kcal : 0;

'''+anchor)
    else:
        old='m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::W)'
        if name=='Cardiovascular':
            after=once(after,'m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::kcal_Per_day)',
                'ihm_signed::effective(m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::kcal_Per_day), Convert(1.0, PowerUnit::W, PowerUnit::kcal_Per_day), 0)')
        elif name=='Endocrine':after=once(after,old,'ihm_signed::effective('+old+', 1.0, 1)')
        elif name=='Nervous':
            if after.count(old)!=2:raise RuntimeError('Expected two Nervous demand readers')
            # Replace full expressions once each, so replacement cannot nest.
            after=after.replace(old,'IHM_EFFECTIVE_PLACEHOLDER',2)
            after=after.replace('IHM_EFFECTIVE_PLACEHOLDER','ihm_signed::effective('+old+', 1.0, 2)',1)
            after=after.replace('IHM_EFFECTIVE_PLACEHOLDER','ihm_signed::effective('+old+', 1.0, 3)',1)
    return after

def build(parent_id,variant_id):
    donor=BASE/'data/raw/physiology/biogears';cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    parent=RUNTIME/'variants'/parent_id;pm=json.loads((parent/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256']:raise RuntimeError('Parent library changed')
    objects=shlex.split((parent/'objects.rsp').read_text())
    if set(objects)!=set(pm['object_sha256']) or any(sha(cwd/o)!=pm['object_sha256'][o] for o in objects):raise RuntimeError('Parent inventory changed')
    variant=RUNTIME/'variants'/variant_id
    if variant.exists():raise RuntimeError('Preserve existing signed variant; choose a fresh name')
    sources={};patched={}
    for name in ('Energy','Tissue','Cardiovascular','Endocrine','Nervous'):
        original=donor/f'projects/biogears/libBiogears/src/engine/Systems/{name}.cpp'
        if original.read_bytes()!=subprocess.check_output(['git','-C',str(donor),'show',f'{REVISION}:{original.relative_to(donor)}']):raise RuntimeError('Original changed')
        matches=[o for o in objects if o.endswith('/'+name+'.cpp.o')]
        if len(matches)!=1:raise RuntimeError('Ambiguous inherited '+name)
        obj=cwd/matches[0]
        source=obj.with_suffix('') if obj.is_absolute() and 'variants' in obj.parts else original
        if not source.exists():raise RuntimeError('Missing inherited source '+str(source))
        if source != original and sha(source) not in (source.parent/'manifest.json').read_text():raise RuntimeError('Inherited source not bound by its manifest: '+str(source))
        sources[name]=(source,matches[0]);patched[name]=patch_source(name,source.read_text())
    variant.mkdir(parents=True)
    header=variant/'native_signed_muscle_port.h';header.write_bytes((BASE/'scripts/native_signed_muscle_port.h').read_bytes())
    cmds=[];source_manifest={}
    for name,(source,inherited_obj) in sources.items():
        path=variant/(name+'.cpp');path.write_text(patched[name]);obj=variant/(name+'.cpp.o')
        patch=variant/(name+'.patch');patch.write_text(''.join(difflib.unified_diff(source.read_text().splitlines(True),patched[name].splitlines(True),fromfile=str(source),tofile=str(path))))
        cmd=['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(path),'-o',str(obj)]
        subprocess.run(cmd,cwd=cwd,check=True);cmds.append(cmd)
        objects[objects.index(inherited_obj)]=str(obj)
        source_manifest[name]={'inherited_source':str(source),'inherited_source_sha256':sha(source),'patched_source_sha256':sha(path),'patch_sha256':sha(patch),'inherited_object':inherited_obj}
    response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    cmd=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());lib=variant/'libbiogears.so.8.0.0';cmd[cmd.index('-o')+1]=str(lib)
    cmd=['@'+str(response) if x=='@CMakeFiles/libbiogears.dir/objects1.rsp' else x for x in cmd];subprocess.run(cmd,cwd=cwd,check=True)
    manifest={'variant':variant_id,'source_revision':REVISION,'parent_variant':parent_id,'parent_manifest_sha256':sha(parent/'manifest.json'),'parent_library_sha256':pm['library_sha256'],'library_sha256':sha(lib),'header_sha256':sha(header),'builder_sha256':sha(Path(__file__)),'sources':source_manifest,'compile_commands':cmds,'link_command':cmd,'object_sha256':{o:sha(cwd/o) for o in objects},'scope':'Diagnostic signed muscle-only demand and independent heat; immutable basal/stress writers; fail-close on rejection, no native rollback claim.'}
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:v for k,v in manifest.items() if k!='object_sha256'},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',default='whole_body_integrity_shared_donor');p.add_argument('--variant',default='whole_body_integrity_signed_muscle_v2');a=p.parse_args();build(a.parent,a.variant)
