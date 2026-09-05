"""Paired native reaction and energy ledger for remaining glucose availability."""
import json,os,subprocess,tempfile,shutil
from pathlib import Path
from audit_native_substrate_integrity import ROOT,RUNTIME,SOURCE,PROBE,sha
PARENT='whole_body_integrity_evaporation_humidity';FIXED='whole_body_integrity_substrate_availability'

def main():
 out=Path(tempfile.mkdtemp(prefix='substrate-integrity-',dir=ROOT/'data/derived/audits'));build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
 probe=PROBE.replace('#include <cassert>','#include <cassert>\n#include <biogears/engine/Systems/Energy.h>\n#include <biogears/engine/Controller/BioGearsConfiguration.h>')
 probe=probe.replace('double glucose,double oxygen)', 'double glucose,double oxygen,double aa,double tag)')
 probe=probe.replace('set(t.m_Triacylglycerol,0);set(t.m_AminoAcids,0);','set(t.m_Triacylglycerol,tag);set(t.m_AminoAcids,aa);')
 probe=probe.replace(' t.CalculateMetabolicConsumptionAndProduction(.02);','''
 auto* vascular=t.m_TissueToVascular[t.m_LiverTissue];
 auto* insulin=vascular->GetSubstanceQuantity(*t.m_Insulin);insulin->GetMass().SetValue(0,MassUnit::g);insulin->Balance(BalanceLiquidBy::Mass);
 auto* glucagon=vascular->GetSubstanceQuantity(*t.m_Glucagon);glucagon->GetMass().SetValue(1e-5,MassUnit::g);glucagon->Balance(BalanceLiquidBy::Mass);
 t.m_energy->GetTotalMetabolicRate().Set(t.m_Patient->GetBasalMetabolicRate());
 t.m_energy->GetExerciseEnergyDemand().SetValue(0,PowerUnit::W);
 t.m_energy->GetEnergyDeficit().SetValue(0,PowerUnit::W);
 auto* urea=t.m_LiverExtracellular->GetSubstanceQuantity(*t.m_Urea);
 double ureaBefore=urea->GetMass(MassUnit::g)/t.m_Urea->GetMolarMass(MassPerAmountUnit::g_Per_mol);
 t.m_FatigueRunningAverage.Reset();
 t.CalculateMetabolicConsumptionAndProduction(.02);''')
 probe=probe.replace('t.m_CO2,t.m_Lactate}', 't.m_CO2,t.m_Lactate,t.m_AminoAcids,t.m_Triacylglycerol}')
 probe=probe.replace(' std::cout<<std::endl;', ''' std::cout<<','<<t.m_Patient->GetBasalMetabolicRate(PowerUnit::kcal_Per_s)*.02<<','<<t.m_energy->GetFatigueLevel().GetValue()<<','<<t.m_data.GetConfiguration().GetEnergyPerATP(EnergyPerAmountUnit::kcal_Per_mol)<<','<<(urea->GetMass(MassUnit::g)/t.m_Urea->GetMolarMass(MassPerAmountUnit::g_Per_mol)-ureaBefore)<<std::endl;''')
 probe=probe.replace('std::stod(argv[3]));','std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]));')
 cpp=out/'probe.cpp';cpp.write_text(probe);binary=out/'probe';cmd=['c++','-std=c++20','-O2',str(cpp)]
 for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:cmd+=['-I',str(p)]
 cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
 with (out/'compile.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 state=out/'frozen_initial.xml';shutil.copyfile(ROOT/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml',state)
 final=out/'frozen_exercise_3600.xml';shutil.copyfile(ROOT/'data/derived/systemic/exertion_v3/exercise/native/states/native_session.xml',final)
 cases=[('aerobic_depletion',1e-8,.01,0,0),('mixed_depletion',1e-8,3e-8,0,0),('anaerobic_only',1e-8,0,0,0),('adequate_glucose',.01,.01,0,0),('scarce_aa_tag',1e-8,.01,1e-9,1e-8),('scarce_oxygen_aa_tag',1e-8,1e-9,1e-9,1e-8)]
 results={};receipts={};checks={};multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip()
 for variant in [PARENT,FIXED]:
  work=out/variant;work.mkdir();selected=RUNTIME/'variants'/variant;env={**os.environ,'OPENBLAS_NUM_THREADS':'1','LD_LIBRARY_PATH':f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'}
  for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(work/name).symlink_to(build/'runtime'/name)
  linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(work/'ldd.txt').write_text(linkage);assert str(selected/'libbiogears.so.8.0.0') in linkage
  receipts[variant]={'library_sha256':sha(selected/'libbiogears.so.8.0.0'),'manifest_sha256':sha(selected/'manifest.json'),'dependency_sha256':{str(Path(s.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(s.split('=>')[1].split(' (')[0].strip())) for s in linkage.splitlines() if '=>' in s}}
  results[variant]={}
  for name,g,o,aa,tag in cases:
   r=subprocess.run([str(binary),str(state),str(g),str(o),str(aa),str(tag)],cwd=work,env=env,capture_output=True,text=True,timeout=120);(work/(name+'.log')).write_text(r.stdout+r.stderr);assert r.returncode==0
   line=next(l for l in r.stdout.splitlines() if l.startswith('SUBSTRATE,'));v=list(map(float,line.split(',')[1:]));initial,remaining,o_final,co2,lac,aa_final,tag_final,basal,fatigue,atp,urea_produced=v
   aa_used=aa-aa_final;tag_used=tag-tag_final;glucose_aerobic=(co2-1.5*aa_used-51*tag_used)/6;glucose_anaerobic=lac/2
   expected_o2=6*glucose_aerobic+1.875*aa_used+72.5*tag_used
   ledger_energy=glucose_aerobic*686+glucose_anaerobic*2*atp+aa_used*387.189+tag_used*7554
   native_energy=.8*basal-fatigue*basal
   row={'initial_mol':{'glucose':g,'oxygen':o,'aa':aa,'tag':tag},'final_mol':{'glucose':remaining,'oxygen':o_final,'co2':co2,'lactate':lac,'aa':aa_final,'tag':tag_final},'glucose_aerobic_mol':glucose_aerobic,'glucose_anaerobic_mol':glucose_anaerobic,'glucose_carbon_error_mol':6*(remaining+glucose_aerobic+glucose_anaerobic-g),'oxygen_stoichiometry_error_mol':o-o_final-expected_o2,'reaction_energy_kcal':ledger_energy,'native_accounted_energy_kcal':native_energy,'energy_error_kcal':native_energy-ledger_energy,'urea_produced_mol':urea_produced,'whole_reaction_carbon_error_mol':6*remaining+3*aa_final+51*tag_final+co2+3*lac+urea_produced-(6*g+3*aa+51*tag),'aa_carbon_model_omission_mol':aa_used,'atp_kcal_per_mol':atp}
   results[variant][name]=row;checks[variant+'/'+name+'/oxygen']=abs(row['oxygen_stoichiometry_error_mol'])<1e-12;checks[variant+'/'+name+'/energy']=abs(row['energy_error_kcal'])<1e-12
   if variant==FIXED:checks[variant+'/'+name+'/glucose_carbon']=abs(row['glucose_carbon_error_mol'])<1e-12
 checks['predecessor_carbon_defect_reproduced']=results[PARENT]['aerobic_depletion']['glucose_carbon_error_mol']>5.9e-8
 for name in ['anaerobic_only','adequate_glucose']:checks[name+'/unchanged']=results[PARENT][name]==results[FIXED][name]
 report={'passed':all(checks.values()),'checks':checks,'results':results,'receipts':receipts,'state_sha256':sha(state),'frozen_final_state_sha256':sha(final),'probe_sha256':sha(binary),'probe_source_sha256':sha(cpp),'compile_command':cmd,'scope':'Glucose reaction carbon and source stoichiometric O2/energy ledger. AA carbon model omission separately disclosed; not whole-body validation.'};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(out);print(json.dumps(checks,indent=2));assert report['passed']
if __name__=='__main__':main()
