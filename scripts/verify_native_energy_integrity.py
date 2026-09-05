"""Actual native cold/neutral branch and matched whole-body energy probes."""
import hashlib,json,math,os,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUNTIME=ROOT/'data/runtime/physiology';SOURCE=ROOT/'data/raw/physiology/biogears'
PROBE=r'''
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Systems/Energy.h>
#include <biogears/engine/Systems/Tissue.h>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/cdm/circuit/thermal/SEThermalCircuitNode.h>
#include <biogears/cdm/circuit/thermal/SEThermalCircuitPath.h>
#include <biogears/cdm/patient/actions/SEExercise.h>
#include <biogears/cdm/scenario/SEPatientActionCollection.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iostream>
#include <cmath>
#include <iomanip>
namespace biogears {
class BioGearsEngineTest {
public:
 static void reset(Energy& e,double temp) {
   e.m_coreNode->GetTemperature().SetReadOnly(false);
   e.m_coreNode->GetTemperature().SetValue(temp,TemperatureUnit::C);
   e.m_coreNode->GetTemperature().SetReadOnly(true);
   e.GetExerciseEnergyDemand().SetValue(0,PowerUnit::W);
   e.GetTotalMetabolicRate().Set(e.m_Patient->GetBasalMetabolicRate());
 }
 static void emit(Energy& e,const SETissueSystem& t,double seconds) {
   std::cout<<"ENERGY,"<<std::setprecision(17)<<seconds<<','<<e.GetExerciseEnergyDemand(PowerUnit::W)<<','<<e.GetTotalMetabolicRate(PowerUnit::W)<<','<<e.m_temperatureGroundToCorePath->GetNextHeatSource(PowerUnit::W)<<','<<e.m_PatientActions->HasExercise()<<','<<t.GetOxygenConsumptionRate(VolumePerTimeUnit::mL_Per_min)<<','<<t.GetMuscleGlycogen(MassUnit::g)<<','<<t.GetLiverGlycogen(MassUnit::g)<<','<<t.GetStoredFat(MassUnit::g)<<','<<t.GetStoredProtein(MassUnit::g)<<','<<e.GetFatigueLevel()<<'\n';
 }
};}
int main(int argc,char**argv){
 using namespace biogears;
 auto bg=CreateBioGearsEngine("probe.log");if(!bg->LoadState(argv[1]))return 2;
 auto* e=const_cast<Energy*>(dynamic_cast<const Energy*>(bg->GetEnergySystem()));
 auto* t=bg->GetTissueSystem();
 auto action=[&](double intensity){SEExercise::SEGeneric g;g.Intensity.SetValue(intensity);SEExercise a{g};if(!bg->ProcessAction(a))throw std::runtime_error("Action rejected");};
 std::string mode=argv[2];
 if(mode=="legacy"){
  BioGearsEngineTest::reset(*e,37.0);e->GetExerciseEnergyDemand().SetValue(1000,PowerUnit::W);
  try{e->PreProcess();std::cout<<"LEGACY,accepted\n";}catch(const CommonDataModelException& ex){std::cout<<"LEGACY,rejected\n";}
  return 0;
 }
 if(mode=="cold"||mode=="neutral"){
  BioGearsEngineTest::reset(*e,mode=="cold"?36.5:37.0);action(0);e->PreProcess();BioGearsEngineTest::emit(*e,*t,0);
  action(.15);for(int k=1;k<=500;k++){e->PreProcess();if(k==1||k==500)BioGearsEngineTest::emit(*e,*t,k*.02);}
  action(0);e->PreProcess();BioGearsEngineTest::emit(*e,*t,10.02);
 }else{
  action(mode=="rest"?0:.15);BioGearsEngineTest::emit(*e,*t,0);
  double previous=0;for(double now:{.02,1.,5.,10.,30.,30.02,31.,60.,120.}){for(int k=0;k<std::llround((now-previous)/.02);++k)if(!bg->AdvanceModelTime())return 7;BioGearsEngineTest::emit(*e,*t,now);if(now==30&&mode=="stop")action(0);previous=now;}
 }
}
'''
FIELDS=['time_s','exercise_demand_w','total_metabolic_w','thermal_next_source_w','exercise_action_active','oxygen_ml_min','muscle_glycogen_g','liver_glycogen_g','stored_fat_g','stored_protein_g','fatigue']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=Path(tempfile.mkdtemp(prefix='energy-integrity-',dir=ROOT/'data/derived/audits'));cpp=out/'probe.cpp';cpp.write_text(PROBE);binary=out/'probe'
 build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
 command=['c++','-std=c++20','-O2',str(cpp)]
 for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:command+=['-I',str(p)]
 command+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)];subprocess.run(command,check=True)
 state=ROOT/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml';multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip()
 results={};checks={};receipts={}
 for variant in ['whole_body_integrity_gi_water','whole_body_integrity_energy']:
  selected=RUNTIME/'variants'/variant;env={**os.environ,'OPENBLAS_NUM_THREADS':'1','LD_LIBRARY_PATH':f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'}
  linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True)
  deps={str(Path(line.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(line.split('=>')[1].split(' (')[0].strip())) for line in linkage.splitlines() if '=>' in line and 'not found' not in line}
  assert str((selected/'libbiogears.so.8.0.0').resolve()) in deps
  receipts[variant]={'dependency_sha256':deps,'manifest_sha256':sha(selected/'manifest.json'),'state_sha256':sha(state),'probe_sha256':sha(binary),'source_sha256':sha(cpp),'compile_command':command}
  results[variant]={}
  for mode in ['cold','neutral','rest','continuous','stop']:
   work=out/variant/mode;work.mkdir(parents=True)
   for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(work/name).symlink_to(build/'runtime'/name)
   run=subprocess.run([str(binary),str(state),mode],cwd=work,env=env,capture_output=True,text=True,timeout=300)
   (work/'stdout.log').write_text(run.stdout);(work/'stderr.log').write_text(run.stderr)
   rows=[dict(zip(FIELDS,map(float,line.split(',')[1:]))) for line in run.stdout.splitlines() if line.startswith('ENERGY,')]
   if run.returncode or len(rows)!=(4 if mode in ['cold','neutral'] else 10):raise RuntimeError(f'Native probe failed {variant}/{mode}; {work}')
   results[variant][mode]=rows
   checks[variant+'_'+mode+'_finite']=all(math.isfinite(v) for row in rows for v in row.values())
  fixed=variant.endswith('_energy')
  work=out/variant/'cold'
  legacy=subprocess.run([str(binary),str(state),'legacy'],cwd=work,env=env,capture_output=True,text=True)
  (out/variant/'legacy.log').write_text(legacy.stdout+legacy.stderr)
  checks[variant+'_legacy_partition_guard']=legacy.returncode==0 and ('LEGACY,rejected' in legacy.stdout if fixed else 'LEGACY,accepted' in legacy.stdout)
  for mode in ['cold','neutral']:
   r=results[variant][mode];baseline=r[0]['total_metabolic_w']
   good=all(abs(row['exercise_demand_w']-150.12*(1-.98**(round(row['time_s']/.02))))<1e-6 and abs(row['total_metabolic_w']-baseline-row['exercise_demand_w'])<1e-6 and abs(row['thermal_next_source_w']-row['total_metabolic_w'])<1e-9 for row in r[1:3]) and r[-1]['exercise_demand_w']==0 and abs(r[-1]['total_metabolic_w']-baseline)<1e-6
   checks[variant+'_'+mode+'_partition_expected']=good if fixed else not good
  r=results[variant]
  checks[variant+'_matched_initial']=r['rest'][0]==r['continuous'][0] # action flag differs; handled below
  checks[variant+'_matched_initial']=all(r['rest'][0][key]==r['continuous'][0][key]==r['stop'][0][key] for key in FIELDS if key!='exercise_action_active')
  checks[variant+'_matched_before_stop']=r['continuous'][:6]==r['stop'][:6]
  if fixed:
   checks['fixed_bound']=all(0<=v['exercise_demand_w']<=150.1200001 for mode in r.values() for v in mode)
   checks['fixed_stop_zero']=all(v['exercise_demand_w']==0 and v['exercise_action_active']==0 for v in r['stop'][6:])
   checks['fixed_stores_nonnegative']=all(v[k]>=0 for mode in r.values() for v in mode for k in ['muscle_glycogen_g','liver_glycogen_g','stored_fat_g','stored_protein_g'])
   checks['fixed_exercise_causal_oxygen']=r['continuous'][-1]['oxygen_ml_min']>r['rest'][-1]['oxygen_ml_min']
   checks['fixed_stop_less_glycogen_use']=r['stop'][-1]['muscle_glycogen_g']>r['continuous'][-1]['muscle_glycogen_g']
 checks['rest_parity']=results['whole_body_integrity_gi_water']['rest']==results['whole_body_integrity_energy']['rest']
 report={'passed':all(checks.values()),'checks':checks,'receipts':receipts,'results':results};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(out);print(json.dumps(checks,indent=2))
 if not report['passed']:raise SystemExit(1)
if __name__=='__main__':main()
