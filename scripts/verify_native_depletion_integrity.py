"""Actual mixed-unit paired depletion and long native stomach-water replay."""
import argparse,hashlib,json,math,os,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUNTIME=ROOT/'data/runtime/physiology';SOURCE=ROOT/'data/raw/physiology/biogears'
PROBE=r'''
#include <cassert>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Systems/Gastrointestinal.h>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/cdm/patient/SENutrition.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iomanip>
#include <iostream>
namespace biogears {class BioGearsEngineTest {public:
 static void state(Gastrointestinal& g){std::cout<<','<<g.GetStomachContents().GetWater(VolumeUnit::mL)<<','<<g.m_WaterDigestionRate.GetValue(VolumePerTimeUnit::mL_Per_s)<<','<<g.m_dT_s;}
 static void direct(Gastrointestinal& g,const char* label,double water,const VolumeUnit& unit){
  auto& s=g.GetStomachContents();s.GetWater().SetValue(water,unit);s.GetSodium().SetValue(1,MassUnit::g);
  s.GetCarbohydrate().SetValue(0,MassUnit::g);s.GetProtein().SetValue(0,MassUnit::g);s.GetFat().SetValue(0,MassUnit::g);s.GetCalcium().SetValue(0,MassUnit::mg);
  g.m_DecrementNutrients=true;g.m_WaterDigestionRate.SetValue(.417,VolumePerTimeUnit::mL_Per_s);
  double initial=s.GetWater(VolumeUnit::mL),w=g.m_SmallIntestineChyme->GetVolume(VolumeUnit::mL),na=g.m_SmallIntestineChymeSodium->GetMass(MassUnit::g);bool failed=false;
  try{g.DigestNutrient();}catch(const std::exception& e){failed=true;}
  std::cout<<"DIRECT,"<<label<<','<<initial<<','<<s.GetWater(VolumeUnit::mL)<<','<<g.m_SmallIntestineChyme->GetVolume(VolumeUnit::mL)-w<<','<<s.GetSodium(MassUnit::g)<<','<<g.m_SmallIntestineChymeSodium->GetMass(MassUnit::g)-na<<','<<failed;
  bool second=false;try{g.DigestNutrient();}catch(const std::exception& e){second=true;}std::cout<<','<<second<<'\n';
 }
};}
int main(int argc,char**argv){using namespace biogears;std::cout<<std::setprecision(17);
 auto bg=CreateBioGearsEngine("probe.log");if(!bg->LoadState(argv[1]))return 2;
 auto* g=const_cast<Gastrointestinal*>(dynamic_cast<const Gastrointestinal*>(bg->GetGastrointestinalSystem()));
 if(std::string(argv[2])=="direct"){
  BioGearsEngineTest::direct(*g,"litre_remainder",1e-9,VolumeUnit::L);
  BioGearsEngineTest::direct(*g,"ml_remainder",1e-6,VolumeUnit::mL);
  BioGearsEngineTest::direct(*g,"litre_partial",.5,VolumeUnit::L);
  BioGearsEngineTest::direct(*g,"invalid_negative",-1e-12,VolumeUnit::L);
 }else{
  for(int k=0;k<60250;k++){
   try{if(!bg->AdvanceModelTime())return 3;}catch(const std::exception& e){std::cout<<"LONG,exception,"<<bg->GetSimulationTime(TimeUnit::s);BioGearsEngineTest::state(*g);std::cout<<','<<e.what()<<'\n';return 0;}
  }
  std::cout<<"LONG,complete,"<<bg->GetSimulationTime(TimeUnit::s);BioGearsEngineTest::state(*g);std::cout<<'\n';
 }
}
'''
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--parent-only',action='store_true');p.add_argument('--fixed-only',action='store_true');args=p.parse_args()
 if args.parent_only and args.fixed_only:p.error('choose at most one variant filter')
 out=Path(tempfile.mkdtemp(prefix='depletion-integrity-',dir=ROOT/'data/derived/audits'));cpp=out/'probe.cpp';cpp.write_text(PROBE);binary=out/'probe';build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
 cmd=['c++','-std=c++20','-O2',str(cpp)]
 for x in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:cmd+=['-I',str(x)]
 cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)];subprocess.run(cmd,check=True)
 state=ROOT/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml';multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip();results={};checks={};receipts={}
 for variant in ([] if args.fixed_only else ['whole_body_integrity_energy'])+([] if args.parent_only else ['whole_body_integrity_depletion']):
  selected=RUNTIME/'variants'/variant;env={**os.environ,'LD_LIBRARY_PATH':f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}','OPENBLAS_NUM_THREADS':'1'}
  linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);assert 'not found' not in linkage
  deps={str(Path(s.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(s.split('=>')[1].split(' (')[0].strip())) for s in linkage.splitlines() if '=>' in s};assert str((selected/'libbiogears.so.8.0.0').resolve()) in deps
  receipts[variant]={'dependency_sha256':deps,'state_sha256':sha(state),'probe_sha256':sha(cpp),'binary_sha256':sha(binary),'compile_command':cmd};results[variant]={}
  for mode in ['direct','long']:
   work=out/variant/mode;work.mkdir(parents=True)
   for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(work/name).symlink_to(build/'runtime'/name)
   run=subprocess.run([str(binary),str(state),mode],cwd=work,env=env,capture_output=True,text=True,timeout=600);(work/'stdout.log').write_text(run.stdout);(work/'stderr.log').write_text(run.stderr)
   assert run.returncode==0,work
   rows=[s.split(',') for s in run.stdout.splitlines() if s.startswith(('DIRECT,','LONG,'))];results[variant][mode]=rows
  direct={r[1]:list(map(float,r[2:])) for r in results[variant]['direct']};fixed=variant.endswith('_depletion')
  for label,(initial,remaining,credit,na,na_credit,failed,second_failed) in direct.items():
   if label=='invalid_negative':good=failed==second_failed==1 and remaining==initial and credit==na_credit==0 and na==1
   elif label=='litre_remainder' and not fixed:good=remaining<0 and second_failed==1 and math.isclose(credit,initial,abs_tol=1e-10)
   else:
    transfer=min(initial,.417*.02);good=failed==second_failed==0 and remaining>=0 and math.isclose(initial-remaining,transfer,abs_tol=1e-10) and math.isclose(credit,transfer,abs_tol=1e-10) and math.isclose(1-na,na_credit,abs_tol=1e-12)
    if transfer==initial:good=good and remaining==0
   checks[variant+'_'+label]=good
  r=results[variant]['long'][0];checks[variant+'_long_replay']=r[1]=='complete' and abs(float(r[2])-1205)<1e-6 and float(r[3])==0 if fixed else r[1]=='exception' and 1195<float(r[2])<1200 and float(r[3])<0 and float(r[4])==.417 and float(r[5])==.02
 report={'passed':all(checks.values()),'checks':checks,'receipts':receipts,'results':results};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(out);print(json.dumps(checks,indent=2));print(json.dumps({k:v['long'] for k,v in results.items()}))
 if not report['passed']:raise SystemExit(1)
if __name__=='__main__':main()
