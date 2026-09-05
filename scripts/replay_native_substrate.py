"""Short paired continuation of frozen native state; no registry/default changes."""
import argparse,json,os,subprocess,tempfile,shutil
from pathlib import Path
from audit_native_substrate_integrity import ROOT,RUNTIME,SOURCE,sha
from verify_native_substrate_integrity import PARENT,FIXED
PROBE=r'''
#include <cassert>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Systems/Tissue.h>
#include <biogears/engine/Systems/Energy.h>
#include <biogears/engine/Systems/BloodChemistry.h>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/engine/Controller/BioGearsCompartments.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iomanip>
#include <iostream>
namespace biogears {class BioGearsEngineTest {public:static void emit(Tissue&t,int time){auto*a=t.m_data.GetCompartments().GetLiquidCompartment("Aorta");std::cout<<"REPLAY,"<<std::setprecision(17)<<time<<','<<a->GetSubstanceQuantity(*t.m_Glucose)->GetConcentration(MassPerVolumeUnit::mg_Per_dL)<<','<<a->GetSubstanceQuantity(*t.m_Lactate)->GetConcentration(MassPerVolumeUnit::mg_Per_dL)<<','<<t.m_data.GetBloodChemistry().GetArterialBloodPH().GetValue()<<','<<t.m_data.GetEnergy().GetTotalMetabolicRate(PowerUnit::W)<<','<<t.GetMuscleGlycogen(MassUnit::g)<<std::endl;}};}
int main(int argc,char**argv){auto bg=biogears::CreateBioGearsEngine("probe.log");if(!bg->LoadState(argv[1]))return 2;auto*t=const_cast<biogears::Tissue*>(dynamic_cast<const biogears::Tissue*>(bg->GetTissueSystem()));biogears::BioGearsEngineTest::emit(*t,0);for(int i=1;i<=3000;++i){if(!bg->AdvanceModelTime())return 7;if(i%500==0)biogears::BioGearsEngineTest::emit(*t,i/50);}bg->SaveStateToFile("final.xml");}
'''
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--state',type=Path,default=ROOT/'data/derived/systemic/exertion_v3/exercise/native/states/native_session.xml');args=parser.parse_args()
 out=Path(tempfile.mkdtemp(prefix='substrate-replay-',dir=ROOT/'data/derived/audits'));state=out/'frozen_state.xml';shutil.copyfile(args.state,state);cpp=out/'probe.cpp';cpp.write_text(PROBE);binary=out/'probe';build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
 cmd=['c++','-std=c++20','-O2',str(cpp)]
 for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:cmd+=['-I',str(p)]
 cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
 with (out/'compile.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 results={};receipts={};multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip()
 for variant in [PARENT,FIXED]:
  work=out/variant;work.mkdir();selected=RUNTIME/'variants'/variant;env={**os.environ,'OPENBLAS_NUM_THREADS':'1','LD_LIBRARY_PATH':f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'}
  for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(work/name).symlink_to(build/'runtime'/name)
  linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(work/'ldd.txt').write_text(linkage)
  r=subprocess.run([str(binary),str(state)],cwd=work,env=env,capture_output=True,text=True,timeout=180);(work/'stdout.log').write_text(r.stdout);(work/'stderr.log').write_text(r.stderr);assert r.returncode==0
  results[variant]=[list(map(float,l.split(',')[1:])) for l in r.stdout.splitlines() if l.startswith('REPLAY,')];assert len(results[variant])==7
  receipts[variant]={'exit_code':r.returncode,'library_sha256':sha(selected/'libbiogears.so.8.0.0'),'final_state_sha256':sha(work/'states/final.xml'),'dependency_sha256':{str(Path(s.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(s.split('=>')[1].split(' (')[0].strip())) for s in linkage.splitlines() if '=>' in s}}
 assert results[PARENT][0]==results[FIXED][0]
 report={'scope':'60s native continuation of identical frozen3600s exercise state; not correction of prior history or one-hour validation','fields':['relative_time_s','aortic_glucose_mg_dl','aortic_lactate_mg_dl','arterial_ph','metabolic_rate_w','muscle_glycogen_g'],'results':results,'receipts':receipts,'state_sha256':sha(state),'probe_sha256':sha(binary),'source_sha256':sha(cpp),'compile_command':cmd};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(out);print(json.dumps(results,indent=2))
if __name__=='__main__':main()
