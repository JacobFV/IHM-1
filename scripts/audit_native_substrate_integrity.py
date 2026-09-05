"""Read-only native substrate branch reproduction; does not patch engine equations."""
import hashlib,json,os,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUNTIME=ROOT/'data/runtime/physiology';SOURCE=ROOT/'data/raw/physiology/biogears'
PROBE=r'''
#include <cassert>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Systems/Tissue.h>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/engine/Controller/BioGearsCompartments.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/compartment/fluid/SELiquidCompartment.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iostream>
#include <iomanip>
namespace biogears { class BioGearsEngineTest { public:
static void run(Tissue& t,double glucose,double oxygen) {
 t.m_ConsumptionProdutionTissues={t.m_LiverTissue};
 auto& c=t.m_data.GetCompartments().GetIntracellularFluid(*t.m_LiverTissue);
 auto set=[&](SESubstance* s,double mol){auto*q=c.GetSubstanceQuantity(*s);q->GetMass().SetValue(mol*s->GetMolarMass(MassPerAmountUnit::g_Per_mol),MassUnit::g);q->Balance(BalanceLiquidBy::Mass);};
 set(t.m_Glucose,glucose);set(t.m_O2,oxygen);set(t.m_CO2,0);set(t.m_Lactate,0);set(t.m_Triacylglycerol,0);set(t.m_AminoAcids,0);set(t.m_Ketones,0);
 t.CalculateMetabolicConsumptionAndProduction(.02);
 std::cout<<"SUBSTRATE,"<<std::setprecision(17)<<glucose;
 for(auto*s:{t.m_Glucose,t.m_O2,t.m_CO2,t.m_Lactate})std::cout<<','<<c.GetSubstanceQuantity(*s)->GetMass(MassUnit::g)/s->GetMolarMass(MassPerAmountUnit::g_Per_mol);
 std::cout<<std::endl;
}};}
int main(int argc,char**argv){auto bg=biogears::CreateBioGearsEngine("probe.log");if(!bg->LoadState(argv[1]))return 2;auto*t=const_cast<biogears::Tissue*>(dynamic_cast<const biogears::Tissue*>(bg->GetTissueSystem()));biogears::BioGearsEngineTest::run(*t,std::stod(argv[2]),std::stod(argv[3]));}
'''
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=Path(tempfile.mkdtemp(prefix='substrate-source-',dir=ROOT/'data/derived/audits'));cpp=out/'probe.cpp';cpp.write_text(PROBE);binary=out/'probe';build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
 cmd=['c++','-std=c++20','-O2',str(cpp)]
 for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:cmd+=['-I',str(p)]
 cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
 with (out/'compile.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 reference=ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json';manifest=json.loads(reference.read_text());selected=Path(manifest['library_path']).parent;multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip();env={**os.environ,'LD_LIBRARY_PATH':f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}','OPENBLAS_NUM_THREADS':'1'}
 state=ROOT/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'
 for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(out/name).symlink_to(build/'runtime'/name)
 results={}
 for name,glucose,oxygen in [('aerobic_depletion',1e-8,.01),('mixed_depletion',1e-8,3e-8),('anaerobic_only',1e-8,0),('adequate_glucose',.01,.01)]:
  r=subprocess.run([str(binary),str(state),str(glucose),str(oxygen)],cwd=out,env=env,capture_output=True,text=True,timeout=120);(out/(name+'_stdout.log')).write_text(r.stdout);(out/(name+'_stderr.log')).write_text(r.stderr)
  rows=[line for line in r.stdout.splitlines() if line.startswith('SUBSTRATE,')];assert r.returncode==0 and len(rows)==1
  v=list(map(float,rows[0].split(',')[1:]));results[name]={'glucose_initial_mol':v[0],'glucose_final_mol':v[1],'oxygen_initial_mol':oxygen,'oxygen_final_mol':v[2],'co2_final_mol':v[3],'lactate_final_mol':v[4],'carbon_initial_mol':6*v[0],'carbon_final_mol':6*v[1]+v[3]+3*v[4],'carbon_created_mol':6*v[1]+v[3]+3*v[4]-6*v[0]}
 linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
 record={'exit_code':r.returncode,'compile_command':cmd,'source_sha256':{p.name:sha(p) for p in (SOURCE/'projects/biogears/libBiogears/src/engine/Systems').glob('*.cpp') if p.name in ['Tissue.cpp','Hepatic.cpp','BloodChemistry.cpp','Saturation.cpp']},'state_sha256':sha(state),'reference_manifest_sha256':sha(reference),'probe_sha256':sha(binary),'probe_source_sha256':sha(cpp),'dependency_sha256':{str(Path(s.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(s.split('=>')[1].split(' (')[0].strip())) for s in linkage.splitlines() if '=>' in s}}
 record['cases']=results
 record['known_defect_reproduced']=results['aerobic_depletion']['carbon_created_mol']>5.9e-8 and results['mixed_depletion']['carbon_created_mol']>2.9e-8 and abs(results['anaerobic_only']['carbon_created_mol'])<1e-15 and abs(results['adequate_glucose']['carbon_created_mol'])<1e-12
 (out/'report.json').write_text(json.dumps(record,indent=2)+'\n');print(out);print(json.dumps(record,indent=2));assert record['known_defect_reproduced']
if __name__=='__main__':main()
