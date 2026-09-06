#!/usr/bin/env python3
"""Small compiled port boundary fixtures; native-source fixtures are separate."""
from pathlib import Path
import subprocess, tempfile
BASE=Path(__file__).resolve().parents[1]
FIXTURE=r'''
#define IHM_SIGNED_MUSCLE_IMPLEMENTATION
#include "native_signed_muscle_port.h"
#include <cassert>
int main(){
 using namespace ihm_signed;
 Record r; r.sequence=1;r.start_s=0;r.end_s=.02;r.reference_id="frozen";
 r.delta_m_W=-1;r.delta_w_W=-2;r.delta_h_W=1;
 validate(r); Scope scope(&r);
 assert(effective(80,1,0)==79);
 assert(heat(80)==81);
 double delta=budget(.001,.0001,.000028,.02,1.0/4184);
 assert(delta==-1*.02/4184); assert(r.allowed_decrement_W>1);
 finish(r);assert(r.heat_count==1&&r.tissue_count==1);
 bool rejected=false;try {heat(80);}catch(const std::runtime_error&){rejected=true;}assert(rejected);
 Record bad;bad.sequence=2;bad.start_s=.02;bad.end_s=.04;bad.reference_id="frozen";
 bad.delta_m_W=-1000;bad.delta_h_W=-1000;validate(bad);
 {Scope other(&bad);rejected=false;try{budget(.001,.0001,.000028,.02,1.0/4184);}catch(const std::runtime_error&){rejected=true;}assert(rejected);assert(bad.tissue_count==0);}
 Record zero;zero.sequence=3;zero.start_s=.04;zero.end_s=.06;zero.reference_id="frozen";validate(zero);
 {Scope other(&zero);assert(heat(80)==80);assert(budget(.001,.0001,.000028,.02,1.0/4184)==0);finish(zero);}
}
'''
def main():
    assert (BASE/'scripts/native_signed_muscle_port.h').exists(), 'signed muscle port missing'
    with tempfile.TemporaryDirectory(prefix='ihm-signed-fixture-') as tmp:
        src=Path(tmp)/'fixture.cpp';src.write_text(FIXTURE)
        exe=Path(tmp)/'fixture'
        subprocess.run(['c++','-std=c++20','-Wall','-Wextra','-Werror','-I'+str(BASE/'scripts'),str(src),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)
    print('PASS compiled signed port boundary fixtures (not native Tissue reaction closure)')

NATIVE_FIXTURE=r'''
#include <cassert>
#include "native_signed_muscle_port.h"
#include "native_body_ports.h"
#include <biogears/engine/Systems/Tissue.h>
#include <biogears/cdm/substance/SESubstance.h>
#include <biogears/engine/Systems/Energy.h>
#include <biogears/cdm/circuit/thermal/SEThermalCircuit.h>
#include <biogears/cdm/circuit/thermal/SEThermalCircuitPath.h>
#include <iostream>
#include <iomanip>
#include <memory>
struct TissueAccess : biogears::Tissue {
 static auto reaction(){return &TissueAccess::CalculateMetabolicConsumptionAndProduction;}
};
int main(int argc,char** argv){
 if(argc!=5 && argc!=6)return 2;
 auto bg=std::make_unique<BioGearsEngine>("fixture.log");
 if(!bg->LoadState(argv[1]))return 3;
 bg->SetAutoTrackFlag(false);
 if(argc==6 && std::string(argv[5])=="scarce") {
  auto* muscle=bg->GetCompartments().GetLiquidCompartment("MuscleTissueIntracellular");
  for(const std::string sub:{"Glucose","Triacylglycerol","Oxygen"}) {
   auto* q=muscle->GetSubstanceQuantity(*bg->GetSubstanceManager().GetSubstance(sub));
   q->GetMass().SetValue(0,MassUnit::g);q->Balance(BalanceLiquidBy::Mass);
  }
  bg->GetTissue().GetMuscleGlycogen().SetValue(0,MassUnit::g);
 }
 ihm_signed::Record r;r.sequence=1;r.start_s=0;r.end_s=.02;r.reference_id="native-frozen-fixture";
 r.delta_m_W=std::stod(argv[2]);r.delta_w_W=std::stod(argv[3]);r.delta_h_W=r.delta_m_W-r.delta_w_W;
 const bool enabled=std::string(argv[4])=="1";
 std::map<std::string,double> out;
 auto observe=[&](const std::string& phase){
  out[phase+".tmr"]=bg->GetEnergy().GetTotalMetabolicRate(PowerUnit::W);
  out[phase+".glycogen_g"]=bg->GetTissue().GetMuscleGlycogen(MassUnit::g);
  for(const std::string organ:{"Muscle","Brain","Fat"}){
   auto* c=bg->GetCompartments().GetLiquidCompartment(organ+"TissueIntracellular");
   for(const std::string sub:{"Glucose","AminoAcids","Triacylglycerol","Oxygen","CarbonDioxide","Lactate"})
    out[phase+"."+organ+"."+sub]=c->GetSubstanceQuantity(*bg->GetSubstanceManager().GetSubstance(sub))->GetMass(MassUnit::g);
  }
 };
 for(const std::string sub:{"Glucose","AminoAcids","Triacylglycerol","Oxygen","CarbonDioxide","Lactate"})
  out["molar."+sub]=bg->GetSubstanceManager().GetSubstance(sub)->GetMolarMass(MassPerAmountUnit::g_Per_mol);
 observe("before");bool rejected=false;std::string error;
 try{
  ihm_signed::validate(r);ihm_signed::Scope scope(enabled?&r:nullptr);
  dynamic_cast<Energy&>(bg->GetEnergy()).PreProcess();
  out["heat"]=bg->GetCircuits().GetInternalTemperatureCircuit().GetPath("GroundToInternalCore")->GetNextHeatSource(PowerUnit::W);
  auto& tissue=dynamic_cast<Tissue&>(bg->GetTissue());
  (tissue.*TissueAccess::reaction())(.02);
  if(enabled)ihm_signed::finish(r);
 }catch(const std::exception& e){rejected=true;error=e.what();}
 observe("after");
 out["rejected"]=rejected;out["allowed_W"]=r.allowed_decrement_W;
 out["base_kcal"]=r.muscle_base_kcal;out["aa_kcal"]=r.obligatory_aa_kcal;
 out["mandatory_kcal"]=r.mandatory_anaerobic_kcal;out["requested_kcal"]=r.muscle_requested_kcal;
 out["unmet_kcal"]=r.unmet_muscle_kcal;out["heat_count"]=r.heat_count;out["tissue_count"]=r.tissue_count;
 std::cout<<"FIXTURE\t{"<<std::setprecision(17);bool first=true;
 for(const auto& [k,v]:out){if(!first)std::cout<<',';first=false;std::cout<<'"'<<k<<"\":"<<v;}
 std::cout<<"}"<<std::endl;if(rejected)std::cerr<<error<<std::endl;
}
'''

def native(variant,state,output):
    import json,os,math,hashlib
    from build_biogears_substrate_variant import RUNTIME,sha
    directory=RUNTIME/'variants'/variant;manifest=json.loads((directory/'manifest.json').read_text())
    assert sha(directory/'libbiogears.so.8.0.0')==manifest['library_sha256']
    output.mkdir(parents=True,exist_ok=False)
    source=output/'native_fixture.cpp';source.write_text(NATIVE_FIXTURE)
    exe=output/'native_fixture';build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';donor=BASE/'data/raw/physiology/biogears';sysroot=RUNTIME/'sysroot'
    command=['c++','-std=c++20','-O1',str(source)]
    for path in [directory,BASE/'scripts',donor/'projects/biogears/libBiogears/include',donor/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(path)]
    command += [str(directory/'libbiogears.so.8.0.0'),'-L',str(lib),f'-Wl,-rpath,{directory}',f'-Wl,-rpath,{lib}','-lbiogears_cdm','-o',str(exe)]
    subprocess.run(command,check=True)
    # SONAME resolution must load this precise isolated variant.
    (output/'libbiogears.so.8.0.0').symlink_to(directory/'libbiogears.so.8.0.0')
    env=dict(os.environ,LD_LIBRARY_PATH=os.pathsep.join(map(str,[output,lib,sysroot/'lib'])))
    cwd=output
    runtime=build/'runtime'
    for name in ('BioGearsConfiguration.xml','config','ecg','environments','nutrition','patients','states','substances','xsd','UCEDefs.conf'):
        if (runtime/name).exists():(output/name).symlink_to(runtime/name,target_is_directory=(runtime/name).is_dir())
    results={}
    def run(name,m,w,enabled=True,scarce=False):
        proc=subprocess.run([str(exe),str(state),str(m),str(w),'1' if enabled else '0']+(['scarce'] if scarce else []),cwd=cwd,env=env,text=True,capture_output=True)
        (output/(name+'.log')).write_text(proc.stdout+proc.stderr)
        proc.check_returncode()
        lines=[l.split('\t',1)[1] for l in proc.stdout.splitlines() if l.startswith('FIXTURE\t')]
        assert len(lines)==1;results[name]=json.loads(lines[0]);return results[name]
    disabled=run('disabled',0,0,False);zero=run('zero',0,0)
    assert not disabled['rejected'] and not zero['rejected']
    for key in disabled:
        if key.startswith(('before.','after.')) or key=='heat':assert disabled[key]==zero[key],key
    allowed=zero['allowed_W'];assert allowed>0
    eps=min(.01,allowed*.01)
    positive=run('positive',eps,0);negative=run('negative',-eps,0);eccentric=run('eccentric',eps,-eps)
    for name,row in [('positive',positive),('negative',negative),('eccentric',eccentric)]:
        assert not row['rejected'],name
        assert row['unmet_kcal']==0,name
        for organ in ('Brain','Fat'):
            for key in zero:
                if key.startswith('after.'+organ):assert row[key]==zero[key],(name,key)
        assert row['after.tmr']==zero['after.tmr']
    # This frozen branch supplies the small increment through native TAG oxidation.
    # Require observed branch identity; do not invent a mixed-branch closure.
    for name,row,delta in [('positive',positive,eps),('negative',negative,-eps),('eccentric',eccentric,eps)]:
        for key in ('after.Muscle.Glucose','after.Muscle.AminoAcids','after.glycogen_g','after.Muscle.Lactate'):
            assert row[key]==zero[key],(name,'changed branch',key)
        tag=(zero['after.Muscle.Triacylglycerol']-row['after.Muscle.Triacylglycerol'])/zero['molar.Triacylglycerol']
        oxygen=(zero['after.Muscle.Oxygen']-row['after.Muscle.Oxygen'])/zero['molar.Oxygen']
        co2=(row['after.Muscle.CarbonDioxide']-zero['after.Muscle.CarbonDioxide'])/zero['molar.CarbonDioxide']
        row['increment_chemical_residual_kcal']=tag*7554-delta*.02/4184
        row['increment_oxygen_residual_mol']=oxygen-72.5*tag
        row['increment_carbon_residual_mol']=co2-51*tag
        assert abs(row['increment_chemical_residual_kcal'])<2e-13
        assert abs(row['increment_oxygen_residual_mol'])<2e-14
        assert abs(row['increment_carbon_residual_mol'])<2e-14
        for sub in ('Glucose','AminoAcids','Triacylglycerol','Oxygen'):
            assert row['before.Muscle.'+sub]>=row['after.Muscle.'+sub],(name,sub,'negative extent')
        assert row['before.glycogen_g']>=row['after.glycogen_g']
    assert math.isclose(positive['heat']-zero['heat'],eps,abs_tol=1e-12)
    assert math.isclose(negative['heat']-zero['heat'],-eps,abs_tol=1e-12)
    assert math.isclose(eccentric['heat']-positive['heat'],eps,abs_tol=1e-12)
    for key in positive:
        if key.startswith('after.'):assert positive[key]==eccentric[key],key
    inside=run('inside',-allowed*(1-1e-7),0)
    outside=run('outside',-allowed*(1+1e-7),0)
    assert not inside['rejected'] and outside['rejected']
    assert outside['tissue_count']==0
    for key in outside:
        if key.startswith('before.') and key!='before.tmr':assert outside[key]==outside[key.replace('before.','after.',1)],key
    scarce=run('scarce',eps,0,scarce=True)
    assert not scarce['rejected'] and scarce['unmet_kcal']>0
    for key,value in scarce.items():
        if key.startswith('after.Muscle.') or key=='after.glycogen_g':assert value>=0,key
    report={'variant':variant,'variant_manifest_sha256':sha(directory/'manifest.json'),'library_sha256':sha(directory/'libbiogears.so.8.0.0'),'state':str(state),'state_sha256':sha(state),'fixture_sha256':sha(source),'executable_sha256':sha(exe),'compile_command':command,'results':results,'scope':'Actual one-call Energy PreProcess and Tissue reactions. Does not establish full native energy/carbon closure or whole-engine rollback.'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'allowed_decrement_W':allowed,'report':str(output/'report.json')}))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--native',action='store_true');p.add_argument('--variant',default='whole_body_integrity_signed_muscle_v2');p.add_argument('--state',type=Path,default=BASE/'data/derived/audits/embodied-native-bigr07yz/body/physiology/input-state.xml');p.add_argument('--output',type=Path,default=BASE/'data/derived/audits/signed-muscle-native-local');a=p.parse_args()
    if a.native:native(a.variant,a.state.resolve(),a.output.resolve())
    else:main()
