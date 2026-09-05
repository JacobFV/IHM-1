#!/usr/bin/env python3
"""Probe real native aqueous sodium transfers, preserving negative regressions."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

BASE=Path(__file__).resolve().parents[1]
RUNTIME=BASE/'data/runtime/physiology'
SOURCE=BASE/'data/raw/physiology/biogears'
AUDIT=BASE/'data/derived/audits/dry-gi-integrity'
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
#include <limits>
#include <exception>
namespace biogears {
class BioGearsEngineTest {
public:
 static void probe(Gastrointestinal& g,const char* label,double sodium,double water,bool missingNa=false,bool missingWater=false,double rate=10.) {
   auto& stomach=g.GetStomachContents();
   // A nonempty carbohydrate pool detects mutation before unknown-water rejection.
   stomach.GetCarbohydrate().SetValue(missingWater?10.:0.,MassUnit::g);
   stomach.GetProtein().SetValue(0.,MassUnit::g);
   stomach.GetFat().SetValue(0.,MassUnit::g);
   stomach.GetCalcium().SetValue(0.,MassUnit::mg);
   stomach.GetSodium().SetValue(sodium,MassUnit::g);
   stomach.GetWater().SetValue(water,VolumeUnit::mL);
   if(missingNa)stomach.GetSodium().Invalidate();
   if(missingWater)stomach.GetWater().Invalidate();
   g.m_DecrementNutrients=true;
   g.m_WaterDigestionRate.SetValue(rate,VolumePerTimeUnit::mL_Per_s);
   const double beforeNa=g.m_SmallIntestineChymeSodium->GetMass(MassUnit::g);
   const double beforeWater=g.m_SmallIntestineChyme->GetVolume(VolumeUnit::mL);
   bool failed=false;
   try {g.DigestNutrient();}
   catch(const std::exception& e){failed=true;std::cout<<"ERROR,"<<label<<','<<e.what()<<'\n';}
   catch(...){failed=true;std::cout<<"ERROR,"<<label<<",nonstandard native exception\n";}
   const double nan=std::numeric_limits<double>::quiet_NaN();
   const bool validNa=stomach.HasSodium(),validWater=stomach.HasWater();
   const bool chymeNa=g.m_SmallIntestineChymeSodium->GetMass().IsValid();
   const bool chymeWater=g.m_SmallIntestineChyme->GetVolume().IsValid();
   std::cout<<"DRYGI,"<<label<<','<<std::setprecision(17)<<sodium<<','<<water<<','
     <<validNa<<','<<(validNa?stomach.GetSodium(MassUnit::g):nan)<<','
     <<validWater<<','<<(validWater?stomach.GetWater(VolumeUnit::mL):nan)<<','
     <<chymeNa<<','<<(chymeNa?g.m_SmallIntestineChymeSodium->GetMass(MassUnit::g)-beforeNa:nan)<<','
     <<chymeWater<<','<<(chymeWater?g.m_SmallIntestineChyme->GetVolume(VolumeUnit::mL)-beforeWater:nan)<<','
     <<failed<<','<<g.m_dT_s<<','<<rate<<','<<stomach.GetCarbohydrate(MassUnit::g)<<'\n';
 }
};
}
int main(int argc,char**argv){
 if(argc!=2)return 2;
 struct Case {const char* label;double sodium;double water;bool missingNa;bool missingWater;double rate=10.;};
 const Case cases[]={
   {"dry_empty",0,0,false,false},{"dry_sodium",1,0,false,false},
   {"depletion",.001,.1,false,false},{"hydrated",1,500,false,false},
   {"hydrated_empty",0,500,false,false},{"unknown_sodium",0,500,true,false},
   {"unknown_water",1,0,false,true},
   {"zero_flow_dry",1,0,false,false,0},
   {"zero_flow_hydrated",1,500,false,false,0},
   {"zero_flow_tiny_water",1,1e-310,false,false,0},
   {"tiny_water_fractional_flow",1,1e-310,false,false,2.5e-309},
 };
 for(const auto& c:cases){
   auto bg=biogears::CreateBioGearsEngine(std::string("probe-")+c.label+".log");
   if(!bg->LoadState(argv[1]))return 3;
   auto* g=const_cast<biogears::Gastrointestinal*>(dynamic_cast<const biogears::Gastrointestinal*>(bg->GetGastrointestinalSystem()));
   if(!g)return 4;
   biogears::BioGearsEngineTest::probe(*g,c.label,c.sodium,c.water,c.missingNa,c.missingWater,c.rate);
 }
}
'''

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify(variant,state):
    out=AUDIT/variant;out.mkdir(parents=True,exist_ok=True)
    cpp=out/'dry_gi_probe.cpp';cpp.write_text(PROBE)
    binary=out/'dry_gi_probe';build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
    command=['c++','-std=c++20','-O2',str(cpp)]
    for path in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(path)]
    command+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
    subprocess.run(command,check=True)
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
        target=out/name
        if not target.exists():target.symlink_to(build/'runtime'/name)
    selected=lib if variant=='upstream' else RUNTIME/'variants'/variant
    multiarch=subprocess.check_output(['c++','-print-multiarch'],text=True).strip()
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','LD_LIBRARY_PATH':f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multiarch}'}
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True)
    (out/'ldd.txt').write_text(linkage)
    resolved=[Path(line.split('=>',1)[1].split(' (',1)[0].strip()) for line in linkage.splitlines() if '=>' in line and 'not found' not in line]
    if selected/'libbiogears.so.8.0.0' not in resolved:raise RuntimeError('Wrong native library resolved')
    result=subprocess.run([str(binary),str(state)],cwd=out,env=env,text=True,capture_output=True,timeout=180)
    (out/'stdout.log').write_text(result.stdout);(out/'stderr.log').write_text(result.stderr)
    rows=[]
    names=['input_sodium_g','input_water_ml','stomach_sodium_valid','stomach_sodium_g','stomach_water_valid','stomach_water_ml','chyme_sodium_valid','chyme_sodium_credit_g','chyme_water_valid','chyme_water_credit_ml','native_exception','dt_s','water_rate_ml_per_s','stomach_carbohydrate_g']
    for line in result.stdout.splitlines():
        if not line.startswith('DRYGI,'):continue
        cells=line.split(',');row=dict(zip(names,map(float,cells[2:])));row['case']=cells[1]
        for name,value in list(row.items()):
            if isinstance(value,float) and not math.isfinite(value):row[name]=None
        expected_water=0 if row['case']=='unknown_water' else min(row['input_water_ml'],row['water_rate_ml_per_s']*row['dt_s'])
        expected_na=(row['input_sodium_g']*(expected_water/row['input_water_ml'])) if row['input_water_ml']>0 and row['case']!='unknown_sodium' else 0.
        if row['case']=='unknown_water':
            passed=bool(row['native_exception']) and bool(row['chyme_sodium_valid']) and row['stomach_sodium_g']==1 and row['chyme_sodium_credit_g']==0 and row['chyme_water_credit_ml']==0 and not row['stomach_water_valid'] and row['stomach_carbohydrate_g']==10
        else:
            passed=not row['native_exception'] and bool(row['chyme_sodium_valid']) and bool(row['chyme_water_valid'])
            passed=passed and math.isclose(row['chyme_sodium_credit_g'] or 0.,expected_na,abs_tol=1e-12) and math.isclose(row['chyme_water_credit_ml'] or 0.,expected_water,abs_tol=1e-10)
            passed=passed and bool(row['stomach_water_valid']) and math.isclose(row['stomach_water_ml'],row['input_water_ml']-expected_water,rel_tol=1e-12,abs_tol=0.)
            if row['case']=='unknown_sodium':passed=passed and not row['stomach_sodium_valid']
            else:passed=passed and bool(row['stomach_sodium_valid']) and math.isclose(row['stomach_sodium_g'] or 0.,row['input_sodium_g']-expected_na,abs_tol=1e-12)
        row.update(expected_sodium_transfer_g=expected_na,expected_water_transfer_ml=expected_water,passed=bool(passed));rows.append(row)
    report={'variant':variant,'native_returncode':result.returncode,'state_path':str(state),'state_sha256':sha(state),
            'probe_sha256':sha(cpp),'binary_sha256':sha(binary),'library_sha256':sha(selected/'libbiogears.so.8.0.0'),
            'dependency_sha256':{str(p):sha(p) for p in resolved if p.is_file()},'compile_command':command,
            'rows':rows,'passed':result.returncode==0 and len(rows)==11 and all(row['passed'] for row in rows)}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['dependency_sha256','compile_command']},indent=2))
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--variant',default='whole_body_integrity_gi_water')
    parser.add_argument('--state',type=Path,default=BASE/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml')
    args=parser.parse_args();raise SystemExit(0 if verify(args.variant,args.state.resolve())['passed'] else 1)
