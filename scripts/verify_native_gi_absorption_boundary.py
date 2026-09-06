#!/usr/bin/env python3
"""Native one-call GI remainder gates plus a bounded retained meal incidence audit."""
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from build_biogears_shared_donor_variant import ROOT, RUNTIME, SOURCE, VARIANT, sha
from verify_gi_shared_donor import limits

GI=RUNTIME/'variants/whole_body_integrity_depletion/Gastrointestinal.cpp'
RETAINED=ROOT/'data/derived/audits/meal-electrolyte-native-02rzs870'
PROBE=r'''
#include <cassert>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/engine/Systems/Gastrointestinal.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iomanip>
#include <iostream>
namespace biogears {
class BioGearsEngineTest {
public:
 static int run(const std::string& mode) {
  BioGears bg("probe.log");
  if(!bg.GetSubstances().LoadSubstanceDirectory()) return 2;
  if(bg.GetState()==EngineState::Active) return 3;
  auto& mgr=bg.GetCompartments();
  const std::vector<std::string> names={"Glucose","Sodium","AminoAcids","Triacylglycerol","Calcium","Chloride"};
  for(const auto& n:names) static_cast<SECompartmentManager&>(mgr).AddLiquidCompartmentSubstance(*bg.GetSubstances().GetSubstance(n));
  auto& c=mgr.CreateLiquidCompartment("ProbeChyme");auto& v=mgr.CreateLiquidCompartment("ProbeVascular");
  c.GetVolume().SetValue(10,VolumeUnit::mL);v.GetVolume().SetValue(10,VolumeUnit::mL);
  for(auto* owner:{&c,&v})for(const auto& n:names){
   auto* q=owner->GetSubstanceQuantity(*bg.GetSubstances().GetSubstance(n));q->GetMass().SetValue(0,MassUnit::g);q->Balance(BalanceLiquidBy::Mass);
  }
  auto q=[&](SELiquidCompartment& owner,const char* name){return owner.GetSubstanceQuantity(*bg.GetSubstances().GetSubstance(name));};
  q(c,"Sodium")->GetMass().SetValue(mode=="sodium_tail"?1.e-7:1.,MassUnit::g);
  double exactGlucose=2.*(9.*(1./(2.+1.)))*(1./3600.)*.02;
  if(mode=="glucose_tail"||mode=="glucose_exact"){
   q(c,"Glucose")->GetMass().SetValue(mode=="glucose_tail"?1.e-6:exactGlucose,MassUnit::g);
   q(c,"Triacylglycerol")->GetMass().SetValue(1.,MassUnit::g);
  } else if(mode=="sodium_tail"||mode=="ample") q(c,"Glucose")->GetMass().SetValue(1.,MassUnit::g);
  else q(c,"Triacylglycerol")->GetMass().SetValue(mode=="fat_tail"?.001:(.019*100.)*.02,MassUnit::mg);
  for(auto* s:c.GetSubstanceQuantities())s->Balance(BalanceLiquidBy::Mass);
  Gastrointestinal g(bg);g.m_dT_s=.02;g.m_SmallIntestineChyme=&c;g.m_vSmallIntestine=&v;
  g.m_SmallIntestineChymeGlucose=q(c,"Glucose");g.m_smallIntestineVascularGlucose=q(v,"Glucose");
  g.m_SmallIntestineChymeSodium=q(c,"Sodium");g.m_SmallIntestineVascularSodium=q(v,"Sodium");
  g.m_SmallIntestineChymeAminoAcids=q(c,"AminoAcids");g.m_smallIntestineVascularAminoAcids=q(v,"AminoAcids");
  g.m_SmallIntestineChymeTriacylglycerol=q(c,"Triacylglycerol");g.m_smallintestineVAscularTriacylglycerol=q(v,"Triacylglycerol");
  g.m_SmallIntestineChymeCalcium=q(c,"Calcium");g.m_SmallIntestineVascularCalcium=q(v,"Calcium");
  std::cout<<"INITIAL"<<std::setprecision(17);for(const auto& n:names)std::cout<<','<<c.GetSubstanceQuantity(*bg.GetSubstances().GetSubstance(n))->GetMass(MassUnit::g);std::cout<<'\n';
  g.AbsorbNutrients();
  std::cout<<"FINAL";for(auto* owner:{&c,&v})for(const auto& n:names)std::cout<<','<<owner->GetSubstanceQuantity(*bg.GetSubstances().GetSubstance(n))->GetMass(MassUnit::g);std::cout<<'\n';
  return 0;
 }
};
}
int main(int argc,char** argv){return argc==2?biogears::BioGearsEngineTest::run(argv[1]):4;}
'''


def retained_audit():
    result={}
    for mode in ['rest','nutrients_only']:
        path=RETAINED/mode/'electrolytes.csv'
        rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(path.open())]
        initial=rows[0];samples=[]
        for row in rows:
            addition=60. if mode=='nutrients_only' and row['time_s']>0 else 0.
            glucose=initial['stomach_carbohydrate_g']+initial['SmallIntestineChyme.Glucose.mass_g']+addition-row['stomach_carbohydrate_g']-row['SmallIntestineChyme.Glucose.mass_g']
            sodium=initial['stomach_sodium_g']+initial['SmallIntestineChyme.Sodium.mass_g']-row['stomach_sodium_g']-row['SmallIntestineChyme.Sodium.mass_g']
            samples.append({'time_s':row['time_s'],'source_mass_equivalent_glucose_net_export_g':glucose,'sodium_net_export_g':sodium,
                            'chyme_chloride_g':row['SmallIntestineChyme.Chloride.mass_g'],'native_sid_mmol_l':row['native_sid_mmol_l'],
                            'arterial_ph':row['Aorta.ph'],'stomach_carbohydrate_g':row['stomach_carbohydrate_g'],'chyme_glucose_g':row['SmallIntestineChyme.Glucose.mass_g']})
        result[mode]={'csv_sha256':sha(path),'samples':samples}
    return result


def main():
    out=Path(tempfile.mkdtemp(prefix='gi-absorption-boundary-',dir=ROOT/'data/derived/audits'))
    print(out,flush=True)
    manifest=json.loads((VARIANT/'manifest.json').read_text());parent=json.loads((RUNTIME/'variants/whole_body_integrity_substrate_availability/manifest.json').read_text())
    gi_obj=str(GI)+'.o'
    if manifest['object_sha256'][gi_obj]!=parent['object_sha256'][gi_obj] or sha(Path(gi_obj))!=manifest['object_sha256'][gi_obj]:raise RuntimeError('GI inherited object changed')
    if sha(VARIANT/'libbiogears.so.8.0.0')!=manifest['library_sha256']:raise RuntimeError('Loaded variant changed')
    cpp=out/'probe.cpp';cpp.write_text(PROBE);binary=out/'probe'
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';cmd=['c++','-std=c++20','-O0',str(cpp)]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:cmd+=['-I',str(p)]
    cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    def run(command,name):
        with (out/(name+'.stdout')).open('w') as stdout,(out/(name+'.stderr')).open('w') as stderr:
            r=subprocess.run(['/usr/bin/time','-v','-o',str(out/(name+'.resources'))]+command,cwd=out,env=env,stdout=stdout,stderr=stderr,timeout=150,preexec_fn=limits)
        if r.returncode:raise RuntimeError(f'{name} failed; receipt {out}')
        return (out/(name+'.stdout')).read_text()
    run(cmd,'compile')
    multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip();env['LD_LIBRARY_PATH']=f'{VARIANT}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
    dependencies=[Path(line.split('=>')[1].split(' (')[0].strip()) for line in linkage.splitlines() if '=>' in line]
    if VARIANT/'libbiogears.so.8.0.0' not in dependencies:raise RuntimeError('Wrong loaded library')
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(out/name).symlink_to(build/'runtime'/name)
    rows={};checks={};species=['glucose','sodium','amino_acids','TAG','calcium','chloride']
    for mode in ['glucose_tail','glucose_exact','sodium_tail','fat_tail','fat_exact','ample']:
        text=run([str(binary),mode],mode)
        initial=list(map(float,next(l for l in text.splitlines() if l.startswith('INITIAL,')).split(',')[1:]))
        end=list(map(float,next(l for l in text.splitlines() if l.startswith('FINAL,')).split(',')[1:]));chyme,vascular=end[:6],end[6:]
        checks[mode+'/finite_nonnegative']=all(math.isfinite(x) and x>=0 for x in end)
        checks[mode+'/paired_species_conservation']=all(abs(a-b-c)<1e-12 for a,b,c in zip(initial,chyme,vascular))
        rows[mode]={'initial_chyme_g':dict(zip(species,initial)),'final_chyme_g':dict(zip(species,chyme)),'vascular_credit_g':dict(zip(species,vascular))}
    for mode in ['glucose_tail','glucose_exact','sodium_tail']:checks[mode+'/finite_glucose_transfer_blocked']=rows[mode]['vascular_credit_g']['glucose']==0 and rows[mode]['initial_chyme_g']['glucose']>0
    for mode in ['fat_tail','fat_exact']:checks[mode+'/finite_TAG_transfer_blocked']=rows[mode]['vascular_credit_g']['TAG']==0 and rows[mode]['initial_chyme_g']['TAG']>0
    checks['ample/transfers_glucose']=rows['ample']['vascular_credit_g']['glucose']>0
    audit=retained_audit()
    checks['retained_chloride_absent']=all(s['chyme_chloride_g']==0 for r in audit.values() for s in r['samples'])
    report={'passed':all(checks.values()),'checks':checks,'native_cases':rows,'retained_meal':audit,
            'scope':'Native AbsorbNutrients solid/ion transfer call; engine inactive so active-only fluid circuit update is skipped. No patient state loaded/advanced. Retained CSV export is net lumen incidence, not gross transporter or portal flux.',
            'gi_source_sha256':sha(GI),'inherited_gi_object_sha256':sha(Path(gi_obj)),'library_sha256':sha(VARIANT/'libbiogears.so.8.0.0'),
            'probe_source_sha256':sha(cpp),'binary_sha256':sha(binary),'compile_command':cmd,'dependency_sha256':{str(p):sha(p) for p in dependencies},
            'resources':{p.name:p.read_text() for p in out.glob('*.resources')}}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({'passed':report['passed'],'checks':len(checks),'failures':[k for k,v in checks.items() if not v]},indent=2))
    return 0 if report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
