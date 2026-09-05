#!/usr/bin/env python3
"""Probe real bilateral renal transfer, transport cap, and glucose bookkeeping."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

BASE = Path(__file__).resolve().parents[1]
RUNTIME = BASE / 'data/runtime/physiology'
SOURCE = BASE / 'data/raw/physiology/biogears'
AUDIT = BASE / 'data/derived/audits/renal-integrity'
PROBE = r'''
#include <cassert>
#include <cmath>
#include <limits>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Systems/Renal.h>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/cdm/circuit/fluid/SEFluidCircuitPath.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iomanip>
#include <iostream>
namespace biogears {
class BioGearsEngineTest {
public:
 static void probe(Renal& r, const char* name, double mass, double maximum, double ratio, double flow) {
   auto& sub=*r.m_Glucose;
   sub.GetClearance().GetRenalTransportMaximum().SetReadOnly(false);
   sub.GetClearance().GetRenalReabsorptionRatio().SetReadOnly(false);
   sub.GetClearance().GetRenalTransportMaximum().SetValue(maximum,MassPerTimeUnit::mg_Per_s);
   sub.GetClearance().GetRenalReabsorptionRatio().SetValue(ratio);
   r.m_LeftReabsorptionPermeabilityModificationFactor=2.;
   r.m_RightReabsorptionPermeabilityModificationFactor=2.;
   r.m_SubstanceTransport.leftGlucoseReabsorptionMass_mg=0;
   r.m_SubstanceTransport.rightGlucoseReabsorptionMass_mg=0;
   SELiquidSubstanceQuantity* tubes[]={r.m_LeftTubules->GetSubstanceQuantity(sub),r.m_RightTubules->GetSubstanceQuantity(sub)};
   SELiquidSubstanceQuantity* blood[]={r.m_LeftPeritubular->GetSubstanceQuantity(sub),r.m_RightPeritubular->GetSubstanceQuantity(sub)};
   SEFluidCircuitPath* paths[]={r.m_LeftReabsorptionResistancePath,r.m_RightReabsorptionResistancePath};
   double initial[2],before[2],proposed[2];
   for(int k=0;k<2;++k) {
     initial[k]=mass/(k+1);
     tubes[k]->GetMass().SetValue(initial[k],MassUnit::mg);
     tubes[k]->Balance(BalanceLiquidBy::Mass);
     paths[k]->GetNextFlow().SetReadOnly(false);
     paths[k]->GetNextFlow().SetValue(flow,VolumePerTimeUnit::mL_Per_s);
     before[k]=blood[k]->GetMass(MassUnit::mg);
     proposed[k]=flow<0 ? 0 : std::isinf(ratio) ? initial[k] : std::min(initial[k],tubes[k]->GetConcentration(MassPerVolumeUnit::mg_Per_mL)*flow*r.m_dt*ratio*.5);
   }
   r.CalculateReabsorptionTransport(sub);
   double book[]={r.m_SubstanceTransport.leftGlucoseReabsorptionMass_mg,r.m_SubstanceTransport.rightGlucoseReabsorptionMass_mg};
   for(int k=0;k<2;++k)
     std::cout << "RENAL," << name << ',' << k << ',' << std::setprecision(17) << initial[k] << ',' << tubes[k]->GetMass(MassUnit::mg) << ',' << blood[k]->GetMass(MassUnit::mg)-before[k] << ',' << proposed[k] << ',' << maximum << ',' << r.m_dt << ',' << book[k] << ',' << sub.GetClearance().GetRenalReabsorptionRate(MassPerTimeUnit::mg_Per_s) << '\n';
 }
};
}
int main(int argc,char** argv) {
 auto bg=biogears::CreateBioGearsEngine("probe.log");
 if(!bg->LoadState(argv[1])) return 2;
 auto* r=const_cast<biogears::Renal*>(dynamic_cast<const biogears::Renal*>(bg->GetRenalSystem()));
 if(!r) return 3;
 double inf=std::numeric_limits<double>::infinity();
 biogears::BioGearsEngineTest::probe(*r,"zero_mass",0,1,inf,2);
 biogears::BioGearsEngineTest::probe(*r,"below_cap",.01,1,inf,2);
 biogears::BioGearsEngineTest::probe(*r,"above_cap",10,1,inf,2);
 biogears::BioGearsEngineTest::probe(*r,"zero_cap",10,0,inf,2);
 biogears::BioGearsEngineTest::probe(*r,"unbounded",10,inf,inf,2);
 biogears::BioGearsEngineTest::probe(*r,"finite_ratio",10,.001,.5,2);
 biogears::BioGearsEngineTest::probe(*r,"zero_ratio",10,1,0,2);
 biogears::BioGearsEngineTest::probe(*r,"zero_flow",10,1,.5,0);
 biogears::BioGearsEngineTest::probe(*r,"backflow",10,1,inf,-2);
 biogears::BioGearsEngineTest::probe(*r,"available_bound",.0001,1,1000000,2);
}

'''

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(variant, state):
    out = AUDIT / variant
    out.mkdir(parents=True, exist_ok=True)
    cpp = out / 'renal_probe.cpp'
    cpp.write_text(PROBE)
    binary = out / 'renal_probe'
    build = RUNTIME / 'biogears-build'
    lib = build / 'outputs/Release/lib'
    command = ['c++', '-std=c++20', '-O2', str(cpp)]
    for path in [SOURCE/'projects/biogears/libBiogears/include', SOURCE/'projects/biogears-common/include', RUNTIME/'sysroot/usr/include', RUNTIME/'sysroot/usr/include/eigen3', build/'projects/biogears/generated/Release']:
        command += ['-I', str(path)]
    command += ['-L', str(lib), '-lbiogears', '-lbiogears_cdm', '-o', str(binary)]
    subprocess.run(command, check=True)
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
        target = out/name
        if not target.exists():
            target.symlink_to(build/'runtime'/name)
    selected = lib if variant == 'upstream' else RUNTIME/'variants'/variant
    multiarch = subprocess.check_output(['c++', '-print-multiarch'], text=True).strip()
    env = {**os.environ, 'OPENBLAS_NUM_THREADS':'1', 'LD_LIBRARY_PATH':f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multiarch}'}
    result = subprocess.run([str(binary),str(state)], cwd=out, env=env, capture_output=True, text=True)
    (out/'stdout.log').write_text(result.stdout)
    (out/'stderr.log').write_text(result.stderr)
    linkage = subprocess.check_output(['ldd',str(binary)],env=env,text=True)
    (out/'ldd.txt').write_text(linkage)
    resolved = [Path(line.split('=>',1)[1].split(' (',1)[0].strip()) for line in linkage.splitlines() if '=>' in line and 'not found' not in line]
    if (selected/'libbiogears.so.8.0.0') not in resolved:
        raise RuntimeError('Probe resolved the wrong BioGears library')
    dependency_hashes = {str(path):sha(path) for path in resolved if path.is_file()}
    rows=[]
    for line in result.stdout.splitlines():
        if not line.startswith('RENAL,'): continue
        _,case,kidney,*values=line.split(',')
        initial,final,credit,proposed,maximum,dt,book,rate=map(float,values)
        expected=min(proposed,maximum*dt)
        passed=all(map(math.isfinite,[initial,final,credit,proposed,dt,book,rate])) and final>=0 and 0<=expected<=initial and all(math.isclose(value,expected,abs_tol=1e-10) for value in [initial-final,credit,book])
        rows.append(dict(case=case,kidney=int(kidney),initial_mg=initial,final_mg=final,peritubular_credit_mg=credit,proposed_mg=proposed,maximum_mg_per_s=maximum if math.isfinite(maximum) else 'infinity',dt_s=dt,bookkeeping_mg=book,total_rate_mg_per_s=rate,expected_transfer_mg=expected,passed=passed))
    for row in rows:
        paired=[r for r in rows if r['case']==row['case']]
        row['passed']=row['passed'] and len(paired)==2 and math.isclose(row['total_rate_mg_per_s']*row['dt_s'],sum(r['expected_transfer_mg'] for r in paired),abs_tol=1e-10)
    report=dict(variant=variant,dependency_sha256=dependency_hashes,library_sha256=sha(selected/'libbiogears.so.8.0.0'),state_path=str(state),state_sha256=sha(state),probe_sha256=sha(cpp),binary_sha256=sha(binary),command=command,ld_library_path=env['LD_LIBRARY_PATH'],native_returncode=result.returncode,rows=rows,passed=result.returncode==0 and len(rows)==20 and all(r['passed'] for r in rows))
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report['passed']

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--variant',default='whole_body_integrity_renal')
    parser.add_argument('--state',type=Path,default=BASE/'data/derived/physiology/native_baseline_v2/states/native_stabilized.xml')
    args=parser.parse_args()
    raise SystemExit(0 if verify(args.variant,args.state.resolve()) else 1)
