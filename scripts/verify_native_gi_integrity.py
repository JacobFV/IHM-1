#!/usr/bin/env python3
"""Compile a friend-access native branch probe; retain both passing and failing evidence."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parents[1]
RUNTIME = BASE / 'data/runtime/physiology'
SOURCE = BASE / 'data/raw/physiology/biogears'
AUDIT = BASE / 'data/derived/audits/gi-integrity'
PROBE = r'''
#include <cassert>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/engine/Systems/Gastrointestinal.h>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/cdm/patient/SENutrition.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iomanip>
#include <iostream>
namespace biogears {
// Existing upstream friend grants test-only access without altering production headers.
class BioGearsEngineTest {
public:
 static void probe(Gastrointestinal& g, double input) {
   g.GetStomachContents().GetCalcium().SetValue(input, MassUnit::mg);
   double before = g.m_SmallIntestineChymeCalcium->GetMass(MassUnit::mg);
   double proposed = g.m_CalciumDigestionRate.GetValue(MassPerTimeUnit::mg_Per_s)
       * g.m_dT_s * g.m_data.GetConfiguration().GetCalciumAbsorptionFraction();
   g.DigestNutrient();
   std::cout << "GI," << std::setprecision(17) << input << ','
     << g.GetStomachContents().GetCalcium(MassUnit::mg) << ','
     << g.m_SmallIntestineChymeCalcium->GetMass(MassUnit::mg)-before << ','
     << proposed << '\n';
 }
};
}
int main(int argc, char** argv) {
 auto bg = biogears::CreateBioGearsEngine("probe.log");
 if (!bg->LoadState(argv[1])) return 2;
 auto* g = const_cast<biogears::Gastrointestinal*>(dynamic_cast<const biogears::Gastrointestinal*>(bg->GetGastrointestinalSystem()));
 if (!g) return 3;
 for(double mass : {500., .001, .0001, 1.e-14, 0., -0.175000000020085})
   biogears::BioGearsEngineTest::probe(*g, mass);
}
'''

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(variant, state):
    out = AUDIT / variant
    out.mkdir(parents=True, exist_ok=True)
    cpp = out / 'calcium_probe.cpp'
    cpp.write_text(PROBE)
    binary = out / 'calcium_probe'
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
        if not line.startswith('GI,'): continue
        initial, final, credit, proposed = map(float,line.split(',')[1:])
        expected = min(max(initial,0), max(proposed,0))
        passed = all(map(math.isfinite,[final,credit,proposed])) and math.isclose(credit,expected,abs_tol=1e-12) and math.isclose(initial-final,expected,abs_tol=1e-12) and (final>=0 if initial>=0 else final==initial) and (final == 0 if 0 <= initial <= proposed else True)
        rows.append(dict(initial_mg=initial,final_mg=final,chyme_credit_mg=credit,proposed_mg=proposed,expected_transfer_mg=expected,passed=passed))
    retained = BASE/'data/derived/physiology/native_hour_rest/states/native_final.xml'
    tree = ET.parse(retained)
    negatives=[float(e.attrib['value']) for e in tree.iter() if e.tag.split('}')[-1]=='Calcium' and float(e.attrib.get('value',0))<0]
    report=dict(variant=variant,dependency_sha256=dependency_hashes,library_sha256=sha(selected/'libbiogears.so.8.0.0'),state_path=str(state),state_sha256=sha(state),probe_sha256=sha(cpp),binary_sha256=sha(binary),command=command,ld_library_path=env['LD_LIBRARY_PATH'],native_returncode=result.returncode,rows=rows,retained_original_state_sha256=sha(retained),retained_negative_calcium_mg=negatives,passed=result.returncode==0 and len(rows)==6 and all(r['passed'] for r in rows))
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report['passed']

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--variant',default='whole_body_integrity')
    parser.add_argument('--state',type=Path,default=BASE/'data/derived/physiology/native_baseline_v2/states/native_stabilized.xml')
    args=parser.parse_args()
    raise SystemExit(0 if verify(args.variant,args.state.resolve()) else 1)
