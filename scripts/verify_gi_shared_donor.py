#!/usr/bin/env python3
"""One-call held-source diffusion probes; no patient initialization or advancement.

Compiles exact extracted diffusion methods against real native compartment/CDM
classes. Raw/capped fluxes in the report are separately labeled mirrored algebra;
ending mass/concentration comes from source-executed native calls.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/raw/physiology/biogears'
RUNTIME = ROOT / 'data/runtime/physiology'
DIFFUSION = SOURCE / 'projects/biogears/libBiogears/src/engine/Systems/Diffusion.cpp'
PIN = '9f90e429a87ae1ba3b7071fc14525905e768d23a27a0f21d15a2581a49b7d48b'
GLUCOSE = SOURCE / 'share/data/substances/Glucose.xml'
CASES = [
    ('zero', [1000, 1000, 1000], .02, 1, False),
    ('small_positive', [1100, 1000, 900], .02, 1, False),
    ('small_negative', [900, 1000, 1100], .02, 1, False),
    ('e_high_ample', [0, 1000, 0], .02, 1, False),
    ('shared_scarce_leaf', [0, 1000, 0], 1, 1, False),
    ('shared_scarce_parent', [0, 1000, 0], 1, 1, True),
    ('adversarial_exact_minus_km', [0, 800000, 800000], .02, 1, False),
    ('adversarial_below_minus_km', [1000, 900000, 900000], .02, 1, False),
]
HEADER = r'''
#include <biogears/engine/Systems/Diffusion.h>
#include <biogears/engine/Controller/BioGears.h>
#include <biogears/cdm/compartment/substances/SELiquidSubstanceQuantity.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <iomanip>
#include <iostream>
namespace biogears {
'''
HARNESS = r'''
class BioGearsEngineTest {
public:
 static int run(int argc, char** argv) {
   BioGears bg("diffusion.log");
   if (!bg.GetSubstances().LoadSubstanceDirectory()) return 2;
   auto* sub = bg.GetSubstances().GetSubstance("Glucose");
   if (!sub) return 3;
   auto& mgr=bg.GetCompartments();
   static_cast<SECompartmentManager&>(mgr).AddLiquidCompartmentSubstance(*sub);
   auto& v=mgr.CreateLiquidCompartment("ProbeV");
   auto& e=mgr.CreateLiquidCompartment("ProbeE");
   auto& i=mgr.CreateLiquidCompartment("ProbeI");
   auto set=[&](SELiquidCompartment& c, double mass, double volume) {
     c.GetVolume().SetValue(volume, VolumeUnit::mL);
     auto& q=*c.GetSubstanceQuantity(*sub);
     q.GetMass().SetValue(mass, MassUnit::ug);
     q.Balance(BalanceLiquidBy::Mass);
   };
   set(v,std::stod(argv[1]),1);set(i,std::stod(argv[3]),1);
   bool parent=std::stoi(argv[6]);
   if (parent) {
     auto& a=mgr.CreateLiquidCompartment("ProbeEa");
     auto& b=mgr.CreateLiquidCompartment("ProbeEb");
     set(a,std::stod(argv[2])*.5,.5);set(b,std::stod(argv[2])*.5,.5);
     e.AddChild(a);e.AddChild(b);e.StateChange();
   } else set(e,std::stod(argv[2]),1);
   DiffusionCalculator d(bg);d.m_dt_s=std::stod(argv[4]);
   DiffusionCalculator::DiffusionCompartmentSet s={nullptr,&v,&e,&i};
   d.CalculateFacilitatedDiffusion(s,*sub,std::stod(argv[5]));
   std::cout<<"RESULT"<<std::setprecision(17);
   for(auto* c:{&v,&e,&i}) std::cout<<','<<c->GetSubstanceQuantity(*sub)->GetMass(MassUnit::ug);
   for(auto* c:{&v,&e,&i}) std::cout<<','<<c->GetSubstanceQuantity(*sub)->GetConcentration(MassPerVolumeUnit::g_Per_mL);
   std::cout<<','<<sub->GetMichaelisCoefficient()<<','<<sub->GetMaximumDiffusionFlux(MassPerAreaTimeUnit::g_Per_cm2_s)<<'\n';
   return 0;
 }
};
}
int main(int argc,char** argv) {return biogears::BioGearsEngineTest::run(argc,argv);}
'''


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))


def extract(text, name):
    start = text.index('void DiffusionCalculator::' + name + '(')
    end = text.index('\n//-------------------------------------------------------------------------------', start)
    return text[start:end]


def mirror(masses, dt, coefficient, km):
    concentration = [m / 1e6 for m in masses]  # every total volume is 1 mL
    gradients = [concentration[0]-concentration[1], concentration[1]-concentration[2]]
    raw = [math.copysign(math.inf, d) if km+d == 0 else 1e6*coefficient*d/(km+d)*dt for d in gradients]
    ve = min(raw[0], masses[0]) if raw[0] > 0 else -min(-raw[0], masses[1])
    ei = min(raw[1], masses[1]) if raw[1] > 0 else -min(-raw[1], masses[2])
    return {'gradient_g_per_mL': gradients, 'raw_ug': [v if math.isfinite(v) else str(v) for v in raw],
            'independently_capped_ug': [ve, ei], 'incidence_delta_ug': [-ve, ve-ei, ei],
            'shared_E_requested_ug': max(-ve, 0)+max(ei, 0)}


def main():
    assert sha(DIFFUSION) == PIN, 'Held diffusion source changed; review before repinning'
    out = Path(tempfile.mkdtemp(prefix='gi-shared-donor-', dir=ROOT/'data/derived/audits'))
    methods = '\n'.join(extract(DIFFUSION.read_text(), name) for name in [
        'CalculateFacilitatedDiffusion', 'DistributeMassbyMassWeighted', 'DistributeMassbyVolumeWeighted'])
    cpp = out/'probe.cpp'; cpp.write_text(HEADER+methods+HARNESS)
    (out/'extracted_methods.cpp').write_text(methods)
    binary = out/'probe'
    build = RUNTIME/'biogears-build'; lib = build/'outputs/Release/lib'
    selected = RUNTIME/'variants/whole_body_integrity_substrate_availability'
    cmd = ['c++', '-std=c++20', '-O0', str(cpp)]
    for p in [SOURCE/'projects/biogears/libBiogears/include', SOURCE/'projects/biogears-common/include',
              RUNTIME/'sysroot/usr/include', RUNTIME/'sysroot/usr/include/eigen3', build/'projects/biogears/generated/Release']:
        cmd += ['-I', str(p)]
    cmd += ['-L', str(lib), '-lbiogears', '-lbiogears_cdm', '-o', str(binary)]
    env = {**os.environ, 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1'}
    def run(command, name, cwd=out):
        with (out/(name+'.stdout')).open('w') as stdout, (out/(name+'.stderr')).open('w') as stderr:
            result = subprocess.run(['/usr/bin/time', '-v', '-o', str(out/(name+'.resources'))]+command,
                                    cwd=cwd, env=env, stdout=stdout, stderr=stderr,
                                    timeout=150, preexec_fn=limits)
        if result.returncode:
            raise RuntimeError(f'{name} failed ({result.returncode}); retained evidence: {out}')
        return (out/(name+'.stdout')).read_text()
    print(out, flush=True)
    run(cmd, 'compile')
    multi = subprocess.check_output(['c++','-print-multiarch'],text=True).strip()
    env['LD_LIBRARY_PATH']=f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'
    linkage = subprocess.check_output(['ldd', str(binary)], env=env, text=True)
    (out/'ldd.txt').write_text(linkage)
    dependencies = [Path(line.split('=>')[1].split(' (')[0].strip()) for line in linkage.splitlines() if '=>' in line]
    assert selected/'libbiogears.so.8.0.0' in dependencies
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
        (out/name).symlink_to(build/'runtime'/name)
    assert sha(out/'substances/Glucose.xml') == sha(GLUCOSE)
    km = float(next(el for el in ET.parse(GLUCOSE).iter() if el.tag.endswith('}MichaelisCoefficient')).attrib['value'])
    rows = {}; checks = {}
    for name, masses, dt, coefficient, parent in CASES:
        stdout = run([str(binary), *map(str,masses),str(dt),str(coefficient),str(int(parent))],name)
        fields = list(map(float,next(line for line in stdout.splitlines() if line.startswith('RESULT,')).split(',')[1:]))
        final, concentrations = fields[:3],fields[3:6]
        mirrored = mirror(masses,dt,coefficient,km)
        predicted = [m+d for m,d in zip(masses,mirrored['incidence_delta_ug'])]
        if parent: predicted[1]=max(predicted[1],0)
        checks[name+'/native_matches_incidence_with_source_parent_debit_cap'] = all(math.isclose(a,b,abs_tol=1e-8,rel_tol=1e-12) for a,b in zip(final,predicted))
        checks[name+'/held_definition'] = fields[6:] == [km,1]
        rows[name]={'initial_mass_ug':masses,'volume_mL':[1,1,1], 'dt_s':dt,'combined_coefficient_g_per_s':coefficient,
                    'extracellular_has_children':parent, 'native_final_mass_ug':final,'native_final_concentration_g_per_mL':concentrations,
                    'native_mass_residual_ug':sum(final)-sum(masses), 'native_nonnegative':all(v>=0 for v in final),
                    'mirrored_algebra_not_native_counters':mirrored,
                    'adversarial_outside_retained_physiological_range':name.startswith('adversarial')}
    checks['shared_leaf_negative_reproduced'] = rows['shared_scarce_leaf']['native_final_mass_ug'][1] == -1000
    checks['shared_parent_mass_creation_reproduced'] = rows['shared_scarce_parent']['native_mass_residual_ug'] == 1000
    checks['leaf_negative_mass_zero_concentration'] = rows['shared_scarce_leaf']['native_final_concentration_g_per_mL'][1] == 0
    checks['below_minus_km_against_gradient'] = rows['adversarial_below_minus_km']['native_final_mass_ug'][0] < 1000
    checks['zero_unchanged'] = rows['zero']['native_final_mass_ug'] == [1000,1000,1000]
    report={'passed':all(checks.values()),'checks':checks,'cases':rows,
            'scope':'Exact extracted held native diffusion methods, real native compartments and CDM balance. No whole-body state loaded or advanced. Diagnostic success includes retained defects, not physiology acceptance.',
            'source_sha256':{str(DIFFUSION):sha(DIFFUSION),str(GLUCOSE):sha(GLUCOSE)},
            'extracted_methods_sha256':sha(out/'extracted_methods.cpp'),'probe_source_sha256':sha(cpp),'binary_sha256':sha(binary),
            'dependency_sha256':{str(p):sha(p) for p in dependencies},'compile_command':cmd,
            'limits':{'address_space_bytes':1024**3,'cpu_seconds_per_child':120,'threads':1,'wall_seconds_per_child':150},
            'resources':{p.name:p.read_text() for p in out.glob('*.resources')}}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':report['passed'],'checks':checks},indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
