#!/usr/bin/env python3
"""Original/corrected real native compartment fixtures, extracted or linked methods."""
import argparse
import json
import math
import os
from pathlib import Path
import random
import subprocess
import tempfile

from build_biogears_shared_donor_variant import DONOR, PARENT, VARIANT, corrected_source, sha
from verify_gi_shared_donor import ROOT, SOURCE, RUNTIME, HEADER, HARNESS, CASES, GLUCOSE, extract, limits
from verify_gi_shared_donor_allocation import allocate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['extracted','linked'],default='extracted')
    args=parser.parse_args()
    out=Path(tempfile.mkdtemp(prefix='gi-native-donor-correction-',dir=ROOT/'data/derived/audits'))
    print(out,flush=True)
    before=DONOR.read_text();after=corrected_source(before)
    variants={'original':before,'corrected':after}
    if args.mode=='linked':
        manifest=json.loads((VARIANT/'manifest.json').read_text())
        if sha(VARIANT/'libbiogears.so.8.0.0')!=manifest['library_sha256'] or sha(VARIANT/'Diffusion.cpp')!=manifest['patched_source_sha256']:
            raise RuntimeError('Corrected production variant identity changed')
        if (VARIANT/'Diffusion.cpp').read_text()!=after:raise RuntimeError('Fixture patch differs from production variant')
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
    multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip()
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    def run(command,work,name):
        with (work/(name+'.stdout')).open('w') as stdout,(work/(name+'.stderr')).open('w') as stderr:
            process=subprocess.run(['/usr/bin/time','-v','-o',str(work/(name+'.resources'))]+command,cwd=work,env=env,
                                   stdout=stdout,stderr=stderr,timeout=150,preexec_fn=limits)
        if process.returncode:raise RuntimeError(f'{name} failed, see {work}')
        return (work/(name+'.stdout')).read_text()
    cases=list(CASES)
    cases += [('scarce_unequal_leaf',[100,1000,900],2.,1.,False),
              ('scarce_unequal_parent',[100,1000,900],2.,1.,True),
              ('empty_donor_leaf',[0,0,0],1.,1.,False),
              ('empty_donor_parent',[0,0,0],1.,1.,True)]
    rng=random.Random(20260905)
    for n in range(8):
        mass=10**rng.uniform(-4,4)
        v=mass*rng.uniform(0,.9);i=mass*rng.uniform(0,.9)
        for parent in [False,True]:
            cases.append((f'budget_rounding_{n}_'+('parent' if parent else 'leaf'),[v,mass,i],10.,1.,parent))
    results={};receipts={};checks={};shared_binary=None
    for label,source_text in variants.items():
        work=out/label;work.mkdir()
        methods='\n'.join(extract(source_text,name) for name in ['CalculateFacilitatedDiffusion','DistributeMassbyMassWeighted','DistributeMassbyVolumeWeighted'])
        cpp=work/'probe.cpp';cpp.write_text(HEADER+(methods if args.mode=='extracted' else '')+HARNESS)
        binary=work/'probe'
        command=['c++','-std=c++20','-O0',str(cpp)]
        for path in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
            command+=['-I',str(path)]
        command+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
        if args.mode=='linked' and shared_binary is not None:
            binary=shared_binary
        else:
            run(command,work,'compile')
            shared_binary=binary
        selected=VARIANT if args.mode=='linked' and label=='corrected' else PARENT
        env['LD_LIBRARY_PATH']=f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'
        linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True)
        (work/'ldd.txt').write_text(linkage)
        dependencies=[Path(line.split('=>')[1].split(' (')[0].strip()) for line in linkage.splitlines() if '=>' in line]
        if selected/'libbiogears.so.8.0.0' not in dependencies:raise RuntimeError('Wrong loaded library')
        for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
            (work/name).symlink_to(build/'runtime'/name)
        if sha(work/'substances/Glucose.xml')!=sha(GLUCOSE):raise RuntimeError('Held glucose changed')
        receipts[label]={'source_sha256':sha(cpp),'binary_sha256':sha(binary),'command':command,
                         'dependency_sha256':{str(p):sha(p) for p in dependencies},'loaded_variant':selected.name}
        results[label]={}
        for name,initial,dt,coefficient,parent in cases:
            stdout=run([str(binary),*map(str,initial),str(dt),str(coefficient),str(int(parent))],work,name)
            fields=list(map(float,next(line for line in stdout.splitlines() if line.startswith('RESULT,')).split(',')[1:]))
            final=fields[:3]
            # Independent expected signed incidence, using initial concentrations.
            gradients=[(initial[0]-initial[1])/1e6,(initial[1]-initial[2])/1e6]
            raw=[math.copysign(math.inf,g) if .8+g==0 else 1e6*coefficient*g/(.8+g)*dt for g in gradients]
            capped=[min(raw[0],initial[0]) if raw[0]>0 else -min(-raw[0],initial[1]),
                    min(raw[1],initial[1]) if raw[1]>0 else -min(-raw[1],initial[2])]
            ve,ei=allocate(initial,*capped,True)
            expected=[initial[0]-ve,initial[1]+ve-ei,initial[2]+ei]
            ledger=[initial[0]-final[0],final[2]-initial[2]]
            residue=sum(final)-sum(initial)
            results[label][name]={'initial_mass_ug':initial,'native_final_mass_ug':final,'native_final_concentration_g_per_mL':fields[3:6],
                                  'mass_residual_ug':residue,'inferred_transfers_from_native_endpoints_ug':ledger,
                                  'mirrored_individually_capped_transfers_ug':capped,'expected_conservative_mass_ug':expected}
            if label=='corrected':
                checks[name+'/finite_nonnegative']=all(math.isfinite(x) and x>=0 for x in final)
                checks[name+'/conservation']=abs(residue)<=1e-8
                checks[name+'/independent_allocation']=all(math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-8) for a,b in zip(final,expected))
                checks[name+'/owner_incidence']=math.isclose(final[1]-initial[1],ledger[0]-ledger[1],rel_tol=1e-12,abs_tol=1e-8)
        receipts[label]['resources']={p.name:p.read_text() for p in work.glob('*.resources')}
    for name in ['zero','small_positive','small_negative','e_high_ample','empty_donor_leaf','empty_donor_parent','adversarial_exact_minus_km','adversarial_below_minus_km']:
        checks[name+'/exact_unchanged_native_outputs']=results['original'][name]==results['corrected'][name]
    checks['original_leaf_negative_reproduced']=results['original']['shared_scarce_leaf']['native_final_mass_ug'][1]==-1000
    checks['original_parent_creation_reproduced']=results['original']['shared_scarce_parent']['mass_residual_ug']==1000
    for layout in ['leaf','parent']:
        checks['corrected_scarce_'+layout+'_500_0_500']=results['corrected']['shared_scarce_'+layout]['native_final_mass_ug']==[500,0,500]
    for stem in ['scarce_unequal']+[f'budget_rounding_{n}' for n in range(8)]:
        checks[stem+'/native_leaf_parent_equivalence']=results['corrected'][stem+'_leaf']['native_final_mass_ug']==results['corrected'][stem+'_parent']['native_final_mass_ug']
    report={'passed':all(checks.values()),'mode':args.mode,'checks':checks,'results':results,'receipts':receipts,
            'donor_source_sha256':sha(DONOR),'patch_generator_sha256':sha(ROOT/'scripts/build_biogears_shared_donor_variant.py'),
            'fixture_sha256':sha(Path(__file__)),'scope':'Actual native original/corrected methods and CDM compartments. No body state initialization or advancement. Signed gradient law is unchanged, including its adversarial defect.'}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':report['passed'],'checks':len(checks),'cases_per_variant':len(cases),'failures':[k for k,v in checks.items() if not v]},indent=2))
    return 0 if report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
