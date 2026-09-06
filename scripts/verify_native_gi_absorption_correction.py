#!/usr/bin/env python3
"""Same actual linked GI fixture and independent acceptance gate for both libraries."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from verify_native_gi_absorption_boundary import PROBE as ORIGINAL_PROBE
from verify_gi_shared_donor import limits
from build_biogears_shared_donor_variant import ROOT,RUNTIME,SOURCE,sha

PARENT=RUNTIME/'variants/whole_body_integrity_signed_muscle_v2'
FIXED=RUNTIME/'variants/whole_body_integrity_gi_absorption'
CASES=['glucose_tail','glucose_exact','sodium_tail','fat_tail','fat_exact','ample','aa_tail','aa_exact','aa_sodium_tail','sequential_sodium','zero_all','aa_ample','fat_ample','multi_ample']
PROBE=ORIGINAL_PROBE.replace('mode=="sodium_tail"?1.e-7:1.', '(mode=="sodium_tail"||mode=="aa_sodium_tail")?1.e-7:(mode=="sequential_sodium"?2.5e-5:(mode=="zero_all"?0.:1.))')
PROBE=PROBE.replace('  for(auto* s:c.GetSubstanceQuantities())', '''  if(mode=="aa_tail"||mode=="aa_exact") {
   q(c,"Glucose")->GetMass().SetValue(0,MassUnit::g);
   q(c,"AminoAcids")->GetMass().SetValue(mode=="aa_tail"?1.e-6:exactGlucose/2.,MassUnit::g);
   q(c,"Triacylglycerol")->GetMass().SetValue(1,MassUnit::g);
  }
  if(mode=="aa_sodium_tail"||mode=="sequential_sodium"||mode=="aa_ample"||mode=="multi_ample") {
   q(c,"Triacylglycerol")->GetMass().SetValue(0,MassUnit::g);
   q(c,"AminoAcids")->GetMass().SetValue(1,MassUnit::g);
   q(c,"Glucose")->GetMass().SetValue((mode=="sequential_sodium"||mode=="multi_ample")?1.:0.,MassUnit::g);
  }
  if(mode=="fat_ample")q(c,"Triacylglycerol")->GetMass().SetValue(1.,MassUnit::g);
  if(mode=="zero_all")q(c,"Triacylglycerol")->GetMass().SetValue(0,MassUnit::g);
  for(auto* s:c.GetSubstanceQuantities())''',1)


def expected_credits(initial):
    glucose,na,aa,tag,calcium,chloride=initial
    requested=(9.*(max(glucose,aa,tag)/(2.+max(glucose,aa,tag))))/3600.*.02
    glucose_extent=min(requested,na,glucose/2.)
    aa_extent=min(requested,na-glucose_extent,aa)
    return [2.*glucose_extent,glucose_extent+aa_extent,aa_extent,min(tag,1.9*.02/1000.),0.,0.]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant',choices=['original','corrected'],default='corrected')
    parser.add_argument('--probe',type=Path)
    parser.add_argument('--original-report',type=Path)
    args=parser.parse_args();selected=PARENT if args.variant=='original' else FIXED
    out=Path(tempfile.mkdtemp(prefix='gi-absorption-correction-',dir=ROOT/'data/derived/audits'));print(out,flush=True)
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';cpp=out/'probe.cpp';cpp.write_text(PROBE)
    binary=args.probe.resolve() if args.probe else out/'probe'
    cmd=['c++','-std=c++20','-O0',str(cpp)]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:cmd+=['-I',str(p)]
    cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    def run(command,name):
        with (out/(name+'.stdout')).open('w') as stdout,(out/(name+'.stderr')).open('w') as stderr:
            r=subprocess.run(['/usr/bin/time','-v','-o',str(out/(name+'.resources'))]+command,cwd=out,env=env,stdout=stdout,stderr=stderr,timeout=150,preexec_fn=limits)
        if r.returncode:raise RuntimeError(f'{name} failed; retained at {out}')
        return (out/(name+'.stdout')).read_text()
    if args.probe:
        if not args.original_report:raise ValueError('Reused probe requires original report identity')
        original=json.loads(args.original_report.read_text())
        if sha(binary)!=original['binary_sha256'] or sha(cpp)!=original['probe_source_sha256']:raise ValueError('Reused probe mismatch')
    else:run(cmd,'compile')
    manifest=json.loads((selected/'manifest.json').read_text())
    if sha(selected/'libbiogears.so.8.0.0')!=manifest['library_sha256']:raise ValueError('Selected library mismatch')
    multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip();env['LD_LIBRARY_PATH']=f'{selected}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
    dependencies=[Path(line.split('=>')[1].split(' (')[0].strip()) for line in linkage.splitlines() if '=>' in line]
    if selected/'libbiogears.so.8.0.0' not in dependencies:raise ValueError('Wrong loaded library')
    for n in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(out/n).symlink_to(build/'runtime'/n)
    rows={};checks={};species=['glucose','sodium','amino_acids','TAG','calcium','chloride']
    for case in CASES:
        text=run([str(binary),case],case)
        initial=list(map(float,next(l for l in text.splitlines() if l.startswith('INITIAL,')).split(',')[1:]))
        final=list(map(float,next(l for l in text.splitlines() if l.startswith('FINAL,')).split(',')[1:]));gut,credit=final[:6],final[6:]
        expected=expected_credits(initial)
        # fat_tail preserves legacy independent sodium/water uptake; it is not
        # part of this paired nutrient correction and is checked by parity.
        indices=[0,2,3,4,5] if case=='fat_tail' else range(6)
        checks[case+'/expected_bounded_transfer']=all(math.isclose(credit[i],expected[i],abs_tol=1e-12,rel_tol=1e-12) for i in indices)
        checks[case+'/nonnegative_finite']=all(math.isfinite(v) and v>=0 for v in final)
        checks[case+'/paired_conservation']=all(abs(a-b-c)<=1e-12 for a,b,c in zip(initial,gut,credit))
        rows[case]={'initial_chyme_g':initial,'final_chyme_g':gut,'native_vascular_credit_g':credit,'expected_credit_g':expected,'species_order':species}
    if args.original_report:
        for case in ['ample','zero_all','aa_ample','fat_ample','multi_ample']:
            checks[case+'/exact_native_parity']=rows[case]==original['cases'][case]
    report={'passed':all(checks.values()),'variant':args.variant,'cases':rows,'checks':checks,'binary_path':str(binary),
            'probe_source_sha256':sha(cpp),'binary_sha256':sha(binary),'library_sha256':sha(selected/'libbiogears.so.8.0.0'),
            'manifest_sha256':sha(selected/'manifest.json'),'dependency_sha256':{str(p):sha(p) for p in dependencies},
            'compile_command':None if args.probe else cmd,'reused_original_report':str(args.original_report) if args.original_report else None,
            'resources':{p.name:p.read_text() for p in out.glob('*.resources')},
            'scope':'Actual linked GI absorption; inactive engine skips active fluid circuit update. All required scalar/mass branches native. Original expected to fail bounded transfer gates.'}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({'passed':report['passed'],'checks':len(checks),'failures':[k for k,v in checks.items() if not v]},indent=2))
    return 0 if report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
