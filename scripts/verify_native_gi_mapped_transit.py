#!/usr/bin/env python3
"""One bounded compile and native-CDM transit fixture; no patient advance."""
import json,os,subprocess,tempfile
from pathlib import Path
from build_biogears_shared_donor_variant import ROOT,RUNTIME,SOURCE,sha
from verify_gi_shared_donor import limits

def main():
    out=Path(tempfile.mkdtemp(prefix='native-gi-mapped-',dir=ROOT/'data/derived/audits'));print(out,flush=True)
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';variant=RUNTIME/'variants/whole_body_integrity_gi_absorption'
    inputs=[ROOT/'scripts/native_gi_mapped_transit.h',ROOT/'scripts/native_gi_mapped_transit_probe.cpp',variant/'libbiogears.so.8.0.0']
    hashes={str(p):sha(p) for p in inputs}
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    cmd=['c++','-std=c++20','-O0',str(inputs[1])]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:cmd+=['-I',str(p)]
    cmd+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(out/'probe')]
    def run(args,name):
        with (out/(name+'.stdout')).open('w') as stdout,(out/(name+'.stderr')).open('w') as stderr:
            r=subprocess.run(['/usr/bin/time','-v','-o',str(out/(name+'.resources'))]+args,cwd=out,env=env,stdout=stdout,stderr=stderr,timeout=90,preexec_fn=limits)
        if r.returncode:raise RuntimeError(f'{name} failed; retained {out}')
    run(cmd,'compile')
    multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip();env['LD_LIBRARY_PATH']=f'{variant}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'
    linkage=subprocess.check_output(['ldd',str(out/'probe')],env=env,text=True);(out/'ldd.txt').write_text(linkage)
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(out/name).symlink_to(build/'runtime'/name)
    run([str(out/'probe')],'native')
    if hashes!={str(p):sha(p) for p in inputs}:raise RuntimeError('Inputs changed during fixture')
    stdout=(out/'native.stdout').read_text();assert 'PASS mapped_species=' in stdout
    report=dict(passed=True,scope='Actual mapped native CDM stores and native calculator PostProcess; synthetic post-absorption state. No patient load, circuit Process solve or physiology advance.',inputs_sha256=hashes,binary_sha256=sha(out/'probe'),compile_command=cmd,stdout=stdout,resources={p.name:p.read_text() for p in out.glob('*.resources')})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(stdout)
if __name__=='__main__':main()
