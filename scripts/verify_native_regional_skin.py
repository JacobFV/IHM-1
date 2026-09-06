#!/usr/bin/env python3
"""Build/run the tiny native regional Skin circuit fixture; no patient advance.

Explicit --compile-and-run is required; default only checks fixture/source pins.
Schedule compilation with the root native slot before invoking that option.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT/'data/runtime/physiology'
SOURCE = ROOT/'data/raw/physiology/biogears'
PARENT = RUNTIME/'variants/whole_body_integrity_signed_muscle_v2'
LIBRARY_PIN = '5a02e2942e6ac5ff3cdce352246e5073de7871b500a572b4ab2d8f50eb87b695'
MANIFEST_PIN = 'ef7d19eb63d038222d60a7e6c9fe43e41e3a3145ad2cd85f76adcec6976f0ac4'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def limits():
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(120,120))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compile-and-run',action='store_true')
    parser.add_argument('--reuse-built',type=Path)
    args=parser.parse_args()
    assert sha(PARENT/'manifest.json')==MANIFEST_PIN, 'Parent manifest changed'
    assert sha(PARENT/'libbiogears.so.8.0.0')==LIBRARY_PIN, 'Parent library changed'
    inputs=[ROOT/'scripts/native_regional_skin.h',ROOT/'scripts/native_regional_skin_fixture.cpp']
    receipts={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    if not args.compile_and_run and args.reuse_built is None:
        print(json.dumps({'prepared':True,'compiled':False,'source_sha256':receipts}));return
    out=args.reuse_built.resolve() if args.reuse_built else Path(tempfile.mkdtemp(prefix='native-regional-skin-',dir=ROOT/'data/derived/audits'))
    for p in inputs:
        if args.reuse_built:
            assert sha(out/p.name)==sha(p), 'Reused fixture source changed'
        else:(out/p.name).write_bytes(p.read_bytes())
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';binary=out/'fixture'
    command=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_regional_skin_fixture.cpp')]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(p)]
    command+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(PARENT)+':'+str(lib)}
    def run(cmd,label,cwd):
        with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
            subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources'))]+cmd,cwd=cwd,env=env,stdout=stdout,stderr=stderr,preexec_fn=limits,check=True,timeout=150)
    print(out,flush=True)
    if not args.reuse_built:run(command,'compile',out)
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
        target=out/name
        if not target.exists():target.symlink_to(build/'runtime'/name)
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True)
    (out/'ldd.txt').write_text(linkage)
    assert str(PARENT/'libbiogears.so.8.0.0') in linkage, 'Fixture linked wrong variant'
    run([str(binary)],'fixture',out)
    line=next(line for line in (out/'fixture.stdout').read_text().splitlines() if line.startswith('RESULT '))
    report=json.loads(line[len('RESULT '):])
    assert report['passed'] and not report['patient_initialized']
    assert sha(PARENT/'manifest.json')==MANIFEST_PIN and sha(PARENT/'libbiogears.so.8.0.0')==LIBRARY_PIN
    assert receipts=={str(p.relative_to(ROOT)):sha(p) for p in inputs}, 'Fixture inputs changed during execution'
    report.update(source_sha256=receipts,parent_variant=PARENT.name,parent_library_sha256=LIBRARY_PIN,parent_manifest_sha256=MANIFEST_PIN,compile_command=command,fixture_executable_sha256=sha(binary),scope='Native graph/compartment fixture only; no production patient integration or chemical transport validation',resources={name:(out/(name+'.resources')).read_text() for name in ['compile','fixture']})
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
