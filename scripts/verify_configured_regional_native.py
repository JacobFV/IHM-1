#!/usr/bin/env python3
"""Opt-in bounded configured native circuit fixture; requires root heavy-slot grant."""
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile
from prepare_configured_regional_skin import ROOT,PARENT,PARENT_MANIFEST_SHA,sha
from ihm.assembly.regional_skin_configuration import validate_configuration
RUNTIME=ROOT/'data/runtime/physiology'
SOURCE=ROOT/'data/raw/physiology/biogears'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared',type=Path,default=ROOT/'data/research/configured_regional_skin/prepared_v1')
    parser.add_argument('--compile-and-run',action='store_true')
    args=parser.parse_args();prepared=args.prepared.resolve()
    receipt=json.loads((prepared/'preparation.json').read_text())
    config=json.loads((prepared/'configuration.json').read_text());validate_configuration(config)
    if sha(PARENT/'manifest.json')!=PARENT_MANIFEST_SHA:raise ValueError('Changed accepted parent')
    for name,key in [('native_configured_regional_skin.h','configured_header_sha256'),
                     ('native_configured_regional_skin_fixture.cpp','fixture_sha256')]:
        if sha(prepared/name)!=receipt[key]:raise ValueError('Changed prepared source')
    if config['configuration_sha256']!=receipt['configuration_sha256']:raise ValueError('Configuration/source mismatch')
    if not args.compile_and_run:
        print(json.dumps({'source_verified':True,'configuration_sha256':config['configuration_sha256'],'native_executed':False}));return
    if sha(PARENT/'libbiogears.so.8.0.0')!=receipt['parent_library_sha256']:raise ValueError('Changed native library')
    out=Path(tempfile.mkdtemp(prefix='configured-regional-native-',dir=ROOT/'data/derived/audits'))
    for name in ['configuration.json','preparation.json','native_configured_regional_skin.h','native_configured_regional_skin_fixture.cpp']:
        shutil.copyfile(prepared/name,out/name)
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';binary=out/'fixture'
    command=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_configured_regional_skin_fixture.cpp')]
    for path in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(path)]
    command+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(PARENT)+':'+str(lib)}
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    def run(cmd,label,timeout):
        with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
            completed=subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources'))]+cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,preexec_fn=limits,timeout=timeout)
        if completed.returncode:raise RuntimeError(label+' failed; retained '+str(out))
    print(str(out),flush=True)
    run(command,'compile',75)
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
        (out/name).symlink_to(build/'runtime'/name)
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True)
    (out/'ldd.txt').write_text(linkage)
    if str(PARENT/'libbiogears.so.8.0.0') not in linkage:raise ValueError('Wrong native runtime lineage')
    run([str(binary)],'fixture',10)
    line=next(x for x in (out/'fixture.stdout').read_text().splitlines() if x.startswith('RESULT '))
    result=json.loads(line[7:])
    if not result['passed']:raise ValueError('Native fixture failed')
    result.update(configuration_sha256=config['configuration_sha256'],native_region_names=[r['name'] for r in config['regions']],
                  parent_variant=PARENT.name,parent_manifest_sha256=PARENT_MANIFEST_SHA,
                  parent_library_sha256=receipt['parent_library_sha256'],compile_command=command,
                  fixture_sha256=sha(binary),production_activation_allowed=False,
                  scope='Actual native circuit solver installation, zero-load and local-pressure/release fixture; no full patient acceptance')
    if sha(PARENT/'manifest.json')!=PARENT_MANIFEST_SHA:raise ValueError('Parent changed during fixture')
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
