#!/usr/bin/env python3
"""Source preparation by default; explicit native slot required for --compile-and-run."""
from pathlib import Path
import argparse,hashlib,json,os,resource,subprocess,tempfile
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/raw/physiology/biogears'
RUNTIME=ROOT/'data/runtime/physiology'
VARIANT=RUNTIME/'variants/whole_body_integrity_gi_absorption'
REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'
SOURCE_FILES=['src/engine/Controller/BioGears.cpp','src/engine/Systems/Respiratory.cpp','src/engine/Systems/Nervous.cpp','include/biogears/cdm/circuit/SECircuitCalculator.inl']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
    pins={}
    for name in SOURCE_FILES:
        p=SOURCE/'projects/biogears/libBiogears'/name
        original=subprocess.check_output(['git','-C',str(SOURCE),'show',f'{REVISION}:{p.relative_to(SOURCE)}'])
        assert original==p.read_bytes(),'Held source changed: '+name
        pins[str(p.relative_to(ROOT))]=sha(p)
    for name in ['native_respiratory_work.h','native_respiratory_work_probe.cpp']:
        p=ROOT/'scripts'/name
        assert p.is_file(),'Missing respiratory work implementation: '+name
        pins[str(p.relative_to(ROOT))]=sha(p)
    return pins
def limits():
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
def native(pins):
    manifest=json.loads((VARIANT/'manifest.json').read_text())
    library=VARIANT/'libbiogears.so.8.0.0'
    assert sha(library)==manifest['library_sha256']
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
    cwd=build/'projects/biogears/libBiogears'
    objects={name:digest for name,digest in manifest['object_sha256'].items() if 'CircuitCalculator.cpp.o' in name or name.endswith('/Respiratory.cpp.o') or name.endswith('/Controller/BioGears.cpp.o')}
    assert len(objects)>=6
    for name,digest in objects.items():assert sha(cwd/name)==digest,'Native object changed'
    out=Path(tempfile.mkdtemp(prefix='native-respiratory-work-',dir=ROOT/'data/derived/audits'))
    for name in ['native_respiratory_work.h','native_respiratory_work_probe.cpp']:(out/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    binary=out/'fixture'
    command=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_respiratory_work_probe.cpp')]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(p)]
    command+=['-L',str(lib),'-lbiogears','-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(VARIANT)+':'+str(lib)}
    def run(cmd,label):
        with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
            subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources')),'nice','-n','10']+cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,preexec_fn=limits,check=True,timeout=70)
    print(out,flush=True);run(command,'compile')
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
    assert str(library) in linkage,'Wrong native library loaded'
    run([str(binary)],'fixture')
    result=json.loads(next(line[7:] for line in (out/'fixture.stdout').read_text().splitlines() if line.startswith('RESULT ')))
    assert result['passed'] and not result['patient_initialized']
    result['observations']=[json.loads(line[5:]) for line in (out/'fixture.stdout').read_text().splitlines() if line.startswith('CASE ')]
    assert len(result['observations'])==13
    assert any(x['source_work_j']>0 for x in result['observations'])
    assert any(x['source_work_j']<0 for x in result['observations'])
    assert prepare()==pins and sha(library)==manifest['library_sha256']
    result.update(source_sha256=pins,variant=VARIANT.name,manifest_sha256=sha(VARIANT/'manifest.json'),library_sha256=sha(library),native_object_sha256=objects,cdm_library_sha256=sha(lib/'libbiogears_cdm.so'),executable_sha256=sha(binary),compile_command=command,resources={name:(out/(name+'.resources')).read_text() for name in ['compile','fixture']},scope='Isolated native respiratory circuit only; no patient, chemistry, shared adapter integration or physical thorax')
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compile-and-run',action='store_true');a=p.parse_args()
    pins=prepare()
    if a.compile_and_run:native(pins)
    else:print(json.dumps({'prepared':True,'compiled':False,'source_sha256':pins}))
