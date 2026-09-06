#!/usr/bin/env python3
"""One full native TU per granted stage; isolated objects, no shared library edits."""
import argparse,json,os,resource,subprocess
from pathlib import Path
from prepare_gi_serialization_repair import prepare,BASE,FILES
from build_biogears_shared_donor_variant import ROOT,RUNTIME,sha

UNITS={'outer':FILES[0],'nested':FILES[1],'ownership':FILES[2]}
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--unit',choices=UNITS,required=True);a=p.parse_args()
    out=a.output.resolve();prepare(out);source=out/UNITS[a.unit];obj=out/(a.unit+'.o')
    if obj.exists():raise ValueError('Retain immutable existing object; use new output directory')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';flags=cwd/'CMakeFiles/libbiogears.dir/flags.make';includes=cwd/'CMakeFiles/libbiogears.dir/includes_CXX.rsp'
    import shlex
    raw={line.split(' = ',1)[0]:line.split(' = ',1)[1] for line in flags.read_text().splitlines() if ' = ' in line}
    command=['c++',*shlex.split(raw['CXX_DEFINES']),*shlex.split(raw['CXX_INCLUDES']),*shlex.split(raw['CXX_FLAGS']),'-I',str((BASE/UNITS[a.unit]).parent),'-c',str(source),'-o',str(obj)]
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_CPU,(90,90));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    with (out/(a.unit+'.stdout')).open('w') as stdout,(out/(a.unit+'.stderr')).open('w') as stderr:
        result=subprocess.run(['/usr/bin/time','-v','-o',str(out/(a.unit+'.resources'))]+command,cwd=cwd,env=env,stdout=stdout,stderr=stderr,timeout=90,preexec_fn=limits)
    if result.returncode:raise RuntimeError(f'Full TU compile failed; receipt {out}')
    path=out/'objects_manifest.json';manifest=json.loads(path.read_text()) if path.exists() else {'objects':{}}
    key='CMakeFiles/libbiogears.dir/'+UNITS[a.unit]+'.o'
    manifest['objects'][key]=dict(object=str(obj),sha256=sha(obj),patched_source=str(source),patched_source_sha256=sha(source),original_source_sha256=sha(BASE/UNITS[a.unit]),command=command,flags_sha256=sha(flags),include_response_sha256=sha(includes),resources=(out/(a.unit+'.resources')).read_text())
    path.write_text(json.dumps(manifest,indent=2)+'\n');print(path)
if __name__=='__main__':main()
