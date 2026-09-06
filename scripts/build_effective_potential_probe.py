"""Isolated diagnostic build; never updates mechanical latest or frozen archives."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,tempfile,signal
from frozen_static_stream import validate
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def prepare():
    casepath=ROOT/'data/derived/lumbar-supine-reference-1qex2x9i/manifest.json';case=json.loads(casepath.read_text());archive=ROOT/case['frozen_engine'];validate(archive)
    for name,digest in case['files'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Reference changed: '+name)
    output=Path(tempfile.mkdtemp(prefix='effective-potential-build-',dir=ROOT/'data/derived'));source=archive.parent/'sources/native_mechanical_stream.cpp';text=source.read_text()
    include='#include "native_static_pose.h"';anchor='   if(command=="evaluate_static_pose")'
    if text.count(include)!=1 or text.count(anchor)!=1:raise ValueError('Frozen source insertion boundary changed')
    text=text.replace(include,include+'\n#include "native_effective_potential_probe.h"').replace(anchor,'   if(command=="evaluate_effective_pose"){const auto payload=ihm_effective::evaluate(model,state,in,surface_foundation);std::cout<<"@IHM "<<payload<<std::endl;continue;}\n'+anchor)
    copied=output/'native_probe.cpp';copied.write_text(text)
    for p in (archive.parent/'sources').glob('*.h'):shutil.copy2(p,output/p.name)
    header=ROOT/'scripts/native_effective_potential_probe.h';shutil.copy2(header,output/header.name)
    inputsource=ROOT/'data/derived/lumbar-supine-static-1g1q08u2/native/inputs';shutil.copytree(inputsource,output/'inputs')
    if sha(output/'inputs/subject_walk_scaled.osim')!=json.loads((ROOT/case['variant_registration']).read_text())['model_sha256']:raise ValueError('Actual98-muscle input mismatch')
    surface=json.loads((ROOT/case['surface_manifest']).read_text())
    if sha(output/'inputs/supine_surface_foundation.txt')!=surface['native_input_sha256']:raise ValueError('Actual seam-omitted input mismatch')
    origin=json.loads((archive.parent/'origin_build_manifest.json').read_text());command=[]
    for token in origin['command']:
        if token.endswith('/native_mechanical_stream.cpp'):token=str(copied)
        elif token.endswith('/native_mechanical_stream'):token=str(output/'probe')
        elif token.startswith('-Wl,-rpath,'):token='-Wl,-rpath,'+str(archive.parent/'libraries')
        elif '.so' in token and token.startswith('/'):token=str(archive.parent/'libraries'/Path(token).name)
        command.append(token)
    include_roots=[ROOT/'data/runtime/opensim/install/opensim/include',ROOT/'data/runtime/opensim/install/simbody/include/simbody']
    sdk={str(p.relative_to(ROOT)):sha(p) for directory in include_roots for p in directory.rglob('*') if p.is_file()}
    files={str(p.relative_to(ROOT)):sha(p) for p in [casepath,header,source,Path(__file__),ROOT/'scripts/verify_native_effective_potential.py',*output.glob('*.h'),copied,*(output/'inputs').iterdir()]}
    manifest=dict(schema='ihm.effective-potential-probe-build.v1',reference_manifest=str(casepath.relative_to(ROOT)),archive_manifest=str(archive.relative_to(ROOT)),source_changes='Only new diagnostic include/command; physical model builder and continuing engine source otherwise unchanged',files=files,sdk_files=sdk,command=command,executable=str((output/'probe').relative_to(ROOT)),native_run=False)
    path=output/'manifest.json';path.write_text(json.dumps(manifest,indent=2)+'\n');return path

def verify(path):
    m=json.loads(path.read_text());validate(ROOT/m['archive_manifest'])
    for name,digest in (m['files']|m['sdk_files']).items():
        if sha(ROOT/name)!=digest:raise ValueError('Diagnostic build source changed: '+name)
    return m

def build(path):
    m=verify(path)
    with (path.parent/'compile.log').open('w') as log:
        process=subprocess.Popen(m['command'],stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),start_new_session=True)
        try:
            code=process.wait(timeout=60)
            if code:raise subprocess.CalledProcessError(code,m['command'])
        except BaseException:
            if process.poll() is None:os.killpg(process.pid,signal.SIGKILL)
            process.wait();raise
    verify(path);m['executable_sha256']=sha(ROOT/m['executable']);path.write_text(json.dumps(m,indent=2)+'\n');print(path)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--build',type=Path);a=p.parse_args()
    if a.prepare:print(prepare())
    elif a.build:build(a.build.resolve())
    else:raise SystemExit('Choose source preparation or coordinated build')
