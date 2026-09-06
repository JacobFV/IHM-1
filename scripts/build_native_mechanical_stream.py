"""One low-priority translation-unit build against retained OpenSim libraries."""
from pathlib import Path
import hashlib,json,os,resource,shlex,shutil,subprocess,tempfile,time
ROOT=Path(__file__).resolve().parents[1]

def build():
    base=ROOT/'data/runtime/mechanical-stream';base.mkdir(parents=True,exist_ok=True)
    out=Path(tempfile.mkdtemp(prefix='build-',dir=base));source=ROOT/'scripts/native_mechanical_stream.cpp';retained=out/source.name;retained.write_bytes(source.read_bytes())
    header=ROOT/'scripts/native_muscle_metabolism.h';(out/header.name).write_bytes(header.read_bytes())
    pose_header=ROOT/'scripts/native_static_pose.h';(out/pose_header.name).write_bytes(pose_header.read_bytes())
    limiter=shutil.which('prlimit')
    if not limiter:raise RuntimeError('prlimit required for bounded compilation')
    runtime=ROOT/'data/runtime/opensim';flags=['-std=c++20','-O0','-DSWIG_PYTHON']
    for p in [runtime/'install/opensim/include',runtime/'install/opensim/include/OpenSim',runtime/'install/simbody/include/simbody']:flags+=['-isystem',str(p)]
    previous=shlex.split((runtime/'dynamics-adapter/build/CMakeFiles/native_opensim_dynamics.dir/link.txt').read_text())
    libraries=[v for v in previous if v.endswith('.so') or '.so.' in v or v.startswith('-Wl,') or v.startswith('-l')]
    library_paths=[Path(v) for v in libraries if Path(v).is_file()]
    library_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in library_paths}
    command=[limiter,'--as=4294967296','--','nice','-n','10','c++',*flags,str(retained),'-o',str(out/'native_mechanical_stream'),*libraries]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    started=time.monotonic()
    with (out/'compile.log').open('w') as log:result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError('Native stream compile failed: '+str(out/'compile.log'))
    if source.read_bytes()!=retained.read_bytes() or header.read_bytes()!=(out/header.name).read_bytes():raise RuntimeError('Source changed during native compilation; unpublished build retained')
    if pose_header.read_bytes()!=(out/pose_header.name).read_bytes():raise RuntimeError('Static pose header changed during compilation')
    if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=digest for p,digest in library_hashes.items()):raise RuntimeError('Library changed during native compilation')
    paths=[source,header,pose_header,out/pose_header.name,out/'native_mechanical_stream',retained,out/header.name]+library_paths
    manifest={'schema':'ihm.mechanical-stream-build.v1','path':str(out.relative_to(ROOT)),
              'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              'command':command,'compile_wall_s':time.monotonic()-started,'maximum_compiler_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'address_space_limit_bytes':4294967296,'header_snapshot_sha256':hashlib.sha256((out/header.name).read_bytes()).hexdigest(),'source_snapshot_sha256':hashlib.sha256(retained.read_bytes()).hexdigest()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with tempfile.NamedTemporaryFile(mode='w',dir=base,delete=False) as pointer:
        pointer.write(json.dumps({'build':str(out.relative_to(ROOT))})+'\n');pointer_name=pointer.name
    os.replace(pointer_name,base/'latest.json');print(json.dumps({'passed':True,'output_dir':str(out)},indent=2))
if __name__=='__main__':build()
