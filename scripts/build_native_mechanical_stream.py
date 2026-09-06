"""One low-priority translation-unit build against retained OpenSim libraries."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shlex,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.native.instance_mass import variant_identity,pointer_directory

def build(instance_mass_variant=None):
    mass_identity=None if instance_mass_variant is None else variant_identity(ROOT,instance_mass_variant)
    base=ROOT/'data/runtime/mechanical-stream';base.mkdir(parents=True,exist_ok=True)
    out=Path(tempfile.mkdtemp(prefix='build-',dir=base));source=ROOT/'scripts/native_mechanical_stream.cpp';retained=out/source.name;retained.write_bytes(source.read_bytes())
    header=ROOT/'scripts/native_muscle_metabolism.h';(out/header.name).write_bytes(header.read_bytes())
    pose_header=ROOT/'scripts/native_static_pose.h';(out/pose_header.name).write_bytes(pose_header.read_bytes())
    surface_header=ROOT/'scripts/native_surface_foundation.h';(out/surface_header.name).write_bytes(surface_header.read_bytes())
    bed_header=ROOT/'scripts/native_bed_compression.h';(out/bed_header.name).write_bytes(bed_header.read_bytes())
    mass_header=ROOT/'scripts/native_local_mass_port.h';(out/mass_header.name).write_bytes(mass_header.read_bytes())
    limiter=shutil.which('prlimit')
    if not limiter:raise RuntimeError('prlimit required for bounded compilation')
    runtime=ROOT/'data/runtime/opensim';flags=['-std=c++20','-O0','-DSWIG_PYTHON']
    for p in [runtime/'install/opensim/include',runtime/'install/opensim/include/OpenSim',runtime/'install/simbody/include/simbody']:flags+=['-isystem',str(p)]
    previous=shlex.split((runtime/'dynamics-adapter/build/CMakeFiles/native_opensim_dynamics.dir/link.txt').read_text())
    libraries=[v for v in previous if v.endswith('.so') or '.so.' in v or v.startswith('-Wl,') or v.startswith('-l')]
    library_paths=[Path(v) for v in libraries if Path(v).is_file()]
    if mass_identity is not None:library_paths += [ROOT/mass_identity['library_path'],ROOT/mass_identity['path']/'manifest.json']
    library_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in library_paths}
    command=[limiter,'--as=4294967296','--','nice','-n','10','c++',*flags,str(retained),'-o',str(out/'native_mechanical_stream'),*libraries,'-ldl']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    started=time.monotonic()
    with (out/'compile.log').open('w') as log:result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError('Native stream compile failed: '+str(out/'compile.log'))
    if source.read_bytes()!=retained.read_bytes() or header.read_bytes()!=(out/header.name).read_bytes():raise RuntimeError('Source changed during native compilation; unpublished build retained')
    if bed_header.read_bytes()!=(out/bed_header.name).read_bytes():raise RuntimeError('Bed compression header changed during compilation')
    if surface_header.read_bytes()!=(out/surface_header.name).read_bytes():raise RuntimeError('Surface foundation header changed during compilation')
    if pose_header.read_bytes()!=(out/pose_header.name).read_bytes():raise RuntimeError('Static pose header changed during compilation')
    if mass_header.read_bytes()!=(out/mass_header.name).read_bytes():raise RuntimeError('Mass port header changed during compilation')
    if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=digest for p,digest in library_hashes.items()):raise RuntimeError('Library changed during native compilation')
    paths=[source,header,mass_header,out/mass_header.name,bed_header,out/bed_header.name,surface_header,out/surface_header.name,pose_header,out/pose_header.name,out/'native_mechanical_stream',retained,out/header.name]+library_paths
    manifest={'schema':'ihm.mechanical-stream-build.v1','instance_mass_variant':mass_identity,'path':str(out.relative_to(ROOT)),
              'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              'command':command,'compile_wall_s':time.monotonic()-started,'maximum_compiler_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'address_space_limit_bytes':4294967296,'header_snapshot_sha256':hashlib.sha256((out/header.name).read_bytes()).hexdigest(),'source_snapshot_sha256':hashlib.sha256(retained.read_bytes()).hexdigest()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    pointer_base=base if mass_identity is None else pointer_directory(ROOT,mass_identity);pointer_base.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w',dir=pointer_base,delete=False) as pointer:
        pointer.write(json.dumps({'build':str(out.relative_to(ROOT))})+'\n');pointer_name=pointer.name
    os.replace(pointer_name,pointer_base/'latest.json');print(json.dumps({'passed':True,'output_dir':str(out)},indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--instance-mass-variant');args=parser.parse_args();build(args.instance_mass_variant)
