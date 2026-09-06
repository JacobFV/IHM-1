"""Bounded positive mass/flux acceptance for an explicitly linked variant."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shutil,subprocess,tempfile,time
from build_simbody_instance_mass_variant import ROOT,validate

def main():
    p=argparse.ArgumentParser();p.add_argument('--variant',type=Path,required=True);p.add_argument('--run-native',action='store_true');a=p.parse_args()
    if not a.run_native:raise SystemExit('Requires coordinated --run-native slot')
    variant=(ROOT/a.variant).resolve();manifest=json.loads((variant/'manifest.json').read_text());assert manifest['complete'],'library not linked'
    validate(variant,manifest);library=variant/'libSimTKsimbody.so.3.9';assert hashlib.sha256(library.read_bytes()).hexdigest()==manifest['library_sha256']
    out=Path(tempfile.mkdtemp(prefix='instance-mass-positive-',dir=ROOT/'data/derived'))
    for name in ('native_simbody_instance_mass_positive_fixture.cpp','native_simbody_instance_mass.h'):(out/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    source_hashes={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ('native_simbody_instance_mass_positive_fixture.cpp','native_simbody_instance_mass.h')}
    runtime=ROOT/'data/runtime/opensim';lib=runtime/'install/simbody/lib';limiter=shutil.which('prlimit');assert limiter
    command=[limiter,'--as=4294967296','--','nice','-n','10','c++','-std=c++20','-O1','-isystem',str(runtime/'install/simbody/include/simbody'),str(out/'native_simbody_instance_mass_positive_fixture.cpp'),str(library),'-L',str(lib),'-lSimTKmath','-lSimTKcommon','-o',str(out/'fixture')]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LD_LIBRARY_PATH=':'.join(map(str,[variant,lib,*[runtime/'sysroot/usr/lib/aarch64-linux-gnu'/x for x in ('lapack','blas','')]])))
    (out/'command.json').write_text(json.dumps(command,indent=2)+'\n');start=time.monotonic()
    with (out/'compile.log').open('w') as log:result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
    if result.returncode:raise RuntimeError('Compile failed: '+str(out/'compile.log'))
    validate(variant,manifest)
    assert all(hashlib.sha256((ROOT/'scripts'/name).read_bytes()).hexdigest()==digest for name,digest in source_hashes.items()),'fixture source changed during build'
    loaded=subprocess.run(['ldd',str(out/'fixture')],capture_output=True,text=True,env=env);(out/'loaded_libraries.txt').write_text(loaded.stdout)
    assert str(library) in loaded.stdout,'variant library not selected'
    loaded_paths=[Path(token) for token in loaded.stdout.split() if token.startswith('/') and Path(token).is_file()]
    loaded_hashes={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in loaded_paths}
    run=subprocess.run([limiter,'--as=4294967296','--','nice','-n','10',str(out/'fixture')],capture_output=True,text=True,env=env)
    (out/'stdout.log').write_text(run.stdout);(out/'stderr.log').write_text(run.stderr)
    if run.returncode:raise RuntimeError('Positive native mass fixture failed: '+str(out/'stderr.log'))
    assert all(hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest for name,digest in loaded_hashes.items()),'loaded library changed during fixture'
    report={**json.loads(run.stdout),'loaded_library_sha256':loaded_hashes,'fixture_source_sha256':source_hashes,'variant':str(variant.relative_to(ROOT)),'variant_library_sha256':manifest['library_sha256'],'wall_s':time.monotonic()-start,'max_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'scope':'State-instance mass, local free-body mass matrix/COM/gravity and momentum/kinetic capture/outflow, concurrent State/replay. Not an OpenSim muscle or whole-body constraint acceptance.'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({**report,'output':str(out)},indent=2))
if __name__=='__main__':main()
