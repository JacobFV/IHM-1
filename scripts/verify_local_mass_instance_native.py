"""Build/run isolated private-ABI local mass fixture in coordinated native slot."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shutil,subprocess,tempfile,time
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');args=parser.parse_args()
    if not args.run_native:raise SystemExit('Requires coordinated --run-native compile/run slot')
    out=Path(tempfile.mkdtemp(prefix='local-mass-instance-',dir=ROOT/'data/derived'))
    runtime=ROOT/'data/runtime/opensim';private=ROOT/'data/raw/mechanics/simbody/Simbody/src'
    source=ROOT/'scripts/native_local_mass_instance_fixture.cpp';retained=out/source.name;retained.write_bytes(source.read_bytes())
    paths=[source,*sorted(private.glob('*.h')),*sorted((runtime/'install/simbody/lib').glob('*.so.3.9'))]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    limiter=shutil.which('prlimit');assert limiter,'prlimit required'
    lib=runtime/'install/simbody/lib'
    command=[limiter,'--as=4294967296','--','nice','-n','10','c++','-std=c++20','-O0','-isystem',str(runtime/'install/simbody/include/simbody'),'-I',str(private),str(retained),'-L',str(lib),'-Wl,-rpath,'+str(lib),'-lSimTKsimbody','-lSimTKmath','-lSimTKcommon','-o',str(out/'fixture')]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LD_LIBRARY_PATH=':'.join(map(str,[lib,*[runtime/'sysroot/usr/lib/aarch64-linux-gnu'/p for p in ('lapack','blas','')]])))
    start=time.monotonic();(out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    with (out/'compile.log').open('w') as log:compiled=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
    if compiled.returncode:raise RuntimeError('Native local mass compile failed: '+str(out/'compile.log'))
    if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=digest for p,digest in hashes.items()):raise RuntimeError('Source/library changed during compile')
    run=subprocess.run([limiter,'--as=4294967296','--','nice','-n','10',str(out/'fixture')],capture_output=True,text=True,env=env)
    (out/'stdout.log').write_text(run.stdout);(out/'stderr.log').write_text(run.stderr)
    if run.returncode:raise RuntimeError('Native local mass fixture failed: '+str(out/'stderr.log'))
    result=json.loads(run.stdout);assert result['passed']
    report={**result,'source_and_library_sha256':hashes,'executable_sha256':hashlib.sha256((out/'fixture').read_bytes()).hexdigest(),'wall_s':time.monotonic()-start,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'scope':'Falsification of private Simbody instance update: reported local mass changes but dynamic mass matrix and kinetic energy do not; state-copy replay, expected flux mismatch. No runtime mass port.'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({**result,'output':str(out),'maximum_child_rss_kib':report['maximum_child_rss_kib']}))
if __name__=='__main__':main()
