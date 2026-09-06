"""Scheduled single-unit native operator extraction; caller must own native slot."""
import argparse,hashlib,json,os,shlex,shutil,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def receipt(path,raw):return {'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--run',action='store_true');args=parser.parse_args()
    out=args.output.resolve()
    if out.exists():raise ValueError('Require fresh build directory')
    runtime=ROOT/'data/runtime/opensim';plan=json.loads((ROOT/'data/research/thoracic_mechanism/native_composition_v1/plan.json').read_bytes());model=ROOT/plan['target_native_model_identity']['path'];raw_model=model.read_bytes()
    if hashlib.sha256(raw_model).hexdigest()!=plan['target_native_model_identity']['sha256']:raise ValueError('Frozen target model changed')
    source=ROOT/'scripts/native_thoracic_operator_snapshot.cpp';raw_source=source.read_bytes()
    link_file=runtime/'dynamics-adapter/build/CMakeFiles/native_opensim_dynamics.dir/link.txt';raw_link=link_file.read_bytes();previous=shlex.split(raw_link.decode())
    flags=[v for v in previous if v.endswith('.so') or '.so.' in v or v.startswith('-Wl,') or v.startswith('-l')]
    inputs={source:raw_source,model:raw_model,link_file:raw_link}
    for name in ['SimbodyMatterSubsystem.h','MobilizedBody_Custom.h','Body.h']:
        header=runtime/'install/simbody/include/simbody/simbody/internal'/name;inputs[header]=header.read_bytes()
    for value in flags:
        p=Path(value)
        if p.is_file():inputs[p]=p.read_bytes()
    out.mkdir(parents=True);retained=out/source.name;retained.write_bytes(raw_source);(out/'assembled_model.osim').write_bytes(raw_model)
    compiler=shutil.which('c++');limiter=shutil.which('prlimit')
    if not compiler or not limiter:raise ValueError('Compiler and prlimit required')
    command=[limiter,'--as=2147483648','--cpu=120','nice','-n','10',compiler,'-std=c++20','-O0','-DSWIG_PYTHON']
    for include in ['install/opensim/include','install/opensim/include/OpenSim','install/simbody/include/simbody']:command+=['-isystem',str(runtime/include)]
    executable=out/'snapshot';command+=[str(retained),'-o',str(executable),*flags]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    evidence={'schema':'thoracic-native-snapshot-build-v1','status':'prepared','inputs':[receipt(p,v) for p,v in inputs.items()],'compile_command':command,'native_integrated':False,'native_snapshot_run':False}
    def persist():
        (out/'manifest.json').write_text(json.dumps(evidence,indent=2)+'\n')
    persist()
    try:
        with (out/'compile.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180,env=env)
        evidence['status']='compiled';evidence['executable']=receipt(executable,executable.read_bytes());persist()
        if args.run:
            run=[limiter,'--as=2147483648','--cpu=60','nice','-n','10',str(executable),str(out/'assembled_model.osim'),str(out/'snapshot.json')]
            evidence['run_command']=run;persist()
            with (out/'run.log').open('w') as log:subprocess.run(run,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90,env=env)
            evidence['native_snapshot_run']=True;evidence['snapshot']=receipt(out/'snapshot.json',(out/'snapshot.json').read_bytes())
        if any(p.read_bytes()!=raw for p,raw in inputs.items()):raise ValueError('Bound input changed during extraction')
        evidence['status']='complete';persist()
    except Exception as error:
        evidence['status']='failed';evidence['error']=str(error);persist();raise
    (out/'manifest.json').write_text(json.dumps(evidence,indent=2)+'\n');print(json.dumps({'build':str(out),'native_snapshot_run':args.run}))

if __name__=='__main__':main()
