"""Replay identical retained native MM inputs against two receipt-bound binaries.

This is a historical binary comparison, not a current-adapter manifest bypass.
Both compiled sources and binaries remain retained and hashed in the output.
"""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, tempfile, time
ROOT=Path(__file__).resolve().parents[1]


def run(execution_path,new_build,duration):
    execution=json.loads(execution_path.read_text())
    inputs=Path(execution['engine_command'][1])
    for name,digest in execution['source_sha256'].items():
        assert hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest,name
    out=Path(tempfile.mkdtemp(prefix='native-bed-solver-parity-',dir=ROOT/'data/derived'))
    report={'passed':False,'scope':'Identical retained initialized MM body, no external ports, held source excitation defaults',
            'execution_source':str(execution_path),'duration_s':duration,'arms':{}}
    for label,build in [('bisection',ROOT/execution['build']['path']),('newton',new_build)]:
        destination=out/label;destination.mkdir()
        binary=build/'native_mechanical_stream';manifest=json.loads((build/'manifest.json').read_text())
        # Historical source paths may intentionally differ now. Verify all retained
        # build-local bytes, including binary/header/source, without retagging them.
        for path,digest in manifest['files'].items():
            target=ROOT/path
            if target.parent==build:
                assert hashlib.sha256(target.read_bytes()).hexdigest()==digest,path
        command=[str(binary),execution['engine_command'][1],str(destination),*execution['engine_command'][3:]]
        started=time.monotonic();rows=[]
        with (destination/'engine.log').open('w') as log:
            process=subprocess.Popen(['prlimit','--as=4294967296','--',*command],stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=log,text=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            def read():
                while True:
                    line=process.stdout.readline()
                    if not line:raise RuntimeError('Native terminated; see '+str(destination/'engine.log'))
                    if line.startswith('@IHM '):
                        value=json.loads(line[5:])
                        if 'error' in value:raise RuntimeError(value['error'])
                        return value
                    log.write(line)
            try:
                rows.append(read())
                for _ in range(round(duration/.005)):
                    process.stdin.write('advance .005 0 0\n');process.stdin.flush();rows.append(read())
            finally:
                if process.poll() is None:
                    process.stdin.write('quit\n');process.stdin.flush()
                    try:process.wait(timeout=5)
                    except subprocess.TimeoutExpired:process.kill();process.wait()
        (destination/'states.json').write_text(json.dumps(rows))
        report['arms'][label]={'build':str(build),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'wall_s':time.monotonic()-started,'states':str(destination/'states.json')}
        print(json.dumps({'arm':label,**report['arms'][label]}),flush=True)
    a=json.loads(Path(report['arms']['bisection']['states']).read_text());b=json.loads(Path(report['arms']['newton']['states']).read_text())
    errors={'coordinate':0.,'kinetic_energy_j':0.,'potential_energy_j':0.,'surface_numeric':0.}
    def numeric_error(x,y):
        if isinstance(x,dict):return max([0.,*[numeric_error(v,y[k]) for k,v in x.items() if k in y]])
        if isinstance(x,list):return max([0.,*[numeric_error(v,w) for v,w in zip(x,y)]])
        if type(x) in (int,float) and type(y) in (int,float):return abs(x-y)
        return 0.
    for x,y in zip(a,b):
        assert x['time_s']==y['time_s']
        errors['coordinate']=max(errors['coordinate'],max(abs(v['value']-y['coordinates'][k]['value']) for k,v in x['coordinates'].items()))
        for k in ('kinetic_energy_j','potential_energy_j'):errors[k]=max(errors[k],abs(x[k]-y[k]))
        errors['surface_numeric']=max(errors['surface_numeric'],numeric_error(x['surface_foundation'],y['surface_foundation']))
    report['maximum_absolute_errors']=errors
    report['thresholds']={'coordinate':1e-7,'kinetic_energy_j':1e-5,'potential_energy_j':1e-5,'surface_numeric':1e-4}
    report['passed']=all(errors[k]<=v for k,v in report['thresholds'].items())
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'output':str(out),**report},indent=2))
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('execution',type=Path);parser.add_argument('--new-build',type=Path);parser.add_argument('--duration',type=float,default=.02);args=parser.parse_args()
    build=args.new_build or ROOT/json.loads((ROOT/'data/runtime/mechanical-stream/latest.json').read_text())['build']
    raise SystemExit(0 if run(args.execution,build,args.duration)['passed'] else 1)
