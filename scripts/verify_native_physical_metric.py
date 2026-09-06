"""Targeted copied-state sided/step-size work audit; no integration or optimizer."""
from pathlib import Path
import argparse, json, os, signal, subprocess, tempfile, time
import xml.etree.ElementTree as ET
import numpy as np
from build_effective_potential_probe import verify, ROOT, sha
from effective_passive_energy_observation import PassiveEnergyObservation

def audit_metric(response):
    inverse=np.asarray(response['inverse_mass_matrix']);mass=np.asarray(response['mass_matrix'])
    r=np.asarray(response['constrained_residual']);tree=np.asarray(response['tree_residual']);a=np.asarray(response['udot'])
    if inverse.shape!=(len(a),len(a)) or not np.all(np.isfinite(np.r_[inverse.ravel(),mass.ravel(),r,a])):raise ValueError('Invalid native mass observation')
    reconstructed=-inverse@r
    error=float(np.max(np.abs(reconstructed-a)))
    symmetry=float(np.max(np.abs(inverse-inverse.T)))
    identity=float(np.max(np.abs(mass@inverse-np.eye(len(a)))))
    minimum_eigenvalue=float(np.linalg.eigvalsh(mass)[0])
    return dict(maximum_acceleration_identity_error=error,inverse_symmetry_error=symmetry,mass_inverse_identity_error=identity,
        minimum_mass_eigenvalue=minimum_eigenvalue,tree_without_constraint_acceleration_error=float(np.max(np.abs(-inverse@tree-a))),
        passed=bool(error<=1e-8 and symmetry<=1e-8 and identity<=1e-8 and minimum_eigenvalue>0),
        acceleration_units='Each translational component m/s^2 and angular component rad/s^2; diagnostic identity error 1e-8 each, physical acceptance remains1e-4')


def run(build_path):
    build=verify(build_path);exe=ROOT/build['executable']
    if sha(exe)!=build.get('executable_sha256'):raise ValueError('Unverified executable')
    folder=build_path.parent;archive=ROOT/build['archive_manifest'];ar=json.loads(archive.read_text())
    case=json.loads((ROOT/build['reference_manifest']).read_text())
    bounds={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.parse(folder/'inputs/subject_walk_scaled.osim').getroot().iter('Coordinate')}
    port=PassiveEnergyObservation(folder/'inputs/subject_walk_scaled.osim')
    output=Path(tempfile.mkdtemp(prefix='native-physical-metric-',dir=ROOT/'data/derived'))
    records=[];process=None;started=time.monotonic()
    report=dict(maximum_evaluations=2,maximum_wall_s=15,build_manifest=str(build_path.relative_to(ROOT)),build_manifest_sha256=sha(build_path),physical_time_advanced_s=0,accepted_equilibrium=False)
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(*args):
        if process is not None and process.poll() is None:os.killpg(process.pid,signal.SIGKILL)
        raise TimeoutError('15 s metric probe cap')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,15)
    try:
        command=['prlimit','--as=4294967296','--','nice','-n','10',str(archive.parent/ar['dynamic_loader']),'--library-path',str(archive.parent/'libraries'),str(exe),str(folder/'inputs'),str(output),'supine',str(case['mass_kg'])]
        log=(output/'engine.log').open('w')
        process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=log,text=True,bufsize=1,start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',LD_LIBRARY_PATH=str(archive.parent/'libraries')))
        def read():
            while True:
                line=process.stdout.readline()
                if not line:raise RuntimeError('Native probe terminated')
                if line.startswith('@IHM '):
                    value=json.loads(line[5:])
                    if 'error' in value:raise ValueError(value['error'])
                    return value
                log.write(line);log.flush()
        initial=read();write('initial.json',initial)
        def evaluate(q,label):
            if len(records)>=2:raise RuntimeError('2-call cap')
            write('pending_request.json',dict(q=q,label=label))
            fields=['evaluate_effective_pose',str(len(q))]
            for n,v in q.items():fields.extend([n,str(v)])
            process.stdin.write(' '.join(fields)+'\n');process.stdin.flush();value=read()
            record=dict(label=label,requested=q,response=value,passive_energy_observation=port.observe(value));records.append(record)
            with (output/'observations.jsonl').open('a') as stream:stream.write(json.dumps(record,allow_nan=False)+'\n')
            return value
        q=json.loads((ROOT/case['old_q_seed']).read_text())['coordinates']
        base=evaluate(q,'supported_old_q_seed')
        failed=json.loads((ROOT/'data/derived/native-force-root-919h2o1p/last_accepted.json').read_text())['requested']
        trial=evaluate(failed,'unsupported_failed_force_root_q')
        process.stdin.write('observe\n');process.stdin.flush();live=read()
        if any(initial[k]!=live[k] for k in initial if k!='kind'):raise ValueError('Continuing state changed')
        audits={name:audit_metric(value) for name,value in [('supported_seed',base),('unsupported_force_root',trial)]}
        report.update(status='diagnostic_complete',continuing_state_unchanged=True,audits=audits,
            metric_gate_passed=all(a['passed'] for a in audits.values()),interpretation='Native mass/sign/constraint identity only; no equilibrium or optimizer acceptance')
    except Exception as error:report.update(status='diagnostic_failed',error=type(error).__name__+': '+str(error)[:65536])
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if process is not None:
            try:
                if process.poll() is None:process.stdin.write('close\n');process.stdin.flush();process.wait(timeout=2)
            except Exception:
                if process.poll() is None:os.killpg(process.pid,signal.SIGKILL)
                process.wait(timeout=2)
        report.update(evaluations=len(records),wall_s=time.monotonic()-started,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'scripts/effective_passive_energy_observation.py')})
        write('report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='rows'}|{'output':str(output)},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');parser.add_argument('build',type=Path);args=parser.parse_args()
    if not args.run_native:raise SystemExit('Coordinated native grant required')
    run(args.build.resolve())
