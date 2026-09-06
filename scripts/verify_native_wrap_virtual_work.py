"""Targeted copied-state sided/step-size work audit; no integration or optimizer."""
from pathlib import Path
import argparse, json, os, signal, subprocess, tempfile, time
import xml.etree.ElementTree as ET
import numpy as np
from build_effective_potential_probe import verify, ROOT, sha
from effective_passive_energy_observation import PassiveEnergyObservation

COORDINATES = [stem + side for side in ('r','l') for stem in ('arm_flex_','arm_add_','arm_rot_','elbow_flex_')]
STEPS = (1e-3, 1e-4, 1e-5)


def directional(base, shifted, coordinate, port):
    names=base['independent_names'];i=names.index(coordinate)
    dq=np.array(shifted['independent_q'])-base['independent_q'];h=dq[i]
    if h==0:raise ValueError('Assembly removed probe direction')
    result=[]
    for b,s in zip(base['muscles'],shifted['muscles']):
        if b['name']!=s['name']:raise ValueError('Muscle ordering changed')
        if not b['name'].startswith(('arm26_BIClong_', 'arm26_BRA_')):continue
        ma=np.array(b['moment_arms_m']);sm=np.array(s['moment_arms_m'])
        dl=s['length_m']-b['length_m']
        eb=port.muscle(b)['corrected_passive_energy_j']+b['active_effective_energy_j']
        es=port.muscle(s)['corrected_passive_energy_j']+s['active_effective_energy_j']
        result.append(dict(name=b['name'],length_difference_m=dl,
            length_derivative_m_per_rad=dl/h,
            base_virtual_work_error_m_per_rad=(dl+ma@dq)/h,
            trapezoid_virtual_work_error_m_per_rad=(dl+.5*(ma+sm)@dq)/h,
            effective_energy_length_work_error_nm=((es-eb)-.5*(b['tendon_force_n']+s['tendon_force_n'])*dl)/h,
            base_moment_arm_m=float(ma[i]), endpoint_moment_arm_m=float(sm[i]),
            minimum_fiber_margin_m=min(b['fiber_length_m']-b['minimum_fiber_length_m'],s['fiber_length_m']-s['minimum_fiber_length_m']),
            minimum_internal_stiffness_n_m=min(b['internal_stiffness_n_m'],s['internal_stiffness_n_m'])))
    offaxis=dq.copy();offaxis[i]=0
    return dict(actual_step_rad=h,maximum_off_axis_q_change=float(np.max(np.abs(offaxis))),
        muscles=result,branch_scope='Sided force/path continuity observations; no explicit native wrap branch identifier emitted')


def run(build_path):
    build=verify(build_path);exe=ROOT/build['executable']
    if sha(exe)!=build.get('executable_sha256'):raise ValueError('Unverified executable')
    folder=build_path.parent;archive=ROOT/build['archive_manifest'];ar=json.loads(archive.read_text())
    case=json.loads((ROOT/build['reference_manifest']).read_text())
    bounds={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.parse(folder/'inputs/subject_walk_scaled.osim').getroot().iter('Coordinate')}
    port=PassiveEnergyObservation(folder/'inputs/subject_walk_scaled.osim')
    output=Path(tempfile.mkdtemp(prefix='native-wrap-work-',dir=ROOT/'data/derived'))
    records=[];process=None;started=time.monotonic()
    report=dict(maximum_evaluations=50,maximum_wall_s=45,build_manifest=str(build_path.relative_to(ROOT)),build_manifest_sha256=sha(build_path),physical_time_advanced_s=0,accepted_equilibrium=False)
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(*args):
        if process is not None and process.poll() is None:os.killpg(process.pid,signal.SIGKILL)
        raise TimeoutError('45 s native audit cap')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,45)
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
            if len(records)>=50:raise RuntimeError('50-call cap')
            write('pending_request.json',dict(q=q,label=label))
            fields=['evaluate_effective_pose',str(len(q))]
            for n,v in q.items():fields.extend([n,str(v)])
            process.stdin.write(' '.join(fields)+'\n');process.stdin.flush();value=read()
            record=dict(label=label,requested=q,response=value,passive_energy_observation=port.observe(value));records.append(record)
            with (output/'observations.jsonl').open('a') as stream:stream.write(json.dumps(record,allow_nan=False)+'\n')
            return value
        # Start from the exact original requested seed, not re-requested assembled q.
        q=json.loads((ROOT/case['old_q_seed']).read_text())['coordinates'];base=evaluate(q,'base');rows=[]
        for coordinate in COORDINATES:
            for step in STEPS:
                for sign in (-1,1):
                    shifted=dict(q);shifted[coordinate]+=sign*step
                    if not bounds[coordinate][0]<=shifted[coordinate]<=bounds[coordinate][1]:raise ValueError('Probe exceeds source bound')
                    value=evaluate(shifted,f'{coordinate}:{sign*step:.8g}')
                    rows.append(dict(coordinate=coordinate,requested_step_rad=sign*step,**directional(base,value,coordinate,port)))
        repeat=evaluate(q,'repeat_base');process.stdin.write('observe\n');process.stdin.flush();live=read()
        if any(initial[k]!=live[k] for k in initial if k!='kind'):raise ValueError('Continuing state changed')
        report.update(status='diagnostic_complete',continuing_state_unchanged=True,rows=rows,
            repeated_request_maximum_actual_q_difference=float(np.max(np.abs(np.array(base['independent_q'])-repeat['independent_q']))),
            repeated_request_energy_difference_j=repeat['effective_energy_j']-base['effective_energy_j'],
            interpretation='Diagnostic only: step convergence and sided differences must be reviewed before any gradient gate acceptance')
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
