"""Gated actual-force trust region on separately attested98-muscle physics."""
from pathlib import Path
import argparse, json, os, signal, subprocess, tempfile, time
import xml.etree.ElementTree as ET
import numpy as np
from build_effective_potential_probe import verify, ROOT, sha
from effective_passive_energy_observation import PassiveEnergyObservation

from prepare_native_force_scaling import model_info, force, GAUGES
from native_force_trust_region import sensitivity_scales, bounded_step, assess_trial


def physical_gate(response,weight,height):
    names=response['mobility_names'];r=np.asarray(response['constrained_residual'])
    if weight<=0 or height<=0 or not np.all(np.isfinite(np.r_[r,response['udot'],weight,height])):raise ValueError('Nonfinite physical gate input')
    def normalized(n):return float(r[names.index(n)]/(weight if n in ('pelvis_tx','pelvis_ty','pelvis_tz') else weight*height))
    support=[normalized(n) for n in ('pelvis_tx','pelvis_tilt','pelvis_rotation')]
    gauges=[normalized(n) for n in GAUGES]
    maximum=max(map(abs,response['udot']))
    constraint=max(abs(response[k]) for k in ('constraint_position_error','constraint_velocity_error','constraint_acceleration_error'))
    return dict(maximum_abs_acceleration=maximum,support_constraints=support,gauge_residual=gauges,
        maximum_constraint_error=constraint,passed=bool(maximum<=1e-4 and max(map(abs,support+gauges))<=1e-4 and constraint<=1e-5))


def run(build_path):
    build=verify(build_path);exe=ROOT/build['executable']
    if sha(exe)!=build.get('executable_sha256'):raise ValueError('Unverified executable')
    folder=build_path.parent;archive=ROOT/build['archive_manifest'];ar=json.loads(archive.read_text())
    case=json.loads((ROOT/build['reference_manifest']).read_text())
    bounds={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.parse(folder/'inputs/subject_walk_scaled.osim').getroot().iter('Coordinate')}
    port=PassiveEnergyObservation(folder/'inputs/subject_walk_scaled.osim')
    output=Path(tempfile.mkdtemp(prefix='native-force-root-',dir=ROOT/'data/derived'))
    records=[];process=None;started=time.monotonic()
    report=dict(maximum_evaluations=200,maximum_wall_s=60,build_manifest=str(build_path.relative_to(ROOT)),build_manifest_sha256=sha(build_path),physical_time_advanced_s=0,accepted_equilibrium=False)
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(*args):
        if process is not None and process.poll() is None:os.killpg(process.pid,signal.SIGKILL)
        raise TimeoutError('60 s native force-root cap')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,60)
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
            if len(records)>=200:raise RuntimeError('200-call cap')
            write('pending_request.json',dict(q=q,label=label))
            fields=['evaluate_effective_pose',str(len(q))]
            for n,v in q.items():fields.extend([n,str(v)])
            process.stdin.write(' '.join(fields)+'\n');process.stdin.flush()
            try:value=read()
            except Exception as error:
                rejection=dict(label=label,requested=q,rejected=True,error=type(error).__name__+': '+str(error));records.append(rejection)
                with (output/'observations.jsonl').open('a') as stream:stream.write(json.dumps(rejection,allow_nan=False)+'\n')
                raise
            record=dict(label=label,requested=q,response=value,passive_energy_observation=port.observe(value));records.append(record)
            with (output/'observations.jsonl').open('a') as stream:stream.write(json.dumps(record,allow_nan=False)+'\n')
            return value
        q=json.loads((ROOT/case['old_q_seed']).read_text())['coordinates']
        current=evaluate(q,'base');names=current['independent_names'];free=[n for n in names if n not in GAUGES]
        indices=[names.index(n) for n in free];B=np.array([bounds[n] for n in free]);_,couplings=model_info(folder/'inputs/subject_walk_scaled.osim')
        surface=json.loads((ROOT/case['surface_manifest']).read_text())
        height=float(np.ptp(np.load(ROOT/surface['arrays_path'])['reference_points_source_m'][:,1]))
        gravity=current['gravity_generalized_forces'];weight=float(np.linalg.norm([gravity[current['mobility_names'].index(n)] for n in ('pelvis_tx','pelvis_ty','pelvis_tz')]))
        coordinate_scale=residual_scale=None;radius=.01;iterations=[]
        write('initial_gate.json',physical_gate(current,weight,height))
        for iteration in range(6):
            gate=physical_gate(current,weight,height)
            if gate['passed']:report['accepted_equilibrium']=True;break
            q=dict(zip(names,current['independent_q']));x=np.array([q[n] for n in free]);r=force(current,couplings)
            D=[];Y=[];maximum_gauge_drift=0.
            for i,n in enumerate(free):
                shifted=dict(q);h=1e-5 if q[n]+1e-5<=bounds[n][1] else -1e-5
                shifted[n]+=h
                if shifted[n]<bounds[n][0]:raise ValueError('No bounded Jacobian direction')
                sample=evaluate(shifted,f'jacobian:{iteration}:{n}')
                D.append(np.array(sample['independent_q'])[indices]-x);Y.append(force(sample,couplings)-r)
                maximum_gauge_drift=max(maximum_gauge_drift,max(abs(sample['independent_q'][names.index(g)]-q[g]) for g in GAUGES))
            condition=float(np.linalg.cond(D))
            if condition>100 or maximum_gauge_drift>1e-8:raise ValueError('Unresolved assembled Jacobian coordinate map')
            J=np.linalg.solve(D,Y).T
            if coordinate_scale is None:
                coordinate_scale,residual_scale,basis=sensitivity_scales(J,B)
                write('scales.json',dict(coordinate_names=free,residual_names=names,coordinate_scale=coordinate_scale.tolist(),residual_scale=residual_scale.tolist(),basis=basis))
            write('last_jacobian.json',dict(iteration=iteration,q=q,force=r.tolist(),jacobian=J.tolist(),condition=condition,maximum_gauge_drift=maximum_gauge_drift))
            accepted=False
            for attempt in range(8):
                step,diagnostic=bounded_step(J,r,x,B,coordinate_scale,residual_scale,radius)
                if diagnostic['predicted_reduction']<=0:break
                trialq=dict(q);trialq.update(zip(free,(x+step).tolist()))
                try:trial=evaluate(trialq,f'trial:{iteration}:{attempt}')
                except ValueError as error:
                    if str(error)!='equal-pressure skin/bed solution exceeds retained domains':raise
                    iterations.append(dict(iteration=iteration,attempt=attempt,radius=radius,accepted=False,domain_rejection=str(error)));radius*=.5;continue
                actual_step=np.array(trial['independent_q'])[indices]-x
                assessment=assess_trial(r,force(trial,couplings),J,actual_step,residual_scale)
                entry=dict(iteration=iteration,attempt=attempt,diagnostic=diagnostic,assessment=assessment,physical_gate=physical_gate(trial,weight,height))
                iterations.append(entry);write('iterations.json',iterations)
                if assessment['accepted']:
                    current=trial;accepted=True
                    write('last_accepted.json',dict(response=current,requested=trialq,physical_gate=entry['physical_gate']))
                    if assessment['reduction_ratio']>.75:radius=min(.1,radius*2)
                    elif assessment['reduction_ratio']<.25:radius*=.5
                    break
                radius*=.5
            if not accepted:report['stop_reason']='No accepted bounded model step';break
        process.stdin.write('observe\n');process.stdin.flush();live=read()
        if any(initial[k]!=live[k] for k in initial if k!='kind'):raise ValueError('Continuing state changed')
        final_gate=physical_gate(current,weight,height)
        report.update(status='bounded_force_root_complete',continuing_state_unchanged=True,physical_gate=final_gate,accepted_equilibrium=final_gate['passed'],iterations=iterations,
            interpretation='Numerical force-residual merit only; source-bound KKT stationarity is never physical equilibrium; forward/reference promotion remains separate')
    except Exception as error:report.update(status='diagnostic_failed',error=type(error).__name__+': '+str(error)[:65536])
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if process is not None:
            try:
                if process.poll() is None:process.stdin.write('close\n');process.stdin.flush();process.wait(timeout=2)
            except Exception:
                if process.poll() is None:os.killpg(process.pid,signal.SIGKILL)
                process.wait(timeout=2)
        report.update(evaluations=len(records),wall_s=time.monotonic()-started,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'scripts/effective_passive_energy_observation.py',ROOT/'scripts/native_force_trust_region.py',ROOT/'scripts/prepare_native_force_scaling.py')})
        write('report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='iterations'}|{'output':str(output)},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');parser.add_argument('build',type=Path);args=parser.parse_args()
    if not args.run_native:raise SystemExit('Coordinated native grant required')
    run(args.build.resolve())
