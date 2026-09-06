"""Gated actual-force trust region on separately attested98-muscle physics."""
from pathlib import Path
import argparse, json, os, signal, subprocess, tempfile, time
import xml.etree.ElementTree as ET
import numpy as np
from build_effective_potential_probe import verify, ROOT, sha
from effective_passive_energy_observation import PassiveEnergyObservation

from prepare_native_force_scaling import model_info, force, GAUGES
from support_feasible_physical_step import native_metric, physical_step, feasible_trial


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


def run(build_path,resume_path=None,initial_radius=None):
    if initial_radius is not None and not 0<initial_radius<=.03:raise ValueError("Invalid explicit native trust radius")
    build=verify(build_path);exe=ROOT/build['executable']
    prerequisite=ROOT/'data/derived/native-physical-metric-04jdqqbn/report.json'
    prerequisite_record=json.loads(prerequisite.read_text())
    if not prerequisite_record['metric_gate_passed'] or prerequisite_record['build_manifest_sha256']!=sha(build_path):raise ValueError('Matching native physical-metric gate required')
    if sha(exe)!=build.get('executable_sha256'):raise ValueError('Unverified executable')
    folder=build_path.parent;archive=ROOT/build['archive_manifest'];ar=json.loads(archive.read_text())
    case=json.loads((ROOT/build['reference_manifest']).read_text())
    bounds={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.parse(folder/'inputs/subject_walk_scaled.osim').getroot().iter('Coordinate')}
    port=PassiveEnergyObservation(folder/'inputs/subject_walk_scaled.osim')
    output=Path(tempfile.mkdtemp(prefix='support-physical-root-',dir=ROOT/'data/derived'))
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
            native_metric(value)
            return value
        q=json.loads((ROOT/case['old_q_seed']).read_text())['coordinates'];resume=None;resume_report=None
        if resume_path is not None:
            resume=json.loads(resume_path.read_text());resume_report=json.loads((resume_path.parent/'report.json').read_text())
            if resume_report['build_manifest_sha256']!=sha(build_path) or not resume_report['continuing_state_unchanged']:raise ValueError('Resume physics identity or continuing-state receipt differs')
            if max(map(abs,resume['physical_gate']['support_constraints']+resume['physical_gate']['gauge_residual']))>1e-4:raise ValueError('Unsupported resume candidate')
            q=dict(resume['requested'])
            report['resume']=dict(candidate=str(resume_path.relative_to(ROOT)),candidate_sha256=sha(resume_path),report_sha256=sha(resume_path.parent/'report.json'),observations_sha256=sha(resume_path.parent/'observations.jsonl'),responses_reused=0,basis='Exact requested coordinates; fresh native response and Jacobian, prior records immutable')
        current=evaluate(q,'base');names=current['independent_names'];free=[n for n in names if n not in GAUGES]
        indices=[names.index(n) for n in free];B=np.array([bounds[n] for n in free]);_,couplings=model_info(folder/'inputs/subject_walk_scaled.osim')
        surface=json.loads((ROOT/case['surface_manifest']).read_text())
        height=float(np.ptp(np.load(ROOT/surface['arrays_path'])['reference_points_source_m'][:,1]))
        gravity=current['gravity_generalized_forces'];weight=float(np.linalg.norm([gravity[current['mobility_names'].index(n)] for n in ('pelvis_tx','pelvis_ty','pelvis_tz')]))
        radius=.03;iterations=[]
        if resume is not None:
            previous=next(v for v in reversed(resume_report['iterations']) if v.get('assessment',{}).get('accepted'))
            radius=previous['diagnostic']['radius'];ratio=previous['assessment']['reduction_ratio']
            if ratio>.75:radius=min(.03,radius*2)
            elif ratio<.25:radius*=.5
            replay_gate=physical_gate(current,weight,height)
            if max(map(abs,replay_gate['support_constraints']+replay_gate['gauge_residual']))>1e-4:raise ValueError('Fresh resume response lost physical support')
            report['resume'].update(radius=radius,maximum_actual_q_difference=float(np.max(np.abs(np.array(current['independent_q'])-resume['response']['independent_q']))),maximum_acceleration_difference=float(np.max(np.abs(np.array(current['udot'])-resume['response']['udot']))))
        if initial_radius is not None:
            report['explicit_initial_radius']=dict(value=initial_radius,basis='Explicit diagnostic radius after derivative formulation change; physical bounds unchanged');radius=initial_radius
        rootnames=('pelvis_tx','pelvis_tilt','pelvis_rotation');rootindices=[current['mobility_names'].index(n) for n in rootnames];rootcolumns=[free.index(n) for n in rootnames]
        rootscale=np.array([weight,weight*height,weight*height])
        report['native_metric_prerequisite']=dict(path=str(prerequisite.relative_to(ROOT)),sha256=sha(prerequisite))
        write('initial_gate.json',physical_gate(current,weight,height))
        for iteration in range(6):
            gate=physical_gate(current,weight,height)
            if gate['passed']:report['accepted_equilibrium']=True;break
            q=dict(zip(names,current['independent_q']));x=np.array([q[n] for n in free]);r=np.asarray(current['constrained_residual']);acceleration=native_metric(current)
            D=[];Y=[];Ya=[];maximum_gauge_drift=0.
            for i,n in enumerate(free):
                shifted=dict(q);h=1e-5 if q[n]+1e-5<=bounds[n][1] else -1e-5
                shifted[n]+=h
                if shifted[n]<bounds[n][0]:raise ValueError('No bounded Jacobian direction')
                sample=evaluate(shifted,f'jacobian:{iteration}:{n}')
                D.append(np.array(sample['independent_q'])[indices]-x);Y.append(np.asarray(sample['constrained_residual'])-r);Ya.append(np.asarray(sample['udot'])-current['udot'])
                maximum_gauge_drift=max(maximum_gauge_drift,max(abs(sample['independent_q'][names.index(g)]-q[g]) for g in GAUGES))
            condition=float(np.linalg.cond(D))
            if condition>100 or maximum_gauge_drift>1e-8:raise ValueError('Unresolved assembled Jacobian coordinate map')
            J=np.linalg.solve(D,Y).T
            frozen_mass_jacobian=-np.asarray(current['inverse_mass_matrix'])@J;full_acceleration_jacobian=np.linalg.solve(D,Ya).T
            Bmetric=full_acceleration_jacobian
            A=J[rootindices]/rootscale[:,None];rootblock=A[:,rootcolumns]
            if np.linalg.cond(rootblock)>1e10:raise ValueError('Ill-conditioned native support correction')
            write('last_jacobian.json',dict(iteration=iteration,q=q,force=r.tolist(),jacobian=J.tolist(),physical_jacobian=Bmetric.tolist(),derivative_basis='Full observed native acceleration Jacobian; includes configuration dependence of MInv and native constraint reactions',condition=condition,maximum_gauge_drift=maximum_gauge_drift,
                frozen_mass_vs_full_acceleration_jacobian_relative_difference=float(np.linalg.norm(frozen_mass_jacobian-full_acceleration_jacobian)/np.linalg.norm(full_acceleration_jacobian))))
            accepted=False
            for attempt in range(8):
                step,diagnostic=physical_step(acceleration,Bmetric,x,B,gate['support_constraints'],A,radius)
                trialq=dict(q);trialq.update(zip(free,(x+step).tolist()))
                try:
                    for correction in range(4):
                        trial=evaluate(trialq,f'trial:{iteration}:{attempt}:correction:{correction}')
                        trialgate=physical_gate(trial,weight,height)
                        if max(map(abs,trialgate['support_constraints']+trialgate['gauge_residual']))<=1e-4:break
                        if correction==3:break
                        change=np.linalg.solve(rootblock,-np.asarray(trialgate['support_constraints']))
                        corrected=np.array(trial['independent_q'])[indices]
                        corrected[rootcolumns]+=change
                        corrected=np.clip(corrected,np.maximum(B[:,0],x-radius),np.minimum(B[:,1],x+radius))
                        trialq=dict(q);trialq.update(zip(free,corrected.tolist()))
                except ValueError as error:
                    if str(error)!='equal-pressure skin/bed solution exceeds retained domains':raise
                    iterations.append(dict(iteration=iteration,attempt=attempt,radius=radius,accepted=False,domain_rejection=str(error)));radius*=.5;continue
                actual_step=np.array(trial['independent_q'])[indices]-x
                prediction=.5*float(acceleration@acceleration-np.sum((acceleration+Bmetric@actual_step)**2))
                assessment=feasible_trial(current,dict(udot=trial['udot'],**trialgate),prediction)
                entry=dict(iteration=iteration,attempt=attempt,support_corrections=correction,diagnostic=diagnostic,assessment=assessment,physical_gate=trialgate)
                iterations.append(entry);write('iterations.json',iterations)
                if assessment['accepted']:
                    current=trial;accepted=True
                    write('last_accepted.json',dict(response=current,requested=trialq,physical_gate=trialgate))
                    if assessment['reduction_ratio']>.75:radius=min(.03,radius*2)
                    elif assessment['reduction_ratio']<.25:radius*=.5
                    break
                radius*=.5
            if not accepted:report['stop_reason']='No accepted bounded model step';break
        process.stdin.write('observe\n');process.stdin.flush();live=read()
        if any(initial[k]!=live[k] for k in initial if k!='kind'):raise ValueError('Continuing state changed')
        final_gate=physical_gate(current,weight,height)
        report.update(status='bounded_force_root_complete',continuing_state_unchanged=True,physical_gate=final_gate,accepted_equilibrium=final_gate['passed'],iterations=iterations,
            interpretation='Native inverse-mass acceleration metric with hard nonlinear support/gauge acceptance; bound KKT reactions never count as equilibrium; no forward/reference promotion')
    except Exception as error:report.update(status='diagnostic_failed',error=type(error).__name__+': '+str(error)[:65536])
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if process is not None:
            try:
                if process.poll() is None:process.stdin.write('close\n');process.stdin.flush();process.wait(timeout=2)
            except Exception:
                if process.poll() is None:os.killpg(process.pid,signal.SIGKILL)
                process.wait(timeout=2)
        report.update(evaluations=len(records),wall_s=time.monotonic()-started,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'scripts/effective_passive_energy_observation.py',ROOT/'scripts/support_feasible_physical_step.py',ROOT/'scripts/prepare_native_force_scaling.py')})
        write('report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='iterations'}|{'output':str(output)},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-native',action='store_true');parser.add_argument('build',type=Path);parser.add_argument('--resume',type=Path);parser.add_argument('--initial-radius',type=float);args=parser.parse_args()
    if not args.run_native:raise SystemExit('Coordinated native grant required')
    run(args.build.resolve(),None if args.resume is None else args.resume.resolve(),args.initial_radius)
