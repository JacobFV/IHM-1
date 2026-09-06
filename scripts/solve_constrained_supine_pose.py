"""Gated static generalized-force solve with mandatory native support constraints."""
from pathlib import Path
import argparse,hashlib,json,signal,sys,tempfile,time
import numpy as np
from scipy.optimize import minimize,NonlinearConstraint,Bounds,least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from build_supine_initial_state import analyze,recipe
from static_pose_journal import PoseJournal,pose_key
from bounded_static_root import interior_origin,local_step,StaticDomainRejection,backtracked_trial,constrained_local_step,balanced_backtracked_trial,recomputed_box_trial


def dynamic_metric(native):
    residual=np.asarray(native['constrained_zero_acceleration_residual_mobility_force'],float)
    acceleration=np.asarray(native['udot'],float);mass_acceleration=np.asarray(native['mass_times_udot'],float)
    identity=float(np.linalg.norm(residual+mass_acceleration));scale=max(1.,float(np.linalg.norm(residual)),float(np.linalg.norm(mass_acceleration)))
    if identity>1e-7*scale:raise ValueError('Native r=-M udot identity failed')
    value=float(acceleration@mass_acceleration)
    if not np.isfinite(value) or value<0:raise ValueError('Nonfinite or negative dynamic mass metric')
    dual=float(-residual@acceleration)
    if abs(value-dual)>1e-7*max(1.,abs(value)):raise ValueError('Primal/dual dynamic metric mismatch')
    normalization=native['mass_kg']*np.linalg.norm(native['gravity_m_s2'])**2
    if not np.isfinite(normalization) or normalization<=0:raise ValueError('Invalid dynamic metric normalization')
    return value/normalization


def acceleration_residual(native):
    """All mobilities, each divided by 1 rad/s² or 1 m/s²; no inertia weights."""
    acceleration=np.asarray(native['udot'],float)
    if acceleration.ndim!=1 or len(acceleration)!=len(native['mobility_rotational']) or not np.all(np.isfinite(acceleration)):
        raise ValueError('Invalid all-mobility acceleration residual')
    return acceleration


def static_converged(entry):
    native=entry['native']
    return bool(np.max(np.abs(acceleration_residual(native)))<=1e-4
        and max(map(abs,entry['support_constraints']+entry['gauge_residual']))<=1e-4
        and all(np.isfinite(native[key]) and abs(native[key])<=1e-5 for key in
            ('constraint_position_error','constraint_velocity_error','constraint_acceleration_error')))


def run(seed_path,material,resume_path=None,resume_cache=None,mode='constrained',frozen_stream=None,equivalence_receipt=None,recompute_box=False,local_radius=.03):
    if not 0<local_radius<=.03:raise ValueError('Invalid initial local radius')
    if recompute_box and mode!='balanced-root':raise ValueError('Recomputed box requires balanced-root mode')
    if frozen_stream is None:
        from ihm.native.mechanical_stream import NativeMechanicalStream
    else:
        from frozen_static_stream import load
        if not frozen_stream.resolve().is_relative_to(ROOT):raise ValueError('Owned frozen stream manifest required')
        NativeMechanicalStream=load(frozen_stream)
    seed_record=json.loads(seed_path.read_text())
    if not seed_record['passed'] or seed_record['material']!=material:raise ValueError('Matching force/moment-supported rigid seed required')
    manifest_path=ROOT/'data/derived/supine-surface-contact-exmzq9pq/manifest.json';manifest=json.loads(manifest_path.read_text())
    initial_path=ROOT/'data/derived/supine-support-5ma720yd/initial_native.json';initial=json.loads(initial_path.read_text())
    model_path=ROOT/'data/derived/supine-support-5ma720yd/plant/native/inputs/subject_walk_scaled.osim'
    analysis=analyze(initial,model_path.read_bytes());base=recipe(analysis,initial,{})
    # Remove only plane-translation/heading gauges from optimization, without
    # locking native joints. Their full native force residuals are still checked.
    gauges=('pelvis_ty','pelvis_tz','pelvis_list')
    all_names=list(base['seed_coordinates']);names=[n for n in all_names if n not in gauges]
    seed=dict(seed_record['best']['seed_coordinates'])
    if resume_path is not None:
        resumed=json.loads(resume_path.read_text())
        if mode in ('constrained','balanced-root') and (max(abs(v) for v in resumed['support_constraints'])>1e-4 or max(abs(v) for v in resumed['gauge_residual'])>1e-4):raise ValueError('Resume candidate lost support constraints')
        seed.update(resumed['coordinates'])
        if resume_cache is None:seed['mtp_angle_r']=0.;seed['mtp_angle_l']=0.
    bounds=[base['coordinate_bounds'][n] for n in names]
    x0=np.array([seed[n] for n in names]);length=float(np.ptp(np.load(ROOT/manifest['arrays_path'])['reference_points_source_m'][:,1]))
    if any(not lo<=value<=hi for value,(lo,hi) in zip(x0,bounds)):raise ValueError('Supported seed outside source bounds')
    output=Path(tempfile.mkdtemp(prefix='constrained-supine-',dir=ROOT/'data/derived'));started=time.monotonic();stream=None;evaluations=[];cache={};best=None;best_feasible=None;journal=None;cache_hits=0
    report=dict(passed=False,accepted_equilibrium=False,physical_time_advanced_s=0,material=material,
                maximum_evaluations=200,maximum_wall_s=60,mode=mode,objective=('all mobility accelerations divided by 1 rad/s^2 or 1 m/s^2, no inertia weighting' if mode in ('acceleration-root','balanced-root') else 'udot^T M udot/(mass*g^2), cross-checked against -r dot udot; no cost floor'),
                toe_seed=('exact resumed coordinates retained for cache reuse' if resume_cache is not None else 'held passive law neutral0 on resume; pure ankle damping has no preferred static angle so retained ankle q'),held_gauge_coordinates={n:seed[n] for n in gauges},
                scope='Native generalized-force static solve with explicit support force/pitch/roll balance; unchanged forward acceptance remains mandatory')
    if frozen_stream is not None:
        if equivalence_receipt is None:raise ValueError('Frozen continuation requires numerical-equivalence receipt')
        equivalence=json.loads(equivalence_receipt.read_text())
        if not equivalence.get('passed') or not equivalence.get('archive_intact') or not equivalence.get('continuing_state_unchanged') or equivalence.get('assembly_accuracy')!=1e-10:
            raise ValueError('Unaccepted frozen numerical-equivalence receipt')
        if equivalence['source_sha256'].get(str(frozen_stream.resolve()))!=hashlib.sha256(frozen_stream.read_bytes()).hexdigest():raise ValueError('Frozen equivalence archive differs')
        report['frozen_numerical_equivalence']=dict(receipt=str(equivalence_receipt),sha256=hashlib.sha256(equivalence_receipt.read_bytes()).hexdigest(),
            scope=equivalence['comparison_basis'],limitation=equivalence['limitation'],strict_generic_failure_retained='data/derived/frozen-static-acceptance-57u8erpr/report.json')
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,stack):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Constrained statics reached60s wall cap')
    def evaluate(values):
        nonlocal best,best_feasible,cache_hits
        key=pose_key(values)
        if key in cache:
            cache_hits+=1
            return cache[key]
        if len(evaluations)>=200:raise RuntimeError('Constrained statics reached200 actual evaluations')
        candidate={**seed,**dict(zip(names,map(float,values)))}
        command=['evaluate_static_pose',str(len(all_names))]
        for name in all_names:command.extend([name,str(candidate[name])])
        write('last_native_request.json',dict(status='pending',new_evaluation=len(evaluations)+1,
            coordinates=candidate,command=' '.join(command),wall_s=time.monotonic()-started))
        try:native=stream._request(' '.join(command))
        except ValueError as error:
            if type(error) is not ValueError or str(error)!='equal-pressure skin/bed solution exceeds retained domains':raise
            rejection=dict(evaluation=len(evaluations)+1,wall_s=time.monotonic()-started,coordinates=candidate,
                command=' '.join(command),classification='copied-state-material-domain-rejection',error=str(error))
            with (output/'rejected_trials.jsonl').open('a') as destination:destination.write(json.dumps(rejection,allow_nan=False)+'\n')
            evaluations.append(dict(evaluation=len(evaluations)+1,wall_s=time.monotonic()-started,rejected=True,error=str(error)))
            write('evaluations.json',evaluations)
            raise StaticDomainRejection(str(error)) from error
        write('last_native_request.json',dict(status='completed',new_evaluation=len(evaluations)+1,
            coordinates=candidate,wall_s=time.monotonic()-started))
        weight=native['mass_kg']*np.linalg.norm(native['gravity_m_s2'])
        scales=np.where(native['mobility_rotational'],weight*length,weight)
        residual=np.asarray(native['constrained_zero_acceleration_residual_mobility_force'])/scales
        mapping=dict(zip(native['mobility_coordinate_names'],range(len(residual))))
        support=residual[[mapping[n] for n in ('pelvis_tx','pelvis_tilt','pelvis_rotation')]]
        root=residual[[mapping[n] for n in gauges]]
        mass_metric=dynamic_metric(native)
        cost=mass_metric if mode=='constrained' else float(acceleration_residual(native)@acceleration_residual(native))
        entry=dict(native=native,coordinates=candidate,cost=cost,dynamic_mass_metric=mass_metric,support_constraints=support.tolist(),gauge_residual=root.tolist())
        journal.append(values,entry)
        cache[key]=entry;evaluations.append(dict(evaluation=len(evaluations)+1,wall_s=time.monotonic()-started,cost=cost,
            support_constraint_norm=float(np.linalg.norm(support)),full_root_residual_norm=float(np.linalg.norm(np.r_[support,root])),
            maximum_abs_udot=float(np.max(np.abs(native['udot']))),maximum_skin_indentation_m=native['maximum_penetration_m']))
        if best is None or cost<best['cost']:best=entry;write('best_candidate.json',best)
        if max(np.max(np.abs(support)),np.max(np.abs(root)))<=1e-4 and (best_feasible is None or cost<best_feasible['cost']):
            best_feasible=entry;write('best_supported_candidate.json',best_feasible)
        write('evaluations.json',evaluations)
        return entry
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,60)
    try:
        stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=77.6122029,
            augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',surface_contact_manifest=manifest_path,bed_material=material)
        execution=json.loads((output/'native/execution.json').read_text())
        report['native_execution_basis']='live worktree validated stream' if frozen_stream is None else 'immutable archived compiled physics; current worktree native source is not used'
        if frozen_stream is not None:
            report['frozen_stream_attestation']=execution['frozen_archive_attestation']
            report['frozen_loader_sha256']=hashlib.sha256((ROOT/'scripts/frozen_static_stream.py').read_bytes()).hexdigest()
        identity=dict(schema='ihm.static-pose-cache.v1',protocol_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            journal_sha256=hashlib.sha256((ROOT/'scripts/static_pose_journal.py').read_bytes()).hexdigest(),
            solver_sha256=hashlib.sha256((ROOT/'scripts/bounded_static_root.py').read_bytes()).hexdigest(),
            source_sha256=execution['source_sha256'],build_files=execution['build']['files'],
            material=material,mass_kg=77.6122029,environment='supine',mode=('acceleration-root' if mode=='balanced-root' else mode),coordinate_order=names,
            held_gauges={n:seed[n] for n in gauges},bounds=bounds,
            surface_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest())
        # JSON normalization keeps tuples/lists identical after disk round-trip.
        identity=json.loads(json.dumps(identity))
        journal=PoseJournal(output,identity,resume_cache,reuse_native_only=resume_cache is not None,framing_root=ROOT);cache=journal.cache
        report['recovered_evaluations']=len(cache)
        for entry in cache.values():
            dynamic_metric(entry['native'])
            if best is None or entry['cost']<best['cost']:best=entry
            if max(map(abs,entry['support_constraints']+entry['gauge_residual']))<=1e-4 and (best_feasible is None or entry['cost']<best_feasible['cost']):best_feasible=entry
        if best is not None:write('best_candidate.json',best)
        if best_feasible is not None:write('best_supported_candidate.json',best_feasible)
        checkpoint=stream.checkpoint();before=stream.snapshot()
        first=evaluate(x0);write('initial_supported_seed_native.json',first)
        def derivatives(q):
            base=evaluate(q);gradient=np.zeros(len(q));jacobian=np.zeros((3,len(q)))
            for index,(lo,hi) in enumerate(bounds):
                step=1e-5 if q[index]+1e-5<=hi else -1e-5
                shifted=np.array(q,float);shifted[index]+=step;value=evaluate(shifted)
                gradient[index]=(value['cost']-base['cost'])/step
                jacobian[:,index]=(np.array(value['support_constraints'])-base['support_constraints'])/step
            return gradient,jacobian
        constraint=NonlinearConstraint(lambda q:np.array(evaluate(q)['support_constraints']),0.,0.,jac=lambda q:derivatives(q)[1])
        def retain_iterate(q,state=None):
            entry=evaluate(q)
            record=dict(iteration=None if state is None else int(state.nit),coordinates=entry['coordinates'],
                values=list(map(float,q)),cost=entry['cost'],support_constraints=entry['support_constraints'],
                gauge_residual=entry['gauge_residual'],trust_radius=None if state is None else float(state.tr_radius))
            with (output/'optimizer_iterates.jsonl').open('a') as destination:destination.write(json.dumps(record,allow_nan=False)+'\n')
            write('last_optimizer_iterate.json',entry)
            return False
        if mode in ('acceleration-root','balanced-root'):
            q=interior_origin(x0,bounds,preserve=mode=='balanced-root')
            report['local_solver']=dict(method=('support-constrained linear Newton plus bounded nonlinear support correction' if mode=='balanced-root' else 'bounded variable least-squares Newton step plus backtracking'),
                maximum_support_corrections=(3 if mode=='balanced-root' else 0),
                maximum_coordinate_step=local_radius,recompute_box=recompute_box,interior_margin=(0. if mode=='balanced-root' else 1e-12),maximum_iterations=10,
                maximum_backtracks=6,origin_coordinate_changes={name:float(v-x0[i]) for i,(name,v) in enumerate(zip(names,q)) if v!=x0[i]})
            success=False;message='Local Newton iteration cap';active_radius=local_radius
            for iteration in range(10):
                current=evaluate(q)
                if static_converged(current):success=True;message='All static checks passed';break
                acceleration=acceleration_residual(current['native']);jac=np.zeros((len(acceleration),len(q)));support_jac=np.zeros((3,len(q)))
                for index,(lo,hi) in enumerate(bounds):
                    step=1e-5 if q[index]+1e-5<=hi else -1e-5
                    shifted=q.copy();shifted[index]+=step
                    shifted_entry=evaluate(shifted)
                    jac[:,index]=(acceleration_residual(shifted_entry['native'])-acceleration)/step
                    support_jac[:,index]=(np.asarray(shifted_entry['support_constraints'])-np.asarray(current['support_constraints']))/step
                if mode=='balanced-root':
                    root_indices=[names.index(name) for name in ('pelvis_tx','pelvis_tilt','pelvis_rotation')]
                    if recompute_box:
                        candidate,trial,active_radius,attempts=recomputed_box_trial(jac,acceleration,q,np.array(bounds),support_jac,
                            current['support_constraints'],current['cost'],evaluate,root_indices,active_radius)
                        diagnostic=dict(recomputed_attempts=attempts,retained_radius=active_radius)
                    else:
                        direction,diagnostic=constrained_local_step(jac,acceleration,q,bounds,support_jac,current['support_constraints'],radius=active_radius)
                        candidate,trial=balanced_backtracked_trial(q,direction,np.array(bounds),current['cost'],evaluate,support_jac,root_indices,radius=active_radius)
                    with (output/'linear_support_steps.jsonl').open('a') as destination:destination.write(json.dumps(dict(iteration=iteration,**diagnostic))+'\n')
                else:
                    direction=local_step(jac,acceleration,q,bounds)
                    candidate,trial=backtracked_trial(q,direction,np.array(bounds),current['cost'],evaluate)
                if candidate is None:message='No decreasing valid local Newton step';break
                q=candidate;retain_iterate(q)
            final=evaluate(q)
            from types import SimpleNamespace
            result=SimpleNamespace(success=success,message=message)
        else:
            result=minimize(lambda q:evaluate(q)['cost'],x0,method='trust-constr',jac=lambda q:derivatives(q)[0],
                 bounds=Bounds(*np.array(bounds).T,keep_feasible=True),constraints=[constraint],callback=retain_iterate,
                 options={'maxiter':30,'gtol':1e-8,'xtol':1e-10,'initial_tr_radius':.03,'verbose':0})
            final=evaluate(result.x)
        write('final_candidate.json',final)
        report.update(optimizer_success=bool(result.success),optimizer_message=str(result.message),
             static_candidate_converged=static_converged(final))
        report['status']='static_candidate_converged' if report['static_candidate_converged'] else 'static_residual_unresolved'
        observed=stream._request('observe');report['continuing_state_unchanged']=all(before[k]==observed[k] for k in before if k!='kind')
        if not report['continuing_state_unchanged']:raise AssertionError('Static solve mutated continuing state')
        stream.restore(checkpoint);stream.release(checkpoint)
    except Exception as error:
        report.update(status='stopped',error=f'{type(error).__name__}: {error}')
        failure=dict(exception_type=type(error).__name__,physical_or_protocol_failure_text=str(error)[:65536])
        if isinstance(error,json.JSONDecodeError):
            failure.update(response_payload=error.doc[:65536],payload_truncated=len(error.doc)>65536,parse_position=error.pos)
        write('failure.json',failure)
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        if frozen_stream is not None:
            try:
                from frozen_static_stream import validate
                validate(frozen_stream);report['frozen_archive_intact_after_run']=True
            except Exception as error:
                report.update(status='frozen_archive_identity_failed',frozen_archive_intact_after_run=False,error=str(error))
        report.update(wall_s=time.monotonic()-started,evaluations=len(evaluations),cache_hits=cache_hits,best_cost=None if best is None else best['cost'],
             best_supported_cost=None if best_feasible is None else best_feasible['cost'],
             source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (seed_path,manifest_path,Path(__file__).resolve(),ROOT/'scripts/bounded_static_root.py',*(() if resume_path is None else (resume_path,)))})
        write('report.json',report);print(json.dumps({**report,'output':str(output)},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-native',action='store_true');parser.add_argument('--seed',type=Path);parser.add_argument('--resume',type=Path);parser.add_argument('--resume-cache',type=Path);parser.add_argument('--frozen-stream',type=Path);parser.add_argument('--equivalence-receipt',type=Path);parser.add_argument('--mode',choices=('constrained','acceleration-root','balanced-root'),default='constrained');parser.add_argument('--material',choices=('MM','HM'),default='MM');parser.add_argument('--recompute-box',action='store_true');parser.add_argument('--local-radius',type=float,default=.03);args=parser.parse_args()
    if not args.run_native or args.seed is None:raise SystemExit('Coordinated --run-native slot and --seed required')
    run(args.seed.resolve(),args.material,None if args.resume is None else args.resume.resolve(),None if args.resume_cache is None else args.resume_cache.resolve(),args.mode,None if args.frozen_stream is None else args.frozen_stream.resolve(),None if args.equivalence_receipt is None else args.equivalence_receipt.resolve(),args.recompute_box,args.local_radius)
