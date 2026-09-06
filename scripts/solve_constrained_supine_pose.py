"""Gated static generalized-force solve with mandatory native support constraints."""
from pathlib import Path
import argparse,hashlib,json,signal,sys,tempfile,time
import numpy as np
from scipy.optimize import minimize,NonlinearConstraint,Bounds
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from build_supine_initial_state import analyze,recipe


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


def run(seed_path,material,resume_path=None):
    from ihm.native.mechanical_stream import NativeMechanicalStream
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
        if max(abs(v) for v in resumed['support_constraints'])>1e-4 or max(abs(v) for v in resumed['gauge_residual'])>1e-4:raise ValueError('Resume candidate lost support constraints')
        seed.update(resumed['coordinates'])
        seed['mtp_angle_r']=0.;seed['mtp_angle_l']=0.
    bounds=[base['coordinate_bounds'][n] for n in names]
    x0=np.array([seed[n] for n in names]);length=float(np.ptp(np.load(ROOT/manifest['arrays_path'])['reference_points_source_m'][:,1]))
    if any(not lo<=value<=hi for value,(lo,hi) in zip(x0,bounds)):raise ValueError('Supported seed outside source bounds')
    output=Path(tempfile.mkdtemp(prefix='constrained-supine-',dir=ROOT/'data/derived'));started=time.monotonic();stream=None;evaluations=[];cache={};best=None;best_feasible=None
    report=dict(passed=False,accepted_equilibrium=False,physical_time_advanced_s=0,material=material,
                maximum_evaluations=200,maximum_wall_s=60,objective='udot^T M udot/(mass*g^2), cross-checked against -r dot udot; no cost floor',
                toe_seed='held passive law neutral0; pure ankle damping has no preferred static angle so retained ankle q',held_gauge_coordinates={n:seed[n] for n in gauges},
                scope='Native generalized-force static solve with explicit support force/pitch/roll balance; unchanged forward acceptance remains mandatory')
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,stack):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Constrained statics reached60s wall cap')
    def evaluate(values):
        nonlocal best,best_feasible
        key=np.asarray(values,dtype=float).tobytes()
        if key in cache:return cache[key]
        if len(evaluations)>=200:raise RuntimeError('Constrained statics reached200 actual evaluations')
        candidate={**seed,**dict(zip(names,map(float,values)))}
        command=['evaluate_static_pose',str(len(all_names))]
        for name in all_names:command.extend([name,str(candidate[name])])
        native=stream._request(' '.join(command))
        weight=native['mass_kg']*np.linalg.norm(native['gravity_m_s2'])
        scales=np.where(native['mobility_rotational'],weight*length,weight)
        residual=np.asarray(native['constrained_zero_acceleration_residual_mobility_force'])/scales
        mapping=dict(zip(native['mobility_coordinate_names'],range(len(residual))))
        support=residual[[mapping[n] for n in ('pelvis_tx','pelvis_tilt','pelvis_rotation')]]
        root=residual[[mapping[n] for n in gauges]]
        cost=dynamic_metric(native);entry=dict(native=native,coordinates=candidate,cost=cost,support_constraints=support.tolist(),gauge_residual=root.tolist())
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
        result=minimize(lambda q:evaluate(q)['cost'],x0,method='trust-constr',jac=lambda q:derivatives(q)[0],
             bounds=Bounds(*np.array(bounds).T,keep_feasible=True),constraints=[constraint],
             options={'maxiter':30,'gtol':1e-8,'xtol':1e-10,'initial_tr_radius':.03,'verbose':0})
        final=evaluate(result.x);write('final_candidate.json',final)
        report.update(optimizer_success=bool(result.success),optimizer_message=str(result.message),
             static_candidate_converged=bool(np.max(np.abs(final['native']['udot']))<=1e-4 and np.max(np.abs(final['support_constraints']))<=1e-4))
        report['status']='static_candidate_converged' if report['static_candidate_converged'] else 'static_residual_unresolved'
        observed=stream._request('observe');report['continuing_state_unchanged']=all(before[k]==observed[k] for k in before if k!='kind')
        if not report['continuing_state_unchanged']:raise AssertionError('Static solve mutated continuing state')
        stream.restore(checkpoint);stream.release(checkpoint)
    except Exception as error:report.update(status='stopped',error=f'{type(error).__name__}: {error}')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        report.update(wall_s=time.monotonic()-started,evaluations=len(evaluations),best_cost=None if best is None else best['cost'],
             best_supported_cost=None if best_feasible is None else best_feasible['cost'],
             source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (seed_path,manifest_path,Path(__file__).resolve(),*(() if resume_path is None else (resume_path,)))})
        write('report.json',report);print(json.dumps({**report,'output':str(output)},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-native',action='store_true');parser.add_argument('--seed',type=Path);parser.add_argument('--resume',type=Path);parser.add_argument('--material',choices=('MM','HM'),default='MM');args=parser.parse_args()
    if not args.run_native or args.seed is None:raise SystemExit('Coordinated --run-native slot and --seed required')
    run(args.seed.resolve(),args.material,None if args.resume is None else args.resume.resolve())
