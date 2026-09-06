"""Gated all-muscle/component static-gradient audit; no optimizer or integration."""
from pathlib import Path
import argparse,json,os,selectors,signal,subprocess,tempfile,time,hashlib
import xml.etree.ElementTree as ET
import numpy as np
from build_effective_potential_probe import verify,ROOT,sha


def pullback(names,independent,couplings,values):
    result=np.array([values[names.index(n)] for n in independent],float)
    for dependent,(parent,coefficient) in couplings.items():result[independent.index(parent)]+=coefficient*values[names.index(dependent)]
    return result


def run(build_path):
    build=verify(build_path);exe=ROOT/build['executable']
    if sha(exe)!=build.get('executable_sha256'):raise ValueError('Unverified probe executable')
    case=json.loads((ROOT/build['reference_manifest']).read_text());archive=ROOT/build['archive_manifest'];ar=json.loads(archive.read_text());folder=build_path.parent
    model=ET.parse(folder/'inputs/subject_walk_scaled.osim').getroot();coordinates={c.get('name'):list(map(float,c.findtext('range').split())) for c in model.iter('Coordinate')};couplings={}
    for c in model.findall('.//ConstraintSet/objects/*'):
        coeff=list(map(float,c.findtext('.//LinearFunction/coefficients').split()))
        if c.tag!='CoordinateCouplerConstraint' or len(coeff)!=2 or coeff[1]!=0:raise ValueError('Unsupported actual kinematic constraint')
        couplings[c.findtext('dependent_coordinate_name')]=(c.findtext('independent_coordinate_names'),coeff[0])
    output=Path(tempfile.mkdtemp(prefix='native-effective-potential-',dir=ROOT/'data/derived'));process=None;started=time.monotonic();records=[]
    report=dict(build_manifest=str(build_path.relative_to(ROOT)),build_manifest_sha256=sha(build_path),maximum_evaluations=65,maximum_wall_s=45,physical_time_advanced_s=0,accepted_equilibrium=False,numerical_static_merit_only=True)
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,frame):
        if process is not None and process.poll() is None:process.kill()
        raise TimeoutError('Effective-potential native45s cap')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,45)
    try:
        command=['prlimit','--as=4294967296','--','nice','-n','10',str(archive.parent/ar['dynamic_loader']),'--library-path',str(archive.parent/'libraries'),str(exe),str(folder/'inputs'),str(output),'supine',str(case['mass_kg'])]
        log=(output/'engine.log').open('w');process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=log,text=True,bufsize=1,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LD_LIBRARY_PATH=str(archive.parent/'libraries')))
        def read():
            while True:
                line=process.stdout.readline()
                if not line:raise RuntimeError('Native effective probe terminated')
                if line.startswith('@IHM '):
                    value=json.loads(line[5:])
                    if 'error' in value:raise ValueError(value['error'])
                    return value
                log.write(line);log.flush()
        before=read();write('initial.json',before)
        def evaluate(q,label):
            if len(records)>=65:raise RuntimeError('Effective-potential65-call cap')
            request=['evaluate_effective_pose',str(len(q))]
            for n,v in q.items():request.extend([n,str(v)])
            write('pending_request.json',dict(label=label,q=q));process.stdin.write(' '.join(request)+'\n');process.stdin.flush();value=read()
            record=dict(label=label,requested=q,response=value);records.append(record)
            with (output/'observations.jsonl').open('a') as f:f.write(json.dumps(record,allow_nan=False)+'\n')
            return value
        q=json.loads((ROOT/case['old_q_seed']).read_text())['coordinates'];base=evaluate(q,'base');write('base.json',base);names=base['independent_names'];q=dict(zip(names,base['independent_q']));pairs=[]
        for n in names:
            h=1e-5 if q[n]+2e-5<=coordinates[n][1] else -1e-5
            if q[n]+2*h<coordinates[n][0]:raise ValueError('No bounded probe direction')
            pair=[]
            for factor in (1,2):
                shifted=dict(q);shifted[n]+=factor*h;pair.append(evaluate(shifted,n+'_'+str(factor)))
            pairs.append(pair)
        D=np.array([2*(np.array(a['independent_q'])-base['independent_q'])-.5*(np.array(b['independent_q'])-base['independent_q']) for a,b in pairs]);condition=float(np.linalg.cond(D))
        if condition>100:raise ValueError('Ill-conditioned actual assembled finite-difference map')
        def gradient(get):
            y=np.array([2*(get(a)-get(base))-.5*(get(b)-get(base)) for a,b in pairs]);return np.linalg.solve(D,y)
        g=gradient(lambda a:a['effective_energy_j']);physical=pullback(base['mobility_names'],names,couplings,base['tree_residual']);constrained=pullback(base['mobility_names'],names,couplings,base['constrained_residual'])
        muscles=[]
        for i,m in enumerate(base['muscles']):
            gl=gradient(lambda a:a['muscles'][i]['length_m']);ge=gradient(lambda a:a['muscles'][i]['passive_energy_j']+a['muscles'][i]['active_effective_energy_j']);ma=np.array(m['moment_arms_m'])
            muscles.append(dict(name=m['name'],type=m['type'],length_moment_arm_error_m=float(np.max(np.abs(gl+ma))),energy_virtual_work_error=float(np.max(np.abs(ge+m['tendon_force_n']*ma))),fiber_minimum_margin_m=m['fiber_length_m']-m['minimum_fiber_length_m'],fiber_velocity_m_s=m['fiber_velocity_m_s'],force_balance_error_n=m['tendon_force_n']-m['fiber_force_along_tendon_n'],native_initialization_tolerance_n=1e-8*m['fiso_n'],activation=m['activation'],internal_stiffness_n_m=m['internal_stiffness_n_m']))
        repeat=evaluate(q,'repeat_base');process.stdin.write('observe\n');process.stdin.flush();live=read();unchanged=all(before[k]==live[k] for k in before if k!='kind')
        if not unchanged:raise ValueError('Probe changed continuing state')
        component_errors={}
        for label,energy in [('gravity','gravity_energy_j'),('joint','joint_energy_j'),('contact',None)]:
            derivative=gradient((lambda a:a[energy]) if energy else (lambda a:a['skin_energy_j']+a['bed_energy_j']))
            expected=-pullback(base['mobility_names'],names,couplings,base[label+'_generalized_forces'])
            component_errors[label]=float(np.max(np.abs(derivative-expected)))
        all_muscles=[m for record in records for m in record['response']['muscles']]
        all_branches_valid=all(m['fiber_length_m']-m['minimum_fiber_length_m']>1e-10 and m['internal_stiffness_n_m']>0 and abs(m['fiber_velocity_m_s'])<=1e-5 for m in all_muscles)
        errors=(g-physical).tolist();report.update(status='diagnostic_complete',continuing_state_unchanged=True,actual_difference_matrix_condition=condition,
            coordinate_names=names,energy_gradient=g.tolist(),physical_gradient=physical.tolist(),gradient_errors=errors,
            maximum_gradient_error=max(map(abs,errors)),component_gradient_errors=component_errors,constraint_reaction_virtual_work_error=float(np.max(np.abs(physical-constrained))),muscle_audit=muscles,
            repeat_effective_energy_difference_j=repeat['effective_energy_j']-base['effective_energy_j'],repeat_maximum_actual_q_difference=float(np.max(np.abs(np.array(repeat['independent_q'])-base['independent_q']))),all_sampled_branches_valid=all_branches_valid,
            component_gradients={k:gradient(lambda a:a[k]).tolist() for k in ('gravity_energy_j','joint_energy_j','skin_energy_j','bed_energy_j','muscle_passive_energy_j','muscle_active_effective_energy_j')},
            gradient_gate_passed=all_branches_valid and max(map(abs,errors))<=1e-4 and max(component_errors.values())<=1e-4 and all(m['length_moment_arm_error_m']<=1e-6 and m['energy_virtual_work_error']<=1e-4 and m['fiber_minimum_margin_m']>1e-10 and m['internal_stiffness_n_m']>0 and abs(m['fiber_velocity_m_s'])<=1e-5 for m in muscles),
            gate_scope='Predeclared numerical gradient fixture tolerances in native coordinate units; not a change to physical equilibrium criteria')
    except Exception as error:report.update(status='diagnostic_failed',error=type(error).__name__+': '+str(error)[:65536])
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if process is not None:
            try:
                if process.poll() is None:process.stdin.write('close\n');process.stdin.flush();process.wait(timeout=2)
            except Exception:process.kill();process.wait(timeout=2)
        report.update(evaluations=len(records),wall_s=time.monotonic()-started,harness_sha256=sha(Path(__file__)));write('report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('muscle_audit','component_gradients','coordinate_names','energy_gradient','physical_gradient','gradient_errors')}|{'output':str(output)},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-native',action='store_true');p.add_argument('build',type=Path);a=p.parse_args()
    if not a.run_native:raise SystemExit('Coordinated native fixture grant required')
    run(a.build.resolve())
