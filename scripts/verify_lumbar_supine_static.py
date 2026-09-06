"""Two gated copied-state diagnostics for separately pinned98-muscle physics."""
from pathlib import Path
import argparse,hashlib,json,signal,sys,tempfile,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from frozen_static_stream import load,validate
from solve_constrained_supine_pose import dynamic_metric


def run(path):
    case=json.loads(path.read_text())
    if case['schema']!='ihm.lumbar-supine-static-reference.v1' or case['cached_responses_reused']!=0:raise ValueError('Fresh98-muscle reference required')
    for name,digest in case['files'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Reference source/input changed: '+name)
    Native=load(ROOT/case['frozen_engine']);output=Path(tempfile.mkdtemp(prefix='lumbar-supine-static-',dir=ROOT/'data/derived'));stream=None;started=time.monotonic()
    report=dict(reference_manifest=str(path.relative_to(ROOT)),reference_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),native_evaluations=0,maximum_evaluations=2,maximum_wall_s=15,physical_time_advanced_s=0,cached_responses_reused=0,accepted_equilibrium=False)
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,frame):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('98-muscle diagnostic exceeded15s')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,15)
    try:
        stream=Native(ROOT,output/'native',environment='supine',target_mass_kg=case['mass_kg'],augmented_registration=case['variant_registration'],surface_contact_manifest=case['surface_manifest'],bed_material=case['material'])
        before=stream.snapshot();write('initial_native.json',before)
        if len(before['muscles'])!=98:raise ValueError('Native did not instantiate98 muscles')
        seed=json.loads((ROOT/case['old_q_seed']).read_text())['coordinates']
        neutral={'pelvis_tx':before['coordinates']['pelvis_tx']['value']}
        summary={}
        for label,q in [('neutral',neutral),('old_q_seed',seed)]:
            command=['evaluate_static_pose',str(len(q))]
            for name,value in q.items():command.extend([name,str(value)])
            write('pending_request.json',dict(label=label,coordinates=q));response=stream._request(' '.join(command));report['native_evaluations']+=1;write(label+'.json',response)
            dynamic_metric(response)
            names=response['mobility_coordinate_names'];residual=response['constrained_zero_acceleration_residual_mobility_force'];hip=names.index('hip_rotation_r')
            summary[label]=dict(maximum_abs_acceleration=max(map(abs,response['udot'])),contact_force_n=response['contact_force_n'],dominant_residuals=sorted(zip(names,residual),key=lambda v:-abs(v[1]))[:10],hip_rotation_r_acceleration=response['udot'][hip],hip_rotation_r_residual_nm=residual[hip],hip_rotation_r_q=response['coordinates']['hip_rotation_r']['value'],hip_direction_scope='Instantaneous zero-speed acceleration only; not a perturbation stability test',constraints={key:response[key] for key in ('constraint_position_error','constraint_velocity_error','constraint_acceleration_error')})
            if not response['continuing_state_unchanged'] or response['physical_time_advanced_s']!=0:raise ValueError('Static protocol mutated live state/time')
        live=stream._request('observe');report['continuing_state_unchanged']=all(before[k]==live[k] for k in before if k!='kind')
        if not report['continuing_state_unchanged']:raise ValueError('Diagnostic altered continuing state')
        report.update(status='two_pose_diagnostic_complete',diagnostic_passed=True,poses=summary)
    except Exception as error:
        report.update(status='diagnostic_failed',diagnostic_passed=False,error=type(error).__name__+': '+str(error)[:65536])
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        validate(ROOT/case['frozen_engine']);report['frozen_archive_intact']=True;report['wall_s']=time.monotonic()-started;report['harness_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();write('report.json',report)
        print(json.dumps({**report,'output':str(output)},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-native',action='store_true');p.add_argument('reference',type=Path);a=p.parse_args()
    if not a.run_native:raise SystemExit('Coordinated --run-native grant required')
    run(a.reference.resolve())
