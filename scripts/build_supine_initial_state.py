"""Build a hashed static-pose seed and diagnosis, never label a transient equilibrated.

This source-only builder does not load a native library or start a simulation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import signal
import time
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from verify_supine_support import LIMITS, metrics


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def analyze(snapshot, model):
    observed = metrics(snapshot)
    source = ET.fromstring(model)
    dependent = {c.findtext('dependent_coordinate_name') for c in source.iter('CoordinateCouplerConstraint')
                 if c.findtext('isEnforced', 'true') == 'true'}
    coordinates = []
    for coordinate in source.iter('Coordinate'):
        name = coordinate.attrib['name']
        coordinates.append(dict(name=name, default_value=float(coordinate.findtext('default_value')),
                                range=list(map(float, coordinate.findtext('range').split())),
                                observed_value=snapshot['coordinates'][name]['value'],
                                dependent=name in dependent,
                                locked=coordinate.findtext('locked', 'false') == 'true',
                                prescribed=coordinate.findtext('prescribed', 'false') == 'true'))
    contacts = []
    normal = -np.asarray(snapshot['gravity_m_s2']) / np.linalg.norm(snapshot['gravity_m_s2'])
    if not np.allclose(normal, [1,0,0], atol=1e-12):
        raise ValueError('Current posterior proxy contract requires source +X normal')
    for contact in snapshot['contacts']:
        gap = contact['center_m'][0] - contact['radius_m'] - snapshot['support_plane_source_x_m']
        contacts.append(dict(body=contact['body_frame'], radius_m=contact['radius_m'], gap_m=gap,
                             normal_force_n=float(np.dot(contact['force_n'],normal))))
    contacts.sort(key=lambda c:c['gap_m'])
    muscle_elements = {m.attrib['name']:m for m in source.iter() if m.tag.endswith('Muscle') and m.attrib.get('name') in snapshot['muscles']}
    muscles = []
    for name, muscle in snapshot['muscles'].items():
        element = muscle_elements[name]
        muscles.append(dict(name=name, model=element.tag, activation=muscle['activation'],
                            excitation=muscle['excitation'], default_activation=None if element.findtext('default_activation') is None else float(element.findtext('default_activation')),
                            fiber_length_over_optimal=muscle['fiber_length_m']/muscle['optimal_fiber_length_m'],
                            tendon_force_over_max_isometric=muscle['tendon_force_n']/muscle['max_isometric_force_n'],
                            source_tendon_slack_length_m=float(element.findtext('tendon_slack_length')),
                            fiber_velocity_m_s=muscle['fiber_velocity_m_s']))
    muscles.sort(key=lambda m:m['tendon_force_over_max_isometric'], reverse=True)
    pelvis = next((c for c in contacts if c['body']=='pelvis'), None)
    contact_diagnosis = dict(plane_source_x_m=snapshot['support_plane_source_x_m'], contacts=contacts,
        positive_force_bodies=[c['body'] for c in contacts if c['normal_force_n']>1e-6],
        largest_gap_m=max(c['gap_m'] for c in contacts),
        rigid_translation_preserves_relative_gaps=True,
        plane_shift_to_touch_pelvis_m=None if pelvis is None else max(0.,pelvis['gap_m']),
        minimum_proxy_penetration_after_that_shift_m=None if pelvis is None else max(0.,pelvis['gap_m']-contacts[0]['gap_m']))
    return dict(schema='ihm.supine-initialization-analysis.v1', observed_time_s=snapshot['time_s'],
                accepted_initial_state=False, acceptance_reason='A pose observation is not a sustained equilibrium receipt',
                metrics=observed, contact_geometry=contact_diagnosis, coordinates=coordinates, muscles=muscles,
                limitations=['Fiber length/optimal is not passive tissue strain or evidence of muscle injury',
                             'Tendon force includes active and passive contributions; these observations do not separate them',
                             'Source joint ranges are optimization bounds, not anatomical validation',
                             'COM/inertia spheres are engineering geometry and do not reproduce posterior skin'])


def recipe(analysis, snapshot, identity):
    return dict(schema='ihm.supine-static-seed.v1', status='seed_only_requires_native_static_solve',
                source_identity=identity, accepted_initial_state=False,
                seed_coordinates={c['name']:c['observed_value'] for c in analysis['coordinates']
                                  if not any(c[k] for k in ('dependent','locked','prescribed'))},
                coordinate_bounds={c['name']:c['range'] for c in analysis['coordinates']
                                   if not any(c[k] for k in ('dependent','locked','prescribed'))},
                fixed_excitation={name:m['excitation'] for name,m in snapshot['muscles'].items()},
                gravity_m_s2=snapshot['gravity_m_s2'],
                support_plane_source_x_m=snapshot['support_plane_source_x_m'],
                static_evaluation={'velocity':'zero only inside isolated static residual evaluation',
                    'muscle_state':'equilibrate muscles at unchanged excitation for each candidate q',
                    'constraints':'assemble dependent coordinates; do not independently optimize them',
                    'residual':'full native generalized accelerations or generalized force balance, including unconstrained root DOFs',
                    'original_state':'restore exactly after every candidate; no published advancement or metabolic integration'},
                bounded_solver_proposal={'method':'bounded nonlinear least squares of static residuals',
                    'maximum_evaluations':200,'maximum_wall_s':60,
                    'accepted_if':'native full residual converges, then unchanged forward acceptance criteria pass',
                    'forbidden_shortcuts':['gravity cancellation','contact preload forces','muscle disabling',
                                          'contact/solver coefficient changes','labeling a zeroed dynamic snapshot settled']},
                forward_acceptance_criteria=LIMITS,
                native_api_required=['isolated evaluate_static_pose with per-coordinate acceleration/force residual',
                                     'accepted full-state export/import bound to exact model/contact/muscle identities'])


def fixture_check():
    run = ROOT/'data/derived/supine-support-5ma720yd'
    snapshot = json.loads((run/'initial_native.json').read_text())
    source = (run/'plant/native/inputs/subject_walk_scaled.osim').read_bytes()
    report = analyze(snapshot, source)
    assert report['accepted_initial_state'] is False
    assert report['contact_geometry']['positive_force_bodies'] == ['torso']
    assert report['contact_geometry']['minimum_proxy_penetration_after_that_shift_m'] > .23
    assert report['metrics']['balance_relative_weight'] > .97
    plan = recipe(report, snapshot, {})
    assert not plan['accepted_initial_state'] and len(plan['fixed_excitation']) == 92
    assert not {'knee_angle_r_beta','knee_angle_l_beta'} & plan['seed_coordinates'].keys()
    shifted = json.loads(json.dumps(snapshot))
    shifted['support_plane_source_x_m'] += 2.
    for contact in shifted['contacts']:
        contact['center_m'][0] += 2.
    shifted_report = analyze(shifted, source)
    assert np.allclose([c['gap_m'] for c in report['contact_geometry']['contacts']],
                       [c['gap_m'] for c in shifted_report['contact_geometry']['contacts']],atol=1e-12)
    assert all(m['activation']==m['default_activation'] for m in report['muscles'] if m['default_activation'] is not None)
    assert all(m['activation']==.05 for m in report['muscles'] if m['model']=='Thelen2003Muscle')
    return dict(passed=True, native_run=False, checks=['retained torso-only support','large pelvis gap',
                'equilibrium rejection','dependent coordinate exclusion','translation invariance','unchanged excitation'])


def solve_native(args, seed):
    """At most 200 actual candidate evaluations; never advance physical time."""
    sys.path.insert(0,str(ROOT))
    from ihm.native.mechanical_stream import NativeMechanicalStream
    from scipy.optimize import least_squares
    output=Path(tempfile.mkdtemp(prefix='supine-static-solve-',dir=ROOT/'data/derived'))
    started=time.monotonic(); stream=None; evaluations=[]; best=None
    report=dict(status='incomplete',accepted_initial_state=False,physical_time_advanced_s=0,
                maximum_evaluations=200,maximum_wall_s=60,
                scope='Static residual search only; no forward-equilibrium acceptance')
    names=list(seed['seed_coordinates'])
    lower=np.array([seed['coordinate_bounds'][n][0] for n in names])
    upper=np.array([seed['coordinate_bounds'][n][1] for n in names])
    initial=np.array([seed['seed_coordinates'][n] for n in names])
    if not np.isfinite(initial).all() or (initial<lower).any() or (initial>upper).any():
        raise ValueError('Static seed outside source bounds')
    def write(name,value):
        (output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,stack):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Static residual search exceeded 60 second budget')
    def evaluate(values):
        nonlocal best
        if len(evaluations)>=200:raise RuntimeError('Static residual search reached 200 actual evaluations')
        if not np.isfinite(values).all() or (values<lower).any() or (values>upper).any():
            raise ValueError('Optimizer candidate outside source bounds')
        command=['evaluate_static_pose',str(len(names))]
        for name,value in zip(names,values):command.extend([name,str(float(value))])
        answer=stream._request(' '.join(command))
        if answer['physical_time_advanced_s']!=0 or not answer['continuing_state_unchanged'] or np.linalg.norm(answer['u'])>1e-12:
            raise AssertionError('Static evaluator changed time/state or retained speed')
        # Translational acceleration scaled by g, angular by g/0.3 m.
        # This changes optimization conditioning only, never native force laws.
        residual=np.array([answer['coordinates'][n]['acceleration']/(9.81/.3 if answer['coordinates'][n]['rotational'] else 9.81) for n in names])
        if not np.isfinite(residual).all():raise ValueError('Nonfinite static residual')
        norm=float(np.linalg.norm(residual))
        row=dict(evaluation=len(evaluations)+1,wall_s=time.monotonic()-started,
                 residual_norm=norm,maximum_abs_udot=float(np.max(np.abs(answer['udot']))),
                 maximum_penetration_m=answer['maximum_penetration_m'])
        evaluations.append(row);write('evaluations.json',evaluations)
        if best is None or norm<best['residual_norm']:
            best=dict(residual_norm=norm,requested_coordinates=dict(zip(names,map(float,values))),native=answer)
            write('best_candidate.json',best)
        return residual
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,60)
    try:
        write('seed.json',seed)
        stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=77.6122029,
             augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
        checkpoint=stream.checkpoint();before=stream.snapshot()
        evaluate(initial)
        observed=stream._request('observe')
        report['candidate_left_live_state_unchanged']=all(before[k]==observed[k] for k in before if k!='kind')
        if not report['candidate_left_live_state_unchanged']:raise AssertionError('Evaluation mutated continuing state')
        restored=stream.restore(checkpoint)
        report['candidate_evaluation_restore_matches']=all(before[k]==restored[k] for k in before if k!='kind')
        if not report['candidate_evaluation_restore_matches']:raise AssertionError('Evaluation changed continuing state')
        # Rejection behavior is part of the native contract, not an optimizer retry.
        for invalid in ('evaluate_static_pose 1 knee_angle_r_beta 0','evaluate_static_pose 1 pelvis_tx 1000'):
            try:stream._request(invalid)
            except ValueError:pass
            else:raise AssertionError('Invalid static-pose coordinate accepted')
        result=least_squares(evaluate,initial,bounds=(lower,upper),max_nfev=200,
                             ftol=1e-8,xtol=1e-8,gtol=1e-8,diff_step=1e-4,x_scale='jac')
        report.update(optimizer_success=bool(result.success),optimizer_message=str(result.message))
        report['status']='static_candidate_converged' if best and max(abs(v) for v in best['native']['udot'])<=1e-4 and best['native']['constraint_position_error']<=1e-5 and best['native']['maximum_penetration_m']<=LIMITS['penetration_m'] else 'static_residual_unresolved'
        restored=stream.restore(checkpoint)
        report['final_restore_matches']=all(before[k]==restored[k] for k in before if k!='kind')
        if not report['final_restore_matches']:raise AssertionError('Static search changed continuing state')
        stream.release(checkpoint)
    except Exception as error:
        report.update(status='stopped',error=f'{type(error).__name__}: {error}')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        report.update(wall_s=time.monotonic()-started,evaluations=len(evaluations),
                      best_residual_norm=None if best is None else best['residual_norm'])
        write('report.json',report);print(json.dumps({**report,'output':str(output)},indent=2))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=ROOT/'data/derived/supine-support-5ma720yd')
    parser.add_argument('--fixture-check',action='store_true')
    parser.add_argument('--solve-native',action='store_true',help='Requires explicit coordinated heavy-resource slot')
    args=parser.parse_args()
    if args.fixture_check:
        print(json.dumps(fixture_check(),indent=2));return
    run=args.run.resolve()
    if not run.is_relative_to(ROOT):
        raise ValueError('Owned retained run required')
    model=run/'plant/native/inputs/subject_walk_scaled.osim'
    snapshot_path=run/'initial_native.json'
    execution_path=run/'plant/native/execution.json'
    execution=json.loads(execution_path.read_text())
    if execution['source_sha256']['subject_walk_scaled.osim'] != sha(model):
        raise ValueError('Retained native model identity mismatch')
    snapshot=json.loads(snapshot_path.read_text())
    analysis=analyze(snapshot,model.read_bytes())
    identity={str(p.relative_to(ROOT)):sha(p) for p in (model,snapshot_path,execution_path,Path(__file__).resolve())}
    if args.solve_native:
        solve_native(args,recipe(analysis,snapshot,identity));return
    output=Path(tempfile.mkdtemp(prefix='supine-initialization-',dir=ROOT/'data/derived'))
    for name,value in [('analysis.json',analysis),('static_seed.json',recipe(analysis,snapshot,identity))]:
        (output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(output),native_run=False,accepted_initial_state=False,
                         metrics=analysis['metrics'],contact_geometry=analysis['contact_geometry']),indent=2))


if __name__=='__main__':
    main()
