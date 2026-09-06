"""Bounded equilibrium diagnostic. Default execution is native-free fixture checking."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import math
import signal
import sys
import tempfile
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
STAGES = (.02, .1, .3, .6, 1.)
LIMITS = dict(balance_relative_weight=.02, acceleration_g=.02,
              kinetic_j_kg=5e-5, kinetic_slope_j_kg_s=1e-4,
              com_speed_m_s=.01, angular_speed_rad_s=.05,
              penetration_m=.01, constraint_error=1e-5,
              audit_relative_weight=1e-7, window_s=.2, earliest_s=.4)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def metrics(frame, previous=None):
    """Native source frame only; no extra inertia from canonical attachments."""
    mass = float(frame['mass_kg'])
    gravity = np.asarray(frame['gravity_m_s2'], float)
    g = float(np.linalg.norm(gravity))
    if frame['environment'] != 'supine' or mass <= 0 or g <= 0:
        raise ValueError('Positive mass/gravity and supine source frame required')
    support = np.asarray(frame['contact_force_n'], float)
    residual = np.asarray(frame['momentum_balance_residual_n'], float)
    # No applied loads in this experiment. The native audit is M(a-g)-contact.
    acceleration = gravity + (support + residual) / mass
    velocity = np.zeros(3)
    max_angular = 0.
    for body in frame['bodies'].values():
        rotation = np.asarray(body['transform_ground'], float)[:3, :3]
        offset = rotation @ np.asarray(body['mass_center_local_m'], float)
        omega = np.asarray(body['angular_velocity_rad_s'], float)
        velocity += body['mass_kg'] * (np.asarray(body['origin_velocity_m_s']) + np.cross(omega, offset)) / mass
        max_angular = max(max_angular, float(np.linalg.norm(omega)))
    penetration = max((max(0., frame['support_plane_source_x_m'] + c['radius_m'] - c['center_m'][0]) for c in frame['contacts']), default=0.)
    row = dict(time_s=frame['time_s'], support_relative_weight=float(np.dot(support, -gravity/g)/(mass*g)),
               balance_relative_weight=float(np.linalg.norm(support+mass*gravity)/(mass*g)),
               acceleration_g=float(np.linalg.norm(acceleration)/g),
               audit_relative_weight=float(np.linalg.norm(residual)/(mass*g)),
               kinetic_j_kg=frame['kinetic_energy_j']/mass,
               com_velocity_m_s=velocity.tolist(), com_speed_m_s=float(np.linalg.norm(velocity)),
               angular_speed_rad_s=max_angular, penetration_m=penetration,
               constraint_error=max(abs(frame['constraint_position_error']), abs(frame['constraint_velocity_error'])))
    if previous is not None:
        dt = row['time_s'] - previous['time_s']
        if dt <= 0:
            raise ValueError('Samples must advance time')
        row['finite_difference_com_acceleration_g'] = float(np.linalg.norm(velocity - previous['com_velocity_m_s']) / (dt*g))
    if not all(math.isfinite(float(x)) for k, x in row.items() if k != 'com_velocity_m_s') or not np.isfinite(velocity).all():
        raise ValueError('Nonfinite diagnostic')
    return row


def classify(rows):
    last = rows[-1]
    divergent = []
    for key, ceiling in dict(support_relative_weight=10., acceleration_g=20., kinetic_j_kg=2.5,
                             penetration_m=.05, constraint_error=.01, audit_relative_weight=.01).items():
        if abs(last[key]) > ceiling:
            divergent.append(key)
    window = [r for r in rows if r['time_s'] >= last['time_s'] - LIMITS['window_s'] - 1e-9]
    span = window[-1]['time_s'] - window[0]['time_s']
    if span >= LIMITS['window_s'] - 1e-9 and last['kinetic_j_kg'] > .1 and last['kinetic_j_kg'] > 4*max(window[0]['kinetic_j_kg'], .001):
        divergent.append('kinetic_growth')
    if divergent:
        return dict(status='diverged', reasons=divergent)
    if last['time_s'] < LIMITS['earliest_s'] or span < LIMITS['window_s'] - 1e-9:
        return dict(status='transient', reasons=['insufficient sustained observation'])
    slope = float(np.polyfit([r['time_s'] for r in window], [r['kinetic_j_kg'] for r in window], 1)[0])
    checks = {key: max(r[key] for r in window) <= LIMITS[key] for key in (
        'balance_relative_weight', 'acceleration_g', 'kinetic_j_kg', 'com_speed_m_s',
        'angular_speed_rad_s', 'penetration_m', 'constraint_error', 'audit_relative_weight')}
    checks['kinetic_slope'] = abs(slope) <= LIMITS['kinetic_slope_j_kg_s']
    checks['finite_difference_acceleration'] = all(r.get('finite_difference_com_acceleration_g', math.inf) <= LIMITS['acceleration_g'] for r in window)
    return dict(status='converged' if all(checks.values()) else 'unsettled', checks=checks,
                kinetic_slope_j_kg_s=slope, window_observed_s=span)


def fixture_check():
    body = dict(mass_kg=10., transform_ground=np.eye(4).tolist(), mass_center_local_m=[0,0,0],
                angular_velocity_rad_s=[0,0,0], origin_velocity_m_s=[0,0,0])
    frame = dict(time_s=0., mass_kg=10., gravity_m_s2=[-10.,0,0], contact_force_n=[0,0,0],
                 momentum_balance_residual_n=[0,0,0], environment='supine', kinetic_energy_j=0.,
                 bodies={'body':body}, support_plane_source_x_m=0.,
                 contacts=[dict(radius_m=.1, center_m=[.1,0,0])],
                 constraint_position_error=0., constraint_velocity_error=0.)
    falling = metrics(frame)
    assert falling['audit_relative_weight'] == 0 and falling['balance_relative_weight'] == 1 and falling['acceleration_g'] == 1
    frame['contact_force_n'] = [100,0,0]
    rows = []
    for i in range(81):
        frame['time_s'] = i*.005
        rows.append(metrics(frame, rows[-1] if rows else None))
    assert classify(rows[:2])['status'] == 'transient'
    assert classify(rows)['status'] == 'converged'
    moving = copy.deepcopy(rows)
    moving[-1]['kinetic_j_kg'] = .01
    assert classify(moving)['status'] == 'unsettled'
    broken = copy.deepcopy(rows)
    broken[-1]['penetration_m'] = .06
    assert classify(broken)['status'] == 'diverged'
    growth = copy.deepcopy(rows)
    growth[-1]['kinetic_j_kg'] = .2
    assert 'kinetic_growth' in classify(growth)['reasons']
    frame['contacts'][0]['center_m'][0] = .08
    assert abs(metrics(frame)['penetration_m']-.02) < 1e-12
    frame['kinetic_energy_j'] = float('nan')
    try:
        metrics(frame)
    except ValueError:
        pass
    else:
        raise AssertionError('Nonfinite input accepted')
    retained = json.loads((ROOT/'data/derived/articulated-acceptance-noo34jds/baseline.json').read_text())
    force = np.array(retained['body_environment']['contact_force_source_n'])
    mass = sum(b['mass_kg'] for b in retained['native_bodies'].values())
    # Source gravity points along -X; retained canonical gravity must not be mixed with source forces.
    retained_balance = float(np.linalg.norm(force + np.array([-9.81*mass,0,0]))/(9.81*mass))
    assert retained_balance > .02
    return dict(fixture_checks_passed=True, retained_initial_transient_balance_relative_weight=retained_balance,
                native_run=False, scope='Classifier fixtures and retained 2 ms transient; no equilibrium acceptance')


def run_native(args):
    from ihm.assembly.articulated import ArticulatedBodyPlant
    output = Path(tempfile.mkdtemp(prefix='supine-support-', dir=ROOT/'data/derived'))
    started = time.monotonic()
    plant = None
    rows, receipts = [], []
    report = dict(passed=False, status='incomplete', criteria=LIMITS, stages_s=STAGES,
                  wall_budget_s=args.wall_budget_s, sample_dt_s=.005,
                  scope='Engineering proxy supine static-support diagnostic; no anatomical or biological validation')

    def write(name, value):
        (output/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')

    def deadline(signum, stack):
        if plant is not None and plant.native.process.poll() is None:
            plant.native.process.kill()
        raise TimeoutError('Supine support wall budget exhausted')

    old_handler = signal.signal(signal.SIGALRM, deadline)
    signal.setitimer(signal.ITIMER_REAL, args.wall_budget_s)
    try:
        write('protocol.json', report)
        plant = ArticulatedBodyPlant(ROOT, output/'plant', environment='supine', target_mass_kg=77.6122029,
                                     augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
        native = plant.native.snapshot()
        assert len(native['muscles']) == 92
        rows.append(metrics(native))
        write('initial_native.json', native)
        for endpoint in STAGES:
            checkpoint = plant.checkpoint()
            before = plant.native.snapshot()
            receipt = dict(stage_target_s=endpoint, start_time_s=before['time_s'], checkpoint_snapshot_sha256=digest(before),
                           checkpoint_scope='Opaque live SimTK State; saved JSON is observation evidence, not a restartable checkpoint')
            receipts.append(receipt)
            write('checkpoint_receipts.json', receipts)
            # Verify an actual restore at each boundary without replaying integration.
            plant.restore(checkpoint)
            restored = plant.native.snapshot()
            receipt['restore_matches'] = digest({k:v for k,v in before.items() if k != 'kind'}) == digest({k:v for k,v in restored.items() if k != 'kind'})
            if not receipt['restore_matches']:
                raise AssertionError('Native checkpoint restore changed state')
            while rows[-1]['time_s'] < endpoint - 1e-10:
                plant.advance(min(.005, endpoint-rows[-1]['time_s']))
                native = plant.native.snapshot()
                rows.append(metrics(native, rows[-1]))
                result = classify(rows)
                write('samples.json', rows)
                if result['status'] == 'diverged':
                    break
            write(f'stage_{endpoint:g}_native.json', native)
            receipt.update(end_time_s=native['time_s'], end_snapshot_sha256=digest(native), outcome=result,
                           wall_s=time.monotonic()-started)
            plant.release(checkpoint)
            receipt['released'] = True
            write('checkpoint_receipts.json', receipts)
            if result['status'] in ('diverged', 'converged'):
                break
        report.update(result, passed=result['status']=='converged', final_time_s=rows[-1]['time_s'])
    except Exception as error:
        report.update(status='stopped', error=f'{type(error).__name__}: {error}',
                      last_successful_time_s=rows[-1]['time_s'] if rows else None,
                      failed_stage_target_s=receipts[-1]['stage_target_s'] if receipts else None)
        if receipts:
            receipts[-1].update(outcome={'status':'stopped', 'error':report['error']},
                                last_successful_time_s=report['last_successful_time_s'])
            write('checkpoint_receipts.json', receipts)
        write('failure.json', {**report, 'traceback':traceback.format_exc()})
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)
        if plant is not None:
            plant.close()
        report.update(wall_s=time.monotonic()-started, samples=len(rows),
                      harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        write('report.json', report)
        print(json.dumps({**report, 'output':str(output)}, indent=2))
    return 0 if report['passed'] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-native', action='store_true', help='Requires a coordinated heavy-resource slot')
    parser.add_argument('--wall-budget-s', type=float, default=60.)
    args = parser.parse_args()
    if not math.isfinite(args.wall_budget_s) or not 1 <= args.wall_budget_s <= 60:
        parser.error('wall budget must be within [1,60] seconds')
    if args.run_native:
        return run_native(args)
    print(json.dumps(fixture_check(), indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
