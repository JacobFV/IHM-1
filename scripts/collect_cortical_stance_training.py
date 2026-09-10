#!/usr/bin/env python3
"""Collect four axis-only native engineering-teacher trajectories for cortex training."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def state_digest(state):
    # Native message kind changes from initialized to restored; physics must agree.
    return hashlib.sha256(json.dumps({k:v for k,v in state.items() if k != 'kind'}, sort_keys=True).encode()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, allow_nan=False) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='data/runtime/motor-learning/patient-stance-axis-data-20260908')
    args = parser.parse_args()
    out = (ROOT / args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    registration_path = ROOT / 'data/models/engineering_stance_v1/registration.json'
    artifact_path = registration_path.parent / 'linearization.npz'
    registration = json.loads(registration_path.read_text())
    mass, dt, seconds = MECHANICAL_TARGET_MASS_KG, .01, 3.
    teacher = NativeStanceLQR(artifact_path, model_sha256=registration['model_sha256'], dt_s=dt, target_mass_kg=mass)
    sources = [Path(__file__), ROOT / 'ihm/native/mechanical_stream.py', ROOT / 'ihm/native/stance_lqr.py', registration_path, artifact_path]
    source_hashes = {str(p.relative_to(ROOT)): digest(p) for p in sources}
    retained = out / 'sources'
    retained.mkdir()
    for path in sources:
        (retained / path.name).write_bytes(path.read_bytes())
    report = {
        'schema': 'ihm.cortical-stance-axis-training.v1',
        'model_sha256': registration['model_sha256'], 'artifact_sha256': digest(artifact_path),
        'source_sha256': source_hashes, 'target_mass_kg': mass, 'dt_s': dt,
        'seconds_requested_per_arm': seconds,
        'teacher': 'NativeStanceLQR, gain=1, exact artifact baseline preserved',
        'muscle_names': teacher.muscle_names, 'state_names': teacher.state_names,
        'heldout_excluded': 'patient98_margin_negative_com_push; diagonal (-5,0,-5) N is never collected or read',
        'force_basis': 'GROUND; point is current pelvis local COM transformed by body_point each step',
        'frame_semantics': 'Each frame is complete PRE-STEP native snapshot, with commands computed from that frame; next_state_time_s is after advance. Trainer must recompute targets from state.',
        'arms': {},
    }
    stream = NativeMechanicalStream(ROOT, out / 'plant', environment='upright', target_mass_kg=mass,
        augmented_registration=str(registration_path), initial_pose=registration.get('initial_pose'))
    initial = stream.snapshot()
    token = stream.checkpoint()
    write_json(out / 'initial_state.json', initial)
    initial_hash = state_digest(initial)
    report['initial_state_sha256'] = initial_hash
    started = time.monotonic()
    try:
        for arm, force in [('positive_x', [5.,0.,0.]), ('negative_x', [-5.,0.,0.]),
                           ('positive_z', [0.,0.,5.]), ('negative_z', [0.,0.,-5.])]:
            stream.restore(token)
            restored = stream.snapshot()
            assert state_digest(restored) == initial_hash, [k for k in restored if k != 'kind' and restored[k] != initial.get(k)]
            frames, error, pulse_steps = [], None, 0
            try:
                for index in range(round(seconds / dt)):
                    state = stream.snapshot()
                    commands, diagnostic = teacher.commands(state)
                    forces, point_receipt = [], None
                    if 1. - 1e-9 <= state['time_s'] < 1.1 - 1e-9:
                        point_receipt = stream.body_point(body='pelvis', station_m=initial['bodies']['pelvis']['mass_center_local_m'])
                        forces = [{'body':'pelvis', 'point_m':point_receipt['point_source_m'], 'force_n':force}]
                        pulse_steps += 1
                    following = stream.advance(dt, actuation=commands, forces=forces)
                    frame = dict(state)
                    frame.update(commands=commands, diagnostics=diagnostic, external_forces=forces,
                        force_point_receipt=point_receipt, next_state_time_s=following['time_s'])
                    frames.append(frame)
                    if (index + 1) % 50 == 0:
                        print(json.dumps({'arm':arm,'time_s':following['time_s']}), flush=True)
            except Exception as exc:
                error = type(exc).__name__ + ': ' + str(exc)
            final = stream.snapshot()
            write_json(out / (arm + '.json'), frames)
            write_json(out / (arm + '-final.json'), final)
            report['arms'][arm] = {'frames_path':str(out / (arm + '.json')), 'frame_count':len(frames),
                'force_n':force, 'pulse_start_s':1., 'pulse_end_s':1.1, 'pulse_steps':pulse_steps,
                'completed_horizon': final['time_s'] >= seconds - 1e-9, 'final_time_s':final['time_s'], 'error':error,
                'frames_sha256':digest(out / (arm + '.json'))}
            write_json(out / 'report.json', report)
    finally:
        stream.release(token)
        stream.close()
    report['wall_seconds'] = time.monotonic() - started
    report['all_four_complete'] = len(report['arms']) == 4 and all(a['completed_horizon'] and a['error'] is None and a['pulse_steps'] == 10 for a in report['arms'].values())
    write_json(out / 'report.json', report)
    print(json.dumps(report, indent=2), flush=True)
    if not report['all_four_complete']:
        raise SystemExit('Incomplete training collection: inspect retained report')


if __name__ == '__main__':
    main()
