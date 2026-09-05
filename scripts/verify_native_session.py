"""Exercise persistent native ownership, action timing, and invalid-command isolation."""
import json
import hashlib
from pathlib import Path
import tempfile
from ihm.native.session import SessionConfig, NativeSession, Meal
from ihm.native import NativeConfig, run_native, load_trajectory


def main():
    for value in (False, float('nan'), -.02, .03, 86400.02):
        try:
            SessionConfig(horizon_s=value)
        except ValueError:
            pass
        else:
            raise AssertionError(f'accepted invalid horizon {value}')
    for kwargs in ({'carbohydrate_g': -1}, {'water_ml': float('inf')}, {'name': 'bad name'}, {}):
        try:
            Meal(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError('accepted invalid meal')
    root = Path(__file__).resolve().parents[1]
    state = root/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'
    out = Path(tempfile.mkdtemp(prefix='session-', dir=root/'data/derived/audits'))
    config = SessionConfig(state_path=state, engine_variant='whole_body_integrity_energy', horizon_s=180)
    batch_dir = out/'batch'
    stream_dir = out/'stream'
    run_native(NativeConfig(seconds=1, state_path=state, engine_variant=config.engine_variant, sample_hz=50), batch_dir)
    batch = load_trajectory(batch_dir)
    with NativeSession(config, stream_dir) as parity:
        launch=json.loads((stream_dir/'execution-inputs.json').read_text())
        environment=json.loads((stream_dir/'environment-inputs/manifest.json').read_text())
        assert environment['capture']['stage']=='pre_start'
        detached=Path(launch['selected_resource_tree'])
        assert detached==stream_dir/'runtime-resources'
        for name in ('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml'):
            assert (stream_dir/name).resolve()==detached/name
        selected_state=Path(launch['selected_patient_input'])
        assert selected_state==stream_dir/'input-state.xml' and not selected_state.is_symlink()
        assert selected_state.stat().st_ino!=state.stat().st_ino
        assert hashlib.sha256(selected_state.read_bytes()).hexdigest()==launch['selected_patient_input_sha256']
        command=json.loads((stream_dir/'manifest.json').read_text())['command']
        assert command[2]==str(selected_state)
        last = parity.step(1)
        errors = {a: abs(last['values'][a]-batch['values'][b][-1]) for a,b in
                  [('lung_volume_ml','TotalLungVolume(mL)'),('heart_rate_per_min','HeartRate(1/min)'),
                   ('mean_arterial_pressure_mmhg','MeanArterialPressure(mmHg)')]}
        assert max(errors.values()) < 1e-5, errors
    with NativeSession(SessionConfig(state_path=state, horizon_s=.3), out/'fractional-horizon') as fractional:
        for _ in range(3):last=fractional.step(.1)
        assert last['elapsed_s']==.3
    # An acknowledged action must produce native demand, not merely carry an
    # intensity scalar on an SEExercise whose mode is still NONE.
    exercise_observations = {}
    for label, intensity in [('rest', 0), ('exercise', .15)]:
        with NativeSession(config, out/label) as demand:
            demand.exercise(intensity)
            exercise_observations[label] = demand.step(120)['values']
            values=exercise_observations[label]
            target=intensity*values['maximum_work_rate_w']
            assert abs(values['exercise_energy_demand_w']-target)<1e-5
            demand.exercise(0)
            assert demand.step(.02)['values']['exercise_energy_demand_w']==0, 'stopped demand remained active'
    demand_delta = {key: exercise_observations['exercise'][key]-exercise_observations['rest'][key]
                    for key in ('metabolic_rate_w', 'oxygen_consumption_ml_per_min')}
    assert all(value > 1 for value in demand_delta.values()), demand_delta
    pending=NativeSession(config,out/'pending-close')
    pending.meal(Meal(carbohydrate_g=1))
    try:pending.close()
    except ValueError:pass
    else:raise AssertionError('graceful close discarded an acknowledged pending meal')
    assert pending.process.poll() is not None
    raw=NativeSession(config,out/'raw-eof')
    raw.meal(Meal(carbohydrate_g=1))
    raw.process.stdin.close()
    assert raw.process.wait(timeout=10)==9, 'raw EOF accepted an unconsumed action'
    raw.close(graceful=False)
    out = out/'actions'
    with NativeSession(config, out) as body:
        initial = body.snapshot()
        for invalid in (.03, -1, True, float('nan'), 181):
            try:
                body.step(invalid)
            except ValueError:
                pass
            else:
                raise AssertionError('accepted invalid step')
        assert body.snapshot()['time_s'] == initial['time_s']
        body.meal(Meal(carbohydrate_g=10, protein_g=2, fat_g=1, water_ml=20))
        try:
            body.meal(Meal(carbohydrate_g=20))
        except ValueError:
            pass
        else:
            raise AssertionError('pending native meal was overwritten')
        try:
            body.save_state()
        except ValueError:
            pass
        else:
            raise AssertionError('unreloadable pending native action was serialized')
        # Native actions are consumed in PreProcess, not a Python shadow stomach.
        after = body.step(.02)
        assert abs(after['time_s']-initial['time_s']-.02) < 1e-7
        assert 9.9 < after['values']['stomach_carbohydrate_g']-initial['values']['stomach_carbohydrate_g'] <= 10
        body.step(.02)
        assert body.snapshot()['values']['stomach_carbohydrate_g'] <= after['values']['stomach_carbohydrate_g']
        body.apnea(1)
        apneic = [body.step(1) for _ in range(35)]
        body.apnea(0)
        recovery = [body.step(1) for _ in range(40)]
        a = [r['values']['lung_volume_ml'] for r in apneic[-10:]]
        b = [r['values']['lung_volume_ml'] for r in recovery[-15:]]
        assert max(b)-min(b) > 50, 'breathing did not resume'
        assert max(a)-min(a) < max(b)-min(b), 'apnea did not suppress excursion'
        saved = body.save_state()
        assert saved.is_file()
        final = body.snapshot()
    report = {'passed': True, 'output': str(out), 'initial_time_s': initial['time_s'],
              'final_time_s': final['time_s'], 'apnea_excursion_ml': max(a)-min(a),
              'recovery_excursion_ml': max(b)-min(b), 'ports': len(final['values']), 'batch_parity_errors': errors,
              'prestart_detached_resources_and_copied_state_selected': True,
              'exercise_minus_rest_120s': demand_delta}
    (out/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
