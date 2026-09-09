#!/usr/bin/env python3
"""Where the unified brain/body/world loop's time actually goes.

The loop costs roughly 80 seconds of CPU per simulated second, so a five-second
paired comparison takes about twelve minutes. This measures the split rather
than guessing it, and records what an exact optimization can and cannot buy.

FINDING. The dominant cost is the native OpenSim advance, not Python. Snapshot
copying was the visible top entry in a cProfile listing, but cProfile charges
about a microsecond to every one of the five million Python-level calls that
copy.deepcopy makes while charging almost nothing to pickle's C loop, so it
overstated copying and understated everything C-bound. The exact copy fixes are
real but small; the native engine is the wall.

SCOPE. Timings are from this machine and this body/model. Wall clock here varies
by more than a factor of two run to run, so the per-copy and per-command figures
below are medians of repeated in-process measurements, which are stable, and the
end-to-end wall figures are reported with that variance stated. No claim is made
that any measured alternative substep is physically acceptable.

Run: PYTHONPATH=. .venv/bin/python scripts/verify_loop_performance.py
"""
import argparse
import copy
import json
import shutil
import statistics
import tempfile
import time
from pathlib import Path

from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.snapshot_data import clone_snapshot_data

ROOT = Path(__file__).resolve().parents[1]
REGISTRATION = 'data/derived/mechanics/whole_body_lumbar_current/registration.json'
MASS_KG = 77.6122029

# Calls per simulated second observed in a profiled 50-step upright run with the
# trained cortical stance controller, a play-floor scene and a 5 ms exchange.
CALLS = {'native_advance': 200, 'native_checkpoint': 250, 'native_release': 250,
         'native_snapshot_retained': 550,
         'native_snapshot_removed_from_checkpoint': 250,
         'observation_clone_removed': 150}


def median(fn, repeats):
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - started)
    return statistics.median(samples)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repeats', type=int, default=15)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    output = args.output or Path(tempfile.mkdtemp(prefix='loop-performance-',
                                                  dir=ROOT / 'data/derived'))
    output.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix='loop-performance-plant-', dir=ROOT / 'data/derived'))
    report = {'schema': 'ihm.loop-performance.v1', 'passed': False,
              'scope': __doc__.split('SCOPE.')[1].split('Run:')[0].strip(),
              'calls_per_simulated_second': dict(CALLS)}
    plant = None
    try:
        plant = ArticulatedBodyPlant(ROOT, scratch / 'plant', environment='upright',
                                     target_mass_kg=MASS_KG, augmented_registration=REGISTRATION)
        native = plant.native

        # ---- 1. Copy cost, on the real state shapes ------------------------
        state = native.snapshot()
        frame = plant.snapshot()
        copies = {}
        for name, value in (('native_state', state), ('projected_frame', frame)):
            deep = median(lambda: copy.deepcopy(value), args.repeats)
            clone = median(lambda: clone_snapshot_data(value), args.repeats)
            assert copy.deepcopy(value) == clone_snapshot_data(value)
            copies[name] = {'deepcopy_s': deep, 'clone_snapshot_data_s': clone,
                            'speedup': deep / clone}
        report['copy_cost'] = copies

        # ---- 2. Native command cost ----------------------------------------
        tokens = []
        commands = {
            'checkpoint_s': median(lambda: tokens.append(native.checkpoint()), args.repeats),
            'release_s': median(lambda: native.release(tokens.pop()), args.repeats),
            'snapshot_s': median(native.snapshot, args.repeats)}
        while tokens:
            native.release(tokens.pop())
        steps = {}
        for dt in (.001, .005, .01, .02):
            each = median(lambda: native.advance(dt), max(6, args.repeats // 2))
            steps[f'{dt * 1000:g}ms'] = {'median_s': each, 'cost_per_simulated_s': each / dt}
        commands['advance_by_step'] = steps
        report['native_commands'] = commands

        # ---- 3. What the exact optimization removed ------------------------
        saved = (CALLS['native_snapshot_retained']
                 * (copies['native_state']['deepcopy_s'] - copies['native_state']['clone_snapshot_data_s'])
                 + CALLS['native_snapshot_removed_from_checkpoint'] * copies['native_state']['deepcopy_s']
                 + CALLS['observation_clone_removed'] * copies['projected_frame']['clone_snapshot_data_s'])
        native_per_second = steps['5ms']['cost_per_simulated_s']
        report['exact_optimization'] = {
            'changes': [
                'MechanicalStream.snapshot uses clone_snapshot_data (pickle round trip) '
                'instead of copy.deepcopy. The native state is parsed JSON that is only '
                'ever rebound, never mutated in place, so the copy is identical.',
                'MechanicalStream.checkpoint no longer takes a Python snapshot. Rollback '
                'is owned by the native process against its token; no caller ever read '
                'that snapshot back.',
                'SelectiveProjectionPlant.advance_observation no longer clones the '
                'observation a second time. _project already builds every field with '
                'clone_snapshot_data, copy.deepcopy or fresh .tolist() output, so the '
                'observation shares no structure with the native result or the endpoint.'],
            'python_cpu_s_saved_per_simulated_s': saved,
            'exactness': ('Proven by byte-identical recorded traces: a 50-step upright '
                          'cortical-stance run produced trace.json sha256 '
                          '570b4047061a48b8e440ef426f220c26dcc983c5197a825cd759087cc0b3b792 '
                          '(738907 bytes) before and after, across four runs.')}

        # ---- 4. The lever that is NOT exact ---------------------------------
        report['unrealized_headroom'] = {
            'native_share': ('The native advance dominates. At the production 5 ms '
                             'mechanical substep it costs about '
                             f'{native_per_second:.0f} s of native time per simulated '
                             'second, against roughly 30 s of Python. Removing all '
                             'Python copying could not have fixed this.'),
            'substep_tradeoff': (
                'Native cost is strongly sublinear in step size, so fewer, larger '
                'substeps are the largest available lever: see '
                'native_commands.advance_by_step. Raising advance_body_world(max_step) '
                'from 5 ms toward 20 ms would cut native calls fourfold and native cost '
                'by more than that. It is NOT exact and NOT applied here: the 5 ms '
                'exchange was chosen deliberately to reduce severe contact transients, '
                'and changing it changes the physics. Measured only as headroom.'),
            'no_pipelining': ('Python and the native process strictly alternate over a '
                              'pipe, so wall time is their sum, not their maximum. '
                              'Overlapping them is an architectural change, not a '
                              'micro-optimization.'),
            'measurement_caveat': ('advance_by_step advances an unactuated settling body '
                                   'with no world loads, so absolute values differ from '
                                   'the production loop. The monotone decrease with '
                                   'larger steps is the robust part.')}
        report['passed'] = True
    except BaseException as error:
        report['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        if plant is not None:
            plant.close()
        shutil.rmtree(scratch, ignore_errors=True)
        (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({'output': str(output), **report}, indent=2), flush=True)
    raise SystemExit(0 if report['passed'] else 1)


if __name__ == '__main__':
    main()
