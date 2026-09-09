#!/usr/bin/env python3
"""Chain movement programs through one continuing native process.

Each stage is run by the retained evaluator in --keep-open mode, so every stage
after the first continues the same native plant instance and the same physical
time. This script sequences stages and retains a chain receipt; it does not
apply any control of its own and never prescribes coordinates or forces.
"""
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVALUATOR = ROOT / 'scripts/evaluate_native_movement_program.py'
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--program', action='append', required=True,
                   help='program json; repeat once per stage, in order')
    p.add_argument('--output', required=True, help='fresh chain output directory')
    p.add_argument('--scope', required=True, help='explicit claim boundary for the receipt')
    p.add_argument('--continue-after-failure', action='store_true')
    a = p.parse_args()

    out = Path(a.output)
    if out.exists():
        raise ValueError('Refusing to overwrite an existing chain output directory')
    out.mkdir(parents=True)
    stages = [str(Path(x)) for x in a.program]
    started = time.monotonic()

    process = subprocess.Popen(
        [sys.executable, str(EVALUATOR), '--keep-open',
         '--program', stages[0], '--output', str(out / 'stage_0')],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=None,
        text=True, bufsize=1, cwd=str(ROOT))

    def wait_ready():
        while True:
            line = process.stdout.readline()
            if not line:
                return False
            if line.strip() == 'READY_FOR_CONTINUATION':
                return True
            print(line.rstrip(), flush=True)

    reports = []
    aborted = None
    try:
        for index, program in enumerate(stages):
            directory = out / f'stage_{index}'
            if index:
                process.stdin.write(json.dumps({'program': program, 'output': str(directory)}) + '\n')
                process.stdin.flush()
            if not wait_ready():
                aborted = f'evaluator exited before stage {index} completed'
                break
            report = json.loads((directory / 'report.json').read_text())
            reports.append({'stage': index, 'program': program, 'directory': str(directory),
                            'program_sha256': sha(program), **report})
            print(json.dumps({'stage': index, 'elapsed_s': report['elapsed_s'],
                              'completed_horizon': report['completed_horizon'],
                              'error': report['error'],
                              'final_per_foot_clearance_m': report['final_per_foot_clearance_m'],
                              'max_clipped_count': report['max_clipped_count']}), flush=True)
            if not report['completed_horizon'] and not a.continue_after_failure:
                aborted = f'stage {index} did not complete its horizon'
                break
    finally:
        try:
            process.stdin.write(json.dumps({'close': True}) + '\n')
            process.stdin.flush()
        except (BrokenPipeError, ValueError):
            pass
        process.wait(timeout=120)

    identities = sorted({r['native_stream_identity'] for r in reports})
    receipt = {
        'schema': 'ihm.native-movement-chain.v1',
        'stages': reports,
        'stage_count': len(reports),
        'requested_stage_count': len(stages),
        'aborted': aborted,
        'all_stages_completed_horizon': bool(reports) and all(r['completed_horizon'] for r in reports) and len(reports) == len(stages),
        'single_native_process': len(identities) == 1,
        'native_stream_identity': identities,
        'chain_final_time_s': reports[-1]['elapsed_s'] if reports else 0.0,
        'driver_sha256': sha(__file__),
        'evaluator_sha256': sha(EVALUATOR),
        'wall_seconds': time.monotonic() - started,
        'walking_demonstrated': False,
        'scope': a.scope,
    }
    (out / 'chain_report.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'stages'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
