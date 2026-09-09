"""Sustained real native fall/contact check; no balance controller or tether."""
from pathlib import Path
import argparse
import json
import sys
import tempfile
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream


def verify(duration_s=5.0, augmented_registration=None):
    output = Path(tempfile.mkdtemp(prefix='native-fall-contact-', dir=ROOT/'data/derived'))
    mass = sum(e['mass_kg'] for e in json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())['entities'])
    started = time.monotonic()
    plant = NativeMechanicalStream(ROOT, output/'native', environment='upright', target_mass_kg=mass,
                                   augmented_registration=augmented_registration)
    samples = []
    try:
        initial = plant.snapshot()
        assert all(abs(c['force_n'][1]) < 1e-8 for c in initial['contacts'] if c['name'].startswith('fall_support_')), 'Body proxies must not overlap the floor in the source standing pose'
        for index in range(round(duration_s/.02)):
            state = plant.advance(.02)
            samples.append({'time_s': state['time_s'], 'pelvis_height_m': state['coordinates']['pelvis_ty']['value'],
                            'kinetic_energy_j': state['kinetic_energy_j'],
                            'momentum_residual_n': float(np.linalg.norm(state['momentum_balance_residual_n'])),
                            'body_floor_force_n': sum(c['force_n'][1] for c in state['contacts'] if c['name'].startswith('fall_support_'))})
            if index % 50 == 0:
                print(json.dumps(samples[-1]), flush=True)
        assert state['contact_model'] == 'source_feet_and_inertia_inscribed_body_spheres'
        assert abs(state['time_s']-duration_s) < 1e-8
        assert min(s['pelvis_height_m'] for s in samples) > 0, 'Pelvis passed through floor'
        assert max(s['body_floor_force_n'] for s in samples) > mass*9.81*.2, 'No non-foot floor reaction'
        assert max(s['momentum_residual_n'] for s in samples) < 1e-5, 'Unbalanced external/contact force'
        report = {'passed': True, 'duration_s': duration_s, 'wall_s': time.monotonic()-started,
                  'scope': 'Uncontrolled fall caught by uncalibrated inertia-derived body contact spheres; not standing or walking',
                  'target_mass_kg': mass, 'augmented_registration': augmented_registration, 'samples': samples}
        (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k!='samples'} | {'output': str(output)}, indent=2))
        return report
    except BaseException as error:
        (output/'failure.json').write_text(json.dumps({'passed': False, 'error': str(error),
            'requested_duration_s': duration_s, 'wall_s': time.monotonic()-started,
            'augmented_registration': augmented_registration, 'samples': samples}, indent=2)+'\n')
        raise
    finally:
        plant.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--duration', type=float, default=5.0)
    parser.add_argument('--augmented-registration')
    args = parser.parse_args()
    if args.duration < 2 or abs(args.duration/.02-round(args.duration/.02)) > 1e-8:
        parser.error('duration must be at least 2 s and an integer multiple of .02 s')
    verify(args.duration, args.augmented_registration)
