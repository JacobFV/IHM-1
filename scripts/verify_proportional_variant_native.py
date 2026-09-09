"""Does the ENGINE accept a re-proportioned body, or only the XML gates?

    python scripts/verify_proportional_variant_native.py \
        --variant data/derived/mechanics/proportional_female/registration.json

Every number in `docs/BODY_PARAMETERS.md` section 3 is measured on XML and on
OpenSim's own path evaluator.  That is the right place to measure most of them,
and it is not the same claim as "the plant loaded it and integrated".  Stature
variants have carried `native_acceptance_complete: false` since they existed.

This runs the variant and the base body through `NativeMechanicalStream` under
the same settings and prints both, because one body's numbers alone say nothing.
It is deliberately short: load, read the state, take a few 2 ms steps, close.  It
is an ACCEPTANCE check -- the plant instantiated this geometry with these muscle
paths and integrated without the solver collapsing -- and it is not a claim about
behaviour.  Nothing here walks, stands, or is driven by anything.

Two things it would catch that no XML gate can.  A refitted path set whose
polynomials are evaluated outside the coordinate box they were fitted on shows up
as a musculotendon at an absurd normalised fibre length.  And a body whose
inertia was scaled inconsistently with its mass shows up as an integration that
will not advance in bounded time -- `advance()` is an error-controlled Simbody
step whose cost is not bounded by `dt`, and this repository has watched it go
from 0.112 s to 24.3 s per 10 ms step for exactly that kind of reason.
"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402

#: Millard's active force-length curve is defined on roughly [0.4, 1.6] normalised
#: fibre lengths; outside it the muscle is in the flat tails and the force is a
#: modelling artefact rather than a measurement. Wider than the 0.2-20
#: path-over-fibre corpus band and for a different purpose.
FIBRE_BAND = (0.3, 1.7)
STEPS = 5
DT_S = 0.002


def run_one(label, registration, mass_kg, environment):
    out = ROOT / ('data/derived/.native-acceptance-%s' % label)
    shutil.rmtree(out, ignore_errors=True)
    started = time.time()
    stream = NativeMechanicalStream(
        ROOT, out, environment=environment, target_mass_kg=mass_kg,
        augmented_registration=str(Path(registration).relative_to(ROOT)))
    record = {'label': label, 'registration': str(Path(registration).relative_to(ROOT)),
              'construct_seconds': round(time.time() - started, 3)}
    try:
        state = stream.state
        muscles = state['muscles']
        ratios = {name: m['fiber_length_m'] / m['optimal_fiber_length_m']
                  for name, m in muscles.items()}
        record.update({
            'muscles': len(muscles),
            'coordinates': len(state['coordinates']),
            'bodies': len(state['bodies']),
            'plant_mass_kg': sum(b['mass_kg'] for b in state['bodies'].values()),
            'normalised_fibre_length_min': min(ratios.values()),
            'normalised_fibre_length_max': max(ratios.values()),
            'worst_short': min(ratios, key=ratios.get),
            'worst_long': max(ratios, key=ratios.get),
            'fibre_band': list(FIBRE_BAND),
            'in_fibre_band': all(FIBRE_BAND[0] <= r <= FIBRE_BAND[1]
                                 for r in ratios.values()),
        })
        step_times = []
        for _ in range(STEPS):
            started = time.time()
            stream.advance(DT_S)
            step_times.append(time.time() - started)
        record.update({
            'steps': STEPS, 'dt_s': DT_S,
            'simulated_time_s': stream.state['time_s'],
            'max_step_seconds': max(step_times),
            'mean_step_seconds': sum(step_times) / len(step_times),
            'advanced': abs(stream.state['time_s'] - STEPS * DT_S) < 1e-9,
        })
    finally:
        stream.close()
        shutil.rmtree(out, ignore_errors=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--variant', required=True)
    parser.add_argument('--base', default='data/derived/mechanics/whole_body_lumbar_current/registration.json')
    parser.add_argument('--environment', default='supine')
    parser.add_argument('--out', default=None)
    args = parser.parse_args()

    variant_record = json.loads((ROOT / args.variant).read_text())
    mass = variant_record['body_parameters']['parameters']['mass_kg']

    rows = [run_one('variant', ROOT / args.variant, mass, args.environment)]
    base_path = ROOT / args.base
    if base_path.exists():
        # Same requested mass, so the comparison is of GEOMETRY and muscle paths
        # and not of how heavy the two bodies were asked to be.
        rows.append(run_one('base', base_path, mass, args.environment))

    summary = {
        'schema': 'ihm.proportional-variant-native-acceptance.v1',
        'note': ('acceptance only: the plant instantiated this geometry with these '
                 'refitted muscle paths and integrated. Not behaviour, not '
                 'locomotion, not a result about the body.'),
        'runs': rows,
        'passed': all(r['advanced'] and r['in_fibre_band'] for r in rows),
    }
    if len(rows) == 2:
        variant, base = rows
        summary['comparison'] = {
            'note': ('both bodies asked for the same target mass, so what differs '
                     'is geometry and muscle path REPRESENTATION -- the muscle '
                     'count is the same 98 either way. The base runs the shipped '
                     '80-path set, so 18 of its muscles keep their GeometryPath; '
                     'the variant runs a 98-path refit and none do.'),
            'muscles_variant_vs_base': [variant['muscles'], base['muscles']],
            'fibre_band_variant': [variant['normalised_fibre_length_min'],
                                   variant['normalised_fibre_length_max']],
            'fibre_band_base': [base['normalised_fibre_length_min'],
                                base['normalised_fibre_length_max']],
            'max_step_seconds_variant_vs_base': [variant['max_step_seconds'],
                                                 base['max_step_seconds']],
        }
    print(json.dumps(summary, indent=2))
    # A variant that has been accepted by the plant should stop saying it has
    # not been. The registration is written before this runs, so the field it
    # ships with is false and only a passing run may change it.
    if summary['passed']:
        record = json.loads((ROOT / args.variant).read_text())
        record['native_acceptance_complete'] = True
        record['native_acceptance'] = {
            'scope': ('the plant instantiated this geometry with these refitted '
                      'muscle paths and integrated %d steps of %g s. Not '
                      'behaviour, not locomotion, not standing.' % (STEPS, DT_S)),
            'environment': args.environment,
            'target_mass_kg': mass,
            'normalised_fibre_length': [rows[0]['normalised_fibre_length_min'],
                                        rows[0]['normalised_fibre_length_max']],
            'compared_against_base': len(rows) == 2,
        }
        (ROOT / args.variant).write_text(json.dumps(record, indent=2) + '\n')
        print('marked native_acceptance_complete in ' + args.variant)
    if args.out:
        out = ROOT / args.out
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, indent=2) + '\n')
        print('written to ' + args.out)
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
