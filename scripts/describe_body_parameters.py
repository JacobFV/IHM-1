"""Print the declared body-parameter schema, and resolve a request against it.

    python scripts/describe_body_parameters.py
    python scripts/describe_body_parameters.py --set stature_m=1.90 --set mass_kg=95
    python scripts/describe_body_parameters.py --json

Every run first recomputes the schema's measured constants from the artifacts
they came from, so a stale constant fails here rather than downstream.
"""
from pathlib import Path
import argparse, json, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.body_parameters import (PARAMETERS, MECHANICAL_SOURCE_MASS_KG,   # noqa: E402
                                 MECHANICAL_TARGET_MASS_KG, MECHANICAL_STATURE_M,
                                 ANATOMICAL_STATURE_M, ANATOMICAL_MASS_KG,
                                 ANATOMY_REGISTRATION_SCALE, STATURE_DISAGREEMENT,
                                 MASS_DISAGREEMENT, assert_measured_defaults, resolve)

STATUS_NOTE = {
    'surfaced': 'already a runtime knob; this schema only names and bounds it',
    'implemented': 'made a knob by this work',
    'declared': 'written down, and reaching the running body through nothing',
    'implemented (proportions only; anatomy unchanged)':
        'made a knob by this work, and the parenthesis is the honest half: it '
        'moves measured PROPORTIONS and changes no anatomical entity',
}


def wrap(text, width, indent):
    out, line = [], ''
    for word in text.split():
        if len(line) + len(word) + 1 > width:
            out.append(line)
            line = word
        else:
            line = (line + ' ' + word) if line else word
    out.append(line)
    return ('\n' + ' ' * indent).join(out)


def parse_assignment(text):
    if '=' not in text:
        raise argparse.ArgumentTypeError('expected name=value, got ' + text)
    name, _, value = text.partition('=')
    try:
        return name.strip(), float(value)
    except ValueError:
        return name.strip(), value.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--set', action='append', type=parse_assignment, default=[],
                        metavar='NAME=VALUE', help='a body parameter to request')
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()

    measured = assert_measured_defaults(ROOT)
    resolved = resolve(dict(args.set))

    if args.json:
        print(json.dumps({'measured': measured, 'schema': [
            {k: (list(v) if isinstance(v, tuple) else v) for k, v in p.items()}
            for p in PARAMETERS], 'resolved': resolved}, indent=2))
        return

    print('the two bodies, measured')
    print('  mechanical (OpenSim Rajagopal subject)')
    print('    source model mass      %12.6f kg   (sum of 22 <Body><mass>)'
          % MECHANICAL_SOURCE_MASS_KG)
    print('    mass actually asked for%12.6f kg   (unattributed literal, %d call sites)'
          % (MECHANICAL_TARGET_MASS_KG, 67))
    print('    stature proxy          %12.6f m    (floor markers to Head marker)'
          % MECHANICAL_STATURE_M)
    print('  anatomical (BodyParts3D atlas + BioGears)')
    print('    composed mass          %12.6f kg' % ANATOMICAL_MASS_KG)
    print('    skin extent            %12.6f m' % ANATOMICAL_STATURE_M)
    print('  they are different objects')
    print('    registration scale     %12.3f      (OpenSim geometry -> atlas)'
          % ANATOMY_REGISTRATION_SCALE)
    print('    stature disagreement   %12.2f %%    (%.4f m predicted vs %.4f m measured)'
          % (100 * STATURE_DISAGREEMENT,
             MECHANICAL_STATURE_M * ANATOMY_REGISTRATION_SCALE, ANATOMICAL_STATURE_M))
    print('    mass disagreement      %12.1f %%' % (100 * MASS_DISAGREEMENT))

    print('\nschema  (ihm.body-parameters.v1)')
    for spec in PARAMETERS:
        domain = ('one of %s' % (spec['domain'],)) if spec['kind'] == 'enum' \
            else ('[%g, %g]' % spec['range'])
        print('\n  %-20s %-14s %s' % (spec['name'], spec['unit'], domain))
        print('    status     %s  -- %s' % (spec['status'], STATUS_NOTE[spec['status']]))
        print('    body       %s' % spec['body'])
        print('    default    %s' % (spec['default'],))
        print('    basis      %s' % wrap(spec['basis'], 66, 15))
        print('    range      %s' % wrap(spec['range_basis'], 66, 15))
        consumers = spec.get('consumers') or ('-- nothing --',)
        print('    consumed by %s' % wrap(', '.join(consumers), 65, 16))
        if spec.get('limitation'):
            print('    LIMITATION %s' % wrap(spec['limitation'], 66, 15))

    print('\nresolved request %s' % (resolved['requested'] or '(all defaults)'))
    for name, value in resolved['parameters'].items():
        print('    %-22s %s' % (name, value))
    print('  derived')
    for name, value in resolved['derived'].items():
        if isinstance(value, (int, float)):
            print('    %-22s %.9g' % (name, value))
    factors = resolved['derived'].get('anisotropic_factors') or {}
    if factors:
        print('  anisotropic per-body factors (fx, fy, fz in the body\'s own frame)')
        for body, triple in sorted(factors.items()):
            print('    %-22s %.6f %.6f %.6f' % (body, *triple))
    realisation = resolved['derived'].get('sex_realisation') or {}
    if realisation.get('what_changed'):
        print('  what sex=%s changed' % realisation['sex'])
        for name, row in realisation['what_changed'].items():
            print('    %-22s %-14.9g -> %-14.9g (x%.6f)'
                  % (name, row['from'], row['to'], row['ratio']))
    if realisation:
        print('  what it did NOT change')
        for key, note in realisation['what_did_not_change'].items():
            if isinstance(note, dict):
                note = note['note']
            print('    %-22s %s' % (key, wrap(note, 60, 27)))
        print('  honest summary: %s' % realisation['honest_summary'])
    print('  limitations of this resolution')
    for note in resolved['limitations']:
        print('    - %s' % wrap(note, 70, 6))


if __name__ == '__main__':
    main()
