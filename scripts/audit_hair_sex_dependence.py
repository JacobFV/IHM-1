"""What a sex-dependent hair distribution would actually need, measured.

    python scripts/audit_hair_sex_dependence.py

``docs/BODY_PARAMETERS.md`` records that seven of the 21 hair fields carry an
``androgen_dependent`` label, that the label is a string in a metadata tuple,
and that there is no androgen variable to drive it.  All true, and it points at
the wrong missing piece.

Read the built field evidence and the picture is sharper.  Every field already
carries the parameter that sexual dimorphism in body hair actually runs through:
``shaft_bearing_fraction``, the fraction of counted follicles that produce a
visible shaft.  The evidence file's own note says why -- follicle NUMBER is
roughly sex-invariant (seago1985 found no sex difference in follicle density on
thigh or upper arm); what differs is how many of those follicles make a terminal
hair.  So a female hair distribution is not a hormone signal wired into a field.
It is one number per field, and this script reports the tier of that number as
the evidence file itself declares it.

The finding this produces is a labelling one and it matters for the ``sex``
parameter: six of the seven androgen-dependent fields set that fraction to
``1.0`` at tier ``assumed``, and 1.0 is specifically an ADULT MALE assumption.
It is not marked as one.  A female-proportioned body built today would inherit
it silently, which is exactly why ``sex`` is declared to reach proportions and
not hair.
"""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

FIELDS = 'data/derived/hair-field-candidate-v1/fields.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fields', default=FIELDS)
    parser.add_argument('--out', default='data/derived/hair-field-candidate-v1/sex_dependence.json')
    args = parser.parse_args()

    path = ROOT / args.fields
    if not path.exists():
        raise SystemExit('%s does not exist; build the hair fields first' % args.fields)
    fields = json.loads(path.read_text())['fields']

    rows = []
    for field in fields:
        morphology = field.get('morphology') or {}
        fraction = morphology.get('shaft_bearing_fraction') or {}
        density = field.get('follicle_density_cm2') or {}
        rows.append({
            'id': field['id'],
            'hair_kind': field['hair_kind'],
            'androgen_dependent': 'androgen' in field['hair_kind'],
            'follicle_density_cm2': density.get('value'),
            'follicle_density_tier': density.get('tier'),
            'shaft_bearing_fraction': fraction.get('value'),
            'shaft_bearing_fraction_tier': fraction.get('tier'),
            'shaft_bearing_fraction_basis': fraction.get('basis'),
            'footprint_cm2': (field.get('footprint') or {}).get('area_cm2'),
        })

    androgen = [r for r in rows if r['androgen_dependent']]
    assumed_one = [r for r in androgen if r['shaft_bearing_fraction'] == 1.0
                   and r['shaft_bearing_fraction_tier'] == 'assumed']
    unmeasured = [r for r in androgen if r['shaft_bearing_fraction'] is None]
    measured = [r for r in androgen if r['shaft_bearing_fraction_tier'] == 'measured']

    report = {
        'schema': 'ihm.hair-sex-dependence.v1',
        'fields': len(rows),
        'androgen_dependent_fields': len(androgen),
        'rows': rows,
        'shaft_bearing_fraction_over_androgen_fields': {
            'assumed_1.0': [r['id'] for r in assumed_one],
            'declared_unmeasured': [r['id'] for r in unmeasured],
            'measured': [r['id'] for r in measured],
        },
        'findings': [
            'The parameter a sex-dependent hair distribution runs through already '
            'exists on every field: shaft_bearing_fraction. What is missing is a '
            'measured value for it, not a hormone variable and not a coupling. '
            'Wiring an androgen signal into these fields today would drive a '
            'quantity whose male value is itself an assumption.',
            'Of the %d androgen-dependent fields, %d set shaft_bearing_fraction '
            'to 1.0 at tier "assumed" (%s), %d declare it unmeasured (%s), and %d '
            'have a measured value. A fraction of 1.0 -- every counted follicle '
            'bears a shaft -- is an ADULT MALE assumption and is not labelled as '
            'one anywhere in the field record.'
            % (len(androgen), len(assumed_one), ', '.join(r['id'] for r in assumed_one),
               len(unmeasured), ', '.join(r['id'] for r in unmeasured) or 'none',
               len(measured)),
            'The one field that refuses the assumption is beard, and the refusal '
            'is correct rather than inconsistent: its density is transferred from '
            'FEMALE CHEEK VELLUS hair, where assuming every follicle bears a '
            'terminal shaft would be badly wrong. The other six are adult male '
            'axillary, pubic, perineal, chest, abdomen and back, where it is '
            'nearly right -- for a man.',
            'Follicle density itself is largely sex-invariant, which is what makes '
            'this the right place to look: seago1985 measured no sex difference on '
            'thigh or upper arm, and the field records cite it.',
        ],
        'what_would_close_it': [
            'A measured female shaft_bearing_fraction for beard, axillary, pubic, '
            'perineal, chest, abdomen and back. Seven numbers. Ferriman-Gallwey '
            'scoring is the standard female body-hair instrument and is ORDINAL, '
            'so it does not supply a fraction directly.',
            'Until those exist, sex reaching hair would move a rendered surface '
            'and no measured quantity, which is the state the sex parameter was '
            'raising over in the first place.',
        ],
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'rows'}, indent=2))
    print('\n%-14s %-28s %-10s %-10s' % ('field', 'kind', 'density', 'shaft frac'))
    for row in rows:
        print('%-14s %-28s %-10s %-10s%s'
              % (row['id'], row['hair_kind'],
                 row['follicle_density_cm2'], row['shaft_bearing_fraction'],
                 '  <- androgen-dependent' if row['androgen_dependent'] else ''))
    print('\nwritten to ' + args.out)


if __name__ == '__main__':
    main()
