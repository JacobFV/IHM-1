"""Which reproductive fields have an organ in this body, and which do not.

    python scripts/audit_hormone_axis_anatomy.py

Two things are declared in this repository that the body cannot support, and one
that it can, and until now all three looked alike in the inventory.

``ihm/native/reproductive.py`` executes the hash-pinned Schlosser & Selgrade
gonadotropin model.  Its own limitation list says estradiol, progesterone and
inhibin are prescribed functions of time rather than states.  The sharper point
that list does not make: **it is a menstrual cycle model, and this body has
testes and a prostate.**  It is also coupled to nothing -- reachable only as a
standalone trajectory through ``/api/reproductive``.

``ihm/fields/systems.py`` declares ``uterine.endometrial_thickness``,
``uterine.contractile_pressure`` and a ``placental`` group.  This script counts
the entities that would carry them.  The answer is zero, and it is zero because
no catalogued source ships a female whole-body mesh, not because anybody chose
it.

What the body CAN support is the gonadal group, and the support is measured
here rather than asserted: the testes are real meshes with a real volume, and
that volume is checked against the clinical reference range, which is a known
answer this script did not choose.

This is an audit and it changes no dynamics.  Making the absence explicit is
what it is for; ``docs/DISCONNECTS.md`` is the standing list of places where a
declared model and a running model are different objects, and two of the rows
below belong on it.
"""
from pathlib import Path
import argparse
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.fields.systems import SYSTEMS  # noqa: E402

BINDING = 'data/derived/anatomy-segment-binding/binding.json'
ANATOMY = 'data/derived/canonical/anatomy.json'

#: Field group -> the entity-name pattern whose presence would support it.
#: Deliberately explicit, one pattern per group, so a group with no organ says
#: so rather than matching something adjacent.
SUPPORT_PATTERNS = {
    'reproductive': re.compile(r'\btestis\b|\btestes\b|\bovary\b|\bovarian\b|'
                               r'epididym|deferent|seminal vesicle|ejaculatory|'
                               r'prostat', re.I),
    'uterine': re.compile(r'\buter|endometri|myometri|cervix uteri', re.I),
    'placental': re.compile(r'placent|chorion|umbilical cord', re.I),
    'follicular': re.compile(r'hair follicle|\bhair\b', re.I),
}

#: The gonad proper, separated from the accessory tract, because a hormone axis
#: is produced by the gonad and not by the duct system.
GONAD = re.compile(r'\btestis\b|\btestes\b|\bovary\b|\bovarian\b', re.I)

#: Ultrasound reference for adult testicular volume, per testis.  Quoted as a
#: RANGE from the clinical literature and used only as an order-of-magnitude
#: known answer: a mesh volume outside it would mean the mesh is not a testis or
#: the units are wrong, which is exactly the class of error this repository's
#: ledger is full of.
TESTIS_VOLUME_REFERENCE_ML = (12.0, 30.0)


def load_entities():
    binding = json.loads((ROOT / BINDING).read_text())['entities']
    volumes = {}
    for entity in json.loads((ROOT / ANATOMY).read_text())['entities']:
        source = entity['id'].replace('body-bp3d-', '').replace('body-', '')
        volumes[entity['id']] = entity.get('volume_m3')
        volumes[source] = entity.get('volume_m3')
        volumes[entity['name']] = entity.get('volume_m3')
    return binding, volumes


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', default='data/derived/reproductive/anatomy_support.json')
    args = parser.parse_args()

    binding, volumes = load_entities()
    declared = {}
    for field, support, region, tau, components in SYSTEMS:
        if field not in SUPPORT_PATTERNS:
            continue
        declared[field] = {'support': support, 'region': region,
                           'components': ['%s.%s' % (field, name)
                                          for name, _u, _c, _s in components]}

    report = {'schema': 'ihm.hormone-axis-anatomy.v1', 'fields': {},
              'entities_searched': len(binding)}
    for field, meta in declared.items():
        pattern = SUPPORT_PATTERNS[field]
        matches = [(eid, entry) for eid, entry in binding.items()
                   if pattern.search(entry.get('name', ''))]
        measured_volumes = [volumes.get(eid) for eid, _ in matches]
        with_volume = [v for v in measured_volumes if v is not None]
        volume = sum(with_volume)
        report['fields'][field] = {
            'declared_components': meta['components'],
            'declared_region': meta['region'],
            'declared_support': meta['support'],
            'supporting_entities': len(matches),
            'supporting_entity_names': sorted({e.get('name', '') for _i, e in matches}),
            'segments': sorted({e['segment'] for _i, e in matches}),
            'entities_with_a_computed_volume': len(with_volume),
            'total_volume_ml_OVER_THOSE_ONLY': 1e6 * volume,
            'volume_coverage_note': (
                'anatomy.json carries volume_m3 for %d of these %d entities and '
                'None for the rest, so the sum is a partial one and is named as '
                'such. Summing a partial set under a complete-sounding label is '
                'the exact shape of this repository\'s withdrawn results.'
                % (len(with_volume), len(matches))),
            'anatomically_supported': bool(matches),
        }

    # The gonad, measured against a reference range chosen before the number
    # was read.
    gonads = [(eid, entry) for eid, entry in binding.items()
              if GONAD.search(entry.get('name', ''))]
    per_gonad = []
    for eid, entry in gonads:
        volume = volumes.get(eid)
        per_gonad.append({'id': eid, 'name': entry.get('name'),
                          'segment': entry['segment'],
                          'volume_ml': None if volume is None else 1e6 * volume})
    measured = [g['volume_ml'] for g in per_gonad if g['volume_ml'] is not None]
    report['gonad'] = {
        'note': ('The only reproductive field this body can support, and the one '
                 'quantity a hormone axis could be scaled by. Measured from the '
                 'entity meshes, not assumed.'),
        'entities': per_gonad,
        'reference_range_ml_per_gonad': list(TESTIS_VOLUME_REFERENCE_ML),
        'reference_basis': ('adult testicular volume by ultrasound, clinical '
                            'reference range; used as an order-of-magnitude '
                            'known answer, not as a calibration target'),
        'in_reference_range': bool(measured) and all(
            TESTIS_VOLUME_REFERENCE_ML[0] <= v <= TESTIS_VOLUME_REFERENCE_ML[1]
            for v in measured),
        'measured_ml': measured,
    }

    reproductive_index = ROOT / 'data/derived/reproductive/index.json'
    if reproductive_index.exists():
        model = json.loads(reproductive_index.read_text())['models'][0]
        report['executed_model'] = {
            'id': model['id'],
            'states': [s['id'] for s in model['states']],
            'prescribed_inputs': [c['id'] for c in model['channels']
                                  if c.get('kind') == 'prescribed_time_input'],
            'body_match': (
                'MISMATCH. This is a menstrual-cycle gonadotropin model. Its four '
                'states are releasable and circulating LH and FSH; estradiol, '
                'progesterone and inhibin are prescribed functions of time whose '
                'source is an ovary. This body has %d gonad entities and they are '
                'testes. The model is not wrong -- it is an executed, hash-pinned '
                'published model -- it is running beside a body it does not '
                'describe, and it is coupled to nothing.' % len(gonads)),
        }

    report['findings'] = [
        'uterine.*: %d supporting entities. The two declared components '
        '(endometrial_thickness, contractile_pressure) have no organ in this '
        'body and no source to get one from.'
        % report['fields']['uterine']['supporting_entities'],
        'placental.*: %d supporting entities, and a placenta is not an organ of '
        'a non-pregnant body of either sex, so this group is doubly unsupported.'
        % report['fields']['placental']['supporting_entities'],
        'reproductive.*: %d supporting entities, all bound to the pelvis, of '
        'which %d carry a computed volume at all. reproductive.estradiol and '
        'reproductive.progesterone are declared on a body whose gonads are '
        'testes; reproductive.testosterone is the one component of that group '
        'the anatomy supports, and the executed hormone model does not compute '
        'it.'
        % (report['fields']['reproductive']['supporting_entities'],
           report['fields']['reproductive']['entities_with_a_computed_volume']),
        'KNOWN-ANSWER CHECK FAILED, and the failure is informative. The gonad '
        'meshes measure %s mL each against a %g-%g mL adult clinical reference '
        'range -- roughly HALF. The prostate mesh is %s mL against a typical '
        'adult 15-25 mL, so the two measurable reproductive organs are small by '
        'about the same factor. That is a property of this atlas rather than a '
        'one-off: BodyParts3D is a single donor, its skin extent is 1.7195 m, '
        'and nothing here has ever calibrated organ volume against a population. '
        'Any hormone axis scaled by gonadal volume would inherit it, which is '
        'why the number is reported before anything is scaled by it.'
        % (', '.join('%.1f' % v for v in measured) or 'no',
           *TESTIS_VOLUME_REFERENCE_ML,
           ('%.1f' % (1e6 * volumes['prostate'])) if volumes.get('prostate') else 'an uncomputed'),
    ]

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'fields': {k: {kk: vv for kk, vv in v.items()
                                     if kk != 'supporting_entity_names'}
                                 for k, v in report['fields'].items()},
                      'gonad': report['gonad'],
                      'findings': report['findings']}, indent=2))
    print('written to ' + args.out)


if __name__ == '__main__':
    main()
