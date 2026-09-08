"""Materialization fidelity: how much of the implicit model one materialization draws.

The implicit model holds every structure the atlas carries. A materialization
selects a subset of it, and which subset is a choice with a measured cost. The
tiers here are derived from the manifest and from the two colocation audits, at
serve time, so nothing downstream has to carry a structure count or an ID list.
"""
from pathlib import Path
import json

MODEL_ID = 'ihm-body'
COLOCATION = 'data/derived/display-structure-promotion-v1/surface-colocation.json'
NAME_RESOLUTION = 'data/derived/tissue-material-assignment-candidate-v1/colocation.json'

# What each tier actually costs in the shipped app, measured through the app
# rather than estimated from the counts: Playwright against the local server,
# `#scene-status` clearing as the readiness signal and
# performance.memory.usedJSHeapSize read three seconds later. A tier is labelled
# with its own number so the reader chooses knowing the price, and the structure
# count the reading was taken at is recorded so a changed atlas shows up as a
# stale measurement rather than as a quietly wrong label.
MEASURED = {
    'body': {'time_to_ready_s': 28, 'js_heap_mb': 350, 'structures_drawn': 2050, 'measured_of': 2229},
    'body-fidelity-registered': {'time_to_ready_s': 39, 'js_heap_mb': 500, 'structures_drawn': 3706, 'measured_of': 3885},
    'body-fidelity-complete': {'time_to_ready_s': 40, 'js_heap_mb': 510, 'structures_drawn': 3818, 'measured_of': 3997},
}
MEASUREMENT = ('Chrome, software rasteriser, 1440x1000, local server, warm asset cache; time from '
               'navigation to the scene status clearing, and performance.memory.usedJSHeapSize after '
               'a forced collection. A cold cache has measured half again as long.')


def _confirmed_duplicates(root):
    """Promoted structures that naming settled as one structure under two names.

    surface-colocation.json measures 291 promoted structures that a geometry-only
    key would have called aliases of an already-drawn structure, and says in its
    own finding that geometry cannot separate the true duplicates in that set from
    true distinct neighbours. colocation.json applies a name key to exactly that
    list and resolves a subset of it. Only the resolved subset is defensible to
    drop: the rest may be real anatomy standing close to something else.
    """
    root = Path(root)
    resolution = json.loads((root / NAME_RESOLUTION).read_bytes())
    if resolution.get('input') != COLOCATION:
        raise ValueError('Colocation name resolution no longer reads the promotion audit it is keyed to')
    measured = json.loads((root / COLOCATION).read_bytes())
    considered = {row['promoted_id'] for row in measured['best_match_per_structure']}
    duplicates = {row['promoted_id'] for row in resolution['resolved']
                  if row.get('verdict') == 'same_structure_under_two_names'}
    if not duplicates or not duplicates <= considered:
        raise ValueError('Resolved duplicates are not a subset of the measured colocation set; rebuild the audits')
    return duplicates, measured, resolution


def materialization_tiers(root, manifest):
    structures = [s for s in manifest['structures'] if s.get('model_id') == MODEL_ID]
    if not structures:
        raise ValueError('The manifest carries no structures for ' + MODEL_ID)
    everything = [s['id'] for s in structures]
    acquired = {s['id'] for s in structures if s.get('evidence_kind') == 'source_geometry'}
    duplicates, measured, resolution = _confirmed_duplicates(root)
    # An ID that no longer names a structure is a stale audit, not a tier.
    known = set(everything)
    duplicates = {i for i in duplicates if i in known}
    sources = {
        'manifest': 'data/derived/app/manifest.json',
        'colocation': COLOCATION,
        'name_resolution': NAME_RESOLUTION,
    }
    tiers = [
        {
            'value': 'body',
            'label': 'Whole body · acquired atlas',
            'omits': sorted(known - acquired),
            'basis': "manifest evidence_kind == 'source_geometry'",
            'note': (f'{len(acquired):,} surfaces acquired from BodyParts3D 4.0 and nothing fitted on top of '
                     'them. This is the atlas as measured, before the display promotion registered another '
                     f'{len(known) - len(acquired):,} structures from other datasets onto it.'),
        },
        {
            'value': 'body-fidelity-registered',
            'label': 'Whole body · registered anatomy, less confirmed duplicates',
            'omits': sorted(duplicates),
            'basis': ("every structure except the promoted ones a name key resolved as "
                      "same_structure_under_two_names in " + NAME_RESOLUTION),
            'note': (f'The acquired atlas plus every registered structure, less the {len(duplicates):,} promoted '
                     f'surfaces that naming settled as an already-drawn structure under a second name. '
                     f"{resolution['unresolved']:,} more sit at alias bounds on an acquired surface and are kept: "
                     'the audit says geometry cannot tell those from true neighbours, so dropping them would '
                     'drop real anatomy.'),
        },
        {
            'value': 'body-fidelity-complete',
            'label': 'Whole body · complete implicit model',
            'omits': [],
            'basis': 'every structure the manifest carries for this model',
            'note': (f'Every one of the {len(known):,} structures in the implicit model, including the '
                     f'{len(duplicates):,} the audits call duplicates of something already drawn. Nothing is '
                     'held back and nothing is cheap about it.'),
        },
    ]
    for tier in tiers:
        tier['structures'] = len(known) - len(tier['omits'])
        cost = dict(MEASURED.get(tier['value'], {}))
        cost['measurement'] = MEASUREMENT
        tier['cost'] = cost
    return {
        'schema': 'ihm.materialization-fidelity.v1',
        'model_id': MODEL_ID,
        'default': 'body',
        'structures': len(known),
        'sources': sources,
        'colocation_finding': measured['finding'],
        'tiers': tiers,
    }
