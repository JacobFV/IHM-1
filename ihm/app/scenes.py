"""Environment tiles read this: identity, label, slot, thumbnail and dependencies.

Slot carries the exclusivity model, exactly as the garment catalogue does. Entries
sharing a slot are mutually exclusive; different slots combine. `requires` names the
slot and ids an entry depends on, so a mattress cannot be offered without its bed.

The live ENVIRONMENTS dict stays authoritative for gravity, axis, plane and supports.
A derived record is attached only where it agrees with the live values field for field;
a disagreeing record is dropped and reported, never served stale.
"""
import hashlib
import json
from pathlib import Path

CATALOGUE = 'data/derived/environment-catalogue-v1/catalogue.json'
MANIFEST = 'data/derived/environment-catalogue-v1/manifest.json'
PHYSICS_FIELDS = ('gravity', 'axis', 'plane', 'supports')


def _derived(root):
    path = root / CATALOGUE
    if not path.is_file():
        return None
    return json.loads(path.read_bytes())


def _thumbnail(root, entry):
    relative = entry.get('thumbnail')
    if not relative:
        return None
    file = (root / relative).resolve()
    if not file.is_relative_to(root) or not file.is_file():
        return None
    return str(file.relative_to(root))


def scene_catalog(root, live):
    """`live` is SceneSessions.catalog(): the running engine's own environment values."""
    root = Path(root).resolve()
    catalog = dict(live)
    environments = []
    notes = []
    derived = None
    try:
        derived = _derived(root)
    except (OSError, ValueError) as error:
        notes.append('Derived environment catalogue unreadable: ' + str(error))
    records = {} if derived is None else {r['id']: r for r in derived['environments']}
    for entry in live['environments']:
        record = records.get(entry['id'])
        if record is not None and any(record.get(f) != entry.get(f) for f in PHYSICS_FIELDS):
            notes.append('Derived record for ' + entry['id'] + ' disagrees with the running engine; rebuild scripts/build_environment_catalogue.py')
            record = None
        merged = dict(entry, slot='environment', kind='environment', thumbnail=None, thumbnail_url=None)
        if record is not None:
            thumbnail = _thumbnail(root, record)
            if thumbnail is None and record.get('thumbnail'):
                notes.append('Thumbnail missing for ' + entry['id'])
            merged = {**record, **{k: v for k, v in entry.items() if k in PHYSICS_FIELDS or k == 'id'},
                      'engine_description': entry.get('description'),
                      'slot': 'environment', 'kind': 'environment', 'thumbnail': thumbnail,
                      'thumbnail_url': None if thumbnail is None else '/api/scene/thumbnail/' + record['id']}
        environments.append(merged)
    catalog['environments'] = environments
    catalog['schema'] = 'ihm.environment-catalog.v1'
    components = []
    if derived is not None:
        for record in derived['components']:
            thumbnail = _thumbnail(root, record)
            if thumbnail is None and record.get('thumbnail'):
                notes.append('Thumbnail missing for ' + record['id'])
            components.append({**record, 'kind': 'component', 'thumbnail': thumbnail,
                               'thumbnail_url': None if thumbnail is None else '/api/scene/thumbnail/' + record['id']})
        catalog['slots'] = derived['slots']
        catalog['exclusivity_model'] = derived['exclusivity_model']
        catalog['camera'] = derived['camera']
        catalog['verified'] = derived['verified']
        catalog['not_selectable'] = derived['not_selectable']
        catalog['catalogue_commit'] = derived['commit']
        manifest = root / MANIFEST
        catalog['catalogue_manifest_sha256'] = hashlib.sha256(manifest.read_bytes()).hexdigest() if manifest.is_file() else None
    else:
        catalog['slots'] = [{'id': 'environment', 'label': 'Environment', 'exclusive': True,
                             'required': True, 'default': None, 'requires': [],
                             'note': 'The body accepts exactly one environment.'}]
        catalog['exclusivity_model'] = 'Slots carry exclusivity: entries sharing a slot are mutually exclusive, different slots combine.'
        notes.append('Derived environment catalogue absent; serving engine environments only. Run scripts/build_environment_catalogue.py')
    catalog['components'] = components
    # One flat array for a tile grid; slot and requires carry the rules, so the
    # front end reads them rather than encoding which tile excludes which.
    catalog['tiles'] = environments + components
    catalog['notes'] = notes
    return catalog


def thumbnail(root, ident, live):
    root = Path(root).resolve()
    catalog = scene_catalog(root, live)
    entry = next((e for e in catalog['environments'] + catalog['components'] if e['id'] == ident), None)
    if entry is None or not entry.get('thumbnail'):
        return None
    file = (root / entry['thumbnail']).resolve()
    if not file.is_relative_to(root) or not file.is_file():
        return None
    return file
