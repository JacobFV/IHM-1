"""Environment tiles read this: identity, label, slot, thumbnail and dependencies.

Slot carries the exclusivity model, exactly as the garment catalogue does. Entries
sharing a slot are mutually exclusive; different slots combine. `requires` names the
slot and ids an entry depends on, so a mattress cannot be offered without its bed.

The live ENVIRONMENTS dict stays authoritative for gravity, axis, plane and supports.
A derived record is attached only where it agrees with the live values field for field;
a disagreeing record is dropped and reported, never served stale.
"""
import hashlib
from functools import lru_cache
import json
from pathlib import Path

CATALOGUE = 'data/derived/environment-catalogue-v1/catalogue.json'
MANIFEST = 'data/derived/environment-catalogue-v1/manifest.json'
PHYSICS_FIELDS = ('gravity', 'axis', 'plane', 'supports')


@lru_cache(maxsize=8)
def _initial_cloth(root, skin_mtime_ns, low, high, mass):
    import gzip,numpy as np
    from ihm.assembly.environment_dynamics import SpringMesh,prepare_cloth
    mesh=SpringMesh('blanket','cloth',low,high,mass)
    skin=json.loads(gzip.decompress((Path(root)/'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz').read_bytes()))
    points=np.array(skin['positions']).reshape(-1,3)
    _,indices=np.unique(np.floor(points/.045).astype(int),axis=0,return_index=True)
    prepare_cloth(mesh,points[indices])
    return mesh.frame()


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
    components, objects, scenes = [], [], []
    if derived is not None:
        for key, kind, target in (('components', 'component', components),
                                  ('objects', 'object', objects),
                                  ('scenes', 'scene', scenes)):
            for record in derived.get(key, []):
                thumbnail = _thumbnail(root, record)
                if thumbnail is None and record.get('thumbnail'):
                    notes.append('Thumbnail missing for ' + record['id'])
                entry = {**record, 'kind': kind, 'thumbnail': thumbnail,
                         'thumbnail_url': None if thumbnail is None else '/api/scene/thumbnail/' + record['id']}
                if kind == 'object':
                    geometry = record.get('geometry')
                    present = geometry is not None and (root / geometry).is_file()
                    if not present:
                        notes.append('Object geometry missing for ' + record['id'])
                    entry['geometry_url'] = '/api/scene/object/' + record['id'] if present else None
                target.append(entry)
        catalog['slots'] = derived['slots']
        catalog['surrounds'] = derived.get('surrounds', [])
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
    # Live embodied capabilities supersede historical display-only catalogue claims.
    # Retain original descriptions as provenance, without presenting them as current behavior.
    from ihm.assembly.environment_dynamics import SCOPE, SpringMesh
    catalog['embodied_environment_physics'] = {'schema':'ihm.environment-state.v1','scope':SCOPE,
        'configuration_field':'environment_selection','state_field':'environment_state'}
    for entry in objects:
        ident=entry['id'];entry['legacy_description']=entry['description']
        mode='cloth' if ident=='blanket' else 'soft' if ident=='pillow' else 'rigid' if ident in ('ball-small','ball-large','block') else 'native_support' if ident=='bed-mattress' else 'fixed'
        entry['embodied_physics']=mode
        entry['description']=entry['label']+' · '+{'cloth':'deforming spring cloth with body reaction forces',
            'soft':'compressible spring lattice with body reaction forces','rigid':'movable rigid body with contact and angular motion',
            'native_support':'native mattress support; firmness requires skin quadrature',
            'fixed':'fixed compound contact geometry'}[mode]+'. Engineering contact parameters; not calibrated.'
        if mode in ('cloth','soft'):
            bounds=entry['bounds_m'];mesh=SpringMesh(ident,mode,bounds['min'],bounds['max'],entry['mass_kg'])
            entry['initial_mesh']=mesh.frame()
            if mode=='cloth':
                skin=root/'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz'
                entry['initial_mesh']=_initial_cloth(str(root),skin.stat().st_mtime_ns,tuple(bounds['min']),tuple(bounds['max']),entry['mass_kg'])
        # Finite deformable objects are supplied in bed scenes. Fixed furniture
        # stays placed by the scene; free balls and blocks are repeatable inserts.
        if mode=='rigid':entry['insertable']=True;entry['insert_label']='Add '+entry['label']
    for entry in scenes:
        entry['legacy_physics']=entry.get('physics')
        entry['physics']=SCOPE;entry['world']['physics']=SCOPE
        entry['world']['light']['kind']='Soft shadowed directional daylight with hemispherical fill'
    for entry in components:
        if entry['slot']=='ambient_thermal':
            entry['live_supported']=entry['id']=='ambient-22c'
            entry['description'] += ' Live embodied adapter: temperature changes are not yet supported.'
    for slot in catalog['slots']:
        if slot['id']=='objects':
            slot['note']='Movable balls and blocks participate in server-owned body contact. Furniture and deformables are placed by scenes.'
            slot['options']=[{'value':e['id'],'label':e['label'],'thumbnail_url':e['thumbnail_url'],'requires':e.get('requires',[])} for e in objects if e.get('insertable')]
            slot['option_count']=len(slot['options'])
            slot['not_insertable']=[e for e in slot.get('not_insertable',[]) if e['value']!='block']
    catalog['components'] = components
    catalog['objects'] = objects
    catalog['scenes'] = scenes
    # One flat array for a tile grid; slot and requires carry the rules, so the
    # front end reads them rather than encoding which tile excludes which.
    catalog['tiles'] = environments + scenes + components + objects
    catalog['notes'] = notes
    return catalog


def _resolve(root, catalog, ident, field):
    entry = next((e for e in catalog['tiles'] if e['id'] == ident), None)
    if entry is None or not entry.get(field):
        return None
    file = (root / entry[field]).resolve()
    if not file.is_relative_to(root) or not file.is_file():
        return None
    return file


def thumbnail(root, ident, live):
    root = Path(root).resolve()
    return _resolve(root, scene_catalog(root, live), ident, 'thumbnail')


def surround_geometry(root, ident, live):
    """World surround geometry: walls, floors, ceilings, ground sheets. Visual only."""
    root = Path(root).resolve()
    catalog = scene_catalog(root, live)
    entry = next((s for s in catalog.get('surrounds', []) if s['id'] == ident), None)
    if entry is None or not entry.get('geometry'):
        return None
    file = (root / entry['geometry']).resolve()
    if not file.is_relative_to(root) or not file.is_file():
        return None
    return file


def object_geometry(root, ident, live):
    """Constructed object geometry, so the viewer can place the object it sees on the tile."""
    root = Path(root).resolve()
    return _resolve(root, scene_catalog(root, live), ident, 'geometry')
