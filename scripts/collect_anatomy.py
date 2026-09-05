#!/usr/bin/env python3
"""Rebuild auditable anatomy indexes from official downloaded archives.

Uses only Python's standard library. No inferred world positions or subject
calibration: OpenSim local frames/functions remain explicit source parameters.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/anatomy'
OUT = ROOT / 'data/derived/anatomy'
REVISION = 'd9b05d470b1a481c222372c85b75772faf8f7792'
UNITS = {
    'mass': 'kg', 'mass_center': 'm', 'inertia': 'kg*m^2',
    'inertia_xx': 'kg*m^2', 'inertia_yy': 'kg*m^2', 'inertia_zz': 'kg*m^2',
    'inertia_xy': 'kg*m^2', 'inertia_xz': 'kg*m^2', 'inertia_yz': 'kg*m^2',
    'max_isometric_force': 'N', 'optimal_fiber_length': 'm',
    'tendon_slack_length': 'm', 'pennation_angle_at_optimal': 'rad',
    'max_contraction_velocity': 'optimal_fiber_lengths/s',
    'activation_time_constant': 's', 'deactivation_time_constant': 's',
    'location': 'm', 'translation': 'm', 'orientation': 'rad',
    'scale_factors': 'dimensionless', 'gravity': 'm/s^2',
}


def rel(path):
    return str(path.relative_to(ROOT))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write(name, data):
    path = OUT / name
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def xml_record(node):
    """Preserve XML children in order, including repeat names/functions/frames."""
    return {'type': node.tag, 'attributes': dict(node.attrib),
            'text': (node.text or '').strip(),
            'children': [xml_record(c) for c in node]}


def leaves(node):
    return {c.tag: {'value': (c.text or '').strip(),
                    'unit': UNITS.get(c.tag, 'unspecified; consult source schema'),
                    'basis': 'source model parameter; no independent measurement or calibration performed by IHM'}
            for c in node if not len(c)}


def parse_model(path, base, geometry_by_name):
    # Legacy OpenSim 3 emits C++ class names with ::, invalid XML namespaces.
    # Escape only tag names during parsing and restore their original spelling.
    raw = path.read_text()
    normalized = re.sub(r'(<\/?)([A-Za-z_][\w:]*)',
                        lambda m: m[1] + m[2].replace('::', '__OPENSIM_SCOPE__'), raw)
    root = ET.fromstring(normalized)
    for element in root.iter():
        element.tag = element.tag.replace('__OPENSIM_SCOPE__', '::')
    model = root.find('Model')
    if model is None:
        raise ValueError(f'No Model in {path}')
    bodies = list(model.findall('./BodySet/objects/Body'))
    joints = list(model.findall('./JointSet/objects/*'))
    # OpenSim 3 serialized joints inside the child body.
    for body in bodies:
        joints.extend(body.findall('./Joint/*'))
    forces = model.findall('./ForceSet/objects/*')
    muscles = [f for f in forces if 'Muscle' in f.tag]
    meshes = []
    for element in model.iter():
        if element.tag not in ('mesh_file', 'geometry_file') or not element.text:
            continue
        ref = element.text.strip()
        candidates = [path.parent/ref, path.parent/'Geometry'/ref, base/'Geometry'/ref]
        resolved = next((p for p in candidates if p.is_file()), None)
        fallback = geometry_by_name.get(Path(ref).name, []) if resolved is None else []
        meshes.append({'reference': ref, 'resolved_path': rel(resolved) if resolved else None,
                       'other_filename_matches': [rel(p) for p in fallback],
                       'coordinate_frame': 'Defined by owning geometry and its parent frame in source XML'})
    muscle_records = []
    for muscle in muscles:
        pts = muscle.findall('./GeometryPath/PathPointSet/objects/*')
        muscle_records.append({'name': muscle.get('name'), 'type': muscle.tag,
            'parameters': leaves(muscle),
            'path_points': [dict(xml_record(p),
                role=('path_origin' if i == 0 else 'path_insertion' if i == len(pts)-1 else 'via_point'),
                coordinate_frame=p.findtext('socket_parent_frame') or p.findtext('body'),
                spatial_unit=model.findtext('length_units') or 'OpenSim default SI meters',
                location_semantics='Body/frame-local model path definition; moving/conditional functions are preserved, not evaluated')
                for i, p in enumerate(pts)],
            'path_wraps': [xml_record(p) for p in muscle.findall('./GeometryPath/PathWrapSet/objects/*')],
            'source_definition': xml_record(muscle)})
    return {'schema_version': 1, 'name': model.get('name'), 'source_path': rel(path),
        'sha256': digest(path), 'source_revision': REVISION,
        'evidence_class': 'published_model_parameters; not raw subject measurements',
        'subject_specific_calibration': False,
        'units': {'length': model.findtext('length_units'), 'force': model.findtext('force_units'),
                  'angles': 'radians in OpenSim model parameters; motion-file conventions are separate'},
        'credits': model.findtext('credits'), 'publications': model.findtext('publications'),
        'model_properties': leaves(model), 'defaults': [xml_record(d) for d in model.findall('./defaults/*')],
        'bodies': [{'name': b.get('name'), 'parameters': leaves(b), 'source_definition': xml_record(b)} for b in bodies],
        'joints': [xml_record(j) for j in joints], 'muscles': muscle_records,
        'geometry_references': meshes,
        'counts': {'bodies': len(bodies), 'joints': len(joints), 'muscles': len(muscles),
                   'path_points': sum(len(m['path_points']) for m in muscle_records),
                   'geometry_references': len(meshes)}}


def extract(archive, destination, strip_first=False):
    with zipfile.ZipFile(archive) as z:
        for member in z.infolist():
            parts = Path(member.filename).parts[1:] if strip_first else Path(member.filename).parts
            if member.is_dir() or not parts:
                continue
            if '..' in parts or Path(*parts).is_absolute():
                raise ValueError('Unsafe archive path')
            target = destination.joinpath(*parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists() or target.stat().st_size != member.file_size:
                target.write_bytes(z.read(member))


def parse_atlas(base):
    tables = {}
    for file in sorted(base.glob('*.txt')):
        with file.open(encoding='utf-8-sig') as f:
            tables[file.stem] = list(csv.DictReader(f, delimiter='\t'))
    id_names = {}
    for row in tables.get('isa_element_parts', []):
        id_names.setdefault(row['element file id'], []).append({'concept_id': row['concept id'], 'name': row['name']})
    mesh_records = []
    for file in sorted(base.rglob('*.obj')):
        vertices = faces = 0
        minimum = [float('inf')]*3
        maximum = [float('-inf')]*3
        with file.open(errors='replace') as f:
            for line in f:
                if line.startswith('v '):
                    xyz = list(map(float, line.split()[1:4]))
                    if len(xyz) != 3:
                        raise ValueError(f'Invalid vertex: {file}')
                    vertices += 1
                    minimum = [min(a,b) for a,b in zip(minimum, xyz)]
                    maximum = [max(a,b) for a,b in zip(maximum, xyz)]
                elif line.startswith('f '):
                    faces += 1
        mesh_records.append({'element_id': file.stem, 'source_path': rel(file),
            'sha256': digest(file), 'bytes': file.stat().st_size, 'vertices': vertices, 'faces': faces,
            'bounds_in_source_coordinates': [minimum, maximum] if vertices else None,
            'concepts': id_names.get(file.stem, [])})
    return {'schema_version': 1, 'source': 'BodyParts3D 4.0, DBCLS',
            'license': 'CC-BY-4.0 per official archive license page updated 2025-02-27',
            'attribution': 'BodyParts3D, © The Database Center for Life Science licensed under CC Attribution 4.0 International',
            'subject_specific_calibration': False,
            'coordinate_semantics': 'Original atlas coordinate system; no transform to OpenSim or runtime applied',
            'coordinate_unit': 'unverified source units; no inferred conversion applied',
            'mesh_detail': 'Official 99% polygon reduction distribution; not full-resolution anatomy',
            'counts': {'meshes': len(mesh_records), 'vertices': sum(m['vertices'] for m in mesh_records),
                       'faces': sum(m['faces'] for m in mesh_records),
                       'table_rows': {k:len(v) for k,v in tables.items()}},
            'meshes': mesh_records, 'tables': tables}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--skip-extract', action='store_true')
    ap.add_argument('--download', action='store_true', help='Download absent archives and metadata from their official sources')
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    opensim = RAW/'opensim-models'
    atlas = RAW/'bodyparts3d'
    if args.download:
        sources = {opensim/'source.zip': 'https://codeload.github.com/opensim-org/opensim-models/zip/'+REVISION}
        for name in ('isa_BP3D_4.0_obj_99.zip', 'isa_parts_list_e.txt', 'partof_parts_list_e.txt',
                     'isa_inclusion_relation_list.txt', 'partof_inclusion_relation_list.txt',
                     'isa_element_parts.txt', 'partof_element_parts.txt', 'README_e.html'):
            sources[atlas/name] = 'https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/'+name
        for name in ('lic.html', 'download.html'):
            sources[atlas/name] = 'https://dbarchive.biosciencedbc.jp/en/bodyparts3d/'+name
        for destination, url in sources.items():
            if not destination.exists():
                destination.parent.mkdir(parents=True, exist_ok=True)
                temporary = destination.with_suffix(destination.suffix+'.partial')
                urllib.request.urlretrieve(url, temporary)
                temporary.replace(destination)
    if not args.skip_extract:
        extract(opensim/'source.zip', opensim/'source', strip_first=True)
        extract(atlas/'isa_BP3D_4.0_obj_99.zip', atlas/'meshes')
    base = opensim/'source'
    geometry_by_name = {}
    for path in base.rglob('*'):
        if path.is_file() and path.suffix.lower() in ('.vtp', '.obj', '.stl'):
            geometry_by_name.setdefault(path.name, []).append(path)
    models = []
    for path in sorted((base/'Models').rglob('*.osim')):
        record = parse_model(path, base, geometry_by_name)
        filename = 'opensim__' + '__'.join(path.relative_to(base/'Models').with_suffix('').parts) + '.json'
        write(filename, record)
        models.append({'name': record['name'], 'source_path': record['source_path'],
                       'derived_path': rel(OUT/filename), 'counts': record['counts']})
    write('opensim_index.json', {'revision': REVISION, 'models': models,
        'counts': {'models': len(models), 'mesh_files': sum(len(v) for v in geometry_by_name.values()),
                   **{k:sum(m['counts'][k] for m in models) for k in ('bodies','joints','muscles','path_points')}}})
    atlas_record = parse_atlas(atlas)
    write('bodyparts3d_index.json', atlas_record)
    inventory = []
    for path in sorted(RAW.rglob('*')):
        if not path.is_file() or path.name == 'inventory.json':
            continue
        if 'opensim-models' in path.parts:
            source_url = ('https://raw.githubusercontent.com/opensim-org/opensim-models/'+REVISION+'/'+str(path.relative_to(base))) if path.is_relative_to(base) else 'https://codeload.github.com/opensim-org/opensim-models/zip/'+REVISION
            license_status = 'No repository-level license found; retain embedded credits/publications; asset redistribution terms unverified'
        else:
            name = path.name if not path.is_relative_to(atlas/'meshes') else 'isa_BP3D_4.0_obj_99.zip'
            source_url = 'https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/'+name
            if name in ('lic.html','download.html'):
                source_url = 'https://dbarchive.biosciencedbc.jp/en/bodyparts3d/'+name
            license_status = 'CC-BY-4.0; see locally retained lic.html'
        inventory.append({'path':rel(path), 'bytes':path.stat().st_size, 'sha256':digest(path),
                          'source_url':source_url, 'license_status':license_status})
    write('inventory.json', {'schema_version':1, 'download_date':'2026-09-05',
        'opensim_revision':REVISION, 'bodyparts3d_version':'4.0 (archive Last-Modified 2013-05-22)',
        'raw_bytes':sum(f['bytes'] for f in inventory), 'files':inventory})
    print(json.dumps({'opensim_models':len(models), 'opensim_counts':json.loads((OUT/'opensim_index.json').read_text())['counts'],
                     'atlas_counts':atlas_record['counts'], 'raw_files':len(inventory),
                     'raw_bytes':sum(f['bytes'] for f in inventory)}, indent=2))


if __name__ == '__main__':
    main()
