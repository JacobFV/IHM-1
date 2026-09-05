"""One canonical body, with source surfaces and explicitly inferred registrations.

All coordinates are meters, X left / Y superior / Z anterior. Anatomical meshes
are reference surfaces, not automatically valid volumetric finite elements.
"""
from __future__ import annotations
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RBFInterpolator

MODEL_ID = 'ihm-body'
FRAME = 'bodyparts3d-display-m'
ROTATION = np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]])


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(value, separators=(',', ':'), allow_nan=False).encode()
    temporary = path.with_name(path.name + '.writing')
    if path.suffix == '.gz':
        with temporary.open('wb') as stream:
            with gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as zipped:
                zipped.write(raw)
    else:
        temporary.write_bytes(raw)
    temporary.replace(path)


def read_geometry(path):
    with gzip.open(path, 'rt') as stream:
        return json.load(stream)


class LandmarkRegistration:
    """Orientation-preserving affine fit plus smooth, regularized residual field.

    Bone bounding-box centers are correspondences, not measured fiducials. The
    0.002 smoothing prior trades exact interpolation for stable inter-bone fit.
    """
    def __init__(self, source, target, smoothing=.002):
        self.source = np.asarray(source, dtype=float)
        self.target = np.asarray(target, dtype=float)
        if self.source.shape != self.target.shape or self.source.ndim != 2 or self.source.shape[1] != 3:
            raise ValueError('Expected matching N by 3 landmarks')
        design = np.c_[self.source, np.ones(len(self.source))]
        if len(self.source) < 5 or np.linalg.matrix_rank(design) < 4 or not np.isfinite(design).all() or not np.isfinite(self.target).all():
            raise ValueError('Landmarks must be finite and span three dimensions')
        self.affine = np.linalg.lstsq(design, self.target, rcond=None)[0]
        if np.linalg.det(self.affine[:3]) <= 0:
            raise ValueError('Landmarks imply a reflected or collapsed body')
        self.smoothing = smoothing
        residual = self.target - design @ self.affine
        self.residual = RBFInterpolator(self.source, residual, smoothing=smoothing, kernel='thin_plate_spline')

    def transform(self, points):
        points = np.asarray(points, dtype=float)
        return np.c_[points, np.ones(len(points))] @ self.affine + self.residual(points)

    def report(self):
        distances = np.linalg.norm(self.transform(self.source) - self.target, axis=1)
        affine_distances = np.linalg.norm(np.c_[self.source, np.ones(len(self.source))] @ self.affine - self.target, axis=1)
        return {'method': 'affine plus regularized thin-plate-spline residual', 'landmark_count': len(self.source),
                'affine_4x3': self.affine.tolist(), 'affine_determinant': float(np.linalg.det(self.affine[:3])),
                'smoothing_prior': self.smoothing, 'affine_rms_m': float(np.sqrt(np.mean(affine_distances**2))),
                'fit_rms_m': float(np.sqrt(np.mean(distances**2))), 'fit_max_m': float(distances.max()),
                'interpretation': 'Geometric correspondence residual, not population variation or measurement confidence.'}

    def jacobian_determinants(self, points, step=1e-5):
        p = np.asarray(points, dtype=float)
        jacobian = np.stack([(self.transform(p + np.eye(3)[i]*step) - self.transform(p - np.eye(3)[i]*step))/(2*step) for i in range(3)], axis=-1)
        return np.linalg.det(jacobian)


def normalized_name(name):
    name = name.lower().replace('(', '').replace(')', '').strip()
    if name.endswith('.l'):
        return 'left ' + name[:-2]
    if name.endswith('.r'):
        return 'right ' + name[:-2]
    return name


def physical_role(name, system):
    name = name.lower()
    if name.startswith('cavity of '):
        return 'fluid_cavity'
    if system == 'lymphatic' and 'node' in name:
        return 'lymph_node_group'
    if system in ('arterial', 'venous'):
        return 'vascular'
    if 'ligament' in name or 'retinacul' in name:
        return 'ligament'
    if 'tendon' in name or 'aponeurosis' in name:
        return 'tendon'
    if 'cartilage' in name or 'intervertebral disc' in name:
        return 'cartilage'
    if system == 'skeletal':
        return 'rigid_bone'
    if system == 'muscular':
        return 'muscle'
    if system == 'nervous':
        return 'nerve'
    if system in ('arterial', 'venous'):
        return 'vascular'
    if name == 'skin':
        return 'skin'
    if system == 'lymphatic' and 'node' in name:
        return 'lymph_node_group'
    if system == 'connective':
        return 'connective_tissue'
    return 'soft_organ'


def mesh_properties(vertices, faces):
    vertices, faces = np.asarray(vertices, dtype=float), np.asarray(faces, dtype=int)
    low, high = vertices.min(0), vertices.max(0)
    centroid = (low + high)/2
    triangles = vertices[faces]
    area = .5 * np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]), axis=1).sum()
    signed = np.einsum('ij,ij->i', triangles[:, 0]-centroid, np.cross(triangles[:, 1]-centroid, triangles[:, 2]-centroid)).sum()/6
    edges = np.sort(np.concatenate((faces[:, :2], faces[:, 1:], faces[:, [2, 0]])), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    watertight = bool(np.all(counts == 2))
    _, _, axes = np.linalg.svd(vertices-centroid, full_matrices=False)
    return {'bounds_m': {'min': low.tolist(), 'max': high.tolist()}, 'centroid_m': centroid.tolist(),
            'centroid_definition': 'axis-aligned surface bounding-box center; not center of mass',
            'surface_area_m2': float(area), 'volume_m3': float(abs(signed)) if watertight else None,
            'volume_method': 'absolute signed surface integral; edge incidence closed, self-intersection not certified' if watertight else 'unknown: surface has boundary or nonmanifold edges',
            'watertight_edge_incidence': watertight, 'principal_axis': axes[0].tolist(),
            'source_vertex_count': len(vertices), 'source_face_count': len(faces)}


def verify_assembly(assembly, root):
    root = Path(root)
    entities = assembly['entities']
    by_id = {e['id']: e for e in entities}
    assert len(by_id) == len(entities), 'Duplicate canonical identity'
    assert assembly['model_id'] == MODEL_ID and assembly['frame']['id'] == FRAME
    assert assembly['frame']['units'] == 'm'
    assert sum(e['evidence_kind'] == 'source_geometry' for e in entities) == 2234
    assert not any(e['evidence_kind'] == 'registered_geometry' and e['role'] == 'rigid_bone' for e in entities)
    checked = set()
    for e in entities:
        assert e['model_id'] == MODEL_ID
        assert e['evidence_kind'] in ('source_geometry', 'registered_geometry', 'synthesized_layer', 'registered_structural_graph')
        assert len(e['centroid_m']) == 3 and np.isfinite(e['centroid_m']).all()
        assert np.all(np.array(e['bounds_m']['max']) >= np.array(e['bounds_m']['min']))
        assert e['uncertainty']['biological']['confidence_percent'] is None
        geometry = e['reference_geometry']
        assert geometry['units'] == 'm' and geometry['frame'] == FRAME
        for record in [geometry, *e['provenance'].get('files', []), *([e['graph']] if 'graph' in e else [])]:
            if record['path'] not in checked:
                assert sha256(root / record['path']) == record['sha256'], record['path']
                checked.add(record['path'])
        if e['evidence_kind'] != 'source_geometry':
            assert e['assumptions'] and e['provenance']['source_ids']
        for connection in e['connections']:
            assert connection['entity_id'] in by_id and connection['entity_id'] != e['id']
    lung_names = {'inferior lobe of left lung', 'inferior lobe of right lung', 'middle lobe of right lung', 'superior lobe of left lung', 'superior lobe of right lung'}
    lungs = [e for e in entities if e['name'] in lung_names]
    assert {e['name'] for e in lungs} == lung_names, 'Missing lung parenchymal lobes'
    for e in lungs:
        assert e['system'] == 'respiratory' and e['role'] == 'soft_organ'
        assert e['centroid_m'][0] > 0 if 'left' in e['name'] else e['centroid_m'][0] < 0
    network = next(e for e in entities if e['role'] == 'lymphatic_network')
    graph = json.loads((root / network['graph']['path']).read_text())
    original_graph = json.loads((root / graph['source_graph']['path']).read_text())
    assert len(graph['nodes']) == 996 and len(graph['edges']) == 1117
    assert sum(n['is_lymph_node'] for n in graph['nodes']) == 272
    assert not graph['directed']
    assert all(g['source_from'] == o['source_from'] and g['source_to'] == o['source_to'] and g['source_length_mm'] == o['source_length_mm'] for g,o in zip(graph['edges'], original_graph['edges']))
    assert all(np.isfinite(n['position_m']).all() for n in graph['nodes'])
    for a in graph['node_group_associations']:
        assert a['candidate_entity_id'] in by_id
        assert a['accepted'] == (a['distance_m'] <= .06)
    reg = assembly['registrations']['z_anatomy']
    assert sha256(root / reg['source_index']['path']) == reg['source_index']['sha256']
    assert reg['jacobian_determinant_min_at_vertices'] > 0, 'Locally folded imported anatomy'
    assert reg['held_out_rms_m'] < .04, 'Registration held-out error above 4 cm engineering acceptance bound'
    assert reg['fit_rms_m'] < .015
    imported = [e for e in entities if e['evidence_kind'] == 'registered_geometry']
    for e in imported:
        name = e['name'].lower()
        if name.startswith('left '):
            assert e['centroid_m'][0] > -.005, e['name']
        if name.startswith('right '):
            assert e['centroid_m'][0] < .005, e['name']
        assert -.87 < e['centroid_m'][1] < .88
    roles = Counter(e['role'] for e in entities)
    for role in ('rigid_bone', 'muscle', 'ligament', 'nerve', 'vascular', 'soft_organ', 'skin', 'skin_layer', 'lymph_node_group'):
        assert roles[role], role
    return {'status': 'passed', 'entities': len(entities), 'verified_files': len(checked), 'roles': dict(roles),
            'source_geometry_entities': 2234, 'registered_additions': len(imported),
            'registration_held_out_rms_m': reg['held_out_rms_m'],
            'checks': ['stable unique identities', 'source and derived byte hashes', 'physical units and frame', 'domain roles', 'connection referential integrity', 'source versus inferred lineage', 'held-out registration residual', 'positive local Jacobians at imported vertices', 'laterality and vertical envelope', 'exact registered lymphatic graph topology'],
            'limits': ['Geometric and numerical checks are not empirical calibration.', 'Positive sampled Jacobians do not prove global injectivity or collision-free tissue interfaces.']}
