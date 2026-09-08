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


def is_cavity_name(name):
    """Whether a name declares a cavity, in either the source or normalized form.

    normalized_name rewrites the Z-Anatomy side suffix as a prefix, so a cavity
    reads 'Cavity of concha.l' before it and 'left cavity of concha' after. A
    bare startswith test sees only the first and silently re-roles the second.
    """
    low = name.lower()
    return any(low.startswith(prefix + 'cavity of ') for prefix in ('', 'left ', 'right '))


def physical_role(name, system):
    name = name.lower()
    if is_cavity_name(name):
        return 'fluid_cavity'
    if system == 'lymphatic' and 'node' in name:
        return 'lymph_node_group'
    if system in ('arterial', 'venous'):
        return 'vascular'
    if 'ligament' in name or 'retinacul' in name:
        return 'ligament'
    if 'tendon' in name or 'aponeurosis' in name:
        return 'tendon'
    # BodyParts3D spells the fibrocartilaginous disc "intervertebral disk"; the
    # spelling variant fell through to rigid_bone and made every disc infinitely
    # stiff. The role follows the record's own FMA is-a anchor FMA55107
    # (cartilage organ). The vocabulary has no fibrocartilage role, so cartilage
    # is nearest, not exact. See DISC-ROLE-FROM-FMA in the assumption ledger.
    if 'cartilage' in name or 'intervertebral disc' in name or 'intervertebral disk' in name:
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


# Integumentary structures whose name says they are a named patch of the body
# wall rather than a tissue volume. Their enclosed volume is the skin they lie
# on, so counting it as tissue would double-count that skin. Owned here, not in
# the promotion lane, so the canonical assembly and every candidate builder
# reach the same role for the same name.
SURFACE_REGION_SYSTEMS = {'integumentary'}
SURFACE_REGION_TOKENS = ('region', 'triangle', 'fossa', 'dorsum', 'surface', 'arch of foot', 'helix',
                         'tragus', 'concha', 'auricle', 'lobule', 'angle of mouth', 'notch of auricle',
                         'philtrum', 'groove', 'crus', 'crura', 'apex of', 'tubercle of auricle')


NAME_STOPWORDS = {'of', 'the', 'a', 'muscle', 'bone', 'part'}


def name_key(name):
    """Side-preserving token multiset for cross-source name matching.

    Z-Anatomy writes 'X of hand.l' where BodyParts3D writes 'X of left hand';
    the ordering differs, the tokens do not. Owned here so the canonical
    assembly and the promotion lane cannot drift apart on what a name match is.
    """
    text = normalized_name(name).replace('-', ' ').replace(',', ' ')
    tokens = [t for t in text.replace('(', ' ').replace(')', ' ').split() if t not in NAME_STOPWORDS]
    return tuple(sorted(tokens))


def refine_role(name, system, base=None):
    """physical_role plus the topographic-region refinement.

    Returns (role, refinement_note). No BodyParts3D entity matches the rule, so
    applying it to the whole assembly re-roles nothing that was already there;
    it exists for the Z-Anatomy body-wall regions the display promotion adds.
    """
    base = physical_role(name, system) if base is None else base
    low = name.lower()
    if system in SURFACE_REGION_SYSTEMS and not is_cavity_name(low) \
       and any(token in low for token in SURFACE_REGION_TOKENS):
        return 'surface_region', 'topographic surface region of the body wall; not a tissue volume'
    return base, None


# Named structures whose source labels its two sides the wrong way round. Each
# entry is a measured mirror inversion in the acquired source, not a defect of
# any transform this repository applies: both sides exist, both are placed on a
# real side of the body, and the two labels are swapped. Listing them by name
# keeps the bilateral mirror check live for every other pair; a seventh swap
# fails the build rather than joining the list silently.
SOURCE_LATERALITY_DEFECTS = {
    'flexor pollicis brevis': 'BodyParts3D; the two thumb muscles sit 0.56 m apart on opposite hands under swapped labels',
    'middle pharyngeal constrictor': 'BodyParts3D; paraxial pair, 18 mm apart under swapped labels',
    'oblique arytenoid': 'BodyParts3D; midline laryngeal pair, 0.3 mm apart under swapped labels',
    'ascending lumbar vein': 'BodyParts3D; paravertebral pair, 36 mm apart under swapped labels',
    'eyelashes': 'Z-Anatomy; the .l and .r lash meshes are swapped, 66 mm apart',
    'lateral temporomandibular ligament': 'Z-Anatomy; the .l and .r ligaments are swapped, 108 mm apart',
}


def bilateral_mirror_census(entities):
    """Every 'left X'/'right X' name pair, and which of them are mirror-inconsistent.

    X is left-positive in this frame, so a genuine bilateral pair must place its
    left-named group at greater mean x than its right-named group. This replaces
    a per-entity sign test that could never have been right: 'right coronary
    artery', 'left portal vein' and 'left medial segment of liver IV' name a
    side of an organ, not a side of the body, and 20 acquired entities fail a
    sign test while being correctly placed -- it never saw them, because it only
    ran on imported rows. A pair, by contrast, is a real mirror claim.

    BodyParts3D authors several meshes under one name ('right fibular vein' is
    three), so the comparison is over the group mean rather than an arbitrary
    representative. A mesh sitting on the far side of the midline from its own
    group is reported as a within-group outlier and not asserted on: it is a
    single mislabelled branch inside a correctly mirrored pair, a different
    defect from a wholly swapped pair.
    """
    groups = {}
    for e in entities:
        groups.setdefault(e['name'].lower(), []).append(e)
    pairs, violations, outliers = 0, [], []
    for name in sorted(groups):
        if not name.startswith('left '):
            continue
        right = groups.get('right ' + name[5:])
        if right is None:
            continue
        left = groups[name]
        pairs += 1
        lx = float(np.mean([e['centroid_m'][0] for e in left]))
        rx = float(np.mean([e['centroid_m'][0] for e in right]))
        for group, mean in ((left, lx), (right, rx)):
            for e in group:
                if len(group) > 1 and e['centroid_m'][0]*mean < 0:
                    outliers.append({'id': e['id'], 'name': e['name'], 'centroid_x_m': e['centroid_m'][0],
                                     'group_mean_x_m': mean, 'group_size': len(group)})
        if not lx > rx:
            violations.append({'stem': name[5:], 'left_ids': [e['id'] for e in left],
                               'right_ids': [e['id'] for e in right],
                               'left_mean_x_m': lx, 'right_mean_x_m': rx, 'separation_m': abs(lx-rx),
                               'known_source_defect': name[5:] in SOURCE_LATERALITY_DEFECTS,
                               'source_defect_note': SOURCE_LATERALITY_DEFECTS.get(name[5:])})
    return {'bilateral_pairs': pairs, 'mirror_violations': violations,
            'unlisted_violations': [v for v in violations if not v['known_source_defect']],
            'within_group_outliers': outliers,
            'rule': "x is left-positive; a 'left X'/'right X' pair must place the left-named "
                    "group at greater mean x than the right-named group",
            'ledger': 'ihm.assembly.anatomy.SOURCE_LATERALITY_DEFECTS'}


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


BP_SOURCE_SURFACE_COUNT = 2234


def surface_identity_key(vertices, faces):
    """Decoded-array identity of one reference surface.

    The stored .json.gz bytes differ on gzip framing alone for surfaces the
    publisher authored twice, so the stored digest cannot see the duplication.
    This hashes the decoded arrays themselves, plus an order-independent facet
    set, so two rows collide only when they are the same surface, vertex for
    vertex and facet for facet. It is not a shape comparison: a re-indexed or
    translated copy of the same anatomy is deliberately NOT matched here.
    """
    vertices = np.ascontiguousarray(np.asarray(vertices, dtype=float))
    faces = np.ascontiguousarray(np.asarray(faces, dtype=np.int64))
    if vertices.ndim != 2 or vertices.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError('Surface identity requires N by 3 vertices and M by 3 faces')
    facets = np.unique(np.sort(faces, axis=1), axis=0)
    digest = hashlib.sha256()
    for block in (vertices.tobytes(), faces.tobytes(), np.ascontiguousarray(facets).tobytes()):
        digest.update(hashlib.sha256(block).digest())
    return digest.hexdigest()


def duplicate_surface_survivors(keys):
    """Map every duplicate-authored id to the id that survives the collapse.

    ``keys`` is {entity_id: surface_identity_key}. Nothing anatomical separates
    two rows that carry one surface, so the survivor is the lexicographically
    earlier canonical id; the choice is recorded, not inferred from evidence.
    Returns {dropped_id: survivor_id}.
    """
    groups = {}
    for identity, key in keys.items():
        groups.setdefault(key, []).append(identity)
    survivors = {}
    for members in groups.values():
        if len(members) < 2:
            continue
        keep = min(members)
        for member in members:
            if member != keep:
                survivors[member] = keep
    return survivors


def attach_regional_support(entities, replace=False):
    """Attach one nearest-bone regional support relation per non-layer entity.

    A support relation establishes a spatial substrate without claiming a
    histological attachment or a physiological flow connection. Bones are the
    ``rigid_bone`` rows only, so re-roling a structure out of ``rigid_bone``
    both removes it from the candidate set and gives it a support of its own.
    ``replace`` drops any existing relation first, so the rule can be re-applied
    to an already-built assembly and reach the same answer as a fresh build.
    """
    bones = [e for e in entities if e['role'] == 'rigid_bone']
    if not bones:
        raise ValueError('Regional support requires at least one rigid bone')
    centers = np.array([e['centroid_m'] for e in bones])
    index_by_id = {e['id']: i for i, e in enumerate(bones)}
    attached = 0
    for e in entities:
        if e['role'] in ('skin_layer', 'lymphatic_network'):
            continue
        if replace:
            e['connections'] = [c for c in e['connections'] if c.get('relation') != 'regional_spatial_support']
        distance = np.linalg.norm(centers - np.array(e['centroid_m']), axis=1)
        if e['id'] in index_by_id:
            distance[index_by_id[e['id']]] = np.inf
        index = int(np.argmin(distance))
        e['connections'].append({'entity_id': bones[index]['id'], 'relation': 'regional_spatial_support',
            'distance_m': float(distance[index]), 'evidence': 'inferred nearest bone bounding-box center; not an attachment or contact constraint'})
        attached += 1
    return attached


def verify_assembly(assembly, root):
    root = Path(root)
    entities = assembly['entities']
    by_id = {e['id']: e for e in entities}
    assert len(by_id) == len(entities), 'Duplicate canonical identity'
    assert assembly['model_id'] == MODEL_ID and assembly['frame']['id'] == FRAME
    assert assembly['frame']['units'] == 'm'
    # The BP scaffold is 2234 authored surfaces. Rows the publisher authored
    # twice are collapsed to one entity, so the entity count is the scaffold
    # minus the recorded collapse; a silent drop still fails this check.
    collapse = assembly.get('duplicate_surface_collapse', {'dropped': []})
    dropped = collapse['dropped']
    assert all(d['survivor'] in by_id and d['id'] not in by_id for d in dropped)
    source_geometry = sum(e['evidence_kind'] == 'source_geometry' for e in entities)
    assert source_geometry == BP_SOURCE_SURFACE_COUNT - len(dropped), source_geometry
    # Registered bones used to be forbidden outright. That gate was a proxy for
    # the real constraint, which is that the acquired BodyParts3D atlas alone
    # defines this body and its size: the registration landmarks are BP bone
    # centres, and a registered bone must therefore not restate a bone the
    # atlas already carries. Enforce that directly, so genuinely absent bones
    # -- ossicles, teeth, phalanges, sesamoids -- can be added without loosening
    # anything, while a registered bone that duplicates an acquired one by name
    # still fails.
    acquired_bone_names = {e['name'].lower() for e in entities
                           if e['evidence_kind'] == 'source_geometry' and e['role'] == 'rigid_bone'}
    registered_bones = [e for e in entities
                        if e['evidence_kind'] == 'registered_geometry' and e['role'] == 'rigid_bone']
    restated = [e['name'] for e in registered_bones if e['name'].lower() in acquired_bone_names]
    assert not restated, restated
    landmark_ids = {l['target_id'] for l in assembly['registrations']['z_anatomy']['landmarks']}
    assert not (landmark_ids & {e['id'] for e in registered_bones}), 'A registered bone anchors the registration'
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
        assert -.87 < e['centroid_m'][1] < .88, e['name']
    # Laterality is checked as a mirror relation over 'left X'/'right X' pairs
    # across the whole assembly, acquired rows included, rather than as a sign
    # test on imported rows only. The sign test could never have been right:
    # 'right coronary artery', 'left portal vein' and 'left medial segment of
    # liver IV' name a side of an organ, not a side of the body, and 20 acquired
    # entities fail it while being correctly placed -- it never saw them because
    # it only ran on imported rows. Every genuine mirror inversion left in the
    # sources is named in SOURCE_LATERALITY_DEFECTS; an unlisted one fails here.
    laterality = bilateral_mirror_census(entities)
    assert not laterality['unlisted_violations'], laterality['unlisted_violations']
    roles = Counter(e['role'] for e in entities)
    for role in ('rigid_bone', 'muscle', 'ligament', 'nerve', 'vascular', 'soft_organ', 'skin', 'skin_layer', 'lymph_node_group'):
        assert roles[role], role
    return {'status': 'passed', 'entities': len(entities), 'verified_files': len(checked), 'roles': dict(roles),
            'source_geometry_entities': source_geometry, 'bp_source_surface_count': BP_SOURCE_SURFACE_COUNT,
            'collapsed_duplicate_surfaces': len(dropped), 'registered_additions': len(imported),
            'registered_bones': len(registered_bones),
            'registration_held_out_rms_m': reg['held_out_rms_m'],
            'laterality': laterality,
            'checks': ['stable unique identities', 'source and derived byte hashes', 'physical units and frame', 'domain roles', 'connection referential integrity', 'source versus inferred lineage', 'held-out registration residual', 'positive local Jacobians at imported vertices', 'bilateral mirror consistency against a named source-defect ledger', 'registered bones do not restate an acquired bone', 'vertical envelope', 'exact registered lymphatic graph topology'],
            'limits': ['Geometric and numerical checks are not empirical calibration.', 'Positive sampled Jacobians do not prove global injectivity or collision-free tissue interfaces.']}
