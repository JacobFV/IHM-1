#!/usr/bin/env python3
"""Geometry + naming primitives for the articular pose-defect audit.

Read-only over `data/derived/canonical/anatomy.json`. Nothing here writes.
"""
import gzip
import hashlib
import json
import re

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

try:
    import igl
except ImportError:  # pragma: no cover - the numpy fallback is exact but quadratic
    igl = None

WELD_DECIMALS = 12


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_geometry(path):
    payload = json.loads(gzip.open(path).read())
    return (np.asarray(payload['positions'], dtype=float).reshape(-1, 3),
            np.asarray(payload['indices'], dtype=np.int64).reshape(-1, 3))


def face_normals(vertices, faces):
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    cross = np.cross(b - a, c - a)
    area = np.linalg.norm(cross, axis=1)
    safe = np.where(area > 0, area, 1.)
    return cross / safe[:, None], area / 2


def signed_volume(vertices, faces):
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    return float(np.einsum('ij,ij->', a, np.cross(b - a, c - a)) / 6)


def topology(vertices, faces):
    """Welded edge-incidence diagnostics; tells us how far from closed a bone is."""
    _, inverse = np.unique(np.round(vertices, WELD_DECIMALS), axis=0, return_inverse=True)
    welded = inverse[faces]
    kept = ((welded[:, 0] != welded[:, 1]) & (welded[:, 1] != welded[:, 2]) & (welded[:, 0] != welded[:, 2]))
    welded = welded[kept]
    directed = np.stack([welded[:, [0, 1]], welded[:, [1, 2]], welded[:, [2, 0]]], 1).reshape(-1, 2)
    _, directed_counts = np.unique(directed, axis=0, return_counts=True)
    _, counts = np.unique(np.sort(directed, axis=1), axis=0, return_counts=True)
    return {'boundary_edges': int((counts == 1).sum()), 'nonmanifold_edges': int((counts > 2).sum()),
            'orientation_consistent': bool((directed_counts == 1).all()),
            'closed': bool((counts == 2).all())}


def closest_on_triangles(points, a, b, c, block=512):
    """Exact closest point on a triangle soup. Returns (distance, closest, face_index)."""
    ab, ac = b - a, c - a
    best_d = np.full(len(points), np.inf)
    best_p = np.zeros((len(points), 3))
    best_i = np.zeros(len(points), dtype=np.int64)
    for start in range(0, len(points), block):
        q = points[start:start + block]
        ap = q[:, None, :] - a[None]
        d1 = np.einsum('kj,ikj->ik', ab, ap)
        d2 = np.einsum('kj,ikj->ik', ac, ap)
        bp = q[:, None, :] - b[None]
        d3 = np.einsum('kj,ikj->ik', ab, bp)
        d4 = np.einsum('kj,ikj->ik', ac, bp)
        cp_ = q[:, None, :] - c[None]
        d5 = np.einsum('kj,ikj->ik', ab, cp_)
        d6 = np.einsum('kj,ikj->ik', ac, cp_)
        va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
        total = va + vb + vc
        safe = np.where(total != 0, total, 1.)
        closest = a[None] + (vb / safe)[..., None] * ab[None] + (vc / safe)[..., None] * ac[None]
        for mask, corner in (((d1 <= 0) & (d2 <= 0), a), ((d3 >= 0) & (d4 <= d3), b), ((d6 >= 0) & (d5 <= d6), c)):
            closest = np.where(mask[..., None], np.broadcast_to(corner[None], closest.shape), closest)
        for mask, origin, edge, num, den in (
                ((vc <= 0) & (d1 >= 0) & (d3 <= 0), a, ab, d1, d1 - d3),
                ((vb <= 0) & (d2 >= 0) & (d6 <= 0), a, ac, d2, d2 - d6),
                ((va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0), b, c - b, d4 - d3, (d4 - d3) + (d5 - d6))):
            t = np.where(den != 0, num / np.where(den != 0, den, 1.), 0.)
            closest = np.where(mask[..., None], origin[None] + t[..., None] * edge[None], closest)
        dist = np.linalg.norm(q[:, None, :] - closest, axis=2)
        pick = dist.argmin(1)
        rows = np.arange(len(q))
        best_d[start:start + block] = dist[rows, pick]
        best_p[start:start + block] = closest[rows, pick]
        best_i[start:start + block] = pick
    return best_d, best_p, best_i


def closest_on_mesh(points, vertices, faces):
    """(distance, closest point, index into `faces`) - libigl AABB when available."""
    if igl is not None:
        squared, which, closest = igl.point_mesh_squared_distance(
            np.ascontiguousarray(points, dtype=np.float64),
            np.ascontiguousarray(vertices, dtype=np.float64),
            np.ascontiguousarray(faces, dtype=np.int32))
        return np.sqrt(np.maximum(squared, 0)), closest, np.asarray(which, dtype=np.int64)
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    return closest_on_triangles(points, a, b, c)


def link_clusters(points, radius):
    """Single-linkage clusters at `radius` via connected components of the radius graph."""
    pairs = cKDTree(points).query_pairs(radius, output_type='ndarray')
    if not len(pairs):
        return np.arange(len(points))
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                       shape=(len(points), len(points)))
    return connected_components(graph, directed=False)[1]


ORD = ('first second third fourth fifth sixth seventh eighth ninth tenth eleventh twelfth').split()
DIGIT_HAND = {'thumb': 1, 'index finger': 2, 'middle finger': 3, 'ring finger': 4, 'little finger': 5}
DIGIT_FOOT = {'big toe': 1, 'second toe': 2, 'third toe': 3, 'fourth toe': 4, 'little toe': 5}
CARPALS = ('scaphoid', 'lunate', 'triquetral', 'pisiform', 'trapezium', 'trapezoid', 'capitate', 'hamate')
TARSALS = ('talus', 'calcaneus', 'navicular', 'cuboid', 'medial cuneiform', 'intermediate cuneiform',
           'lateral cuneiform')
CRANIAL = ('frontal bone', 'parietal bone', 'occipital bone', 'temporal bone', 'sphenoid bone', 'ethmoid',
           'vomer', 'maxilla', 'palatine bone', 'zygomatic bone', 'nasal bone', 'lacrimal bone',
           'inferior nasal concha')


def token(name):
    """(kind, side, index) for a canonical bone name; kind is the anatomical class."""
    lowered = name.lower()
    side = 'l' if 'left' in lowered else ('r' if 'right' in lowered else '')
    if 'tooth' in lowered:
        return ('tooth', side, 0)
    for kind, pattern in (('phalanx', r'(distal|middle|proximal) phalanx of (?:left |right )?(.+)'),):
        m = re.fullmatch(pattern, lowered)
        if m:
            level, digit = m.group(1), m.group(2)
            if digit in DIGIT_HAND:
                return (f'{level}_phalanx_hand', side, DIGIT_HAND[digit])
            if digit in DIGIT_FOOT:
                return (f'{level}_phalanx_foot', side, DIGIT_FOOT[digit])
    m = re.fullmatch(r'(?:left |right )?(\w+) metacarpal bone', lowered)
    if m and m.group(1) in ORD:
        return ('metacarpal', side, ORD.index(m.group(1)) + 1)
    m = re.fullmatch(r'(?:left |right )?(\w+) metatarsal bone', lowered)
    if m and m.group(1) in ORD:
        return ('metatarsal', side, ORD.index(m.group(1)) + 1)
    m = re.fullmatch(r'(?:left |right )?(\w+) rib', lowered)
    if m and m.group(1) in ORD:
        return ('rib', side, ORD.index(m.group(1)) + 1)
    m = re.fullmatch(r'(\w+) (cervical|thoracic|lumbar) vertebra', lowered)
    if m:
        return (f'vertebra_{m.group(2)}', '', ORD.index(m.group(1)) + 1)
    if lowered == 'atlas':
        return ('vertebra_cervical', '', 1)
    if lowered == 'axis':
        return ('vertebra_cervical', '', 2)
    if lowered.startswith('intervertebral disk'):
        return ('disc', '', 0)
    if 'navicular bone of' in lowered:
        return ('navicular', side, 0)
    if 'sesamoid' in lowered:
        return ('sesamoid', side, 0)
    for kind in TARSALS:
        if re.search(r'\b' + kind + r'\b', lowered):
            return (kind.replace(' ', '_'), side, 0)
    for kind in CARPALS:
        if re.search(r'\b' + kind + r'\b', lowered):
            return (kind, side, 0)
    for kind in ('femur', 'tibia', 'fibula', 'patella', 'humerus', 'radius', 'ulna', 'scapula', 'clavicle',
                 'sacrum', 'coccyx', 'mandible', 'manubrium', 'hyoid bone', 'xiphoid process',
                 'body of sternum', 'hip bone'):
        if lowered.endswith(kind):
            return (kind.replace(' ', '_'), side, 0)
    for kind in CRANIAL:
        if lowered.endswith(kind):
            return ('cranial:' + kind.replace(' ', '_'), side, 0)
    return ('other:' + lowered, side, 0)


HAND_CARPAL = set(CARPALS)
FOOT_TARSAL = {t.replace(' ', '_') for t in TARSALS}
CUNEIFORMS = {'medial_cuneiform', 'intermediate_cuneiform', 'lateral_cuneiform'}


SPINE_BASE = {'vertebra_cervical': 0, 'vertebra_thoracic': 7, 'vertebra_lumbar': 19}


def spinal_level(kind, index):
    return SPINE_BASE[kind] + index if kind in SPINE_BASE else None


def joint_of(ta, tb):
    """(label, joint_class) for an ordered-agnostic pair of bone tokens, or None."""
    (ka, sa, ia), (kb, sb, ib) = sorted([ta, tb])
    pair = {ka, kb}
    same_side = sa == sb and sa != ''
    same_ray = ia == ib

    if 'tooth' in pair:
        return ('dental', 'dental')
    if ka.startswith('cranial:') and kb.startswith('cranial:'):
        return ('cranial_suture', 'suture')
    if ka.startswith('cranial:') or kb.startswith('cranial:'):
        other = kb if ka.startswith('cranial:') else ka
        skull = ka if ka.startswith('cranial:') else kb
        if other == 'mandible' and 'temporal' in skull:
            return ('temporomandibular', 'synovial_fibrocartilaginous')
        if other == 'vertebra_cervical' and 'occipital' in skull:
            return ('atlanto_occipital', 'synovial')
        return ('craniofacial_other', 'suture')
    if 'disc' in pair:
        return ('discovertebral' if 'vertebra' in ka or 'vertebra' in kb else 'disc_other', 'amphiarthrosis')
    if ka.startswith('vertebra') and kb.startswith('vertebra'):
        if {ia, ib} == {1, 2} and ka == kb == 'vertebra_cervical':
            return ('atlantoaxial', 'synovial')
        if abs(spinal_level(ka, ia) - spinal_level(kb, ib)) != 1:
            return ('nonadjacent_vertebral', 'other')
        return ('zygapophyseal', 'synovial')
    if 'rib' in pair:
        if not (ka.startswith('vertebra') or kb.startswith('vertebra')):
            return ('costal_other', 'other')
        vertebra, rib = ((ka, ia), ib) if ka.startswith('vertebra') else ((kb, ib), ia)
        if vertebra[0] != 'vertebra_thoracic' or not 0 <= vertebra[1] - rib <= 1:
            return ('rib_nonarticular_contact', 'other')
        return ('costovertebral', 'synovial')
    if pair == {'manubrium', 'body of sternum'.replace(' ', '_')}:
        return ('manubriosternal', 'symphysis')
    if pair == {'body_of_sternum', 'xiphoid_process'}:
        return ('xiphisternal', 'symphysis')
    if pair == {'hip_bone', 'sacrum'}:
        return ('sacroiliac', 'synovial_fibrocartilaginous')
    if pair == {'hip_bone'} and sa != sb:
        return ('pubic_symphysis', 'symphysis')
    if pair == {'hip_bone', 'femur'} and same_side:
        return ('coxofemoral', 'synovial')
    if pair == {'femur', 'tibia'} and same_side:
        return ('tibiofemoral', 'synovial')
    if pair == {'femur', 'patella'} and same_side:
        return ('patellofemoral', 'synovial')
    if pair == {'tibia', 'patella'}:
        return ('patellar_ligament_span', 'non_articular')
    if pair == {'tibia', 'fibula'} and same_side:
        return ('tibiofibular', 'synovial')
    if pair in ({'tibia', 'talus'}, {'fibula', 'talus'}) and same_side:
        return ('talocrural', 'synovial')
    if pair == {'talus', 'calcaneus'} and same_side:
        return ('subtalar', 'synovial')
    if pair == {'talus', 'navicular'} and same_side:
        return ('talonavicular', 'synovial')
    if pair == {'calcaneus', 'cuboid'} and same_side:
        return ('calcaneocuboid', 'synovial')
    if pair == {'navicular', 'cuboid'} and same_side:
        return ('cubonavicular', 'synovial')
    if 'navicular' in pair and (pair - {'navicular'}) <= CUNEIFORMS and same_side:
        return ('cuneonavicular', 'synovial')
    if pair <= CUNEIFORMS | {'cuboid'} and same_side:
        return ('intertarsal_distal', 'synovial')
    if 'metatarsal' in pair and (pair - {'metatarsal'}) <= CUNEIFORMS | {'cuboid'} and same_side:
        return ('tarsometatarsal', 'synovial')
    if ka == kb == 'metatarsal' and same_side:
        return ('intermetatarsal', 'synovial') if abs(ia - ib) == 1 else ('nonadjacent_ray', 'other')
    if pair == {'metatarsal', 'proximal_phalanx_foot'} and same_side and same_ray:
        return ('metatarsophalangeal', 'synovial')
    if pair == {'metatarsal', 'sesamoid'} and same_side:
        return ('metatarsosesamoid', 'synovial')
    if pair == {'proximal_phalanx_foot', 'middle_phalanx_foot'} and same_side and same_ray:
        return ('interphalangeal_foot_proximal', 'synovial')
    if pair == {'middle_phalanx_foot', 'distal_phalanx_foot'} and same_side and same_ray:
        return ('interphalangeal_foot_distal', 'synovial')
    if pair == {'proximal_phalanx_foot', 'distal_phalanx_foot'} and same_side and same_ray and ia == 1:
        return ('interphalangeal_hallux', 'synovial')
    if pair == {'scapula', 'humerus'} and same_side:
        return ('glenohumeral', 'synovial')
    if pair == {'scapula', 'clavicle'} and same_side:
        return ('acromioclavicular', 'synovial')
    if pair == {'clavicle', 'manubrium'}:
        return ('sternoclavicular', 'synovial_fibrocartilaginous')
    if pair == {'humerus', 'ulna'} and same_side:
        return ('humeroulnar', 'synovial')
    if pair == {'humerus', 'radius'} and same_side:
        return ('humeroradial', 'synovial')
    if pair == {'radius', 'ulna'} and same_side:
        return ('radioulnar', 'synovial')
    if 'radius' in pair and (pair - {'radius'}) <= HAND_CARPAL and same_side:
        return ('radiocarpal', 'synovial')
    if 'ulna' in pair and (pair - {'ulna'}) <= HAND_CARPAL and same_side:
        return ('ulnocarpal', 'synovial_fibrocartilaginous')
    if pair <= HAND_CARPAL and same_side:
        proximal = {'scaphoid', 'lunate', 'triquetral', 'pisiform'}
        return ('intercarpal_proximal' if pair <= proximal else
                ('midcarpal' if len(pair & proximal) == 1 else 'intercarpal_distal'), 'synovial')
    if 'metacarpal' in pair and (pair - {'metacarpal'}) <= HAND_CARPAL and same_side:
        return ('carpometacarpal_thumb' if max(ia, ib) == 1 and 'trapezium' in pair else 'carpometacarpal',
                'synovial')
    if ka == kb == 'metacarpal' and same_side:
        return ('intermetacarpal', 'synovial') if abs(ia - ib) == 1 else ('nonadjacent_ray', 'other')
    if pair == {'metacarpal', 'proximal_phalanx_hand'} and same_side and same_ray:
        return ('metacarpophalangeal', 'synovial')
    if pair == {'proximal_phalanx_hand', 'middle_phalanx_hand'} and same_side and same_ray:
        return ('interphalangeal_hand_proximal', 'synovial')
    if pair == {'middle_phalanx_hand', 'distal_phalanx_hand'} and same_side and same_ray:
        return ('interphalangeal_hand_distal', 'synovial')
    if pair == {'proximal_phalanx_hand', 'distal_phalanx_hand'} and same_side and same_ray and ia == 1:
        return ('interphalangeal_thumb', 'synovial')
    if 'hyoid_bone' in pair:
        return ('hyoid_self', 'non_articular')
    return ('unclassified', 'other')
