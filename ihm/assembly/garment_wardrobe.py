"""Wardrobe of externally sourced garment meshes fitted to this body.

Garments are not synthesized here. Each is an acquired CC0 mesh authored against
the MakeHuman base mesh, registered onto this body's watertight outer envelope by
a recorded similarity transform, a piecewise-rigid limb pose correction and a
Laplacian-regularised shrinkwrap with a standoff. Fit residuals, penetration and
edge-length distortion are measured, not asserted.

Only numpy and scipy are imported so this module loads in both the main venv and
the out-of-process libigl venv.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

REFERENCE_HEIGHT_M = 1.7194712
SLAB_AREA_REJECT_M2 = 2.5
SOURCE_UNITS_PER_M = 10.0  # MakeHuman base mesh is authored in decimetres.


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def load_obj(path):
    """Wavefront OBJ vertices and fan-triangulated faces. Positions only."""
    positions = []
    triangles = []
    with open(path, 'r', errors='replace') as handle:
        for line in handle:
            if line.startswith('v '):
                parts = line.split()
                positions.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('f '):
                index = [int(token.split('/')[0]) - 1 for token in line.split()[1:]]
                if len(index) < 3:
                    raise ValueError(f'{path}: face with {len(index)} corners')
                for k in range(1, len(index) - 1):
                    triangles.append((index[0], index[k], index[k + 1]))
    if not positions or not triangles:
        raise ValueError(f'{path}: no vertices or no faces')
    v = np.asarray(positions, float)
    t = np.asarray(triangles, np.int64)
    if t.min() < 0 or t.max() >= len(v):
        raise ValueError(f'{path}: face index outside the vertex table')
    return v, t


def face_components(vertex_count, triangles):
    graph = coo_matrix((np.ones(3 * len(triangles)),
                        (np.concatenate([triangles[:, 0], triangles[:, 1], triangles[:, 2]]),
                         np.concatenate([triangles[:, 1], triangles[:, 2], triangles[:, 0]]))),
                       shape=(vertex_count, vertex_count))
    count, label = connected_components(graph, directed=False)
    return count, label


def keep_components(positions, triangles, wanted):
    """Retain only the listed face components, ordered by descending vertex count."""
    _, label = face_components(len(positions), triangles)
    order = sorted(np.unique(label), key=lambda k: -int((label == k).sum()))
    selected = {order[i] for i in wanted}
    triangles = triangles[np.isin(label[triangles[:, 0]], list(selected))]
    if not len(triangles):
        raise ValueError('component selection removed every face')
    return compact(positions, triangles)


def compact(positions, triangles):
    p = positions[triangles]
    area = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    keep = (area > 1e-14) & (triangles[:, 0] != triangles[:, 1]) & (triangles[:, 1] != triangles[:, 2]) \
        & (triangles[:, 0] != triangles[:, 2])
    triangles = triangles[keep]
    if not len(triangles):
        raise ValueError('every face was degenerate')
    used, inverse = np.unique(triangles, return_inverse=True)
    return positions[used], inverse.reshape(-1, 3).astype(np.int64)


def weld(positions, triangles, tolerance_m=1e-6):
    key = np.round(positions / tolerance_m).astype(np.int64)
    _, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
    return compact(positions[first], inverse[triangles])


def mesh_area_m2(positions, triangles):
    p = positions[triangles]
    return float(0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum())


def edge_table(triangles):
    edges = np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]])
    return np.unique(np.sort(edges, axis=1), axis=0)


def boundary_edge_count(triangles):
    edges = np.sort(np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]), axis=1)
    _, count = np.unique(edges, axis=0, return_counts=True)
    return int((count == 1).sum()), int((count > 2).sum())


def orient_outward(positions, triangles):
    """Flip face components whose winding points at their own component axis."""
    _, label = face_components(len(positions), triangles)
    triangles = triangles.copy()
    flips = 0
    for k in np.unique(label[triangles[:, 0]]):
        mask = label[triangles[:, 0]] == k
        p = positions[triangles[mask]]
        normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        radial = p.mean(1) - positions[label == k].mean(0)
        radial[:, 1] = 0.0
        if float(np.sum(normal * radial)) < 0:
            triangles[mask] = triangles[mask][:, [0, 2, 1]]
            flips += 1
    return triangles, flips


SLOTS = [
    {'id': 'underwear_bottom', 'label': 'Underwear (bottom)', 'region': 'pelvis', 'layer': 0},
    {'id': 'underwear_top', 'label': 'Underwear (top)', 'region': 'torso', 'layer': 0},
    {'id': 'torso_base', 'label': 'Torso base layer', 'region': 'torso', 'layer': 1},
    {'id': 'torso_outer', 'label': 'Torso outer layer', 'region': 'torso', 'layer': 2},
    {'id': 'legs', 'label': 'Legs', 'region': 'pelvis_legs', 'layer': 1},
    {'id': 'feet_inner', 'label': 'Feet inner layer', 'region': 'feet', 'layer': 0},
    {'id': 'feet_outer', 'label': 'Feet outer layer', 'region': 'feet', 'layer': 1},
    {'id': 'hands', 'label': 'Hands', 'region': 'hands', 'layer': 1},
    {'id': 'head', 'label': 'Head', 'region': 'head', 'layer': 1},
    {'id': 'neck', 'label': 'Neck', 'region': 'neck', 'layer': 1},
]

# One catalogue entry per retained source mesh. `components` selects face
# components by descending vertex count where one source file holds several
# garments; `standoff_mm` is the clearance the shrinkwrap targets.
CATALOGUE = [
    {'id': 'briefs', 'name': 'Briefs', 'slots': ['underwear_bottom'], 'colour': '#d7dade',
     'pack': 'underwear01', 'asset': 'wolgade_female_panties_01', 'obj': 'f_panties_01.obj', 'standoff_mm': 2.0},
    {'id': 'bra-top', 'name': 'Bra top', 'slots': ['underwear_top'], 'colour': '#d7dade',
     'pack': 'underwear01', 'asset': 'wolgade_female_top_01', 'obj': 'f_top_01.obj', 'standoff_mm': 2.0},
    {'id': 't-shirt', 'name': 'T-shirt', 'slots': ['torso_base'], 'colour': '#7d9ea6',
     'pack': 'shirts01', 'asset': 'elvs_crude_t-shirt_male', 'obj': 'crude_male_shirt.obj', 'standoff_mm': 4.0},
    {'id': 'polo-shirt', 'fit_mode': 'drape', 'name': 'Polo shirt', 'slots': ['torso_base'], 'colour': '#5f8a7d',
     'pack': 'shirts01', 'asset': 'namuhekam_male_polo_shirt', 'obj': 'Polo_t-shirt.obj', 'standoff_mm': 4.0},
    {'id': 'tucked-t-shirt', 'name': 'Tucked T-shirt', 'slots': ['torso_base'], 'colour': '#93887a',
     'pack': 'shirts01', 'asset': 'toigo_basic_tucked_t-shirt', 'obj': 't_shirt_basic_tucked.obj', 'standoff_mm': 4.0},
    {'id': 'fisherman-sweater', 'fit_mode': 'drape', 'name': 'Fisherman sweater', 'slots': ['torso_base'], 'colour': '#b9b0a0',
     'pack': 'shirts01', 'asset': 'toigo_fisherman_sweater', 'obj': 'sweater_fisherman.obj', 'standoff_mm': 7.0},
    {'id': 'turtleneck', 'name': 'Turtleneck top', 'slots': ['torso_base'], 'colour': '#4c5560',
     'pack': 'shirts01', 'asset': 'toigo_turtleneck_halter_top', 'obj': 'turtleneck_halter.obj', 'standoff_mm': 4.0},
    {'id': 'tank-top', 'name': 'Tank top', 'slots': ['torso_base'], 'colour': '#8fa9a8',
     'pack': 'shirts01', 'asset': 'toigo_keyhole_tank_top', 'obj': 'tank_keyhole_neck.obj', 'standoff_mm': 3.0},
    {'id': 'camisole', 'name': 'Camisole', 'slots': ['torso_base'], 'colour': '#c9b6bd',
     'pack': 'shirts01', 'asset': 'toigo_camisole_top', 'obj': 'camisole_top.obj', 'standoff_mm': 3.0},
    {'id': 'tube-top', 'name': 'Tube top', 'slots': ['torso_base'], 'colour': '#9c6f86',
     'pack': 'shirts01', 'asset': 'skalldyrssuppe_tube_top_funky_colors', 'obj': 'tube_top.obj', 'standoff_mm': 3.0},
    {'id': 'suit-jacket', 'fit_mode': 'drape', 'name': 'Suit jacket', 'slots': ['torso_outer'], 'colour': '#3a3f4a',
     'pack': 'system', 'asset': 'male_elegantsuit01', 'obj': 'male_elegantsuit01.obj', 'components': [0],
     'standoff_mm': 9.0},
    {'id': 'suit-trousers', 'name': 'Suit trousers', 'slots': ['legs'], 'colour': '#3a3f4a',
     'pack': 'system', 'asset': 'male_elegantsuit01', 'obj': 'male_elegantsuit01.obj', 'components': [1],
     'standoff_mm': 6.0},
    {'id': 'casual-shirt', 'fit_mode': 'drape', 'name': 'Casual shirt', 'slots': ['torso_base'], 'colour': '#6d7f96',
     'pack': 'system', 'asset': 'male_casualsuit03', 'obj': 'male_casualsuit03.obj', 'components': [0],
     'standoff_mm': 6.0},
    {'id': 'casual-trousers', 'name': 'Casual trousers', 'slots': ['legs'], 'colour': '#4a5162',
     'pack': 'system', 'asset': 'male_casualsuit03', 'obj': 'male_casualsuit03.obj', 'components': [1],
     'standoff_mm': 6.0},
    {'id': 'wool-trousers', 'name': 'Wool trousers', 'slots': ['legs'], 'colour': '#5a5548',
     'pack': 'pants01', 'asset': 'toigo_wool_pants', 'obj': 'pants_wool.obj', 'standoff_mm': 6.0},
    {'id': 'cargo-trousers', 'fit_mode': 'drape', 'name': 'Cargo trousers', 'slots': ['legs'], 'colour': '#6d6b52',
     'pack': 'pants01', 'asset': 'cortu_cargo_pants', 'obj': 'cargo_pants.obj', 'standoff_mm': 6.0},
    {'id': 'denim-shorts', 'name': 'Denim shorts', 'slots': ['legs'], 'colour': '#3f5f82',
     'pack': 'pants01', 'asset': 'cortu_jeans_shorts', 'obj': 'jean_shorts.obj', 'standoff_mm': 6.0},
    {'id': 'harem-trousers', 'fit_mode': 'drape', 'name': 'Harem trousers', 'slots': ['legs'], 'colour': '#7a5b74',
     'pack': 'pants01', 'asset': 'toigo_harem_pants', 'obj': 'pants_harem.obj', 'standoff_mm': 8.0},
    {'id': 'long-skirt', 'fit_mode': 'drape', 'name': 'Long skirt', 'slots': ['legs'], 'colour': '#6b4a63',
     'pack': 'skirts01', 'asset': 'toigo_long_full_skirt', 'obj': 'skirt_full_long.obj', 'standoff_mm': 8.0},
    {'id': 'mini-skirt', 'fit_mode': 'drape', 'name': 'Mini skirt', 'slots': ['legs'], 'colour': '#8a4f5e',
     'pack': 'skirts01', 'asset': 'frankyaye_mini_skirt_02', 'obj': 'mini_skirt_02.obj', 'standoff_mm': 6.0},
    {'id': 'halter-dress', 'fit_mode': 'drape', 'name': 'Halter dress', 'slots': ['torso_base', 'legs'], 'colour': '#7a4b5c',
     'pack': 'dress01', 'asset': 'toigo_halter_dress_knee_length', 'obj': 'dress_knee_halter.obj', 'standoff_mm': 6.0},
    {'id': 'midi-dress', 'fit_mode': 'drape', 'name': 'Midi dress', 'slots': ['torso_base', 'legs'], 'colour': '#57506e',
     'pack': 'dress01', 'asset': 'toigo_halter_dress_midi', 'obj': 'dress_midi_halter.obj', 'standoff_mm': 6.0},
    {'id': 'tunic', 'fit_mode': 'drape', 'name': 'Tunic', 'slots': ['torso_base', 'legs'], 'colour': '#b3a07f',
     'pack': 'dress01', 'asset': 'wdg_mycenaean_tunic', 'obj': 'mycenaean_tunic.obj', 'standoff_mm': 7.0},
    {'id': 'ankle-socks', 'name': 'Ankle socks', 'slots': ['feet_inner'], 'colour': '#e2e5e9',
     'pack': 'underwear04', 'asset': 'joepal_crude_low_socks', 'obj': 'crudelowsocks.obj', 'standoff_mm': 2.0},
    {'id': 'crew-socks', 'name': 'Crew socks', 'slots': ['feet_inner'], 'colour': '#dadfe6',
     'pack': 'underwear04', 'asset': 'joepal_crude_high_socks', 'obj': 'crudehighsocks.obj', 'standoff_mm': 2.0},
    {'id': 'stockings', 'name': 'Stockings', 'slots': ['feet_inner'], 'colour': '#cdbfb6',
     'pack': 'underwear01', 'asset': 'marco_105_stocking01', 'obj': 'stocking01.obj', 'standoff_mm': 2.0},
    {'id': 'shoes', 'name': 'Shoes', 'slots': ['feet_outer'], 'colour': '#33302e',
     'pack': 'system', 'asset': 'shoes01', 'obj': 'shoes01.obj', 'standoff_mm': 3.0},
    {'id': 'dress-shoes', 'name': 'Dress shoes', 'slots': ['feet_outer'], 'colour': '#241f1c',
     'pack': 'system', 'asset': 'shoes02', 'obj': 'shoes02.obj', 'standoff_mm': 3.0},
    {'id': 'gloves', 'name': 'Gloves', 'slots': ['hands'], 'colour': '#5a4636',
     'pack': 'gloves01', 'asset': 'toigo_gloves_medium', 'obj': 'gloves_medium.obj', 'standoff_mm': 2.0},
    {'id': 'long-gloves', 'name': 'Long gloves', 'slots': ['hands'], 'colour': '#4c3c52',
     'pack': 'gloves01', 'asset': 'toigo_gloves_long', 'obj': 'gloves_long.obj', 'standoff_mm': 2.0},
    {'id': 'newsboy-cap', 'name': 'Newsboy cap', 'slots': ['head'], 'colour': '#5f5a4e',
     'pack': 'hats01', 'asset': 'jujube_newsboy_cap', 'obj': 'newsboy_cap.obj', 'standoff_mm': 4.0},
    {'id': 'cloche-hat', 'name': 'Cloche hat', 'slots': ['head'], 'colour': '#7b4f4f',
     'pack': 'hats01', 'asset': 'aethelraed_unraed_cloche_hat', 'obj': 'cloche_hat.obj', 'standoff_mm': 4.0},
    {'id': 'fedora', 'name': 'Fedora', 'slots': ['head'], 'colour': '#4b4139',
     'pack': 'system', 'asset': 'fedora01', 'obj': 'fedora.obj', 'standoff_mm': 4.0},
]


def slot_model(catalogue=CATALOGUE, slots=SLOTS):
    known = {s['id'] for s in slots}
    members = {}
    for entry in catalogue:
        for slot in entry['slots']:
            if slot not in known:
                raise ValueError(f'{entry["id"]} claims undeclared slot {slot!r}')
            members.setdefault(slot, []).append(entry['id'])
    rows = []
    for entry in catalogue:
        blocked = sorted({other['id'] for other in catalogue if other['id'] != entry['id']
                          and set(other['slots']) & set(entry['slots'])})
        rows.append({'garment': entry['id'], 'occupies': list(entry['slots']), 'excludes': blocked})
    return {
        'schema': 'ihm.garment-slot-model.v1',
        'rule': 'Two garments are mutually exclusive if and only if their occupied slot sets intersect. '
                'Garments with disjoint slot sets are freely combinable. A garment may occupy several slots, '
                'in which case it excludes every garment in any of them.',
        'slots': slots,
        'slot_members': {slot: sorted(ids) for slot, ids in sorted(members.items())},
        'empty_slots': sorted(known - set(members)),
        'garments': rows,
        'draw_order': [s['id'] for s in sorted(slots, key=lambda s: (s['region'], s['layer']))],
        'draw_order_basis': 'ascending declared layer within a body region; a display convention, not a contact result',
    }


def load_envelope(npz_path):
    data = np.load(npz_path, allow_pickle=False)
    positions = np.asarray(data['positions'], float)
    triangles = np.asarray(data['indices'], np.int64)
    area = mesh_area_m2(positions, triangles)
    if area > SLAB_AREA_REJECT_M2:
        raise ValueError(f'target surface area {area:.4f} m2 exceeds the one-sided limit; this is the skin slab')
    height = float(positions[:, 1].max() - positions[:, 1].min())
    if abs(height - REFERENCE_HEIGHT_M) > 1e-3:
        raise ValueError(f'target surface height {height:.6f} m is not the canonical reference height')
    return positions, triangles, area, height
