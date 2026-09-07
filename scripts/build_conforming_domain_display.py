"""Export per-entity boundary surfaces of the conforming tetrahedral domains for the workbench.

Interior tetrahedra are not renderable. What a viewer needs is, for every material owner in a
conforming volume, the faces of that owner's tets that no other tet of the same owner shares:
the owner's own closed boundary inside the shared vertex array. Extracting it changes no
volume, moves no vertex, and invents no label -- every triangle here is a face of a tet the
source domain already wrote, and every owner name/role is copied from the domain's own member
record.

The exported display keeps the domain's single vertex array, so interfaces stay coincident and
the owners still meet exactly where the mesher made them meet.

  .venv/bin/python scripts/build_conforming_domain_display.py --out data/derived/conforming-domain-display-v1
"""
from pathlib import Path
import argparse, gzip, hashlib, json
import numpy as np

BASE = Path(__file__).resolve().parents[1]
FACES = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def verify(relative, expected):
    path = BASE / relative
    digest = sha(path)
    if digest != expected:
        raise ValueError(f'{relative} changed on disk ({digest} != {expected}); rebuild the domain')
    return digest


def owner_boundary(tets, labels):
    """Per-label boundary triangles: the label's tet faces incident to exactly one of its tets.

    Vertex indices stay in the domain's own numbering, so the shared nodes stay shared.
    """
    surfaces = {}
    for label in np.unique(labels):
        block = tets[labels == label]
        faces = block[:, FACES].reshape(-1, 3)
        keys = np.sort(faces, axis=1)
        _, first, counts = np.unique(keys, axis=0, return_index=True, return_counts=True)
        surfaces[int(label)] = faces[first[counts == 1]]
    return surfaces


def pack(vertices, surfaces, structures):
    """One shared position array plus one index run per role; ranges locate each structure."""
    used = np.unique(np.concatenate([surfaces[s['label']].reshape(-1) for s in structures]))
    remap = np.full(len(vertices), -1, np.int64)
    remap[used] = np.arange(len(used))
    positions = np.asarray(vertices[used], np.float32)
    quantization = float(np.max(np.abs(vertices[used] - positions.astype(np.float64))))
    roles, order = {}, {}
    for structure in structures:
        roles.setdefault(structure['role'], []).append(structure)
    groups, records = [], []
    for role in sorted(roles):
        indices, ranges = [], []
        offset = 0
        for structure in sorted(roles[role], key=lambda s: -len(surfaces[s['label']])):
            triangles = remap[surfaces[structure['label']]]
            indices.append(triangles.reshape(-1))
            ranges.append({'entity_id': structure['entity_id'], 'name': structure['name'],
                           'system': structure.get('system'), 'tets': int(structure['tets']),
                           'triangle_start': offset, 'triangle_count': len(triangles)})
            offset += len(triangles)
            order[structure['entity_id']] = role
        flat = np.concatenate(indices) if indices else np.zeros(0, np.int64)
        groups.append({'role': role, 'triangle_count': offset,
                       'triangles': flat.astype(np.int32).tolist(), 'structures': ranges})
        records.append(offset)
    return positions, groups, quantization, sum(records)


def hand_domain(source='data/derived/hand-reflexive-grip-domain-v1', resolution='coarse'):
    directory = Path(source)
    manifest = json.loads((BASE / directory / 'manifest.json').read_bytes())
    inputs = {str(directory / name): verify(directory / name, digest)
              for name, digest in manifest['artifacts_sha256'].items()
              if name in ('tetmesh-coarse.npz', 'members.jsonl', 'domain.json')}
    inputs[str(directory / 'manifest.json')] = sha(BASE / directory / 'manifest.json')
    domain = json.loads((BASE / directory / 'domain.json').read_bytes())[
        'resolutions'][resolution]
    members = [json.loads(line) for line in (BASE / directory / 'members.jsonl').read_text().splitlines()]
    if len(members) != domain['member_count']:
        raise ValueError('Member roll does not match the meshed resolution')
    volume = np.load(BASE / directory / f'tetmesh-{resolution}.npz', allow_pickle=False)
    vertices, tets, owner = volume['TV'], volume['TT'], volume['owner']
    if domain['conforming_mesh']['tets'] != len(tets) or domain['conforming_mesh']['tet_vertices'] != len(vertices):
        raise ValueError('Recorded conforming mesh size does not match the stored volume')
    counted = {row['entity_id']: row['owned_tets'] for row in domain['conforming_mesh']['per_structure']}
    surfaces = owner_boundary(tets, owner)
    structures = []
    for index, member in enumerate(members):
        if index not in surfaces:
            continue
        tet_count = int((owner == index).sum())
        if counted.get(member['entity_id']) != tet_count:
            raise ValueError('Owner column disagrees with the domain per-structure tet census')
        structures.append({'label': index, 'entity_id': member['entity_id'], 'name': member['name'],
                           'role': member['role'], 'system': member['system'], 'tets': tet_count})
    unclaimed = int((owner == -1).sum())
    if -1 in surfaces and unclaimed:
        if domain['conforming_mesh']['complement_tets'] != unclaimed:
            raise ValueError('Complement tet count disagrees with the domain record')
        structures.append({'label': -1, 'entity_id': 'unclaimed-complement', 'role': 'unclaimed_complement',
                           'name': 'unclaimed complement (inside the meshed box, outside every member)',
                           'system': 'none', 'tets': unclaimed})
    accept = domain['acceptance']
    payload = {
        'id': f'hand-reflexive-grip-{resolution}', 'kind': 'conforming_tetrahedral_domain',
        'title': 'Left hand · conforming tetrahedral domain',
        'summary': (f"{domain['conforming_mesh']['tets']:,} tetrahedra, {len(structures)} material owners, "
                    f"one shared vertex array of {domain['conforming_mesh']['tet_vertices']:,} nodes."),
        'frame': {'id': 'bodyparts3d-display-m', 'units': 'm'},
        'mesher': {'name': domain['mesher'], 'epsilon_relative': domain['epsilon_relative'],
                   'edge_length_relative': domain['edge_length_relative'],
                   'seconds': domain['ftetwild']['seconds']},
        'volume': {'tets': domain['conforming_mesh']['tets'], 'vertices': domain['conforming_mesh']['tet_vertices'],
                   'structure_tets': domain['conforming_mesh']['structure_tets'],
                   'complement_tets': domain['conforming_mesh']['complement_tets'],
                   'total_volume_m3': domain['conforming_mesh']['total_volume_m3'],
                   'min_dihedral_deg': accept['min_dihedral_deg'],
                   'interface_faces_between_two_structures': domain['shared_node_audit']['interface_faces_between_two_structures'],
                   'shared_nodes_at_interfaces': domain['shared_node_audit']['shared_nodes_at_interfaces']},
        'limitations': [
            'Boundary extraction only. Interior tetrahedra are not shipped and cannot be inspected here.',
            f"fTetWild displaced the interfaces: maximum surface deviation {accept['max_interface_deviation_m']:.3e} m "
            f"against a {accept['surface_deviation_budget_m']:.2e} m budget.",
            f"Per-structure volume error reaches {accept['max_abs_per_structure_volume_error']:.3f} against a "
            f"{accept['volume_error_budget']} budget; all {domain['member_count']} members exceed it.",
            'The body-envelope owner surrounds the tissue owners, so it hides them until it is switched off.',
            'Colour encodes the recorded material role only. No mechanical state is displayed.'],
        'source_manifest': str(directory / 'manifest.json'),
    }
    return payload, vertices, surfaces, structures, inputs


def whole_body_domain(spacing='0.01m'):
    directory = Path('data/derived/material-domains') / f'whole-body-{spacing}'
    manifest = json.loads((BASE / directory / 'manifest.json').read_bytes())
    inputs = {str(directory / name): verify(directory / name, digest)
              for name, digest in manifest['artifacts_sha256'].items()}
    inputs[str(directory / 'manifest.json')] = sha(BASE / directory / 'manifest.json')
    rows = {}
    for line in (BASE / directory / 'sources.jsonl').read_text().splitlines():
        row = json.loads(line)
        rows[row['entity_id']] = row
    volume = np.load(BASE / directory / 'whole-body-domain.npz', allow_pickle=False)
    vertices, tets, labels = volume['vertices_m'], volume['tetrahedra'], volume['material_index']
    ids = volume['source_ids']
    if manifest['mesh']['tetrahedra'] != len(tets) or manifest['mesh']['nodes'] != len(vertices):
        raise ValueError('Recorded partition mesh size does not match the stored volume')
    surfaces = owner_boundary(tets, labels)
    tets_per_cell = len(tets) // len(volume['cell_indices'])
    structures = []
    for label in sorted(surfaces):
        entity = str(ids[label])
        row = rows[entity]
        count = int((labels == label).sum())
        if row['owned_cells'] * tets_per_cell != count:
            raise ValueError('Material column disagrees with the recorded owned-cell census')
        structures.append({'label': label, 'entity_id': entity, 'name': row['name'],
                           'role': row['role'], 'system': row['system'], 'tets': count})
    payload = {
        'id': f'whole-body-{spacing}', 'kind': 'voxel_material_partition',
        'title': f'Whole body · material partition ({float(volume["spacing_m"]) * 1000:g} mm)',
        'summary': (f"{len(tets):,} tetrahedra over {len(volume['cell_indices']):,} occupied cells, "
                    f"{len(structures)} owning entities of {len(ids)} sources."),
        'frame': manifest['frame'],
        'mesher': {'name': 'voxel partition + 6-tet cell split', 'spacing_m': float(volume['spacing_m']),
                   'seconds': manifest['wall_seconds']},
        'volume': {'tets': len(tets), 'vertices': len(vertices), 'cells': len(volume['cell_indices']),
                   'entities_claiming_no_cell': manifest['entities_claiming_no_cell'],
                   'entities_owning_no_cell': manifest['entities_owning_no_cell'],
                   'winding_trees': manifest['winding_trees']},
        'limitations': [
            'Boundary extraction only. Interior tetrahedra are not shipped and cannot be inspected here.',
            f"Cell spacing is {float(volume['spacing_m']) * 1000:g} mm; every surface here is that voxel staircase, "
            'not an anatomical boundary.',
            f"{manifest['entities_owning_no_cell']} of {len(ids)} source entities own no cell at this spacing "
            'and are absent from the display.',
            'Ownership is a priority-ranked winding-number partition, not a mechanical assembly.'],
        'source_manifest': str(directory / 'manifest.json'),
    }
    return payload, vertices, surfaces, structures, inputs


BUILDERS = {'hand-reflexive-grip-coarse': hand_domain,
            'whole-body-0.01m': lambda: whole_body_domain('0.01m'),
            'whole-body-0.005m': lambda: whole_body_domain('0.005m')}


def build(out, domains):
    out = Path(out) if Path(out).is_absolute() else BASE / out
    out.mkdir(parents=True, exist_ok=True)
    index = {'schema': 'ihm.conforming-domain-display.v1',
             'builder': 'scripts/build_conforming_domain_display.py',
             'builder_sha256': sha(__file__),
             'projection': ('per-owner boundary faces of the source tetrahedral volume; shared vertex array kept; '
                            'float32 display positions; interior tetrahedra omitted'),
             'domains': []}
    for name in domains:
        payload, vertices, surfaces, structures, inputs = BUILDERS[name]()
        positions, groups, quantization, triangles = pack(vertices, surfaces, structures)
        if quantization > 1e-6:
            raise ValueError('Display quantization exceeds 1 micrometre')
        lo, hi = positions.min(axis=0), positions.max(axis=0)
        payload.update({
            'schema_version': 1,
            'anchor': {'origin_m': ((lo + hi) / 2).tolist(),
                       'extent_m': (hi - lo).tolist(),
                       'local_axes': {'0': [1, 0, 0], '1': [0, 1, 0], '2': [0, 0, 1]}},
            'geometry': {'positions_m': positions.reshape(-1).tolist(),
                         'vertex_count': len(positions), 'triangle_count': triangles,
                         'position_quantization_max_error_m': quantization},
            'groups': groups,
            'structure_count': len(structures),
            'frames': [{'time_s': 0.0}],
            'source_hashes': inputs})
        path = out / f'{name}.json.gz'
        path.write_bytes(gzip.compress(json.dumps(payload, separators=(',', ':'), allow_nan=False).encode(), mtime=0))
        index['domains'].append({
            'id': payload['id'], 'title': payload['title'], 'summary': payload['summary'],
            'display_path': path.name, 'display_sha256': sha(path), 'bytes': path.stat().st_size,
            'structure_count': len(structures), 'triangle_count': triangles,
            'vertex_count': len(positions), 'source_tets': payload['volume']['tets'],
            'roles': [g['role'] for g in groups], 'source_hashes': inputs})
    (out / 'index.json').write_text(json.dumps(index, indent=2) + '\n')
    return index


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='data/derived/conforming-domain-display-v1')
    parser.add_argument('--domains', default='hand-reflexive-grip-coarse,whole-body-0.01m')
    args = parser.parse_args()
    print(json.dumps(build(args.out, [d for d in args.domains.split(',') if d]), indent=2))
