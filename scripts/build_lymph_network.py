#!/usr/bin/env python3
"""Compile publisher coordinates/topology into an isolated source-frame model."""
import argparse
import gzip
import json
import math
from pathlib import Path
from collect_lymph_network import ROOT, RAW, SOURCES, sha

OUT = ROOT / 'data/derived/lymphatic'
MODEL = 'published-lymphatic-network'
LIMITATIONS = [
    'PlasticBoy-derived authored anatomical estimate; not a measured dataset or independent human subject.',
    'Source graph includes author-modeled connections and output vertices; see 2020 paper section 2.2.',
    'Undirected structural topology: From/To columns retain source ordering, not established lymph flow direction.',
    'Display edges are straight endpoint chords, not lumen surfaces or solved vessel centerlines.',
    'Radii, pressures, flows, valve states and subject calibration are unknown; no transport solver is run.',
    'Separate source specimen/frame; no anatomical registration to BodyParts3D, Z-Anatomy or OpenSim.',
]

def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(value, separators=(',', ':'), allow_nan=False).encode()
    if path.suffix == '.gz':
        with path.open('wb') as f:
            with gzip.GzipFile(fileobj=f, mode='wb', mtime=0) as z:
                z.write(content)
    else:
        path.write_bytes(content)


def rows(path):
    return [line.split() for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]


def build():
    provenance = json.loads((RAW / 'provenance.json').read_text())
    for name, (_, expected) in SOURCES.items():
        if sha(RAW / name) != expected:
            raise ValueError(f'Source identity changed: {name}')
    for member in provenance['archive_members']:
        if sha(ROOT / member['path']) != member['sha256']:
            raise ValueError('Extracted source identity changed')
    directory = RAW / 'extracted/supplementary files/Data-based HLS graph'
    nodes = []
    for i, row in enumerate(rows(directory / 'graph_vertices.txt'), 1):
        if len(row) != 4 or row[3] not in ('0', '1'):
            raise ValueError('Invalid source vertex')
        p = list(map(float, row[:3]))
        if not all(math.isfinite(x) for x in p):
            raise ValueError('Nonfinite coordinates')
        nodes.append({'id': i, 'position_mm': p, 'is_lymph_node': row[3] == '1', 'anatomical_name': None})
    edges = []
    for i, row in enumerate(rows(directory / 'graph_edges.txt'), 1):
        if len(row) != 3:
            raise ValueError('Invalid source edge')
        a, b = map(int, row[:2]); length = float(row[2])
        if not (1 <= a <= len(nodes) and 1 <= b <= len(nodes)) or a == b or not math.isfinite(length) or length <= 0:
            raise ValueError('Invalid source edge reference/length')
        chord = math.dist(nodes[a-1]['position_mm'], nodes[b-1]['position_mm'])
        edges.append({'id': i, 'source_from': a, 'source_to': b, 'source_length_mm': length,
                      'endpoint_chord_length_mm': chord, 'radius_mm': None, 'flow_ml_s': None,
                      'flow_direction': None})
    if (len(nodes), len(edges), sum(n['is_lymph_node'] for n in nodes)) != (996, 1117, 272):
        raise ValueError('Unexpected publisher graph counts')
    graph = {'schema_version': 1, 'model_id': MODEL, 'directed': False,
             'edge_order_interpretation': 'source From/To column order retained; no physiological flow direction asserted',
             'source_frame': 'savinkov2020-xyz-mm', 'coordinate_units': 'mm',
             'evidence_kind': provenance['evidence_kind'], 'measured_dataset': False,
             'nodes': nodes, 'edges': edges, 'provenance': provenance, 'limitations': LIMITATIONS}
    dump(OUT / 'graph.json', graph)
    # Rigid display rotation maps source Z (height) to screen Y. Unit conversion only.
    rotated = [[n['position_mm'][0], n['position_mm'][2], -n['position_mm'][1]] for n in nodes]
    low = [min(p[i] for p in rotated) for i in range(3)]
    high = [max(p[i] for p in rotated) for i in range(3)]
    center = [(a+b)/2 for a, b in zip(low, high)]
    positions = [[(p[i]-center[i])/1000 for i in range(3)] for p in rotated]
    line_positions = [c for e in edges for v in (e['source_from'], e['source_to']) for c in positions[v-1]]
    payload = {'positions': line_positions, 'units': 'm', 'primitive': 'line_segments', 'source_edge_ids': [e['id'] for e in edges],
               'source_graph_sha256': sha(OUT / 'graph.json'), 'path_interpretation': LIMITATIONS[3],
               'radius_mm': None, 'flow_ml_s': None, 'directed': False}
    geometry_path = OUT / 'geometry' / (MODEL + '.json.gz')
    dump(geometry_path, payload)
    source = {'label': 'Savinkov et al. 2020 · published lymphatic graph',
              'url': 'https://doi.org/10.3390/math8122236', 'sha256': sha(RAW / 'mathematics-08-02236-s001.zip'),
              'source_revision': provenance['source_revision'], 'units': 'mm', 'frame': graph['source_frame'],
              'specimen': 'PlasticBoy-derived anatomical estimate; no measured subject', 'status': 'published structural model; uncalibrated',
              'license': provenance['license'], 'evidence_kind': provenance['evidence_kind']}
    model = {'id': MODEL, 'name': 'Published lymphatic network · anatomical estimate',
             'description': '996 source vertices and 1,117 structural edges; 272 lymph-node flags. PlasticBoy-derived model, not measured anatomy. No flow or radius data.',
             'frame': MODEL + '-display-m', 'source_frame': graph['source_frame'], 'source_units': 'mm', 'display_units': 'm',
             'calibration_status': 'published model geometry; not measured or subject calibrated',
             'bounds': {'min': [(x-center[i])/1000 for i, x in enumerate(low)], 'max': [(x-center[i])/1000 for i, x in enumerate(high)]},
             'display_transform': {'rotation': [[1,0,0],[0,0,1],[0,-1,0]], 'scale': .001, 'translation': [-x/1000 for x in center]},
             'source': source, 'attribution': [provenance['attribution']], 'limitations': LIMITATIONS,
             'graph_path': str((OUT / 'graph.json').relative_to(ROOT)), 'graph_sha256': sha(OUT / 'graph.json')}
    structure = {'id': MODEL, 'model_id': MODEL, 'name': 'Lymphatic vessels · published structural graph',
                 'system': 'lymphatic', 'kind': 'lines', 'geometry_url': '/api/geometry/' + MODEL,
                 'geometry_sha256': sha(geometry_path), 'color': '#7bc58e', 'default_visible': True,
                 'source': source, 'calibration_status': model['calibration_status'], 'limitations': LIMITATIONS,
                 'source_vertex_count': len(nodes), 'source_edge_count': len(edges), 'lymph_node_flag_count': 272,
                 'directed': False, 'flow_ml_s': None, 'radius_mm': None}
    fragment = {'models': [model], 'structures': [structure]}
    dump(OUT / 'manifest_fragment.json', fragment)
    print(json.dumps({'vertices': len(nodes), 'edges': len(edges), 'lymph_node_flags': 272, 'manifest_fragment': str(OUT / 'manifest_fragment.json')}))
    return fragment


def append_to_manifest(manifest_path=ROOT / 'data/derived/app/manifest.json'):
    """Explicit integration hook; build() does not edit the shared application manifest."""
    import shutil
    fragment = json.loads((OUT / 'manifest_fragment.json').read_text())
    source_geometry=OUT/'geometry'/(MODEL+'.json.gz')
    if sha(source_geometry)!=fragment['structures'][0]['geometry_sha256'] or sha(OUT/'graph.json')!=fragment['models'][0]['graph_sha256']:
        raise ValueError('Lymphatic graph/display identity changed before integration')
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    for key in ('models', 'structures'):
        manifest[key] = [x for x in manifest[key] if x.get('model_id', x['id']) != MODEL] + fragment[key]
    destination = manifest_path.parent / 'geometry' / (MODEL + '.json.gz')
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(OUT / 'geometry' / destination.name, destination)
    temporary = manifest_path.with_suffix('.lymphatic.tmp')
    dump(temporary, manifest)
    temporary.replace(manifest_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--append', action='store_true', help='explicitly integrate built model into app manifest')
    args = parser.parse_args()
    build()
    if args.append:
        append_to_manifest()
