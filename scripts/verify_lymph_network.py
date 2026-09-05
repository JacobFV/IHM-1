#!/usr/bin/env python3
"""Independently verify graph extraction and rendered endpoints against the pinned ZIP."""
import gzip
import json
import math
import zipfile
from build_lymph_network import ROOT, RAW, OUT, MODEL, sha
from collect_lymph_network import SOURCES


def verify():
    archive = RAW / 'mathematics-08-02236-s001.zip'
    assert sha(archive) == SOURCES[archive.name][1]
    with zipfile.ZipFile(archive) as z:
        def read(name):
            return [r.split() for r in z.read('supplementary files/Data-based HLS graph/' + name).decode().splitlines() if r and not r.startswith('#')]
        vertices = read('graph_vertices.txt'); edges = read('graph_edges.txt')
    graph_path = OUT / 'graph.json'
    graph = json.loads(graph_path.read_text())
    assert len(vertices) == len(graph['nodes']) == 996
    assert len(edges) == len(graph['edges']) == 1117
    assert sum(n['is_lymph_node'] for n in graph['nodes']) == 272
    assert graph['directed'] is False and graph['measured_dataset'] is False
    neighbors = {i: set() for i in range(1, 997)}
    unique_edges = set()
    max_length_error = 0
    for i, (row, node) in enumerate(zip(vertices, graph['nodes']), 1):
        assert node['id'] == i and node['position_mm'] == list(map(float, row[:3]))
        assert node['is_lymph_node'] == bool(int(row[3]))
    for i, (row, edge) in enumerate(zip(edges, graph['edges']), 1):
        a, b = map(int, row[:2]); length = float(row[2])
        assert (edge['id'], edge['source_from'], edge['source_to'], edge['source_length_mm']) == (i, a, b, length)
        assert all(edge[k] is None for k in ('radius_mm', 'flow_ml_s', 'flow_direction'))
        neighbors[a].add(b); neighbors[b].add(a); unique_edges.add(tuple(sorted((a,b))))
        max_length_error = max(max_length_error, abs(length - math.dist(list(map(float, vertices[a-1][:3])), list(map(float, vertices[b-1][:3])))))
    assert len(unique_edges) == 1117
    # Rounded source lengths agree with independent Euclidean endpoint calculation.
    assert max_length_error < .001
    visited = set(); stack = [1]
    while stack:
        v = stack.pop()
        if v in visited: continue
        visited.add(v); stack.extend(neighbors[v] - visited)
    assert len(visited) == 996
    geometry_path = OUT / 'geometry' / (MODEL + '.json.gz')
    geometry = json.loads(gzip.decompress(geometry_path.read_bytes()))
    fragment = json.loads((OUT / 'manifest_fragment.json').read_text())
    assert fragment['structures'][0]['geometry_sha256'] == sha(geometry_path)
    assert geometry['source_graph_sha256'] == fragment['models'][0]['graph_sha256'] == sha(graph_path)
    positions = geometry['positions']
    assert len(positions) == 1117 * 6 and geometry['source_edge_ids'] == list(range(1,1118))
    assert geometry['directed'] is False
    t = fragment['models'][0]['display_transform']
    assert t['scale'] == .001 and t['rotation'] == [[1,0,0],[0,0,1],[0,-1,0]]
    for i, row in enumerate(edges):
        for endpoint, source_id in enumerate(map(int, row[:2])):
            p = list(map(float, vertices[source_id-1][:3]))
            offset = i*6 + endpoint*3
            for axis in range(3):
                expected = sum(t['rotation'][axis][j]*p[j] for j in range(3))*.001+t['translation'][axis]
                assert math.isclose(positions[offset+axis], expected, abs_tol=1e-12)
    import tempfile
    from pathlib import Path
    from unittest.mock import patch
    import build_lymph_network as builder
    with tempfile.TemporaryDirectory() as td:
        fixture=Path(td);(fixture/'geometry').mkdir()
        (fixture/'manifest_fragment.json').write_text(json.dumps(fragment))
        (fixture/'geometry'/(MODEL+'.json.gz')).write_bytes(b'changed geometry')
        (fixture/'graph.json').write_bytes(graph_path.read_bytes())
        destination=fixture/'app.json';destination.write_text('{"models":[],"structures":[]}')
        before=destination.read_bytes()
        with patch.object(builder,'OUT',fixture):
            try:builder.append_to_manifest(destination)
            except ValueError:pass
            else:raise AssertionError('changed lymphatic geometry accepted')
        assert destination.read_bytes()==before
    report = {'passed': True, 'source_vertices': 996, 'source_edges': 1117, 'lymph_node_flags': 272,
              'connected_components': 1, 'rendered_segments': 1117,
              'max_source_length_vs_chord_error_mm': max_length_error,
              'source_archive_sha256': sha(archive), 'graph_sha256': sha(graph_path),
              'checks': ['pinned publisher archive', 'every source coordinate/flag/edge/length preserved',
                         'full source graph connected without duplicate undirected edges', 'no assumed radius/flow/direction',
                         'all rendered endpoints match rigid source-to-display transform', 'artifact hashes'],
              'measurement_validation': False, 'flow_validation': False}
    (OUT / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return report

if __name__ == '__main__':
    verify()
