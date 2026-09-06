#!/usr/bin/env python3
"""Bounded retained-byte geometry/provenance audit; no network or native process."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.lymph_registration import SkinTerritoryInventory


def audit(root=ROOT):
    def receipt(path, expected=None):
        raw = (root / path).read_bytes()
        sha = hashlib.sha256(raw).hexdigest()
        if expected is not None and sha != expected:
            raise ValueError('Changed retained source: ' + path)
        return {'path': path, 'sha256': sha, 'bytes': len(raw)}
    ap = 'data/derived/canonical/anatomy.json'
    gp = 'data/derived/canonical/lymphatic_graph.json'
    anatomy = json.loads((root / ap).read_text())
    graph = json.loads((root / gp).read_text())
    skin = [e for e in anatomy['entities'] if e['role'] == 'skin']
    if len(skin) != 1:
        raise ValueError('Explicit skin inventory required for changed canonical assembly')
    skin = skin[0]; ref = skin['reference_geometry']
    inventory = SkinTerritoryInventory(root / ref['path'], ref['sha256']).report()
    groups = []
    for e in anatomy['entities']:
        if e['role'] != 'lymph_node_group':
            continue
        g = e['reference_geometry']
        groups.append({'id': e['id'], 'name': e['name'],
                       'geometry': receipt(g['path'], g['sha256']),
                       'provenance': e['provenance'], 'drainage_identity_validated': False})
    provenance = json.loads((root / 'data/raw/lymphatic/provenance.json').read_text())
    source_receipts = [receipt(x['path'], x['sha256']) for x in provenance['files'] + provenance['archive_members']]
    return {'schema': 'lymph_registration_source_audit_v1',
            'canonical_receipts': [receipt(ap), receipt(gp), receipt(ref['path'], ref['sha256'])],
            'frame': anatomy['frame'], 'skin': {'id': skin['id'], 'provenance': skin['provenance'],
                                               'inventory': inventory},
            'lymph_node_groups': groups, 'publisher_source_receipts': source_receipts,
            'publisher_license': {k: provenance[k] for k in ['license', 'license_evidence', 'third_party_basis']},
            'graph': {'nodes': len(graph['nodes']), 'edges': len(graph['edges']),
                      'directed': graph['directed'],
                      'named_nodes': sum(n.get('anatomical_name') is not None for n in graph['nodes']),
                      'lymph_nodes': sum(n['is_lymph_node'] for n in graph['nodes']),
                      'group_associations': len(graph['node_group_associations']),
                      'accepted_spatial_candidates': sum(x['accepted'] for x in graph['node_group_associations']),
                      'validated_drainage_associations': 0,
                      'validated_skin_triangle_territories': 0},
            'native_engineering_regions_rebound': False,
            'registration_status': 'unresolved_no_source_skin_masks_or_node_identity_crosswalk',
            'native_commands': []}


if __name__ == '__main__':
    target = ROOT / 'data/research/lymph_registration/canonical_audit.json'
    result = audit()
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': str(target.relative_to(ROOT)), 'graph': result['graph'],
                      'skin': result['skin']['inventory']['triangle_count'],
                      'named_node_groups': len(result['lymph_node_groups'])}))
