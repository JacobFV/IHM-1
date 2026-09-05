#!/usr/bin/env python3
"""Verify provenance, source geometry, display geometry and source-frame transform."""
import gzip
import json
from pathlib import Path
import tempfile
import numpy as np
from build_extended_anatomy import ROOT, RAW, OUT, MODEL_ID, sha256, append_manifest, validate_provenance, validate_cached_sources, validate_fragment
from collect_extended_anatomy import REVISION, EXPECTED_SHA256

def verify():
    provenance=validate_provenance()
    assert provenance['revision']==REVISION
    for record in provenance['files']:
        assert sha256(ROOT/record['path'])==record['sha256']==EXPECTED_SHA256[Path(record['path']).name]
        assert (ROOT/record['path']).stat().st_size==record['bytes']
    assert sha256(ROOT/provenance['blend_path'])==provenance['blend_sha256']
    index=validate_cached_sources(provenance)
    fragment=json.loads((OUT/'manifest_fragment.json').read_text())
    coverage=json.loads((OUT/'coverage.json').read_text())
    validate_fragment(fragment,provenance,index)
    assert index['source_blend_sha256']==provenance['blend_sha256']
    assert len(index['meshes'])==len(fragment['structures'])==coverage['structure_count']
    assert len({s['id'] for s in fragment['structures']})==len(fragment['structures'])
    model=fragment['models'][0];assert model['id']==MODEL_ID
    assert not model['cross_family_registration'] and model['independent_subject_count']==0
    t=model['display_transform'];rotation=np.asarray(t['rotation']);scale=t['scale'];translation=np.asarray(t['translation'])
    global_min=np.full(3,np.inf);global_max=-global_min
    node_objects=[];curve_count=0
    for entry,structure in zip(index['meshes'],fragment['structures']):
        assert entry['id']==structure['id'] and structure['source']['object_name']==entry['name']
        assert '7: Nervous system & Sense organs' not in entry['source_collections']
        assert not entry['name'].endswith(('.j','.g','.t','.i')) and 'kidney' not in entry['name'].lower()
        assert sha256(ROOT/entry['source_geometry_path'])==entry['source_geometry_sha256']
        with np.load(ROOT/entry['source_geometry_path'],allow_pickle=False) as d:
            v=d['vertices'];f=d['faces'];local=d['local_vertices'];local_faces=d['local_faces']
        matrix=np.array(entry['matrix_world'])
        np.testing.assert_allclose(v,local@matrix[:3,:3].T+matrix[:3,3],atol=1e-12)
        np.testing.assert_array_equal(f,local_faces[:,::-1] if np.linalg.det(matrix[:3,:3])<0 else local_faces)
        assert sha256(ROOT/entry['base_geometry_path'])==entry['base_geometry_sha256']
        with np.load(ROOT/entry['base_geometry_path'],allow_pickle=False) as base:
            assert len(base['faces'])==entry['base_triangles'] and len(base['vertices'])==entry['base_vertices']
            original_matrix=np.array(entry['original_matrix_world'])
            np.testing.assert_allclose(base['vertices'],base['local_vertices']@original_matrix[:3,:3].T+original_matrix[:3,3],atol=1e-12)
            np.testing.assert_array_equal(base['faces'],base['local_faces'][:,::-1] if np.linalg.det(original_matrix[:3,:3])<0 else base['local_faces'])
        assert v.shape==(entry['source_vertices'],3) and f.shape==(entry['source_triangles'],3)
        assert np.isfinite(v).all() and f.min()>=0 and f.max()<len(v)
        np.testing.assert_allclose([v.min(0),v.max(0)],entry['bounds'])
        transformed=v@rotation.T*scale+translation
        global_min=np.minimum(global_min,transformed.min(0));global_max=np.maximum(global_max,transformed.max(0))
        path=ROOT/structure['geometry_path'];assert sha256(path)==structure['geometry_sha256']
        with gzip.open(path,'rt') as stream:p=json.load(stream)
        positions=np.array(p['positions']).reshape(-1,3);faces=np.array(p['indices']).reshape(-1,3)
        normals=np.array(p['normals']).reshape(-1,3)
        assert len(positions)>0 and len(faces)>0 and positions.shape==normals.shape
        assert np.isfinite(positions).all() and np.isfinite(normals).all()
        assert faces.min()>=0 and faces.max()<len(positions)
        assert not p['display_decimation'] and p['full_resolution']
        assert p['original_faces']==p['source_triangles']==p['output_triangles']==len(f)==len(faces)
        np.testing.assert_array_equal(faces,f)
        assert structure['source']['license']=='CC-BY-SA-4.0'
        if not p['display_decimation']:
            np.testing.assert_allclose(positions,transformed,atol=1e-8)
        if structure['system']=='lymphatic' and 'node' in entry['name'].lower():node_objects.append(entry['name'])
        curve_count+=entry['source_type']=='CURVE'
    np.testing.assert_allclose(global_min,model['bounds']['min'],atol=1e-7)
    np.testing.assert_allclose(global_max,model['bounds']['max'],atol=1e-7)
    assert len(node_objects)==coverage['lymph_node_group_objects'] and len(node_objects)>100
    assert coverage['lymphatic_vessel_objects']==0
    assert coverage['source_triangles']==coverage['output_triangles']==sum(e['source_triangles'] for e in index['meshes'])
    assert len(json.loads((OUT/'raw_object_inventory.json').read_text())['objects'])==7184
    # Source-authored procedural nodes and mirrored gland must actually be evaluated.
    named={e['name']:e for e in index['meshes']}
    assert named['Median sacral nodes']['base_triangles']==0 < named['Median sacral nodes']['source_triangles']
    for name in ['Thyroid gland','Laryngopharynx']:
        e=named[name]
        assert any(m['type']=='MIRROR' for m in e['evaluated_modifiers'])
        assert e['source_triangles']>e['base_triangles']
        assert np.ptp(np.array(e['bounds']),axis=0)[0]>np.ptp(np.array(e['base_bounds']),axis=0)[0]*1.5
    for term in ['axillary','inguinal','jugular','mesenteric']:
        assert any(term in n.lower() for n in node_objects),term
    # Verify integration is idempotent and preserves unrelated records without touching real manifest.
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp)/'manifest.json'
        sentinel={'id':'existing','model_id':'existing'}
        path.write_text(json.dumps({'models':[{'id':'existing'}],'structures':[sentinel],'limitations':[]}))
        append_manifest(path);append_manifest(path)
        merged=json.loads(path.read_text())
        assert merged['structures'][0]==sentinel and len(merged['structures'])==len(fragment['structures'])+1
        assert len(merged['models'])==2
    print(json.dumps({'status':'passed','source_and_display_objects':len(fragment['structures']),
        'lymph_node_group_objects':len(node_objects),'tessellated_source_curves':curve_count,
        'counts_by_system':coverage['counts_by_system']},indent=2))

if __name__=='__main__':verify()
