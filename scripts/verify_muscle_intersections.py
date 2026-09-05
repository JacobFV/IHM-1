"""Verify retained intersection evidence and original face identities, no GUI."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import numpy as np


def verify(directory):
    directory=Path(directory).resolve();manifest=json.loads((directory/'manifest.json').read_text())
    assert manifest['schema']=='ihm.muscle-intersection-audit.v1'
    assert manifest['source_geometry_modified'] is False
    for name,digest in manifest['artifacts_sha256'].items():
        path=(directory/name).resolve()
        assert path.is_relative_to(directory) and hashlib.sha256(path.read_bytes()).hexdigest()==digest
    records=[json.loads(line) for line in (directory/'records.jsonl').read_text().splitlines()]
    assert len({r['entity_id'] for r in records})==len(records)==425
    examined=passed=intersecting=0
    for record in records:
        base=directory/'entities'/record['entity_id']
        assert json.loads((base/'report.json').read_text())==record
        assert hashlib.sha256((base/'input-bulk.json.gz').read_bytes()).hexdigest()==record['input_geometry_sha256']
        assert hashlib.sha256((base/'input-mapping.json').read_bytes()).hexdigest()==record['input_mapping_sha256']
        if not record['prior_bulk_candidate']:
            assert not record['geometric_screen_passed'] and not (base/'intersection-evidence.npz').exists()
            continue
        examined+=1
        geometry=json.loads(gzip.decompress((base/'input-bulk.json.gz').read_bytes()))
        x=np.asarray(geometry['positions'],float).reshape(-1,3);f=np.asarray(geometry['indices'],int).reshape(-1,3)
        mapping=json.loads((base/'input-mapping.json').read_text())
        original=np.asarray(mapping['bulk_source_face_indices'])
        with np.load(base/'intersection-evidence.npz',allow_pickle=False) as data:
            pairs=data['intersecting_bulk_face_pairs'];inverse=data['source_to_welded_vertex'];welded=data['welded_vertices_m']
            assert np.array_equal(x,welded[inverse])
            assert pairs.ndim==2 and pairs.shape[1]==2
            assert not len(pairs) or (pairs.min()>=0 and pairs.max()<len(f) and np.all(pairs[:,0]<pairs[:,1]))
            assert len(pairs)==len(np.unique(pairs,axis=0))==record['intersection_pair_count']
            assert np.array_equal(original[pairs],data['intersecting_original_face_pairs'])
            assert len(data['nonmanifold_welded_vertices'])==record['nonmanifold_vertex_count']
            assert record['self_intersection_free']==(len(pairs)==0)
            assert record['vertex_links_manifold']==(record['nonmanifold_vertex_count']==0)
            assert record['geometric_screen_passed']==(record['self_intersection_free'] and record['vertex_links_manifold'])
        assert record['mass_assigned_kg'] is None and record['canonical_handoff_applied'] is False
        passed+=record['geometric_screen_passed'];intersecting+=not record['self_intersection_free']
    summary=json.loads((directory/'summary.json').read_text())
    assert (examined,passed,intersecting)==(summary['prior_bulk_candidates'],summary['screened_candidates'],summary['self_intersecting_entities'])
    return dict(passed=True,entities=len(records),examined=examined,intersection_free_candidates=passed,self_intersecting=intersecting,
                claim='Artifact integrity, exact coordinate welding and complete original face-pair mapping; CGAL numerical detection has a separate native fixture check')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args()
    result=verify(a.directory);print(json.dumps(result,indent=2))
