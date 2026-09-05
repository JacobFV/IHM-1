"""Independent complete provenance checks for the seven-mesh derived result."""
from pathlib import Path
import argparse,gzip,json,tempfile
import numpy as np
from inspect_remaining_muscle_bulk import self_test,paired_disks
from decompose_muscle_dimensions import ROOT,cancellation
from audit_surface_volume_readiness import sha,topology

def verify(directory):
    out=Path(directory);manifest=json.loads((out/'manifest.json').read_text());parent=ROOT/manifest['parent_manifest_path']
    assert sha(parent)==manifest['parent_manifest_sha256']
    for name,digest in manifest['artifacts_sha256'].items():assert sha(out/name)==digest,name
    reports=[];faces=0;components=0
    for path in sorted((out/'entities').glob('*/report.json')):
        r=json.loads(path.read_text());reports.append(r);directory=path.parent
        source=r['source_geometry'];assert sha(ROOT/source['path'])==source['sha256']
        raw=r['raw_source_evidence'];assert sha(ROOT/raw['path'])==raw['sha256']
        canonical=json.loads(gzip.decompress((ROOT/source['path']).read_bytes()));x=np.asarray(canonical['positions'],np.float64).reshape(-1,3);t=np.asarray(canonical['indices']).reshape(-1,3)
        mapping=json.loads((directory/'mapping.json').read_text());bulk=np.asarray(mapping['bulk_source_face_indices'],int);strata=np.asarray(mapping['strata_source_face_indices'],int)
        assert np.array_equal(np.sort(np.r_[bulk,strata]),np.arange(len(t))) and not np.intersect1d(bulk,strata).size
        for name,ids in [('bulk',bulk),('strata',strata)]:
            data=json.loads(gzip.decompress((directory/(name+'.json.gz')).read_bytes()))
            assert np.asarray(data['positions'],np.float64).tobytes()==x.reshape(-1).tobytes()
            assert np.array_equal(np.asarray(data['indices'],int).reshape(-1,3),t[ids])
        assert topology(x,t[bulk])==r['remaining_bulk_topology']
        lines=(ROOT/raw['path']).read_text().splitlines();all_moved=[];pairs=[]
        for disk in mapping['new_paired_disks']:
            ids=np.asarray(disk['canonical_source_face_indices'],int);local=t[ids];selected,recomputed=paired_disks(x,local)
            assert np.array_equal(selected,np.arange(len(ids))) and len(recomputed)==1
            for face,line_number in zip(ids,disk['raw_obj_face_lines']):
                line=lines[line_number-1].split();assert line[0]=='f'
                assert np.array_equal([int(s.split('/')[0])-1 for s in line[1:]],t[face])
            for pair in recomputed[0]['opposite_pairs']:pairs.append({'source_face_indices':ids[pair['source_face_indices']].tolist()})
            all_moved.extend(ids.tolist());components+=1
        assert np.array_equal(np.sort(all_moved),np.sort(mapping['new_disk_source_face_indices']))
        proof=cancellation(x,t,pairs)
        assert proof['maximum_pair_winding_residual']<1e-10 and proof['maximum_scale_normalized_volume_residual']<1e-10
        for proposal in r['unapplied_proposals']:
            assert proposal['applied'] is False and set(proposal['source_face_indices'])<=set(bulk)
        faces+=len(t)
    assert len(reports)==7 and components==9
    result={'passed':True,'entities_verified':len(reports),'source_faces_all_retained':faces,'new_disk_components_verified':components,
            'new_bulk_candidates':sum(r['bulk_candidate'] for r in reports),'whole_corpus_bulk_candidates':418+sum(r['bulk_candidate'] for r in reports),
            'manifest_sha256':sha(out/'manifest.json'),'verifier_sha256':sha(__file__),
            'checks':['artifact_and_source_hashes','original_face_exhaustive_partition','coordinate_bytes','exact_retained_triangle_indices',
                      'original_OBJ_face_lines','entire_paired_disk_support','signed_cancellation','unapplied_proposals_remain_retained','remaining_bulk_reaudit']}
    destination=Path(tempfile.mkdtemp(prefix='remaining-muscle-audit-',dir=ROOT/'artifacts/verification'))
    (destination/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(result,verification_path=str(destination/'verification.json')),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--existing',type=Path);args=parser.parse_args();self_test()
    if args.existing:verify(args.existing)
