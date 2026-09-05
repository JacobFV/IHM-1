"""Read-only verification of every retained bulk/stratum face and coordinate."""
from pathlib import Path
import argparse,gzip,hashlib,json,tempfile
import numpy as np
from decompose_muscle_dimensions import components,cancellation,self_test,ROOT
from audit_surface_volume_readiness import sha

def verify(directory):
    out=Path(directory);manifest=json.loads((out/'manifest.json').read_text())
    for name,digest in manifest['artifacts_sha256'].items():
        if sha(out/name)!=digest:raise ValueError('Changed decomposition artifact: '+name)
    sources={r['source_sha256']:r for r in map(json.loads,(out/'sources.jsonl').read_text().splitlines())}
    entities=list(map(json.loads,(out/'entities.jsonl').read_text().splitlines()));checked=set();total_faces=0;total_pairs=0
    max_winding=0.;max_volume=0.
    for entity in entities:
        digest=entity['source_record_sha256']
        if digest in checked:continue
        checked.add(digest);source=sources[digest];path=ROOT/entity['reference_geometry']['path']
        if sha(path)!=digest:raise ValueError('Original canonical source changed; use its retained matching version')
        raw=json.loads(gzip.decompress(path.read_bytes()));x=np.asarray(raw['positions'],np.float64).reshape(-1,3);t=np.asarray(raw['indices']).reshape(-1,3)
        mapping=json.loads((out/source['mapping']).read_text());bulk=np.asarray(mapping['bulk_source_face_indices'],int);side=np.asarray(mapping['strata_source_face_indices'],int)
        assert np.array_equal(np.sort(np.r_[bulk,side]),np.arange(len(t)))
        assert not np.intersect1d(bulk,side).size and mapping['faces_discarded']==0
        for name,ids in [('bulk',bulk),('strata',side)]:
            geometry=json.loads(gzip.decompress((out/source['representations'][name]['path']).read_bytes()))
            assert np.asarray(geometry['positions'],np.float64).tobytes()==x.reshape(-1).tobytes()
            assert np.array_equal(np.asarray(geometry['indices'],int).reshape(-1,3),t[ids])
        _,labels,count=components(x,t);pairs=mapping['isolated_opposite_pairs']
        assert np.array_equal(np.sort([i for pair in pairs for i in pair['source_face_indices']]),side)
        for pair in pairs:
            i,j=pair['source_face_indices'];assert labels[i]==labels[j] and count[labels[i]]==2
            a,b=x[t[i]],x[t[j]]
            assert any(a.tobytes()==np.roll(b[::-1],k,axis=0).tobytes() for k in range(3))
        proof=cancellation(x,t,pairs)
        assert proof==source['signed_cancellation']
        assert proof['maximum_pair_winding_residual']<1e-10 and proof['maximum_scale_normalized_volume_residual']<1e-10
        max_winding=max(max_winding,proof['maximum_pair_winding_residual']);max_volume=max(max_volume,proof['maximum_pair_signed_volume_residual_m3'])
        assert source['no_bulk_faces']==(not len(bulk))
        if not len(bulk):assert source['remaining_bulk_topology'] is None and not source['bulk_candidate']
        else:assert source['remaining_bulk_topology']['faces']==len(bulk)
        total_faces+=len(t);total_pairs+=len(pairs)
    report={'passed':True,'manifest_sha256':sha(out/'manifest.json'),'entities':len(entities),'unique_sources':len(checked),
            'source_faces_all_accounted_for':total_faces,'classified_pairs_checked':total_pairs,
            'maximum_pair_winding_residual':max_winding,'maximum_pair_signed_volume_residual_m3':max_volume,
            'checks':['artifact_hashes','original_source_identity','byte_identical_coordinate_arrays','disjoint_exhaustive_source_face_partition',
                      'exact_retained_triangle_indices','isolated_two_face_components_only','opposite_cyclic_coordinate_identity',
                      'off_surface_winding_and_signed_volume_cancellation','pair_only_has_no_bulk'],
            'verifier_sha256':sha(__file__),'qualification':'No physical thickness, mass, occupied volume or canonical activation established.'}
    destination=Path(tempfile.mkdtemp(prefix='muscle-dimensions-audit-',dir=ROOT/'artifacts/verification'))
    (destination/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(report,verification_path=str(destination/'verification.json')),indent=2));return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--existing',type=Path);args=parser.parse_args();self_test()
    if args.existing:verify(args.existing)
