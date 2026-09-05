"""Conservative derived muscle surfaces; opposite windings remain unresolved."""
from pathlib import Path
import argparse,gzip,hashlib,json,shutil,time
import numpy as np
from audit_surface_volume_readiness import ROOT,sha,topology

def clean(vertices,triangles):
    x=np.asarray(vertices,dtype=np.float64);t=np.asarray(triangles)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():raise ValueError('Finite 3D coordinates required')
    if t.ndim!=2 or t.shape[1]!=3 or not len(t) or t.dtype.kind not in 'iu' or t.min()<0 or t.max()>=len(x):raise ValueError('Valid triangles required')
    # Byte-level coordinate identity is stricter than proximity and numeric
    # equality: even distinct +0/-0 bit patterns are not silently coalesced.
    _,vertex_ids=np.unique(np.ascontiguousarray(x).view(np.dtype((np.void,24))).ravel(),return_inverse=True)
    keys=vertex_ids[t];start=np.argmin(keys,axis=1)
    cyclic=np.take_along_axis(keys,(np.arange(3)[None,:]+start[:,None])%3,axis=1)
    zero=np.all(np.cross(x[t[:,1]]-x[t[:,0]],x[t[:,2]]-x[t[:,0]])==0,axis=1)
    groups={}
    for i,key in enumerate(np.sort(keys,axis=1)):
        if not zero[i]:groups.setdefault(tuple(key),[]).append(i)
    retained=[];duplicate_of={};ambiguous=[]
    for ids in groups.values():
        orientations={tuple(cyclic[i]) for i in ids}
        if len(orientations)>1:
            ambiguous.append(ids);retained.extend(ids)
        else:
            retained.append(ids[0]);duplicate_of.update({i:ids[0] for i in ids[1:]})
    retained=np.array(sorted(retained),dtype=np.int64)
    mapping=np.full(len(t),-1,dtype=np.int64);mapping[retained]=np.arange(len(retained))
    for source,owner in duplicate_of.items():mapping[source]=mapping[owner]
    return t[retained].copy(),{'source_face_to_retained_face':mapping.tolist(),
              'retained_source_face_indices':retained.tolist(),'excluded_zero_area_source_faces':np.flatnonzero(zero).tolist(),
              'removed_duplicate_source_faces':sorted(duplicate_of),'opposite_winding_source_face_groups':ambiguous,
              'coordinate_identity':'float64 coordinate bytes; cyclic orientation preserved; no tolerance welding',
              'coordinates_changed':False,'opposite_winding_policy':'retain entire ambiguous group'}

def self_test():
    x=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    t=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    duplicated=np.r_[x,x];faces=np.r_[t,np.roll(t[0]+4,1)[None,:]]
    result,mapping=clean(duplicated,faces)
    assert np.array_equal(result,t) and mapping['source_face_to_retained_face']==[0,1,2,3,0]
    assert topology(duplicated,result)['topological_candidate']
    # Exactly the same oriented coordinate triangle set, without multiplicity.
    def oriented_set(v,f):
        keys=[]
        for face in v[f]:
            rotations=[np.roll(face,i,axis=0).tobytes() for i in range(3)];keys.append(min(rotations))
        return set(keys)
    assert oriented_set(duplicated,faces)==oriented_set(duplicated,result)
    ambiguous=np.r_[faces,t[0,::-1][None,:]];unchanged,m=clean(duplicated,ambiguous)
    assert np.array_equal(unchanged,ambiguous) and len(m['opposite_winding_source_face_groups'])==1
    assert not m['removed_duplicate_source_faces'] and not topology(duplicated,unchanged)['topological_candidate']
    result,m=clean(x,np.r_[t,[[0,0,1]]]);assert np.array_equal(result,t) and m['source_face_to_retained_face'][-1]==-1
    near=x.copy();near[3]=[1e-14,0,0];_,m=clean(near,np.r_[t,[[0,1,3]]]);assert 4 in m['excluded_zero_area_source_faces']
    # Distinct coordinate bytes cannot become a duplicate through rounding.
    shifted=np.r_[x,x.copy()];shifted[5,0]=np.nextafter(1.,2.)
    _,m=clean(shifted,np.r_[t,(t[0]+4)[None,:]]);assert not m['removed_duplicate_source_faces']
    slender=np.array([[0.,0,0],[1.,0,0],[0,1e-170,0]])
    retained,m=clean(slender,np.array([[0,1,2]]));assert len(retained)==1,'A small nonzero cross product is not exact zero area'
    print('PASS cyclic coordinate-duplicate equivalence, opposite-winding ambiguity, explicit zero-area mapping and precision retention')

def build(output):
    out=Path(output)
    if out.exists():raise ValueError('Choose a fresh derived cleanup directory')
    out.mkdir(parents=True);(out/'inputs').mkdir();(out/'geometry').mkdir();(out/'mappings').mkdir()
    anatomy_path=ROOT/'data/derived/canonical/anatomy.json';original=anatomy_path.read_bytes();anatomy=json.loads(original)
    (out/'inputs/anatomy.json').write_bytes(original)
    for p in (Path(__file__),Path(__file__).with_name('audit_surface_volume_readiness.py')):shutil.copyfile(p,out/'inputs'/p.name)
    cache={};records=[];started=time.monotonic();source_hashes={}
    for entity in anatomy['entities']:
        if entity.get('system')!='muscular':continue
        ref=entity['reference_geometry'];path=ROOT/ref['path']
        if ref['representation']!='triangular_surface' or ref['units']!='m':raise ValueError('Unsupported muscle reference')
        data_bytes=path.read_bytes()
        if hashlib.sha256(data_bytes).hexdigest()!=ref['sha256']:raise ValueError('Canonical muscle geometry identity changed')
        source_hashes[ref['path']]=ref['sha256']
        for source in entity['provenance']['files']:
            if source['path'] not in source_hashes:source_hashes[source['path']]=sha(ROOT/source['path'])
            if source_hashes[source['path']]!=source['sha256']:raise ValueError('Muscle raw source identity changed')
        digest=ref['sha256']
        if digest not in cache:
            data=json.loads(gzip.decompress(data_bytes));x=np.asarray(data['positions'],float).reshape(-1,3);t=np.asarray(data['indices']).reshape(-1,3)
            derived,mapping=clean(x,t);before=topology(x,t);after=topology(x,derived)
            data['indices']=derived.reshape(-1).tolist()
            data['display_faces']=len(derived)
            data['derived_surface_cleanup']={'source_sha256':digest,'source_faces':len(t),'retained_faces':len(derived),
                                             'vertex_positions_changed':False,'vertex_normals':'retained from source',
                                             'physical_volume_validated':False}
            geometry=out/'geometry'/(digest+'.json.gz');geometry.write_bytes(gzip.compress(json.dumps(data,separators=(',',':'),allow_nan=False).encode(),mtime=0))
            restored=json.loads(gzip.decompress(geometry.read_bytes()))
            assert np.asarray(restored['positions'],np.float64).tobytes()==x.reshape(-1).tobytes()
            mapping_path=out/'mappings'/(digest+'.json');mapping_path.write_text(json.dumps(mapping,separators=(',',':'))+'\n')
            cache[digest]={'source_sha256':digest,'derived_geometry':str(geometry.relative_to(out)),
                           'derived_geometry_sha256':sha(geometry),'mapping':str(mapping_path.relative_to(out)),
                           'mapping_sha256':sha(mapping_path),'before':before,'after':after,
                           'removed_duplicate_faces':len(mapping['removed_duplicate_source_faces']),
                           'removed_zero_area_faces':len(mapping['excluded_zero_area_source_faces']),
                           'opposite_winding_groups':len(mapping['opposite_winding_source_face_groups'])}
        m=cache[digest];records.append({'entity_id':entity['id'],'name':entity['name'],'reference_geometry':ref,
                 'source_files':entity['provenance']['files'],'source_record_sha256':digest,
                 'before_candidate':m['before']['topological_candidate'],'after_candidate':m['after']['topological_candidate'],
                 'became_candidate':not m['before']['topological_candidate'] and m['after']['topological_candidate'],
                 'remaining_topology_blockers':m['after']['topology_blockers'],
                 'removed_duplicate_faces':m['removed_duplicate_faces'],'removed_zero_area_faces':m['removed_zero_area_faces'],
                 'opposite_winding_groups':m['opposite_winding_groups'],'canonical_activation_applied':False})
        if len(records)%100==0:print('Audited derived muscles',len(records),flush=True)
    if anatomy_path.read_bytes()!=original:raise ValueError('Canonical anatomy changed during cleanup')
    for name,values in [('sources',cache.values()),('entities',records)]:
        with (out/(name+'.jsonl')).open('x') as stream:
            for value in values:stream.write(json.dumps(value,allow_nan=False)+'\n')
    summary={'muscles':len(records),'unique_sources':len(cache),
             **{key:sum(int(r[key]) for r in records) for key in ('before_candidate','after_candidate','became_candidate','removed_duplicate_faces','removed_zero_area_faces','opposite_winding_groups')},
             'still_blocked':sum(not r['after_candidate'] for r in records),'wall_seconds':time.monotonic()-started,
             'status':'derived_geometry_only; opposite-winding ambiguity retained; no physical volume or activation claim'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    manifest={'schema':'ihm.muscle-duplicate-cleanup.v1','anatomy_sha256':hashlib.sha256(original).hexdigest(),
              'source_files_sha256':source_hashes,'artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},
              'retention':'Source anatomy and code copied, full coordinate arrays retained in derived geometry, all face mappings retained; raw source files remain at original hashed paths.',
              'original_files_modified':False}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(summary,output_dir=str(out)),indent=2));return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.self_test:self_test()
    if args.output:build(args.output)
    elif not args.self_test:parser.error('Supply --output or --self-test')
