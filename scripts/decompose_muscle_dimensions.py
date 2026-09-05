"""Preserve isolated opposite-triangle components as explicit surface strata."""
from pathlib import Path
import argparse,gzip,hashlib,json,math,shutil,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from audit_surface_volume_readiness import ROOT,sha,topology

def components(vertices,triangles):
    x=np.asarray(vertices,np.float64);t=np.asarray(triangles)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():raise ValueError('Finite vertex coordinates required')
    if t.ndim!=2 or t.shape[1]!=3 or not len(t) or t.dtype.kind not in 'iu' or t.min()<0 or t.max()>=len(x):raise ValueError('Valid triangle connectivity required')
    _,inverse=np.unique(np.ascontiguousarray(x).view(np.dtype((np.void,24))).ravel(),return_inverse=True)
    welded=inverse[t];edge=np.concatenate((welded[:,[0,1]],welded[:,[1,2]],welded[:,[2,0]]))
    _,ids,count=np.unique(np.sort(edge,axis=1),axis=0,return_inverse=True,return_counts=True)
    order=np.argsort(ids,kind='stable');face_ids=np.tile(np.arange(len(t)),3)[order]
    starts=np.r_[0,np.cumsum(count)[:-1]];first=np.repeat(face_ids[starts],count)
    graph=coo_matrix((np.ones(len(face_ids)),(first,face_ids)),shape=(len(t),len(t))).tocsr()
    n,labels=connected_components(graph,directed=False)
    return welded,labels,np.bincount(labels,minlength=n)

def split_strata(vertices,triangles):
    x=np.asarray(vertices,np.float64);t=np.asarray(triangles);welded,labels,counts=components(x,t)
    order=np.argsort(labels,kind='stable');starts=np.r_[0,np.cumsum(counts)[:-1]];pairs=[]
    for component in np.flatnonzero(counts==2):
        ids=order[starts[component]:starts[component]+2];a,b=welded[ids]
        if len(set(a))!=3 or not np.array_equal(np.sort(a),np.sort(b)):continue
        a=np.roll(a,-int(np.argmin(a)));b=np.roll(b,-int(np.argmin(b)))
        if not np.array_equal(a,b[[0,2,1]]):continue
        tri=x[t[ids[0]]]
        if np.all(np.cross(tri[1]-tri[0],tri[2]-tri[0])==0):continue
        pairs.append({'edge_component_id':int(component),'source_face_indices':sorted(ids.tolist())})
    sheet=np.array(sorted(i for p in pairs for i in p['source_face_indices']),dtype=np.int64)
    bulk=np.setdiff1d(np.arange(len(t)),sheet,assume_unique=True)
    return bulk,sheet,pairs

def solid_angle(point,triangle):
    a,b,c=np.asarray(triangle,float)-point;la,lb,lc=(np.linalg.norm(v) for v in (a,b,c))
    numerator=float(a@np.cross(b,c))
    denominator=la*lb*lc+(a@b)*lc+(b@c)*la+(c@a)*lb
    return float(2*np.arctan2(numerator,denominator))

def cancellation(vertices,triangles,pairs):
    x=np.asarray(vertices,float);t=np.asarray(triangles);max_winding=0.;max_volume=0.;max_naive=0.;max_relative_volume=0.;checks=0
    for pair in pairs:
        a,b=x[t[pair['source_face_indices']]];normal=np.cross(a[1]-a[0],a[2]-a[0]);normal/=np.linalg.norm(normal)
        center=a.mean(axis=0);scale=max(np.linalg.norm(a[1]-a[0]),np.linalg.norm(a[2]-a[0]),np.linalg.norm(a[2]-a[1]))
        # Normal offsets guarantee samples off the source triangle, on both sides.
        for distance in (.1,-.1,1.,-1.,3.):
            point=center+distance*scale*normal
            residual=abs(solid_angle(point,a)+solid_angle(point,b))/(4*np.pi)
            max_winding=max(max_winding,residual);checks+=1
        for origin in (np.zeros(3),center,center+scale*np.array([2.,3.,5.])):
            aa=a-origin;bb=b-origin
            naive_a=float(aa[0]@np.cross(aa[1],aa[2]))/6;naive_b=float(bb[0]@np.cross(bb[1],bb[2]))/6
            max_naive=max(max_naive,abs(naive_a+naive_b))
            # Equivalent determinant with local edges avoids subtracting large
            # products when a tiny surface lies far from the canonical origin.
            va=float(aa[0]@np.cross(aa[1]-aa[0],aa[2]-aa[0]))/6
            vb=float(bb[0]@np.cross(bb[1]-bb[0],bb[2]-bb[0]))/6
            residual=abs(va+vb);max_volume=max(max_volume,residual)
            max_relative_volume=max(max_relative_volume,residual/max(abs(va),abs(vb),scale**3,np.finfo(float).tiny))
    return {'classified_pairs':len(pairs),'off_surface_winding_samples':checks,
            'maximum_pair_winding_residual':max_winding,'maximum_pair_signed_volume_residual_m3':max_volume,
            'maximum_naive_global_product_volume_residual_m3':max_naive,
            'maximum_scale_normalized_volume_residual':max_relative_volume,
            'test_points':'centroid plus normal offsets of ±0.1, ±1 and +3 maximum edge lengths',
            'volume_evaluation':'anchor dot cross(local edges)/6; algebraically identical to triple absolute coordinates, with less cancellation error',
            'qualification':'Numerical cancellation checks at off-surface samples support the analytic orientation-reversal identity; this is not a self-intersection or physical-thickness test.'}

def self_test():
    x=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    tetra=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    sheet=np.array([[3.,0,0],[4.,0,0],[3.,1.,0]])
    vertices=np.r_[x,sheet,sheet];faces=np.r_[tetra,[[4,5,6],[7,9,8]]]
    bulk,side,pairs=split_strata(vertices,faces)
    assert np.array_equal(bulk,np.arange(4)) and np.array_equal(side,[4,5])
    assert topology(vertices,faces[bulk])['topological_candidate']
    check=cancellation(vertices,faces,pairs);assert check['maximum_pair_winding_residual']<1e-12 and check['maximum_scale_normalized_volume_residual']<1e-12
    # A pair sharing an edge with a larger component must stay unresolved there.
    _,side,pairs=split_strata(vertices,np.r_[faces,[[4,5,0]]]);assert not len(side) and not pairs
    _,side,pairs=split_strata(sheet,np.array([[0,1,2],[1,2,0]]));assert not len(side)
    bulk,side,pairs=split_strata(sheet,np.array([[0,1,2],[0,2,1]]));assert not len(bulk) and len(side)==2
    translated=np.array([[-.111778,-.846894,.0901013],[-.111766,-.8468797,.0900903],[-.111760,-.8468868,.0901043]])
    faces=np.array([[0,1,2],[2,1,0]]);_,_,pairs=split_strata(translated,faces)
    assert cancellation(translated,faces,pairs)['maximum_scale_normalized_volume_residual']<1e-10,'Volume arithmetic must handle small triangles far from the origin'
    print('PASS isolated-pair separation, larger-component retention, same-winding rejection, pair-only no-bulk and signed cancellation')

def build(output):
    out=Path(output)
    if out.exists():raise ValueError('Choose a fresh dimensional decomposition directory')
    out.mkdir(parents=True)
    for name in ('inputs','bulk','strata','mappings'):(out/name).mkdir()
    anatomy_path=ROOT/'data/derived/canonical/anatomy.json';original=anatomy_path.read_bytes();anatomy=json.loads(original)
    (out/'inputs/anatomy.json').write_bytes(original)
    for path in (Path(__file__),Path(__file__).with_name('audit_surface_volume_readiness.py')):shutil.copyfile(path,out/'inputs'/path.name)
    sources={};records=[];hashes={};started=time.monotonic()
    for entity in anatomy['entities']:
        if entity.get('system')!='muscular':continue
        ref=entity['reference_geometry'];path=ROOT/ref['path'];content=path.read_bytes();digest=hashlib.sha256(content).hexdigest()
        if digest!=ref['sha256'] or ref['units']!='m' or ref['representation']!='triangular_surface':raise ValueError('Unsupported or changed canonical muscle source')
        hashes[ref['path']]=digest
        for source in entity['provenance']['files']:
            if source['path'] not in hashes:hashes[source['path']]=sha(ROOT/source['path'])
            if hashes[source['path']]!=source['sha256']:raise ValueError('Changed raw source identity')
        if digest not in sources:
            data=json.loads(gzip.decompress(content));x=np.asarray(data['positions'],float).reshape(-1,3);t=np.asarray(data['indices']).reshape(-1,3)
            bulk,sheet,pairs=split_strata(x,t);before=topology(x,t);after=topology(x,t[bulk]) if len(bulk) else None
            check=cancellation(x,t,pairs)
            if check['maximum_pair_winding_residual']>1e-10 or check['maximum_scale_normalized_volume_residual']>1e-10:raise ValueError('Classified pair fails numerical cancellation check')
            mapping={'source_face_count':len(t),'bulk_source_face_indices':bulk.tolist(),'strata_source_face_indices':sheet.tolist(),
                     'isolated_opposite_pairs':pairs,'strata_role':'unresolved zero-thickness surface; no assigned material, thickness or mass',
                     'coordinates_changed':False,'faces_discarded':0,'canonical_activation_applied':False}
            (out/'mappings'/(digest+'.json')).write_text(json.dumps(mapping,separators=(',',':'))+'\n')
            representations={}
            for name,ids in [('bulk',bulk),('strata',sheet)]:
                derived=dict(data);derived['indices']=t[ids].reshape(-1).tolist();derived['display_faces']=len(ids)
                derived['dimensional_representation']={'source_sha256':digest,'role':name,'source_faces':len(t),'retained_faces':len(ids),
                                                      'physical_volume_or_mass_assigned':False,'original_coordinate_array_retained':True}
                target=out/name/(digest+'.json.gz');target.write_bytes(gzip.compress(json.dumps(derived,separators=(',',':'),allow_nan=False).encode(),mtime=0))
                representations[name]={'path':str(target.relative_to(out)),'sha256':sha(target),'faces':len(ids)}
            sources[digest]={'source_sha256':digest,'representations':representations,'mapping':f'mappings/{digest}.json',
                             'before':before,'remaining_bulk_topology':after,'signed_cancellation':check,
                             'bulk_candidate':bool(after and after['topological_candidate']),
                             'classified_strata_pairs':len(pairs),'no_bulk_faces':not len(bulk)}
        source=sources[digest];records.append({'entity_id':entity['id'],'name':entity['name'],'reference_geometry':ref,
                   'source_files':entity['provenance']['files'],'source_record_sha256':digest,
                   'before_candidate':source['before']['topological_candidate'],'bulk_candidate':source['bulk_candidate'],
                   'has_unresolved_strata':source['classified_strata_pairs']>0,'strata_pairs':source['classified_strata_pairs'],
                   'no_bulk_faces':source['no_bulk_faces'],
                   'remaining_bulk_blockers':source['remaining_bulk_topology']['topology_blockers'] if source['remaining_bulk_topology'] else ['no_volumetric_bulk'],
                   'canonical_activation_applied':False,'mass_assigned_kg':None})
        if len(records)%100==0:print('Decomposed muscle sources',len(records),flush=True)
    if anatomy_path.read_bytes()!=original:raise ValueError('Canonical anatomy changed during decomposition')
    for name,values in [('sources',sources.values()),('entities',records)]:
        with (out/(name+'.jsonl')).open('x') as stream:
            for value in values:stream.write(json.dumps(value,allow_nan=False)+'\n')
    summary={'muscles':len(records),'unique_sources':len(sources),'before_candidates':sum(r['before_candidate'] for r in records),
             'bulk_candidates':sum(r['bulk_candidate'] for r in records),
             'bulk_candidates_with_unresolved_strata':sum(r['bulk_candidate'] and r['has_unresolved_strata'] for r in records),
             'new_bulk_candidates':sum(r['bulk_candidate'] and not r['before_candidate'] for r in records),
             'pair_only_no_bulk_muscles':sum(r['no_bulk_faces'] for r in records),'strata_pairs':sum(r['strata_pairs'] for r in records),
             'still_blocked_bulk':sum(not r['bulk_candidate'] and not r['no_bulk_faces'] for r in records),
             'maximum_pair_winding_residual':max(s['signed_cancellation']['maximum_pair_winding_residual'] for s in sources.values()),
             'maximum_pair_signed_volume_residual_m3':max(s['signed_cancellation']['maximum_pair_signed_volume_residual_m3'] for s in sources.values()),
             'maximum_scale_normalized_volume_residual':max(s['signed_cancellation']['maximum_scale_normalized_volume_residual'] for s in sources.values()),
             'faces_discarded':0,'wall_seconds':time.monotonic()-started,
             'qualification':'Geometric dimensional decomposition, not physical deletion, inferred thickness, occupied volume, material assignment or canonical activation.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    manifest={'schema':'ihm.muscle-dimensional-decomposition.v1','anatomy_sha256':hashlib.sha256(original).hexdigest(),
              'source_files_sha256':hashes,'artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},
              'original_data_modified':False,'mass_claims_created':False,'source_retention':'Original source hashes plus copied anatomy/code; both derived coordinate arrays and disjoint face mappings retain every canonical face.'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(summary,output_dir=str(out)),indent=2));return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.self_test:self_test()
    if args.output:build(args.output)
    elif not args.self_test:parser.error('Supply --output or --self-test')
