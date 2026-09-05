"""Read-only canonical topology audit; closed surfaces are not occupied volumes."""
from pathlib import Path
import argparse,collections,gzip,hashlib,json,os,shutil,tempfile,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parents[1]
UNRESOLVED=('self_intersections_not_checked','vertex_link_manifoldness_not_checked',
            'outward_orientation_and_cavity_convention_not_established',
            'shell_nesting_and_intercomponent_overlap_not_checked','interentity_overlap_and_exclusive_material_ownership_not_checked',
            'physical_tissue_occupancy_and_material_parameters_not_established')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def edges(triangles):
    directed=np.concatenate((triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]))
    unique,inverse,count=np.unique(np.sort(directed,axis=1),axis=0,return_inverse=True,return_counts=True)
    balance=np.bincount(inverse,weights=np.where(directed[:,0]<directed[:,1],1,np.where(directed[:,0]>directed[:,1],-1,0)),minlength=len(unique))
    report={'unique_edges':len(unique),'boundary_edges':int(np.count_nonzero(count==1)),
            'nonmanifold_edges':int(np.count_nonzero(count>2)),'collapsed_edges':int(np.count_nonzero(unique[:,0]==unique[:,1])),
            'two_face_winding_conflicts':int(np.count_nonzero((count==2)&(balance!=0))),
            'edge_incidence_histogram':{str(k):int(v) for k,v in zip(*np.unique(count,return_counts=True))}}
    return report,unique,inverse,count

def topology(vertices,triangles):
    raw=np.asarray(vertices,float);tri=np.asarray(triangles)
    if raw.ndim!=2 or raw.shape[1]!=3 or not len(raw) or not np.isfinite(raw).all():raise ValueError('Invalid or nonfinite vertex coordinates')
    if tri.ndim!=2 or tri.shape[1]!=3 or not len(tri) or tri.dtype.kind not in 'iu' or tri.min()<0 or tri.max()>=len(raw):raise ValueError('Invalid triangular connectivity')
    raw_edges=edges(tri)[0]
    x,reverse=np.unique(raw,axis=0,return_inverse=True);t=reverse[tri]
    edge_report,unique,inverse,count=edges(t)
    used=np.unique(t);unused=len(x)-len(used)
    vgraph=coo_matrix((np.ones(len(unique)),(unique[:,0],unique[:,1])),shape=(len(x),len(x))).tocsr()
    _,vlabels=connected_components(vgraph,directed=False)
    vertex_components=len(np.unique(vlabels[used]))
    # Shared-edge face components use a star per edge, avoiding dense cliques
    # at a nonmanifold edge. No Python adjacency list per source vertex.
    order=np.argsort(inverse,kind='stable');face_index=np.tile(np.arange(len(t)),3)[order]
    starts=np.r_[0,np.cumsum(count)[:-1]];first=np.repeat(face_index[starts],count)
    adjacency=coo_matrix((np.ones(len(face_index)),(first,face_index)),shape=(len(t),len(t))).tocsr()
    ncomponents,labels=connected_components(adjacency,directed=False)
    a,b,c=x[t[:,0]],x[t[:,1]],x[t[:,2]]
    area_twice=np.linalg.norm(np.cross(b-a,c-a),axis=1)
    scale=np.maximum.reduce((np.sum((b-a)**2,axis=1),np.sum((c-a)**2,axis=1),np.sum((c-b)**2,axis=1)))
    repeated=np.any(np.diff(np.sort(t,axis=1),axis=1)==0,axis=1)
    degenerate=area_twice==0
    near=(area_twice<=64*np.finfo(float).eps*scale)&~degenerate
    _,duplicates=np.unique(np.sort(t,axis=1),axis=0,return_counts=True)
    integral=np.einsum('ij,ij->i',a,np.cross(b,c))/6
    components=[]
    face_order=np.argsort(labels,kind='stable');component_counts=np.bincount(labels)
    offset=0
    for component,n in enumerate(component_counts):
        ids=face_order[offset:offset+n];offset+=n;nodes=np.unique(t[ids]);center=(x[nodes].min(axis=0)+x[nodes].max(axis=0))/2
        centered=float(np.einsum('ij,ij->i',a[ids]-center,np.cross(b[ids]-center,c[ids]-center)).sum()/6)
        components.append({'id':component,'faces':int(n),'vertices':len(nodes),
                           'signed_integral_canonical_origin_m3':float(integral[ids].sum()),
                           'signed_integral_local_origin_m3':centered,'local_origin_m':center.tolist()})
    blockers=[]
    for key,value in [('boundary_edges',edge_report['boundary_edges']),('nonmanifold_edges',edge_report['nonmanifold_edges']),
                      ('winding_conflicts',edge_report['two_face_winding_conflicts']),('repeated_index_faces',repeated.sum()),
                      ('zero_area_faces',degenerate.sum()),('near_degenerate_faces',near.sum()),('duplicate_faces',np.sum(duplicates-1))]:
        if value:blockers.append(key)
    return {'raw_vertices':len(raw),'welded_vertices':len(x),'welded_duplicate_vertices':len(raw)-len(x),
            'unused_welded_vertices':unused,'faces':len(t),'raw_edges':raw_edges,'welded_edges':edge_report,
            'repeated_index_faces':int(repeated.sum()),'zero_area_faces':int(degenerate.sum()),
            'near_degenerate_faces':int(near.sum()),'near_degenerate_rule':'twice_area <= 64*float64_epsilon*maximum_squared_edge; exact-zero reported separately',
            'duplicate_faces_ignoring_winding':int(np.sum(duplicates-1)),
            'vertex_connected_components':vertex_components,'edge_connected_face_components':int(ncomponents),
            'components':components,'surface_area_m2':float(area_twice.sum()/2),
            'signed_surface_integral_m3':float(integral.sum()),'occupied_or_union_volume_m3':None,
            'closed_edge_incidence':bool(np.all(count==2)),'consistent_two_face_winding':edge_report['two_face_winding_conflicts']==0,
            'topological_candidate':not blockers,'topology_blockers':blockers,
            'unresolved_volume_checks':list(UNRESOLVED),
            'integral_interpretation':'Signed triangle surface integral, not occupied/union volume; open surfaces are origin-dependent and closed overlapping/nested components require further interpretation.'}

def audit(output):
    out=Path(output)
    if out.exists():raise ValueError('Choose a fresh audit output directory')
    out.mkdir(parents=True);anatomy_path=ROOT/'data/derived/canonical/anatomy.json'
    anatomy_bytes=anatomy_path.read_bytes();anatomy=json.loads(anatomy_bytes);began=time.monotonic()
    (out/'inputs').mkdir();(out/'inputs/anatomy.json').write_bytes(anatomy_bytes);shutil.copyfile(__file__,out/'inputs/audit_surface_volume_readiness.py')
    identity_cache={};mesh_cache={};systems={};total=collections.Counter();reasons=collections.Counter();all_files={}
    def check(reference):
        name=reference.get('path');expected=reference.get('sha256')
        if not isinstance(name,str):return 'missing_source_path'
        path=ROOT/name
        if Path(name).is_absolute() or '..' in Path(name).parts or not path.resolve().is_relative_to(ROOT):return 'unsafe_source_path'
        if name not in identity_cache:
            try:
                before=path.stat();digest=sha(path);after=path.stat()
                if (before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns):raise ValueError('Source changed while hashing')
                identity_cache[name]=digest;all_files[name]=digest
            except (OSError,ValueError) as error:identity_cache[name]={'error':str(error)}
        actual=identity_cache[name]
        if isinstance(actual,dict):return 'missing_or_unstable_source_file'
        if actual!=expected:return 'source_sha256_mismatch'
        return None
    with (out/'sources.jsonl').open('x') as source_log,(out/'entities.jsonl').open('x') as entity_log:
        for index,e in enumerate(anatomy['entities']):
            ref=e.get('reference_geometry',{});failures=[]
            for reference in [ref,*e.get('provenance',{}).get('files',[])]:
                issue=check(reference)
                if issue:failures.append({'path':reference.get('path'),'expected_sha256':reference.get('sha256'),'reason':issue})
            representation=ref.get('representation');digest=ref.get('sha256');metrics=None
            if not failures and representation in ('triangular_surface','surface_shell_quadrature_layer') and ref.get('units')=='m' and ref.get('frame')==anatomy['frame']['id']:
                if digest not in mesh_cache:
                    try:
                        path=ROOT/ref['path'];content=path.read_bytes()
                        if hashlib.sha256(content).hexdigest()!=digest:raise ValueError('Geometry changed after identity check')
                        data=json.loads(gzip.decompress(content) if str(path).endswith('.gz') else content)
                        positions=np.asarray(data.get('positions',data.get('vertices')),float).reshape(-1,3)
                        indices=np.asarray(data.get('indices',data.get('faces')))
                        if indices.size%3:raise ValueError('Nontriangular index count')
                        metrics=topology(positions,indices.reshape(-1,3))
                        mesh_cache[digest]=metrics
                    except (OSError,ValueError,TypeError) as error:mesh_cache[digest]={'error':str(error)}
                    source_log.write(json.dumps({'sha256':digest,'first_reference_path':ref['path'],'topology':mesh_cache[digest]},allow_nan=False)+'\n');source_log.flush()
                metrics=mesh_cache[digest]
            topology_blockers=(metrics.get('topology_blockers',[]) if metrics else [])
            blockers=[f['reason'] for f in failures]+topology_blockers
            if representation!='triangular_surface':blockers.append('reference_is_'+str(representation))
            if ref.get('units')!='m' or ref.get('frame')!=anatomy['frame']['id']:blockers.append('unsupported_units_or_frame')
            if metrics and 'error' in metrics:blockers.append('invalid_geometry')
            candidate=not blockers and metrics is not None and metrics.get('topological_candidate',False)
            record={'entity_id':e['id'],'name':e['name'],'system':e.get('system','unclassified'),
                    'source_id':e.get('source_id'),'evidence_kind':e.get('evidence_kind'),'reference_geometry':ref,
                    'provenance_files':e.get('provenance',{}).get('files',[]),'source_identity_failures':failures,
                    'topology_source_sha256':digest if metrics else None,'topological_candidate':candidate,
                    'raw_vertices':metrics.get('raw_vertices') if metrics else None,
                    'welded_vertices':metrics.get('welded_vertices') if metrics else None,
                    'faces':metrics.get('faces') if metrics else None,
                    'blocking_reasons':sorted(set(blockers)),'remaining_volume_checks':list(UNRESOLVED),
                    'previous_watertight_edge_incidence':e.get('watertight_edge_incidence'),
                    'closed_after_exact_welding':bool(metrics and metrics.get('closed_edge_incidence')),
                    'previously_open_now_closed':bool(e.get('watertight_edge_incidence') is False and metrics and metrics.get('closed_edge_incidence')),
                    'occupied_or_union_volume_m3':None}
            entity_log.write(json.dumps(record,allow_nan=False)+'\n')
            system=systems.setdefault(record['system'],collections.Counter());system['entities']+=1;total['entities']+=1
            for key,value in [('topological_candidates',candidate),('previously_open_now_closed',record['previously_open_now_closed']),
                              ('closed_after_exact_welding',record['closed_after_exact_welding']),('source_identity_failures',bool(failures)),('blocked',bool(blockers))]:
                system[key]+=int(value);total[key]+=int(value)
            reasons.update(record['blocking_reasons'])
            if (index+1)%200==0:print(f'Audited {index+1}/{len(anatomy["entities"])} entities; {total["topological_candidates"]} topological candidates',flush=True)
    if anatomy_path.read_bytes()!=anatomy_bytes:raise ValueError('Canonical anatomy changed during audit')
    summary={'counts':dict(total),'unique_triangle_sources_audited':len(mesh_cache),'source_files_hashed_once':len(identity_cache),
             'unique_mesh_totals':{key:sum(m.get(key,0) for m in mesh_cache.values()) for key in ('raw_vertices','welded_vertices','faces','zero_area_faces','near_degenerate_faces','duplicate_faces_ignoring_winding')},
             'by_system':{k:dict(v) for k,v in sorted(systems.items())},'blocking_reason_counts':dict(reasons),
             'wall_seconds':time.monotonic()-began,'unresolved_checks':list(UNRESOLVED),
             'status':'topological_readiness_audit_not_volume_validation','original_data_modified':False}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    manifest={'schema':'ihm.surface-volume-readiness.v1','anatomy_sha256':hashlib.sha256(anatomy_bytes).hexdigest(),
              'source_files_sha256':all_files,'auditor_sha256':sha(out/'inputs/audit_surface_volume_readiness.py'),
              'artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},
              'source_retention':'Anatomy and audit code copied; geometry/raw inputs remain at recorded original paths with hashes. This is not a detached full-source archive.'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(summary,output_dir=str(out)),indent=2));return summary

def self_test():
    x=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    t=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    closed=topology(x,t);assert closed['topological_candidate'] and np.isclose(closed['signed_surface_integral_m3'],1/6)
    split=topology(x[t].reshape(-1,3),np.arange(12).reshape(-1,3))
    assert split['raw_edges']['boundary_edges']==12 and split['welded_edges']['boundary_edges']==0 and split['welded_duplicate_vertices']==8
    assert topology(x,t[:-1])['welded_edges']['boundary_edges']==3
    flipped=t.copy();flipped[0]=flipped[0,::-1];assert topology(x,flipped)['welded_edges']['two_face_winding_conflicts']==3
    assert topology(x,np.r_[t,t[:1]])['duplicate_faces_ignoring_winding']==1
    assert topology(x,np.r_[t,[[0,0,1]]])['repeated_index_faces']==1
    collinear=np.r_[x,[[.5,0,0]]];assert topology(collinear,np.r_[t,[[0,4,1]]])['zero_area_faces']==1
    two=topology(np.r_[x,x+3],np.r_[t,t+4]);assert two['vertex_connected_components']==two['edge_connected_face_components']==2
    touching=topology(np.r_[x,-x],np.r_[t,t+4]);assert touching['vertex_connected_components']==1 and touching['edge_connected_face_components']==2
    reversed_mesh=topology(x,t[:,::-1]);assert reversed_mesh['topological_candidate'] and reversed_mesh['signed_surface_integral_m3']<0
    assert reversed_mesh['occupied_or_union_volume_m3'] is None and 'self_intersections_not_checked' in reversed_mesh['unresolved_volume_checks']
    print('PASS exact weld, open boundaries, inconsistent winding, degeneracy, duplicates, components and signed-integral qualification')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.self_test:self_test()
    if args.output:audit(args.output)
    elif not args.self_test:parser.error('Supply --output with a fresh directory or --self-test')
