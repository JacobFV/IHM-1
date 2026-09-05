"""Source-preserving disk-stratum separation for the remaining seven muscles."""
from pathlib import Path
import argparse,gzip,hashlib,json,shutil,tempfile,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from decompose_muscle_dimensions import components,cancellation
from audit_surface_volume_readiness import ROOT,sha,topology

def graph_connected(edges):
    vertices,inverse=np.unique(edges,return_inverse=True);pairs=inverse.reshape(-1,2)
    graph=coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(len(vertices),len(vertices))).tocsr()
    return connected_components(graph,directed=False,return_labels=False)==1

def disk_support(triangles):
    t=np.asarray(triangles,int);v=np.unique(t)
    all_edges=np.concatenate((t[:,[0,1]],t[:,[1,2]],t[:,[2,0]]));edge,count=np.unique(np.sort(all_edges,axis=1),axis=0,return_counts=True)
    boundary=edge[count==1];euler=len(v)-len(edge)+len(t)
    if np.any(count>2) or not len(boundary) or euler!=1:return False
    vertices,degree=np.unique(boundary,return_counts=True)
    if np.any(degree!=2) or not graph_connected(boundary):return False
    boundary_set=set(vertices)
    for vertex in v:
        incident=t[np.any(t==vertex,axis=1)];link=incident[incident!=vertex].reshape(-1,2)
        _,degree=np.unique(link,return_counts=True)
        if not graph_connected(link) or np.any(degree>2):return False
        if vertex in boundary_set:
            if np.count_nonzero(degree==1)!=2:return False
        elif np.any(degree!=2):return False
    return True

def paired_disks(vertices,triangles):
    x=np.asarray(vertices,float);t=np.asarray(triangles);w,labels,counts=components(x,t);selected=[];reports=[]
    for component in np.flatnonzero((counts>2)&(counts%2==0)):
        ids=np.flatnonzero(labels==component);_,inverse,n=np.unique(np.sort(w[ids],axis=1),axis=0,return_inverse=True,return_counts=True)
        if not np.all(n==2):continue
        pairs=[];support=[];valid=True
        for key in range(len(n)):
            a,b=ids[inverse==key];wa,wb=w[a],w[b]
            if len(set(wa))!=3 or not any(np.array_equal(wa,np.roll(wb[::-1],k)) for k in range(3)):
                valid=False;break
            if np.all(np.cross(x[t[a,1]]-x[t[a,0]],x[t[a,2]]-x[t[a,0]])==0):valid=False;break
            pairs.append({'source_face_indices':[int(a),int(b)]});support.append(wa)
        if not valid or not disk_support(np.asarray(support)):continue
        nodes=np.unique(t[ids]);coordinates=np.unique(x[nodes],axis=0)
        reports.append({'parent_bulk_edge_component_id':int(component),'parent_bulk_face_indices':ids.tolist(),
                        'opposite_pairs':pairs,'unoriented_support_faces':len(support),'coordinate_vertices':len(coordinates),
                        'coordinate_extent_m':np.ptp(coordinates,axis=0).tolist(),
                        'center_m':coordinates.mean(axis=0).tolist(),
                        'singular_values_m':np.linalg.svd(coordinates-coordinates.mean(axis=0),compute_uv=False).tolist(),
                        'classification':'entire all-paired component with manifold-disk support; physical sheet role/thickness unresolved'})
        selected.extend(ids.tolist())
    return np.array(sorted(selected),int),reports

def self_test():
    x=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],[.5,.5,0],[0,0,1]])
    fan=np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
    doubled=np.r_[fan,fan[:,::-1]];selected,record=paired_disks(x,doubled)
    assert np.array_equal(selected,np.arange(8)) and len(record)==1
    tetra=np.array([[1,2,5],[0,5,2],[0,1,5],[0,2,1]])
    selected,_=paired_disks(x,np.r_[tetra,tetra[:,::-1]]);assert not len(selected),'A doubled closed shell is not an open disk stratum'
    selected,_=paired_disks(x,np.r_[doubled,[[0,1,5]]]);assert not len(selected),'Do not extract pairs from unpaired bulk'
    assert not disk_support(np.r_[fan,[[0,4,5]]]),'Nonmanifold vertex/edge support must fail'
    print('PASS paired-disk retention, doubled-closed-shell rejection and connected unpaired bulk rejection')

def raw_obj_evidence(entity,canonical):
    source=entity['provenance']['files'][0];path=ROOT/source['path']
    if sha(path)!=source['sha256']:raise ValueError('Changed raw OBJ source')
    vertices=[];faces=[];lines=[];header=[]
    for line_number,line in enumerate(path.read_text().splitlines(),1):
        parts=line.split()
        if not parts:continue
        if parts[0]=='#' and len(header)<12:header.append(line)
        elif parts[0]=='v':vertices.append(list(map(float,parts[1:4])))
        elif parts[0]=='f':
            if len(parts)!=4:raise ValueError('Raw OBJ contains nontriangle faces requiring a separate triangulation correspondence')
            faces.append([int(z.split('/')[0])-1 for z in parts[1:]]);lines.append(line_number)
    transform=entity['provenance']['source_to_canonical']
    transformed=np.asarray(vertices)@np.asarray(transform['rotation']).T*transform['scale']+np.asarray(transform['translation'])
    if not np.array_equal(transformed,np.asarray(canonical['positions']).reshape(-1,3)) or not np.array_equal(faces,np.asarray(canonical['indices']).reshape(-1,3)):
        raise ValueError('Raw OBJ correspondence is not exact; do not infer exporter provenance')
    return {'path':source['path'],'sha256':source['sha256'],'vertex_transform_exact':True,'face_order_exact':True,
            'source_header':header,'face_line_numbers':lines,
            'canonical_license_label':entity['provenance'].get('license'),
            'source_role':'authored atlas, not anatomical imaging or microscopic validation; duplicate geometry already exists in raw OBJ'}

def build(output):
    out=Path(output)
    if out.exists():raise ValueError('Choose a fresh remaining-muscle result directory')
    out.mkdir(parents=True);(out/'inputs').mkdir();(out/'entities').mkdir()
    parent=ROOT/'data/derived/muscle-dimensional-decomposition-v2';parent_manifest=json.loads((parent/'manifest.json').read_text())
    parent_sources={s['source_sha256']:s for s in map(json.loads,(parent/'sources.jsonl').read_text().splitlines())}
    entity_records=list(map(json.loads,(parent/'entities.jsonl').read_text().splitlines()));anatomy_path=ROOT/'data/derived/canonical/anatomy.json'
    anatomy_bytes=anatomy_path.read_bytes();anatomy=json.loads(anatomy_bytes);by_id={e['id']:e for e in anatomy['entities']}
    if hashlib.sha256(anatomy_bytes).hexdigest()!=parent_manifest['anatomy_sha256']:raise ValueError('Canonical anatomy differs from parent decomposition')
    (out/'inputs/anatomy.json').write_bytes(anatomy_bytes)
    for p in (Path(__file__),Path(__file__).with_name('decompose_muscle_dimensions.py'),Path(__file__).with_name('audit_surface_volume_readiness.py')):shutil.copyfile(p,out/'inputs'/p.name)
    shutil.copyfile(parent/'manifest.json',out/'inputs/parent-manifest.json');records=[];started=time.monotonic()
    for old in entity_records:
        if old['bulk_candidate']:continue
        entity=by_id[old['entity_id']];source=parent_sources[old['source_record_sha256']]
        for ref in (source['representations']['bulk'],source['representations']['strata']):
            if sha(parent/ref['path'])!=ref['sha256']:raise ValueError('Parent derived geometry changed')
        mapping_path=parent/source['mapping']
        if sha(mapping_path)!=parent_manifest['artifacts_sha256'][source['mapping']]:raise ValueError('Parent face mapping changed')
        mapping=json.loads(mapping_path.read_text());canonical_path=ROOT/entity['reference_geometry']['path']
        if sha(canonical_path)!=entity['reference_geometry']['sha256']:raise ValueError('Canonical geometry changed')
        data=json.loads(gzip.decompress(canonical_path.read_bytes()));x=np.asarray(data['positions']).reshape(-1,3);t=np.asarray(data['indices']).reshape(-1,3)
        raw=raw_obj_evidence(entity,data);parent_bulk=np.asarray(mapping['bulk_source_face_indices']);local=t[parent_bulk]
        selected,disks=paired_disks(x,local);moved=parent_bulk[selected]
        bulk=np.setdiff1d(parent_bulk,moved);strata=np.sort(np.r_[mapping['strata_source_face_indices'],moved]).astype(int)
        after=topology(x,t[bulk]);pairs=[]
        for disk in disks:
            disk['canonical_source_face_indices']=parent_bulk[disk['parent_bulk_face_indices']].tolist()
            disk['raw_obj_face_lines']=[raw['face_line_numbers'][i] for i in disk['canonical_source_face_indices']]
            for pair in disk['opposite_pairs']:pairs.append({'source_face_indices':parent_bulk[pair['source_face_indices']].tolist()})
        proof=cancellation(x,t,pairs)
        if proof['maximum_pair_winding_residual']>1e-10 or proof['maximum_scale_normalized_volume_residual']>1e-10:raise ValueError('Disk pair cancellation failed')
        zero=np.flatnonzero(np.all(np.cross(x[t[bulk,1]]-x[t[bulk,0]],x[t[bulk,2]]-x[t[bulk,0]])==0,axis=1))
        proposals=[]
        if len(zero):
            without=topology(x,t[np.delete(bulk,zero)])
            proposals.append({'kind':'exclude_degenerate_triangles','applied':False,'source_face_indices':bulk[zero].tolist(),
                              'raw_obj_face_lines':[raw['face_line_numbers'][i] for i in bulk[zero]],
                              'resulting_topology_blockers':without['topology_blockers'],'boundary_edges_if_excluded':without['welded_edges']['boundary_edges'],
                              'reason':'Exclusion changes boundary topology. Keep source faces; any local retriangulation or hole closure requires explicit geometry constraints and renewed validation.'})
        directory=out/'entities'/entity['id'];directory.mkdir()
        for name,ids in [('bulk',bulk),('strata',strata)]:
            derived=dict(data);derived['indices']=t[ids].reshape(-1).tolist();derived['display_faces']=len(ids)
            derived['dimensional_representation']={'role':name,'source_sha256':entity['reference_geometry']['sha256'],'physical_volume_or_mass_assigned':False}
            (directory/(name+'.json.gz')).write_bytes(gzip.compress(json.dumps(derived,separators=(',',':'),allow_nan=False).encode(),mtime=0))
        face_map={'bulk_source_face_indices':bulk.tolist(),'strata_source_face_indices':strata.tolist(),
                  'new_disk_source_face_indices':moved.tolist(),'new_paired_disks':disks,'faces_discarded':0}
        (directory/'mapping.json').write_text(json.dumps(face_map,separators=(',',':'))+'\n')
        raw.pop('face_line_numbers')
        record={'entity_id':entity['id'],'name':entity['name'],'source_geometry':entity['reference_geometry'],
                'raw_source_evidence':raw,'parent_remaining_bulk_topology':source['remaining_bulk_topology'],
                'remaining_bulk_topology':after,'bulk_candidate':after['topological_candidate'],'new_disk_components':len(disks),
                'new_surface_faces':len(moved),'new_pair_cancellation':proof,'unapplied_proposals':proposals,
                'geometry_cause':'Exact paired sheets already exist in authored raw OBJ; an exporter/authoring artifact is plausible but its cause and anatomical role are not established.',
                'synthesis_scope':'No missing anatomical volume inferred. For any unresolved gap, conditioned synthesis needs local topology, anatomical boundary and attachment evidence; no arbitrary filling applied.',
                'canonical_activation_applied':False,'mass_assigned_kg':None}
        (directory/'report.json').write_text(json.dumps(record,indent=2)+'\n');records.append(record)
    summary={'reviewed_muscles':len(records),'parent_bulk_candidates':418,'additional_bulk_candidates':sum(r['bulk_candidate'] for r in records),
             'whole_muscle_corpus_bulk_candidates':418+sum(r['bulk_candidate'] for r in records),
             'new_disk_components':sum(r['new_disk_components'] for r in records),'new_surface_faces':sum(r['new_surface_faces'] for r in records),
             'remaining_blocked':[r['entity_id'] for r in records if not r['bulk_candidate']],
             'faces_discarded':0,'coordinates_changed':False,'wall_seconds':time.monotonic()-started,
             'qualification':'Geometric only; no physical role/thickness, occupied tissue volume, self-intersection, mass ownership or calibration established.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    manifest={'schema':'ihm.remaining-muscle-bulk.v1','parent_manifest_path':str((parent/'manifest.json').relative_to(ROOT)),
              'parent_manifest_sha256':sha(parent/'manifest.json'),'anatomy_sha256':hashlib.sha256(anatomy_bytes).hexdigest(),
              'source_files_sha256':{r['source_geometry']['path']:r['source_geometry']['sha256'] for r in records},
              'artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},'original_files_modified':False}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(summary,output_dir=str(out)),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.self_test:self_test()
    if args.output:build(args.output)
    elif not args.self_test:parser.error('Supply --output or --self-test')
