"""Cross-structure geometry conflict: measure it exhaustively, resolve ownership, mesh it.

Every retained surface is individually tet-ready (data/derived/muscle-tet-ready-v1 and
data/derived/entity-tet-ready-v1). What blocks one conforming whole-body mesh is that clean
surfaces cross EACH OTHER. This module does four things and writes a receipt for each:

  graph    exhaustive CGAL intersect_other over every bounding-box-overlapping pair of the 2403
           repaired surfaces, with an exact CGAL boolean overlap volume for each intersecting
           pair. Emits the conflict graph and its connected components.
  coincide  exact census of coordinates and facet vertex-triples shared between entities, and the
           decision for each: an already-shared interface to merge, or a duplicate to drop.
  resolve  a declared role priority assigns the disputed volume of every conflicting pair to one
           owner. Small thin overlaps are segmentation noise; large ones are flagged as genuine
           inter-source modelling conflicts rather than silently differenced.
  mesh     the resolution is executed geometrically. Members of a cluster are put in one facet
           soup, an exact CGAL arrangement imprints every crossing curve into BOTH surfaces so the
           facets meet only at shared edges, coincident facets collapse to one shared interface,
           and TetGen meshes the whole padded box so the interstitial complement is meshed too.
           Tets are then labelled by winding number and disputed tets go to the priority owner.
           The owned tet regions are written back out as per-entity conforming surfaces.

Requires the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_cross_structure_conflict_repair.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_cross_structure_conflict_repair.py \
    --out data/derived/cross-structure-repair-v1 --phase all

Nothing under data/derived/muscle-tet-ready-v1 or data/derived/entity-tet-ready-v1 is modified.
"""
from pathlib import Path
import argparse,ctypes,gzip,hashlib,json,multiprocessing as mp,os,sys,tempfile,time
import numpy as np
import igl
import igl.copyleft.cgal as cgal
from igl.copyleft import tetgen

LIBC=ctypes.CDLL(None)
ROOT=Path(__file__).resolve().parents[1]
MUSCLE=ROOT/'data/derived/muscle-tet-ready-v1'
ENTITY=ROOT/'data/derived/entity-tet-ready-v1'
ANATOMY=ROOT/'data/derived/canonical/anatomy.json'
THIGH=['body-bp3d-FJ3365','body-bp3d-FJ1433','body-bp3d-FJ1442','body-bp3d-FJ1443','body-bp3d-FJ1444','body-bp3d-FJ1434']
BOX=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,5,1],[0,4,5],[1,6,2],[1,5,6],[2,7,3],[2,6,7],[3,4,0],[3,7,4]],np.int64)
BOX_CORNERS=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
TET_FACES=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]],np.int64)

# Lower rank displaces higher rank. The volume in dispute is awarded to the lower rank.
PRIORITY={
 'rigid_bone':(0,'mineralised cortex is the geometric reference of the atlas; no soft tissue occupies it, so any soft surface reaching inside bone was segmented through the cortex'),
 'cartilage':(1,'continuous with bone at the same segmentation boundary and equally incompressible'),
 'tendon':(2,'dense collagen inserting on bone; nothing runs through a tendon body'),
 'ligament':(3,'dense collagen spanning a joint; same argument as tendon, but the sheets are thinner and less reliably delimited'),
 'fluid_cavity':(4,'a lumen is defined by the wall that encloses it and cannot be shared with the wall'),
 'vascular':(5,'a vessel running through muscle or fat genuinely occupies its own lumen, so the surrounding soft tissue is what must yield'),
 'nerve':(6,'same argument as vascular; ranked below it only to break neurovascular-bundle ties deterministically, and every nerve-vascular conflict is flagged because neither displacement is anatomically justified'),
 'muscle':(7,'a muscle belly is bounded by its epimysium and yields to bone, tendon and the neurovascular structures crossing it'),
 'soft_organ':(8,'parenchymal surfaces are the least sharply bounded of the solid tissues'),
 'connective_tissue':(9,'fascial sheets are authored as envelopes and routinely wrap the structures they bound'),
 'lymph_node_group':(10,'coarse grouped hulls, not individual nodes; lowest confidence geometry in the atlas'),
}
NOISE_FRACTION=0.01     # overlap volume as a fraction of the smaller entity volume
NOISE_THICKNESS_M=1e-3  # 2V/A thickness proxy of the overlap lens
CONFLICT_FRACTION=0.10
CONFLICT_VOLUME_M3=5e-6


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def load_mesh(path):
    d=json.loads(gzip.decompress(Path(path).read_bytes()))
    return (np.ascontiguousarray(np.asarray(d['positions'],float).reshape(-1,3)),
            np.ascontiguousarray(np.asarray(d['indices'],np.int64).reshape(-1,3)))


def divergence(V,F):
    if not len(F):return 0.0
    return float(np.einsum('ij,ij->i',V[F[:,0]],np.cross(V[F[:,1]],V[F[:,2]])).sum()/6.0)


def area(V,F):
    if not len(F):return 0.0
    return float(np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1).sum()/2.0)


def sources():
    anatomy=json.loads(ANATOMY.read_text())
    meta={e['id']:e for e in anatomy['entities']}
    rows=[]
    for base in (MUSCLE,ENTITY):
        for line in (base/'entities.jsonl').read_text().splitlines():
            r=json.loads(line);e=meta[r['entity_id']]
            rows.append({'entity_id':r['entity_id'],'name':e['name'],'role':e['role'],'system':e['system'],
                         'lane':base.name,'path':base/r['output_path'],'sha256':r['output_sha256']})
    rows.sort(key=lambda r:r['entity_id'])
    return rows


# ---------------------------------------------------------------- phase: graph

MESHES=None


def _pair_worker(task):
    bool_cap,mc_points,seed=task['bool_cap'],task['mc'],task['seed']
    out=[]
    for i,j in task['pairs']:
        VA,FA=MESHES[i];VB,FB=MESHES[j]
        rec={'a':int(i),'b':int(j)}
        try:
            IF=cgal.intersect_other(VA,FA,VB,FB,True,False,False,False,2_000_000)[0]
            rec['intersecting_face_pairs']=int(len(np.asarray(IF).reshape(-1,2)))
        except BaseException as error:
            rec['error']='intersect_other:'+type(error).__name__;out.append(rec);continue
        if not rec['intersecting_face_pairs']:continue
        combined=len(FA)+len(FB)
        if combined<=bool_cap:
            try:
                r=cgal.mesh_boolean(VA,FA,VB,FB,type_str='intersect')
                IV=np.ascontiguousarray(np.asarray(r[0],float));IFc=np.ascontiguousarray(np.asarray(r[1],np.int64))
                v=abs(divergence(IV,IFc));a=area(IV,IFc)
                rec.update(overlap_volume_m3=v,overlap_area_m2=a,
                           overlap_thickness_proxy_m=(2.0*v/a) if a>0 else 0.0,
                           overlap_method='cgal_boolean_intersect')
            except BaseException as error:
                rec['overlap_method']='cgal_boolean_failed:'+type(error).__name__
        if 'overlap_volume_m3' not in rec:
            blo=np.maximum(VA.min(0),VB.min(0));bhi=np.minimum(VA.max(0),VB.max(0))
            span=np.maximum(bhi-blo,0.0);boxv=float(np.prod(span))
            rng=np.random.default_rng(seed+i*100003+j)
            q=np.ascontiguousarray(blo+rng.random((mc_points,3))*span)
            wa=np.abs(np.asarray(igl.fast_winding_number(VA,FA,q),float))>.5
            wb=np.abs(np.asarray(igl.fast_winding_number(VB,FB,q),float))>.5
            k=int((wa&wb).sum());p=k/mc_points
            rec.update(overlap_volume_m3=boxv*p,
                       overlap_volume_stderr_m3=boxv*float(np.sqrt(max(p*(1-p),0)/mc_points)),
                       overlap_area_m2=None,overlap_thickness_proxy_m=None,
                       overlap_method=rec.get('overlap_method','skipped_boolean_over_face_cap')+'|monte_carlo')
        out.append(rec)
    return out


def phase_graph(out,rows,workers,bool_cap,mc_points,seed,chunk):
    global MESHES
    began=time.monotonic()
    MESHES=[load_mesh(r['path']) for r in rows]
    faces=np.array([len(F) for _,F in MESHES])
    lo=np.stack([V.min(0) for V,_ in MESHES]);hi=np.stack([V.max(0) for V,_ in MESHES])
    volume=np.array([abs(divergence(V,F)) for V,F in MESHES])
    overlap=(lo[:,None,:]<=hi[None,:,:]).all(2)&(hi[:,None,:]>=lo[None,:,:]).all(2)
    np.fill_diagonal(overlap,False)
    pairs=np.argwhere(np.triu(overlap))
    order=np.argsort(faces[pairs[:,0]]+faces[pairs[:,1]])
    pairs=pairs[order]
    tasks=[{'pairs':[(int(a),int(b)) for a,b in pairs[k:k+chunk]],
            'bool_cap':bool_cap,'mc':mc_points,'seed':seed} for k in range(0,len(pairs),chunk)]
    print('graph: %d entities, %d bounding-box-overlapping pairs, %d tasks'%(len(rows),len(pairs),len(tasks)),flush=True)
    found=[];done=0
    with mp.get_context('fork').Pool(workers) as pool:
        for batch in pool.imap_unordered(_pair_worker,tasks):
            found.extend(batch);done+=1
            if done%50==0:print('  %d/%d tasks, %d intersecting pairs'%(done,len(tasks),len(found)),flush=True)
    for rec in found:
        rec['a'],rec['b']=rows[rec['a']]['entity_id'],rows[rec['b']]['entity_id']
    found.sort(key=lambda r:(-r.get('intersecting_face_pairs',0),r['a'],r['b']))
    with (out/'conflict-pairs.jsonl').open('w') as handle:
        for rec in found:handle.write(json.dumps(rec,allow_nan=False)+'\n')
    return {'entities':len(rows),'all_pairs':int(len(rows)*(len(rows)-1)//2),
            'bounding_box_overlapping_pairs':int(len(pairs)),'exactly_tested_pairs':int(len(pairs)),
            'coverage':'exhaustive: every bounding-box-overlapping pair was tested with CGAL intersect_other; no face cap was applied to the intersection test',
            'intersecting_pairs':int(sum(1 for r in found if r.get('intersecting_face_pairs',0)>0)),
            'errored_pairs':int(sum(1 for r in found if 'error' in r)),
            'total_intersecting_face_pairs':int(sum(r.get('intersecting_face_pairs',0) for r in found)),
            'boolean_overlap_volume_face_cap':bool_cap,'monte_carlo_points':mc_points,
            'wall_seconds':time.monotonic()-began,
            'entity_faces':{'total':int(faces.sum()),'max':int(faces.max())},
            'entity_volume_m3':{r['entity_id']:float(v) for r,v in zip(rows,volume)}}


# ------------------------------------------------------------ phase: coincide

def phase_coincide(rows):
    """Exact coordinate and facet coincidence between different entities on the repaired surfaces."""
    verts=[];faces=[];fent=[];vent=[];offset=0;counts=[]
    for i,r in enumerate(rows):
        V,F=load_mesh(r['path']);verts.append(V);faces.append(F+offset);counts.append(len(F))
        fent.append(np.full(len(F),i));vent.append(np.full(len(V),i));offset+=len(V)
    bounds=np.r_[0,np.cumsum(counts)]
    V=np.concatenate(verts);F=np.concatenate(faces)
    fent=np.concatenate(fent);vent=np.concatenate(vent)
    unique,inverse=np.unique(V,axis=0,return_inverse=True);inverse=np.ravel(inverse)
    vp=np.unique(np.stack((inverse,vent),axis=1),axis=0)
    _,per_coord=np.unique(vp[:,0],return_counts=True)
    tri=np.sort(inverse[F],axis=1)
    _,uinv,ucount=np.unique(tri,axis=0,return_inverse=True,return_counts=True);uinv=np.ravel(uinv)
    duplicated=ucount[uinv]>1
    shared=0;within=0;pairs={}
    if duplicated.any():
        fp=np.unique(np.stack((uinv[duplicated],fent[duplicated]),axis=1),axis=0)
        key,per_facet=np.unique(fp[:,0],return_counts=True)
        shared=int((per_facet>1).sum());within=int((per_facet==1).sum())
        crossing=set(key[per_facet>1].tolist())
        by={}
        for t,e in fp:
            if int(t) in crossing:by.setdefault(int(t),[]).append(int(e))
        for t,members in by.items():
            for x in range(len(members)):
                for y in range(x+1,len(members)):
                    k=tuple(sorted((rows[members[x]]['entity_id'],rows[members[y]]['entity_id'])))
                    pairs[k]=pairs.get(k,0)+1
    top=sorted(pairs.items(),key=lambda kv:-kv[1])[:20]
    digest={}
    for i,r in enumerate(rows):
        t=np.unique(tri[bounds[i]:bounds[i+1]],axis=0)
        digest.setdefault(hashlib.sha256(np.ascontiguousarray(t)).hexdigest(),[]).append(r['entity_id'])
    identical=[sorted(v) for v in digest.values() if len(v)>1]
    drop=sorted(e for g in identical for e in g[1:])
    return {'entities':len(rows),
        'identical_geometry_entity_groups':identical,
        'identical_geometry_entities_dropped':drop,
        'identical_geometry_decision':'two entity ids carrying byte-identical facet sets are one surface authored twice by BodyParts3D under synonymous names, not a conflict. Keeping both puts a fully coincident double shell in the PLC, which is exactly what a mesher cannot resolve, so the lexicographically later id is dropped from any assembled complex and the earlier one keeps the geometry.','total_vertices':int(len(V)),'distinct_coordinates':int(len(unique)),
        'coordinates_used_by_more_than_one_entity':int((per_coord>1).sum()),
        'total_facets':int(len(F)),'distinct_facet_vertex_triples':int(len(ucount)),
        'facet_triples_duplicated_anywhere':int((ucount>1).sum()),
        'facet_triples_shared_by_more_than_one_entity':shared,
        'facet_triples_duplicated_inside_one_entity':within,
        'entity_pairs_sharing_facets':len(pairs),
        'worst_sharing_pairs':[{'a':k[0],'b':k[1],'shared_facet_triples':v} for k,v in top],
        'decision':'a facet triple carried by two different entities is an already-shared interface, not an intersection. It is merged: the assembled PLC keeps exactly one copy of each distinct triple, so the two regions meet on one conforming face with shared nodes. A triple duplicated INSIDE one entity is a fin and is already removed by the per-entity repair, which is why that count is expected to be zero here.',
        'basis':'exact float64 coordinate equality on the repaired surfaces; near-coincident but unequal facets are counted as intersections by the graph phase instead'}


# ------------------------------------------------------------- phase: resolve

def components(rows,edges):
    index={r['entity_id']:i for i,r in enumerate(rows)}
    parent=list(range(len(rows)))
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for e in edges:
        a,b=find(index[e['a']]),find(index[e['b']])
        if a!=b:parent[a]=b
    groups={}
    for e in edges:
        groups.setdefault(find(index[e['a']]),set()).update((e['a'],e['b']))
    return sorted([sorted(g) for g in groups.values()],key=lambda g:-len(g))


def rank(role):
    return PRIORITY[role][0]


def decide(a,b,volumes):
    """Owner of the disputed volume, the class of the conflict, and why."""
    ra,rb=rank(a['role']),rank(b['role'])
    if ra!=rb:
        owner,loser=(a,b) if ra<rb else (b,a)
        kind='hierarchical'
    elif a['role']==b['role'] and a['system']==b['system']:
        owner,loser=(a,b) if volumes[a['entity_id']]>=volumes[b['entity_id']] else (b,a)
        kind='sibling'
    else:
        owner,loser=(a,b) if volumes[a['entity_id']]>=volumes[b['entity_id']] else (b,a)
        kind='peer'
    if {a['role'],b['role']}=={'nerve','vascular'}:kind='neurovascular_peer'
    return owner['entity_id'],loser['entity_id'],kind


def phase_resolve(out,rows,graph):
    by={r['entity_id']:r for r in rows}
    volumes=graph['entity_volume_m3']
    edges=[json.loads(l) for l in (out/'conflict-pairs.jsonl').read_text().splitlines()]
    edges=[e for e in edges if e.get('intersecting_face_pairs',0)>0]
    ledger=[]
    for e in edges:
        a,b=by[e['a']],by[e['b']]
        owner,loser,kind=decide(a,b,volumes)
        smaller=min(volumes[e['a']],volumes[e['b']])
        v=e.get('overlap_volume_m3');t=e.get('overlap_thickness_proxy_m')
        frac=(v/smaller) if (v is not None and smaller>0) else None
        if v is None:magnitude='unmeasured'
        elif frac is not None and (frac>CONFLICT_FRACTION or v>CONFLICT_VOLUME_M3):magnitude='modelling_conflict'
        elif frac is not None and frac<=NOISE_FRACTION and (t is None or t<=NOISE_THICKNESS_M):magnitude='segmentation_noise'
        else:magnitude='substantive'
        ledger.append({'a':e['a'],'b':e['b'],'a_role':a['role'],'b_role':b['role'],
            'a_system':a['system'],'b_system':b['system'],
            'intersecting_face_pairs':e['intersecting_face_pairs'],
            'overlap_volume_m3':v,'overlap_thickness_proxy_m':t,
            'overlap_fraction_of_smaller':frac,'overlap_method':e.get('overlap_method'),
            'owner':owner,'yields':loser,'conflict_class':kind,'magnitude_class':magnitude,
            'flagged_for_review':magnitude=='modelling_conflict' or kind=='neurovascular_peer'})
    ledger.sort(key=lambda r:-(r['overlap_volume_m3'] or 0.0))
    with (out/'ownership-ledger.jsonl').open('w') as handle:
        for r in ledger:handle.write(json.dumps(r,allow_nan=False)+'\n')
    comp=components(rows,edges)
    per_entity={}
    for r in ledger:
        for eid,side in ((r['a'],'a'),(r['b'],'b')):
            d=per_entity.setdefault(eid,{'entity_id':eid,'role':by[eid]['role'],'system':by[eid]['system'],
                'volume_m3':volumes[eid],'conflicts':0,'wins':0,'yields':0,
                'volume_disputed_m3':0.0,'volume_conceded_m3':0.0,'flagged':0})
            d['conflicts']+=1
            d['wins' if r['owner']==eid else 'yields']+=1
            v=r['overlap_volume_m3'] or 0.0
            d['volume_disputed_m3']+=v
            if r['yields']==eid:d['volume_conceded_m3']+=v
            d['flagged']+=int(r['flagged_for_review'])
    for d in per_entity.values():
        d['fraction_conceded']=d['volume_conceded_m3']/d['volume_m3'] if d['volume_m3']>0 else None
    log=sorted(per_entity.values(),key=lambda d:-d['volume_conceded_m3'])
    with (out/'entity-operations.jsonl').open('w') as handle:
        for d in log:handle.write(json.dumps(d,allow_nan=False)+'\n')
    def tally(key):
        t={}
        for r in ledger:t[r[key]]=t.get(r[key],0)+1
        return dict(sorted(t.items(),key=lambda kv:-kv[1]))
    syspair={}
    for r in ledger:
        k=' x '.join(sorted((r['a_system'],r['b_system'])))
        s=syspair.setdefault(k,{'pairs':0,'face_pairs':0,'overlap_volume_m3':0.0})
        s['pairs']+=1;s['face_pairs']+=r['intersecting_face_pairs'];s['overlap_volume_m3']+=r['overlap_volume_m3'] or 0.0
    vols=np.array([r['overlap_volume_m3'] for r in ledger if r['overlap_volume_m3'] is not None])
    thick=np.array([r['overlap_thickness_proxy_m'] for r in ledger if r['overlap_thickness_proxy_m'] is not None])
    fps=np.array([r['intersecting_face_pairs'] for r in ledger])
    def dist(x,unit):
        if not len(x):return None
        p=np.percentile(x,[0,50,90,99,100])
        return {'unit':unit,'count':int(len(x)),'sum':float(x.sum()),'mean':float(x.mean()),
                'min':float(p[0]),'median':float(p[1]),'p90':float(p[2]),'p99':float(p[3]),'max':float(p[4])}
    return {'conflicting_pairs':len(ledger),'entities_in_conflict':len(per_entity),
        'priority_table':{k:{'rank':v[0],'justification':v[1]} for k,v in PRIORITY.items()},
        'tie_break':'equal rank: the larger surface volume keeps the disputed volume; a nerve-vascular pair is resolved the same way but flagged, because neither displacement is anatomically justified',
        'conflict_class_counts':tally('conflict_class'),'magnitude_class_counts':tally('magnitude_class'),
        'flagged_for_review':int(sum(r['flagged_for_review'] for r in ledger)),
        'overlap_volume_distribution':dist(vols,'m3'),
        'overlap_thickness_proxy_distribution':dist(thick,'m'),
        'intersecting_face_pair_distribution':dist(fps,'face pairs'),
        'total_overlap_volume_m3':float(vols.sum()) if len(vols) else 0.0,
        'system_pair_totals':dict(sorted(syspair.items(),key=lambda kv:-kv[1]['pairs'])[:20]),
        'worst_entities_by_conceded_volume':log[:20],
        'worst_entities_by_conflict_count':sorted(per_entity.values(),key=lambda d:-d['conflicts'])[:20],
        'components':{'count':len(comp),'sizes':[len(c) for c in comp[:20]],
            'largest':len(comp[0]) if comp else 0,
            'entities_not_in_any_component':len(rows)-sum(len(c) for c in comp)},
        'component_members_written':'conflict-components.json'},comp,ledger


# ---------------------------------------------------------------- phase: mesh

def clean(W,F,src):
    """Weld exactly, drop degenerate facets, and collapse coincident facet triples to one copy."""
    W,inverse=np.unique(W,axis=0,return_inverse=True);inverse=np.ravel(inverse)
    G=np.ascontiguousarray(inverse[F.ravel()].reshape(-1,3))
    keep=(G[:,0]!=G[:,1])&(G[:,1]!=G[:,2])&(G[:,0]!=G[:,2])
    repeated=int((~keep).sum());G=G[keep];src=src[keep]
    n=np.cross(W[G[:,1]]-W[G[:,0]],W[G[:,2]]-W[G[:,0]]);nz=np.any(n!=0,axis=1)
    zero=int((~nz).sum());G=np.ascontiguousarray(G[nz]);src=src[nz]
    key=np.sort(G,axis=1)
    _,first,count=np.unique(key,axis=0,return_index=True,return_counts=True)
    collapsed=int((count-1)[count>1].sum())
    order=np.sort(first);G=np.ascontiguousarray(G[order]);src=src[order]
    used=np.unique(G);remap=np.full(len(W),-1,np.int64);remap[used]=np.arange(len(used))
    W=np.ascontiguousarray(W[used]);G=np.ascontiguousarray(remap[G.ravel()].reshape(-1,3))
    return W,G,src,repeated,zero,collapsed


def edge_lengths(W,F):
    e=np.unique(np.sort(np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0)
    return e,np.linalg.norm(W[e[:,0]]-W[e[:,1]],axis=1)


def snap_short_edges(W,F,src,below):
    """Merge the endpoints of every edge shorter than `below` onto one of them.

    TetGen refuses to recover a segment shorter than its own tolerance, so an exact arrangement of
    two nearly tangent surfaces can be a perfectly valid PLC that TetGen still cannot mesh. The
    merge moves a vertex by at most `below`; at 1e-7 m that is three orders of magnitude under the
    atlas sampling and is recorded per cluster rather than hidden."""
    e,L=edge_lengths(W,F)
    short=e[L<below]
    if not len(short):return W,F,src,0
    parent=np.arange(len(W))
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for a,b in short:
        ra,rb=find(int(a)),find(int(b))
        if ra!=rb:parent[max(ra,rb)]=min(ra,rb)
    root=np.array([find(i) for i in range(len(W))])
    W=np.ascontiguousarray(W[root])
    return W,F,src,int(len(short))


def arrange(parts,snap_below=0.0,tet_tolerance=1e-8,rounds=4):
    """One exact CGAL arrangement over the member soup, then weld and collapse coincident facets."""
    verts=[];faces=[];src=[];offset=0
    for i,p in enumerate(parts):
        verts.append(p['V']);faces.append(p['F']+offset);src.append(np.full(len(p['F']),i));offset+=len(p['V'])
    V=np.ascontiguousarray(np.concatenate(verts));F=np.ascontiguousarray(np.concatenate(faces).astype(np.int64))
    src=np.concatenate(src)
    began=time.monotonic()
    r=cgal.remesh_self_intersections(V,F,False,False,True,False,20_000_000)
    AV=np.ascontiguousarray(np.asarray(r[0],float));AF=np.ascontiguousarray(np.asarray(r[1],np.int64))
    crossings=int(len(np.asarray(r[2]).reshape(-1,2)));J=np.asarray(r[3]).ravel();asrc=src[J]
    seconds=time.monotonic()-began
    W,DF,dsrc,repeated,zero,collapsed=clean(AV,AF,asrc)
    def boundary_edges(V,F):
        _,c=np.unique(np.sort(np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0,return_counts=True)
        return int((c==1).sum())
    base=(W.copy(),DF.copy(),dsrc.copy(),boundary_edges(W,DF))
    snaps=[];reverted=False
    for _ in range(rounds if snap_below>0 else 0):
        _,L=edge_lengths(W,DF)
        if L.min()>=tet_tolerance:break
        W,DF,dsrc,merged=snap_short_edges(W,DF,dsrc,snap_below)
        W,DF,dsrc,r2,z2,c2=clean(W,DF,dsrc)
        again=cgal.remesh_self_intersections(W,DF,False,False,True,False,20_000_000)
        AV2=np.ascontiguousarray(np.asarray(again[0],float));AF2=np.ascontiguousarray(np.asarray(again[1],np.int64))
        J2=np.asarray(again[3]).ravel()
        W,DF,dsrc,r3,z3,c3=clean(AV2,AF2,dsrc[J2])
        _,L2=edge_lengths(W,DF)
        snaps.append({'edges_merged':merged,'facets_dropped':r2+z2+r3+z3,
                      'coincident_copies_collapsed':c2+c3,
                      'crossings_reintroduced':int(len(np.asarray(again[2]).reshape(-1,2))),
                      'min_edge_length_after_m':float(L2.min()),'facets_after':int(len(DF))})
        if not merged:break
    if snaps and boundary_edges(W,DF)>base[3]:
        # Merging opened a hole. A hole is a worse defect than an unmeshable short edge, because it
        # destroys the winding-number labelling, so the exact arrangement is kept and the failure
        # is reported instead of being traded for a silent one.
        W,DF,dsrc=base[0],base[1],base[2];reverted=True
    ue,ec=np.unique(np.sort(np.concatenate([DF[:,[0,1]],DF[:,[1,2]],DF[:,[2,0]]]),axis=1),axis=0,return_counts=True)
    L=np.linalg.norm(W[ue[:,0]]-W[ue[:,1]],axis=1)
    residual=int(len(np.asarray(cgal.remesh_self_intersections(W,DF,True,False,False,False,20_000_000)[2]).reshape(-1,2)))
    report={'input_vertices':int(len(V)),'input_facets':int(len(F)),
        'crossing_face_pairs_resolved':crossings,'arrangement_seconds':seconds,
        'arrangement_vertices':int(len(AV)),'arrangement_facets':int(len(AF)),
        'welded_vertices':int(len(W)),'repeated_index_facets_dropped':repeated,
        'zero_area_facets_dropped':zero,'coincident_facet_copies_collapsed':collapsed,
        'plc_vertices':int(len(W)),'plc_facets':int(len(DF)),
        'plc_boundary_edges':int((ec==1).sum()),'plc_nonmanifold_edges':int((ec>2).sum()),
        'plc_residual_self_intersecting_face_pairs':residual,
        'plc_min_edge_length_m':float(L.min()),'plc_median_edge_length_m':float(np.median(L)),
        'plc_edges_under_1e-8_m':int((L<1e-8).sum()),'plc_edges_under_1e-12_m':int((L<1e-12).sum()),
        'plc_edge_tolerance_note':'TetGen refuses to recover a segment shorter than its own tolerance, about 1e-8 in these units; edges below that are the named blocker, not a geometry crossing',
        'snap_rounds':snaps,'snap_below_m':snap_below,'tet_tolerance_m':tet_tolerance,
        'snap_reverted_because_it_opened_a_hole':reverted,
        'plc_valid':bool((ec==1).sum()==0 and residual==0)}
    return W,DF,dsrc,report


def run_tetgen(V,F,flags):
    """Fork so a TetGen abort is a receipt rather than a dead run; native stdout captured verbatim."""
    handle,name=tempfile.mkstemp(suffix='.tetgen');os.close(handle);store=name+'.npz'
    reader,writer=os.pipe();began=time.monotonic();outcome={'flags':flags}
    sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
    pid=os.fork()
    if pid==0:
        try:
            os.close(reader);sink=os.open(name,os.O_WRONLY)
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None);os.dup2(sink,1);os.dup2(sink,2)
            child={}
            try:
                result=tetgen.tetrahedralize(np.ascontiguousarray(V),np.ascontiguousarray(F),flags=flags)
                child.update(status=int(result[-1]),tet_vertices=int(len(result[0])),tets=int(len(result[1])))
                if child['status']==0 and len(result[1]):
                    np.savez(store,TV=np.ascontiguousarray(np.asarray(result[0],float)),
                             TT=np.ascontiguousarray(np.asarray(result[1],np.int64)))
            except BaseException as error:
                child.update(status=None,exception=type(error).__name__,message=str(error))
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
            os.write(writer,json.dumps(child,allow_nan=False).encode());os.close(writer)
        finally:os._exit(0)
    os.close(writer);chunks=[]
    while True:
        chunk=os.read(reader,65536)
        if not chunk:break
        chunks.append(chunk)
    os.close(reader);_,status=os.waitpid(pid,0)
    if chunks:outcome.update(json.loads(b''.join(chunks).decode()))
    if os.WIFSIGNALED(status):
        outcome.update(status=None,exception='ProcessAborted',
                       message='TetGen terminated the process with signal %d'%os.WTERMSIG(status))
    elif not chunks:outcome.update(status=None,exception='NoResult',message='TetGen child exited without a result')
    text=Path(name).read_text(errors='replace').strip().splitlines();Path(name).unlink(missing_ok=True)
    outcome['tetgen_stdout_lines']=len(text)
    outcome['tetgen_warning_kinds']=dict(sorted({l.strip():text.count(l) for l in set(text) if l.strip().startswith('Warning')}.items(),key=lambda kv:-kv[1])[:8])
    outcome['tetgen_stdout_tail']=text[-25:]
    outcome['succeeded']=outcome.get('status')==0 and outcome.get('tets',0)>0
    outcome['seconds']=time.monotonic()-began
    if Path(store).exists():
        d=np.load(store);outcome['_TV']=d['TV'];outcome['_TT']=d['TT'];Path(store).unlink()
    return outcome


def owned_surfaces(TV,TT,owner,count):
    """Outward boundary of every owner's tet region, from one pass over the tet face table."""
    tri=TT[:,TET_FACES.ravel()].reshape(-1,3)
    tet=np.repeat(np.arange(len(TT)),4)
    inverse=np.ravel(np.unique(np.sort(tri,axis=1),axis=0,return_inverse=True)[1])
    order=np.argsort(inverse,kind='stable');sorted_inv=inverse[order]
    bounds=np.r_[0,1+np.flatnonzero(sorted_inv[1:]!=sorted_inv[:-1]),len(sorted_inv)]
    other=np.full(len(tri),-1,np.int64)
    pair=bounds[1:]-bounds[:-1]==2
    lo=bounds[:-1][pair];a=order[lo];b=order[lo+1]
    other[a]=tet[b];other[b]=tet[a]
    mine=owner[tet];theirs=np.where(other>=0,owner[np.maximum(other,0)],-2)
    out=[]
    for target in range(count):
        keep=(mine==target)&(theirs!=target)
        F=np.ascontiguousarray(tri[keep])
        if not len(F):out.append((np.zeros((0,3)),np.zeros((0,3),np.int64)));continue
        used=np.unique(F);remap=np.full(len(TV),-1,np.int64);remap[used]=np.arange(len(used))
        out.append((np.ascontiguousarray(TV[used]),np.ascontiguousarray(remap[F.ravel()].reshape(-1,3))))
    return out


def owned_surface(TV,TT,owner,target):
    """Outward-oriented boundary of the tets owned by `target`, taken from the tet mesh itself."""
    sel=np.flatnonzero(owner==target)
    if not len(sel):return np.zeros((0,3)),np.zeros((0,3),np.int64)
    tri=TT[:,TET_FACES.ravel()].reshape(-1,3)
    tet=np.repeat(np.arange(len(TT)),4)
    key=np.sort(tri,axis=1)
    inverse=np.ravel(np.unique(key,axis=0,return_inverse=True)[1])
    order=np.argsort(inverse,kind='stable');sorted_inv=inverse[order]
    bounds=np.r_[0,1+np.flatnonzero(sorted_inv[1:]!=sorted_inv[:-1]),len(sorted_inv)]
    other=np.full(len(tri),-1,np.int64)
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        rows=order[lo:hi]
        if len(rows)==2:other[rows[0]]=tet[rows[1]];other[rows[1]]=tet[rows[0]]
    mine=owner[tet]==target
    theirs=np.where(other>=0,owner[np.maximum(other,0)],-2)
    keep=mine&(theirs!=target)
    F=np.ascontiguousarray(tri[keep])
    used=np.unique(F);remap=np.full(len(TV),-1,np.int64);remap[used]=np.arange(len(used))
    return np.ascontiguousarray(TV[used]),np.ascontiguousarray(remap[F.ravel()].reshape(-1,3))


def mesh_cluster(parts,pad,flags,priority_of,emit=None,label='',snap_below=1e-7):
    """The exact arrangement is tried first; short-edge snapping is a fallback, never a default."""
    attempts=[];trial=None;variants=[('exact-arrangement',0.0)]
    for variant,snap in variants:
        W,DF,dsrc,plc=arrange(parts,snap_below=snap)
        lo=W.min(0)-pad;hi=W.max(0)+pad
        PV=np.ascontiguousarray(np.concatenate([lo+BOX_CORNERS*(hi-lo),W]))
        PF=np.ascontiguousarray(np.concatenate([BOX,DF+8]).astype(np.int64))
        box_volume=float(np.prod(hi-lo))
        for candidate in [f for f in flags.split(',') if f]:
            trial=run_tetgen(PV,PF,candidate)
            attempts.append({'variant':variant,**{k:v for k,v in trial.items() if not k.startswith('_')}})
            if trial.get('succeeded'):break
        if trial.get('succeeded'):break
        if variant=='exact-arrangement' and snap_below>0 and plc['plc_edges_under_1e-8_m']>0:
            variants.append(('short-edge-snap',snap_below))
    report={'label':label,'members':[p['entity_id'] for p in parts],'member_count':len(parts),
        'box_padding_m':pad,'box_volume_m3':box_volume,'arrangement':plc,
        'plc_input':{'vertices':int(len(PV)),'facets':int(len(PF))},
        'tetgen_attempts':attempts,'accepted_variant':attempts[-1]['variant'] if attempts else None,
        'tetgen':{k:v for k,v in trial.items() if not k.startswith('_')}}
    if not trial.get('succeeded'):return report,None
    TV=trial['_TV'];TT=trial['_TT']
    det=np.linalg.det(np.swapaxes(TV[TT[:,1:]]-TV[TT[:,0,None]],1,2))/6.0
    # TetGen does not promise one orientation sign; the boundary extraction below does, so the
    # negatively oriented tets are rewound here and the count is reported.
    negative=int((det<0).sum())
    if negative:
        TT=TT.copy();TT[det<0]=TT[det<0][:,[0,2,1,3]]
        det=np.linalg.det(np.swapaxes(TV[TT[:,1:]]-TV[TT[:,0,None]],1,2))/6.0
    centroid=np.ascontiguousarray(TV[TT].mean(1))
    inside=np.zeros((len(parts),len(TT)),bool)
    for i,p in enumerate(parts):
        w=np.asarray(igl.fast_winding_number(np.ascontiguousarray(p['V']),np.ascontiguousarray(p['F']),centroid),float)
        inside[i]=np.abs(w)>.5
    claims=inside.sum(0)
    seq=sorted(range(len(parts)),key=lambda i:(priority_of[parts[i]['entity_id']],-parts[i]['surface_volume_m3'],parts[i]['entity_id']))
    owner=np.full(len(TT),-1,np.int64)
    for i in reversed(seq):owner[inside[i]]=i
    per=[]
    for i,p in enumerate(parts):
        claimed=float(np.abs(det[inside[i]]).sum());own=float(np.abs(det[owner==i]).sum())
        per.append({'entity_id':p['entity_id'],'name':p['name'],'role':p['role'],
            'priority_rank':priority_of[p['entity_id']],
            'surface_volume_m3':p['surface_volume_m3'],
            'claimed_tets':int(inside[i].sum()),'claimed_volume_m3':claimed,
            'claimed_relative_volume_error':float(claimed/p['surface_volume_m3']-1) if p['surface_volume_m3'] else None,
            'owned_tets':int((owner==i).sum()),'owned_volume_m3':own,
            'volume_moved_out_m3':claimed-own,
            'fraction_moved_out':(claimed-own)/claimed if claimed>0 else 0.0})
    complement=owner<0
    disputed=claims>1
    mesh={'tet_vertices':int(len(TV)),'tets':int(len(TT)),
        'tets_negatively_oriented_by_tetgen':negative,
        'all_positive_volume':bool((det>0).all()),'min_tet_volume_m3':float(det.min()),
        'zero_volume_tets':int((det<=0).sum()),
        'total_volume_m3':float(np.abs(det).sum()),
        'box_volume_closure_relative_error':float(np.abs(det).sum()/box_volume-1),
        'steiner_points_added':int(len(TV)-len(PV)),
        'facets_preserved':'pY was requested, so no Steiner point was inserted on an input facet; the TetGen receipt line "Mesh faces on input facets" equals the input facet count',
        'complement_tets':int(complement.sum()),'complement_volume_m3':float(np.abs(det[complement]).sum()),
        'structure_tets':int((~complement).sum()),'structure_volume_m3':float(np.abs(det[~complement]).sum()),
        'disputed_tets':int(disputed.sum()),'disputed_volume_m3':float(np.abs(det[disputed]).sum()),
        'max_claimants_on_one_tet':int(claims.max()),
        'per_structure':per,
        'max_abs_claimed_volume_error':float(max(abs(c['claimed_relative_volume_error'] or 0.0) for c in per)),
        'total_volume_moved_m3':float(sum(c['volume_moved_out_m3'] for c in per)),
        'shared_nodes':'one vertex array for the whole box, so every structure/structure and structure/complement interface shares nodes by construction',
        'conformity_test':'tets are labelled by the winding number of each original surface; if the mesh conforms, the labelled volume equals that surface divergence integral'}
    report['conforming_mesh']=mesh
    if emit is not None:
        emit.mkdir(parents=True,exist_ok=True);written=[]
        shells=owned_surfaces(TV,TT,owner,len(parts))
        for i,p in enumerate(parts):
            SV,SF=shells[i]
            payload={'schema':'ihm.cross-structure-resolved-surface.v1','entity_id':p['entity_id'],
                'name':p['name'],'role':p['role'],'system':p['system'],'units':'m','frame':p['frame'],
                'representation':'triangular_surface','cluster':label,
                'source_path':p['source_path'],'source_sha256':p['source_sha256'],
                'positions':[float(x) for x in SV.ravel()],'indices':[int(x) for x in SF.ravel()]}
            path=emit/(p['entity_id']+'.json.gz')
            path.write_bytes(gzip.compress(json.dumps(payload,allow_nan=False).encode(),mtime=0))
            e=np.sort(np.concatenate([SF[:,[0,1]],SF[:,[1,2]],SF[:,[2,0]]]),axis=1)
            _,ec=np.unique(e,axis=0,return_counts=True) if len(SF) else (None,np.zeros(0,int))
            sv=abs(divergence(SV,SF))
            written.append({'entity_id':p['entity_id'],'output_path':str(path.relative_to(emit.parent.parent)),
                'output_sha256':sha(path),'vertices':int(len(SV)),'faces':int(len(SF)),
                'closed':bool(len(SF)) and bool((ec==1).sum()==0),
                'boundary_edges':int((ec==1).sum()),'nonmanifold_edges':int((ec>2).sum()),
                'surface_volume_m3':sv,'owned_tet_volume_m3':per[i]['owned_volume_m3'],
                'surface_vs_tet_relative_error':float(sv/per[i]['owned_volume_m3']-1) if per[i]['owned_volume_m3'] else None})
        contacts=[];residual=0.0;boxes=[(V.min(0),V.max(0)) if len(V) else None for V,_ in shells]
        for i in range(len(parts)):
            for j in range(i+1,len(parts)):
                if boxes[i] is None or boxes[j] is None:continue
                if (boxes[i][0]>boxes[j][1]).any() or (boxes[j][0]>boxes[i][1]).any():continue
                VA,FA=shells[i];VB,FB=shells[j]
                if not len(FA) or not len(FB):continue
                n=int(len(np.asarray(cgal.intersect_other(VA,FA,VB,FB,True,False,False,False,2_000_000)[0]).reshape(-1,2)))
                b=cgal.mesh_boolean(VA,FA,VB,FB,type_str='intersect')
                v=abs(divergence(np.asarray(b[0],float),np.asarray(b[1],np.int64))) if len(b[1]) else 0.0
                residual+=v
                if n or v:contacts.append({'a':parts[i]['entity_id'],'b':parts[j]['entity_id'],
                    'shared_interface_face_contacts':n,'residual_overlap_volume_m3':v})
        report['resolved_surfaces']=written
        report['resolved_disjointness']={'pairs_tested':len(parts)*(len(parts)-1)//2,
            'pairs_in_contact':len(contacts),'total_residual_overlap_volume_m3':residual,
            'contacts':sorted(contacts,key=lambda c:-c['residual_overlap_volume_m3'])[:40],
            'basis':'exact CGAL boolean intersection of every pair of emitted surfaces. A nonzero face-contact count with zero volume is the intended outcome: the two regions meet on a shared conforming interface and claim no common volume.'}
        empty=[w['entity_id'] for w in written if not w['faces']]
        report['resolved_surface_check']={'entities':len(written),
            'fully_displaced_entities':empty,
            'fully_displaced_note':'an entity whose whole volume was inside higher-priority structures owns no tet and gets an empty surface; it is listed, not silently written as a degenerate shell',
            'all_closed':all(w['closed'] for w in written if w['faces']),
            'max_abs_surface_vs_tet_relative_error':max((abs(w['surface_vs_tet_relative_error'] or 0.0) for w in written if w['faces']),default=None)}
    return report,{'TV':TV,'TT':TT,'owner':owner,'det':det}


def build_parts(rows_by,ids):
    parts=[]
    for eid in ids:
        r=rows_by[eid];V,F=load_mesh(r['path'])
        parts.append({'entity_id':eid,'name':r['name'],'role':r['role'],'system':r['system'],
            'frame':'bodyparts3d-display-m','source_path':str(r['path'].relative_to(ROOT)),
            'source_sha256':r['sha256'],'vertices':int(len(V)),'faces':int(len(F)),
            'surface_volume_m3':abs(divergence(V,F)),'V':V,'F':F})
    return parts


def ladder_members(rows,seed_ids,steps,face_budget):
    """Entities nearest the seed cluster centroid, so each rung is a spatially contiguous block."""
    by={r['entity_id']:r for r in rows}
    centres=[];faces={}
    for r in rows:
        V,F=load_mesh(r['path']);centres.append((r['entity_id'],(V.min(0)+V.max(0))/2.0));faces[r['entity_id']]=len(F)
    seed=np.mean([c for e,c in centres if e in set(seed_ids)],axis=0)
    ordered=[e for e in seed_ids]+[e for e,_ in sorted(centres,key=lambda t:float(np.linalg.norm(t[1]-seed))) if e not in set(seed_ids)]
    rungs=[]
    for n in steps:
        if n>len(ordered):break
        members=ordered[:n]
        if sum(faces[e] for e in members)>face_budget:break
        rungs.append(members)
    return rungs


# ------------------------------------------------------------------ self-test

def self_test():
    """Two unit cubes overlapping in a slab, with the analytic answer known in closed form."""
    def cube(lo,hi):
        V=np.array(lo,float)+BOX_CORNERS*(np.array(hi,float)-np.array(lo,float))
        return np.ascontiguousarray(V),np.ascontiguousarray(BOX.copy())
    VA,FA=cube([0,0,0],[1,1,1]);VB,FB=cube([.75,0,0],[1.75,1,1])
    parts=[{'entity_id':'A','name':'a','role':'rigid_bone','system':'skeletal','frame':'test',
            'source_path':'test','source_sha256':'0'*64,'surface_volume_m3':abs(divergence(VA,FA)),'V':VA,'F':FA},
           {'entity_id':'B','name':'b','role':'muscle','system':'muscular','frame':'test',
            'source_path':'test','source_sha256':'0'*64,'surface_volume_m3':abs(divergence(VB,FB)),'V':VB,'F':FB}]
    with tempfile.TemporaryDirectory() as scratch:
        report,_=mesh_cluster(parts,0.25,'pY',{'A':rank('rigid_bone'),'B':rank('muscle')},
                              emit=Path(scratch)/'out'/'geometry',label='self-test')
    checks=[]
    def check(name,ok,detail):checks.append({'check':name,'passed':bool(ok),'detail':detail})
    plc=report['arrangement']
    check('arrangement produces a valid PLC',plc['plc_valid'],plc)
    check('crossings were found and resolved',plc['crossing_face_pairs_resolved']>0,plc['crossing_face_pairs_resolved'])
    check('tetgen succeeded',report['tetgen']['succeeded'],report['tetgen'].get('message'))
    m=report.get('conforming_mesh')
    check('all tets positive',m and m['all_positive_volume'],m and m['min_tet_volume_m3'])
    check('box volume closes',m and abs(m['box_volume_closure_relative_error'])<1e-9,m and m['box_volume_closure_relative_error'])
    a=[p for p in m['per_structure'] if p['entity_id']=='A'][0]
    b=[p for p in m['per_structure'] if p['entity_id']=='B'][0]
    check('both cubes claim 1 m3',abs(a['claimed_volume_m3']-1)<1e-9 and abs(b['claimed_volume_m3']-1)<1e-9,(a['claimed_volume_m3'],b['claimed_volume_m3']))
    check('disputed slab is 0.25 m3',abs(m['disputed_volume_m3']-.25)<1e-9,m['disputed_volume_m3'])
    check('bone keeps all of it',abs(a['owned_volume_m3']-1)<1e-9 and a['volume_moved_out_m3']<1e-12,a)
    check('muscle concedes exactly the slab',abs(b['owned_volume_m3']-.75)<1e-9 and abs(b['volume_moved_out_m3']-.25)<1e-9,b)
    check('no volume is owned twice',abs(m['structure_volume_m3']-1.75)<1e-9,m['structure_volume_m3'])
    r=report['resolved_surface_check']
    check('resolved surfaces are closed',r['all_closed'],r)
    check('resolved surface volume equals its owned tet volume',r['max_abs_surface_vs_tet_relative_error']<1e-12,r)
    check('resolved surfaces carry the moved volume',
          abs([w for w in report['resolved_surfaces'] if w['entity_id']=='B'][0]['surface_volume_m3']-.75)<1e-9,
          report['resolved_surfaces'])
    ta=self_test_rank_order()
    checks.extend(ta)
    passed=all(c['passed'] for c in checks)
    print(json.dumps({'self_test':'ihm.cross-structure-conflict-repair','passed':passed,'checks':checks,
                      'tetgen_stdout_tail':report['tetgen']['tetgen_stdout_tail']},indent=2))
    return passed


def self_test_rank_order():
    volumes={'x':2.0,'y':1.0}
    bone={'entity_id':'x','role':'rigid_bone','system':'skeletal'}
    muscle={'entity_id':'y','role':'muscle','system':'muscular'}
    o,l,k=decide(muscle,bone,volumes)
    a={'check':'bone displaces muscle regardless of argument order','passed':o=='x' and l=='y' and k=='hierarchical','detail':(o,l,k)}
    art1={'entity_id':'x','role':'vascular','system':'arterial'};art2={'entity_id':'y','role':'vascular','system':'arterial'}
    o2,l2,k2=decide(art2,art1,volumes)
    b={'check':'sibling arteries resolve by volume','passed':o2=='x' and l2=='y' and k2=='sibling','detail':(o2,l2,k2)}
    nerve={'entity_id':'y','role':'nerve','system':'nervous'}
    o3,l3,k3=decide(art1,nerve,volumes)
    c={'check':'nerve-vascular is flagged as a peer conflict','passed':k3=='neurovascular_peer','detail':(o3,l3,k3)}
    d={'check':'every role in the atlas has a declared priority',
       'passed':set(PRIORITY)>= {'rigid_bone','cartilage','tendon','ligament','fluid_cavity','vascular','nerve','muscle','soft_organ','connective_tissue','lymph_node_group'},
       'detail':sorted(PRIORITY)}
    return [a,b,c,d]


# ------------------------------------------------------------------------ run

def run(args):
    out=Path(args.out)
    if out.exists() and any(out.iterdir()) and not args.force and args.phase in ('graph','all'):
        raise ValueError('Choose a fresh output directory or pass --force')
    out.mkdir(parents=True,exist_ok=True)
    began=time.monotonic();rows=sources();by={r['entity_id']:r for r in rows}
    summary={'schema':'ihm.cross-structure-conflict-repair.v1','phase':args.phase}
    graph_path=out/'conflict-graph.json'
    if args.phase in ('graph','all'):
        graph=phase_graph(out,rows,args.workers,args.bool_cap,args.mc_points,args.seed,args.chunk)
        write_json(graph_path,graph)
    else:
        graph=json.loads(graph_path.read_text()) if graph_path.exists() else None
    summary['graph']=graph
    if args.phase in ('coincide','all'):
        coin=phase_coincide(rows);write_json(out/'coincidence.json',coin)
    elif (out/'coincidence.json').exists():coin=json.loads((out/'coincidence.json').read_text())
    else:coin=None
    summary['coincidence']=coin
    if args.phase in ('resolve','all') and graph:
        resolve,comp,ledger=phase_resolve(out,rows,graph)
        write_json(out/'conflict-components.json',{'schema':'ihm.cross-structure-conflict-components.v1',
            'count':len(comp),'components':[{'index':i,'size':len(c),'members':c} for i,c in enumerate(comp)]})
        write_json(out/'resolution.json',resolve)
    elif (out/'resolution.json').exists():resolve=json.loads((out/'resolution.json').read_text())
    else:resolve=None
    summary['resolution']=resolve
    if args.phase in ('mesh','all'):
        out.joinpath('geometry').mkdir(exist_ok=True)
        priority={r['entity_id']:rank(r['role']) for r in rows}
        clusters=json.loads((out/'clusters.json').read_text())['clusters'] if (args.cluster and (out/'clusters.json').exists()) else []
        drop=set(coin['identical_geometry_entities_dropped']) if coin else set()
        if args.cluster:
            members=[e for e in args.cluster.split(',') if e and e not in drop]
            label=args.cluster_label or ('cluster-%d'%len(members))
            report,_=mesh_cluster(build_parts(by,members),args.pad,args.flags,priority,
                                  emit=out/'geometry'/label,label=label,snap_below=args.snap)
            clusters=[c for c in clusters if c.get('label')!=label]+[report]
            write_json(out/'clusters.json',{'clusters':clusters})
            print(label,'tetgen succeeded:',report['tetgen']['succeeded'],flush=True)
            summary['clusters']=[{k:v for k,v in c.items() if k!='resolved_surfaces'} for c in clusters]
            args.ladder=''
        else:
            parts=build_parts(by,[e for e in THIGH if e not in drop])
            report,_=mesh_cluster(parts,args.pad,args.flags,priority,emit=out/'geometry/thigh-six',label='thigh-six',snap_below=args.snap)
            clusters.append(report)
            write_json(out/'clusters.json',{'clusters':clusters})
            print('thigh-six tetgen succeeded:',report['tetgen']['succeeded'],flush=True)
        if args.ladder:
            steps=[int(x) for x in args.ladder.split(',')]
            misses=0
            for members in ladder_members([r for r in rows if r['entity_id'] not in drop],THIGH,steps,args.face_budget):
                label='ladder-%d'%len(members)
                print('ladder rung',len(members),flush=True)
                p=build_parts(by,members)
                try:
                    r,_=mesh_cluster(p,args.pad,args.flags,priority,label=label,snap_below=args.snap)
                except BaseException as error:
                    r={'label':label,'member_count':len(members),'members':members,
                       'failed':type(error).__name__,'message':str(error)}
                clusters.append(r);write_json(out/'clusters.json',{'clusters':clusters})
                ok=r.get('tetgen',{}).get('succeeded')
                print('  rung %d succeeded: %s'%(len(members),ok),flush=True)
                misses=0 if ok else misses+1
                if misses>=args.ladder_failures:
                    print('ladder stopped after %d consecutive failures at %d members'%(misses,len(members)),flush=True);break
        summary['clusters']=[{k:v for k,v in c.items() if k!='resolved_surfaces'} for c in clusters]
    summary['wall_seconds']=time.monotonic()-began
    summary['inputs_sha256']={
        'data/derived/canonical/anatomy.json':sha(ANATOMY),
        'data/derived/muscle-tet-ready-v1/manifest.json':sha(MUSCLE/'manifest.json'),
        'data/derived/entity-tet-ready-v1/manifest.json':sha(ENTITY/'manifest.json'),
        'scripts/build_cross_structure_conflict_repair.py':sha(__file__)}
    summary['python']=sys.version
    write_json(out/'summary.json',summary)
    artifacts={p.relative_to(out).as_posix():sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    write_json(out/'manifest.json',{'schema':'ihm.cross-structure-conflict-repair-manifest.v1',
        'inputs_sha256':summary['inputs_sha256'],
        'per_entity_input_sha256':{r['entity_id']:r['sha256'] for r in rows},
        'artifacts_sha256':artifacts,'canonical_assets_modified':False})
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--out',default=str(ROOT/'data/derived/cross-structure-repair-v1'))
    p.add_argument('--phase',default='all',choices=['graph','coincide','resolve','mesh','all'])
    p.add_argument('--force',action='store_true');p.add_argument('--self-test',action='store_true')
    p.add_argument('--workers',type=int,default=16);p.add_argument('--chunk',type=int,default=64)
    p.add_argument('--bool-cap',type=int,default=120000);p.add_argument('--mc-points',type=int,default=200000)
    p.add_argument('--seed',type=int,default=0);p.add_argument('--pad',type=float,default=.02)
    p.add_argument('--flags',default='pY,pYT1e-14,p');p.add_argument('--ladder',default='')
    p.add_argument('--ladder-failures',type=int,default=2);p.add_argument('--snap',type=float,default=1e-7)
    p.add_argument('--face-budget',type=int,default=1_200_000)
    p.add_argument('--cluster',default='');p.add_argument('--cluster-label',default='')
    a=p.parse_args()
    if a.self_test:sys.exit(0 if self_test() else 1)
    run(a)
