"""Execute the recorded cross-structure ownership rule geometrically over the whole atlas.

data/derived/cross-structure-repair-v1 measured 11047 genuinely intersecting pairs among the 2403
repaired surfaces and recorded, for each, which structure owns the disputed volume. The decision was
executed geometrically only inside two meshed clusters (54 entities). This module executes it for
every pair and emits a conflict-free surface set, independent of any whole-body meshing.

  duplicates  five byte-identical id pairs are one surface authored twice; the later id is dropped
              rather than differenced against itself.
  resolve     a strict total order (declared role rank, then larger volume, then id) makes ownership
              globally acyclic. Entities are swept in that order, so when a structure is processed
              every higher-priority neighbour is already final and it is cut once against each.
              CGAL boolean minus first; an exact CGAL arrangement + winding classification is the
              fallback. The method is recorded per pair.
  imprint     the cut facets a loser inherits are subdivisions of its owner's facets. The owner is
              re-triangulated to carry the same subdivision so the interface has shared nodes, not
              merely coincident geometry. Guarded per owner: adopted only if the owner stays closed,
              manifold, self-intersection-free and volume-identical.
  verify      exhaustive CGAL intersect_other over every bounding-box-overlapping pair of the
              emitted set, with an exact boolean intersection volume for every pair that touches;
              plus per-entity closedness, manifoldness, self-intersection, orientation and
              volume-gated TetGen (harness imported from verify_muscle_tet_ready_surfaces).

Requires the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_conflict_free_surface_atlas.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_conflict_free_surface_atlas.py \
    --out data/derived/conflict-free-atlas-v1 --phase all

Nothing under data/derived/muscle-tet-ready-v1, data/derived/entity-tet-ready-v1 or
data/derived/cross-structure-repair-v1 is modified.
"""
from pathlib import Path
import argparse,ctypes,gzip,hashlib,json,multiprocessing as mp,os,resource,signal,sys,tempfile,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import igl
import igl.copyleft.cgal as cgal

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(Path(__file__).resolve().parent))
from build_cross_structure_conflict_repair import (PRIORITY,rank,sha,write_json,load_mesh,divergence,
                                                   sources,decide)
from build_muscle_tet_ready_surfaces import diagnose,signed_volume
from verify_muscle_tet_ready_surfaces import tetrahedralize,VOLUME_GATE

LIBC=ctypes.CDLL(None)
REPAIR=ROOT/'data/derived/cross-structure-repair-v1'
CUTOFF=20_000_000
BOX=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,5,1],[0,4,5],[1,6,2],[1,5,6],[2,7,3],[2,6,7],[3,4,0],[3,7,4]],np.int64)
BOX_CORNERS=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
SELF=-1


# ---------------------------------------------------------------- primitives

FACET_DTYPE=np.dtype([('c%d'%i,'f8') for i in range(9)])


def facet_key(V,F):
    """Exact key per facet, invariant to winding and vertex order: the lexicographically sorted
    coordinate triple. Sortable and hashable, so two surfaces can be intersected facet-wise."""
    if not len(F):return np.zeros(0,dtype=FACET_DTYPE)
    T=V[F]
    order=np.lexsort((T[:,:,2].T,T[:,:,1].T,T[:,:,0].T),axis=0)
    T=np.take_along_axis(T,order.T[:,:,None],axis=1)
    return np.ascontiguousarray(T).reshape(len(F),9).view(FACET_DTYPE).ravel()


def weld(V,F):
    if not len(F):return np.zeros((0,3)),np.zeros((0,3),np.int64)
    U,inv=np.unique(V,axis=0,return_inverse=True);inv=np.ravel(inv)
    G=np.ascontiguousarray(inv[F.ravel()].reshape(-1,3))
    keep=(G[:,0]!=G[:,1])&(G[:,1]!=G[:,2])&(G[:,0]!=G[:,2])
    G=G[keep]
    used=np.unique(G) if len(G) else np.zeros(0,np.int64)
    remap=np.full(len(U),-1,np.int64);remap[used]=np.arange(len(used))
    return np.ascontiguousarray(U[used]),np.ascontiguousarray(remap[G.ravel()].reshape(-1,3))


def compact(V,F,extra=None):
    if not len(F):return np.zeros((0,3)),np.zeros((0,3),np.int64),(np.zeros(0,np.int64) if extra is not None else None)
    used=np.unique(F);remap=np.full(len(V),-1,np.int64);remap[used]=np.arange(len(used))
    return (np.ascontiguousarray(V[used]),np.ascontiguousarray(remap[F.ravel()].reshape(-1,3)),
            None if extra is None else np.ascontiguousarray(extra))


def edge_counts(F):
    if not len(F):return 0,0
    e=np.sort(np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1)
    _,c=np.unique(e,axis=0,return_counts=True)
    return int((c==1).sum()),int((c>2).sum())


def component_volumes(V,F):
    if not len(F):return np.zeros(0)
    n,c=igl.facet_components(np.ascontiguousarray(F));c=np.asarray(c).ravel()
    out=[]
    for k in range(int(n)):
        S=F[c==k]
        out.append(float(np.einsum('ij,ij->i',V[S[:,0]],np.cross(V[S[:,1]],V[S[:,2]])).sum()/6.0))
    return np.array(out)


def boundary_loops(F):
    """Cycles of the directed edges a closed surface is missing, one per hole."""
    d=np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]])
    key=np.sort(d,axis=1)
    _,first,count=np.unique(key,axis=0,return_index=True,return_counts=True)
    lone=first[count==1]
    if not len(lone):return []
    nxt={}
    for a,b in d[lone]:
        if int(b) in nxt:return None
        nxt[int(b)]=int(a)
    loops=[];seen=set()
    for start in list(nxt):
        if start in seen:continue
        loop=[start];seen.add(start);cur=nxt[start]
        while cur!=start:
            if cur in seen or cur not in nxt:return None
            loop.append(cur);seen.add(cur);cur=nxt[cur]
        if len(loop)<3:return None
        loops.append(loop)
    return loops


def close_holes(V,F):
    """Fan-triangulate every hole from its own centroid. A shell with a hole has no interior, so it
    cannot be a boolean operand at all; the repair is minimal, guarded and recorded per entity."""
    loops=boundary_loops(F)
    if loops is None:return None,{'filled':False,'reason':'boundary edges do not form simple cycles'}
    if not loops:return (V,F),{'filled':False,'reason':'already closed'}
    verts=[V];faces=[F];n=len(V);added=0
    for loop in loops:
        c=V[np.asarray(loop)].mean(0)
        verts.append(c[None,:])
        patch=np.array([[n,loop[i],loop[(i+1)%len(loop)]] for i in range(len(loop))],np.int64)
        faces.append(patch);n+=1;added+=len(patch)
    W=np.ascontiguousarray(np.concatenate(verts));G=np.ascontiguousarray(np.concatenate(faces).astype(np.int64))
    b,nm=edge_counts(G)
    si=self_intersections(W,G)
    ok=b==0 and nm==0 and si==0
    return ((W,G) if ok else None),{'filled':bool(ok),'holes':len(loops),'facets_added':added,
        'boundary_edges_after':b,'nonmanifold_edges_after':nm,'self_intersecting_face_pairs_after':si,
        'volume_before_m3':abs(divergence(V,F)),'volume_after_m3':abs(divergence(W,G)),
        'reason':None if ok else 'the fan patch did not produce a clean closed shell'}


def self_intersections(V,F):
    if not len(F):return 0
    return int(len(np.asarray(cgal.remesh_self_intersections(V,F,True,False,False,False,CUTOFF)[2]).reshape(-1,2)))


# ------------------------------------------------------------------ boolean

def boolean_minus(VA,FA,VB,FB):
    """CGAL boolean difference. Returns (V,F,J) with J indexing the concatenated (FA;FB) soup."""
    r=cgal.mesh_boolean(np.ascontiguousarray(VA),np.ascontiguousarray(FA),
                        np.ascontiguousarray(VB),np.ascontiguousarray(FB),type_str='minus')
    return (np.ascontiguousarray(np.asarray(r[0],float)),
            np.ascontiguousarray(np.asarray(r[1],np.int64)),
            np.ascontiguousarray(np.asarray(r[2],np.int64)).ravel())


def _cyclic(F):
    """Winding-preserving canonical form: rotate each triple so the smallest index leads."""
    r=np.argmin(F,axis=1)
    return np.stack([F[np.arange(len(F)),(r+k)%3] for k in range(3)],axis=1)


def arrangement_minus(VA,FA,VB,FB):
    """Fallback: one exact CGAL arrangement of the two-surface soup, then side classification.

    A's facets whose interior side is outside B are kept; B's facets whose exterior side is inside A
    are kept with reversed winding, so the cut boundary A inherits is literally B's own subdivided
    geometry. Coincident facets are settled combinatorially rather than by a winding number, because
    a sample point on a coincident facet sits exactly on the other surface."""
    V=np.ascontiguousarray(np.concatenate([VA,VB]))
    F=np.ascontiguousarray(np.concatenate([FA,FB+len(VA)]).astype(np.int64))
    r=cgal.remesh_self_intersections(V,F,False,False,True,False,CUTOFF)
    AV=np.ascontiguousarray(np.asarray(r[0],float));AF=np.ascontiguousarray(np.asarray(r[1],np.int64))
    J=np.ascontiguousarray(np.asarray(r[3],np.int64)).ravel()
    U,inv=np.unique(AV,axis=0,return_inverse=True);inv=np.ravel(inv)
    G=np.ascontiguousarray(inv[AF.ravel()].reshape(-1,3))
    keep=(G[:,0]!=G[:,1])&(G[:,1]!=G[:,2])&(G[:,0]!=G[:,2])
    G=G[keep];J=J[keep]
    if len(G):
        n=np.cross(U[G[:,1]]-U[G[:,0]],U[G[:,2]]-U[G[:,0]])
        nz=np.any(n!=0,axis=1);G=G[nz];J=J[nz];n=n[nz]
    if not len(G):return np.zeros((0,3)),np.zeros((0,3),np.int64),np.zeros(0,np.int64)
    froma=J<len(FA)
    live=np.ones(len(G),bool)
    _,ik=np.unique(np.sort(G,axis=1),axis=0,return_inverse=True);ik=np.ravel(ik)
    counts=np.bincount(ik)
    if (counts>1).any():
        cyc=_cyclic(G)
        buckets={}
        for i in np.flatnonzero(counts[ik]>1):buckets.setdefault(int(ik[i]),[]).append(int(i))
        for members in buckets.values():
            aside=[i for i in members if froma[i]];bside=[i for i in members if not froma[i]]
            if not aside or not bside:continue
            for j in bside:
                live[j]=False
                for i in aside:
                    if live[i] and np.array_equal(cyc[i],cyc[j]):live[i]=False;break
    length=np.linalg.norm(n,axis=1)
    unit=n/np.maximum(length,np.finfo(float).tiny)[:,None]
    e=np.linalg.norm(U[G[:,1]]-U[G[:,0]],axis=1)+np.linalg.norm(U[G[:,2]]-U[G[:,1]],axis=1)+np.linalg.norm(U[G[:,0]]-U[G[:,2]],axis=1)
    delta=.25*length/np.maximum(e,np.finfo(float).tiny)
    centre=U[G].mean(1)
    ai=np.flatnonzero(live&froma);bi=np.flatnonzero(live&~froma)
    take=np.zeros(len(G),bool)
    if len(ai):
        q=np.ascontiguousarray(centre[ai]-unit[ai]*delta[ai,None])
        take[ai]=np.abs(np.asarray(igl.fast_winding_number(np.ascontiguousarray(VB),np.ascontiguousarray(FB),q),float))<=.5
    if len(bi):
        q=np.ascontiguousarray(centre[bi]+unit[bi]*delta[bi,None])
        take[bi]=np.abs(np.asarray(igl.fast_winding_number(np.ascontiguousarray(VA),np.ascontiguousarray(FA),q),float))>.5
    out=np.concatenate([G[take&froma],G[take&~froma][:,[0,2,1]]])
    src=np.concatenate([J[take&froma],J[take&~froma]])
    V2,F2,_=compact(U,np.ascontiguousarray(out))
    return V2,F2,src


def cut(VA,FA,origin,VB,FB,owner_faces):
    """Difference B out of A, carrying per-facet provenance.

    B is the disjoint union of every higher-priority neighbour, so a structure is cut once against
    all of them rather than once per neighbour; `owner_faces` labels each B facet with its owner."""
    combined=np.concatenate([origin,owner_faces])
    try:
        V,F,J=boolean_minus(VA,FA,VB,FB)
        if len(F) and (J.max()>=len(combined) or J.min()<0):raise ValueError('birth index out of range')
        return V,F,(combined[J] if len(F) else np.zeros(0,np.int64)),'cgal_boolean_minus'
    except BaseException as boolean_error:
        try:
            V,F,J=arrangement_minus(VA,FA,VB,FB)
            return V,F,(combined[J] if len(F) else np.zeros(0,np.int64)),'cgal_arrangement:'+type(boolean_error).__name__
        except BaseException as arrangement_error:
            return None,None,None,'failed:%s|%s'%(type(boolean_error).__name__,type(arrangement_error).__name__)


def isolated_cut(VA,FA,origin,VB,FB,owner_faces,mem_bytes,timeout_s):
    """Run one difference in a forked child under an address-space cap and a wall-clock deadline.

    An exact boolean on two anatomical surfaces can blow up without warning; a child that overruns
    is killed and recorded as a refused operation instead of taking the host down with it."""
    handle,store=tempfile.mkstemp(suffix='.npz');os.close(handle);os.unlink(store)
    handle,noise=tempfile.mkstemp(suffix='.cut');os.close(handle)
    reader,writer=os.pipe();began=time.monotonic()
    sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
    pid=os.fork()
    if pid==0:
        try:
            os.close(reader)
            sink=os.open(noise,os.O_WRONLY)
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None);os.dup2(sink,1);os.dup2(sink,2)
            if mem_bytes:resource.setrlimit(resource.RLIMIT_AS,(mem_bytes,mem_bytes))
            V,F,O,method=cut(VA,FA,origin,VB,FB,owner_faces)
            if V is not None:np.savez(store,V=V,F=F,O=O)
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
            os.write(writer,json.dumps({'method':method,'ok':V is not None}).encode())
            os.close(writer)
        except BaseException as error:
            try:os.write(writer,json.dumps({'method':'failed:'+type(error).__name__,'ok':False}).encode())
            except BaseException:pass
        finally:os._exit(0)
    os.close(writer)
    deadline=began+timeout_s;killed=False
    while True:
        finished,_=os.waitpid(pid,os.WNOHANG)
        if finished:break
        if time.monotonic()>deadline:
            os.kill(pid,signal.SIGKILL);os.waitpid(pid,0);killed=True;break
        time.sleep(.05)
    chunks=[]
    try:
        while True:
            chunk=os.read(reader,65536)
            if not chunk:break
            chunks.append(chunk)
    finally:os.close(reader)
    seconds=time.monotonic()-began
    tail=Path(noise).read_text(errors='replace').strip().splitlines()[-6:]
    Path(noise).unlink(missing_ok=True)
    if killed:
        Path(store).unlink(missing_ok=True)
        return None,None,None,'refused_over_time_cap:%ds'%timeout_s,seconds
    if not chunks:
        Path(store).unlink(missing_ok=True)
        return None,None,None,'refused_child_died:'+('|'.join(tail)[:120] or 'no message; probably the address-space cap'),seconds
    payload=json.loads(b''.join(chunks).decode())
    if not payload['ok']:
        Path(store).unlink(missing_ok=True)
        return None,None,None,payload['method'],seconds
    data=np.load(store);V=data['V'];F=data['F'];O=data['O'];data.close()
    Path(store).unlink(missing_ok=True)
    return (np.ascontiguousarray(V),np.ascontiguousarray(F,dtype=np.int64),
            np.ascontiguousarray(O,dtype=np.int64),payload['method'],seconds)


# ------------------------------------------------------------------- inputs

def entity_digest(V,F):
    U,inv=np.unique(V,axis=0,return_inverse=True);inv=np.ravel(inv)
    T=np.unique(np.sort(inv[F.ravel()].reshape(-1,3),axis=1),axis=0)
    h=hashlib.sha256();h.update(np.ascontiguousarray(U).tobytes());h.update(np.ascontiguousarray(T).tobytes())
    return h.hexdigest()


def load_inputs():
    rows=sources()
    meshes={};volume={};digest={}
    for r in rows:
        V,F=load_mesh(r['path']);meshes[r['entity_id']]=(V,F)
        volume[r['entity_id']]=abs(divergence(V,F));digest[r['entity_id']]=entity_digest(V,F)
    return rows,meshes,volume,digest


def duplicate_groups(rows,digest):
    by={}
    for r in rows:by.setdefault(digest[r['entity_id']],[]).append(r['entity_id'])
    groups=sorted([sorted(v) for v in by.values() if len(v)>1])
    return groups,sorted(e for g in groups for e in g[1:])


def load_edges(drop,keep_of):
    edges=[json.loads(l) for l in (REPAIR/'conflict-pairs.jsonl').read_text().splitlines()]
    edges=[e for e in edges if e.get('intersecting_face_pairs',0)>0]
    kept=[];removed=[];rewired=[];seen=set()
    for e in edges:
        a,b=keep_of.get(e['a'],e['a']),keep_of.get(e['b'],e['b'])
        if a==b:
            removed.append({'a':e['a'],'b':e['b'],'reason':'byte-identical duplicate pair, not a conflict',
                            'intersecting_face_pairs':e['intersecting_face_pairs'],
                            'overlap_volume_m3':e.get('overlap_volume_m3')});continue
        key=tuple(sorted((a,b)))
        if (e['a'],e['b'])!=(a,b) or (e['a'],e['b'])!=key:
            if e['a'] in drop or e['b'] in drop:
                rewired.append({'from':[e['a'],e['b']],'to':list(key),'already_present':key in seen})
        if key in seen:continue
        seen.add(key)
        kept.append({'a':key[0],'b':key[1],'intersecting_face_pairs':e['intersecting_face_pairs'],
                     'overlap_volume_m3':e.get('overlap_volume_m3'),
                     'overlap_thickness_proxy_m':e.get('overlap_thickness_proxy_m'),
                     'overlap_method':e.get('overlap_method')})
    return kept,removed,rewired


def order_keys(rows,volume,active):
    return {r['entity_id']:(rank(r['role']),-volume[r['entity_id']],r['entity_id'])
            for r in rows if r['entity_id'] in active}


def layers(active,edges,key):
    owners={e:[] for e in active}
    for r in edges:
        a,b=r['a'],r['b']
        if key[a]<key[b]:owners[b].append(a)
        else:owners[a].append(b)
    level={};order=sorted(active,key=lambda e:key[e])
    for e in order:
        level[e]=0 if not owners[e] else 1+max(level[o] for o in owners[e])
    buckets={}
    for e in order:buckets.setdefault(level[e],[]).append(e)
    return owners,level,[buckets[k] for k in sorted(buckets)]


# ------------------------------------------------------------------ census

CSTATE={}


def _census_task(task):
    out=[]
    for i,j in task:
        a,b=CSTATE['ids'][i],CSTATE['ids'][j]
        VA,FA=CSTATE['meshes'][a];VB,FB=CSTATE['meshes'][b]
        rec={'a':a,'b':b}
        try:
            n=int(len(np.asarray(cgal.intersect_other(VA,FA,VB,FB,True,False,False,False,2_000_000)[0]).reshape(-1,2)))
        except BaseException as error:
            rec['error']='intersect_other:'+type(error).__name__;out.append(rec);continue
        rec['intersecting_face_pairs']=n
        if not n:
            # No surface crossing means each connected component is wholly in or wholly out, so
            # containment is the only remaining way to share volume, and it is invisible to a
            # surface-crossing test.
            inside=False
            for (V1,F1),(V2,F2) in ((( VA,FA),(VB,FB)),((VB,FB),(VA,FA))):
                if not len(F1) or not len(F2):continue
                count,comp=igl.facet_components(np.ascontiguousarray(F1))
                comp=np.asarray(comp).ravel()
                seeds=[]
                for k in range(int(count)):
                    sel=F1[comp==k]
                    if len(sel):seeds.append(V1[sel[0]].mean(0))
                if not seeds:continue
                q=np.ascontiguousarray(np.asarray(seeds,float))
                w=np.abs(np.asarray(igl.fast_winding_number(np.ascontiguousarray(V2),np.ascontiguousarray(F2),q),float))
                if (w>.5).any():inside=True;break
            rec['containment']=bool(inside)
            if not inside:continue
        try:
            r=cgal.mesh_boolean(VA,FA,VB,FB,type_str='intersect')
            IV=np.asarray(r[0],float);IFc=np.asarray(r[1],np.int64)
            rec['overlap_volume_m3']=abs(divergence(IV,IFc)) if len(IFc) else 0.0
            rec['overlap_method']='cgal_boolean_intersect'
        except BaseException as error:
            rec['overlap_volume_m3']=None;rec['overlap_method']='cgal_boolean_failed:'+type(error).__name__
        out.append(rec)
    return out


def census(meshes,ids,workers,chunk,log):
    """Re-run the exhaustive pairwise conflict test on the surfaces actually being resolved.

    The recorded census gated on CGAL intersect_other, which cannot see a structure lying wholly
    inside another: a nested pair crosses no face. It also predates the duplicate drop and the hole
    fill, both of which change which surfaces exist and what volume they bound."""
    began=time.monotonic()
    lo=np.stack([meshes[e][0].min(0) for e in ids]);hi=np.stack([meshes[e][0].max(0) for e in ids])
    ov=(lo[:,None,:]<=hi[None,:,:]).all(2)&(hi[:,None,:]>=lo[None,:,:]).all(2)
    np.fill_diagonal(ov,False)
    idx=np.argwhere(np.triu(ov))
    faces=np.array([len(meshes[e][1]) for e in ids])
    idx=idx[np.argsort(faces[idx[:,0]]+faces[idx[:,1]])]
    pairs=[(int(i),int(j)) for i,j in idx]
    CSTATE.update(meshes=meshes,ids=ids)
    tasks=[pairs[k:k+chunk] for k in range(0,len(pairs),chunk)]
    log('  census: %d entities, %d bounding-box-overlapping pairs, %d tasks'%(len(ids),len(pairs),len(tasks)))
    found,dead=run_tasks(_census_task,tasks,workers,log,'census tasks',True)
    crossing=[r for r in found if r.get('intersecting_face_pairs',0)>0]
    nested=[r for r in found if r.get('intersecting_face_pairs',0)==0 and r.get('containment')]
    report={'entities':len(ids),'bounding_box_overlapping_pairs':len(pairs),
        'coverage':'exhaustive over every bounding-box-overlapping pair: CGAL intersect_other for surface crossing, then a winding-number containment test per connected component for the pairs that cross no face',
        'crossing_pairs':len(crossing),'nested_pairs':len(nested),
        'conflicting_pairs':len(crossing)+len(nested),
        'errored_pairs':[r for r in found if 'error' in r],
        'tasks_whose_worker_died':len(dead),
        'seconds':time.monotonic()-began}
    return [r for r in found if r.get('intersecting_face_pairs',0)>0 or r.get('containment')],report


# ------------------------------------------------------------ resolve sweep

STATE={}


def _resolve_task(payload):
    began=time.monotonic()
    entity_id=payload['entity_id'];V,F=payload['V'],payload['F']
    operand={o:(OV,OF) for o,OV,OF in payload['owners']}
    index={o:i for o,i,_ in payload['labels']};label_to_id={i:o for o,i,_ in payload['labels']}
    weight={o:w for o,_,w in payload['labels']}
    mem=payload['mem'];timeout=payload['timeout'];budget=payload['budget']
    final=operand
    before=abs(divergence(V,F))
    origin=np.full(len(F),SELF,np.int64)
    ops={};usable=[]
    todo=sorted(operand,key=lambda o:(-(weight.get(o,0.0) or 0.0),o))
    for owner in todo:
        OV,OF=operand[owner]
        ops[owner]={'a':entity_id,'b':owner,'owner':owner,'yields':entity_id,
                    'recorded_overlap_volume_m3':weight.get(owner),
                    'interface_facets_inherited':0}
        if not len(OF):
            ops[owner]['method']='skipped_empty_owner';continue
        boundary=edge_counts(OF)[0]
        if boundary:
            ops[owner]['method']='skipped_open_owner';ops[owner]['owner_boundary_edges']=int(boundary)
            ops[owner]['note']='an open shell bounds no volume, so there is nothing to subtract'
            continue
        usable.append(owner)
    was_closed=edge_counts(F)[0]==0
    attempts=[];spent=0.0
    if usable and len(F):
        lo=V.min(0);hi=V.max(0)
        near=[o for o in usable
              if not ((final[o][0].min(0)>hi).any() or (final[o][0].max(0)<lo).any())]
        for o in usable:
            if o not in near:ops[o]['method']='skipped_disjoint_bounding_boxes'
        if near:
            BV=[];BF=[];lab=[];offset=0
            for o in near:
                OV,OF=final[o]
                BV.append(OV);BF.append(OF+offset);lab.append(np.full(len(OF),index[o],np.int64));offset+=len(OV)
            BV=np.ascontiguousarray(np.concatenate(BV));BF=np.ascontiguousarray(np.concatenate(BF).astype(np.int64))
            lab=np.concatenate(lab)
            NV,NF,NO,method,seconds=isolated_cut(V,F,origin,BV,BF,lab,mem,timeout)
            spent+=seconds
            attempts.append({'strategy':'union-of-owners','owners':len(near),'owner_facets':int(len(BF)),
                             'method':method,'seconds':seconds,'succeeded':NV is not None})
            if NV is not None and was_closed and len(NF) and edge_counts(NF)[0]:
                NV=None;attempts[-1].update(succeeded=False,method=method+'|reverted_open_result')
            if NV is not None:
                V,F,origin=NV,NF,NO
                for o in near:ops[o]['method']=method
            else:
                # one bad owner can poison the union, so fall back to cutting them one at a time,
                # under a shared budget: an entity must not be able to consume the run.
                refused=0
                for o in near:
                    left=budget-spent
                    if not len(F):ops[o]['method']='skipped_empty_operand';continue
                    if left<=5.0:
                        ops[o]['method']='refused_entity_time_budget_exhausted';refused+=1;continue
                    OV,OF=final[o]
                    lab1=np.full(len(OF),index[o],np.int64)
                    NV,NF,NO,m,sec=isolated_cut(V,F,origin,OV,OF,lab1,mem,min(timeout,left))
                    spent+=sec
                    if NV is not None and was_closed and len(NF) and edge_counts(NF)[0]:
                        NV=None;m=m+'|reverted_open_result'
                    if NV is not None:V,F,origin=NV,NF,NO
                    ops[o]['method']='sequential|'+m;ops[o]['seconds']=sec
                attempts.append({'strategy':'sequential-fallback','owners':len(near),
                                 'refused_for_budget':refused,'seconds':spent,'budget_s':budget})
    if len(origin):
        for oi,count in zip(*np.unique(origin[origin!=SELF],return_counts=True)):
            owner=label_to_id.get(int(oi))
            if owner in ops:ops[owner]['interface_facets_inherited']=int(count)
    after=abs(divergence(V,F))
    rows=list(ops.values())
    for r in rows:r.setdefault('method','not_attempted')
    return {'entity_id':entity_id,'V':V,'F':F,'origin':origin,'ops':rows,'attempts':attempts,
            'volume_before_m3':before,'volume_after_m3':after,'faces_after':int(len(F)),
            'seconds':time.monotonic()-began}


def resolve(active,meshes,edges,owners,index,workers,log,mem,timeout,budget):
    """Cut every structure once, against the union of the repaired surfaces of the higher-priority
    neighbours it yields to.

    The cut uses each owner's ORIGINAL surface, not its resolved one, and that is not an
    approximation: a resolved owner is a subset of its original, so a result disjoint from the
    original is disjoint from the resolved surface too. The only volume the two conventions could
    assign differently is volume in E, O and a third structure H that displaced O; but then E and H
    themselves overlap, the exhaustive pairwise test recorded that as a conflict, and E is cut by H
    as well. So the assignment is identical and every entity becomes independent, which matters:
    ownership chains here run 60 deep, and cutting against resolved owners serialises the atlas
    behind them."""
    weight={tuple(sorted((r['a'],r['b']))):(r['overlap_volume_m3'] or 0.0) for r in edges}
    final={};results={};ops=[];work=[]
    for e in active:
        if not owners[e]:
            V,F=meshes[e];v=abs(divergence(V,F))
            final[e]=(V,F)
            results[e]={'entity_id':e,'origin':np.full(len(F),SELF,np.int64),
                        'volume_before_m3':v,'volume_after_m3':v,'faces_after':int(len(F)),
                        'attempts':[],'seconds':0.0}
        else:work.append(e)
    order=sorted(work,key=lambda e:-sum(len(meshes[o][1]) for o in owners[e]))

    def build(e):
        return {'entity_id':e,'V':meshes[e][0],'F':meshes[e][1],
                'owners':[(o,meshes[o][0],meshes[o][1]) for o in owners[e]],
                'labels':[(o,index[o],weight.get(tuple(sorted((e,o))),0.0)) for o in owners[e]],
                'mem':mem,'timeout':timeout,'budget':budget}

    began=time.monotonic();done=0;spoke=0.0
    if workers>1 and len(order)>1:
        with mp.get_context('fork').Pool(workers) as pool:
            for r in pool.imap_unordered(_resolve_task,(build(e) for e in order),chunksize=1):
                final[r['entity_id']]=(r['V'],r['F']);ops.extend(r['ops']);done+=1
                results[r['entity_id']]={k:v for k,v in r.items() if k not in ('V','F','ops')}
                if time.monotonic()-spoke>60:
                    spoke=time.monotonic()
                    log('  sweep: %d/%d cut, %.0f s'%(done,len(order),time.monotonic()-began))
    else:
        for e in order:
            r=_resolve_task(build(e))
            final[r['entity_id']]=(r['V'],r['F']);ops.extend(r['ops']);done+=1
            results[r['entity_id']]={k:v for k,v in r.items() if k not in ('V','F','ops')}
    log('  sweep: %d/%d cut in %.1f s'%(done,len(order),time.monotonic()-began))
    missing=[e for e in active if e not in final]
    if missing:raise RuntimeError('%d entities were never cut: %s'%(len(missing),missing[:5]))
    return final,results,ops


# -------------------------------------------------------------------- imprint

def _imprint_worker(payload):
    owner,OV,OF,patches=payload['owner'],payload['V'],payload['F'],payload['patches']
    origin=payload.get('origin')
    base=abs(divergence(OV,OF))
    faces=sum(len(p[1]) for p in patches)
    record={'entity_id':owner,'patch_facets':faces,'losers':len(patches),
            'faces_before':int(len(OF)),'volume_m3':base}
    if not faces or not len(OF):
        record.update(adopted=False,reason='no patch');return record,None
    t0=time.monotonic()
    try:
        V=[OV];F=[OF];offset=len(OV)
        for PV,PF in patches:
            V.append(PV);F.append(PF[:,[0,2,1]]+offset);offset+=len(PV)
        SV=np.ascontiguousarray(np.concatenate(V));SF=np.ascontiguousarray(np.concatenate(F).astype(np.int64))
        r=cgal.remesh_self_intersections(SV,SF,False,False,True,False,CUTOFF)
        AV=np.ascontiguousarray(np.asarray(r[0],float));AF=np.ascontiguousarray(np.asarray(r[1],np.int64))
        J=np.ascontiguousarray(np.asarray(r[3],np.int64)).ravel()
        U,inv=np.unique(AV,axis=0,return_inverse=True);inv=np.ravel(inv)
        G=np.ascontiguousarray(inv[AF.ravel()].reshape(-1,3))
        alive=(G[:,0]!=G[:,1])&(G[:,1]!=G[:,2])&(G[:,0]!=G[:,2])
        G=G[alive];J=J[alive];W=U
        # the owner's own copies lead the soup, so the first occurrence of a duplicated facet is the
        # owner's and the coincident patch copy is the one dropped.
        _,first=np.unique(np.sort(G,axis=1),axis=0,return_index=True)
        take=np.sort(first);G=np.ascontiguousarray(G[take]);J=J[take]
        newborn=(np.asarray(origin,np.int64)[np.minimum(J,len(OF)-1)] if origin is not None
                 else np.full(len(G),SELF,np.int64))
        newborn=np.where(J<len(OF),newborn,SELF).astype(np.int64)
        used=np.unique(G);remap=np.full(len(W),-1,np.int64);remap[used]=np.arange(len(used))
        W=np.ascontiguousarray(W[used]);G=np.ascontiguousarray(remap[G.ravel()].reshape(-1,3))
        vol=abs(divergence(W,G));boundary,nonmanifold=edge_counts(G)
        residual=self_intersections(W,G)
        ok=(boundary==0 and nonmanifold==0 and residual==0 and base>0
            and abs(vol/base-1)<1e-9 and len(G)>0)
        record.update(faces_after=int(len(G)),volume_after_m3=vol,
                      relative_volume_error=float(vol/base-1) if base>0 else None,
                      boundary_edges=boundary,nonmanifold_edges=nonmanifold,
                      residual_self_intersecting_face_pairs=residual,
                      seconds=time.monotonic()-t0,adopted=bool(ok))
        record['facets_inherited_from_a_higher_priority_owner']=int((newborn!=SELF).sum())
        record['facets_the_arrangement_could_not_attribute']=int((J>=len(OF)).sum())
        if not ok:record['reason']='guard rejected the imprint; the owner is kept unrefined'
        return record,((W,G,newborn) if ok else None)
    except BaseException as error:
        record.update(adopted=False,reason='%s: %s'%(type(error).__name__,str(error)[:200]),
                      seconds=time.monotonic()-t0)
        return record,None


def imprint(final,results,index_to_id,workers,cap,log):
    patches={}
    for eid,r in results.items():
        origin=r['origin'];V,F=final[eid]
        if not len(F) or (origin==SELF).all():continue
        for oi in np.unique(origin[origin!=SELF]):
            owner=index_to_id[int(oi)]
            sel=F[origin==oi]
            PV,PF,_=compact(V,np.ascontiguousarray(sel))
            patches.setdefault(owner,[]).append((PV,PF))
    tasks=[]
    for owner,ps in sorted(patches.items()):
        total=sum(len(p[1]) for p in ps)+len(final[owner][1])
        if total>cap:
            tasks.append({'owner':owner,'V':np.zeros((0,3)),'F':np.zeros((0,3),np.int64),'patches':[],
                          'skipped':{'entity_id':owner,'adopted':False,'patch_facets':int(total),
                                     'reason':'soup over --imprint-cap facets'}})
            continue
        tasks.append({'owner':owner,'V':final[owner][0],'F':final[owner][1],'patches':ps,
                      'origin':results[owner]['origin']})
    live=[t for t in tasks if not t.get('skipped')]
    records=[t['skipped'] for t in tasks if t.get('skipped')]
    log('  imprint: %d owners carry a cut patch, %d attempted'%(len(patches),len(live)))
    if live:
        if workers>1 and len(live)>1:
            with mp.get_context('fork').Pool(min(workers,len(live))) as pool:
                got=list(pool.imap_unordered(_imprint_worker,live,chunksize=1))
        else:
            got=[_imprint_worker(t) for t in live]
        for record,mesh in got:
            records.append(record)
            if mesh is not None:
                final[record['entity_id']]=(mesh[0],mesh[1])
                results[record['entity_id']]['origin']=mesh[2]
    return final,sorted(records,key=lambda r:r['entity_id'])


# ------------------------------------------------------------------- emit

def emit(out,final,rows_by,results,index_to_id):
    geometry=out/'geometry';geometry.mkdir(parents=True,exist_ok=True)
    written=[]
    for eid in sorted(final):
        V,F=final[eid];r=rows_by[eid]
        payload={'schema':'ihm.conflict-free-surface.v1','entity_id':eid,'name':r['name'],'role':r['role'],
                 'system':r['system'],'units':'m','frame':'bodyparts3d-display-m',
                 'representation':'triangular_surface','priority_rank':rank(r['role']),
                 'source_path':str(Path(r['path']).relative_to(ROOT)),'source_sha256':r['sha256'],
                 'positions':[float(x) for x in V.ravel()],'indices':[int(x) for x in F.ravel()]}
        path=geometry/(eid+'.json.gz')
        path.write_bytes(gzip.compress(json.dumps(payload,allow_nan=False).encode(),mtime=0))
        written.append({'entity_id':eid,'output_path':path.relative_to(out).as_posix(),
                        'output_sha256':sha(path),'vertices':int(len(V)),'faces':int(len(F))})
    return written


# ------------------------------------------------------------------ verify

VSTATE={}


def run_tasks(fn,tasks,workers,log,label,expand,depth=0):
    """Fan out over processes without letting a native abort hang the run.

    A CGAL call can terminate its process outright. multiprocessing.Pool then waits forever for a
    result that will never arrive, so the executor is used instead: a dead worker fails every
    outstanding future, the survivors are retried one task at a time, and whatever still kills a
    worker is recorded as an errored task rather than silently lost."""
    collected=[];retry=[]
    if not tasks:return collected,[]
    executor=ProcessPoolExecutor(max_workers=max(1,min(workers,len(tasks))),mp_context=mp.get_context('fork'))
    try:
        futures={executor.submit(fn,t):t for t in tasks}
        done=0
        for future in as_completed(futures):
            try:
                value=future.result()
                collected.extend(value) if expand else collected.append(value)
            except BaseException as error:
                retry.append((futures[future],type(error).__name__))
            done+=1
            if done%50==0:log('    %s %d/%d, %d collected, %d to retry'%(label,done,len(tasks),len(collected),len(retry)))
    finally:
        try:executor.shutdown(wait=False,cancel_futures=True)
        except BaseException:pass
    if not retry:return collected,[]
    log('    %s: %d tasks lost a worker, retrying them singly'%(label,len(retry)))
    if depth>=2:return collected,[t for t,_ in retry]
    singles=[]
    for task,_ in retry:
        singles.extend([[x] for x in task] if isinstance(task,list) else [task])
    more,dead=run_tasks(fn,singles,workers,log,label+'-retry',expand,depth+1)
    collected.extend(more)
    return collected,dead


def isolated_tetgen(V,F,flags,timeout_s):
    """The imported TetGen harness has no wall clock, and a resolved surface can keep it busy for
    hours. It is run in its own session so the whole process group can be killed on the deadline."""
    handle,noise=tempfile.mkstemp(suffix='.tet');os.close(handle)
    reader,writer=os.pipe();began=time.monotonic()
    sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
    pid=os.fork()
    if pid==0:
        try:
            os.close(reader);os.setsid()
            sink=os.open(noise,os.O_WRONLY);os.dup2(sink,1);os.dup2(sink,2)
            outcome=tetrahedralize(np.ascontiguousarray(V),np.ascontiguousarray(F),flags)
            outcome['tetgen_stdout_tail']=outcome.pop('tetgen_stdout',[])[-8:]
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
            os.write(writer,json.dumps(outcome,allow_nan=False).encode());os.close(writer)
        except BaseException as error:
            try:os.write(writer,json.dumps({'succeeded':False,'exception':type(error).__name__}).encode())
            except BaseException:pass
        finally:os._exit(0)
    os.close(writer)
    deadline=began+timeout_s;killed=False
    while True:
        finished,_=os.waitpid(pid,os.WNOHANG)
        if finished:break
        if time.monotonic()>deadline:
            try:os.killpg(pid,signal.SIGKILL)
            except BaseException:
                try:os.kill(pid,signal.SIGKILL)
                except BaseException:pass
            os.waitpid(pid,0);killed=True;break
        time.sleep(.05)
    chunks=[]
    try:
        while True:
            chunk=os.read(reader,1<<20)
            if not chunk:break
            chunks.append(chunk)
    finally:os.close(reader)
    Path(noise).unlink(missing_ok=True)
    seconds=time.monotonic()-began
    if killed:
        return {'succeeded':False,'exception':'TimeCap','flags':flags,'seconds':seconds,
                'message':'TetGen was still running after %d s and was killed; the surface is not proven to tetrahedralise, and it is not proven to fail either'%timeout_s}
    if not chunks:
        return {'succeeded':False,'exception':'NoResult','flags':flags,'seconds':seconds,
                'message':'the TetGen child exited without a result'}
    outcome=json.loads(b''.join(chunks).decode());outcome['seconds']=seconds
    return outcome


def _entity_check(entity_id):
    V,F=VSTATE['final'][entity_id]
    before=VSTATE['before'][entity_id]
    record={'entity_id':entity_id,'role':VSTATE['role'][entity_id],'system':VSTATE['system'][entity_id],
            'priority_rank':rank(VSTATE['role'][entity_id]),
            'volume_before_m3':before,'volume_after_m3':abs(divergence(V,F)),
            'vertices':int(len(V)),'faces':int(len(F))}
    record['volume_conceded_m3']=before-record['volume_after_m3']
    record['fraction_conceded']=(record['volume_conceded_m3']/before) if before>0 else None
    if not len(F):
        record.update(closed=False,manifold=False,self_intersection_free=True,orientation_consistent=False,
                      annihilated=True,components=0,boundary_edges=0,nonmanifold_edges=0,
                      nonmanifold_vertices=0,self_intersecting_face_pairs=0,
                      tetgen={'succeeded':False,'reason':'no geometry'},still_tet_ready=False)
        return record
    d=diagnose(V,F,intersections=True)
    comp=component_volumes(V,F)
    record.update(annihilated=False,closed=bool(d['closed']),
                  boundary_edges=d['boundary_edges'],nonmanifold_edges=d['nonmanifold_edges'],
                  nonmanifold_vertices=d['nonmanifold_vertices'],components=d['face_components'],
                  manifold=bool(d['boundary_edges']==0 and d['nonmanifold_edges']==0 and d['nonmanifold_vertices']==0),
                  self_intersecting_face_pairs=d['self_intersecting_face_pairs'],
                  self_intersection_free=bool(d['self_intersection_free']),
                  orientation_consistent=bool(d['orientation_consistent']),
                  all_components_outward=bool(len(comp)) and bool((comp>0).all() or (comp<0).all()),
                  duplicate_face_copies=d['duplicate_face_copies'])
    t=isolated_tetgen(V,F,VSTATE['flags'],VSTATE['tetgen_timeout'])
    record['tetgen']={k:v for k,v in t.items() if k!='tetgen_stdout'}
    record['still_tet_ready']=bool(t.get('succeeded'))
    return record


def _pair_check(task):
    out=[]
    for a,b in task:
        VA,FA=VSTATE['final'][a];VB,FB=VSTATE['final'][b]
        if not len(FA) or not len(FB):continue
        rec={'a':a,'b':b}
        try:
            IF=cgal.intersect_other(VA,FA,VB,FB,True,False,False,False,2_000_000)[0]
            n=int(len(np.asarray(IF).reshape(-1,2)))
        except BaseException as error:
            rec['error']='intersect_other:'+type(error).__name__;out.append(rec);continue
        if not n:continue
        rec['intersecting_face_pairs']=n
        try:
            r=cgal.mesh_boolean(VA,FA,VB,FB,type_str='intersect')
            IV=np.asarray(r[0],float);IFc=np.asarray(r[1],np.int64)
            rec['residual_overlap_volume_m3']=abs(divergence(IV,IFc)) if len(IFc) else 0.0
            rec['residual_method']='cgal_boolean_intersect'
        except BaseException as error:
            rec['residual_overlap_volume_m3']=None
            rec['residual_method']='cgal_boolean_failed:'+type(error).__name__
        out.append(rec)
    return out


def verify(out,final,before,rows_by,workers,flags,chunk,log,volume_tol,tetgen_timeout):  # noqa: C901
    VSTATE.update(final=final,before=before,flags=flags,tetgen_timeout=tetgen_timeout,
                  role={e:rows_by[e]['role'] for e in final},
                  system={e:rows_by[e]['system'] for e in final})
    ids=sorted(final)
    lo=np.array([final[e][0].min(0) if len(final[e][1]) else np.full(3,np.inf) for e in ids])
    hi=np.array([final[e][0].max(0) if len(final[e][1]) else np.full(3,-np.inf) for e in ids])
    ov=(lo[:,None,:]<=hi[None,:,:]).all(2)&(hi[:,None,:]>=lo[None,:,:]).all(2)
    np.fill_diagonal(ov,False)
    idx=np.argwhere(np.triu(ov))
    pairs=[(ids[int(i)],ids[int(j)]) for i,j in idx]
    log('  verify: %d entities, %d bounding-box-overlapping pairs'%(len(ids),len(pairs)))
    tasks=[pairs[k:k+chunk] for k in range(0,len(pairs),chunk)]
    began=time.monotonic()
    residual,dead_pairs=run_tasks(_pair_check,tasks,workers,log,'pair tasks',True)
    log('  verify: %d touching pairs found in %.1f s, %d pair tasks unrecoverable'
        %(len(residual),time.monotonic()-began,len(dead_pairs)))
    began=time.monotonic()
    entities,dead_entities=run_tasks(_entity_check,ids,workers,log,'entity checks',False)
    entities.sort(key=lambda r:r['entity_id'])
    log('  verify: per-entity checks in %.1f s, %d entities unrecoverable'%(time.monotonic()-began,len(dead_entities)))
    volumetric=[r for r in residual if (r.get('residual_overlap_volume_m3') or 0.0)>volume_tol]
    touching=[r for r in residual if r.get('intersecting_face_pairs') and (r.get('residual_overlap_volume_m3') or 0.0)<=volume_tol]
    report={'entities':len(ids),'bounding_box_overlapping_pairs':len(pairs),
        'coverage':'exhaustive: every bounding-box-overlapping pair of the emitted set tested with CGAL intersect_other, then an exact CGAL boolean intersection volume for every pair that reported a face contact',
        'pairs_with_face_contact':len(residual),
        'pairs_with_residual_overlap_volume':len(volumetric),
        'residual_volume_tolerance_m3':volume_tol,
        'total_residual_overlap_volume_m3':float(sum((r.get('residual_overlap_volume_m3') or 0.0) for r in residual)),
        'pairs_touching_with_zero_volume':len(touching),
        'boolean_failed_pairs':[r for r in residual if r.get('residual_method','').startswith('cgal_boolean_failed')],
        'errored_pairs':[r for r in residual if 'error' in r],
        'pairs_whose_worker_died':[list(x) if isinstance(x,list) else x for x in dead_pairs],
        'entities_whose_worker_died':dead_entities,
        'worst_residual':sorted(volumetric,key=lambda r:-(r['residual_overlap_volume_m3'] or 0))[:40],
        'interpretation':'a face contact with zero boolean volume is the intended outcome: the two surfaces meet on a shared interface and claim no common volume. Only a nonzero boolean intersection volume is an unresolved conflict.'}
    return report,residual,entities


# ---------------------------------------------------------------- conformity

def conformity(final,results,index_to_id,edges):
    wanted={}
    for eid,r in sorted(results.items()):
        origin=r['origin'];V,F=final[eid]
        if not len(F) or not len(origin):continue
        for oi in np.unique(origin[origin!=SELF]):
            wanted.setdefault(index_to_id[int(oi)],[]).append((eid,np.ascontiguousarray(F[origin==oi]),V))
    per=[]
    for owner,claims in sorted(wanted.items()):
        OV,OF=final[owner]
        okeys=np.sort(facet_key(OV,OF))
        for eid,sel,V in claims:
            k=facet_key(V,sel)
            matched=int(len(np.intersect1d(k,okeys,assume_unique=False)))
            per.append({'loser':eid,'owner':owner,'interface_facets':int(len(sel)),
                        'facets_shared_with_owner':matched,
                        'conforming_fraction':matched/len(sel) if len(sel) else None})
    resolved={tuple(sorted((r['a'],r['b']))) for r in edges}
    withiface={tuple(sorted((p['loser'],p['owner']))) for p in per}
    full=[p for p in per if p['interface_facets'] and p['facets_shared_with_owner']==p['interface_facets']]
    partial=[p for p in per if 0<p['facets_shared_with_owner']<p['interface_facets']]
    none=[p for p in per if p['interface_facets'] and p['facets_shared_with_owner']==0]
    facets=sum(p['interface_facets'] for p in per)
    shared=sum(p['facets_shared_with_owner'] for p in per)
    report={'resolved_pairs':len(resolved),
        'pairs_that_produced_an_interface':len(withiface),
        'pairs_with_no_interface_left':len(resolved-withiface),
        'interfaces_fully_conforming':len(full),'interfaces_partially_conforming':len(partial),
        'interfaces_coincident_but_not_conforming':len(none),
        'interface_facets':facets,'interface_facets_shared_with_owner':shared,
        'interface_facet_conforming_fraction':(shared/facets) if facets else None,
        'definition':'an interface facet is genuinely conforming when the exact coordinate triple the loser carries also exists as a facet of the owner, so the two surfaces share nodes there; a facet that is coincident with the owner but differently triangulated is merely non-overlapping',
        'worst':sorted([p for p in per if p['conforming_fraction'] is not None],key=lambda p:(p['conforming_fraction'],-p['interface_facets']))[:20]}
    return report,per


# ---------------------------------------------------------------- self-test

def _cube(lo,hi):
    """Outward-oriented axis-aligned box. The atlas surfaces are outward-oriented, and CGAL reads an
    inward-oriented shell as the complement, so the fixture must match the data."""
    lo=np.asarray(lo,float);hi=np.asarray(hi,float)
    V=np.ascontiguousarray(lo+BOX_CORNERS*(hi-lo));F=np.ascontiguousarray(BOX[:,[0,2,1]].copy())
    assert divergence(V,F)>0
    return V,F


def self_test():
    checks=[]
    def check(name,ok,detail):checks.append({'check':name,'passed':bool(ok),'detail':detail})

    VA,FA=_cube([0,0,0],[1,1,1]);VB,FB=_cube([.75,0,0],[1.75,1,1])
    V,F,J=boolean_minus(VB,FB,VA,FA)
    check('boolean minus removes exactly the analytic slab',abs(abs(divergence(V,F))-.75)<1e-12,abs(divergence(V,F)))
    b,n=edge_counts(F)
    check('boolean minus stays closed and manifold',b==0 and n==0,(b,n))
    V2,F2,_=arrangement_minus(VB,FB,VA,FA)
    check('arrangement fallback agrees with the boolean',abs(abs(divergence(V2,F2))-.75)<1e-12,abs(divergence(V2,F2)))
    check('arrangement fallback stays closed',edge_counts(F2)==(0,0),edge_counts(F2))

    VS,FS=_cube([.4,.4,.4],[.6,.6,.6])
    V3,F3,origin=cut(VA,FA,np.full(len(FA),SELF,np.int64),VS,FS,np.full(len(FS),7,np.int64))[:3]
    check('nested owner leaves a cavity of exactly its own volume',abs(abs(divergence(V3,F3))-(1-.008))<1e-12,abs(divergence(V3,F3)))
    inherited=facet_key(V3,np.ascontiguousarray(F3[origin==7]))
    owned=set(facet_key(VS,FS).tolist())
    check('a fully engulfed owner hands over every one of its facets verbatim',
          len(inherited)==len(FS) and all(k in owned for k in inherited.tolist()),(len(inherited),len(FS)))

    V4,F4,_,_=cut(VS,FS,np.full(len(FS),SELF,np.int64),VA,FA,np.zeros(len(FA),np.int64))
    check('a structure entirely inside a higher-priority one is annihilated, not shrunk',len(F4)==0,int(len(F4)))

    key={'bone':(0,-1.0,'bone'),'muscle':(7,-1.0,'muscle'),'lymph':(10,-1.0,'lymph')}
    edges=[{'a':'bone','b':'muscle'},{'a':'muscle','b':'lymph'},{'a':'bone','b':'lymph'}]
    owners,level,buckets=layers(['bone','muscle','lymph'],edges,key)
    check('the sweep orders a three-way chain by priority',
          [b for b in buckets]==[['bone'],['muscle'],['lymph']] and sorted(owners['lymph'])==['bone','muscle'],
          (buckets,owners))
    check('the total order reproduces the recorded pairwise rule',
          decide({'entity_id':'y','role':'muscle','system':'muscular'},
                 {'entity_id':'x','role':'rigid_bone','system':'skeletal'},{'x':2.0,'y':1.0})[0]=='x'
          and key['bone']<key['muscle']<key['lymph'],None)

    # three mutually overlapping cubes: the resolved volumes must sum to the union, never more.
    parts={'a':_cube([0,0,0],[1,1,1]),'b':_cube([.5,0,0],[1.5,1,1]),'c':_cube([.25,.5,0],[1.25,1.5,1])}
    roles={'a':0,'b':7,'c':10}
    keyed={e:(roles[e],-1.0,e) for e in parts}
    tri=[{'a':'a','b':'b'},{'a':'a','b':'c'},{'a':'b','b':'c'}]
    ow,_,bk=layers(list(parts),tri,keyed)
    fin,fres,fops=resolve(sorted(parts),parts,
                          [{'a':t['a'],'b':t['b'],'overlap_volume_m3':0.0} for t in tri],
                          ow,{e:i for i,e in enumerate(sorted(parts))},1,lambda m:None,8<<30,120.0,300.0)
    chain={}
    for e in sorted(parts,key=lambda x:keyed[x]):
        V,F=parts[e];o=np.full(len(F),SELF,np.int64)
        for owner in sorted(ow[e],key=lambda x:keyed[x]):
            V,F,o,_=cut(V,F,o,*chain[owner],np.zeros(len(chain[owner][1]),np.int64))
        chain[e]=(V,F)
    check('cutting against original owners matches cutting against resolved owners',
          all(abs(abs(divergence(*fin[e]))-abs(divergence(*chain[e])))<1e-12 for e in fin),
          {e:(abs(divergence(*fin[e])),abs(divergence(*chain[e]))) for e in fin})
    total=sum(abs(divergence(*fin[e])) for e in fin)
    union=cgal.mesh_boolean(*parts['a'],*parts['b'],type_str='union')
    union=cgal.mesh_boolean(np.asarray(union[0],float),np.asarray(union[1],np.int64),*parts['c'],type_str='union')
    uvol=abs(divergence(np.asarray(union[0],float),np.asarray(union[1],np.int64)))
    check('a three-way overlap is conceded once, not twice',abs(total-uvol)<1e-12,(total,uvol))
    resid=0.0
    for x,y in (('a','b'),('a','c'),('b','c')):
        if not len(fin[x][1]) or not len(fin[y][1]):continue
        r=cgal.mesh_boolean(*fin[x],*fin[y],type_str='intersect')
        resid+=abs(divergence(np.asarray(r[0],float),np.asarray(r[1],np.int64))) if len(r[1]) else 0.0
    check('the resolved three-way set shares no volume',resid<1e-15,resid)

    VD1,FD1=_cube([0,0,0],[1,1,1]);VD2,FD2=_cube([0,0,0],[1,1,1])
    rows=[{'entity_id':'later','role':'muscle'},{'entity_id':'earlier','role':'muscle'}]
    dig={'later':entity_digest(VD1,FD1),'earlier':entity_digest(VD2,FD2)}
    groups,dropped=duplicate_groups(rows,dig)
    check('byte-identical ids are grouped and the later one dropped',groups==[['earlier','later']] and dropped==['later'],(groups,dropped))
    check('a moved vertex breaks the digest',entity_digest(VD1+1e-12,FD1)!=entity_digest(VD1,FD1),None)

    t=tetrahedralize(*_cube([0,0,0],[1,1,1]),'pYq1.414')
    check('the imported volume-gated TetGen harness passes a unit cube',t['succeeded'] and t['tets']>0,
          {'tets':t.get('tets'),'gate':t.get('volume_gate'),'error':t.get('tet_volume_vs_surface_relative_error')})
    ti=isolated_tetgen(*_cube([0,0,0],[1,1,1]),'pYq1.414',120.0)
    check('the time-capped TetGen wrapper agrees with the harness',
          ti.get('succeeded') and ti.get('tets')==t['tets'],{'tets':ti.get('tets'),'harness':t['tets']})
    tz=isolated_tetgen(*_cube([0,0,0],[1,1,1]),'pYq1.414',0.0)
    check('a TetGen run over the time cap is killed and recorded',
          tz.get('exception')=='TimeCap' and not tz.get('succeeded'),tz.get('exception'))

    PV,PF,_=compact(V3,np.ascontiguousarray(F3[origin==7]))
    rec,mesh=_imprint_worker({'owner':'o','V':VS,'F':FS,'patches':[(PV,PF)],
                              'origin':np.full(len(FS),SELF,np.int64)})
    check('imprinting an owner with the patch it already owns leaves its volume exact',
          rec['adopted'] and abs(rec['relative_volume_error'])<1e-12,rec)

    # a partial cut: the owner is refined so the loser's inherited facets become shared facets
    VP,FP=_cube([.4,.4,.4],[.6,.6,1.4])
    VQ,FQ,oq,_=cut(VA,FA,np.full(len(FA),SELF,np.int64),VP,FP,np.full(len(FP),3,np.int64))
    QP,QF,_=compact(VQ,np.ascontiguousarray(FQ[oq==3]))
    shared_before=len(set(facet_key(QP,QF).tolist())&set(facet_key(VP,FP).tolist()))
    rec2,mesh2=_imprint_worker({'owner':'p','V':VP,'F':FP,'patches':[(QP,QF)],
                                'origin':np.full(len(FP),SELF,np.int64)})
    after_mesh=(mesh2[0],mesh2[1]) if mesh2 else (VP,FP)
    shared_after=len(set(facet_key(QP,QF).tolist())&set(facet_key(*after_mesh).tolist())) if rec2['adopted'] else shared_before
    check('imprinting turns a partial cut into shared facets',
          rec2['adopted'] and shared_after>shared_before and shared_after==len(QF),
          {'before':shared_before,'after':shared_after,'interface_facets':int(len(QF))})

    VH,FH=_cube([0,0,0],[1,1,1])
    holed=np.ascontiguousarray(FH[2:])
    fixed,hrec=close_holes(VH,holed)
    check('a hole is filled back to the exact original volume',
          hrec['filled'] and abs(abs(divergence(*fixed))-1)<1e-12,hrec)
    check('an already closed surface is left alone',close_holes(VH,FH)[1]['reason']=='already closed',None)

    big=_cube([0,0,0],[1,1,1]);small=_cube([.4,.4,.4],[.6,.6,.6])
    cm={'big':big,'small':small};cids=['big','small']
    import igl.copyleft.cgal as _c
    crossing=int(len(np.asarray(_c.intersect_other(big[0],big[1],small[0],small[1],True,False,False,False,2_000_000)[0]).reshape(-1,2)))
    cfound,crep=census(cm,cids,1,8,lambda m:None)
    check('the census sees a nested pair that crosses no face',
          crossing==0 and crep['nested_pairs']==1 and abs(cfound[0]['overlap_volume_m3']-.008)<1e-12,
          {'crossing_face_pairs':crossing,'nested':crep['nested_pairs'],
           'overlap':cfound[0]['overlap_volume_m3'] if cfound else None})
    far={'x':_cube([0,0,0],[1,1,1]),'y':_cube([0,0,0],[1,1,1])}
    far['y']=(far['y'][0]+np.array([0.,0,0.999]),far['y'][1])
    _,frep=census({'x':_cube([0,0,0],[1,1,1]),'y':_cube([5,5,5],[6,6,6])},['x','y'],1,8,lambda m:None)
    check('the census does not invent a conflict between disjoint surfaces',
          frep['conflicting_pairs']==0,frep['conflicting_pairs'])

    key_a=facet_key(*_cube([0,0,0],[1,1,1]));key_b=facet_key(*_cube([0,0,0],[1,1,1]))
    iv,iff,io,im,isec=isolated_cut(VB,FB,np.full(len(FB),SELF,np.int64),VA,FA,
                                   np.zeros(len(FA),np.int64),8<<30,60)
    check('the isolated cut returns the same answer as the in-process one',
          iv is not None and abs(abs(divergence(iv,iff))-.75)<1e-12,(im,None if iv is None else abs(divergence(iv,iff))))
    _,_,_,tm,_=isolated_cut(VB,FB,np.full(len(FB),SELF,np.int64),VA,FA,
                            np.zeros(len(FA),np.int64),8<<30,0.0)
    check('an overrunning cut is refused, not left to run',tm.startswith('refused_over_time_cap'),tm)
    _,_,_,mm,_=isolated_cut(VB,FB,np.full(len(FB),SELF,np.int64),VA,FA,
                            np.zeros(len(FA),np.int64),1<<26,60)
    check('a cut over the address-space cap is refused, not left to swap',mm.startswith('refused_child_died'),mm)

    check('the facet key is winding-invariant and exact',set(key_a.tolist())==set(key_b.tolist())
          and len(set(key_a.tolist()))==12,len(set(key_a.tolist())))

    passed=all(c['passed'] for c in checks)
    print(json.dumps({'self_test':'ihm.conflict-free-surface-atlas','passed':passed,'checks':checks},indent=2))
    return passed


# --------------------------------------------------------------------- run

def run(args):
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    began=time.monotonic()
    def log(message):print(message,flush=True)

    log('loading %s inputs'%'repaired');t0=time.monotonic()
    rows,meshes,volume,digest=load_inputs()
    rows_by={r['entity_id']:r for r in rows}
    log('  %d entities, %d facets in %.1f s'%(len(rows),sum(len(F) for _,F in meshes.values()),time.monotonic()-t0))

    groups,dropped=duplicate_groups(rows,digest)
    recorded=json.loads((REPAIR/'coincidence.json').read_text())
    dup={'groups':groups,'dropped':dropped,
         'recomputed_independently':True,
         'agrees_with_recorded_coincidence_phase':sorted(dropped)==sorted(recorded['identical_geometry_entities_dropped']),
         'recorded_dropped':recorded['identical_geometry_entities_dropped'],
         'rule':'two ids carrying the same welded coordinate set and the same sorted facet triples are one surface authored twice under synonymous BodyParts3D ids. The lexicographically later id is dropped and its conflicts are rewired to the kept id rather than differencing a surface against itself.'}
    recorded_pairs=sum(1 for l in (REPAIR/'conflict-pairs.jsonl').read_text().splitlines()
                       if json.loads(l).get('intersecting_face_pairs',0)>0)
    keep_of={}
    for g in groups:
        for e in g[1:]:keep_of[e]=g[0]
    active=[r['entity_id'] for r in rows if r['entity_id'] not in set(dropped)]
    edges,removed,rewired=load_edges(set(dropped),keep_of)
    dup['self_pair_edges_removed']=removed
    dup['edges_rewired_to_the_kept_id']=rewired
    dup['recorded_conflicting_pairs']=recorded_pairs
    dup['conflicting_pairs_after_duplicate_removal']=len(edges)
    write_json(out/'duplicates.json',dup)
    log('duplicates: %d groups, %d ids dropped, %d conflicting pairs remain'%(len(groups),len(dropped),len(edges)))

    repairs=[]
    for e in active:
        V,F=meshes[e]
        if edge_counts(F)[0]==0:continue
        fixed,record=close_holes(V,F)
        record['entity_id']=e;record['role']=rows_by[e]['role'];record['name']=rows_by[e]['name']
        record['boundary_edges_before']=edge_counts(F)[0]
        if fixed is not None:meshes[e]=fixed;volume[e]=abs(divergence(*fixed))
        repairs.append(record)
    write_json(out/'hole-fill.json',{'schema':'ihm.conflict-free-atlas-hole-fill.v1',
        'open_input_surfaces':len(repairs),'filled':sum(1 for r in repairs if r['filled']),
        'records':repairs,
        'why':'CGAL mesh_boolean requires an operand with a well defined interior. Six of the 2403 repaired surfaces carry boundary edges, so every conflict they take part in would otherwise fall back to the arrangement path or be skipped. The hole is fan-triangulated from its own centroid, and the fill is adopted only if the shell becomes closed, manifold and self-intersection-free.'})
    log('hole fill: %d open input surfaces, %d closed'%(len(repairs),sum(1 for r in repairs if r['filled'])))

    t0=time.monotonic()
    fresh,census_report=census({e:meshes[e] for e in active},sorted(active),args.workers,args.chunk,log)
    recorded_set={tuple(sorted((x['a'],x['b']))) for x in edges}
    fresh_set={tuple(sorted((r['a'],r['b']))) for r in fresh}
    census_report['agreement_with_the_recorded_census']={
        'recorded_conflicting_pairs_after_duplicate_removal':len(edges),
        'pairs_found_again':len(recorded_set&fresh_set),
        'recorded_pairs_not_found_now':sorted(recorded_set-fresh_set)[:50],
        'recorded_pairs_not_found_now_count':len(recorded_set-fresh_set),
        'pairs_the_recorded_census_missed':len(fresh_set-recorded_set),
        'missed_because_nested':sum(1 for r in fresh if r.get('intersecting_face_pairs',0)==0
                                    and tuple(sorted((r['a'],r['b']))) not in recorded_set),
        'note':'a nested pair crosses no face, so the recorded census could not see it; the hole fill also closes five lung lobes that then enclose volume they did not bound before'}
    write_json(out/'conflict-census.json',{'schema':'ihm.conflict-free-atlas-census.v1',**census_report})
    with (out/'conflict-census.jsonl').open('w') as handle:
        for r in sorted(fresh,key=lambda r:-(r.get('overlap_volume_m3') or 0.0)):
            handle.write(json.dumps(r,allow_nan=False)+'\n')
    log('census: %d conflicting pairs (%d crossing, %d nested) in %.1f s; recorded census missed %d'
        %(census_report['conflicting_pairs'],census_report['crossing_pairs'],census_report['nested_pairs'],
          time.monotonic()-t0,census_report['agreement_with_the_recorded_census']['pairs_the_recorded_census_missed']))
    edges=[{'a':r['a'],'b':r['b'],'intersecting_face_pairs':r.get('intersecting_face_pairs',0),
            'overlap_volume_m3':r.get('overlap_volume_m3'),'overlap_method':r.get('overlap_method'),
            'nested':r.get('intersecting_face_pairs',0)==0} for r in fresh]

    key=order_keys(rows,volume,set(active))
    owners,level,buckets=layers(active,edges,key)
    index={e:i for i,e in enumerate(sorted(active))};index_to_id={i:e for e,i in index.items()}
    log('sweep: %d entities in %d priority levels (deepest chain %d)'%(len(active),len(buckets),max(level.values())))

    resume=args.phase=='verify'
    if resume:
        ops=[json.loads(l) for l in (out/'pair-operations.jsonl').read_text().splitlines()]
        imprints=json.loads((out/'imprint.json').read_text())['records'] if (out/'imprint.json').exists() else []
        conform={k:v for k,v in json.loads((out/'conformity.json').read_text()).items() if k!='schema'}
        written=[]
        for path in sorted((out/'geometry').glob('*.json.gz')):
            written.append({'entity_id':path.stem.replace('.json',''),
                            'output_path':path.relative_to(out).as_posix(),'output_sha256':sha(path)})
        log('resume: %d emitted surfaces, %d recorded pair operations'%(len(written),len(ops)))
        edges=[{'a':r['a'],'b':r['b'],'intersecting_face_pairs':r.get('intersecting_face_pairs',0),
                'overlap_volume_m3':r.get('overlap_volume_m3')}
               for r in (json.loads(l) for l in (out/'conflict-census.jsonl').read_text().splitlines())]
        key=order_keys(rows,volume,set(active))
        owners,level,buckets=layers(active,edges,key)
        return finish(args,out,began,log,rows,rows_by,meshes,volume,active,edges,owners,level,buckets,
                      dropped,groups,dup,repairs,recorded_pairs,ops,imprints,conform,written,None)

    t0=time.monotonic()
    final,results,ops=resolve(active,{e:meshes[e] for e in active},edges,owners,index,args.workers,log,
                              int(args.memory_cap_gb*(1<<30)),args.time_cap_s,args.entity_budget_s)
    log('resolve: %d pair operations in %.1f s'%(len(ops),time.monotonic()-t0))
    with (out/'pair-operations.jsonl').open('w') as handle:
        for r in sorted(ops,key=lambda r:(r['yields'],r['owner'])):handle.write(json.dumps(r,allow_nan=False)+'\n')

    imprints=[]
    if not args.no_imprint:
        t0=time.monotonic()
        final,imprints=imprint(final,results,index_to_id,args.workers,args.imprint_cap,log)
        log('imprint: %d owners processed, %d adopted, in %.1f s'%(len(imprints),sum(1 for r in imprints if r.get('adopted')),time.monotonic()-t0))
        write_json(out/'imprint.json',{'schema':'ihm.conflict-free-atlas-imprint.v1',
            'owners_with_a_cut_patch':len(imprints),
            'adopted':sum(1 for r in imprints if r.get('adopted')),
            'rejected':[r for r in imprints if not r.get('adopted')][:200],
            'rejected_count':sum(1 for r in imprints if not r.get('adopted')),
            'guard':'an imprint is adopted only if the owner stays closed, manifold, self-intersection-free and its volume is unchanged to 1e-9 relative; otherwise the owner keeps its unrefined triangulation and the rejection is recorded',
            'records':imprints})

    conform,interfaces=conformity(final,results,index_to_id,edges)
    write_json(out/'conformity.json',{'schema':'ihm.conflict-free-atlas-conformity.v1',**conform})
    with (out/'interfaces.jsonl').open('w') as handle:
        for r in sorted(interfaces,key=lambda r:(r['loser'],r['owner'])):handle.write(json.dumps(r,allow_nan=False)+'\n')
    log('conformity: %d/%d interface facets share nodes with their owner'%(conform['interface_facets_shared_with_owner'],conform['interface_facets']))

    written=emit(out,final,rows_by,results,index_to_id)
    log('emit: %d surfaces written'%len(written))
    return finish(args,out,began,log,rows,rows_by,meshes,volume,active,edges,owners,level,buckets,
                  dropped,groups,dup,repairs,recorded_pairs,ops,imprints,conform,written,final)


def finish(args,out,began,log,rows,rows_by,meshes,volume,active,edges,owners,level,buckets,
           dropped,groups,dup,repairs,recorded_pairs,ops,imprints,conform,written,final):
    """Verification and receipts, from the emitted files rather than from memory."""
    reloaded={w['entity_id']:load_mesh(out/w['output_path']) for w in written}
    drift=[] if final is None else [w['entity_id'] for w in written
           if not np.array_equal(reloaded[w['entity_id']][1],final[w['entity_id']][1])
           or not np.array_equal(reloaded[w['entity_id']][0],final[w['entity_id']][0])]
    log('emit: %d surfaces reloaded, %d differ from memory'%(len(reloaded),len(drift)))
    for w in written:
        V,F=reloaded[w['entity_id']];w['vertices']=int(len(V));w['faces']=int(len(F))
    before={e:volume[e] for e in reloaded}
    t0=time.monotonic()
    residual_report,residual,entities=verify(out,reloaded,before,rows_by,args.workers,args.flags,args.chunk,log,
                                            args.volume_tolerance,args.tetgen_timeout_s)
    residual_report['round_trip_mismatches']=drift
    residual_report['round_trip_note']=('the emitted files were reloaded and compared against the in-memory result'
                                        if final is not None else 'resumed run: verification reads the emitted files only')
    log('verify: %.1f s'%(time.monotonic()-t0))
    write_json(out/'residual-conflicts.json',{'schema':'ihm.conflict-free-atlas-residual.v1',**residual_report})
    with (out/'residual-pairs.jsonl').open('w') as handle:
        for r in sorted(residual,key=lambda r:-(r.get('residual_overlap_volume_m3') or 0.0)):
            handle.write(json.dumps(r,allow_nan=False)+'\n')

    degree={}
    for x in edges:
        degree[x['a']]=degree.get(x['a'],0)+1;degree[x['b']]=degree.get(x['b'],0)+1
    ledger={}
    for line in (REPAIR/'entity-operations.jsonl').read_text().splitlines():
        r=json.loads(line);ledger[r['entity_id']]=r
    by_written={w['entity_id']:w for w in written}
    log_rows=[]
    for r in entities:
        e=r['entity_id']
        r=dict(r)
        r['name']=rows_by[e]['name']
        r['owners']=sorted(owners.get(e,[]))
        r['conflicts']=degree.get(e,0)
        r['yields_to']=len(owners.get(e,[]))
        r['priority_level']=level.get(e)
        r['output_path']=by_written[e]['output_path'];r['output_sha256']=by_written[e]['output_sha256']
        led=ledger.get(e)
        r['ledger_conceded_volume_m3']=led['volume_conceded_m3'] if led else None
        r['ledger_minus_measured_m3']=(led['volume_conceded_m3']-r['volume_conceded_m3']) if led else None
        log_rows.append(r)
    with (out/'entities.jsonl').open('w') as handle:
        for r in sorted(log_rows,key=lambda r:r['entity_id']):handle.write(json.dumps(r,allow_nan=False)+'\n')

    annihilated=[r for r in log_rows if r['annihilated']]
    heavy=[r for r in log_rows if r['fraction_conceded'] is not None and r['fraction_conceded']>.5]
    regressions=[r for r in log_rows if not r['still_tet_ready'] and not r['annihilated']]
    ledger_sum=sum(v['volume_conceded_m3'] for v in ledger.values())
    measured_sum=sum(r['volume_conceded_m3'] for r in log_rows)
    write_json(out/'annihilation.json',{'schema':'ihm.conflict-free-atlas-annihilation.v1',
        'entities_reduced_to_nothing':[{k:r[k] for k in ('entity_id','name','role','system','volume_before_m3','owners')} for r in annihilated],
        'entities_reduced_to_nothing_count':len(annihilated),
        'entities_conceding_over_half_their_volume':sorted(
            [{k:r[k] for k in ('entity_id','name','role','system','volume_before_m3','volume_after_m3','volume_conceded_m3','fraction_conceded','owners')} for r in heavy],
            key=lambda r:-r['fraction_conceded']),
        'entities_conceding_over_half_their_volume_count':len(heavy),
        'measured_total_conceded_volume_m3':measured_sum,
        'ledger_total_conceded_volume_m3':ledger_sum,
        'double_counting_note':'the ledger sums the conceded volume of every pair independently, so a volume disputed by three structures is counted more than once and a per-entity fraction can exceed 1. The measured figure here is volume before minus volume after on the emitted surface, so each cubic metre is conceded exactly once.',
        'finding':'an entity that concedes most of its volume is a source conflict surfaced by the rule, not a successful resolution'})

    summary={'schema':'ihm.conflict-free-surface-atlas.v1',
        'inputs':{'entities_in':len(rows),'duplicate_ids_dropped':len(dropped),'entities_out':len(written),
                  'conflicting_pairs_in':recorded_pairs,'conflicting_pairs_resolved':len(edges),
                  'facets_in':int(sum(len(F) for _,F in meshes.values())),
                  'facets_out':int(sum(w['faces'] for w in written))},
        'ordering':{'rule':'strict total order (declared role rank, then larger volume, then entity id); the lower key owns the disputed volume',
                    'priority_table':{k:v[0] for k,v in PRIORITY.items()},
                    'levels':len(buckets),'deepest_chain':max(level.values()) if level else 0,
                    'why':'a total order makes ownership globally acyclic; each structure is cut once against the union of the higher-priority neighbours it yields to'},
        'methods':{m:sum(1 for o in ops if o['method']==m) for m in sorted({o['method'] for o in ops})},
        'volume':{'before_m3':float(sum(volume[e] for e in by_written)),
                  'after_m3':float(sum(r['volume_after_m3'] for r in log_rows)),
                  'conceded_m3':measured_sum,
                  'ledger_pairwise_sum_m3':ledger_sum,
                  'recorded_total_pairwise_overlap_m3':0.0031367032975771186},
        'duplicates':{'groups':len(groups),'dropped':dropped,'agrees_with_recorded':dup['agrees_with_recorded_coincidence_phase']},
        'hole_fill':{'open_input_surfaces':len(repairs),'filled':sum(1 for r in repairs if r['filled']),
                     'entities':[r['entity_id'] for r in repairs]},
        'census':json.loads((out/'conflict-census.json').read_text()) if (out/'conflict-census.json').exists() else None,
        'residual':residual_report,
        'conformity':conform,
        'imprint':{'owners':len(imprints),'adopted':sum(1 for r in imprints if r.get('adopted'))} if imprints else None,
        'per_entity':{'closed':sum(1 for r in log_rows if r['closed']),
                      'manifold':sum(1 for r in log_rows if r['manifold']),
                      'self_intersection_free':sum(1 for r in log_rows if r['self_intersection_free']),
                      'orientation_consistent':sum(1 for r in log_rows if r['orientation_consistent']),
                      'tetrahedralized':sum(1 for r in log_rows if r['still_tet_ready']),
                      'tetgen_flags':args.flags,'tetgen_volume_gate':VOLUME_GATE,
                      'tetgen_timeout_s':args.tetgen_timeout_s,
                      'tetgen_timed_out':sum(1 for r in log_rows if r['tetgen'].get('exception')=='TimeCap'),
                      'total_tets':int(sum(r['tetgen'].get('tets',0) or 0 for r in log_rows)),
                      'annihilated':len(annihilated),
                      'tetgen_regressions':[{k:r[k] for k in ('entity_id','name','role','faces','volume_after_m3','fraction_conceded')}|{'tetgen':{k:r['tetgen'].get(k) for k in ('status','tets','exception','message','tet_volume_vs_surface_relative_error')}} for r in regressions],
                      'tetgen_regression_count':len(regressions)},
        'annihilation':{'reduced_to_nothing':len(annihilated),'conceding_over_half':len(heavy)},
        'wall_seconds':time.monotonic()-began,'python':sys.version}
    write_json(out/'summary.json',summary)

    inputs={'data/derived/canonical/anatomy.json':sha(ROOT/'data/derived/canonical/anatomy.json'),
            'data/derived/muscle-tet-ready-v1/manifest.json':sha(ROOT/'data/derived/muscle-tet-ready-v1/manifest.json'),
            'data/derived/entity-tet-ready-v1/manifest.json':sha(ROOT/'data/derived/entity-tet-ready-v1/manifest.json'),
            'data/derived/cross-structure-repair-v1/manifest.json':sha(REPAIR/'manifest.json'),
            'data/derived/cross-structure-repair-v1/conflict-pairs.jsonl':sha(REPAIR/'conflict-pairs.jsonl'),
            'data/derived/cross-structure-repair-v1/ownership-ledger.jsonl':sha(REPAIR/'ownership-ledger.jsonl'),
            'data/derived/cross-structure-repair-v1/coincidence.json':sha(REPAIR/'coincidence.json'),
            'scripts/build_conflict_free_surface_atlas.py':sha(__file__),
            'scripts/build_cross_structure_conflict_repair.py':sha(ROOT/'scripts/build_cross_structure_conflict_repair.py'),
            'scripts/verify_muscle_tet_ready_surfaces.py':sha(ROOT/'scripts/verify_muscle_tet_ready_surfaces.py')}
    artifacts={p.relative_to(out).as_posix():sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    write_json(out/'manifest.json',{'schema':'ihm.conflict-free-surface-atlas-manifest.v1',
        'inputs_sha256':inputs,
        'per_entity_input_sha256':{r['entity_id']:r['sha256'] for r in rows},
        'artifacts_sha256':artifacts,
        'canonical_assets_modified':False,
        'native_cgal_sha256':sha(Path(cgal.pyigl_copyleft_cgal.__file__))})
    print(json.dumps({k:v for k,v in summary.items() if k not in ('residual','conformity')},indent=2)[:8000])
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--out',default=str(ROOT/'data/derived/conflict-free-atlas-v1'))
    p.add_argument('--phase',default='all',choices=['all','verify'])
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--workers',type=int,default=10)
    p.add_argument('--chunk',type=int,default=64)
    p.add_argument('--flags',default='pYq1.414')
    p.add_argument('--no-imprint',action='store_true')
    p.add_argument('--imprint-cap',type=int,default=400_000)
    p.add_argument('--volume-tolerance',type=float,default=1e-15)
    p.add_argument('--memory-cap-gb',type=float,default=8.0)
    p.add_argument('--time-cap-s',type=float,default=300.0)
    p.add_argument('--entity-budget-s',type=float,default=900.0)
    p.add_argument('--tetgen-timeout-s',type=float,default=600.0)
    a=p.parse_args()
    if a.self_test:sys.exit(0 if self_test() else 1)
    run(a)
