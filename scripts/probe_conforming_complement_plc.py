"""How far the canonical structure surfaces are from one global piecewise-linear complex.

A single conforming tetrahedralization of the body plus its interstitial complement needs the
whole surface set to be one valid PLC: every facet closed and consistently oriented, no facet
self-intersecting, no two facets from different structures crossing, and no coincident facets.
This measures each of those conditions on the retained surfaces -- exhaustively where that is
cheap, on an explicit random sample where it is not -- and then actually runs TetGen on one
anatomically adjacent cluster inside a padded box so the complement between the structures has
to be meshed too.

Requires the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/probe_conforming_complement_plc.py

Not measured here: any global tetrahedralization is attempted, element quality of the cluster
mesh beyond volume conformity, and whether a conforming mesh would solve better than voxels.
"""
from pathlib import Path
import argparse,ctypes,gzip,hashlib,importlib.util,json,os,sys,tempfile,time
import numpy as np
import igl
import igl.copyleft.cgal as cgal
from igl.copyleft import tetgen

LIBC=ctypes.CDLL(None)
ROOT=Path(__file__).resolve().parents[1]
MUSCLE=ROOT/'data/derived/muscle-tet-ready-v1'
CLUSTER=['body-bp3d-FJ3365','body-bp3d-FJ1433','body-bp3d-FJ1442','body-bp3d-FJ1443','body-bp3d-FJ1444','body-bp3d-FJ1434']
BOX=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,5,1],[0,4,5],[1,6,2],[1,5,6],[2,7,3],[2,6,7],[3,4,0],[3,7,4]])
BOX_CORNERS=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m
MD=load_module(ROOT/'ihm/assembly/material_domains.py','material_domains')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def mesh(path):
    raw=Path(path).read_bytes();d=json.loads(gzip.decompress(raw) if str(path).endswith('.gz') else raw)
    return (np.ascontiguousarray(np.asarray(d.get('positions',d.get('vertices')),float).reshape(-1,3)),
            np.ascontiguousarray(np.asarray(d.get('indices',d.get('faces')),np.int64).reshape(-1,3)))

def divergence(V,F):
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]]
    return float(np.einsum('ij,ij->i',A,np.cross(B,C)).sum()/6)

def sources():
    anatomy=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text());frame=anatomy['frame']['id']
    repaired={json.loads(l)['entity_id']:json.loads(l) for l in (MUSCLE/'entities.jsonl').read_text().splitlines()}
    rows=[]
    for e in anatomy['entities']:
        ref=e['reference_geometry']
        if ref.get('representation')!='triangular_surface' or ref.get('units')!='m' or ref.get('frame')!=frame:continue
        path=MUSCLE/repaired[e['id']]['output_path'] if e['id'] in repaired else ROOT/ref['path']
        rows.append({'entity_id':e['id'],'name':e['name'],'system':e['system'],'path':path,
                     'lo':np.asarray(e['bounds_m']['min']),'hi':np.asarray(e['bounds_m']['max'])})
    return rows

def coincidence(rows):
    """Exact coordinate coincidence across different entities: shared vertices and duplicated facets."""
    verts=[];faces=[];fent=[];vent=[];offset=0
    for i,r in enumerate(rows):
        V,F=mesh(r['path']);verts.append(V);faces.append(F+offset)
        fent.append(np.full(len(F),i));vent.append(np.full(len(V),i));offset+=len(V)
    V=np.concatenate(verts);F=np.concatenate(faces)
    fent=np.concatenate(fent);vent=np.concatenate(vent)
    unique,inverse=np.unique(V,axis=0,return_inverse=True);inverse=np.ravel(inverse)
    vp=np.unique(np.stack((inverse,vent),axis=1),axis=0)
    _,per_coord=np.unique(vp[:,0],return_counts=True)
    tri=np.sort(inverse[F],axis=1)
    _,uinv,ucount=np.unique(tri,axis=0,return_inverse=True,return_counts=True);uinv=np.ravel(uinv)
    duplicated=ucount[uinv]>1
    cross=0
    if duplicated.any():
        fp=np.unique(np.stack((uinv[duplicated],fent[duplicated]),axis=1),axis=0)
        _,per_facet=np.unique(fp[:,0],return_counts=True);cross=int((per_facet>1).sum())
    return {'entities':len(rows),'total_vertices':int(len(V)),'distinct_coordinates':int(len(unique)),
        'coordinates_used_by_more_than_one_entity':int((per_coord>1).sum()),
        'total_facets':int(len(F)),'distinct_facet_vertex_triples':int(len(ucount)),
        'facets_duplicated_anywhere':int(duplicated.sum()),
        'facet_triples_shared_by_more_than_one_entity':cross,
        'basis':'exact float64 coordinate equality; near-coincident but unequal facets are not counted'}

def pair_audit(rows,sample,seed,max_faces):
    lo=np.stack([r['lo'] for r in rows]);hi=np.stack([r['hi'] for r in rows])
    overlap=(lo[:,None,:]<=hi[None,:,:]).all(2)&(hi[:,None,:]>=lo[None,:,:]).all(2)
    np.fill_diagonal(overlap,False)
    pairs=np.argwhere(np.triu(overlap))
    rng=np.random.default_rng(seed)
    cache={}
    def get(i):
        if i not in cache:cache[i]=mesh(rows[i]['path'])
        return cache[i]
    tested=[];intersecting=0;faces_pairs=0;began=time.perf_counter()
    order=rng.permutation(len(pairs))
    for k in order:
        if len(tested)>=sample:break
        i,j=int(pairs[k][0]),int(pairs[k][1])
        VA,FA=get(i);VB,FB=get(j)
        if len(FA)+len(FB)>max_faces:continue
        IF=cgal.intersect_other(VA,FA,VB,FB,True,False,False,False,1000)[0]
        n=int(len(IF));intersecting+=n>0;faces_pairs+=n
        tested.append({'a':rows[i]['entity_id'],'b':rows[j]['entity_id'],'a_system':rows[i]['system'],
                       'b_system':rows[j]['system'],'intersecting_face_pairs':n})
    return {'entities':len(rows),'bounding_box_overlapping_pairs':int(len(pairs)),
        'all_pairs':int(len(rows)*(len(rows)-1)//2),
        'exactly_tested_pairs':len(tested),'pairs_that_actually_intersect':int(intersecting),
        'fraction_of_tested_pairs_intersecting':float(intersecting/max(len(tested),1)),
        'total_intersecting_face_pairs_found':int(faces_pairs),
        'projected_intersecting_pairs_over_all_bbox_overlaps':float(len(pairs)*intersecting/max(len(tested),1)),
        'seconds':time.perf_counter()-began,'max_combined_faces_tested':max_faces,
        'worst':sorted(tested,key=lambda t:-t['intersecting_face_pairs'])[:12],
        'basis':'exact CGAL intersect_other on a uniform random sample of bounding-box-overlapping pairs; pairs above the face cap are skipped, so the projection is a lower bound if large meshes intersect more'}

def tetrahedralize(V,F,flags):
    """Fork so a TetGen abort is recorded rather than killing the probe; native output captured."""
    handle,name=tempfile.mkstemp(suffix='.tetgen');os.close(handle)
    store=name+'.npz';reader,writer=os.pipe();began=time.monotonic();outcome={'flags':flags}
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
    os.close(reader);_,wait_status=os.waitpid(pid,0)
    if chunks:outcome.update(json.loads(b''.join(chunks).decode()))
    if os.WIFSIGNALED(wait_status):
        outcome.update(status=None,exception='ProcessAborted',
                       message='TetGen terminated the process with signal %d'%os.WTERMSIG(wait_status))
    elif not chunks:outcome.update(status=None,exception='NoResult',message='TetGen child exited without a result')
    text=Path(name).read_text(errors='replace').strip().splitlines();Path(name).unlink(missing_ok=True)
    outcome['tetgen_stdout_lines']=len(text)
    outcome['tetgen_warning_kinds']=dict(sorted({l.strip():text.count(l) for l in set(text) if l.startswith('Warning')}.items(),key=lambda kv:-kv[1])[:6])
    outcome['tetgen_stdout_tail']=text[-12:]
    outcome['succeeded']=outcome.get('status')==0 and outcome.get('tets',0)>0
    outcome['seconds']=time.monotonic()-began
    if Path(store).exists():
        d=np.load(store);outcome['_TV']=d['TV'];outcome['_TT']=d['TT'];Path(store).unlink()
    return outcome

def enclose(parts,pad,flags):
    """Mesh a padded box containing the given surfaces, so the interstitial complement is meshed too."""
    lo=np.min([p['V'].min(0) for p in parts],axis=0)-pad;hi=np.max([p['V'].max(0) for p in parts],axis=0)+pad
    box=lo+BOX_CORNERS*(hi-lo)
    V=[box];F=[BOX];offset=8
    for p in parts:V.append(p['V']);F.append(p['F']+offset);offset+=len(p['V'])
    V=np.ascontiguousarray(np.concatenate(V));F=np.ascontiguousarray(np.concatenate(F).astype(np.int64))
    trial=tetrahedralize(V,F,flags)
    report={'members':[p['entity_id'] for p in parts],'box_volume_m3':float(np.prod(hi-lo)),
        'plc_input':{'vertices':int(len(V)),'facets':int(len(F))},
        'tetgen':{k:v for k,v in trial.items() if not k.startswith('_')}}
    if trial.get('succeeded'):
        TV=trial['_TV'];TT=trial['_TT']
        d=np.linalg.det(np.swapaxes(TV[TT[:,1:]]-TV[TT[:,0,None]],1,2))/6
        centroid=np.ascontiguousarray(TV[TT].mean(1));labels=np.full(len(TT),-1)
        conformity=[]
        for i,p in enumerate(parts):
            w=igl.fast_winding_number(np.ascontiguousarray(p['V']),np.ascontiguousarray(p['F']),centroid)
            inside=np.abs(w)>.5;labels[inside]=i;volume=float(np.abs(d[inside]).sum())
            conformity.append({'entity_id':p['entity_id'],'tets':int(inside.sum()),'tet_volume_m3':volume,
                'surface_volume_m3':p['surface_volume_m3'],
                'relative_volume_error':float(volume/p['surface_volume_m3']-1)})
        complement=labels<0
        report['conforming_mesh']={'tet_vertices':int(len(TV)),'tets':int(len(TT)),
            'all_positive_volume':bool((d>0).all()),'min_tet_volume_m3':float(d.min()),
            'total_volume_m3':float(np.abs(d).sum()),
            'box_volume_closure_relative_error':float(np.abs(d).sum()/report['box_volume_m3']-1),
            'complement_tets':int(complement.sum()),'complement_volume_m3':float(np.abs(d[complement]).sum()),
            'structure_tets':int((~complement).sum()),'structure_volume_m3':float(np.abs(d[~complement]).sum()),
            'per_structure':conformity,
            'max_abs_structure_volume_error':float(max(abs(c['relative_volume_error']) for c in conformity)),
            'shared_nodes':'one vertex array for the whole box, so every structure/complement interface shares nodes by construction',
            'conformity_test':'tets are labelled by centroid winding; if the mesh conforms, the labelled volume equals the surface divergence volume'}
    return report

def cluster_trial(rows,pad,flags):
    by={r['entity_id']:r for r in rows};parts=[]
    for eid in CLUSTER:
        V,F=mesh(by[eid]['path'])
        parts.append({'entity_id':eid,'name':by[eid]['name'],'system':by[eid]['system'],
            'vertices':int(len(V)),'faces':int(len(F)),'surface_volume_m3':abs(divergence(V,F)),'V':V,'F':F})
    pairwise=[]
    crossing=np.zeros((len(parts),len(parts)),bool)
    for a in range(len(parts)):
        for b in range(a+1,len(parts)):
            n=int(len(cgal.intersect_other(parts[a]['V'],parts[a]['F'],parts[b]['V'],parts[b]['F'],True,False,False,False,100000)[0]))
            pairwise.append({'a':parts[a]['entity_id'],'b':parts[b]['entity_id'],'intersecting_face_pairs':n})
            crossing[a,b]=crossing[b,a]=n>0
    self_int=[{'entity_id':p['entity_id'],
        'self_intersecting_face_pairs':int(len(cgal.remesh_self_intersections(p['V'],p['F'],True,False,False,False,100000)[2]))} for p in parts]
    keep=[]
    for i in sorted(range(len(parts)),key=lambda i:-parts[i]['surface_volume_m3']):
        if not any(crossing[i,j] for j in keep):keep.append(i)
    return {'cluster':[{k:p[k] for k in ('entity_id','name','system','vertices','faces','surface_volume_m3')} for p in parts],
        'box_padding_m':pad,'pairwise_intersections':pairwise,'self_intersections':self_int,
        'mutually_intersecting_pairs':int(sum(1 for p in pairwise if p['intersecting_face_pairs'])),
        'pairs':len(pairwise),
        'full_cluster':enclose(parts,pad,flags),
        'non_intersecting_subset':enclose([parts[i] for i in sorted(keep)],pad,flags),
        'subset_rule':'greedy maximal set with no mutual intersection, largest surface volume first'}

def run(out,force,sample,seed,max_faces,pad,flags,skip_pairs):
    out=Path(out)
    if out.exists() and any(out.iterdir()) and not force:raise ValueError('Choose a fresh output directory or pass --force')
    began=time.monotonic();rows=sources()
    print('sources',len(rows),flush=True)
    coin=coincidence(rows);print('coincidence done',flush=True)
    pairs=None if skip_pairs else pair_audit(rows,sample,seed,max_faces)
    print('pair audit done',flush=True)
    cluster=cluster_trial(rows,pad,flags);print('cluster done',flush=True)
    conditions=[
        {'condition':'every facet set closed and consistently oriented',
         'status':'2376 of 2404 close under exact welding; 28 do not',
         'evidence':'data/derived/material-domains/whole-body-0.01m/manifest.json sources.not_closed_ids'},
        {'condition':'no facet set self-intersecting',
         'status':'425 repaired muscles are self-intersection-free; the 1979 canonical non-muscle surfaces have never been certified',
         'evidence':'data/derived/muscle-tet-ready-v1/summary.json after.self_intersection_free; this probe self-tests only the cluster'},
        {'condition':'no two facets from different structures crossing',
         'status':'measured by sample here','evidence':'pair_audit'},
        {'condition':'no coincident facets','status':'measured exhaustively here','evidence':'coincidence'},
        {'condition':'an outer bound enclosing the complement',
         'status':'available: watertight genus-0 outer envelope, 69.720 L',
         'evidence':'data/derived/outer-envelope/validation.json'}]
    summary={'schema':'ihm.conforming-complement-plc-probe.v1','wall_seconds':time.monotonic()-began,
        'plc_conditions':conditions,'coincidence':coin,'pair_audit':pairs,'cluster_trial':cluster,
        'inputs_sha256':{'data/derived/canonical/anatomy.json':sha(ROOT/'data/derived/canonical/anatomy.json'),
            'data/derived/muscle-tet-ready-v1/manifest.json':sha(MUSCLE/'manifest.json'),
            'ihm/assembly/material_domains.py':sha(ROOT/'ihm/assembly/material_domains.py'),
            'scripts/probe_conforming_complement_plc.py':sha(__file__)},
        'python':sys.version,'unverified':[
            'no global tetrahedralization of the whole body was attempted',
            'cross-structure intersection is sampled, not exhaustive, and pairs above the face cap are skipped',
            'non-muscle canonical surfaces are not certified self-intersection-free',
            'the cluster complement is bounded by a padded box, not by the outer envelope',
            'element quality of the cluster mesh is not measured beyond positive volume and volume conformity']}
    out.mkdir(parents=True,exist_ok=True)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (out/'manifest.json').write_text(json.dumps({'schema':'ihm.conforming-complement-plc-probe-manifest.v1',
        'inputs_sha256':summary['inputs_sha256'],'artifacts_sha256':{'summary.json':sha(out/'summary.json')},
        'canonical_assets_modified':False},indent=2)+'\n')
    print(json.dumps(summary,indent=2,allow_nan=False))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default=str(ROOT/'data/derived/conforming-complement-probe-v1'))
    p.add_argument('--force',action='store_true');p.add_argument('--sample',type=int,default=250)
    p.add_argument('--seed',type=int,default=0);p.add_argument('--max-faces',type=int,default=60000)
    p.add_argument('--pad',type=float,default=.02);p.add_argument('--flags',default='pY')
    p.add_argument('--skip-pairs',action='store_true');a=p.parse_args()
    run(a.out,a.force,a.sample,a.seed,a.max_faces,a.pad,a.flags,a.skip_pairs)
