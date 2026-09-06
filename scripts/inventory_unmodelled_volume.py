"""Inventory the canonical body interior that no entity mesh claims, against the outer skin envelope."""
from pathlib import Path
import argparse,collections,gzip,hashlib,json,time
import numpy as np
import igl
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
ENVELOPE=ROOT/'data/derived/outer-envelope'
DEPTH_EDGES_MM=(0,2,4,8,16,24,32,48,64,96,128,1e9)
UNVERIFIED=('entity_occupancy_is_generalized_winding_number_over_0.5_not_a_certified_solid',
            'open_and_self_intersecting_source_surfaces_get_a_winding_number_interpretation_only',
            'void_is_geometric_absence_of_a_claiming_entity_not_evidence_of_anatomical_emptiness',
            'nearest_entity_attribution_is_a_voxel_surface_distance_label_not_an_anatomical_assignment',
            'skin_and_shell_layer_entities_claim_a_1.8mm_slab_that_a_voxel_grid_this_coarse_cannot_resolve')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as s:
        for b in iter(lambda:s.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def mesh_of(path):
    b=Path(path).read_bytes()
    d=json.loads(gzip.decompress(b) if str(path).endswith('.gz') else b)
    return (np.asarray(d['positions'],float).reshape(-1,3),
            np.ascontiguousarray(np.asarray(d['indices'],np.int64).reshape(-1,3)),hashlib.sha256(b).hexdigest())

def divergence_volume(V,F):
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]]
    return float(np.einsum('ij,ij->i',A,np.cross(B,C)).sum()/6),float(np.linalg.norm(np.cross(B-A,C-A),axis=1).sum()/2)

def sample(V,F,h):
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]];e1=B-A;e2=C-A
    L=np.maximum.reduce([np.linalg.norm(e1,axis=1),np.linalg.norm(e2,axis=1),np.linalg.norm(C-B,axis=1)])
    n=np.ceil(L/(h/2.)).astype(int)+1;out=[]
    for k in np.unique(n):
        idx=np.where(n==k)[0];g=np.arange(k+1)/k
        uu,vv=np.meshgrid(g,g,indexing='ij');m=(uu+vv)<=1+1e-12;uu=uu[m];vv=vv[m];cs=max(1,int(2e6/len(uu)))
        for s in range(0,len(idx),cs):
            b=idx[s:s+cs]
            out.append((A[b][:,None,:]+uu[None,:,None]*e1[b][:,None,:]+vv[None,:,None]*e2[b][:,None,:]).reshape(-1,3))
    return np.concatenate(out) if out else np.zeros((0,3))

def occupancy(entities,lo,h,dim,interior,collect_owner=True):
    """Per-entity generalized winding number over the entity bounding box; counts stack into the grid."""
    count=np.zeros(tuple(dim),np.uint16);owner=np.full(tuple(dim),-1,np.int32) if collect_owner else None
    rows=[];axes=[lo[i]+h*np.arange(dim[i]) for i in range(3)]
    for k,(e,V,F) in enumerate(entities):
        a=np.maximum(((V.min(0)-h-lo)/h).astype(int),0);b=np.minimum(((V.max(0)+h-lo)/h).astype(int)+1,dim-1)
        vox=0;inside_vox=0
        if (b>=a).all():
            sl=[np.arange(a[i],b[i]+1) for i in range(3)]
            Q=np.ascontiguousarray(np.stack(np.meshgrid(*[axes[i][sl[i]] for i in range(3)],indexing='ij'),-1).reshape(-1,3))
            w=igl.fast_winding_number(V,F,Q)>0.5
            if w.any():
                ii=np.stack(np.meshgrid(*sl,indexing='ij'),-1).reshape(-1,3)[w]
                count[ii[:,0],ii[:,1],ii[:,2]]+=1;vox=int(w.sum())
                inside_vox=int(interior[ii[:,0],ii[:,1],ii[:,2]].sum())
        if collect_owner:
            P=sample(V,F,h)
            if len(P):
                ij=np.clip(((P-lo)/h+0.5).astype(np.int32),0,dim-1);owner[ij[:,0],ij[:,1],ij[:,2]]=k
        vol,area=divergence_volume(V,F)
        rows.append({'entity_id':e['id'],'name':e['name'],'system':e.get('system'),'role':e.get('role'),
                     'watertight_edge_incidence':e.get('watertight_edge_incidence'),
                     'signed_divergence_volume_m3':vol,'surface_area_m2':area,
                     'winding_voxel_volume_m3':vox*h**3,'winding_voxel_volume_inside_envelope_m3':inside_vox*h**3})
        if (k+1)%400==0:print('  occupancy %d/%d'%(k+1,len(entities)),flush=True)
    return count,owner,rows

def components(void,h,lo,depth,owner_lab,names,limit):
    lab,n=ndimage.label(void,ndimage.generate_binary_structure(3,1))
    cnt=np.bincount(lab.ravel());cnt[0]=0;order=np.argsort(-cnt)
    out=[]
    for k in order[:limit]:
        if cnt[k]==0:break
        m=lab==k;idx=np.argwhere(m);p=idx*h+lo;dd=depth[m]
        near=collections.Counter(owner_lab[m].tolist())
        top=[{'entity':names[i] if i>=0 else None,'volume_m3':c*h**3} for i,c in near.most_common(6)]
        out.append({'component':int(k),'voxels':int(cnt[k]),'volume_m3':float(cnt[k]*h**3),
                    'volume_L':float(cnt[k]*h**3*1000),'centroid_m':p.mean(0).tolist(),
                    'min_m':p.min(0).tolist(),'max_m':p.max(0).tolist(),
                    'depth_below_skin_mm':{'mean':float(dd.mean()*1000),'median':float(np.median(dd)*1000),
                                           'p05':float(np.percentile(dd,5)*1000),'p95':float(np.percentile(dd,95)*1000),
                                           'max':float(dd.max()*1000)},
                    'nearest_entity_share':top})
    return out,int(n),cnt

def histogram(depth_mm,weight_m3):
    edges=np.asarray(DEPTH_EDGES_MM,float);b=np.digitize(depth_mm,edges)-1
    v=np.bincount(b,minlength=len(edges)-1)[:len(edges)-1]*weight_m3
    return [{'lo_mm':float(edges[i]),'hi_mm':(None if edges[i+1]>1e8 else float(edges[i+1])),
             'volume_m3':float(v[i]),'fraction':float(v[i]/max(v.sum(),1e-30))} for i in range(len(v))]

def analyse(XV,TT,entities,h,detail):
    lo=XV.min(0)-2*h;dim=np.ceil((XV.max(0)+2*h-lo)/h).astype(int)+1
    axes=[lo[i]+h*np.arange(dim[i]) for i in range(3)]
    Q=np.ascontiguousarray(np.stack(np.meshgrid(*axes,indexing='ij'),-1).reshape(-1,3))
    interior=(igl.fast_winding_number(XV,TT,Q)>0.5).reshape(tuple(dim))
    count,owner,rows=occupancy(entities,lo,h,dim,interior,collect_owner=detail)
    occ=count>0;void=interior&~occ
    res={'grid_h_m':h,'grid_dim':dim.tolist(),'grid_origin_m':lo.tolist(),
         'interior_volume_m3':float(interior.sum()*h**3),
         'occupied_volume_m3':float((interior&occ).sum()*h**3),
         'void_volume_m3':float(void.sum()*h**3),
         'void_fraction_of_interior':float(void.sum()/max(interior.sum(),1)),
         'entity_claims_outside_envelope_m3':float((occ&~interior).sum()*h**3)}
    oc=count[interior&occ]
    res['overlap']={'mean_entities_per_occupied_voxel':float(oc.mean()) if oc.size else 0.0,
                    'fraction_of_occupied_volume_with_3_or_more':float((oc>=3).mean()) if oc.size else 0.0,
                    'max_entities_in_one_voxel':int(oc.max()) if oc.size else 0,
                    'histogram':{str(k):float(v*h**3) for k,v in zip(*[x.tolist() for x in np.unique(oc,return_counts=True)])},
                    'summed_entity_claim_volume_m3':float(count[interior].sum()*h**3)}
    if not detail:return res,None,None,None,None,None
    idx=np.argwhere(void);P=np.ascontiguousarray(idx*h+lo)
    sq,_,_=igl.point_mesh_squared_distance(P,XV,TT)
    depth=np.zeros(tuple(dim));depth[void]=np.sqrt(sq)
    names=[e['name']+' ['+e['id']+']' for e,_,_ in entities]
    _,ind=ndimage.distance_transform_edt(owner<0,sampling=(h,h,h),return_indices=True)
    olab=owner[ind[0],ind[1],ind[2]]
    comps,ncomp,cnt=components(void,h,lo,depth,olab,names,120)
    res['void_components']={'count':ncomp,'largest_volume_m3':float(cnt.max()*h**3),
                            'largest_fraction_of_void':float(cnt.max()/max(void.sum(),1)),
                            'components_over_0.1L':int((cnt*h**3>1e-4).sum())}
    res['depth_histogram']=histogram(depth[void]*1000,h**3)
    sysshare=collections.Counter()
    for i,c in collections.Counter(olab[void].tolist()).items():
        sysshare[(entities[i][0].get('system') if i>=0 else 'none')]+=c*h**3
    res['void_volume_by_nearest_entity_system_m3']={k:float(v) for k,v in sorted(sysshare.items(),key=lambda x:-x[1])}
    top=collections.Counter()
    for i,c in collections.Counter(olab[void].tolist()).items():top[names[i] if i>=0 else 'none']+=c*h**3
    res['void_volume_by_nearest_entity_m3']=dict(sorted(({k:float(v) for k,v in top.items()}).items(),key=lambda x:-x[1])[:40])
    return res,comps,rows,void,depth,(lo,dim,count,interior)

def run(out,force,primary,sensitivity):
    out=Path(out)
    if out.exists() and any(out.iterdir()) and not force:raise ValueError('Choose a fresh output directory or pass --force')
    out.mkdir(parents=True,exist_ok=True);began=time.monotonic()
    env=np.load(ENVELOPE/'outer-envelope.npz');XV=np.ascontiguousarray(env['positions']);TT=np.ascontiguousarray(env['indices'].astype(np.int64))
    ab=(ROOT/'data/derived/canonical/anatomy.json').read_bytes();an=json.loads(ab)
    ents=[];shas={};firsts=[];skipped=collections.Counter()
    for e in an['entities']:
        ref=e['reference_geometry']
        if ref.get('representation')!='triangular_surface' or ref.get('units')!='m' or ref.get('frame')!=an['frame']['id']:
            skipped[str(ref.get('representation'))]+=1;continue
        V,F,s=mesh_of(ROOT/ref['path'])
        if s!=ref['sha256']:raise ValueError('Geometry sha256 mismatch for '+e['id'])
        if s not in shas:firsts.append(len(ents))
        shas.setdefault(s,[]).append(e['id']);ents.append((e,V,F))
    print('entities',len(ents),'unique meshes',len(shas),'skipped',dict(skipped),flush=True)
    envvol,envarea=divergence_volume(XV,TT)
    res,comps,rows,void,depth,grid=analyse(XV,TT,ents,primary,True)
    sens=[res|{'void_components':None,'depth_histogram':None,'void_volume_by_nearest_entity_m3':None,
               'void_volume_by_nearest_entity_system_m3':None}]
    for h in sensitivity:
        if abs(h-primary)<1e-12:continue
        r,_,_,_,_,_=analyse(XV,TT,ents,h,False);sens.append(r);print('sensitivity',h,round(r['void_fraction_of_interior'],4),flush=True)
    lo,dim,count,interior=grid
    summary={'schema':'ihm.unmodelled-volume.v1','primary_grid_h_m':primary,
             'envelope':{'path':'data/derived/outer-envelope/outer-envelope.npz','faces':int(len(TT)),
                         'divergence_volume_m3':envvol,'surface_area_m2':envarea,
                         'sha256':sha(ENVELOPE/'outer-envelope.npz')},
             'entities_considered':len(ents),'unique_geometry_files':len(shas),
             'entities_sharing_geometry':int(sum(len(v) for v in shas.values() if len(v)>1)),
             'skipped_by_representation':dict(skipped),
             'summed_abs_divergence_volume_unique_meshes_m3':float(sum(abs(divergence_volume(ents[i][1],ents[i][2])[0]) for i in firsts)),
             'summed_surface_area_unique_meshes_m2':float(sum(divergence_volume(ents[i][1],ents[i][2])[1] for i in firsts)),
             'primary':res,'resolution_sensitivity':[{k:r[k] for k in ('grid_h_m','interior_volume_m3','occupied_volume_m3',
                 'void_volume_m3','void_fraction_of_interior','overlap')} for r in sens],
             'unverified':list(UNVERIFIED),'canonical_assets_modified':False,'wall_seconds':time.monotonic()-began}
    summary['entities_claiming_no_voxel_at_primary_h']=int(sum(1 for r in rows if r['winding_voxel_volume_m3']==0))
    summary['void_components_listed_in_jsonl']=len(comps)
    summary['void_components_unlisted_volume_m3']=float(res['void_volume_m3']-sum(c['volume_m3'] for c in comps))
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    with (out/'void_components.jsonl').open('w') as f:
        for c in comps:f.write(json.dumps(c,allow_nan=False)+'\n')
    with (out/'entity_occupancy.jsonl').open('w') as f:
        for r in sorted(rows,key=lambda r:-r['winding_voxel_volume_inside_envelope_m3']):f.write(json.dumps(r,allow_nan=False)+'\n')
    np.savez_compressed(out/'void_field.npz',origin_m=lo,h_m=primary,void=np.packbits(void),
                        interior=np.packbits(interior),entity_count=count.astype(np.uint8),dim=np.asarray(dim))
    (out/'manifest.json').write_text(json.dumps({'schema':'ihm.unmodelled-volume-manifest.v1',
        'inputs_sha256':{'data/derived/canonical/anatomy.json':hashlib.sha256(ab).hexdigest(),
                         'data/derived/outer-envelope/outer-envelope.npz':sha(ENVELOPE/'outer-envelope.npz'),
                         'data/derived/outer-envelope/manifest.json':sha(ENVELOPE/'manifest.json')},
        'geometry_sha256_count':len(shas),'builder_sha256':sha(__file__),
        'artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},
        'canonical_assets_modified':False},indent=2)+'\n')
    print(json.dumps({'output_dir':str(out),'primary':{k:res[k] for k in ('grid_h_m','interior_volume_m3','occupied_volume_m3','void_volume_m3','void_fraction_of_interior','overlap','void_components')},
                      'depth_histogram':res['depth_histogram'],'top_components':comps[:8]},indent=2))
    return summary

def self_test():
    o=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
    q=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,5,1],[0,4,5],[1,6,2],[1,5,6],[2,7,3],[2,6,7],[3,4,0],[3,7,4]])[:,::-1]
    v,a=divergence_volume(o,q);assert abs(v-1)<1e-12 and abs(a-6)<1e-12
    cube=lambda c,s:(o*s+c,q)
    e=lambda i:{'id':'e%d'%i,'name':'e%d'%i,'system':'test','role':'test','watertight_edge_incidence':True}
    inner=[(e(0),)+cube(np.array([.1,.1,.1]),.3),(e(1),)+cube(np.array([.15,.15,.15]),.3),(e(2),)+cube(np.array([.6,.6,.6]),.3)]
    inner=[(x[0],np.ascontiguousarray(x[1]),np.ascontiguousarray(x[2])) for x in inner]
    r,comps,rows,void,depth,_=analyse(np.ascontiguousarray(o),np.ascontiguousarray(q),inner,0.02,True)
    assert abs(r['interior_volume_m3']-1)<0.02,r['interior_volume_m3']
    claimed=0.3**3*3-0.25**3  # two overlapping 0.3 cubes offset by 0.05 in each axis, plus a disjoint one
    assert abs(r['occupied_volume_m3']-claimed)<0.004,(r['occupied_volume_m3'],claimed)
    assert abs(r['void_volume_m3']-(1-claimed))<0.02
    assert r['overlap']['max_entities_in_one_voxel']==2 and 1.0<r['overlap']['mean_entities_per_occupied_voxel']<1.6
    assert abs(r['overlap']['summed_entity_claim_volume_m3']-3*0.3**3)<0.004
    assert r['void_components']['count']==1 and sum(h['volume_m3'] for h in r['depth_histogram'])>0
    assert abs(sum(h['volume_m3'] for h in r['depth_histogram'])-r['void_volume_m3'])<1e-12
    assert comps[0]['depth_below_skin_mm']['max']>0 and comps[0]['nearest_entity_share'][0]['entity'] is not None
    assert len(rows)==3 and abs(rows[0]['signed_divergence_volume_m3']-0.027)<1e-9
    far=analyse(np.ascontiguousarray(o),np.ascontiguousarray(q),
                [(e(9),np.ascontiguousarray(o*.3+5.),np.ascontiguousarray(q))],0.02,False)[0]
    assert far['void_fraction_of_interior']==1.0 and far['occupied_volume_m3']==0
    edge=analyse(np.ascontiguousarray(o),np.ascontiguousarray(q),
                 [(e(8),np.ascontiguousarray(o*.4+.8),np.ascontiguousarray(q))],0.02,False)[0]
    assert abs(edge['occupied_volume_m3']-0.008)<0.002 and edge['entity_claims_outside_envelope_m3']>0
    print('PASS divergence volumes, winding occupancy, overlap counting, void components and depth binning')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'data/derived/unmodelled-volume')
    p.add_argument('--h',type=float,default=0.008);p.add_argument('--sensitivity',type=float,nargs='*',default=[0.012,0.008,0.005,0.004])
    p.add_argument('--force',action='store_true');p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:self_test()
    else:run(a.output,a.force,a.h,a.sensitivity)
