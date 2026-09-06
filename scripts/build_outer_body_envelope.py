"""Extract a watertight outward outer-skin envelope from the double-sided canonical skin slab."""
from pathlib import Path
import argparse,gzip,hashlib,json,time
import numpy as np
import igl
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parents[1]
SKIN='body-bp3d-FJ2810'
DETECT_H=0.002
PROBE_DEPTHS_M=(0.0008,0.0012,0.0016,0.0020,0.0026,0.0034,0.0045,0.0060,0.0080,0.0120)
BULK_SEED_M=(0.0,0.05,0.0)
UNVERIFIED=('triangle_self_intersection_not_exhaustively_tested',
            'aperture_caps_are_fan_surfaces_not_anatomical_geometry',
            'sheet_assignment_oracle_is_a_2mm_voxel_flood_fill_not_an_exact_visibility_test',
            'envelope_is_the_outer_face_of_an_authored_1.8mm_skin_slab_not_a_measured_body_surface')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as s:
        for b in iter(lambda:s.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def weld(V,F):
    X,inv=np.unique(np.asarray(V,float),axis=0,return_inverse=True);T=inv[np.asarray(F,np.int64)]
    return X,T[~(np.diff(np.sort(T,axis=1),axis=1)==0).any(1)]

def halfedges(T):
    d=np.concatenate((T[:,[0,1]],T[:,[1,2]],T[:,[2,0]]))
    u,inv,c=np.unique(np.sort(d,axis=1),axis=0,return_inverse=True,return_counts=True)
    return d,u,inv,c

def edge_faces(T):
    """Two incident face ids per edge; only defined when every edge carries exactly two faces."""
    d,u,inv,c=halfedges(T)
    if not np.all(c==2):raise ValueError('Edge-face pairing needs every edge at incidence 2')
    o=np.argsort(inv,kind='stable');fi=np.tile(np.arange(len(T)),3)[o]
    return np.stack((fi[0::2],fi[1::2]),1),u

def fan_counts(T,nv):
    """Faces around each vertex, glued across shared edges: one fan per vertex means vertex-manifold."""
    d,u,inv,c=halfedges(T);pair=np.where(c==2)[0]
    corner=lambda f,v:3*f+np.argmax(T[f]==v[:,None],axis=1)
    rows=[];cols=[]
    if len(pair):
        o=np.argsort(inv,kind='stable');fi=np.tile(np.arange(len(T)),3)[o];pos=np.searchsorted(inv[o],pair)
        f0=fi[pos];f1=fi[pos+1]
        for v in (u[pair,0],u[pair,1]):
            rows.append(corner(f0,v));cols.append(corner(f1,v))
    r=np.concatenate(rows) if rows else np.zeros(0,np.int64)
    cc=np.concatenate(cols) if cols else np.zeros(0,np.int64)
    g=coo_matrix((np.ones(len(r)),(r,cc)),shape=(3*len(T),3*len(T)))
    _,lab=connected_components(g,directed=False)
    vid=np.empty(3*len(T),np.int64)
    for k in range(3):vid[3*np.arange(len(T))+k]=T[:,k]
    keys=np.unique(np.stack((vid,lab),1),axis=0)
    return np.bincount(keys[:,0],minlength=nv)

def topology(X,T):
    d,u,inv,c=halfedges(T)
    bal=np.bincount(inv,weights=np.where(d[:,0]<d[:,1],1,-1),minlength=len(u))
    A,B,C=X[T[:,0]],X[T[:,1]],X[T[:,2]];n=np.cross(B-A,C-A);a2=np.linalg.norm(n,axis=1)
    used=np.unique(T);fans=fan_counts(T,len(X))
    _,vlab=connected_components(coo_matrix((np.ones(len(u)),(u[:,0],u[:,1])),shape=(len(X),len(X))),directed=False)
    return {'vertices_used':int(len(used)),'vertices_stored':int(len(X)),'faces':int(len(T)),'edges':int(len(u)),
            'edge_incidence_histogram':{str(k):int(v) for k,v in zip(*[x.tolist() for x in np.unique(c,return_counts=True)])},
            'watertight_every_edge_incidence_2':bool(np.all(c==2)),
            'nonmanifold_edges':int((c>2).sum()),'boundary_edges':int((c==1).sum()),
            'winding_conflicts':int(((c==2)&(bal!=0)).sum()),'consistent_winding':bool(not ((c==2)&(bal!=0)).any()),
            'nonmanifold_vertices':int((fans[used]!=1).sum()),'vertex_link_single_fan':bool(bool((fans[used]==1).all())),
            'euler_characteristic':int(len(used)-len(u)+len(T)),
            'connected_components':int(len(np.unique(vlab[used]))),
            'zero_area_faces':int((a2==0).sum()),
            'duplicate_faces':int(np.sum(np.unique(np.sort(T,axis=1),axis=0,return_counts=True)[1]-1)),
            'surface_area_m2':float(a2.sum()/2),
            'signed_volume_m3':float(np.einsum('ij,ij->i',A,np.cross(B,C)).sum()/6)}

def rasterize(X,T,lo,h,dim):
    """Barycentric point sampling at h/3 spacing; dense enough to block 6-connected flood fill."""
    O=np.zeros(tuple(dim),bool);A,B,C=X[T[:,0]],X[T[:,1]],X[T[:,2]];e1=B-A;e2=C-A
    L=np.maximum.reduce([np.linalg.norm(e1,axis=1),np.linalg.norm(e2,axis=1),np.linalg.norm(C-B,axis=1)])
    n=np.ceil(L/(h/3.)).astype(int)+1
    for k in np.unique(n):
        idx=np.where(n==k)[0];g=np.arange(k+1)/k
        uu,vv=np.meshgrid(g,g,indexing='ij');m=(uu+vv)<=1+1e-12;uu=uu[m];vv=vv[m];cs=max(1,int(4e6/len(uu)))
        for s in range(0,len(idx),cs):
            b=idx[s:s+cs]
            P=(A[b][:,None,:]+uu[None,:,None]*e1[b][:,None,:]+vv[None,:,None]*e2[b][:,None,:]).reshape(-1,3)
            ij=((P-lo)/h+0.5).astype(np.int32);O[ij[:,0],ij[:,1],ij[:,2]]=True
    return O

def grid_of(X,h,pad=4):
    lo=X.min(0)-pad*h;return lo,np.ceil((X.max(0)+pad*h-lo)/h).astype(int)+1

def apertures(X,T,h=DETECT_H):
    """Smallest free-space erosion that severs body bulk from ambient air localises the slab apertures."""
    lo,dim=grid_of(X,h);O=rasterize(X,T,lo,h,dim);S1=ndimage.generate_binary_structure(3,1);F=~O
    seed=tuple(((np.asarray(BULK_SEED_M)-lo)/h+0.5).astype(int))
    if O[seed]:raise ValueError('Bulk seed landed on the slab')
    E=F.copy();r=0
    while True:
        r+=1;E=ndimage.binary_erosion(E,S1,border_value=1)
        l,_=ndimage.label(E,S1);cnt=np.bincount(l.ravel());cnt[0]=0;b=l[seed];a=int(np.argmax(cnt))
        if b and b!=a:break
        if r>12:raise ValueError('Free space never separated; slab is not a closed two-sheet shell')
    J=ndimage.binary_dilation(l==b,S1,iterations=r+2)&ndimage.binary_dilation(l==a,S1,iterations=r+2)&F
    jl,jn=ndimage.label(J,ndimage.generate_binary_structure(3,3));jc=np.bincount(jl.ravel());jc[0]=0
    boxes=[]
    for k in range(1,jn+1):
        if jc[k]<10:continue
        p=np.argwhere(jl==k)*h+lo
        boxes.append({'voxels':int(jc[k]),'min_m':p.min(0).tolist(),'max_m':p.max(0).tolist(),'centroid_m':p.mean(0).tolist()})
    if not boxes:raise ValueError('Bulk separated but no aperture throat isolated')
    return {'detect_h_m':h,'separating_erosion_radius_voxels':int(r),'throat_voxels':int(J.sum()),
            'discarded_throat_fragments_under_10_voxels':int(jn-len(boxes)),'apertures':boxes}

def solid_field(X,T,boxes,h,margin=0.002):
    lo,dim=grid_of(X,h);O=rasterize(X,T,lo,h,dim);P=np.zeros(tuple(dim),bool)
    for bx in boxes:
        P[tuple(slice(max(0,int((bx['min_m'][i]-margin-lo[i])/h)),min(int(dim[i]),int((bx['max_m'][i]+margin-lo[i])/h)+2)) for i in range(3))]=True
    S1=ndimage.generate_binary_structure(3,1)
    lab,_=ndimage.label(~(O|P),S1);solid=lab!=lab[0,0,0]
    return lo,dim,(O|P),solid,{'grid_h_m':h,'grid_dim':dim.tolist(),'plug_voxels':int(P.sum()),
        'plug_volume_m3':float(P.sum()*h**3),'slab_voxel_volume_m3':float(O.sum()*h**3),
        'sealed_solid_volume_m3':float(solid.sum()*h**3),'sealed_solid_components':int(ndimage.label(solid,S1)[1])}

def classify(cent,N,lo,h,dim,occ,solid,steps=8):
    """March outward one voxel at a time; the first non-slab voxel says air (outer) or bulk (inner)."""
    out=np.zeros(len(cent),np.int8)
    for s in range(1,steps+1):
        ij=np.clip(((cent+(s*h)*N-lo)/h+0.5).astype(np.int32),0,np.asarray(dim)-1)
        o=occ[ij[:,0],ij[:,1],ij[:,2]];q=solid[ij[:,0],ij[:,1],ij[:,2]];new=(out==0)&(~o)
        out[new]=np.where(q[new],-1,1)
    return out

def pair_faces(X,T,cent,N):
    best=np.full(len(T),-1,np.int64);dist=np.full(len(T),np.inf)
    for d in PROBE_DEPTHS_M:
        _,idx,cp=igl.point_mesh_squared_distance(np.ascontiguousarray(cent-d*N),X,T)
        ok=(np.einsum('ij,ij->i',N,N[idx])<-0.6)&(idx!=np.arange(len(T)))
        t=np.linalg.norm(cp-cent,axis=1);up=ok&(t<dist);best[up]=idx[up];dist[up]=t[up]
    mutual=np.where((best>=0)&(best[np.maximum(best,0)]==np.arange(len(T))),best,-1)
    return best,dist,mutual

def smooth(lab0,mutual,ef,rounds=12,prior=1.0):
    """Neighbour majority under the constraint that a mutual anti-parallel pair straddles the slab."""
    f0,f1=ef[:,0],ef[:,1];NF=len(lab0);lab=lab0.copy();hist=[]
    tot=np.bincount(f0,minlength=NF)+np.bincount(f1,minlength=NF)
    for _ in range(rounds):
        o=lab.astype(float)
        sc=2*(np.bincount(f0,weights=o[f1],minlength=NF)+np.bincount(f1,weights=o[f0],minlength=NF))-tot+prior*(2*lab0-1)
        new=sc>0;p=mutual>=0;q=mutual[p];d=sc[p]-sc[q]
        new[p]=np.where(d!=0,d>0,lab[p]);hist.append(int((new!=lab).sum()));lab=new
        if hist[-1]==0:break
    p=mutual>=0
    return lab,{'icm_rounds':len(hist),'icm_changes_per_round':hist,'mutually_paired_faces':int(p.sum()),
                'pair_exclusivity_violations':int((lab[p]==lab[mutual[p]]).sum())}

def largest_component(T,lab,ef):
    idx=np.where(lab)[0];ren=-np.ones(len(T),np.int64);ren[idx]=np.arange(len(idx))
    m=lab[ef[:,0]]&lab[ef[:,1]]
    _,cl=connected_components(coo_matrix((np.ones(int(m.sum())),(ren[ef[m,0]],ren[ef[m,1]])),shape=(len(idx),len(idx))),directed=False)
    cnt=np.bincount(cl);keep=np.zeros(len(T),bool);keep[idx[cl==int(np.argmax(cnt))]]=True
    return keep,{'outer_face_components':int(len(cnt)),'largest_component_faces':int(keep.sum()),
                 'discarded_stray_faces':int(lab.sum()-keep.sum())}

def cap(X,T,keep):
    """Grow over pinch vertices, then cone every remaining boundary loop onto its centroid."""
    VF=coo_matrix((np.ones(3*len(T),np.int8),(T.ravel(),np.repeat(np.arange(len(T)),3))),shape=(len(X),len(T))).tocsr()
    grew=[]
    for _ in range(10):
        d,u,inv,c=halfedges(T[keep]);b=d[np.isin(inv,np.where(c==1)[0])]
        vs,vc=np.unique(b.ravel(),return_counts=True);pin=vs[vc>2]
        if not len(pin) and not (c>2).any():break
        grew.append(int(len(pin)));keep=keep.copy();keep[np.unique(VF[pin].tocoo().col)]=True
    S=T[keep];d,u,inv,c=halfedges(S);b=d[np.isin(inv,np.where(c==1)[0])]
    _,ll=connected_components(coo_matrix((np.ones(len(b)),(b[:,0],b[:,1])),shape=(len(X),len(X))),directed=False)
    ids=np.unique(ll[np.unique(b.ravel())]) if len(b) else np.zeros(0,int)
    caps=[];apex=[];loops=[]
    for k in ids:
        E=b[ll[b[:,0]]==k];vv=np.unique(E.ravel());ai=len(X)+len(apex)
        loops.append({'apex_index':int(ai),'boundary_edges':int(len(E)),'boundary_vertices':int(len(vv)),
                      'centroid_m':X[vv].mean(0).tolist(),'extent_mm':((X[vv].max(0)-X[vv].min(0))*1000).tolist()})
        apex.append(X[vv].mean(0));caps.extend((int(y),int(x),ai) for x,y in E)
    XV=np.vstack([X,np.asarray(apex)]) if apex else X
    TT=np.vstack([S,np.asarray(caps,np.int64).reshape(-1,3)]) if caps else S
    used=np.unique(TT);ren=-np.ones(len(XV),np.int64);ren[used]=np.arange(len(used))
    return XV[used],ren[TT],{'kept_source_faces':int(keep.sum()),'pinch_vertices_absorbed_per_round':grew,
                             'boundary_loops_capped':int(len(ids)),'cap_faces':len(caps),'apex_vertices':len(apex),
                             'loops':sorted(loops,key=lambda x:-x['boundary_edges'])}

def vertex_normals(X,T):
    A,B,C=X[T[:,0]],X[T[:,1]],X[T[:,2]];n=np.cross(B-A,C-A);out=np.zeros_like(X)
    for k in range(3):np.add.at(out,T[:,k],n)
    return out/np.maximum(np.linalg.norm(out,axis=1,keepdims=True),1e-30)

def run(out,force,voxel_check_h=0.004):
    out=Path(out)
    if out.exists() and any(out.iterdir()) and not force:raise ValueError('Choose a fresh output directory or pass --force')
    out.mkdir(parents=True,exist_ok=True);began=time.monotonic()
    ap=ROOT/'data/derived/canonical/anatomy.json';ab=ap.read_bytes();an=json.loads(ab)
    pf=json.loads((ROOT/'data/derived/canonical/profile.json').read_bytes())
    ent=next(e for e in an['entities'] if e['id']==SKIN);ref=ent['reference_geometry']
    gb=(ROOT/ref['path']).read_bytes();gs=hashlib.sha256(gb).hexdigest()
    if gs!=ref['sha256']:raise ValueError('Skin geometry sha256 disagrees with anatomy.json')
    d=json.loads(gzip.decompress(gb))
    X,T=weld(np.asarray(d['positions'],float).reshape(-1,3),np.asarray(d['indices'],np.int64).reshape(-1,3))
    src=topology(X,T);print('slab',{k:src[k] for k in ('faces','watertight_every_edge_incidence_2','euler_characteristic','surface_area_m2','signed_volume_m3')},flush=True)
    if not src['watertight_every_edge_incidence_2'] or not src['consistent_winding']:
        raise ValueError('Skin slab is not a closed consistently wound surface')
    A,B,C=X[T[:,0]],X[T[:,1]],X[T[:,2]];Nf=np.cross(B-A,C-A);area=np.linalg.norm(Nf,axis=1)/2
    Nf=Nf/(2*area[:,None]);cent=np.ascontiguousarray((A+B+C)/3)
    best,dist,mutual=pair_faces(X,T,cent,Nf);pr=best>=0
    pairing={'faces':int(len(T)),'faces_with_antiparallel_partner':int(pr.sum()),'fraction':float(pr.mean()),
             'partner_area_fraction':float(area[pr].sum()/area.sum()),'probe_depths_m':list(PROBE_DEPTHS_M),
             'slab_thickness_mm_percentiles':{str(q):float(np.percentile(dist[pr],q)*1000) for q in (1,5,25,50,75,95,99)}}
    print('pairing',round(pairing['fraction'],5),pairing['slab_thickness_mm_percentiles']['50'],flush=True)
    aps=apertures(X,T);print('apertures',aps['separating_erosion_radius_voxels'],len(aps['apertures']),flush=True)
    lo,dim,occ,solid,sol=solid_field(X,T,aps['apertures'],DETECT_H)
    if sol['sealed_solid_components']!=1 or not 0.03<sol['sealed_solid_volume_m3']<0.20:
        raise ValueError('Sealed flood fill produced no single plausible body solid: %r'%sol)
    print('solid',sol,flush=True)
    o=classify(cent,Nf,lo,DETECT_H,dim,occ,solid);raw=(o==1);p=mutual>=0
    ef,_=edge_faces(T)
    lab,icm=smooth(raw,mutual,ef)
    icm.update({'oracle_outer_faces':int(raw.sum()),'oracle_undetermined_faces':int((o==0).sum()),
                'oracle_pair_consistency':float((raw[p]!=raw[mutual[p]]).mean()),
                'final_outer_faces':int(lab.sum()),'final_outer_area_m2':float(area[lab].sum()),
                'final_inner_area_m2':float(area[~lab].sum())})
    keep,comp=largest_component(T,lab,ef);icm.update(comp)
    XV,TT,capr=cap(X,T,keep)
    env=topology(XV,TT)
    print('envelope',{k:env[k] for k in ('faces','watertight_every_edge_incidence_2','euler_characteristic','signed_volume_m3','surface_area_m2')},flush=True)
    glo,gdim=grid_of(XV,voxel_check_h,pad=2)
    Q=np.stack(np.meshgrid(*[glo[i]+voxel_check_h*np.arange(gdim[i]) for i in range(3)],indexing='ij'),-1).reshape(-1,3)
    vv=float((igl.fast_winding_number(XV,TT,np.ascontiguousarray(Q))>0.5).sum()*voxel_check_h**3)
    kg=float(pf['mass_kg']);hm=float(pf['height_m']);gkg=env['signed_volume_m3']*1010.0
    bsa=lambda m:0.007184*(hm*100)**0.725*m**0.425
    checks={'watertight_every_edge_incidence_2':env['watertight_every_edge_incidence_2'],
            'edge_incidence_histogram':env['edge_incidence_histogram'],
            'edge_manifold_no_edge_over_two_faces':env['nonmanifold_edges']==0,
            'vertex_manifold_single_fan_per_vertex':env['vertex_link_single_fan'],
            'nonmanifold_vertices':env['nonmanifold_vertices'],
            'consistent_winding_zero_conflicts':env['consistent_winding'],'winding_conflicts':env['winding_conflicts'],
            'euler_characteristic':env['euler_characteristic'],'genus':int((2-env['euler_characteristic'])//2),
            'single_connected_component':env['connected_components']==1,
            'zero_area_faces':env['zero_area_faces'],'duplicate_faces':env['duplicate_faces'],
            'outward_oriented_positive_volume':env['signed_volume_m3']>0,
            'enclosed_volume_m3':env['signed_volume_m3'],'enclosed_volume_L':env['signed_volume_m3']*1000,
            'independent_voxel_volume_m3':vv,'independent_voxel_h_m':voxel_check_h,
            'voxel_vs_divergence_relative_difference':abs(vv-env['signed_volume_m3'])/env['signed_volume_m3'],
            'outer_surface_area_m2':env['surface_area_m2'],
            'source_slab_area_m2':src['surface_area_m2'],
            'outer_over_slab_area_ratio':env['surface_area_m2']/src['surface_area_m2'],
            'profile_mass_kg':kg,'profile_height_m':hm,
            'target_volume_m3_profile_mass_over_1010':kg/1010.0,
            'volume_vs_profile_target_relative':env['signed_volume_m3']/(kg/1010.0)-1.0,
            'implied_mass_kg_at_1010':gkg,
            'dubois_bsa_m2_from_profile_mass':bsa(kg),'area_vs_dubois_profile_relative':env['surface_area_m2']/bsa(kg)-1.0,
            'dubois_bsa_m2_from_envelope_mass':bsa(gkg),'area_vs_dubois_envelope_relative':env['surface_area_m2']/bsa(gkg)-1.0}
    body={'positions':XV.reshape(-1).tolist(),'indices':TT.reshape(-1).tolist(),
          'normals':vertex_normals(XV,TT).reshape(-1).tolist(),'frame':an['frame']['id'],'units':'m',
          'representation':'triangular_surface',
          'derivation':'outer sheet of canonical skin slab %s with aperture caps'%SKIN}
    (out/'outer-envelope.json.gz').write_bytes(gzip.compress(json.dumps(body,allow_nan=False).encode(),9))
    np.savez_compressed(out/'outer-envelope.npz',positions=XV,indices=TT)
    rep={'schema':'ihm.outer-body-envelope.v1','source_entity':SKIN,'source_geometry_path':ref['path'],
         'source_geometry_sha256':gs,'anatomy_sha256':hashlib.sha256(ab).hexdigest(),
         'source_slab_topology':src,'antiparallel_pairing':pairing,'aperture_detection':aps,
         'sealed_flood_fill':sol,'sheet_assignment':icm,'capping':capr,'envelope_topology':env,
         'validation':checks,'unverified':list(UNVERIFIED),'wall_seconds':time.monotonic()-began}
    (out/'validation.json').write_text(json.dumps(rep,indent=2,allow_nan=False)+'\n')
    (out/'manifest.json').write_text(json.dumps({'schema':'ihm.outer-body-envelope-manifest.v1',
        'inputs_sha256':{ref['path']:gs,'data/derived/canonical/anatomy.json':hashlib.sha256(ab).hexdigest(),
                         'data/derived/canonical/profile.json':sha(ROOT/'data/derived/canonical/profile.json')},
        'builder_sha256':sha(__file__),
        'artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},
        'canonical_assets_modified':False},indent=2)+'\n')
    print(json.dumps({'output_dir':str(out),'validation':checks},indent=2))
    return rep

def self_test():
    v=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]]);t=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    r=topology(v,t)
    assert r['watertight_every_edge_incidence_2'] and r['consistent_winding'] and r['euler_characteristic']==2
    assert r['vertex_link_single_fan'] and abs(r['signed_volume_m3']-1/6)<1e-12 and r['connected_components']==1
    f=t.copy();f[0]=f[0,::-1];assert topology(v,f)['winding_conflicts']==3
    assert topology(v,t[:,::-1])['signed_volume_m3']<0
    two=topology(np.r_[v,v+5],np.r_[t,t+4]);assert two['euler_characteristic']==4 and two['connected_components']==2
    bow=topology(np.r_[v,-v[1:]],np.r_[t,np.array([[0,4,5],[0,5,6],[0,6,4],[4,6,5]])])
    assert bow['nonmanifold_vertices']==1 and not bow['vertex_link_single_fan'] and bow['watertight_every_edge_incidence_2']
    o=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
    q=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,5,1],[0,4,5],[1,6,2],[1,5,6],[2,7,3],[2,6,7],[3,4,0],[3,7,4]])[:,::-1]
    s=0.05;X=np.r_[o,(o-.5)*(1-2*s)+.5];T=np.r_[q,q[:,::-1]+8]
    sl=topology(X,T);assert sl['watertight_every_edge_incidence_2'] and abs(sl['signed_volume_m3']-(1-(1-2*s)**3))<1e-12
    ef,_=edge_faces(T);assert len(ef)==sl['edges']
    lab,rep=smooth(np.r_[np.ones(12,bool),np.zeros(12,bool)],np.full(24,-1,np.int64),ef)
    assert lab[:12].all() and not lab[12:].any()
    noisy=np.r_[np.ones(12,bool),np.zeros(12,bool)];noisy[0]=False;noisy[13]=True
    assert (smooth(noisy,np.full(24,-1,np.int64),ef)[0]==np.r_[np.ones(12,bool),np.zeros(12,bool)]).all()
    keep,_=largest_component(T,np.r_[np.ones(12,bool),np.zeros(12,bool)],ef)
    XV,TT,cr=cap(X,T,keep)
    ct=topology(XV,TT)
    assert ct['watertight_every_edge_incidence_2'] and ct['euler_characteristic']==2 and cr['boundary_loops_capped']==0
    assert abs(ct['signed_volume_m3']-1.0)<1e-12 and abs(ct['surface_area_m2']-6.0)<1e-12
    open_keep=keep.copy();open_keep[0]=False
    XO,TO,cr2=cap(X,T,open_keep)
    to=topology(XO,TO)
    assert cr2['boundary_loops_capped']==1 and to['watertight_every_edge_incidence_2'] and to['euler_characteristic']==2
    lo,dim=grid_of(X,0.05,pad=2);O=rasterize(X,T,lo,0.05,dim)
    lab2,_=ndimage.label(~O,ndimage.generate_binary_structure(3,1))
    assert lab2.max()>=2 and lab2[0,0,0]!=lab2[tuple(((np.array([.5,.5,.5])-lo)/0.05+.5).astype(int))]
    print('PASS topology receipts, bowtie detection, slab handling, smoothing, capping and watertight rasterization')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'data/derived/outer-envelope')
    p.add_argument('--force',action='store_true');p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:self_test()
    else:run(a.output,a.force)
