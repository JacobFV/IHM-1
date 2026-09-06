"""Representation error of the whole-body voxel partition against conforming muscle tetrahedra.

For every one of the 425 repaired muscle surfaces this measures, on the same grid the
whole-body partition uses, the cells the muscle occupies, its voxel volume and its voxel
boundary area, and compares them with that muscle's own surface divergence volume, surface
area and the retained TetGen conforming mesh. On a stratified sample spanning the size range
it also measures symmetric surface deviation between the voxel boundary and the true surface.

Requires the isolated libigl environment for exact point-triangle distance:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/measure_voxel_vs_conforming_muscle.py

Measured here: occupancy, volume, boundary area, cell counts, sampled symmetric surface
deviation. Not measured here: mechanical response of either representation, conforming tet
quality beyond what the retained TetGen trials already record, and anything outside muscle.
"""
from pathlib import Path
import argparse,gzip,hashlib,importlib.util,json,sys,time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
MUSCLE=ROOT/'data/derived/muscle-tet-ready-v1'
TRIALS=ROOT/'data/derived/muscle-tet-ready-verification-v1/tetgen-trials.jsonl'
DIRS=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
QUAD={0:((1,0,0),(1,1,0),(1,1,1),(1,0,1)),1:((0,0,0),(0,0,1),(0,1,1),(0,1,0)),
      2:((0,1,0),(0,1,1),(1,1,1),(1,1,0)),3:((0,0,0),(1,0,0),(1,0,1),(0,0,1)),
      4:((0,0,1),(1,0,1),(1,1,1),(0,1,1)),5:((0,0,0),(0,1,0),(1,1,0),(1,0,0))}

def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m
MD=load_module(ROOT/'ihm/assembly/material_domains.py','material_domains')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def surface(path):
    d=json.loads(gzip.decompress(Path(path).read_bytes()))
    return (np.ascontiguousarray(np.asarray(d['positions'],float).reshape(-1,3)),
            np.ascontiguousarray(np.asarray(d['indices'],np.int64).reshape(-1,3)))

def divergence(V,F):
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]]
    return float(np.einsum('ij,ij->i',A,np.cross(B,C)).sum()/6),float(np.linalg.norm(np.cross(B-A,C-A),axis=1).sum()/2)

def occupancy(V,F,origin,h):
    """Cell-center winding occupancy on the global lattice, restricted to the surface bounding box."""
    lo=np.floor((V.min(0)-origin)/h).astype(np.int64);hi=np.ceil((V.max(0)-origin)/h).astype(np.int64)
    cells=np.stack(np.meshgrid(*(np.arange(lo[i],hi[i]+1) for i in range(3)),indexing='ij'),-1).reshape(-1,3)
    centers=origin+(cells+.5)*h
    tree=MD.WindingHierarchy(V,F)
    box=np.all((centers>=tree.lower)&(centers<=tree.upper),axis=1)
    w=np.zeros(len(centers));w[box]=tree(centers[box])
    return cells[np.abs(w)>.5]

def boundary(cells,origin,h):
    """Exposed cube faces of the occupied set, as two triangles per face on shared lattice nodes."""
    if not len(cells):return np.zeros((0,3)),np.zeros((0,3),np.int64),0
    key=lambda a:(a[:,0]+(1<<20))*(1<<42)+(a[:,1]+(1<<20))*(1<<21)+(a[:,2]+(1<<20))
    present=set(key(cells).tolist());quads=[]
    for d in range(6):
        step=DIRS[d];mask=np.array([k not in present for k in key(cells+step).tolist()])
        if mask.any():
            corners=np.array(QUAD[d])
            quads.append(cells[mask][:,None,:]+corners[None,:,:])
    if not quads:return np.zeros((0,3)),np.zeros((0,3),np.int64),0
    quads=np.concatenate(quads)
    nodes,index=np.unique(quads.reshape(-1,3),axis=0,return_inverse=True);index=index.reshape(-1,4)
    tri=np.concatenate((index[:,[0,1,2]],index[:,[0,2,3]]))
    return np.ascontiguousarray(origin+nodes*h),np.ascontiguousarray(tri),len(quads)

def sample_surface(V,F,step,rng):
    """Deterministic barycentric lattice sampling at roughly the requested edge spacing."""
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]];e1=B-A;e2=C-A
    L=np.maximum.reduce([np.linalg.norm(e1,axis=1),np.linalg.norm(e2,axis=1),np.linalg.norm(C-B,axis=1)])
    n=np.clip(np.ceil(L/step).astype(int),1,12);out=[]
    for k in np.unique(n):
        idx=np.flatnonzero(n==k);g=(np.arange(k+1))/k
        u,v=np.meshgrid(g,g,indexing='ij');m=(u+v)<=1+1e-12;u=u[m];v=v[m]
        chunk=max(1,int(2e6/max(len(u),1)))
        for s in range(0,len(idx),chunk):
            b=idx[s:s+chunk]
            out.append((A[b][:,None,:]+u[None,:,None]*e1[b][:,None,:]+v[None,:,None]*e2[b][:,None,:]).reshape(-1,3))
    return np.ascontiguousarray(np.concatenate(out)) if out else np.zeros((0,3))

def deviation(V,F,XV,XF,h):
    if not len(XF):return None
    rng=np.random.default_rng(0)
    qs=sample_surface(V,F,h/4,rng);qv=sample_surface(XV,XF,h/4,rng)
    d1=np.sqrt(igl.point_mesh_squared_distance(qs,XV,XF)[0])
    d2=np.sqrt(igl.point_mesh_squared_distance(qv,V,F)[0])
    return {'true_samples':int(len(qs)),'voxel_samples':int(len(qv)),
        'mean_true_to_voxel_m':float(d1.mean()),'p95_true_to_voxel_m':float(np.percentile(d1,95)),
        'max_true_to_voxel_m':float(d1.max()),'mean_voxel_to_true_m':float(d2.mean()),
        'p95_voxel_to_true_m':float(np.percentile(d2,95)),'max_voxel_to_true_m':float(d2.max()),
        'symmetric_hausdorff_m':float(max(d1.max(),d2.max())),
        'symmetric_mean_m':float((d1.mean()+d2.mean())/2),
        'cell_diagonal_m':float(np.sqrt(3)*h),'hausdorff_over_cell_diagonal':float(max(d1.max(),d2.max())/(np.sqrt(3)*h)),
        'basis':'barycentric lattice samples at h/4 on each surface, exact point-triangle distance to the other; a sampled Hausdorff bound, not a certified one'}

def run(spacings,sample,out,force):
    out=Path(out)
    if out.exists() and any(out.iterdir()) and not force:raise ValueError('Choose a fresh output directory or pass --force')
    began=time.monotonic()
    trials={}
    for line in TRIALS.read_text().splitlines():
        r=json.loads(line);trials[r['entity_id']]=r
    entities=[json.loads(l) for l in (MUSCLE/'entities.jsonl').read_text().splitlines()]
    anatomy=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
    frame=anatomy['frame']['id']
    lows=[]
    for e in anatomy['entities']:
        ref=e['reference_geometry']
        if ref.get('representation')!='triangular_surface' or ref.get('units')!='m' or ref.get('frame')!=frame:continue
        lows.append(e['bounds_m']['min'])
    world_min=np.min(lows,axis=0)
    origins={h:np.floor(world_min/h+1e-10)*h for h in spacings}
    for h in spacings:
        built=ROOT/f'data/derived/material-domains/whole-body-{h:g}m/manifest.json'
        if built.exists():
            recorded=np.asarray(json.loads(built.read_text())['grid']['origin_m'])
            if not np.allclose(recorded,origins[h],atol=1e-12):raise ValueError(f'Grid origin disagrees with the built {h} m partition')
    meshes={}
    for e in entities:
        V,F=surface(MUSCLE/e['output_path']);vol,area=divergence(V,F);meshes[e['entity_id']]=(V,F,abs(vol),area)
    order=sorted(meshes,key=lambda k:meshes[k][2])
    picks=set(np.array(order)[np.unique(np.linspace(0,len(order)-1,sample).astype(int))].tolist())
    rows=[]
    for n,e in enumerate(entities):
        eid=e['entity_id'];V,F,vol,area=meshes[eid];t=trials.get(eid,{}).get('repaired',{})
        row={'entity_id':eid,'name':e['name'],'faces':int(len(F)),'vertices':int(len(V)),
            'surface_volume_m3':vol,'surface_area_m2':area,
            'conforming':{'tets':t.get('tets'),'tet_vertices':t.get('tet_vertices'),'tet_volume_m3':t.get('tet_volume_m3'),
                'succeeded':t.get('succeeded'),'flags':t.get('flags'),
                'tet_volume_relative_error':(abs(t['tet_volume_m3']-vol)/vol if t.get('tet_volume_m3') else None)},
            'deviation_sampled':eid in picks,'spacings':{}}
        for h in spacings:
            cells=occupancy(V,F,origins[h],h)
            XV,XF,faces=boundary(cells,origins[h],h)
            r={'cells':int(len(cells)),'voxel_volume_m3':float(len(cells)*h**3),
               'volume_relative_error':float(len(cells)*h**3/vol-1),
               'voxel_boundary_faces':int(faces),'voxel_boundary_area_m2':float(faces*h**2),
               'area_relative_error':float(faces*h**2/area-1),
               'voxel_tets_six_per_cell':int(6*len(cells)),
               'conforming_tets_per_voxel_tet':(t['tets']/(6*len(cells)) if t.get('tets') and len(cells) else None)}
            if eid in picks:r['deviation']=deviation(V,F,XV,XF,h)
            row['spacings'][f'{h:g}']=r
        rows.append(row)
        if (n+1)%50==0:print('  %d/%d'%(n+1,len(entities)),flush=True)
    summary={'schema':'ihm.voxel-vs-conforming-muscle.v1','entities':len(rows),'spacings_m':list(spacings),
        'grid_origins_m':{f'{h:g}':origins[h].tolist() for h in spacings},
        'grid_alignment':'same lattice as data/derived/material-domains/whole-body-<h>m',
        'conforming_reference':{'path':'data/derived/muscle-tet-ready-verification-v1/tetgen-trials.jsonl',
            'sha256':sha(TRIALS),'flags':'pYq1.414',
            'total_tets':int(sum(r['conforming']['tets'] or 0 for r in rows)),
            'median_tets':float(np.median([r['conforming']['tets'] for r in rows if r['conforming']['tets']])),
            'tetrahedralized':int(sum(1 for r in rows if r['conforming']['succeeded'])),
            'max_tet_vs_surface_volume_relative_error':float(max(r['conforming']['tet_volume_relative_error'] or 0 for r in rows))},
        'surface_volume_total_m3':float(sum(r['surface_volume_m3'] for r in rows)),
        'surface_area_total_m2':float(sum(r['surface_area_m2'] for r in rows)),'per_spacing':{}}
    for h in spacings:
        k=f'{h:g}';s=[r['spacings'][k] for r in rows]
        cells=np.array([x['cells'] for x in s]);err=np.array([x['volume_relative_error'] for x in s])
        aerr=np.array([x['area_relative_error'] for x in s])
        dev=[r['spacings'][k]['deviation'] for r in rows if r['deviation_sampled'] and r['spacings'][k].get('deviation')]
        summary['per_spacing'][k]={'spacing_m':h,'cell_volume_m3':h**3,
            'total_voxel_volume_m3':float(cells.sum()*h**3),
            'total_voxel_volume_relative_error':float(cells.sum()*h**3/summary['surface_volume_total_m3']-1),
            'total_voxel_boundary_area_m2':float(sum(x['voxel_boundary_faces'] for x in s)*h**2),
            'total_area_relative_error':float(sum(x['voxel_boundary_faces'] for x in s)*h**2/summary['surface_area_total_m2']-1),
            'total_voxel_tets':int(6*cells.sum()),
            'conforming_tets_per_voxel_tet':float(summary['conforming_reference']['total_tets']/max(6*cells.sum(),1)),
            'muscles_with_no_cell':int((cells==0).sum()),
            'muscles_with_no_cell_ids':[r['entity_id'] for r in rows if r['spacings'][k]['cells']==0],
            'muscles_under_8_cells':int((cells<8).sum()),'muscles_under_64_cells':int((cells<64).sum()),
            'volume_of_muscles_with_no_cell_m3':float(sum(r['surface_volume_m3'] for r in rows if r['spacings'][k]['cells']==0)),
            'volume_error_percentiles':{p:float(np.percentile(err[cells>0],int(p))) for p in ('5','25','50','75','95')},
            'abs_volume_error_over_25pct':int((np.abs(err)>.25).sum()),'abs_volume_error_over_50pct':int((np.abs(err)>.5).sum()),
            'area_error_percentiles':{p:float(np.percentile(aerr[cells>0],int(p))) for p in ('5','25','50','75','95')},
            'deviation_sample':{'muscles':len(dev),
                'median_symmetric_mean_m':float(np.median([d['symmetric_mean_m'] for d in dev])) if dev else None,
                'median_symmetric_hausdorff_m':float(np.median([d['symmetric_hausdorff_m'] for d in dev])) if dev else None,
                'max_symmetric_hausdorff_m':float(max(d['symmetric_hausdorff_m'] for d in dev)) if dev else None,
                'median_hausdorff_over_cell_diagonal':float(np.median([d['hausdorff_over_cell_diagonal'] for d in dev])) if dev else None}}
    summary['wall_seconds']=time.monotonic()-began
    summary['unverified']=['voxel occupancy is cell-center generalized winding over 0.5, not a conservative or volume-preserving voxelization',
        'surface deviation is a sampled bound at h/4 spacing, not a certified Hausdorff distance',
        'the conforming side is the retained pYq1.414 TetGen trial, not a re-run, and its element quality is not re-measured here',
        'no mechanical solve is performed on either representation, so no stress or strain error is measured']
    out.mkdir(parents=True,exist_ok=True)
    with (out/'muscles.jsonl').open('w') as f:
        for r in rows:f.write(json.dumps(r,allow_nan=False)+'\n')
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (out/'manifest.json').write_text(json.dumps({'schema':'ihm.voxel-vs-conforming-muscle-manifest.v1',
        'inputs_sha256':{'data/derived/muscle-tet-ready-v1/manifest.json':sha(MUSCLE/'manifest.json'),
            'data/derived/muscle-tet-ready-verification-v1/tetgen-trials.jsonl':sha(TRIALS),
            'data/derived/canonical/anatomy.json':sha(ROOT/'data/derived/canonical/anatomy.json'),
            'ihm/assembly/material_domains.py':sha(ROOT/'ihm/assembly/material_domains.py'),
            'scripts/measure_voxel_vs_conforming_muscle.py':sha(__file__)},
        'python':sys.version,'igl':getattr(igl,'__version__','2.6.2'),
        'artifacts_sha256':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'},
        'canonical_assets_modified':False},indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ('entities','conforming_reference','surface_volume_total_m3','per_spacing','wall_seconds')},indent=2,allow_nan=False))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--spacings-m',type=float,nargs='*',default=[.01,.005])
    p.add_argument('--sample',type=int,default=36)
    p.add_argument('--out',default=str(ROOT/'data/derived/voxel-vs-conforming-muscle-v1'))
    p.add_argument('--force',action='store_true');a=p.parse_args()
    run(tuple(a.spacings_m),a.sample,a.out,a.force)
