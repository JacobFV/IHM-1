"""Hierarchical winding reproduces the dense kernel exactly on the retained pelvic domains.

Verified here: cell-center occupancy, node/tetrahedron counts, per-source cell
counts, overlap count and material volumes at 0.004 m and 0.002 m are identical
to the retained manifests; the raised grid cap admits a whole-body 5 mm grid and
still rejects 2 mm with a memory-explicit message; the boundary-cap substitution
is algebraic, so agreement is to floating-point round-off, not a tolerance.

Not verified here: any whole-body partition is actually built (only its cell and
point-face work is projected from cached entity extents), tetrahedral quality,
mass activation, and any libigl path -- libigl 2.6.2 publishes no aarch64 abi3
wheel, so the 3.13 main venv has no igl and none is imported.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.material_domains import (VOXEL_CELL_CAP,GRID_BYTES_PER_CELL,WindingHierarchy,
    oriented_boundary,source_surface,voxel_partition)

SOURCES=['body-bp3d-FJ3132','body-bp3d-FJ3133','body-bp3d-FJ3134']
RETAINED={.004:'data/derived/material-domains/pelvis-0.004m',.002:'data/derived/material-domains/pelvis-0.002m'}
OUT=ROOT/'data/derived/material-domains/winding-kernel'

def grid(surfaces,h):
    allx=np.concatenate([np.asarray(s[1],float) for s in surfaces])
    origin=np.floor(allx.min(axis=0)/h+1e-10)*h
    divisions=np.ceil((allx.max(axis=0)-origin)/h-1e-10).astype(int)
    cells=np.stack(np.meshgrid(*(np.arange(n) for n in divisions),indexing='ij'),axis=-1).reshape(-1,3)
    return origin+(cells+.5)*h

def dense(points,x,t,chunk=96):
    """Byte-for-byte the pre-change kernel; not imported, so the baseline cannot drift."""
    p=np.asarray(points,float);tri=np.asarray(x,float)[np.asarray(t,int)];result=[]
    for start in range(0,len(p),chunk):
        a=tri[None,:,0]-p[start:start+chunk,None,:]
        b=tri[None,:,1]-p[start:start+chunk,None,:]
        c=tri[None,:,2]-p[start:start+chunk,None,:]
        la=np.linalg.norm(a,axis=2);lb=np.linalg.norm(b,axis=2);lc=np.linalg.norm(c,axis=2)
        determinant=np.einsum('ijk,ijk->ij',a,np.cross(b,c))
        denominator=la*lb*lc+np.einsum('ijk,ijk->ij',a,b)*lc+np.einsum('ijk,ijk->ij',b,c)*la+np.einsum('ijk,ijk->ij',c,a)*lb
        result.append(np.sum(2*np.arctan2(determinant,denominator),axis=1)/(4*np.pi))
    return np.concatenate(result) if result else np.empty(0)

def entity_extents():
    cache=OUT/'entity-extents.json'
    if cache.exists():return json.loads(cache.read_text())
    rows=[]
    for e in json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())['entities']:
        p=ROOT/e['reference_geometry']['path']
        if not p.exists():continue
        g=json.load(gzip.open(p,'rt') if p.suffix=='.gz' else p.open())
        f=g.get('indices',g.get('faces'));v=g.get('positions',g.get('vertices'))
        if f is None or v is None:continue
        v=np.asarray(v,float).reshape(-1,3)
        rows.append({'id':e['id'],'faces':len(np.asarray(f,int).reshape(-1,3)),'lo':v.min(axis=0).tolist(),'hi':v.max(axis=0).tolist()})
    OUT.mkdir(parents=True,exist_ok=True);cache.write_text(json.dumps(rows,separators=(',',':')))
    return rows

def projection(rows,h):
    lo=np.min([r['lo'] for r in rows],axis=0);hi=np.max([r['hi'] for r in rows],axis=0)
    divisions=np.ceil((hi-lo)/h).astype(int);total=int(np.prod(divisions,dtype=float));faces=sum(r['faces'] for r in rows)
    culled=sum(int(np.prod(np.maximum(np.ceil((np.array(r['hi'])-np.array(r['lo']))/h),1)))*r['faces'] for r in rows)
    return {'spacing_m':h,'divisions':[int(d) for d in divisions],'cells':total,'entities':len(rows),'faces':faces,
            'dense_point_faces':float(total)*faces,'aabb_culled_point_faces':float(culled),
            'aabb_cull_factor':float(total)*faces/culled,'within_cell_cap':total<=VOXEL_CELL_CAP}

def scale_benchmark(entity_id,h,sample,seed=0):
    """Full hierarchy over a whole-body grid; dense only on a random subsample."""
    by={e['id']:e for e in json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())['entities']}
    x,t,_=source_surface(ROOT/by[entity_id]['reference_geometry']['path'])
    rows=entity_extents();lo=np.min([r['lo'] for r in rows],axis=0);hi=np.max([r['hi'] for r in rows],axis=0)
    origin=np.floor(lo/h)*h;divisions=np.ceil((hi-origin)/h).astype(int)
    cells=np.stack(np.meshgrid(*(np.arange(n) for n in divisions),indexing='ij'),axis=-1).reshape(-1,3)
    centers=origin+(cells+.5)*h
    start=time.perf_counter();tree=WindingHierarchy(x,t);build=time.perf_counter()-start
    box=np.all((centers>=tree.lower)&(centers<=tree.upper),axis=1);points=centers[box]
    start=time.perf_counter();fast=tree(points);hierarchy=time.perf_counter()-start
    index=np.random.default_rng(seed).choice(len(points),min(sample,len(points)),replace=False)
    start=time.perf_counter();reference=dense(points[index],x,t);baseline=time.perf_counter()-start
    assert np.all((np.abs(reference)>.5)==(np.abs(fast[index])>.5))
    projected=baseline/len(index)*len(points)
    return {'entity_id':entity_id,'spacing_m':h,'faces':len(t),'grid_cells':int(len(centers)),'evaluated_points':int(len(points)),
            'tree_nodes':len(tree.nodes),'build_s':build,'hierarchy_s':hierarchy,'occupied_cells':int((np.abs(fast)>.5).sum()),
            'dense_sample_points':int(len(index)),'dense_sample_s':baseline,'dense_point_faces_per_s':len(index)*len(t)/baseline,
            'dense_projected_full_s':projected,'sampled_max_winding_deviation':float(np.abs(reference-fast[index]).max()),
            'speedup':projected/(build+hierarchy),
            'scope':'dense side is a sampled rate extrapolated to the same point set, not a full dense run'}

def main(scale_spacing_m,scale_sample):
    anatomy=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text());by={e['id']:e for e in anatomy['entities']}
    surfaces=[]
    for source_id in SOURCES:
        path=ROOT/by[source_id]['reference_geometry']['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==by[source_id]['reference_geometry']['sha256']
        x,t,_=source_surface(path);surfaces.append((source_id,x,t))
        assert len(oriented_boundary(t))==0
    domains=[]
    for h,retained in RETAINED.items():
        manifest=json.loads((ROOT/retained/'manifest.json').read_text())
        centers=grid(surfaces,h);deviation=0.;reference_s=hierarchy_s=0.;cap_only=0
        for _,x,t in surfaces:
            start=time.perf_counter();slow=dense(centers,x,t);reference_s+=time.perf_counter()-start
            start=time.perf_counter()
            tree=WindingHierarchy(x,t);box=np.all((centers>=tree.lower)&(centers<=tree.upper),axis=1)
            fast=np.zeros(len(centers));fast[box]=tree(centers[box]);hierarchy_s+=time.perf_counter()-start
            assert np.all((np.abs(slow)>.5)==(np.abs(fast)>.5)),'occupancy classification diverged'
            deviation=max(deviation,float(np.abs(slow-fast).max()));cap_only+=int((~box).sum())
        partition=voxel_partition(surfaces,spacing_m=h)
        volumes=np.bincount(partition['material_index'],minlength=len(SOURCES))*h**3/6
        retained_volumes=[r['volume_m3'] for r in manifest['material_regions']]
        assert len(partition['vertices_m'])==manifest['vertices'],(len(partition['vertices_m']),manifest['vertices'])
        assert len(partition['tetrahedra'])==manifest['tetrahedra']
        assert partition['overlapping_source_cells']==manifest['overlapping_source_cells']
        assert np.allclose(volumes,retained_volumes,rtol=0,atol=0)
        domains.append({'spacing_m':h,'retained':retained,'grid_cells':int(len(centers)),
            'aabb_rejected_point_source_pairs':cap_only,'max_winding_deviation':deviation,
            'occupancy_identical':True,'vertices':int(len(partition['vertices_m'])),'tetrahedra':int(len(partition['tetrahedra'])),
            'source_occupied_cells':partition['source_occupied_cells'],'overlapping_source_cells':partition['overlapping_source_cells'],
            'volume_ml':float(volumes.sum()*1e6),'retained_volume_ml':float(sum(retained_volumes)*1e6),
            'dense_reference_s':reference_s,'hierarchy_s':hierarchy_s,'speedup':reference_s/hierarchy_s})
    rows=entity_extents();projections=[projection(rows,h) for h in (.01,.005,.002)]
    assert projections[1]['within_cell_cap'] and not projections[2]['within_cell_cap']
    rejected=None
    try:voxel_partition([('probe',np.array([[0.,0,0],[1e6,1e6,1e6]]),np.zeros((0,3),int))],spacing_m=.002)
    except ValueError as error:rejected=str(error)
    assert rejected and 'GiB' in rejected and 'cell cap' in rejected
    scale=scale_benchmark('body-skin-epidermis',scale_spacing_m,scale_sample)
    igl_status={'main_venv_python':'.'.join(map(str,sys.version_info[:3])),'igl_importable':False,'reason':None}
    try:
        import igl;igl_status['igl_importable']=True;igl_status['reason']='present'
    except ImportError as error:igl_status['reason']=str(error)
    receipt={'schema_version':1,'cell_cap':VOXEL_CELL_CAP,'grid_bytes_per_cell':GRID_BYTES_PER_CELL,
        'cell_cap_memory_gib':VOXEL_CELL_CAP*GRID_BYTES_PER_CELL/2**30,'cap_rejection_message':rejected,
        'kernel':'per-entity AABB cull + exact boundary-cap AABB hierarchy (leaf 64)',
        'exactness_basis':'a face subset plus the apex fan over its oriented boundary is a closed surface inside the subset bounding box, so outside that box the subset winding is exactly the negated cap winding',
        'retained_domains':domains,'whole_body_projection':projections,'scale_benchmark':scale,'libigl':igl_status,
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in [ROOT/'ihm/assembly/material_domains.py',Path(__file__)]},
        'not_verified':['no whole-body partition was built; whole-body figures are projected cell and point-face counts',
                        'tetrahedral quality, mass activation and collision interfaces unchanged and unverified',
                        'no libigl/CGAL fast winding number path exists on this machine and none is exercised']}
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(receipt,indent=2,allow_nan=False))
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scale-spacing-m',type=float,default=.01);p.add_argument('--scale-sample',type=int,default=1500)
    a=p.parse_args();main(a.scale_spacing_m,a.scale_sample)
