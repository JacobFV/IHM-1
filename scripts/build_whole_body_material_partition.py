"""Measured whole-body voxel material partition from every canonical closed-surface source.

Runs the retained ihm.assembly.material_domains.voxel_partition unmodified over all
2404 triangular-surface entities at one spacing; the only wrapper is a timing shim
around WindingHierarchy that calls through without touching any arithmetic.

Measured here: cell-center generalized-winding occupancy per entity, explicit
ordered overlap resolution, shared-node hexahedral six-tet output, per-entity and
per-system owned volume, wall time and peak RSS.

Not measured here: any density, mass activation, tissue interface, contact, or
constitutive property; voxel occupancy is a discretization of source surfaces, not
a segmentation of imaging.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,resource,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import ihm.assembly.material_domains as md
from ihm.assembly.material_domains import source_surface,voxel_partition

MUSCLE=ROOT/'data/derived/muscle-tet-ready-v1'
ENVELOPE=ROOT/'data/derived/outer-envelope'
TIMING={'build_s':0.,'query_s':0.,'trees':0,'nodes':0}

class TimedHierarchy(md.WindingHierarchy):
    """Identical arithmetic; records tree construction and traversal time only."""
    def __init__(self,*a,**k):
        s=time.perf_counter();super().__init__(*a,**k)
        TIMING['build_s']+=time.perf_counter()-s;TIMING['trees']+=1;TIMING['nodes']+=len(self.nodes)
    def __call__(self,*a,**k):
        s=time.perf_counter();r=super().__call__(*a,**k);TIMING['query_s']+=time.perf_counter()-s;return r

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def mesh(path):
    raw=Path(path).read_bytes()
    d=json.loads(gzip.decompress(raw) if str(path).endswith('.gz') else raw)
    return (np.asarray(d.get('positions',d.get('vertices')),float).reshape(-1,3),
            np.asarray(d.get('indices',d.get('faces')),int).reshape(-1,3),hashlib.sha256(raw).hexdigest())

def divergence(V,F):
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]]
    return float(np.einsum('ij,ij->i',A,np.cross(B,C)).sum()/6),float(np.linalg.norm(np.cross(B-A,C-A),axis=1).sum()/2)

def collect():
    """Canonical triangular surfaces, with the repaired muscle surface substituted where one exists."""
    anatomy=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
    mechanics={e['id']:e for e in json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())['entities']}
    frame=anatomy['frame']['id']
    repaired={}
    for line in (MUSCLE/'entities.jsonl').read_text().splitlines():
        r=json.loads(line);repaired[r['entity_id']]=r
    rows=[];skipped={}
    for e in anatomy['entities']:
        ref=e['reference_geometry']
        if ref.get('representation')!='triangular_surface' or ref.get('units')!='m' or ref.get('frame')!=frame:
            skipped[str(ref.get('representation'))]=skipped.get(str(ref.get('representation')),0)+1;continue
        canonical=ROOT/ref['path']
        if sha(canonical)!=ref['sha256']:raise ValueError('Stale canonical geometry for '+e['id'])
        if e['id'] in repaired:
            r=repaired[e['id']]
            if r['source_sha256']!=ref['sha256']:raise ValueError('Repaired muscle built from a different canonical source: '+e['id'])
            path=MUSCLE/r['output_path'];origin='muscle-tet-ready-v1'
            if sha(path)!=r['output_sha256']:raise ValueError('Stale repaired muscle geometry for '+e['id'])
        else:path=canonical;origin='canonical'
        V,F,digest=mesh(path)
        try:
            wx,wt,meta=source_surface(path);closure='closed_after_exact_welding';welded=len(wx)
        except ValueError as error:
            closure=str(error);wx,inverse=np.unique(V,axis=0,return_inverse=True);wt=inverse[F];welded=len(wx)
        vol,area=divergence(V,F)
        rows.append({'entity_id':e['id'],'name':e['name'],'system':e['system'],'role':e.get('role'),
            'geometry_origin':origin,'path':str(path.relative_to(ROOT)),'geometry_sha256':digest,
            'faces':int(len(F)),'source_vertices':int(len(V)),'welded_vertices':int(welded),'closure':closure,
            'signed_divergence_volume_m3':vol,'surface_area_m2':area,
            'allocated_mass_kg':mechanics.get(e['id'],{}).get('mass_kg'),
            'lo':V.min(0).tolist(),'hi':V.max(0).tolist(),'_x':wx,'_t':wt})
    return anatomy,rows,skipped

def build(spacing_m,output_dir=None,force=False,smoke=0):
    out=Path(output_dir or ROOT/f'data/derived/material-domains/whole-body-{spacing_m:g}m')
    if out.exists() and any(out.iterdir()) and not force:raise ValueError('Choose a fresh material-domain directory or pass --force')
    began=time.monotonic();t=time.perf_counter()
    anatomy,rows,skipped=collect();load_s=time.perf_counter()-t
    if smoke:rows=rows[::max(1,len(rows)//smoke)]  # stride subset for a wiring smoke test only
    # Descending absolute divergence volume: voxel_partition gives a contested cell to
    # the later source, so the smallest, most spatially specific claimant owns it.
    order=sorted(range(len(rows)),key=lambda i:(-abs(rows[i]['signed_divergence_volume_m3']),rows[i]['entity_id']))
    rows=[rows[i] for i in order]
    for rank,r in enumerate(rows):r['priority_rank']=rank
    surfaces=[(r['entity_id'],r['_x'],r['_t']) for r in rows]
    md.WindingHierarchy=TimedHierarchy
    t=time.perf_counter();partition=voxel_partition(surfaces,spacing_m=spacing_m);partition_s=time.perf_counter()-t
    md.WindingHierarchy=TimedHierarchy.__mro__[1]
    h=partition['spacing_m'];owner=partition['cell_owner']
    owned=np.bincount(owner,minlength=len(rows))
    for i,r in enumerate(rows):
        r.pop('_x');r.pop('_t')
        r['claimed_cells']=int(partition['source_occupied_cells'][i]);r['owned_cells']=int(owned[i])
        r['claimed_volume_m3']=r['claimed_cells']*h**3;r['owned_volume_m3']=r['owned_cells']*h**3
        r['owned_effective_density_kg_m3']=(r['allocated_mass_kg']/r['owned_volume_m3']
            if r['allocated_mass_kg'] and r['owned_cells'] else None)
    systems={}
    for r in rows:
        s=systems.setdefault(r['system'],{'entities':0,'entities_with_no_cell':0,'faces':0,'claimed_cells':0,
            'owned_cells':0,'abs_divergence_volume_m3':0.,'allocated_mass_kg':0.})
        s['entities']+=1;s['faces']+=r['faces'];s['claimed_cells']+=r['claimed_cells'];s['owned_cells']+=r['owned_cells']
        s['entities_with_no_cell']+=int(r['owned_cells']==0)
        s['abs_divergence_volume_m3']+=abs(r['signed_divergence_volume_m3'])
        s['allocated_mass_kg']+=r['allocated_mass_kg'] or 0.
    for s in systems.values():
        s['owned_volume_m3']=s['owned_cells']*h**3;s['claimed_volume_m3']=s['claimed_cells']*h**3
        s['owned_effective_density_kg_m3']=s['allocated_mass_kg']/s['owned_volume_m3'] if s['owned_cells'] else None
    divisions=np.ceil((np.max([r['hi'] for r in rows],axis=0)-partition['origin_m'])/h-1e-10).astype(int)
    x=partition['vertices_m'];tets=partition['tetrahedra']
    det=np.linalg.det(np.swapaxes(x[tets[:,1:]]-x[tets[:,0,None]],1,2))
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20
    out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'whole-body-domain.npz',vertices_m=x,tetrahedra=tets.astype(np.int32),
        material_index=partition['material_index'].astype(np.int32),cell_indices=partition['cell_indices'].astype(np.int32),
        cell_owner=owner.astype(np.int32),origin_m=partition['origin_m'],spacing_m=np.array(h),
        source_ids=np.array(partition['source_ids']))
    with (out/'sources.jsonl').open('w') as f:
        for r in rows:f.write(json.dumps(r,allow_nan=False)+'\n')
    manifest={'schema':'ihm.whole-body-material-partition.v1','spacing_m':h,
        'frame':anatomy['frame'],'wall_seconds':time.monotonic()-began,'peak_rss_gib':peak,
        'stages_s':{'load_and_weld':load_s,'voxel_partition':partition_s,
                    'winding_tree_build':TIMING['build_s'],'winding_traversal':TIMING['query_s'],
                    'grid_and_tetrahedralization':partition_s-TIMING['build_s']-TIMING['query_s']},
        'winding_trees':TIMING['trees'],'winding_tree_nodes':TIMING['nodes'],
        'sources':{'entities':len(rows),'skipped_by_representation':skipped,
            'geometry_origin':{'canonical':sum(1 for r in rows if r['geometry_origin']=='canonical'),
                               'muscle-tet-ready-v1':sum(1 for r in rows if r['geometry_origin']!='canonical')},
            'closed_after_exact_welding':sum(1 for r in rows if r['closure']=='closed_after_exact_welding'),
            'not_closed_generalized_winding_only':sum(1 for r in rows if r['closure']!='closed_after_exact_welding'),
            'not_closed_ids':[r['entity_id'] for r in rows if r['closure']!='closed_after_exact_welding'],
            'faces':sum(r['faces'] for r in rows),'welded_vertices':sum(r['welded_vertices'] for r in rows),
            'summed_abs_divergence_volume_m3':sum(abs(r['signed_divergence_volume_m3']) for r in rows)},
        'ordering':'descending abs signed divergence volume, ties by entity id; voxel_partition gives a contested cell to the later source, so the smallest claimant owns it',
        'grid':{'origin_m':partition['origin_m'].tolist(),'divisions':[int(d) for d in divisions],
            'grid_cells':int(np.prod(divisions,dtype=float)),
            'occupied_bounding_divisions':[int(d) for d in np.ptp(partition['cell_indices'],axis=0)+1],
            'occupied_cells':int(len(owner)),'cell_volume_m3':h**3,
            'boundary_discretization_diagonal_m':partition['boundary_discretization_diagonal_m']},
        'overlap':{'rule':partition['overlap_rule'],'cells_claimed_by_more_than_one_source':partition['overlapping_source_cells'],
            'summed_claimed_cells':int(sum(r['claimed_cells'] for r in rows)),
            'reassigned_cells':int(sum(r['claimed_cells'] for r in rows))-int(len(owner))},
        'mesh':{'nodes':int(len(x)),'tetrahedra':int(len(tets)),'min_tet_volume_m3':float(det.min()/6),
            'max_tet_volume_m3':float(det.max()/6),'all_tetrahedra_positive':bool((det>0).all()),
            'total_tet_volume_m3':float(det.sum()/6)},
        'entities_claiming_no_cell':int(sum(1 for r in rows if r['claimed_cells']==0)),
        'entities_owning_no_cell':int(sum(1 for r in rows if r['owned_cells']==0)),
        'systems':dict(sorted(systems.items(),key=lambda kv:-kv[1]['owned_cells'])),
        'envelope_reference':{'path':'data/derived/outer-envelope/outer-envelope.npz',
            'sha256':sha(ENVELOPE/'outer-envelope.npz'),
            'divergence_volume_m3':json.loads((ROOT/'data/derived/unmodelled-volume/summary.json').read_text())['envelope']['divergence_volume_m3']},
        'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [
            ROOT/'data/derived/canonical/anatomy.json',ROOT/'data/derived/canonical/mechanics.json',
            MUSCLE/'manifest.json',ENVELOPE/'manifest.json',
            ROOT/'ihm/assembly/material_domains.py',Path(__file__)]},
        'geometry_sha256_count':len({r['geometry_sha256'] for r in rows}),
        'instrumentation':'WindingHierarchy subclassed to record build and traversal time; arithmetic and control flow unchanged',
        'not_verified':['occupancy is cell-center generalized winding over 0.5, not a certified solid or a measured segmentation',
            'the 28 non-closed sources get a generalized-winding interpretation only',
            'no density, mass activation, constitutive law, interface or contact is resolved by this partition',
            'skin and other sub-cell-thickness shells cannot be resolved at this spacing',
            'tetrahedral quality beyond positive volume is not assessed'],
        'canonical_assets_modified':False}
    manifest['artifacts_sha256']={p.name:sha(p) for p in [out/'whole-body-domain.npz',out/'sources.jsonl']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:manifest[k] for k in ('spacing_m','wall_seconds','peak_rss_gib','stages_s','grid','overlap','mesh',
        'entities_claiming_no_cell','entities_owning_no_cell')},indent=2))
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--spacing-m',type=float,default=.01)
    p.add_argument('--output-dir');p.add_argument('--force',action='store_true')
    p.add_argument('--smoke',type=int,default=0);a=p.parse_args()
    build(a.spacing_m,a.output_dir,a.force,a.smoke)
