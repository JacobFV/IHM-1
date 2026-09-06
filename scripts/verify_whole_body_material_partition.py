"""Independent validation of a built whole-body voxel material partition.

Verified here, re-derived from the stored arrays rather than from the builder's
own state: every tetrahedron has positive volume and the exact cell volume; the
node set is exactly the shared corner lattice of the occupied cells, so touching
cells share nodes; no cell index appears twice and every owner is in range; the
partition's ownership and vacancy decisions reproduce under the dense reference
winding kernel on a random sample; occupied volume against the outer envelope
interior measured on the same grid; per-entity claimed volume against that
entity's own surface divergence integral.

Not verified here: tetrahedral quality beyond positivity, any density, mass,
constitutive or interface property, and whether an occupied cell is anatomically
that tissue -- occupancy is cell-center generalized winding over 0.5.
"""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.material_domains import WindingHierarchy,winding_numbers,source_surface
from scripts.build_whole_body_material_partition import mesh

CORNERS=np.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]])
ENVELOPE=ROOT/'data/derived/outer-envelope'
KEY=lambda a:(a[:,0]+(1<<20))*(1<<42)+(a[:,1]+(1<<20))*(1<<21)+(a[:,2]+(1<<20))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def mesh_checks(d):
    x=d['vertices_m'];tets=d['tetrahedra'].astype(np.int64);cells=d['cell_indices'].astype(np.int64)
    owner=d['cell_owner'].astype(np.int64);h=float(d['spacing_m']);origin=d['origin_m']
    det=np.linalg.det(np.swapaxes(x[tets[:,1:]]-x[tets[:,0,None]],1,2));vol=det/6
    lattice=np.rint((x-origin)/h).astype(np.int64)
    expected=np.unique((cells[:,None,:]+CORNERS).reshape(-1,3),axis=0)
    key=lambda a:a[np.lexsort(a.T[::-1])]
    node_ok=bool(len(lattice)==len(expected) and np.array_equal(key(lattice),key(expected))
                 and np.allclose(x,origin+lattice*h,rtol=0,atol=1e-12))
    _,tet_incidence=np.unique(tets.reshape(-1),return_counts=True)
    _,cell_incidence=np.unique((cells[:,None,:]+CORNERS).reshape(-1,3),axis=0,return_counts=True)
    return {'nodes':int(len(x)),'tetrahedra':int(len(tets)),'occupied_cells':int(len(cells)),
        'all_tetrahedra_positive':bool((vol>0).all()),'min_tet_volume_m3':float(vol.min()),'max_tet_volume_m3':float(vol.max()),
        'exact_cell_volume_m3':h**3/6,'max_tet_volume_deviation_m3':float(np.abs(vol-h**3/6).max()),
        'tets_per_cell':float(len(tets)/len(cells)),'total_tet_volume_m3':float(vol.sum()),
        'occupied_cell_volume_m3':float(len(cells)*h**3),
        'tet_volume_matches_cell_volume':bool(abs(vol.sum()-len(cells)*h**3)<1e-9*len(cells)*h**3),
        'node_set_is_exact_shared_corner_lattice':node_ok,
        'nodes_referenced_by_every_tet_in_range':bool(tets.min()>=0 and tets.max()<len(x)),
        'unreferenced_nodes':int(len(x)-len(np.unique(tets))),
        'nodes_shared_by_more_than_one_cell':int((cell_incidence>1).sum()),
        'max_cells_per_node':int(cell_incidence.max()),'max_tet_incidences_per_node':int(tet_incidence.max()),
        'cell_indices_unique':bool(len(np.unique(KEY(cells)))==len(cells)),
        'owner_in_range':bool(owner.min()>=0 and owner.max()<len(d['source_ids'])),
        'material_index_matches_cell_owner':bool(np.array_equal(d['material_index'].astype(np.int64),np.repeat(owner,6)))}

def envelope_check(d,chunk):
    """Interior of the retained outer envelope measured on this partition's own grid."""
    h=float(d['spacing_m']);origin=d['origin_m'];cells=d['cell_indices'].astype(np.int64)
    env=np.load(ENVELOPE/'outer-envelope.npz')
    V=np.asarray(env['positions'],float);F=np.asarray(env['indices'],np.int64)
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]]
    envelope_volume=float(np.einsum('ij,ij->i',A,np.cross(B,C)).sum()/6)
    lo=np.floor((V.min(0)-origin)/h).astype(int);hi=np.ceil((V.max(0)-origin)/h).astype(int)
    grid=np.stack(np.meshgrid(*(np.arange(lo[i],hi[i]+1) for i in range(3)),indexing='ij'),-1).reshape(-1,3)
    centers=origin+(grid+.5)*h
    t=time.perf_counter();tree=WindingHierarchy(V,F)
    box=np.all((centers>=tree.lower)&(centers<=tree.upper),axis=1)
    w=np.zeros(len(centers));w[box]=tree(centers[box],chunk);interior=np.abs(w)>.5
    seconds=time.perf_counter()-t
    interior_key=KEY(grid[interior])
    inside=np.isin(KEY(cells),interior_key)
    return interior_key,{'envelope_faces':int(len(F)),'envelope_divergence_volume_m3':envelope_volume,
        'envelope_sha256':sha(ENVELOPE/'outer-envelope.npz'),'winding_seconds':seconds,
        'grid_cells_tested':int(len(grid)),'interior_cells':int(interior.sum()),
        'interior_volume_m3':float(interior.sum()*h**3),
        'interior_volume_vs_divergence_ratio':float(interior.sum()*h**3/envelope_volume),
        'occupied_cells_inside_envelope':int(inside.sum()),
        'occupied_volume_inside_envelope_m3':float(inside.sum()*h**3),
        'occupied_cells_outside_envelope':int((~inside).sum()),
        'occupied_volume_outside_envelope_m3':float((~inside).sum()*h**3),
        'fill_fraction_of_interior':float(inside.sum()/max(interior.sum(),1)),
        'void_fraction_of_interior':float(1-inside.sum()/max(interior.sum(),1))}

def dense_check(d,rows,owned_sample,vacant_sample,seed,interior_key):
    """Reproduce ownership and vacancy at sampled cells with the dense O(points x faces) kernel."""
    h=float(d['spacing_m']);origin=d['origin_m'];cells=d['cell_indices'].astype(np.int64)
    owner=d['cell_owner'].astype(np.int64);ids=[str(s) for s in d['source_ids']]
    rank={r['entity_id']:i for i,r in enumerate(rows)}
    lo=np.array([r['lo'] for r in rows]);hi=np.array([r['hi'] for r in rows])
    rng=np.random.default_rng(seed)
    pick=rng.choice(len(cells),min(owned_sample,len(cells)),replace=False)
    occupied=set(KEY(cells).tolist())
    free=np.array(sorted(set(interior_key.tolist())-occupied))
    take=free[rng.choice(len(free),min(vacant_sample,len(free)),replace=False)]
    vacant=np.stack(((take//(1<<42))-(1<<20),((take//(1<<21))%(1<<21))-(1<<20),(take%(1<<21))-(1<<20)),axis=1)
    cache={}
    def surface(i):
        if i not in cache:
            r=rows[i]
            try:x,t,_=source_surface(ROOT/r['path'])
            except ValueError:
                V,F,_=mesh(ROOT/r['path']);x,inverse=np.unique(V,axis=0,return_inverse=True);t=inverse[F]
            cache[i]=(x,t)
        return cache[i]
    owner_confirmed=0;owner_failed=[];priority_confirmed=0;priority_failed=[];pairs=0
    t0=time.perf_counter()
    for j in pick:
        p=(origin+(cells[j]+.5)*h)[None,:]
        candidates=np.flatnonzero(np.all((p>=lo-1e-12)&(p<=hi+1e-12),axis=1)).tolist()
        own=rank[ids[owner[j]]];claims=[]
        for i in candidates:
            x,t=surface(i);pairs+=len(t)
            if abs(float(winding_numbers(p,x,t)[0]))>.5:claims.append(i)
        if own in claims:owner_confirmed+=1
        else:owner_failed.append({'cell':cells[j].tolist(),'owner':ids[owner[j]],'dense_claims':[rows[i]['entity_id'] for i in claims]})
        if claims and max(claims)==own:priority_confirmed+=1
        elif claims:priority_failed.append({'cell':cells[j].tolist(),'owner':ids[owner[j]],'last_dense_claimant':rows[max(claims)]['entity_id']})
    vacancy_confirmed=0;vacancy_failed=[]
    for row in vacant:
        p=(origin+(row+.5)*h)[None,:]
        candidates=np.flatnonzero(np.all((p>=lo-1e-12)&(p<=hi+1e-12),axis=1)).tolist()
        claims=[]
        for i in candidates:
            x,t=surface(i);pairs+=len(t)
            if abs(float(winding_numbers(p,x,t)[0]))>.5:claims.append(rows[i]['entity_id'])
        if claims:vacancy_failed.append({'cell':row.tolist(),'dense_claims':claims})
        else:vacancy_confirmed+=1
    return {'sampled_owned_cells':int(len(pick)),'owner_reproduced_by_dense_kernel':owner_confirmed,
        'owner_failures':owner_failed[:10],'ordering_rule_reproduced':priority_confirmed,
        'ordering_failures':priority_failed[:10],
        'sampled_vacant_cells':int(len(vacant)),'vacancy_reproduced_by_dense_kernel':vacancy_confirmed,
        'vacancy_failures':vacancy_failed[:10],'dense_point_face_pairs':int(pairs),
        'dense_seconds':time.perf_counter()-t0,
        'scope':'dense kernel evaluated at sampled cell centers only, over every source whose bounding box contains the point'}

def entity_check(rows,h):
    claimed=np.array([r['claimed_volume_m3'] for r in rows])
    divergence=np.array([abs(r['signed_divergence_volume_m3']) for r in rows])
    nonzero=(claimed>0)&(divergence>0)
    ratio=claimed[nonzero]/divergence[nonzero]
    order=np.argsort(-np.abs(claimed-divergence))
    return {'entities':len(rows),'entities_claiming_no_cell':int((claimed==0).sum()),
        'divergence_volume_of_entities_claiming_no_cell_m3':float(divergence[claimed==0].sum()),
        'largest_divergence_volume_claiming_no_cell_m3':float(divergence[claimed==0].max()) if (claimed==0).any() else 0.,
        'summed_claimed_volume_m3':float(claimed.sum()),'summed_abs_divergence_volume_m3':float(divergence.sum()),
        'summed_ratio':float(claimed.sum()/divergence.sum()),
        'ratio_percentiles':{p:float(np.percentile(ratio,int(p))) for p in ('5','25','50','75','95')},
        'entities_within_20_percent_of_divergence':int((np.abs(ratio-1)<=.2).sum()),
        'entities_within_50_percent_of_divergence':int((np.abs(ratio-1)<=.5).sum()),
        'entities_compared':int(nonzero.sum()),
        'worst_absolute_deviations':[{'entity_id':rows[i]['entity_id'],'name':rows[i]['name'],'system':rows[i]['system'],
            'faces':rows[i]['faces'],'closure':'closed' if rows[i]['closure']=='closed_after_exact_welding' else 'open',
            'divergence_volume_m3':float(divergence[i]),'claimed_volume_m3':float(claimed[i]),
            'owned_volume_m3':rows[i]['owned_volume_m3']} for i in order[:15]],
        'basis':'claimed volume is the entity own winding occupancy before overlap resolution; owned volume is after'}

def external_cross_check(rows):
    """Independent implementation: libigl fast_winding_number on an 8 mm node grid, from the retained void inventory."""
    path=ROOT/'data/derived/unmodelled-volume/entity_occupancy.jsonl'
    other={};
    for line in path.read_text().splitlines():
        r=json.loads(line);other[r['entity_id']]=r['winding_voxel_volume_m3']
    ids=[r['entity_id'] for r in rows if r['entity_id'] in other]
    a=np.array([next(r for r in rows if r['entity_id']==i)['claimed_volume_m3'] for i in ids])
    b=np.array([other[i] for i in ids])
    m=(a>0)&(b>0);ratio=a[m]/b[m]
    return {'reference':'data/derived/unmodelled-volume/entity_occupancy.jsonl',
        'reference_sha256':sha(path),'reference_kernel':'libigl fast_winding_number, 8 mm node-centred grid, canonical muscle geometry',
        'entities_compared':int(len(ids)),'this_summed_claimed_volume_m3':float(a.sum()),
        'reference_summed_volume_m3':float(b.sum()),'summed_ratio':float(a.sum()/b.sum()),
        'both_nonzero':int(m.sum()),'median_ratio':float(np.median(ratio)),
        'ratio_p05':float(np.percentile(ratio,5)),'ratio_p95':float(np.percentile(ratio,95)),
        'zero_here_nonzero_there':int(((a==0)&(b>0)).sum()),'nonzero_here_zero_there':int(((a>0)&(b==0)).sum()),
        'scope':'different winding implementation, spacing and grid alignment, and canonical rather than repaired muscle surfaces, so per-entity agreement is expected only to discretization'}

def run(domain,out,owned_sample,vacant_sample,chunk,seed):
    domain=Path(domain);out=Path(out)
    manifest=json.loads((domain/'manifest.json').read_text())
    rows=[json.loads(l) for l in (domain/'sources.jsonl').read_text().splitlines()]
    d={k:v for k,v in np.load(domain/'whole-body-domain.npz',allow_pickle=False).items()}
    h=float(d['spacing_m'])
    began=time.monotonic()
    m=mesh_checks(d);interior_key,env=envelope_check(d,chunk)
    dense=dense_check(d,rows,owned_sample,vacant_sample,seed,interior_key)
    ent=entity_check(rows,h);cross=external_cross_check(rows)
    owned=np.bincount(d['cell_owner'].astype(np.int64),minlength=len(rows))
    receipt={'schema':'ihm.whole-body-material-partition-verification.v1','domain':str(domain.relative_to(ROOT)),
        'spacing_m':h,'wall_seconds':time.monotonic()-began,
        'mesh':m,'envelope':env,'dense_reference_sample':dense,'per_entity_volume':ent,'external_cross_check':cross,
        'ownership_consistency':{'summed_owned_volume_m3':float(owned.sum()*h**3),
            'owned_equals_occupied':bool(int(owned.sum())==m['occupied_cells']),
            'owned_volume_vs_envelope_divergence':float(owned.sum()*h**3/env['envelope_divergence_volume_m3'])},
        'inputs_sha256':{str((domain/n).relative_to(ROOT)):sha(domain/n) for n in
            ('manifest.json','sources.jsonl','whole-body-domain.npz')}|
            {'ihm/assembly/material_domains.py':sha(ROOT/'ihm/assembly/material_domains.py'),
             'scripts/build_whole_body_material_partition.py':sha(ROOT/'scripts/build_whole_body_material_partition.py'),
             str(Path(__file__).relative_to(ROOT)):sha(__file__)},
        'build_manifest_agreement':{k:bool(manifest['mesh'][k]==m[k] if k in m else False)
            for k in ('nodes','tetrahedra','all_tetrahedra_positive')},
        'not_verified':['tetrahedral quality beyond positive volume','density, mass activation, constitutive law, interfaces, contact',
            'anatomical correctness of any occupancy decision','sub-cell-thickness shells such as skin at this spacing']}
    passed=(m['all_tetrahedra_positive'] and m['node_set_is_exact_shared_corner_lattice'] and m['cell_indices_unique']
        and m['owner_in_range'] and m['material_index_matches_cell_owner'] and m['tet_volume_matches_cell_volume']
        and not dense['owner_failures'] and not dense['ordering_failures'] and not dense['vacancy_failures']
        and receipt['ownership_consistency']['owned_equals_occupied'])
    receipt['all_structural_checks_passed']=bool(passed)
    out.mkdir(parents=True,exist_ok=True)
    name=f'receipt-{h:g}m.json';(out/name).write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(receipt,indent=2,allow_nan=False))
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--domain',default='data/derived/material-domains/whole-body-0.01m')
    p.add_argument('--out',default='data/derived/material-domains/whole-body-verification')
    p.add_argument('--owned-sample',type=int,default=120);p.add_argument('--vacant-sample',type=int,default=120)
    p.add_argument('--chunk',type=int,default=None);p.add_argument('--seed',type=int,default=0);a=p.parse_args()
    run(ROOT/a.domain,ROOT/a.out,a.owned_sample,a.vacant_sample,a.chunk,a.seed)
