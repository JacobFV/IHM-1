"""Evidence-attached, non-overlapping volume materializations in canonical SI.

Voxel cells discretize source closed surfaces; they are not measured voxels.
Ordered source identities resolve overlap explicitly, without adding two masses.
"""
from pathlib import Path
import gzip, hashlib, itertools, json
import numpy as np

def source_surface(path):
    """Weld only byte-equal vertex coordinates, preserving all source positions."""
    path=Path(path)
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else path.open()) as f:g=json.load(f)
    raw=np.asarray(g.get('positions',g.get('vertices')),float).reshape(-1,3)
    faces=np.asarray(g.get('indices',g.get('faces')),int).reshape(-1,3)
    x,reverse=np.unique(raw,axis=0,return_inverse=True);faces=reverse[faces]
    edges=np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
    _,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
    if np.any(counts!=2):raise ValueError('Source surface is not closed after exact-coordinate welding')
    ordered=np.sort(edges,axis=1)
    _,inv=np.unique(ordered,axis=0,return_inverse=True)
    balance=np.bincount(inv,weights=np.where(edges[:,0]<edges[:,1],1,-1))
    if np.any(balance):raise ValueError('Source surface winding is inconsistent')
    return x,faces,{'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                    'input_vertices':len(raw),'welded_vertices':len(x),'faces':len(faces),
                    'welding':'exact coordinate equality only; no moved coordinates or hole filling',
                    'watertight_after_welding':True}

def winding_numbers(points,vertices,triangles,chunk=96):
    """Generalized winding by oriented solid angle; no ray-jitter assumption."""
    p=np.asarray(points,float);tri=np.asarray(vertices,float)[np.asarray(triangles,int)]
    result=[]
    for start in range(0,len(p),chunk):
        a=tri[None,:,0]-p[start:start+chunk,None,:]
        b=tri[None,:,1]-p[start:start+chunk,None,:]
        c=tri[None,:,2]-p[start:start+chunk,None,:]
        la=np.linalg.norm(a,axis=2);lb=np.linalg.norm(b,axis=2);lc=np.linalg.norm(c,axis=2)
        determinant=np.einsum('ijk,ijk->ij',a,np.cross(b,c))
        denominator=la*lb*lc+np.einsum('ijk,ijk->ij',a,b)*lc+np.einsum('ijk,ijk->ij',b,c)*la+np.einsum('ijk,ijk->ij',c,a)*lb
        result.append(np.sum(2*np.arctan2(determinant,denominator),axis=1)/(4*np.pi))
    return np.concatenate(result) if result else np.empty(0)

def voxel_partition(surfaces,*,spacing_m):
    """Use cell-center winding occupancy; later sources win intersecting cells.

    Every occupied cube is split by the same six-tet pattern. Shared interfaces
    have shared nodes. Boundary position resolution is limited by the cell
    diagonal, irrespective of how many output triangles are generated.
    """
    h=float(spacing_m)
    if isinstance(spacing_m,bool) or not np.isfinite(h) or h<=0:raise ValueError('Positive finite voxel spacing required')
    if not surfaces or len({s[0] for s in surfaces})!=len(surfaces):raise ValueError('Unique nonempty source surfaces required')
    allx=np.concatenate([np.asarray(s[1],float) for s in surfaces])
    origin=np.floor(allx.min(axis=0)/h+1e-10)*h
    divisions=np.ceil((allx.max(axis=0)-origin)/h-1e-10).astype(int)
    if np.any(divisions<1) or np.prod(divisions)>2_000_000:raise ValueError('Unsupported voxel domain extent')
    cells=np.stack(np.meshgrid(*(np.arange(n) for n in divisions),indexing='ij'),axis=-1).reshape(-1,3)
    centers=origin+(cells+.5)*h
    owners=np.full(len(cells),-1,int);coverage=np.zeros(len(cells),int)
    source_counts=[]
    for i,(name,x,t) in enumerate(surfaces):
        winding=winding_numbers(centers,x,t);inside=np.abs(winding)>.5
        coverage+=inside;owners[inside]=i;source_counts.append(int(inside.sum()))
    occupied=owners>=0;cells=cells[occupied];owners=owners[occupied]
    if not len(cells):raise ValueError('No resolved material cells; refine the grid')
    corners=np.array(list(itertools.product((0,1),repeat=3)))
    raw=(cells[:,None,:]+corners).reshape(-1,3)
    unique,index=np.unique(raw,axis=0,return_inverse=True);vertex_map=index.reshape(-1,8)
    pattern=[]
    for order in itertools.permutations(range(3)):
        v=np.zeros(3,int);tet=[0]
        for axis in order:v=v.copy();v[axis]+=1;tet.append(int(np.flatnonzero(np.all(corners==v,axis=1))[0]))
        pattern.append(tet)
    tetrahedra=vertex_map[:,np.array(pattern)].reshape(-1,4)
    vertices=origin+unique*h
    det=np.linalg.det(np.swapaxes(vertices[tetrahedra[:,1:]]-vertices[tetrahedra[:,0,None]],1,2))
    negative=det<0;tetrahedra[negative]=tetrahedra[negative][:,[0,2,1,3]]
    return {'vertices_m':vertices,'tetrahedra':tetrahedra,'material_index':np.repeat(owners,6),
            'cell_indices':cells,'cell_owner':owners,'origin_m':origin,'spacing_m':h,
            'source_ids':[s[0] for s in surfaces],'source_occupied_cells':source_counts,
            'overlapping_source_cells':int(np.sum(coverage>1)),
            'overlap_rule':'later source in explicit ordered list owns the entire intersecting cell',
            'boundary_discretization_diagonal_m':float(np.sqrt(3)*h)}

class MaterialOwnership:
    """Exclusive allocation from an existing body ledger, never additive mass."""
    def __init__(self,source_mass_kg):
        self.source_mass_kg=dict(source_mass_kg);self.claims={};self.source_owner={}
        if any(isinstance(v,bool) or not np.isfinite(v) or v<=0 for v in self.source_mass_kg.values()):raise ValueError('Material ledger masses must be positive finite values')
    def claim(self,domain_id,source_ids):
        ids=list(source_ids)
        if domain_id in self.claims or len(set(ids))!=len(ids) or not ids:raise ValueError('Unique nonempty domain claim required')
        if any(s not in self.source_mass_kg or s in self.source_owner for s in ids):raise ValueError('Unknown or already owned source material')
        self.claims[domain_id]=ids
        self.source_owner.update({s:domain_id for s in ids})
    def mass_kg(self,domain_id):return sum(self.source_mass_kg[s] for s in self.claims[domain_id])
