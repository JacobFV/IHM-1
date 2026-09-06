"""Evidence-attached, non-overlapping volume materializations in canonical SI.

Voxel cells discretize source closed surfaces; they are not measured voxels.
Ordered source identities resolve overlap explicitly, without adding two masses.
"""
from pathlib import Path
import gzip, hashlib, itertools, json
import numpy as np

VOXEL_CELL_CAP=16_000_000
GRID_BYTES_PER_CELL=113

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

def _corner_winding(p,tri,chunk=None):
    """Dense oriented solid angle over explicit triangle corner coordinates."""
    if not len(p):return np.empty(0)
    if not len(tri):return np.zeros(len(p))
    step=int(chunk) if chunk else max(1,min(len(p),(1<<17)//len(tri)))  # keep the (step,faces,3) temporaries cache resident
    result=[]
    for start in range(0,len(p),step):
        a=tri[None,:,0]-p[start:start+step,None,:]
        b=tri[None,:,1]-p[start:start+step,None,:]
        c=tri[None,:,2]-p[start:start+step,None,:]
        la=np.linalg.norm(a,axis=2);lb=np.linalg.norm(b,axis=2);lc=np.linalg.norm(c,axis=2)
        determinant=np.einsum('ijk,ijk->ij',a,np.cross(b,c))
        denominator=la*lb*lc+np.einsum('ijk,ijk->ij',a,b)*lc+np.einsum('ijk,ijk->ij',b,c)*la+np.einsum('ijk,ijk->ij',c,a)*lb
        result.append(np.sum(2*np.arctan2(determinant,denominator),axis=1)/(4*np.pi))
    return np.concatenate(result)

def winding_numbers(points,vertices,triangles,chunk=96):
    """Generalized winding by oriented solid angle; no ray-jitter assumption.

    Reference kernel: dense O(points x faces). Retained as the equivalence
    baseline for WindingHierarchy, not as the production evaluator.
    """
    return _corner_winding(np.asarray(points,float),np.asarray(vertices,float)[np.asarray(triangles,int)],chunk)

def oriented_boundary(faces):
    """Oriented edges the face set does not cancel internally; empty when closed."""
    faces=np.asarray(faces,int)
    e=np.concatenate((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
    unique,inverse=np.unique(np.sort(e,axis=1),axis=0,return_inverse=True)
    net=np.rint(np.bincount(inverse,weights=np.where(e[:,0]<e[:,1],1.,-1.),minlength=len(unique))).astype(int)
    keep=np.flatnonzero(net)
    if not len(keep):return np.empty((0,2),int)
    edges=np.repeat(unique[keep],np.abs(net[keep]),axis=0)
    flip=np.repeat(net[keep]<0,np.abs(net[keep]));edges[flip]=edges[flip][:,::-1]
    return edges

class WindingHierarchy:
    """Exact hierarchical solid angle; far subtrees collapse to their boundary cap.

    A face subset and the apex fan over its oriented boundary form a closed
    surface contained in the subset bounding box, because the apex is the box
    centre. Outside that box the subset winding is therefore exactly the negated
    cap winding: the substitution is algebraic, not an error-bounded far-field
    approximation, so classification cannot drift. A closed subset caps to
    nothing, which is the bounding-volume cull.
    """
    def __init__(self,vertices,triangles,leaf=64,margin=1e-9):
        self.x=np.asarray(vertices,float);self.faces=np.asarray(triangles,int)
        self.leaf=int(leaf);self.margin=float(margin);self.nodes=[]
        self._build(np.arange(len(self.faces)))
        self.lower=self.nodes[0]['lower'];self.upper=self.nodes[0]['upper']
    def _build(self,index):
        corners=self.x[self.faces[index]];flat=corners.reshape(-1,3)
        lower=flat.min(axis=0);upper=flat.max(axis=0);center=(lower+upper)/2
        pad=self.margin*max(float((upper-lower).max()),1e-30)
        boundary=oriented_boundary(self.faces[index])
        cap=None
        if len(boundary)<len(index):
            a=self.x[boundary[:,0]];b=self.x[boundary[:,1]]
            cap=np.stack((np.broadcast_to(center,a.shape),b,a),axis=1)
        node={'lower':lower-pad,'upper':upper+pad,'cap':cap,'children':None,'tri':None}
        self.nodes.append(node);slot=len(self.nodes)-1
        if len(index)<=self.leaf:node['tri']=corners;return slot
        mid=corners.mean(axis=1);axis=int(np.argmax(mid.max(axis=0)-mid.min(axis=0)))
        order=index[np.argsort(mid[:,axis],kind='stable')];half=len(order)//2
        node['children']=(self._build(order[:half]),self._build(order[half:]))
        return slot
    def __call__(self,points,chunk=None):
        p=np.asarray(points,float);out=np.zeros(len(p))
        if not len(p):return out
        stack=[(0,np.arange(len(p)))]
        while stack:
            i,sel=stack.pop()
            node=self.nodes[i];q=p[sel]
            if node['cap'] is not None:
                near=np.all((q>=node['lower'])&(q<=node['upper']),axis=1)
                if not near.all():
                    j=sel[~near]
                    if len(node['cap']):out[j]-=_corner_winding(p[j],node['cap'],chunk)
                    sel=sel[near]
                    if not len(sel):continue
            if node['children'] is None:out[sel]+=_corner_winding(p[sel],node['tri'],chunk)
            else:stack.extend((c,sel) for c in node['children'])
        return out

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
    count=float(np.prod(divisions,dtype=float))
    if np.any(divisions<1) or count>VOXEL_CELL_CAP:raise ValueError(
        f'Voxel domain {tuple(int(d) for d in divisions)} = {count:.3g} cells exceeds the {VOXEL_CELL_CAP:.3g} cell cap; '
        f'the grid phase holds about {GRID_BYTES_PER_CELL} B/cell (int64 cell indices, float64 centers, owner, coverage, '
        f'per-source winding and mask), so the cap bounds pre-occupancy grid memory near '
        f'{VOXEL_CELL_CAP*GRID_BYTES_PER_CELL/2**30:.1f} GiB; the requested grid would need {count*GRID_BYTES_PER_CELL/2**30:.3g} GiB')
    cells=np.stack(np.meshgrid(*(np.arange(n) for n in divisions),indexing='ij'),axis=-1).reshape(-1,3)
    centers=origin+(cells+.5)*h
    owners=np.full(len(cells),-1,int);coverage=np.zeros(len(cells),int)
    source_counts=[]
    for i,(name,x,t) in enumerate(surfaces):
        x=np.asarray(x,float);inside=np.zeros(len(cells),bool)
        tree=WindingHierarchy(x,t)
        box=np.all((centers>=tree.lower)&(centers<=tree.upper),axis=1)
        if box.any():inside[box]=np.abs(tree(centers[box]))>.5
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
