"""Pair discovery, tied/contact partition and self collision for nodal contact.

`ihm.assembly.contact_dynamics` resolves a declared node/triangle pair and
audits its impulses. It has no pair search: every interface must be named by
hand, one direction at a time. This module supplies the missing spatial layer
without touching that kernel, so its conservation receipts are unchanged.

Three separate things live here and are deliberately not conflated:

1. Broad phase. Bounding-volume pair discovery (sweep and prune, or a BVH, or
   the dense O(N^2) matrix) followed by a surface-cell refinement that is the
   same 4 mm co-occupancy test the coupling assessment used to measure the
   adjacency ground truth. Stage one is a superset filter and its recall is
   measured, not assumed; stage two decides.
2. Interface partition. An anatomically bound interface is a CONSTRAINT, not a
   contact. `partition_interfaces` splits discovered pairs by the recorded
   classification and refuses to hand a tied pair to the contact kernel.
3. Self collision. One deformable body against itself, with topological
   neighbours excluded by shared-vertex ring, since a triangle always touches
   the triangles it shares an edge with.

Every geometric statement is a proxy on an occupancy grid or on axis-aligned
bounds, exactly as the assessment says. It selects candidate interfaces; the
narrow phase still decides whether anything is touching.
"""
from dataclasses import dataclass, field
import numpy as np

__all__=['surface_samples','voxel_keys','encode_cells','aabb_bounds','dense_aabb_pairs',
         'sweep_and_prune_pairs','bvh_pairs','BodySurface','SurfaceCellIndex','BroadPhase',
         'BroadPhaseResult','partition_interfaces','InterfacePartition','self_collision_candidates',
         'SymmetricContact','expand_symmetric','merge_symmetric_receipts']

_PACK_LOW=-512
_PACK_SPAN=1<<20


def _positive(value,name):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value) or not np.isfinite(value) or value<=0:
        raise ValueError(f'Positive finite {name} required')
    return float(value)


def _nonnegative(value,name):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value) or not np.isfinite(value) or value<0:
        raise ValueError(f'Nonnegative finite {name} required')
    return float(value)


def encode_cells(index):
    """Pack a signed voxel index triple into one int64, as the assessment packs it."""
    index=np.asarray(index,np.int64)
    if index.ndim!=2 or index.shape[1]!=3:raise ValueError('Cell index triples required')
    if len(index) and (index.min()<_PACK_LOW or index.max()>=_PACK_SPAN+_PACK_LOW):
        raise ValueError('Voxel index outside packing range')
    return ((index[:,0]-_PACK_LOW)<<40)|((index[:,1]-_PACK_LOW)<<20)|(index[:,2]-_PACK_LOW)


def surface_samples(triangles,h):
    """Barycentric lattice dense enough that no h cell crossed by a face is missed."""
    triangles=np.asarray(triangles,float)
    if triangles.ndim!=3 or triangles.shape[1:]!=(3,3):raise ValueError('Triangle corner array required')
    h=_positive(h,'cell size')
    if not len(triangles):return np.zeros((0,3))
    area=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)/2
    per=np.clip(np.ceil(np.sqrt(area/((h/2)**2))).astype(np.int64),1,64)
    out=[]
    for count in np.unique(per):
        picked=triangles[per==count]
        u=(np.arange(count)+.5)/count
        a,b=np.meshgrid(u,u,indexing='ij')
        w=np.stack([a.ravel(),b.ravel()],axis=1)
        w=np.where(w.sum(axis=1,keepdims=True)>1,1-w,w)
        bary=np.concatenate([1-w.sum(axis=1,keepdims=True),w],axis=1)
        out.append((picked[:,None,:,:]*bary[None,:,:,None]).sum(axis=2).reshape(-1,3))
    return np.concatenate(out,axis=0)


def voxel_keys(vertices,triangles,h):
    """Occupied h cells of a surface: its vertices plus a face-dense sample lattice."""
    vertices=np.asarray(vertices,float);triangles=np.asarray(triangles,np.int64)
    if vertices.ndim!=2 or vertices.shape[1]!=3 or not np.isfinite(vertices).all():
        raise ValueError('Finite vertex array required')
    if triangles.ndim!=2 or triangles.shape[1]!=3 or (len(triangles) and (triangles.min()<0 or triangles.max()>=len(vertices))):
        raise ValueError('Valid triangle index array required')
    points=np.concatenate([vertices,surface_samples(vertices[triangles],h)],axis=0)
    return np.unique(encode_cells(np.floor(points/h).astype(np.int64)))


def aabb_bounds(vertices):
    v=np.asarray(vertices,float)
    if v.ndim!=2 or v.shape[1]!=3 or not len(v) or not np.isfinite(v).all():raise ValueError('Finite vertices required')
    return v.min(axis=0),v.max(axis=0)


def _checked_boxes(lo,hi,margin_m):
    lo=np.asarray(lo,float);hi=np.asarray(hi,float)
    if lo.ndim!=2 or lo.shape[1]!=3 or lo.shape!=hi.shape or not np.isfinite(lo).all() or not np.isfinite(hi).all():
        raise ValueError('Finite (n,3) bounds required')
    if np.any(hi<lo):raise ValueError('Inverted bounding box')
    margin=_nonnegative(margin_m,'broad phase margin')
    return lo-margin,hi+margin


def dense_aabb_pairs(lo,hi,margin_m=0.):
    """Exact axis-aligned overlap by the O(N^2) boolean matrix the assessment used."""
    low,high=_checked_boxes(lo,hi,margin_m)
    n=len(low);overlap=np.ones((n,n),bool)
    for k in range(3):
        overlap&=(low[:,k][:,None]<=high[:,k][None,:])&(low[:,k][None,:]<=high[:,k][:,None])
    np.fill_diagonal(overlap,False)
    i,j=np.nonzero(np.triu(overlap,1))
    return np.stack([i,j],axis=1).astype(np.int64),{'algorithm':'dense','overlap_tests':int(n*(n-1)//2)}


def sweep_and_prune_pairs(lo,hi,margin_m=0.,axis=None):
    """Single-axis sweep and prune. Same pair set as dense AABB, fewer tests."""
    low,high=_checked_boxes(lo,hi,margin_m)
    n=len(low)
    if n<2:return np.zeros((0,2),np.int64),{'algorithm':'sweep_and_prune','axis':0,'overlap_tests':0,'interval_tests':0}
    if axis is None:axis=int(np.argmax(high.max(axis=0)-low.min(axis=0)))
    if axis not in (0,1,2):raise ValueError('Sweep axis must be 0, 1 or 2')
    order=np.argsort(low[:,axis],kind='stable')
    slo=low[order];shi=high[order]
    other=[k for k in range(3) if k!=axis]
    out=[];overlap_tests=0;interval_tests=0
    for a in range(n-1):
        limit=shi[a,axis]
        # Sorted starts: the sweep front ends at the first body starting past a's end.
        stop=int(np.searchsorted(slo[a+1:,axis],limit,side='right'))
        if stop<=0:continue
        window=slice(a+1,a+1+stop)
        interval_tests+=stop
        keep=np.ones(stop,bool)
        for k in other:
            keep&=(slo[window,k]<=shi[a,k])&(slo[a,k]<=shi[window,k])
        overlap_tests+=stop
        hit=np.flatnonzero(keep)
        if len(hit):
            partner=order[a+1+hit]
            out.append(np.stack([np.full(len(hit),order[a]),partner],axis=1))
    pairs=np.concatenate(out,axis=0) if out else np.zeros((0,2),np.int64)
    pairs=np.sort(pairs,axis=1)
    if len(pairs):pairs=np.unique(pairs,axis=0)
    return pairs.astype(np.int64),{'algorithm':'sweep_and_prune','axis':axis,
        'overlap_tests':int(overlap_tests),'interval_tests':int(interval_tests)}


@dataclass
class _Bvh:
    lo:np.ndarray
    hi:np.ndarray
    left:np.ndarray
    right:np.ndarray
    start:np.ndarray
    count:np.ndarray
    order:np.ndarray


def _build_bvh(low,high,leaf_size):
    n=len(low);order=np.arange(n)
    centre=(low+high)/2
    node_lo=[];node_hi=[];node_left=[];node_right=[];node_start=[];node_count=[]
    def new_node():
        node_lo.append(np.zeros(3));node_hi.append(np.zeros(3))
        node_left.append(-1);node_right.append(-1);node_start.append(0);node_count.append(0)
        return len(node_lo)-1
    root=new_node()
    work=[(root,0,n)]
    while work:
        node,s,e=work.pop()
        members=order[s:e]
        node_lo[node]=low[members].min(axis=0);node_hi[node]=high[members].max(axis=0)
        if e-s<=leaf_size:
            node_start[node]=s;node_count[node]=e-s;continue
        c=centre[members]
        extent=c.max(axis=0)-c.min(axis=0)
        axis=int(np.argmax(extent))
        if extent[axis]<=0:
            node_start[node]=s;node_count[node]=e-s;continue
        local=np.argsort(c[:,axis],kind='stable')
        order[s:e]=members[local]
        mid=s+(e-s)//2
        l=new_node();r=new_node()
        node_left[node]=l;node_right[node]=r
        work.append((l,s,mid));work.append((r,mid,e))
    return _Bvh(np.array(node_lo),np.array(node_hi),np.array(node_left,int),np.array(node_right,int),
                np.array(node_start,int),np.array(node_count,int),order)


def bvh_pairs(lo,hi,margin_m=0.,leaf_size=8):
    """Median-split AABB BVH with a self-descent. Same pair set as dense AABB."""
    low,high=_checked_boxes(lo,hi,margin_m)
    if leaf_size<1 or int(leaf_size)!=leaf_size:raise ValueError('Positive integer leaf size required')
    n=len(low)
    if n<2:return np.zeros((0,2),np.int64),{'algorithm':'bvh','nodes':0,'node_tests':0,'overlap_tests':0}
    bvh=_build_bvh(low,high,int(leaf_size))
    node_tests=0;overlap_tests=0;out=[]
    def leaf_members(node):
        s=bvh.start[node];return bvh.order[s:s+bvh.count[node]]
    def box_overlap(a,b):
        return bool(np.all(bvh.lo[a]<=bvh.hi[b]) and np.all(bvh.lo[b]<=bvh.hi[a]))
    def brute(mi,mj,same):
        nonlocal overlap_tests
        if same:
            if len(mi)<2:return
            i,j=np.triu_indices(len(mi),1)
            a=mi[i];b=mi[j]
        else:
            a=np.repeat(mi,len(mj));b=np.tile(mj,len(mi))
        overlap_tests+=len(a)
        keep=np.ones(len(a),bool)
        for k in range(3):
            keep&=(low[a,k]<=high[b,k])&(low[b,k]<=high[a,k])
        hit=np.flatnonzero(keep)
        if len(hit):out.append(np.stack([a[hit],b[hit]],axis=1))
    stack=[(0,0,True)]
    while stack:
        x,y,same=stack.pop()
        if same:
            if bvh.count[x]:
                brute(leaf_members(x),None,True);continue
            l,r=bvh.left[x],bvh.right[x]
            stack.append((l,l,True));stack.append((r,r,True));stack.append((l,r,False))
        else:
            node_tests+=1
            if not box_overlap(x,y):continue
            lx,ly=bvh.count[x],bvh.count[y]
            if lx and ly:
                brute(leaf_members(x),leaf_members(y),False);continue
            if not lx and (ly or _volume(bvh,x)>=_volume(bvh,y)):
                stack.append((bvh.left[x],y,False));stack.append((bvh.right[x],y,False))
            else:
                stack.append((x,bvh.left[y],False));stack.append((x,bvh.right[y],False))
    pairs=np.concatenate(out,axis=0) if out else np.zeros((0,2),np.int64)
    pairs=np.sort(pairs,axis=1)
    if len(pairs):pairs=np.unique(pairs,axis=0)
    return pairs.astype(np.int64),{'algorithm':'bvh','nodes':int(len(bvh.lo)),
        'node_tests':int(node_tests),'overlap_tests':int(overlap_tests)}


def _volume(bvh,node):
    e=np.maximum(bvh.hi[node]-bvh.lo[node],0.)
    return float(e[0]*e[1]*e[2])


@dataclass(frozen=True)
class BodySurface:
    """One deformable owner's boundary surface, in the frame the broad phase runs in."""
    body_id:str
    vertices_m:np.ndarray
    triangles:np.ndarray

    def bounds(self):return aabb_bounds(self.vertices_m)


class SurfaceCellIndex:
    """Uniform hash of surface samples: which bodies enter the same h cell."""

    def __init__(self,surfaces,cell_m):
        self.cell_m=_positive(cell_m,'surface cell size')
        keys=[];owners=[]
        for i,s in enumerate(surfaces):
            k=voxel_keys(s.vertices_m,s.triangles,self.cell_m)
            keys.append(k);owners.append(np.full(len(k),i,np.int32))
        self.key=np.concatenate(keys) if keys else np.zeros(0,np.int64)
        self.owner=np.concatenate(owners) if owners else np.zeros(0,np.int32)
        order=np.argsort(self.key,kind='stable')
        self.key=self.key[order];self.owner=self.owner[order]
        cut=np.flatnonzero(np.diff(self.key))+1
        self.starts=np.concatenate(([0],cut)) if len(self.key) else np.zeros(0,np.int64)
        self.ends=np.concatenate((cut,[len(self.key)])) if len(self.key) else np.zeros(0,np.int64)
        self.entries=int(len(self.key))
        self.occupied_cells=int(len(self.starts))

    def shared_cell_counts(self,restrict_to=None):
        """Shared h-cell count per unordered body pair. `restrict_to` prunes the scan."""
        allowed=None
        if restrict_to is not None:
            allowed=set()
            for a,b in np.asarray(restrict_to,np.int64):
                allowed.add((int(a),int(b)) if a<b else (int(b),int(a)))
        counts={}
        occupancy=self.ends-self.starts
        for s,e in zip(self.starts[occupancy>1],self.ends[occupancy>1]):
            members=np.unique(self.owner[s:e])
            for x in range(len(members)):
                for y in range(x+1,len(members)):
                    key=(int(members[x]),int(members[y]))
                    if allowed is not None and key not in allowed:continue
                    counts[key]=counts.get(key,0)+1
        return counts


@dataclass(frozen=True)
class BroadPhaseResult:
    pairs:tuple
    shared_cells:dict
    body_ids:tuple
    stats:dict

    def as_index_pairs(self):
        return np.array(self.pairs,np.int64).reshape(-1,2)

    def as_id_pairs(self):
        return tuple((self.body_ids[a],self.body_ids[b]) for a,b in self.pairs)


class BroadPhase:
    """Two stage discovery: bounding prune, then measured surface co-occupancy.

    Stage one is only a superset filter, and it is a superset only if its margin
    covers the cell quantisation of stage two: two surfaces sharing one h cell
    differ by less than h on every axis, so a margin of h/2 per box is the
    smallest one that cannot produce a false negative. The default enforces it.
    """

    def __init__(self,surfaces,*,cell_m,margin_m=None,prune='sweep_and_prune',leaf_size=8):
        self.surfaces=tuple(surfaces)
        if not self.surfaces:raise ValueError('At least one surface required')
        ids=[s.body_id for s in self.surfaces]
        if len(set(ids))!=len(ids):raise ValueError('Duplicate body id in broad phase input')
        self.body_ids=tuple(ids)
        self.cell_m=_positive(cell_m,'surface cell size')
        self.margin_m=self.cell_m/2 if margin_m is None else _nonnegative(margin_m,'broad phase margin')
        if prune not in ('sweep_and_prune','bvh','dense','none'):raise ValueError('Unknown prune stage')
        self.prune=prune;self.leaf_size=leaf_size
        bounds=[s.bounds() for s in self.surfaces]
        self.lo=np.array([b[0] for b in bounds]);self.hi=np.array([b[1] for b in bounds])

    def candidate_pairs(self):
        if self.prune=='none':
            n=len(self.surfaces);i,j=np.triu_indices(n,1)
            return np.stack([i,j],axis=1).astype(np.int64),{'algorithm':'none','overlap_tests':0}
        if self.prune=='dense':return dense_aabb_pairs(self.lo,self.hi,self.margin_m)
        if self.prune=='bvh':return bvh_pairs(self.lo,self.hi,self.margin_m,self.leaf_size)
        return sweep_and_prune_pairs(self.lo,self.hi,self.margin_m)

    def discover(self):
        candidates,prune_stats=self.candidate_pairs()
        index=SurfaceCellIndex(self.surfaces,self.cell_m)
        counts=index.shared_cell_counts(restrict_to=candidates)
        pairs=tuple(sorted(counts))
        stats={'bodies':len(self.surfaces),'cell_m':self.cell_m,'margin_m':self.margin_m,
               'prune':prune_stats,'candidate_pairs':int(len(candidates)),
               'surface_cell_entries':index.entries,'occupied_surface_cells':index.occupied_cells,
               'discovered_pairs':len(pairs),
               'prune_precision':float(len(pairs)/len(candidates)) if len(candidates) else 0.}
        return BroadPhaseResult(pairs,dict(counts),self.body_ids,stats)


@dataclass(frozen=True)
class InterfacePartition:
    contact_pairs:tuple
    tied_pairs:tuple
    unknown_pairs:tuple
    counts:dict


def partition_interfaces(result,labels,*,unknown_as='contact'):
    """Split discovered pairs into contacts and anatomically tied constraints.

    A tied interface is a constraint. Handing it to the contact kernel welds it
    with a one-sided impulse law and destroys the very sliding the per-entity
    representation exists to express, so tied pairs are returned separately and
    never as contacts. `labels` maps an unordered id pair to a recorded interface
    label; anything the classification did not resolve is reported as unknown and
    routed by an explicit policy, never silently.
    """
    if unknown_as not in ('contact','tied','refuse'):raise ValueError('Unknown-interface policy required')
    tied_labels={'tied','tied_by_partition','tied_by_tree_continuity'}
    contact_labels={'slide'}
    lookup={}
    for key,value in labels.items():
        a,b=key
        lookup[(a,b) if a<b else (b,a)]=value
    contact=[];tied=[];unknown=[];counts={}
    for i,j in result.pairs:
        ai,bi=result.body_ids[i],result.body_ids[j]
        label=lookup.get((ai,bi) if ai<bi else (bi,ai),'unknown')
        counts[label]=counts.get(label,0)+1
        if label in tied_labels:tied.append((i,j))
        elif label in contact_labels:contact.append((i,j))
        else:
            unknown.append((i,j))
            if unknown_as=='contact':contact.append((i,j))
            elif unknown_as=='tied':tied.append((i,j))
    if unknown and unknown_as=='refuse':
        raise ValueError(f'{len(unknown)} discovered interfaces carry no recorded classification')
    return InterfacePartition(tuple(contact),tuple(tied),tuple(unknown),counts)


def self_collision_candidates(vertices_m,triangles,*,cell_m,ring=1):
    """Triangle pairs of ONE body that share an h cell, excluding topological neighbours.

    A triangle always occupies the cells of the triangles it shares an edge or a
    vertex with. Reporting those is not a collision, it is the mesh. `ring`
    excludes pairs within that many shared-vertex hops; ring 1 is shared-vertex.
    """
    vertices=np.asarray(vertices_m,float);tri=np.asarray(triangles,np.int64)
    if tri.ndim!=2 or tri.shape[1]!=3 or not len(tri):raise ValueError('Triangle array required')
    if ring<0 or int(ring)!=ring:raise ValueError('Nonnegative integer ring required')
    h=_positive(cell_m,'self collision cell size')
    corners=vertices[tri]
    lo=corners.min(axis=1);hi=corners.max(axis=1)
    cells=[];owners=[]
    for f in range(len(tri)):
        a=np.floor(lo[f]/h).astype(np.int64);b=np.floor(hi[f]/h).astype(np.int64)
        grid=np.stack(np.meshgrid(*[np.arange(a[k],b[k]+1) for k in range(3)],indexing='ij'),axis=-1).reshape(-1,3)
        cells.append(encode_cells(grid));owners.append(np.full(len(grid),f,np.int64))
    key=np.concatenate(cells);owner=np.concatenate(owners)
    order=np.argsort(key,kind='stable');key=key[order];owner=owner[order]
    cut=np.flatnonzero(np.diff(key))+1
    starts=np.concatenate(([0],cut));ends=np.concatenate((cut,[len(key)]))
    raw=set()
    for s,e in zip(starts,ends):
        if e-s<2:continue
        members=np.unique(owner[s:e])
        for x in range(len(members)):
            for y in range(x+1,len(members)):
                raw.add((int(members[x]),int(members[y])))
    neighbour=_ring_neighbours(tri,int(ring))
    kept=sorted(p for p in raw if p not in neighbour)
    return {'candidate_pairs':np.array(kept,np.int64).reshape(-1,2),
            'cell_pairs_before_topology':len(raw),
            'topological_pairs_excluded':len(raw)-len(kept),
            'ring':int(ring),'cell_m':h,'triangles':int(len(tri)),
            'cell_entries':int(len(key)),'occupied_cells':int(len(starts))}


def _ring_neighbours(tri,ring):
    """Unordered triangle pairs within `ring` shared-vertex hops."""
    if ring<=0:return set()
    by_vertex={}
    for f,face in enumerate(tri):
        for v in face:by_vertex.setdefault(int(v),[]).append(f)
    adjacency={f:set() for f in range(len(tri))}
    for members in by_vertex.values():
        for a in members:
            adjacency[a].update(members)
    for a in adjacency:adjacency[a].discard(a)
    reach={f:set(adjacency[f]) for f in range(len(tri))}
    for _ in range(ring-1):
        grown={}
        for f,s in reach.items():
            g=set(s)
            for n in s:g.update(adjacency[n])
            g.discard(f);grown[f]=g
        reach=grown
    out=set()
    for f,s in reach.items():
        for g in s:
            if f<g:out.add((f,g))
    return out


@dataclass(frozen=True)
class SymmetricContact:
    """One declaration for one interface, expanded into the kernel's two passes.

    `NodeTriangleContact` is one directional: nodes of A against triangles of B.
    An interface therefore needs two declarations today, and nothing checks that
    the second was written or that both carry the same friction and search band.
    This carries the interface once and expands it, so the two passes cannot
    disagree. It changes no impulse: the expansion is exactly the two directional
    passes the kernel already performs, in a fixed order, and every conservation
    receipt is reported per pass and summed.
    """
    a_id:str
    b_id:str
    friction_static:float
    friction_kinetic:float
    search_distance_m:float
    node_ids_a:object=None
    node_ids_b:object=None
    triangle_ids_a:object=None
    triangle_ids_b:object=None

    def __post_init__(self):
        if self.a_id==self.b_id:raise ValueError('A symmetric interface needs two distinct owners')


def expand_symmetric(contact,node_triangle_contact_cls):
    """Return the two directional passes of one symmetric interface, in order."""
    forward=node_triangle_contact_cls(a_id=contact.a_id,b_id=contact.b_id,
        friction_static=contact.friction_static,friction_kinetic=contact.friction_kinetic,
        search_distance_m=contact.search_distance_m,node_ids=contact.node_ids_a,
        triangle_ids=contact.triangle_ids_b)
    reverse=node_triangle_contact_cls(a_id=contact.b_id,b_id=contact.a_id,
        friction_static=contact.friction_static,friction_kinetic=contact.friction_kinetic,
        search_distance_m=contact.search_distance_m,node_ids=contact.node_ids_b,
        triangle_ids=contact.triangle_ids_a)
    return forward,reverse


def merge_symmetric_receipts(forward,reverse):
    """Sum the two directional receipts of one interface without hiding either."""
    out={'passes':(forward,reverse)}
    for key in ('dissipation_j','kinetic_transfer_a_j','kinetic_transfer_b_j'):
        out[key]=float(forward[key]+reverse[key])
    for key in ('contact_count','unresolved_edge_contacts'):
        out[key]=int(forward[key]+reverse[key])
    out['max_preprojection_penetration_m']=float(max(forward['max_preprojection_penetration_m'],
                                                     reverse['max_preprojection_penetration_m']))
    for key in ('paired_impulse_residual_ns','angular_impulse_residual_nms'):
        out[key]=np.asarray(forward[key])+np.asarray(reverse[key])
    out['interface_impulse_residual_ns']=(np.asarray(forward['impulse_a_ns']).sum(axis=0)
        +np.asarray(forward['impulse_b_ns']).sum(axis=0)+np.asarray(reverse['impulse_a_ns']).sum(axis=0)
        +np.asarray(reverse['impulse_b_ns']).sum(axis=0))
    return out


def surface_edges(triangles):
    """Unique undirected edges of a triangle surface, in ascending vertex order."""
    t=np.asarray(triangles,np.int64)
    if t.ndim!=2 or t.shape[1]!=3 or not len(t):raise ValueError('Triangle array required')
    e=np.concatenate((t[:,[0,1]],t[:,[1,2]],t[:,[2,0]]),axis=0)
    return np.unique(np.sort(e,axis=1),axis=0)


def _segment_closest(pa0,pa1,pb0,pb1):
    """Closest points of two segments, returned with their interior parameters."""
    d1=pa1-pa0;d2=pb1-pb0;r=pa0-pb0
    a=np.einsum('ij,ij->i',d1,d1);e=np.einsum('ij,ij->i',d2,d2)
    f=np.einsum('ij,ij->i',d2,r);c=np.einsum('ij,ij->i',d1,r);b=np.einsum('ij,ij->i',d1,d2)
    den=a*e-b*b
    with np.errstate(divide='ignore',invalid='ignore'):
        s=np.where(den>1e-30,(b*f-c*e)/np.where(den>1e-30,den,1.),0.)
    s=np.clip(s,0.,1.)
    with np.errstate(divide='ignore',invalid='ignore'):
        t=np.where(e>1e-30,(b*s+f)/np.where(e>1e-30,e,1.),0.)
    t_clipped=np.clip(t,0.,1.)
    with np.errstate(divide='ignore',invalid='ignore'):
        s=np.where(a>1e-30,np.clip((b*t_clipped-c)/np.where(a>1e-30,a,1.),0.,1.),s)
    p=pa0+s[:,None]*d1;q=pb0+t_clipped[:,None]*d2
    return s,t_clipped,p,q,den


def resolve_edge_edge_contact(position_a_m,velocity_a_m_s,mass_a_kg,
                              position_b_m,velocity_b_m_s,mass_b_kg,
                              *,edges_a,edges_b,friction_static,friction_kinetic,
                              search_distance_m,mobile_a=None,mobile_b=None):
    """Zero-restitution edge/edge impacts with the node/triangle kernel's audit.

    `resolve_node_triangle_contact` refuses an edge case rather than selecting a
    farther plane, and counts it as `unresolved_edge_contacts`. Two crossed
    ridges therefore pass straight through each other. This closes that case with
    the SAME impulse algebra: one contact point carrying two weights on each side
    instead of one and three, so the paired impulse cancels exactly and the
    dissipation identity 0.5*jn^2*w - jt.(u_t + 0.5*w*jt) is unchanged.

    Two deliberate differences from the triangle kernel, stated rather than hidden:

    * No position projection. An edge crossing carries no inside/outside test, so
      the penetration sign is not defined and a projection would be a guess. The
      constraint is applied at velocity level on approach only, which means the
      search distance must exceed the per-step approach for it to bite. The
      minimum gap is reported so that failure is visible.
    * Nearly parallel edges are refused, not resolved: their common normal is
      ill-conditioned. They are counted in `unresolved_parallel_edges`.
    """
    xa=np.array(position_a_m,float,copy=True);xb=np.array(position_b_m,float,copy=True)
    va=np.array(velocity_a_m_s,float,copy=True);vb=np.array(velocity_b_m_s,float,copy=True)
    ma=np.asarray(mass_a_kg,float);mb=np.asarray(mass_b_kg,float)
    ea=np.asarray(edges_a,np.int64);eb=np.asarray(edges_b,np.int64)
    kinetic_a_before=float(.5*np.sum(ma[:,None]*va*va));kinetic_b_before=float(.5*np.sum(mb[:,None]*vb*vb))
    for x,v,m in ((xa,va,ma),(xb,vb,mb)):
        if x.ndim!=2 or x.shape[1]!=3 or x.shape!=v.shape or m.shape!=(len(x),) or not all(np.isfinite(z).all() for z in (x,v,m)) or np.any(m<=0):
            raise ValueError('Finite nodal contact arrays and positive mass required')
    for e,x in ((ea,xa),(eb,xb)):
        if e.ndim!=2 or e.shape[1]!=2 or not len(e) or e.min()<0 or e.max()>=len(x) or np.any(e[:,0]==e[:,1]):
            raise ValueError('Valid distinct edge index pairs required')
    for value in (friction_static,friction_kinetic):
        if isinstance(value,bool) or not np.isfinite(value) or value<0:raise ValueError('Nonnegative finite Coulomb coefficient required')
    if friction_kinetic>friction_static:raise ValueError('Kinetic friction exceeds static friction')
    distance=_positive(search_distance_m,'contact search distance')
    move_a=np.ones(len(xa),bool) if mobile_a is None else np.asarray(mobile_a,bool)
    move_b=np.ones(len(xb),bool) if mobile_b is None else np.asarray(mobile_b,bool)
    if move_a.shape!=ma.shape or move_b.shape!=mb.shape:raise ValueError('Invalid mobility masks')
    from scipy.spatial import cKDTree
    wa=move_a/ma;wb=move_b/mb
    ja=np.zeros_like(xa);jb=np.zeros_like(xb)
    b0=xb[eb[:,0]];b1=xb[eb[:,1]]
    mid=(b0+b1)/2;half=float(np.linalg.norm(b1-b0,axis=1).max()/2)
    tree=cKDTree(mid)
    loss=0.;contacts=0;parallel=0;endpoint=0;angular=np.zeros(3)
    minimum_gap=np.inf
    for k in range(len(ea)):
        a0=xa[ea[k,0]];a1=xa[ea[k,1]]
        centre=(a0+a1)/2;reach=float(np.linalg.norm(a1-a0)/2)
        candidates=np.asarray(tree.query_ball_point(centre,half+reach+distance),int)
        if not len(candidates):continue
        s,t,p,q,den=_segment_closest(np.repeat(a0[None],len(candidates),0),np.repeat(a1[None],len(candidates),0),
                                     b0[candidates],b1[candidates])
        d=np.linalg.norm(p-q,axis=1)
        j=int(np.argmin(d))
        gap=float(d[j])
        if gap>distance:continue
        minimum_gap=min(minimum_gap,gap)
        # An endpoint solution is a vertex/edge or vertex/face case and belongs to
        # the node/triangle kernel, which owns the inside/outside test.
        if not (1e-9<s[j]<1-1e-9 and 1e-9<t[j]<1-1e-9):endpoint+=1;continue
        d1=a1-a0;d2=b1[candidates[j]]-b0[candidates[j]]
        axis=np.cross(d1,d2);axis_norm=float(np.linalg.norm(axis))
        if axis_norm<=1e-9*float(np.linalg.norm(d1)*np.linalg.norm(d2)):parallel+=1;continue
        normal=axis/axis_norm
        separation=p[j]-q[j]
        if float(separation@normal)<0:normal=-normal
        na=ea[k];nb=eb[candidates[j]]
        alpha=np.array([1-s[j],s[j]]);beta=np.array([1-t[j],t[j]])
        inverse_mass=float(np.sum(alpha**2*wa[na])+np.sum(beta**2*wb[nb]))
        if inverse_mass==0:continue
        relative=alpha@va[na]-beta@vb[nb];vn=float(relative@normal)
        if vn>=0:continue
        contacts+=1
        jn=-vn/inverse_mass
        tangential=relative-vn*normal;speed=float(np.linalg.norm(tangential))
        needed=speed/inverse_mass
        magnitude=needed if needed<=friction_static*jn else min(needed,friction_kinetic*jn)
        jt=-magnitude*tangential/speed if speed>0 else np.zeros(3)
        impulse=jn*normal+jt
        loss+=.5*jn*jn*inverse_mass-float(jt@(tangential+.5*inverse_mass*jt))
        contribution_a=alpha[:,None]*impulse[None,:]
        contribution_b=-beta[:,None]*impulse[None,:]
        va[na]+=wa[na,None]*contribution_a;vb[nb]+=wb[nb,None]*contribution_b
        np.add.at(ja,na,contribution_a);np.add.at(jb,nb,contribution_b)
        angular+=np.cross(xa[na],contribution_a).sum(axis=0)+np.cross(xb[nb],contribution_b).sum(axis=0)
    return {'position_a_m':xa,'position_b_m':xb,'velocity_a_m_s':va,'velocity_b_m_s':vb,
            'impulse_a_ns':ja,'impulse_b_ns':jb,
            'position_correction_a_m':np.zeros_like(xa),'position_correction_b_m':np.zeros_like(xb),
            'dissipation_j':float(loss),'contact_count':contacts,
            'unresolved_parallel_edges':parallel,'endpoint_solutions_deferred':endpoint,
            'minimum_edge_gap_m':float(minimum_gap) if np.isfinite(minimum_gap) else None,
            'paired_impulse_residual_ns':ja.sum(axis=0)+jb.sum(axis=0),
            'angular_impulse_residual_nms':angular,
            'kinetic_transfer_a_j':float(.5*np.sum(ma[:,None]*va*va)-kinetic_a_before),
            'kinetic_transfer_b_j':float(.5*np.sum(mb[:,None]*vb*vb)-kinetic_b_before)}


def contact_triangle_ids(surface_a,surface_b,*,cell_m,dilation=1):
    """Triangles of B whose cells lie within `dilation` cells of A's occupied cells.

    The kernel accepts a triangle subset per declared pair. Passing the whole
    target mesh makes its single global enclosing-sphere radius the worst triangle
    in the entire body, and makes its KD tree the whole body. This narrows both
    without changing which triangle is nearest, provided the dilation covers the
    search band.
    """
    h=_positive(cell_m,'contact cell size')
    keys_a=voxel_keys(surface_a.vertices_m,surface_a.triangles,h)
    idx=np.stack([((keys_a>>40)&0xFFFFF),((keys_a>>20)&0xFFFFF),(keys_a&0xFFFFF)],axis=1)
    offsets=np.array([[i,j,k] for i in range(-dilation,dilation+1)
                      for j in range(-dilation,dilation+1) for k in range(-dilation,dilation+1)],np.int64)
    grown=np.unique((idx[:,None,:]+offsets[None]).reshape(-1,3),axis=0)
    grown_keys=np.unique((grown[:,0]<<40)|(grown[:,1]<<20)|grown[:,2])
    corners=surface_b.vertices_m[surface_b.triangles]
    lo=np.floor(corners.min(axis=1)/h).astype(np.int64);hi=np.floor(corners.max(axis=1)/h).astype(np.int64)
    keep=np.zeros(len(surface_b.triangles),bool)
    for f in range(len(surface_b.triangles)):
        a=lo[f];b=hi[f]
        grid=np.stack(np.meshgrid(*[np.arange(a[k],b[k]+1) for k in range(3)],indexing='ij'),axis=-1).reshape(-1,3)
        if np.isin(encode_cells(grid),grown_keys).any():keep[f]=True
    return np.flatnonzero(keep)
