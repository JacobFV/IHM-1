"""Objective elastic tethers and exact source-triangle garment/body pairings.

Tethers are explicit engineering supports, not anatomical attachments or friction.
No forces are applied to a scene implicitly.
"""
from dataclasses import dataclass
from pathlib import Path
import gzip,json
import numpy as np
from scipy.spatial import cKDTree
from .contact_dynamics import _triangle_projection
from .whole_garment import materialize_garments

@dataclass(frozen=True)
class BarycentricAttachment:
    garment_node:int
    body_nodes:tuple
    barycentric:tuple
    rest_length_m:float
    stiffness_n_m:float
    def __post_init__(self):
        ids=np.asarray(self.body_nodes);w=np.asarray(self.barycentric,float)
        if isinstance(self.garment_node,(bool,np.bool_)) or not isinstance(self.garment_node,(int,np.integer)) or self.garment_node<0:raise ValueError('Nonnegative garment node required')
        if ids.shape!=(3,) or ids.dtype.kind not in 'iu' or np.any(ids<0) or len(set(ids))!=3:raise ValueError('Three distinct body nodes required')
        if w.shape!=(3,) or not np.isfinite(w).all() or np.any(w<0) or abs(w.sum()-1)>1e-12:raise ValueError('Convex barycentric weights required')
        for v,zero in ((self.rest_length_m,True),(self.stiffness_n_m,False)):
            if isinstance(v,(bool,np.bool_)) or not np.isscalar(v) or not np.isfinite(v) or v<0 or (not zero and v==0):raise ValueError('Finite physical spring parameters required')


def attachment_forces(garment_positions_m,body_positions_m,attachments,*,wrench_origin_m=(0.,0.,0.)):
    a=np.asarray(garment_positions_m,float);b=np.asarray(body_positions_m,float);origin=np.asarray(wrench_origin_m,float)
    if any(x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all() for x in (a,b)) or origin.shape!=(3,) or not np.isfinite(origin).all():raise ValueError('Finite positions and wrench origin required')
    fa=np.zeros_like(a);fb=np.zeros_like(b);energy=0.;records=[]
    for p in attachments:
        ids=np.asarray(p.body_nodes);w=np.asarray(p.barycentric)
        if p.garment_node>=len(a) or ids.max()>=len(b):raise ValueError('Attachment outside owner')
        q=w@b[ids];delta=a[p.garment_node]-q;length=float(np.linalg.norm(delta))
        if length<=1e-14 and p.rest_length_m>0:raise ValueError('Collapsed nonzero-rest attachment has unresolved force direction')
        extension=length-p.rest_length_m
        force=-p.stiffness_n_m*delta if p.rest_length_m==0 else -p.stiffness_n_m*extension*delta/length
        fa[p.garment_node]+=force;np.add.at(fb,ids,-w[:,None]*force);energy+=.5*p.stiffness_n_m*extension**2
        records.append({'length_m':length,'extension_m':extension,'force_n':force.tolist(),'body_point_m':q.tolist()})
    torque_a=np.cross(a-origin,fa).sum(0);torque_b=np.cross(b-origin,fb).sum(0)
    return {'garment_force_n':fa,'body_force_n':fb,'energy_j':float(energy),'body_wrench':{'force_n':fb.sum(0),'torque_nm':torque_b,'origin_m':origin},
            'paired_force_residual_n':fa.sum(0)+fb.sum(0),'paired_torque_residual_nm':torque_a+torque_b,'attachments':records}


class SourceSurface:
    """Unmodified source triangles with conservative exact-nearest candidate search."""
    def __init__(self,x,tri):
        self.positions=np.asarray(x,float);self.triangles=np.asarray(tri,int);self.vertices=self.positions[self.triangles]
        self.centers=self.vertices.mean(1);self.radius=float(np.linalg.norm(self.vertices-self.centers[:,None],axis=2).max());self.tree=cKDTree(self.centers)
    def closest(self,point):
        point=np.asarray(point,float);_,seed=self.tree.query(point)
        # Any point on this seed bounds nearest distance; enclosing triangle
        # radii then include every triangle that could be nearer.
        upper=np.linalg.norm(point-self.vertices[seed,0]);candidates=np.asarray(self.tree.query_ball_point(point,upper+self.radius+1e-12),int)
        tri=self.vertices[candidates];gap,normals,beta,inside,d2=_triangle_projection(point,tri)
        k=int(np.argmin(d2));t=tri[k]
        if inside[k]:weights=np.maximum(beta[k],0);weights/=weights.sum();kind='face'
        else:
            alternatives=[]
            for i,j in ((0,1),(1,2),(2,0)):
                edge=t[j]-t[i];u=float(np.clip((point-t[i])@edge/(edge@edge),0,1));w=np.zeros(3);w[i]=1-u;w[j]=u
                alternatives.append((np.linalg.norm(point-w@t),w))
            _,weights=min(alternatives,key=lambda v:v[0]);kind='vertex' if np.max(weights)==1 else 'edge'
        q=weights@t;index=int(candidates[k])
        return {'triangle_index':index,'body_nodes':self.triangles[index].tolist(),'barycentric':weights.tolist(),'body_point_m':q.tolist(),
                'distance_m':float(np.linalg.norm(point-q)),'closest_feature':kind,'oriented_face_gap_m':float(gap[k]),
                'reconstruction_error_m':float(np.linalg.norm(q-weights@self.positions[self.triangles[index]]))}


def build_source_pairings(root):
    root=Path(root);bodies,identity=materialize_garments(root,areal_density_kg_m2=.18,edge_stiffness_n_m=12.)
    skin=json.loads(gzip.decompress((root/identity['skin']['path']).read_bytes()));x=np.asarray(skin['positions']).reshape(-1,3);tri=np.asarray(skin['indices']).reshape(-1,3);surface=SourceSurface(x,tri)
    records={}
    for name,body in bodies.items():
        audit=body.pattern_audit;supports=[];samples=[]
        for i,loop in enumerate(audit['boundary_loops']):
            label=audit['opening_labels'][str(i)];loop=np.array(loop)
            if name=='shorts' and label=='waistband':selected=loop[::4]
            elif name=='shirt' and 'arm_opening' in label:
                order=np.argsort(body.position_m[loop,1],kind='stable');selected=loop[order[-8:]]
            else:selected=[]
            for node in selected:supports.append({'garment_node':int(node),'region':label,'role':'optional_engineered_elastic_support',**surface.closest(body.position_m[node])})
        # Sparse body-interface samples are metadata only, never implicit tethers.
        for label,mask in ([('torso',body.position_m[:,1]>.15)] if name=='shirt' else [('hips',body.position_m[:,1]>-.1),('thighs',body.position_m[:,1]<-.2)]):
            ids=np.flatnonzero(mask)
            for node in ids[np.linspace(0,len(ids)-1,12,dtype=int)]:samples.append({'garment_node':int(node),'region':label,'role':'contact_query_metadata_only',**surface.closest(body.position_m[node])})
        records[name]={'support_pairings':supports,'contact_samples':samples}
    return {'schema':'ihm.garment-body-pairings.v1','identity':identity,'garments':records,
            'frame':{'units':'m','x':'left','y':'superior','z':'anterior'},
            'limitations':['Source triangle correspondence does not establish biological attachment','Tether rest lengths/stiffness must be explicit in any experiment','Nearest source face winding is retained; global exterior/cavity and containment semantics are not established','Contact metadata applies no constraint or force','Source coordinates and topology are unchanged']}
