"""Dynamic tetrahedral solids and finite-mass node/triangle Coulomb contact.

One clock advances every owner, including compatible Cloth objects. Position
projection and symplectic integration defects are measured, never hidden as
physical dissipation. This is nodal contact, not continuous/self collision.
"""
from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree
from .mechanics_backend import DeformableRegion

def _finite_positive(value,name):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value) or not np.isfinite(value) or value<=0:raise ValueError(f'Positive finite {name} required')
    return float(value)

def boundary_triangles(tets):
    t=np.asarray(tets,int)
    faces=np.concatenate((t[:,[1,2,3]],t[:,[0,3,2]],t[:,[0,1,3]],t[:,[0,2,1]]))
    _,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
    if np.any(count>2):raise ValueError('Nonmanifold volume interface')
    return faces[count[inv]==1]

class DynamicTetrahedra:
    """Finite-strain inertial solid sharing the independently checked NH law."""
    def __init__(self,vertices_m,tetrahedra,*,mu_pa,lambda_pa,density_kg_m3,fixed_nodes=(),material_owner_ids=()):
        self.region=DeformableRegion(vertices_m,tetrahedra,mu_pa=mu_pa,lambda_pa=lambda_pa,density_kg_m3=density_kg_m3)
        self.reference_position_m=self.region.reference.copy();self.position_m=self.reference_position_m.copy()
        self.velocity_m_s=np.zeros_like(self.position_m);self.mass_kg=self.region.nodal_mass.copy()
        if np.any(self.mass_kg<=0):raise ValueError('Every solid node must own positive material mass')
        self.triangles=boundary_triangles(self.region.tets);self.surface_nodes=np.unique(self.triangles)
        raw_fixed=np.asarray(fixed_nodes)
        if raw_fixed.ndim!=1 or (len(raw_fixed) and raw_fixed.dtype.kind not in 'iu'):raise ValueError('Integer fixed nodes required')
        self.fixed_nodes=raw_fixed.astype(int)
        if len(np.unique(self.fixed_nodes))!=len(self.fixed_nodes) or np.any(self.fixed_nodes<0) or np.any(self.fixed_nodes>=len(self.position_m)):raise ValueError('Invalid fixed nodes')
        self.material_owner_ids=tuple(material_owner_ids)
        if len(set(self.material_owner_ids))!=len(self.material_owner_ids):raise ValueError('Duplicate material owner within a body')
        self.time_s=0.
        x=self.reference_position_m[self.region.tets]
        areas=np.stack([np.linalg.norm(np.cross(x[:,b]-x[:,a],x[:,c]-x[:,a]),axis=1)/2 for a,b,c in ((1,2,3),(0,3,2),(0,1,3),(0,2,1))],axis=1)
        height=3*self.region.volumes/areas.max(axis=1)
        rho=np.broadcast_to(np.asarray(density_kg_m3,float),self.region.volumes.shape)
        wave=np.sqrt((self.region.lam+2*self.region.mu)/rho)
        self.max_explicit_dt_s=float(.1*np.min(height/wave))

    def elastic_forces(self,positions_m):
        energy,gradient=self.region.energy_gradient(positions_m)
        return -gradient,energy

    def minimum_jacobian(self,positions_m=None):
        return float(np.linalg.det(self.region.deformation(self.position_m if positions_m is None else positions_m)).min())

    def checkpoint(self):return (self.position_m.copy(),self.velocity_m_s.copy(),self.time_s)

    def restore(self,checkpoint):
        x,v,time=checkpoint
        x=np.asarray(x,float);v=np.asarray(v,float)
        self.elastic_forces(x)
        if v.shape!=x.shape or not np.isfinite(v).all() or not np.isfinite(time) or time<0:raise ValueError('Invalid dynamic checkpoint')
        self.position_m=x.copy();self.velocity_m_s=v.copy();self.time_s=float(time)

def _triangle_projection(point,triangles):
    e0=triangles[:,1]-triangles[:,0];e1=triangles[:,2]-triangles[:,0]
    normal=np.cross(e0,e1);norm=np.linalg.norm(normal,axis=1)
    if np.any(norm<1e-16):raise ValueError('Collapsed contact triangle')
    normal/=norm[:,None];r=point-triangles[:,0]
    gap=np.einsum('ij,ij->i',r,normal)
    aa=np.einsum('ij,ij->i',e0,e0);bb=np.einsum('ij,ij->i',e0,e1);cc=np.einsum('ij,ij->i',e1,e1)
    dd=np.einsum('ij,ij->i',r,e0);ee=np.einsum('ij,ij->i',r,e1);den=aa*cc-bb*bb
    beta1=(cc*dd-bb*ee)/den;beta2=(aa*ee-bb*dd)/den
    beta=np.stack((1-beta1-beta2,beta1,beta2),axis=1)
    inside=np.all(beta>=-1e-10,axis=1)
    # Exact nearest point distance also considers edges. A nearest edge/vertex
    # contact is explicitly reported unresolved instead of selecting a farther
    # plane behind a closed surface.
    distance2=np.where(inside,gap*gap,np.inf)
    for a,b in ((0,1),(1,2),(2,0)):
        edge=triangles[:,b]-triangles[:,a]
        u=np.clip(np.einsum('ij,ij->i',point-triangles[:,a],edge)/np.einsum('ij,ij->i',edge,edge),0,1)
        q=triangles[:,a]+u[:,None]*edge
        distance2=np.minimum(distance2,np.sum((point-q)**2,axis=1))
    return gap,normal,beta,inside,distance2

def resolve_node_triangle_contact(position_a_m,velocity_a_m_s,mass_a_kg,
                                  position_b_m,velocity_b_m_s,mass_b_kg,triangles_b,
                                  *,node_ids,friction_static,friction_kinetic,
                                  search_distance_m,mobile_a=None,mobile_b=None):
    """Sequential zero-restitution impacts with finite target triangle inertia.

    A node must be nearest an interior face projection, within the explicit
    search distance and on its negative side. Contact forces are not inferred
    from overlap appearance. Edges, self-collision and swept collision remain
    outside this kernel's coverage.
    """
    xa=np.array(position_a_m,float,copy=True);xb=np.array(position_b_m,float,copy=True)
    va=np.array(velocity_a_m_s,float,copy=True);vb=np.array(velocity_b_m_s,float,copy=True)
    ma=np.asarray(mass_a_kg,float);mb=np.asarray(mass_b_kg,float);tri=np.asarray(triangles_b)
    kinetic_a_before=float(.5*np.sum(ma[:,None]*va*va));kinetic_b_before=float(.5*np.sum(mb[:,None]*vb*vb))
    ids=np.asarray(node_ids)
    for x,v,m in ((xa,va,ma),(xb,vb,mb)):
        if x.ndim!=2 or x.shape[1]!=3 or x.shape!=v.shape or m.shape!=(len(x),) or not all(np.isfinite(z).all() for z in (x,v,m)) or np.any(m<=0):raise ValueError('Finite nodal contact arrays and positive mass required')
    if tri.ndim!=2 or tri.shape[1]!=3 or tri.dtype.kind not in 'iu' or not len(tri) or tri.min()<0 or tri.max()>=len(xb):raise ValueError('Invalid target triangles')
    if ids.ndim!=1 or ids.dtype.kind not in 'iu' or (len(ids) and (ids.min()<0 or ids.max()>=len(xa))) or len(np.unique(ids))!=len(ids):raise ValueError('Unique valid contact nodes required')
    for value in (friction_static,friction_kinetic):
        if isinstance(value,bool) or not np.isfinite(value) or value<0:raise ValueError('Nonnegative finite Coulomb coefficient required')
    if friction_kinetic>friction_static:raise ValueError('Kinetic friction exceeds static friction')
    distance=_finite_positive(search_distance_m,'contact search distance')
    move_a=np.ones(len(xa),bool) if mobile_a is None else np.asarray(mobile_a,bool)
    move_b=np.ones(len(xb),bool) if mobile_b is None else np.asarray(mobile_b,bool)
    if move_a.shape!=ma.shape or move_b.shape!=mb.shape:raise ValueError('Invalid mobility masks')
    wa=move_a/ma;wb=move_b/mb;ja=np.zeros_like(xa);jb=np.zeros_like(xb)
    correction_a=np.zeros_like(xa);correction_b=np.zeros_like(xb)
    loss=0.;contacts=0;edge_unresolved=0;penetration=0.;angular=np.zeros(3)
    reference_triangles=xb[tri];centers=reference_triangles.mean(axis=1)
    radius=float(np.linalg.norm(reference_triangles-centers[:,None,:],axis=2).max())
    tree=cKDTree(centers);target_motion_bound=0.
    for i in ids:
        # A conservative enclosing-sphere search, expanded by every accumulated
        # target projection, cannot discard a triangle within the search band.
        candidates=np.asarray(tree.query_ball_point(xa[i],radius+distance+target_motion_bound),int)
        if not len(candidates):continue
        gaps,normals,weights,inside,d2=_triangle_projection(xa[i],xb[tri[candidates]])
        k=int(np.argmin(d2))
        if d2[k]>distance**2 or gaps[k]>1e-12:continue
        if not inside[k]:edge_unresolved+=1;continue
        nodes=tri[candidates[k]];beta=weights[k];normal=normals[k];gap=float(gaps[k])
        inverse_mass=wa[i]+float(np.sum(beta**2*wb[nodes]))
        if inverse_mass==0:continue
        penetration=max(penetration,-gap);contacts+=1
        correction=max(0.,-gap)/inverse_mass*normal
        delta_a=wa[i]*correction;delta_b=-wb[nodes,None]*beta[:,None]*correction
        xa[i]+=delta_a;xb[nodes]+=delta_b;correction_a[i]+=delta_a;np.add.at(correction_b,nodes,delta_b)
        target_motion_bound+=float(np.linalg.norm(delta_b,axis=1).max())
        relative=va[i]-beta@vb[nodes];vn=float(relative@normal)
        jn=max(-vn,0.)/inverse_mass
        tangential=relative-vn*normal;speed=np.linalg.norm(tangential)
        needed=speed/inverse_mass
        magnitude=needed if needed<=friction_static*jn else min(needed,friction_kinetic*jn)
        jt=-magnitude*tangential/speed if speed>0 else np.zeros(3)
        impulse=jn*normal+jt;target=-beta[:,None]*impulse
        loss+=.5*jn*jn*inverse_mass-float(jt@(tangential+.5*inverse_mass*jt))
        va[i]+=wa[i]*impulse;vb[nodes]+=wb[nodes,None]*target
        ja[i]+=impulse;np.add.at(jb,nodes,target)
        angular+=np.cross(xa[i],impulse)+np.cross(xb[nodes],target).sum(axis=0)
    return {'position_a_m':xa,'position_b_m':xb,'velocity_a_m_s':va,'velocity_b_m_s':vb,
            'impulse_a_ns':ja,'impulse_b_ns':jb,'position_correction_a_m':correction_a,
            'position_correction_b_m':correction_b,'dissipation_j':float(loss),
            'contact_count':contacts,'unresolved_edge_contacts':edge_unresolved,
            'max_preprojection_penetration_m':float(penetration),
            'paired_impulse_residual_ns':ja.sum(axis=0)+jb.sum(axis=0),
            'angular_impulse_residual_nms':angular,
            'kinetic_transfer_a_j':float(.5*np.sum(ma[:,None]*va*va)-kinetic_a_before),
            'kinetic_transfer_b_j':float(.5*np.sum(mb[:,None]*vb*vb)-kinetic_b_before)}

@dataclass(frozen=True)
class NodeTriangleContact:
    a_id:str
    b_id:str
    friction_static:float
    friction_kinetic:float
    search_distance_m:float
    node_ids:object=None
    triangle_ids:object=None

def step_coupled(bodies,dt_s,*,contacts=(),external_forces_n=None,gravity_m_s2=(0.,-9.81,0.)):
    """Atomic common-clock symplectic step of all elastic mass owners."""
    dt=_finite_positive(dt_s,'step');gravity=np.asarray(gravity_m_s2,float)
    if gravity.shape!=(3,) or not np.isfinite(gravity).all() or not bodies:raise ValueError('Bodies and finite gravity required')
    if len({id(b) for b in bodies.values()})!=len(bodies):raise ValueError('Same material body supplied under multiple owners')
    owners=[owner for b in bodies.values() for owner in getattr(b,'material_owner_ids',())]
    if len(set(owners))!=len(owners):raise ValueError('Duplicate physical material owners')
    times=[b.time_s for b in bodies.values()]
    if not np.isfinite(times).all() or min(times)<0:raise ValueError('Invalid body clock')
    if max(times)-min(times)>1e-12:raise ValueError('Body clocks disagree')
    if any(dt>b.max_explicit_dt_s for b in bodies.values()):raise ValueError('Step exceeds an explicit elastic stability limit')
    ext={} if external_forces_n is None else external_forces_n
    if set(ext)-set(bodies):raise ValueError('External force references unknown owner')
    states={};energy_before=0.;support=np.zeros(3);external_impulse=np.zeros(3)
    for key,b in bodies.items():
        x=b.position_m.copy();v=b.velocity_m_s.copy();mass=b.mass_kg
        force,elastic=b.elastic_forces(x);external=np.asarray(ext.get(key,np.zeros_like(x)),float)
        if external.shape!=x.shape or not np.isfinite(external).all():raise ValueError('Invalid external nodal force')
        fixed=np.asarray(getattr(b,'fixed_nodes',()),int);mobile=np.ones(len(x),bool);mobile[fixed]=False
        total=force+external+mass[:,None]*gravity
        predicted_v=v+dt*total/mass[:,None];support-=np.sum(mass[fixed,None]*predicted_v[fixed],axis=0)
        predicted_v[fixed]=0.;predicted_x=x+dt*predicted_v
        states[key]={'oldx':x,'oldv':v,'x':predicted_x,'v':predicted_v,'mass':mass,'mobile':mobile,'external':external}
        energy_before+=elastic+.5*np.sum(mass[:,None]*v*v)-np.sum(mass[:,None]*x*gravity)
        external_impulse+=dt*np.sum(external+mass[:,None]*gravity,axis=0)
    receipts=[]
    for pair in contacts:
        if pair.a_id==pair.b_id or pair.a_id not in states or pair.b_id not in states:raise ValueError('Contact requires two distinct known owners')
        a,b=states[pair.a_id],states[pair.b_id];ba,bb=bodies[pair.a_id],bodies[pair.b_id]
        nodes=np.unique(ba.triangles) if pair.node_ids is None else np.asarray(pair.node_ids)
        triangles=bb.triangles if pair.triangle_ids is None else bb.triangles[np.asarray(pair.triangle_ids)]
        result=resolve_node_triangle_contact(a['x'],a['v'],a['mass'],b['x'],b['v'],b['mass'],triangles,
                    node_ids=nodes,friction_static=pair.friction_static,friction_kinetic=pair.friction_kinetic,
                    search_distance_m=pair.search_distance_m,mobile_a=a['mobile'],mobile_b=b['mobile'])
        a['x']=result['position_a_m'];a['v']=result['velocity_a_m_s'];b['x']=result['position_b_m'];b['v']=result['velocity_b_m_s']
        support-=result['impulse_a_ns'][~a['mobile']].sum(axis=0)+result['impulse_b_ns'][~b['mobile']].sum(axis=0)
        receipts.append(result)
    energy_after=0.;external_work=0.;momentum=np.zeros(3);elastic_total=0.;kinetic_total=0.
    for key,b in bodies.items():
        s=states[key];x,v,mass=s['x'],s['v'],s['mass']
        _,elastic=b.elastic_forces(x) # Positive J check before any owner's commit.
        kinetic=.5*np.sum(mass[:,None]*v*v)
        energy_after+=elastic+kinetic-np.sum(mass[:,None]*x*gravity)
        external_work+=np.sum(s['external']*(x-s['oldx']))
        momentum+=np.sum(mass[:,None]*(v-s['oldv']),axis=0)
        elastic_total+=elastic;kinetic_total+=kinetic
    loss=sum(r['dissipation_j'] for r in receipts)
    for key,b in bodies.items():b.position_m=states[key]['x'];b.velocity_m_s=states[key]['v'];b.time_s=times[0]+dt
    return {'time_s':times[0]+dt,'elastic_energy_j':float(elastic_total),'kinetic_energy_j':float(kinetic_total),
            'total_energy_j':float(energy_after),'external_work_j':float(external_work),'contact_dissipation_j':loss,
            'numerical_energy_defect_j':float(energy_after-energy_before-external_work+loss),
            'momentum_residual_ns':momentum-external_impulse-support,'support_impulse_ns':support,
            'contacts':receipts,'state_owners':list(bodies)}
