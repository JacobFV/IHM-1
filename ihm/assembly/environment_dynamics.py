"""Server-owned scene mechanics, coupled through canonical material force ports.

Engineering models, not calibrated tissue/contact measurements. Skin samples are
rigidly attached to native bone owners; they are not an independent human mass.
Cloth uses stretch/shear/bend springs, cushions a three-dimensional spring lattice,
and free props use six-degree-of-freedom Newton/Euler integration. Contacts use
compliant normal forces and regularized Coulomb friction. Exchange is explicit
on the embodied runtime's 5 ms mechanical exchange, with 1 ms environment substeps
and equal/opposite contact impulses.
"""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

CATALOGUE = 'data/derived/environment-catalogue-v1/catalogue.json'
# Bending rigidity in N*m. Chosen so the fabric bending length
# (B/(areal weight))**(1/3) is at least the 47 mm cloth vertex spacing, which is
# the coarsest drape this mesh can represent; a smaller value leaves creases the
# mesh cannot resolve. It sits in the heavy-fabric range but is an engineering
# value, not a measured Kawabata modulus.
CLOTH_BENDING_RIGIDITY_N_M = 2.5e-3
# Declared blanket sheet thickness for triangle-level self contact, in metres.
CLOTH_THICKNESS_M = .006
CLOTH_STRETCH_SWEEPS = 12
# Alternating material/contact projection passes per substep. Six left the
# strain ceiling unconverged on 59 of 200 held-body substeps at 1.12039; twelve
# converges every substep at 1.12010, and the loop exits early once it does, so
# a larger cap costs nothing when the constraints are compatible.
CLOTH_COVER_PASSES = 12
SCOPE = ('Server-owned spring cloth, volumetric spring-lattice pillow, rigid props and fixed compound colliders. '
         'Sampled skin contact loads act on the native articulated body at material points; reciprocal impulses act on props. '
         'Engineering coefficients, approximate skin registration, fixed cylinder/canopy bounding boxes and vertex rigid manifolds, 5 ms embodied exchange; approximate vertex cloth self-contact (no triangle collision), '
         'no resolved skin pressure/thermal insulation or clinically validated contact prediction. Mattress support remains native-owned.')


def resolve_selection(root, environment, selection=None):
    """Only catalogue identities are accepted; clients cannot supply physics or paths."""
    if selection is None:return None, {}
    if not isinstance(selection, dict) or set(selection)-{'scene','objects','bed_support_model','mattress_material','ambient_thermal'}:
        raise ValueError('Unknown environment selection')
    root=Path(root);catalog=json.loads((root/CATALOGUE).read_bytes())
    base={'supine':'bed','upright':'floor','free':'studio'}[environment]
    records={r['id']:r for key in ('scenes','objects','components') for r in catalog[key]}
    selected={};options={}
    for slot in ('scene','bed_support_model','mattress_material','ambient_thermal'):
        ident=selection.get(slot)
        if ident is None:continue
        if not isinstance(ident,str) or ident not in records or records[ident]['slot']!=slot:raise ValueError('Unknown '+slot)
        selected[slot]=ident
    ids=selection.get('objects',[])
    if not isinstance(ids,list) or len(ids)>16 or any(not isinstance(i,str) or i not in records or records[i]['slot']!='objects' for i in ids):
        raise ValueError('Expected at most 16 catalogue objects')
    selected['objects']=list(ids)
    for ident in [v for k,v in selected.items() if k!='objects']+ids:
        for requirement in records[ident].get('requires',[]):
            value=base if requirement['slot']=='environment' else selected.get(requirement['slot'])
            if value not in requirement['any_of']:raise ValueError('Unsatisfied environment dependency for '+ident)
    for slot in ('bed_support_model','mattress_material'):
        if slot in selected:
            for setting in records[selected[slot]]['selection']:
                if setting['parameter'] in ('surface_contact_manifest','bed_material'):options[setting['parameter']]=setting['value']
    # Ambient choices target a different native adapter. Do not silently claim
    # that a UI selection has changed this continuing body's thermal boundary.
    if selected.get('ambient_thermal') not in (None,'ambient-22c'):
        raise ValueError('Live environment air-temperature input is unavailable for this native body adapter')
    return selected,options


def box_contact(points, low, high, radius):
    """Sphere/AABB signed distance and outward normal, including interior points."""
    closest=np.clip(points,low,high);delta=points-closest;distance=np.linalg.norm(delta,axis=1)
    normal=delta/np.maximum(distance[:,None],1e-12)
    inside=np.all((points>=low)&(points<=high),axis=1)
    if np.any(inside):
        margins=np.concatenate((points[inside]-low,high-points[inside]),axis=1)
        face=np.argmin(margins,axis=1);normal[inside]=0
        rows=np.flatnonzero(inside);normal[rows,face%3]=np.where(face<3,-1.,1.)
        distance[inside]=-np.min(margins,axis=1)
    return np.maximum(0.,radius-distance),normal


def contact_force(depth, normal, relative, stiffness, damping, friction=.5):
    vn=np.sum(relative*normal,axis=1)
    magnitude=np.where(depth>0,np.maximum(0.,stiffness*depth-damping*vn),0.)
    tangent=relative-vn[:,None]*normal
    speed=np.linalg.norm(tangent,axis=1)
    return magnitude[:,None]*(normal-friction*tangent/np.maximum(speed[:,None],.05))


class SpringMesh:
    def __init__(self, ident, kind, low, high, mass):
        self.id=ident;self.kind=kind
        nx,ny,nz=(19,23,1) if kind=='cloth' else (11,9,3)
        low=np.array(low,float);high=np.array(high,float)
        axes=[np.linspace(low[i],high[i],n) for i,n in enumerate((nx,ny,nz))]
        if kind=='cloth':axes[2]=np.array([high[2]])
        self.x=np.stack(np.meshgrid(*axes,indexing='ij'),-1).reshape(-1,3)
        if kind=='soft':
            u=(self.x[:,0]-low[0])/(high[0]-low[0]);v=(self.x[:,1]-low[1])/(high[1]-low[1])
            layer=(self.x[:,2]-low[2])/(high[2]-low[2])
            dome=.28+.72*np.sqrt(np.maximum(0.,np.sin(np.pi*u)*np.sin(np.pi*v)))
            self.x[:,2]=low[2]+layer*(high[2]-low[2])*dome
        self.v=np.zeros_like(self.x);self.mass=mass/len(self.x);self.rest=self.x.copy()
        self.fixed=np.zeros(len(self.x),bool)
        if kind=='soft':self.fixed=self.x[:,2]==low[2]
        grid=np.arange(len(self.x)).reshape(nx,ny,nz);edges=[]
        # Axis, diagonal shear, and two-edge bending springs. Soft bodies also
        # include through-thickness diagonals to resist shearing/collapse.
        offsets=[(1,0,0),(0,1,0),(1,1,0),(1,-1,0),(2,0,0),(0,2,0)]
        if kind=='soft':offsets += [(0,0,1),(1,0,1),(0,1,1),(1,1,1),(1,-1,1)]
        for dx,dy,dz in offsets:
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        if 0<=i+dx<nx and 0<=j+dy<ny and k+dz<nz:
                            edges.append((grid[i,j,k],grid[i+dx,j+dy,k+dz],.08 if max(dx,dy)>1 else 1.))
        self.edges=np.array(edges);self.a=self.edges[:,0].astype(int);self.b=self.edges[:,1].astype(int)
        self.length=np.linalg.norm(self.x[self.b]-self.x[self.a],axis=1)
        self.stiffness=(45. if kind=='cloth' else 28.)*self.edges[:,2]
        faces=[]
        def quad(a,b,c,d):faces.extend(((a,b,c),(a,c,d)))
        for i in range(nx-1):
            for j in range(ny-1):
                quad(grid[i,j,-1],grid[i+1,j,-1],grid[i+1,j+1,-1],grid[i,j+1,-1])
                if nz>1:quad(grid[i,j,0],grid[i,j+1,0],grid[i+1,j+1,0],grid[i+1,j,0])
        if nz>1:
            for k in range(nz-1):
                for i in range(nx-1):
                    for j in (0,ny-1):quad(grid[i,j,k],grid[i+1,j,k],grid[i+1,j,k+1],grid[i,j,k+1])
                for j in range(ny-1):
                    for i in (0,nx-1):quad(grid[i,j,k],grid[i,j+1,k],grid[i,j+1,k+1],grid[i,j,k+1])
        self.indices=np.array(faces).reshape(-1).tolist()
        self.areal_density_kg_m2=mass/max(1e-12,float((high[0]-low[0])*(high[1]-low[1]))) if kind=='cloth' else None
        self.self_contact=None;self.self_contact_state={};self.skin_barrier_state={}
        self.cover_mode=False;self.stretch_constraint=None;self.stretch_state={};self.constraint_state={};self.preparation={}
        self.bending=None;self.bending_state={};self.triangle_contact=None;self.triangle_contact_state={};self.triangle_reference=None
        if kind=='cloth':
            from .cloth_contact import ClothSelfContact
            from .cloth_bending import ClothBending
            from .cloth_self_collision import ClothSelfCollision
            self.self_contact=ClothSelfContact(self)
            # Distance springs alone leave a crease free: an isometric fold costs
            # no energy, so folds were an unpenalised mesh mode rather than a
            # render artifact. The hinge term is the missing shell physics.
            faces=np.asarray(self.indices,int).reshape(-1,3)
            self.bending=ClothBending(self.rest,faces,rigidity_n_m=CLOTH_BENDING_RIGIDITY_N_M)
            self.triangle_contact=ClothSelfCollision(faces,len(self.x),thickness_m=CLOTH_THICKNESS_M)

    def spring_forces(self):
        d=self.x[self.b]-self.x[self.a];length=np.linalg.norm(d,axis=1)
        n=d/np.maximum(length[:,None],1e-12)
        rate=np.sum((self.v[self.b]-self.v[self.a])*n,axis=1)
        f=(self.stiffness*(length-self.length)+.035*rate)[:,None]*n
        total=np.zeros_like(self.x);np.add.at(total,self.a,f);np.add.at(total,self.b,-f)
        if self.bending is not None:total+=self.bending.forces(self.x)
        return total

    def elastic_energy(self):
        length=np.linalg.norm(self.x[self.b]-self.x[self.a],axis=1)
        stored=float(.5*np.sum(self.stiffness*(length-self.length)**2))
        return stored+(0. if self.bending is None else self.bending.energy(self.x))

    def frame(self):
        return {'id':self.id,'kind':self.kind,'positions':self.x.reshape(-1).tolist(),'indices':self.indices,'self_contact':deepcopy(self.self_contact_state),
                'bending':deepcopy(self.bending_state),'triangle_self_contact':deepcopy(self.triangle_contact_state),'fixed_nodes':np.flatnonzero(self.fixed).tolist(),'skin_barrier_last_substep':deepcopy(self.skin_barrier_state),'material_stretch':deepcopy(self.stretch_state),'constraint_last_substep':deepcopy(self.constraint_state),'preparation':deepcopy(self.preparation)}


def project_skin_barrier(mesh, tree, points, h, velocity=None):
    """Project shell overlaps; remove only inward relative normal velocity.

    Skin velocity is prescribed for this explicit environment substep. The
    opposite velocity impulse is returned to the body's force ports. Projection
    does not become velocity or a fabricated contact impulse. Spring-energy
    changes from projection are reported separately; gravitational projection
    work is unresolved because this API has no gravity argument.
    """
    if not np.isfinite(h) or h <= 0:
        raise ValueError('Positive finite substep required')
    points = np.asarray(points, dtype=float)
    velocity = np.zeros_like(points) if velocity is None else np.asarray(velocity, dtype=float)
    if velocity.shape != points.shape or not np.isfinite(velocity).all():
        raise ValueError('Skin velocity must match finite skin points')
    mass = np.broadcast_to(np.asarray(mesh.mass, dtype=float), (len(mesh.x),))
    impulses = np.zeros_like(points)
    before_kinetic = float(.5 * np.sum(mass[:, None] * mesh.v ** 2))

    def spring_energy():
        if not all(hasattr(mesh, attr) for attr in ('a', 'b', 'length', 'stiffness')):
            return None
        length = np.linalg.norm(mesh.x[mesh.b]-mesh.x[mesh.a], axis=1)
        return float(.5 * np.sum(mesh.stiffness * (length-mesh.length)**2))

    before_spring = spring_energy()
    maximum = 0.
    total = 0.
    boundary_work = 0.
    contact_nodes = set()
    for _ in range(3):
        distance, nearest = tree.query(mesh.x)
        delta = mesh.x - points[nearest]
        normal = delta / np.maximum(distance[:, None], 1e-12)
        coincident = distance < 1e-12
        if np.any(coincident):
            # A sample sphere supplies no normal at its exact centre. Prefer
            # boundary motion, then a deterministic axis; geometric ambiguity
            # remains an acknowledged limitation of this particle shell.
            fallback = velocity[nearest[coincident]].copy()
            norms = np.linalg.norm(fallback, axis=1)
            fallback[norms < 1e-12] = [0., 0., 1.]
            normal[coincident] = fallback / np.linalg.norm(fallback, axis=1)[:, None]
        depth = np.maximum(0., .048-distance)
        depth[mesh.fixed] = 0.
        active = depth > 0
        if not np.any(active):
            break
        maximum = max(maximum, float(depth.max()))
        total += float(depth.sum())
        contact_nodes.update(np.flatnonzero(active).tolist())
        mesh.x += depth[:, None] * normal
        relative_normal = np.sum((mesh.v-velocity[nearest])*normal, axis=1)
        dv = np.where(active, np.maximum(0., -relative_normal), 0.)[:, None] * normal
        particle_impulse = mass[:, None]*dv
        mesh.v += dv
        np.add.at(impulses, nearest, -particle_impulse)
        boundary_work += float(np.sum(particle_impulse*velocity[nearest]))
    after_kinetic = float(.5 * np.sum(mass[:, None] * mesh.v ** 2))
    after_spring = spring_energy()
    kinetic_change = after_kinetic-before_kinetic
    mesh.skin_barrier_state = {
        'vertex_count': len(contact_nodes),
        'max_projection_m': maximum,
        'sum_projection_m': total,
        'spring_energy_change_j': None if before_spring is None else after_spring-before_spring,
        'kinetic_energy_change_j': kinetic_change,
        'prescribed_boundary_work_j': boundary_work,
        'normal_impulse_dissipation_j': boundary_work-kinetic_change,
        'scope': 'Discrete sample-sphere projection; spring projection work reported, gravitational work unresolved; prescribed boundary work is an explicit-coupling estimate, not measured native work.'}
    return impulses



def drape_initial_cloth(mesh, skin_points, plane=-.24):
    """Collision-conforming initial placement; flat fabric rest lengths are retained.

    This is an initial condition, not a claim of static equilibrium. Subsequent
    server steps relax the springs and transfer every contact reaction to the body.
    """
    tree=cKDTree(skin_points[:,:2]);radius=.048
    neighborhoods=tree.query_ball_point(mesh.x[:,:2],radius)
    for i,neighbors in enumerate(neighborhoods):
        height=plane+.006
        if neighbors:
            points=skin_points[neighbors];d2=np.sum((points[:,:2]-mesh.x[i,:2])**2,axis=1)
            height=max(height,float(np.max(points[:,2]+np.sqrt(np.maximum(0.,radius**2-d2)))))
        mesh.x[i,2]=height+.002



def top_support_heights(x, plane, supports=()):
    """Top surfaces available to a loose supine cover; fixed world owns reaction."""
    height=np.full(len(x),plane,dtype=float)
    for low,high,_ident in supports:
        inside=np.all((x[:,:2]>=np.asarray(low)[:2])&(x[:,:2]<=np.asarray(high)[:2]),axis=1)&(x[:,2]>=low[2]-.003)
        height[inside]=np.maximum(height[inside],high[2])
    return height


def constrain_top_cover(mesh,cover,h,gravity,plane,supports=(),passes=CLOTH_COVER_PASSES,swept_from=None):
    """Alternate material and one-sided contact constraints, auditing net work.

    This is restricted loose top-cover physics. It cannot model underside wraps.
    Position repair is nonconservative and is never converted to velocity.
    ``swept_from`` supplies the positions this substep started from so the
    triangle-level pass can run its continuous test; without it only discrete
    proximity is enforced and fast tunnelling is not detected.
    """
    x0=mesh.x.copy();v0=mesh.v.copy();mass=mesh.mass
    spring0=mesh.elastic_energy()
    impulses=np.zeros_like(cover.points);support_impulse=np.zeros(3);boundary_work=0.;dissipation=0.
    if mesh.triangle_contact is not None:
        # Vertex/triangle and edge/edge contact, including the swept test. The
        # particle guard below only keeps vertices apart, so two coarse faces
        # could still cross between their vertices.
        #
        # The pass is skipped only under an exact clearance argument: the
        # distance between any two primitives can shrink by at most twice the
        # largest vertex displacement since the last evaluation, so while
        # 2*travel stays below the measured clearance margin no contact is
        # possible. Folding the sheet shrinks that margin and the pass resumes
        # on every substep. This is a bound, not a heuristic sampling stride.
        reference=getattr(mesh,'triangle_reference',None)
        clearance=mesh.triangle_contact_state.get('residual_min_distance_m',0.) if mesh.triangle_contact_state else 0.
        travel=0. if reference is None else float(np.linalg.norm(mesh.x-reference,axis=1).max(initial=0.))
        if reference is not None and np.isfinite(clearance) and clearance-2*travel>mesh.triangle_contact.thickness:
            mesh.triangle_contact_state=dict(mesh.triangle_contact_state,evaluated=False,
                clearance_margin_m=clearance-2*travel,vertex_travel_since_evaluation_m=travel)
        else:
            mesh.triangle_contact_state=mesh.triangle_contact.resolve(mesh,previous_positions=swept_from)
            mesh.triangle_contact_state.update(evaluated=True,clearance_margin_m=0.,vertex_travel_since_evaluation_m=0.)
            mesh.triangle_reference=mesh.x.copy()
            support_impulse+=np.asarray(mesh.triangle_contact_state['fixed_support_impulse_ns'])
    for iteration in range(passes):
        mesh.stretch_state=mesh.stretch_constraint.resolve(mesh,gravity=gravity)
        support_impulse+=np.asarray(mesh.stretch_state['support_impulse_ns'])
        mesh.self_contact_state=mesh.self_contact.resolve(mesh)
        support_impulse+=np.asarray(mesh.self_contact_state['support_impulse_ns'])
        depth=np.maximum(0.,top_support_heights(mesh.x,plane,supports)+.003-mesh.x[:,2]);depth[mesh.fixed]=0
        active=depth>0;mesh.x[:,2]+=depth
        # Static support contact: passive normal and Coulomb tangential impulses.
        normal_dv=np.where(active,np.maximum(0.,-mesh.v[:,2]),0.)
        dv=np.zeros_like(mesh.v);dv[:,2]=normal_dv
        speed=np.linalg.norm(mesh.v[:,:2],axis=1)
        dv[:,:2]=-np.minimum(speed,.5*normal_dv)[:,None]*mesh.v[:,:2]/np.maximum(speed[:,None],1e-15)
        mesh.v+=dv;support_impulse-=mass*dv.sum(axis=0)
        impulses+=cover.project(mesh,h,friction=.45)
        boundary_work+=mesh.skin_barrier_state['prescribed_boundary_work_j']
        dissipation+=mesh.skin_barrier_state['contact_impulse_dissipation_j']
        solver=mesh.stretch_constraint
        ratios=np.linalg.norm(mesh.x[solver.b]-mesh.x[solver.a],axis=1)/solver.rest
        if np.all(ratios<=solver.ratio+solver.tolerance):break
    solver=mesh.stretch_constraint
    ratios=np.linalg.norm(mesh.x[solver.b]-mesh.x[solver.a],axis=1)/solver.rest
    mesh.stretch_state.update(max_extension_ratio=float(ratios.max(initial=0.)),
        converged=bool(np.all(ratios<=solver.ratio+solver.tolerance)),
        violating_edges=int(np.count_nonzero(ratios>solver.ratio+solver.tolerance)))
    spring_change=mesh.elastic_energy()-spring0
    gravity_change=float(-np.sum(mass*(mesh.x-x0)*gravity))
    kinetic_change=float(.5*mass*np.sum(mesh.v**2-v0**2))
    mesh.constraint_state=dict(passes=iteration+1,spring_projection_energy_change_j=spring_change,
        gravity_projection_energy_change_j=gravity_change,kinetic_energy_change_j=kinetic_change,
        prescribed_skin_work_j=boundary_work,skin_impulse_dissipation_j=dissipation,
        unresolved_energy_change_j=spring_change+gravity_change+kinetic_change-boundary_work,
        fixed_world_support_impulse_ns=support_impulse.tolist(),
        triangle_self_contact=deepcopy(mesh.triangle_contact_state),
        scope='Loose supine top cover, no tucks or underside wrapping; material, hinge-bending, triangle self-contact and sampled upper-envelope constraints. Net nonconservative projection work reported, not declared conserved.')
    return impulses


def prepare_cloth(mesh, skin_points, plane=-.24, object_supports=()):
    """Drop a free flat sheet over held skin before clock zero.

    This sets a loose-blanket initial condition, not a certified equilibrium.
    The old mattress tucks conflict with the one-sided top-cover topology.
    """
    from .cloth_cover import SupineClothCover
    from .cloth_stretch import ClothStretchConstraint
    mesh.cover_mode=True;mesh.fixed[:]=False;mesh.triangle_reference=None;mesh.triangle_contact_state={}
    # Two sweeps left the strain limit unconverged in the runtime (1.72 against
    # a 1.12 ceiling). The projection is cheap; the sweep count is not the place
    # to save time, and the constraint is reported converged or not either way.
    mesh.stretch_constraint=ClothStretchConstraint(mesh,iterations=CLOTH_STRETCH_SWEEPS)
    mesh.x=mesh.rest.copy();mesh.x[:,2]=float(np.max(skin_points[:,2]))+.052;mesh.v[:]=0
    cover=SupineClothCover(skin_points);h=.001;gravity=np.array([0.,0.,-9.81])
    projection_work=0.
    for _ in range(1500):
        force=mesh.spring_forces()+mesh.mass*gravity
        depth,normal,nearest=cover.query(mesh.x)
        force+=contact_force(depth*normal[:,2],normal,mesh.v,60.,.12,.45)
        normal=np.zeros_like(mesh.x);normal[:,2]=1
        depth=np.maximum(0.,top_support_heights(mesh.x,plane,object_supports)+.003-mesh.x[:,2])
        force+=contact_force(depth,normal,mesh.v,100.,.15)
        mesh.v+=force/mesh.mass*h;mesh.v*=np.exp(-3*h);swept=mesh.x.copy();mesh.x+=mesh.v*h
        constrain_top_cover(mesh,cover,h,gravity,plane,object_supports,passes=3,swept_from=swept)
        projection_work+=mesh.constraint_state['unresolved_energy_change_j']
    reset=float(.5*mesh.mass*np.sum(mesh.v**2));mesh.v[:]=0
    mesh.preparation=dict(duration_s=1.5,held_body=True,equilibrium_accepted=False,
        boundary='Loose flat sheet dropped above body; no fixed tucks',kinetic_reset_j=reset,
        unresolved_constraint_energy_change_j=projection_work)



class RigidProp:
    def __init__(self, ident, record, offset):
        self.id=ident;self.kind='rigid';self.radius=record.get('radius_m')
        low=np.array(record['bounds_m']['min']);high=np.array(record['bounds_m']['max'])
        self.origin=(low+high)/2;self.x=self.origin+offset;self.v=np.zeros(3);self.omega=np.zeros(3);self.rotation=np.eye(3)
        self.half=(high-low)/2;self.mass=float(record['mass_kg'])
        from .rigid_contact import BOUNCY_BALL
        self.material=dict(BOUNCY_BALL if self.radius else {'name':'engineering-block','restitution':0.,'friction':.5,'bounce_threshold_m_s':.05})
        self.material.update(record.get('contact_material',{}))
        if not 0<=self.material['restitution']<=1 or not 0<=self.material['friction']:raise ValueError('Invalid rigid contact material')
        self.contact_dissipation_j=0.;self.contact_projection_m=0.
        self.inertia=np.full(3,.4*self.mass*self.radius**2) if self.radius else self.mass/3*(np.sum(self.half**2)-self.half**2)

    def integrate(self,dt,force,torque,*,plane=None,axis=2):
        from .mechanics import _rotation
        if self.radius and plane is not None:
            from .rigid_contact import sphere_plane_step
            loss,_,repair=sphere_plane_step(self.x,self.v,self.omega,dt,force/self.mass,self.radius,self.mass,axis,plane,
                self.material['restitution'],self.material['friction'],self.material['bounce_threshold_m_s'])
            self.contact_dissipation_j+=loss;self.contact_projection_m+=repair
        else:
            self.x+=self.v*dt+.5*force/self.mass*dt**2;self.v+=force/self.mass*dt
        local=self.rotation.T@self.omega
        self.omega+=self.rotation@((self.rotation.T@torque-np.cross(local,self.inertia*local))/self.inertia)*dt
        self.rotation=_rotation(self.omega*dt)@self.rotation

    def frame(self):
        return {'id':self.id,'kind':'rigid','position_m':self.x.tolist(),'origin_m':self.origin.tolist(),'rotation_matrix':self.rotation.tolist(), 'velocity_m_s':self.v.tolist(), 'contact_material':dict(self.material), 'contact_dissipation_j':self.contact_dissipation_j, 'contact_projection_m':self.contact_projection_m}


class EnvironmentDynamics:
    def __init__(self, root, environment, selection, registration):
        self.root=Path(root);self.selection=deepcopy(selection);self.time_s=0.;self.contacts=[]
        self.gravity=np.array({'supine':[0,0,-9.81],'upright':[0,-9.81,0],'free':[0,0,0]}[environment],float)
        self.axis=2 if environment=='supine' else 1
        self.base_plane=-.24 if environment=='supine' else -.96
        from .world_frame import GravityAlignedWorldFrame
        self.native_frame_bound=hasattr(registration,'global_map') and hasattr(registration,'reference')
        if self.native_frame_bound:
            self.world_frame=GravityAlignedWorldFrame.from_native(registration,registration.reference,environment=environment,template_support_height_m=self.base_plane)
            self.gravity=self.world_frame.world_gravity_m_s2.copy()
        else:
            # Explicit held-canonical fixtures have no native gravity/support owner.
            self.world_frame=GravityAlignedWorldFrame(np.eye(4),source_to_world_rotation=np.eye(3),source_gravity_m_s2=self.gravity)

        self.static=[];self.object_static=[];self.soft=[];self.rigid=[];self.sources={}
        catalog=self.read(CATALOGUE);records={r['id']:r for r in catalog['objects']}
        scene=next((s for s in catalog['scenes'] if s['id']==selection.get('scene')),None)
        placements=deepcopy(scene['placements']) if scene else []
        # A bare bed now has an actual frame, cushion and blanket as well.
        if not scene and environment=='supine':placements=[{'object':i,'offset_m':[0,0,0]} for i in ('bed-frame','bed-mattress','pillow','blanket')]
        for ident in selection.get('objects',[]):placements.append({'object':ident,'offset_m':[0,0,0]})
        self.placements=[];counts={}
        for placement in placements:
            ident=placement['object'];counts[ident]=counts.get(ident,0)+1
            instance=ident+'-'+str(counts[ident]);offset=np.array(placement.get('offset_m',[0,0,0]),float)
            record=records[ident];geometry=self.read(record['geometry'])
            self.placements.append({'id':ident,'instance':instance,'offset_m':offset.tolist()})
            if ident=='blanket' or ident=='pillow':
                bounds=record['bounds_m'];self.soft.append(SpringMesh(instance,'cloth' if ident=='blanket' else 'soft',np.array(bounds['min'])+offset,np.array(bounds['max'])+offset,record['mass_kg']))
            elif ident in ('ball-small','ball-large','block'):self.rigid.append(RigidProp(instance,record,offset))
            elif ident=='bed-mattress':
                bounds=record['bounds_m'];self.object_static.append((np.asarray(bounds['min'])+offset,np.asarray(bounds['max'])+offset,instance))
            elif ident!='bed-mattress':
                for part in geometry['parts']:self.add_part(part,offset,instance)
        if scene:
            world=scene['world'];geometry=self.read(world['surround_geometry'])
            for part in geometry['parts']:
                part=deepcopy(part)
                # Surround part metadata is gravity-frame, its mesh is canonical.
                if environment=='supine':
                    for field in ('min_m','max_m'):part[field]=np.array(part[field])[[2,0,1]].tolist()
                self.add_part(part,np.zeros(3),'world:'+part['name'])
        # The native body owns its bed support; cloth/props separately contact
        # the bounded mattress and actual room floor, never an infinite bed top.
        floors=[high[self.axis] for low,high,ident in self.static if ident in ('world:floor','world:sheet')]
        self.object_plane=max(floors) if floors else self.base_plane
        # Contact sampling in canonical skin coordinates, attached to named bones.
        skin=self.read('data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz')
        points=np.array(skin['positions']).reshape(-1,3)
        _,indices=np.unique(np.floor(points/.045).astype(int),axis=0,return_index=True);points=points[indices]
        self.ids=[];self.offsets=[]
        for p in points:
            owner=registration._ranking(p)[0][1]
            ident=registration.groups[owner]['canonical_bones'][0]
            self.ids.append(ident);self.offsets.append(p-np.array(registration.specs[ident]['centroid_m']))
        self.offsets=np.array(self.offsets);self.previous=None;self.last_impulse=np.zeros(3)
        for mesh in self.soft:
            if mesh.kind=='cloth':prepare_cloth(mesh,self.world_frame.points_to_world(points),self.object_plane,self.static+self.object_static)

    def read(self, relative):
        path=(self.root/relative).resolve()
        if not path.is_relative_to(self.root):raise ValueError('Environment asset escaped workspace')
        raw=path.read_bytes();self.sources[relative]=hashlib.sha256(raw).hexdigest()
        return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)

    def add_part(self,part,offset,ident):
        if 'min_m' in part:low=np.array(part['min_m'])+offset;high=np.array(part['max_m'])+offset
        else:
            center=np.array(part['centre_m'])+offset;half=np.full(3,part['radius_m'])
            if part['primitive']=='cylinder':half[part.get('axis',1)]=part['height_m']/2
            low=center-half;high=center+half
        self.static.append((low,high,ident))

    def skin_points(self,entities):
        points=np.array([np.array(entities[i]['centroid_m'])+np.array(entities[i]['rotation_matrix'])@p for i,p in zip(self.ids,self.offsets)])
        frame=getattr(self,'world_frame',None)
        return points if frame is None else frame.points_to_world(points)

    @staticmethod
    def prop_contact(prop,points,radius):
        if prop.radius:
            delta=points-prop.x;distance=np.linalg.norm(delta,axis=1)
            normal=delta/np.maximum(distance[:,None],1e-12)
            normal[distance<1e-12]=[1,0,0]
            return np.maximum(0.,prop.radius+radius-distance),normal
        local=(points-prop.x)@prop.rotation
        depth,normal=box_contact(local,-prop.half,prop.half,radius)
        return depth,normal@prop.rotation.T

    def owns(self, ident):
        return isinstance(ident,str) and any(o.id==ident for o in self.soft+self.rigid)

    def validate_object_forces(self, forces, *, _world_coordinates=False):
        """Bind current world force points to material coordinates before mutation."""
        if not isinstance(forces,(list,tuple)) or len(forces)>32:raise ValueError('At most 32 object force ports required')
        objects={o.id:o for o in self.soft+self.rigid};result=[]
        def vector(value):
            if not isinstance(value,(list,tuple,np.ndarray)) or len(value)!=3 or any(isinstance(v,(bool,np.bool_)) for v in value):raise ValueError('Object force vectors require three finite numbers')
            value=np.asarray(value,dtype=float)
            if value.shape!=(3,) or not np.isfinite(value).all():raise ValueError('Object force vectors require three finite numbers')
            return value
        for port in forces:
            if not isinstance(port,dict) or set(port)!={'id','force_n','point_m'}:raise ValueError('Object force requires id, force_n and point_m')
            if not self.owns(port['id']):raise ValueError('Unknown environment object force owner')
            owner=objects[port['id']];force=vector(port['force_n']);point=vector(port['point_m'])
            frame=getattr(self,'world_frame',None)
            if frame is not None and not _world_coordinates:
                force=frame.vectors_to_world(force);point=frame.points_to_world(point)
            if np.linalg.norm(force)>100:raise ValueError('Object force exceeds 100 N envelope')
            entry={'id':owner.id,'force_n':force.tolist()}
            if owner.kind=='rigid':
                local=owner.rotation.T@(point-owner.x)
                accepted=np.linalg.norm(local)<=owner.radius+.02 if owner.radius else np.all(np.abs(local)<=owner.half+.02)
                if not accepted:raise ValueError('Object force point is stale or outside its material owner')
                entry['local_point_m']=local.tolist()
            else:
                distance=np.linalg.norm(owner.x-point,axis=1);index=int(np.argmin(distance))
                if distance[index]>max(.025,.75*float(np.median(owner.length))):raise ValueError('Cloth force point is stale or outside its material owner')
                if owner.fixed[index]:raise ValueError('Cannot force an anchored cloth vertex')
                entry['node_index']=index
            result.append(entry)
        return result

    def checkpoint(self):
        return deepcopy((self.soft,self.rigid,self.previous,self.time_s,self.contacts,self.last_impulse,getattr(self,'last_object_forces',[])))

    def restore(self,state):
        self.soft,self.rigid,self.previous,self.time_s,self.contacts,self.last_impulse,self.last_object_forces=deepcopy(state)

    def advance(self,dt,entities,object_forces=()):
        # Revalidate material entries before changing any state. The public
        # runtime binds world points once; each substep follows those materials.
        if not isinstance(object_forces,(list,tuple)) or len(object_forces)>32:raise ValueError('At most 32 object force ports required')
        objects={o.id:o for o in self.soft+self.rigid};raw=[]
        for port in object_forces:
            if not isinstance(port,dict) or not self.owns(port.get('id')):raise ValueError('Unknown object material force')
            obj=objects[port['id']]
            if obj.kind=='rigid':
                if set(port)!={'id','force_n','local_point_m'}:raise ValueError('Rigid force requires material point')
                local=np.asarray(port['local_point_m'],dtype=float)
                if local.shape!=(3,) or not np.isfinite(local).all():raise ValueError('Invalid rigid material point')
                point=obj.x+obj.rotation@local
            else:
                if set(port)!={'id','force_n','node_index'}:raise ValueError('Cloth force requires material node')
                node=port['node_index']
                if type(node) is not int or not 0<=node<len(obj.x):raise ValueError('Invalid cloth material node')
                if obj.fixed[node]:raise ValueError('Cannot force an anchored cloth vertex')
                point=obj.x[node]
            raw.append({'id':obj.id,'force_n':port['force_n'],'point_m':point.tolist()})
        # Revalidation is for force bounds and geometry only; retain the exact
        # original cloth node even when two folded vertices occupy one position.
        self.validate_object_forces(raw,_world_coordinates=True)
        points=self.skin_points(entities);velocity=np.zeros_like(points) if self.previous is None else (points-self.previous)/dt
        tree=cKDTree(points);impulses=np.zeros_like(points);h=dt/int(np.ceil(dt/.001));count=round(dt/h)
        object_plane=getattr(self,'object_plane',self.base_plane)
        object_supports=self.static+getattr(self,'object_static',[])
        cover=None
        if any(mesh.cover_mode for mesh in self.soft):
            from .cloth_cover import SupineClothCover
            cover=SupineClothCover(points,velocity)
        self.contacts=[]
        # Fixed walls/furniture act directly on the body; native support planes
        # retain their existing owner, avoiding duplicate body floor support.
        for low,high,ident in self.static:
            if ident in ('world:floor','world:sheet') and abs(high[self.axis]-self.base_plane)<1e-6:continue
            depth,normal=box_contact(points,low,high,.012)
            force=contact_force(depth,normal,velocity,250.,1.)
            impulses+=force*dt
        for _ in range(count):
            rigid_loads={p.id:[np.zeros(3),np.zeros(3)] for p in self.rigid}
            for prop in self.rigid:
                for port in object_forces:
                    if port['id']!=prop.id:continue
                    force=np.asarray(port['force_n']);lever=prop.rotation@np.asarray(port['local_point_m'])
                    rigid_loads[prop.id][0]+=force;rigid_loads[prop.id][1]+=np.cross(lever,force)
            # Contact manifolds exchange forces at the same world stations,
            # preserving linear and angular momentum between free props.
            for i,a in enumerate(self.rigid):
                for b in self.rigid[i+1:]:
                    stations=a.x[None,:] if a.radius else np.array([[x,y,z] for x in (-1,1) for y in (-1,1) for z in (-1,1)])*a.half@a.rotation.T+a.x
                    depth,n=self.prop_contact(b,stations,a.radius or 0.)
                    contact=stations-n*(a.radius or 0.)
                    relative=a.v+np.cross(a.omega,contact-a.x)-b.v-np.cross(b.omega,contact-b.x)
                    f=contact_force(depth,n,relative,1500.,8.)
                    rigid_loads[a.id][0]+=f.sum(axis=0);rigid_loads[b.id][0]-=f.sum(axis=0)
                    rigid_loads[a.id][1]+=np.cross(contact-a.x,f).sum(axis=0)
                    rigid_loads[b.id][1]-=np.cross(contact-b.x,f).sum(axis=0)
            for mesh in self.soft:
                f=mesh.spring_forces()+mesh.mass*self.gravity
                for port in object_forces:
                    if port['id']==mesh.id:f[port['node_index']]+=np.asarray(port['force_n'])
                if mesh.cover_mode:
                    depth,normal,nearest=cover.query(mesh.x)
                    relative=mesh.v.copy();covered=nearest>=0;relative[covered]-=velocity[nearest[covered]]
                    force=contact_force(depth*normal[:,2],normal,relative,60.,.12,.45)
                    nearest=np.maximum(nearest,0)
                else:
                    distance,nearest=tree.query(mesh.x)
                    delta=mesh.x-points[nearest];normal=delta/np.maximum(distance[:,None],1e-12)
                    force=contact_force(np.maximum(0.,.048-distance),normal,mesh.v-velocity[nearest],60.,.12,.45)
                force[mesh.fixed]=0
                f+=force;np.add.at(impulses,nearest,-force*h)
                for prop in self.rigid:
                    depth,n=self.prop_contact(prop,mesh.x,.003)
                    relative=mesh.v-prop.v-np.cross(prop.omega,mesh.x-prop.x)
                    contact=contact_force(depth,n,relative,80.,.15)
                    contact[mesh.fixed]=0
                    f+=contact;rigid_loads[prop.id][0]-=contact.sum(axis=0)
                    rigid_loads[prop.id][1]-=np.cross(mesh.x-prop.x,contact).sum(axis=0)
                for low,high,_ident in object_supports:
                    depth,n=box_contact(mesh.x,low,high,.003)
                    f+=contact_force(depth,n,mesh.v,80.,.15)
                if np.any(self.gravity):
                    depth=np.maximum(0.,object_plane+.003-mesh.x[:,self.axis]);normal=np.zeros_like(mesh.x);normal[:,self.axis]=1
                    f+=contact_force(depth,normal,mesh.v,100.,.15)
                mesh.v+=f/mesh.mass*h;mesh.v*=np.exp(-.3*h);swept=mesh.x.copy();mesh.x+=mesh.v*h
                mesh.x[mesh.fixed]=mesh.rest[mesh.fixed];mesh.v[mesh.fixed]=0
                if mesh.cover_mode:
                    impulses+=constrain_top_cover(mesh,cover,h,self.gravity,object_plane,object_supports,swept_from=swept)
                elif mesh.kind=='cloth':
                    mesh.self_contact_state=mesh.self_contact.resolve(mesh)
                    impulses+=project_skin_barrier(mesh,tree,points,h,velocity=velocity)
            for prop in self.rigid:
                force=prop.mass*self.gravity.copy()+rigid_loads[prop.id][0];torque=rigid_loads[prop.id][1].copy()
                local=(points-prop.x)@prop.rotation
                if prop.radius:
                    distance=np.linalg.norm(local,axis=1);normal=(local/np.maximum(distance[:,None],1e-12))@prop.rotation.T
                    depth=np.maximum(0.,prop.radius+.012-distance)
                else:
                    depth,n=box_contact(local,-prop.half,prop.half,.012);normal=n@prop.rotation.T
                rel=velocity-(prop.v+np.cross(prop.omega,points-prop.x))
                body_force=contact_force(depth,normal,rel,250.,1.,.45)
                impulses+=body_force*h;force-=np.sum(body_force,axis=0);torque-=np.sum(np.cross(points-prop.x,body_force),axis=0)
                corners=np.array([[x,y,z] for x in (-1,1) for y in (-1,1) for z in (-1,1)])*prop.half
                if prop.radius:
                    corners=np.concatenate((np.eye(3),-np.eye(3)))*prop.radius
                support=(corners if prop.radius else corners@prop.rotation.T)+prop.x
                sv=prop.v+np.cross(prop.omega,support-prop.x)
                supports=np.zeros_like(support)
                for low,high,_ident in object_supports:
                    if prop.radius and abs(high[self.axis]-object_plane)<1e-6 and _ident in ('world:floor','world:sheet'):continue
                    depth,n=box_contact(support,low,high,0.)
                    supports+=contact_force(depth,n,sv,1500.,8.)
                if np.any(self.gravity) and not prop.radius:
                    n=np.zeros_like(support);n[:,self.axis]=1
                    supports+=contact_force(np.maximum(0.,object_plane-support[:,self.axis]),n,sv,1500.,8.)
                force+=np.sum(supports,axis=0);torque+=np.sum(np.cross(support-prop.x,supports),axis=0)
                prop.integrate(h,force,torque,plane=object_plane if np.any(self.gravity) else None,axis=self.axis)
        for mesh in self.soft:
            if mesh.bending is not None:mesh.bending_state=mesh.bending.state(mesh.x,mesh.areal_density_kg_m2)
            if not np.isfinite(mesh.x).all() or np.max(np.linalg.norm(mesh.v,axis=1))>30:raise RuntimeError('Environment deformation exceeded integration envelope')
        for prop in self.rigid:
            if not np.isfinite(prop.x).all() or np.linalg.norm(prop.v)>30:raise RuntimeError('Rigid environment exceeded integration envelope')
        forces=impulses/dt;active=np.flatnonzero(np.linalg.norm(forces,axis=1)>1e-9)
        frame=getattr(self,'world_frame',None)
        output_points=points if frame is None else frame.points_to_canonical(points)
        output_forces=forces if frame is None else frame.vectors_to_canonical(forces)
        ports=[{'id':self.ids[i],'point_m':output_points[i].tolist(),'force_n':output_forces[i].tolist()} for i in active]
        if any(np.linalg.norm(forces[i])>1000 for i in active):raise RuntimeError('Environment contact exceeded force envelope')
        self.last_object_forces=deepcopy(list(object_forces))
        self.contacts=ports;self.last_impulse=np.sum(output_forces,axis=0)*dt;self.previous=points.copy();self.time_s+=dt
        return ports

    def frame(self):
        return {'schema':'ihm.environment-state.v1','time_s':self.time_s,'selection':self.selection,'scope':SCOPE,
                'placements':self.placements,'objects':[o.frame() for o in self.soft+self.rigid],
                'world_frame':None if not hasattr(self,'world_frame') else self.world_frame.metadata(),
                'native_frame_bound':getattr(self,'native_frame_bound',False),
                'object_coordinate_frame':'gravity_aligned_world','body_force_coordinate_frame':'canonical',
                'contact_count':len(self.contacts),'body_impulse_ns':self.last_impulse.tolist(),
                'body_force_ports':deepcopy(self.contacts),'object_force_ports':deepcopy(getattr(self,'last_object_forces',[])),'source_sha256':self.sources,
                'parameters':{'substep_s':.001,'skin_sample_spacing_m':.045,'contact_radius_m':.048,
                              'cloth_initialization':'Loose flat top cover dropped for 1.5 s against held skin, zeroed initial velocity; no equilibrium acceptance',
                              'object_support':'Bounded mattress/furniture and room floor, separate from native body support',
                              'cloth_max_extension_ratio':1.12,
                              'cloth_self_contact_separation_m':.012,'cloth_bending_rigidity_n_m':CLOTH_BENDING_RIGIDITY_N_M,
                              'cloth_thickness_m':CLOTH_THICKNESS_M,'cloth_stretch_sweeps':CLOTH_STRETCH_SWEEPS,'cloth_cover_passes':CLOTH_COVER_PASSES,'ball_restitution':.75,'ball_friction':.35,'cloth_spring_n_m':45.,'pillow_spring_n_m':28.,'coefficient_basis':'uncalibrated engineering values'}}
