"""Server-owned scene mechanics, coupled through canonical material force ports.

Engineering models, not calibrated tissue/contact measurements. Skin samples are
rigidly attached to native bone owners; they are not an independent human mass.
Cloth uses stretch/shear/bend springs, cushions a three-dimensional spring lattice,
and free props use six-degree-of-freedom Newton/Euler integration. Contacts use
compliant normal forces and regularized Coulomb friction. Exchange is explicit
at 20 ms, with 1 ms environment substeps and equal/opposite contact impulses.
"""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

CATALOGUE = 'data/derived/environment-catalogue-v1/catalogue.json'
SCOPE = ('Server-owned spring cloth, volumetric spring-lattice pillow, rigid props and fixed compound colliders. '
         'Sampled skin contact loads act on the native articulated body at material points; reciprocal impulses act on props. '
         'Engineering coefficients, approximate skin registration, fixed cylinder/canopy bounding boxes and vertex rigid manifolds, 20 ms explicit exchange; no cloth self-contact, '
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

    def spring_forces(self):
        d=self.x[self.b]-self.x[self.a];length=np.linalg.norm(d,axis=1)
        n=d/np.maximum(length[:,None],1e-12)
        rate=np.sum((self.v[self.b]-self.v[self.a])*n,axis=1)
        f=(self.stiffness*(length-self.length)+.035*rate)[:,None]*n
        total=np.zeros_like(self.x);np.add.at(total,self.a,f);np.add.at(total,self.b,-f)
        return total

    def frame(self):
        return {'id':self.id,'kind':self.kind,'positions':self.x.reshape(-1).tolist(),'indices':self.indices}


def project_skin_barrier(mesh,tree,points,h):
    """Hard shell guard supplements compliant force contact for thin cloth.

    Position corrections become impulses, including their reciprocal body load.
    This prevents stiff textile springs pulling vertices through the sampled skin.
    """
    impulses=np.zeros_like(points)
    for _ in range(3):
        distance,nearest=tree.query(mesh.x);delta=mesh.x-points[nearest]
        normal=delta/np.maximum(distance[:,None],1e-12)
        correction=np.maximum(0.,.048-distance)[:,None]*normal
        correction[mesh.fixed]=0
        mesh.x+=correction;mesh.v+=correction/h
        np.add.at(impulses,nearest,-mesh.mass*correction/h)
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



def prepare_cloth(mesh, skin_points, plane=-.24):
    """Precondition fabric against a held reference body before clock zero.

    This dissipative setup avoids starting a stretched drape with a large impulse.
    It is separate from the continuing body's clock and native force exchange.
    """
    drape_initial_cloth(mesh,skin_points,plane)
    # The foot hem and two side corners are tucked into the mattress: an explicit fixed boundary,
    # rather than a visually placed panel that immediately slides off the bed.
    foot=mesh.rest[:,1]==mesh.rest[:,1].min()
    side_corners=(mesh.rest[:,1]==mesh.rest[:,1].max())&((mesh.rest[:,0]==mesh.rest[:,0].min())|(mesh.rest[:,0]==mesh.rest[:,0].max()))
    mesh.fixed=foot|side_corners
    mesh.x[mesh.fixed,2]=plane+.006;mesh.rest[mesh.fixed]=mesh.x[mesh.fixed]
    tree=cKDTree(skin_points);h=.001
    for _ in range(1500):
        force=mesh.spring_forces();force[:,2]-=mesh.mass*9.81
        distance,nearest=tree.query(mesh.x);delta=mesh.x-skin_points[nearest]
        normal=delta/np.maximum(distance[:,None],1e-12)
        force+=contact_force(np.maximum(0.,.048-distance),normal,mesh.v,60.,.12,.45)
        normal=np.zeros_like(mesh.x);normal[:,2]=1
        force+=contact_force(np.maximum(0.,plane+.003-mesh.x[:,2]),normal,mesh.v,100.,.15)
        mesh.v+=force/mesh.mass*h;mesh.v*=np.exp(-3*h);mesh.x+=mesh.v*h
        project_skin_barrier(mesh,tree,skin_points,h)
        mesh.x[mesh.fixed]=mesh.rest[mesh.fixed];mesh.v[mesh.fixed]=0
    mesh.v[:]=0



class RigidProp:
    def __init__(self, ident, record, offset):
        self.id=ident;self.kind='rigid';self.radius=record.get('radius_m')
        low=np.array(record['bounds_m']['min']);high=np.array(record['bounds_m']['max'])
        self.origin=(low+high)/2;self.x=self.origin+offset;self.v=np.zeros(3);self.omega=np.zeros(3);self.rotation=np.eye(3)
        self.half=(high-low)/2;self.mass=float(record['mass_kg'])
        self.inertia=np.full(3,.4*self.mass*self.radius**2) if self.radius else self.mass/3*(np.sum(self.half**2)-self.half**2)

    def integrate(self,dt,force,torque):
        from .mechanics import _rotation
        self.v+=force/self.mass*dt;self.x+=self.v*dt
        local=self.rotation.T@self.omega
        self.omega+=self.rotation@((self.rotation.T@torque-np.cross(local,self.inertia*local))/self.inertia)*dt
        self.rotation=_rotation(self.omega*dt)@self.rotation

    def frame(self):
        return {'id':self.id,'kind':'rigid','position_m':self.x.tolist(),'origin_m':self.origin.tolist(),'rotation_matrix':self.rotation.tolist()}


class EnvironmentDynamics:
    def __init__(self, root, environment, selection, registration):
        self.root=Path(root);self.selection=deepcopy(selection);self.time_s=0.;self.contacts=[]
        self.gravity=np.array({'supine':[0,0,-9.81],'upright':[0,-9.81,0],'free':[0,0,0]}[environment],float)
        self.axis=2 if environment=='supine' else 1
        self.base_plane=-.24 if environment=='supine' else -.96
        self.static=[];self.soft=[];self.rigid=[];self.sources={}
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
            if mesh.kind=='cloth':prepare_cloth(mesh,points,self.base_plane)

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
        return np.array([np.array(entities[i]['centroid_m'])+np.array(entities[i]['rotation_matrix'])@p for i,p in zip(self.ids,self.offsets)])

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

    def checkpoint(self):
        return deepcopy((self.soft,self.rigid,self.previous,self.time_s,self.contacts,self.last_impulse))

    def restore(self,state):
        self.soft,self.rigid,self.previous,self.time_s,self.contacts,self.last_impulse=deepcopy(state)

    def advance(self,dt,entities):
        points=self.skin_points(entities);velocity=np.zeros_like(points) if self.previous is None else (points-self.previous)/dt
        tree=cKDTree(points);impulses=np.zeros_like(points);h=dt/int(np.ceil(dt/.001));count=round(dt/h)
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
                for low,high,_ident in self.static:
                    depth,n=box_contact(mesh.x,low,high,.003)
                    f+=contact_force(depth,n,mesh.v,80.,.15)
                if np.any(self.gravity):
                    depth=np.maximum(0.,self.base_plane+.003-mesh.x[:,self.axis]);normal=np.zeros_like(mesh.x);normal[:,self.axis]=1
                    f+=contact_force(depth,normal,mesh.v,100.,.15)
                mesh.v+=f/mesh.mass*h;mesh.v*=np.exp(-.3*h);mesh.x+=mesh.v*h
                mesh.x[mesh.fixed]=mesh.rest[mesh.fixed];mesh.v[mesh.fixed]=0
                if mesh.kind=='cloth':impulses+=project_skin_barrier(mesh,tree,points,h)
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
                for low,high,_ident in self.static:
                    depth,n=box_contact(support,low,high,0.)
                    supports+=contact_force(depth,n,sv,1500.,8.)
                if np.any(self.gravity):
                    n=np.zeros_like(support);n[:,self.axis]=1
                    supports+=contact_force(np.maximum(0.,self.base_plane-support[:,self.axis]),n,sv,1500.,8.)
                force+=np.sum(supports,axis=0);torque+=np.sum(np.cross(support-prop.x,supports),axis=0)
                prop.integrate(h,force,torque)
        for mesh in self.soft:
            if not np.isfinite(mesh.x).all() or np.max(np.linalg.norm(mesh.v,axis=1))>30:raise RuntimeError('Environment deformation exceeded integration envelope')
        for prop in self.rigid:
            if not np.isfinite(prop.x).all() or np.linalg.norm(prop.v)>30:raise RuntimeError('Rigid environment exceeded integration envelope')
        forces=impulses/dt;active=np.flatnonzero(np.linalg.norm(forces,axis=1)>1e-9)
        ports=[{'id':self.ids[i],'point_m':points[i].tolist(),'force_n':forces[i].tolist()} for i in active]
        if any(np.linalg.norm(forces[i])>1000 for i in active):raise RuntimeError('Environment contact exceeded force envelope')
        self.contacts=ports;self.last_impulse=np.sum(impulses,axis=0);self.previous=points.copy();self.time_s+=dt
        return ports

    def frame(self):
        return {'schema':'ihm.environment-state.v1','time_s':self.time_s,'selection':self.selection,'scope':SCOPE,
                'placements':self.placements,'objects':[o.frame() for o in self.soft+self.rigid],
                'contact_count':len(self.contacts),'body_impulse_ns':self.last_impulse.tolist(),
                'body_force_ports':deepcopy(self.contacts),'source_sha256':self.sources,
                'parameters':{'substep_s':.001,'skin_sample_spacing_m':.045,'contact_radius_m':.048,
                              'cloth_initialization':'Foot hem and two side corners tucked at the mattress; 1.5 s damped settling against held canonical skin before time zero',
                              'cloth_spring_n_m':45.,'pillow_spring_n_m':28.,'coefficient_basis':'uncalibrated engineering values'}}
