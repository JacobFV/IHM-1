"""Coupled cervical/thoracic source metric with a shared parent body twist.

Six independent neck coordinates expand into 24 source coordinates through
18 explicitly affine couplers. This is an opt-in material system, not a force
that can be added to the existing native torso without replacing its inertia.
"""
import hashlib,json,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from .cervical_registration import numbers,reference_joint
from .thoracic_free_dynamics import ThoracicFreeDynamics


def skew(x):return np.cross(x,np.eye(3)).T

def vee(a):return np.array([a[2,1]-a[1,2],a[0,2]-a[2,0],a[1,0]-a[0,1]])/2

def constant(t):return (np.asarray(t,float),np.zeros((6,4,4)),np.zeros((4,4)),np.zeros((4,4)))

def product(a,b):
    x,dx,xt,xtt=a;y,dy,yt,ytt=b
    return x@y,np.einsum('iab,bc->iac',dx,y)+np.einsum('ab,ibc->iac',x,dy),xt@y+x@yt,xtt@y+2*xt@yt+x@ytt


def affine(function):
    """Admit explicit affine functions only, never sample an arbitrary spline."""
    if function.tag=='Constant':
        c=numbers(function.findtext('value'))
        if c.shape==(1,):return 0.,float(c[0])
    if function.tag=='LinearFunction':
        c=numbers(function.findtext('coefficients'))
        if c.shape==(1,):return 0.,float(c[0])
        if c.shape==(2,):return tuple(c)
    if function.tag=='SimmSpline':
        x=numbers(function.findtext('x'));y=numbers(function.findtext('y'))
        if x.shape==y.shape==(2,) and x[1]>x[0]:
            slope=(y[1]-y[0])/(x[1]-x[0]);return float(slope),float(y[0]-slope*x[0])
    raise ValueError('Unsupported non-affine cervical function')


class CervicalKinematics:
    def __init__(self,recipe):
        self.recipe=recipe;self.names=list(recipe['coordinates']);self.constraints={}
        for row in recipe['constraints']:
            c=ET.fromstring(row['source_xml']);dep=c.findtext('dependent_coordinate_name');ind=(c.findtext('independent_coordinate_names') or '').split()
            functions=c.find('coupled_coordinates_function')
            if c.tag!='CoordinateCouplerConstraint' or c.findtext('isDisabled')!='false' or len(ind)!=1 or dep in self.constraints or functions is None or len(functions)!=1:raise ValueError('Unsupported cervical constraint')
            self.constraints[dep]=(ind[0],*affine(functions[0]))
        self.independent=[n for n in self.names if n not in self.constraints]
        if len(self.independent)!=6 or len(self.constraints)!=18:raise ValueError('Cervical coordinate inventory changed')
        self.n=np.zeros((24,6));self.offset=np.zeros(24);known={name:(np.eye(6)[i],0.) for i,name in enumerate(self.independent)}
        pending=dict(self.constraints)
        while pending:
            done=[]
            for dep,(ind,a,b) in pending.items():
                if ind in known:known[dep]=(a*known[ind][0],a*known[ind][1]+b);done.append(dep)
            if not done:raise ValueError('Unresolved/cyclic cervical couplers')
            for dep in done:del pending[dep]
        for i,name in enumerate(self.names):self.n[i],self.offset[i]=known[name]
        if np.linalg.matrix_rank(self.n)!=6:raise ValueError('Cervical independent N-map is rank deficient')
        self.joints={name:ET.fromstring(body['source_joint_xml']) for name,body in recipe['bodies'].items()}
        self.limits={name:numbers(ET.fromstring(row['source_xml']).findtext('range')) for name,row in recipe['coordinates'].items()}
        self.maps=known

    def coordinates(self,q):
        q=np.asarray(q,float)
        if q.shape!=(6,) or not np.isfinite(q).all() or np.max(np.abs(q))>.05:raise ValueError('Outside cervical engineering domain +/-0.05rad')
        full=self.n@q+self.offset
        for name,value in zip(self.names,full):
            bounds=self.limits[name]
            if bounds.shape!=(2,) or value<bounds[0] or value>bounds[1]:raise ValueError('Outside source cervical coordinate range')
        return q,dict(zip(self.names,full))

    def frames(self,q,speed=None):
        q,full=self.coordinates(q);speed=np.zeros(6) if speed is None else np.asarray(speed,float)
        if speed.shape!=(6,) or not np.isfinite(speed).all():raise ValueError('Six finite cervical speeds required')
        root=constant(self.recipe['donor_spine_to_target_torso']);frames={'spine':root};pending=dict(self.joints)
        while pending:
            done=[]
            for name,joint in pending.items():
                parent=joint.findtext('parent_body')
                if parent not in frames:continue
                # Existing strict evaluator rejects unsupported ordering, reverse,
                # nonzero translations and unresolved functions before jet use.
                reference_joint(joint,full)
                def fixed(location,orientation):
                    t=np.eye(4);t[:3,:3]=Rotation.from_euler('XYZ',numbers(joint.findtext(orientation))).as_matrix();t[:3,3]=numbers(joint.findtext(location));return t
                local=constant(fixed('location_in_parent','orientation_in_parent'))
                for axis in joint.findall('./SpatialTransform/TransformAxis')[:3]:
                    names=(axis.findtext('coordinates') or '').split();a,b=affine(axis.find('function')[0]);row,offset=self.maps[names[0]] if names else (np.zeros(6),0.)
                    row=a*row;angle=float(row@q+a*offset+b);rate=float(row@speed);vector=numbers(axis.findtext('axis'));k=skew(vector)
                    r=Rotation.from_rotvec(vector*angle).as_matrix();t=np.eye(4);t[:3,:3]=r;der=np.zeros((6,4,4));der[:,:3,:3]=row[:,None,None]*(r@k)
                    td=np.zeros((4,4));tdd=np.zeros((4,4));td[:3,:3]=r@k*rate;tdd[:3,:3]=r@k@k*rate**2
                    local=product(local,(t,der,td,tdd))
                local=product(local,constant(np.linalg.inv(fixed('location','orientation'))))
                frames[name]=product(frames[parent],local);done.append(name)
            if not done:raise ValueError('Unresolved cervical body graph')
            for name in done:del pending[name]
        return {name:frames[name] for name in self.joints}

    def metric(self,q,velocity):
        u=np.asarray(velocity,float)
        if u.shape!=(12,) or not np.isfinite(u).all():raise ValueError('Parent6 + cervical6 speeds required')
        v=u[:3];omega=u[3:6];speed=u[6:];matrix=np.zeros((12,12));bias=np.zeros(12);mass=0.;first=np.zeros(3);second=np.zeros((3,3));direct=0.
        for name,(t,d,td,tdd) in self.frames(q,speed).items():
            body=self.recipe['bodies'][name]['mass_properties_in_new_body_frame'];m=body['mass_kg'];c=np.r_[body['center_m'],1.];x=(t@c)[:3];j=np.einsum('iab,b->ia',d,c)[:,:3].T
            h=(tdd@c)[:3];r=t[:3,:3];rd=td[:3,:3];rdd=tdd[:3,:3];w=np.stack([vee(di[:3,:3]@r.T) for di in d],axis=1)
            relative=w@speed;alpha=vee(rdd@r.T+rd@rd.T);inertia=r@np.asarray(body['inertia_kg_m2'])@r.T
            a=np.c_[np.eye(3),-skew(x),j];b=np.c_[np.zeros((3,3)),np.eye(3),w]
            conv=np.cross(omega,v)+np.cross(omega,np.cross(omega,x))+2*np.cross(omega,j@speed)+h
            angular_conv=np.cross(omega,relative)+alpha;absolute=omega+relative
            matrix+=m*a.T@a+b.T@inertia@b;bias+=m*a.T@conv+b.T@(inertia@angular_conv+np.cross(absolute,inertia@absolute))
            direct+=.5*m*np.sum((a@u)**2)+.5*(b@u)@inertia@(b@u)
            mass+=m;first+=m*x;second+=np.trace(inertia)/2*np.eye(3)-inertia+m*np.outer(x,x)
        return {'mass_matrix':matrix,'inertial_bias':bias,'mass_kg':mass,'H':first,'Q':second,'direct_energy_J':float(direct)}


class CoupledThoracicMetric:
    def __init__(self,thorax,composition_plan):
        self.thorax=thorax;self.thorax_dynamics=ThoracicFreeDynamics(thorax)
        path=Path(composition_plan);raw=path.read_bytes();self.identity=hashlib.sha256(raw).hexdigest();plan=json.loads(raw)
        if plan['native_activation_allowed'] is not False:raise ValueError('Source composition cannot authorize native activation')
        root=Path(__file__).resolve().parents[2];sources={}
        for key,receipt in plan['source_receipts'].items():
            content=(root/receipt['path']).read_bytes()
            if hashlib.sha256(content).hexdigest()!=receipt['sha256']:raise ValueError('Composition source identity changed')
            sources[key]=json.loads(content)
        if hashlib.sha256(thorax.path.read_bytes()).hexdigest()!=plan['source_receipts']['thoracic_mechanism']['sha256']:raise ValueError('Wrong thorax mechanism')
        target=plan['target_native_model_identity']['sha256']
        if any(sources[k]['target_model_identity']['sha256']!=target for k in ['cervical_registration','thoracic_anatomy']):raise ValueError('Composite native targets differ')
        original=sources['cervical_inertia']['current_torso'];mass=original['mass_kg'];center=np.asarray(original['center_m']);inertia=np.asarray(original['inertia_kg_m2'])
        expected={'mass_kg':mass,'H_first_moment_kg_m':mass*center,'Q_second_moment_kg_m2':np.trace(inertia)/2*np.eye(3)-inertia+mass*np.outer(center,center)}
        if any(not np.allclose(plan['current_native_torso'][key],value,atol=2e-12,rtol=0) for key,value in expected.items()):raise ValueError('Original torso removal ledger differs from frozen source')
        self.neck=CervicalKinematics(sources['cervical_registration']);self.active=thorax.active+list(range(32,38));self.neck_slots=list(range(6))+list(range(32,38))
        self.plan=plan

    def coordinates(self,q):
        q=np.asarray(q,float)
        if q.shape!=(32,) or not np.isfinite(q).all():raise ValueError('26thorax +6neck finite coordinates required')
        self.thorax.coordinates(q[:26]);self.neck.coordinates(q[26:]);return q

    def evaluate(self,q,velocity,pressure_pa=0.):
        q=self.coordinates(q);u=np.asarray(velocity,float)
        if u.shape!=(38,) or not np.isfinite(u).all():raise ValueError('Parent6 +thorax26 +neck6 speeds required')
        thorax=self.thorax_dynamics.evaluate(q[:26],u[:32]);neck=self.neck.metric(q[26:],u[self.neck_slots])
        matrix=np.zeros((38,38));matrix[:32,:32]=thorax['mass_matrix'];matrix[np.ix_(self.neck_slots,self.neck_slots)]+=neck['mass_matrix']
        bias=np.zeros(38);bias[:32]=thorax['inertial_bias'];bias[self.neck_slots]+=neck['inertial_bias']
        eigen=np.linalg.eigvalsh(matrix[np.ix_(self.active,self.active)])
        if eigen[0]<=0:raise ValueError('Coupled active metric is not positive definite')
        port=self.pressure_port(q,u,pressure_pa);force=port['generalized_force'];acceleration=np.zeros(38);acceleration[self.active]=np.linalg.solve(matrix[np.ix_(self.active,self.active)],(force-bias)[self.active])
        return {'mass_matrix':matrix,'inertial_bias':bias,'velocity_derivative':acceleration,'active_eigenvalues':eigen,
                'mass_kg':thorax['mass_kg']+neck['mass_kg'],'kinetic_energy_J':float(.5*u@matrix@u),
                'direct_material_energy_J':thorax['direct_material_energy_J']+neck['direct_energy_J'],'pressure_port':port,
                'physical_locked_joint_reaction':None}

    def pressure_port(self,q,velocity,pressure_pa):
        q=self.coordinates(q);u=np.asarray(velocity,float)
        if u.shape!=(38,) or not np.isfinite(u).all() or any(abs(u[6+k])>1e-14 for k in self.thorax.locked):raise ValueError('38finite speeds respecting rib locks required')
        cavity=self.thorax.cavity(q[:26],pressure_pa=pressure_pa);force=np.zeros(38);force[6:32]=cavity['pressure_generalized_force'];j=np.zeros(38);j[6:32]=cavity['volume_jacobian_m3_per_coordinate']
        rate=float(j@u)
        return {'geometric_volume_m3':cavity['geometric_volume_m3'],'volume_rate_m3_per_s':rate,'generalized_force':force,
                'volume_speed_jacobian':j,'pressure_pa':float(pressure_pa),'mechanical_power_W':float(pressure_pa*rate),
                'native_mapping_established':False,'convention':'Positive pressure is outward mechanical effort; no native source-flow or gas-volume equivalence assumed'}
