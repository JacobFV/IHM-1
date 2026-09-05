"""SI reduced whole-body mechanics, distinct from native OpenSim integration.

Bones have three rigid translation DOFs with constrained orientation; tissues share these translations/rotations and have
an affine hyperelastic deformation. Attachment links exchange force and torque.
The pressure/volume solve is quasistatic; rigid translation uses linearly implicit Euler.
"""
from __future__ import annotations
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized


def neo_hookean(f, mu, lam):
    """Compressible neo-Hookean density [J/m3], first Piola stress [Pa]."""
    f=np.asarray(f,dtype=float)
    j=float(np.linalg.det(f))
    if not np.isfinite(f).all() or not math.isfinite(j) or j<=0:
        raise ValueError('Finite positive deformation determinant required')
    if mu<=0 or lam<0:raise ValueError('Invalid Lame moduli')
    logj=math.log(j); invt=np.linalg.inv(f).T
    return (float(mu/2*(np.sum(f*f)-3)-mu*logj+lam/2*logj*logj),
            mu*(f-invt)+lam*logj*invt)


def tetra_force_energy(reference, current, mu, lam):
    """Exact constant-strain tetra internal nodal forces and elastic energy."""
    x=np.asarray(reference,dtype=float);y=np.asarray(current,dtype=float)
    dm=(x[1:]-x[0]).T; ds=(y[1:]-y[0]).T
    vol=abs(float(np.linalg.det(dm)))/6
    if vol<=1e-18:raise ValueError('Degenerate reference tetrahedron')
    inv=np.linalg.inv(dm);w,p=neo_hookean(ds@inv,mu,lam)
    h=-vol*p@inv.T
    forces=np.vstack([-h.sum(axis=1),h.T])
    return forces,vol*w


def _rotation(v):
    angle=float(np.linalg.norm(v))
    if angle<1e-15:return np.eye(3)
    x,y,z=v/angle;k=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return np.eye(3)+math.sin(angle)*k+(1-math.cos(angle))*(k@k)


def _value(material, key):return float(material[key]['value'])


class BodyMechanics:
    """A small-motion structural model. All inputs and outputs use SI units.

    from_dict consumes the reproducible mechanics.json artifact. Per-entity
    affine transforms act about reference centroid; no source vertices change.
    Drivers are explicit boundary inputs, never mislabeled live physiology.
    """
    def __init__(self, payload):
        self.payload=payload;self.specs=payload['entities'];self.ids=[e['id'] for e in self.specs]
        self.index={id:i for i,id in enumerate(self.ids)}
        self.x0=np.array([e['centroid_m'] for e in self.specs],float)
        self.x=self.x0.copy();n=len(self.ids)
        self.v=np.zeros((n,3));self.omega=np.zeros((n,3));self.r=np.tile(np.eye(3),(n,1,1))
        self.deformation=np.tile(np.eye(3),(n,1,1))
        self.mass=np.array([e['mass_kg'] for e in self.specs])
        self.inertia=np.array([e['inertia_diagonal_kg_m2'] for e in self.specs])
        if (self.mass<=0).any() or (self.inertia<=0).any():raise ValueError('Positive mass and inertia required')
        self.time=0.;self.work=0.;self.dissipation=0.;self.affine_work=0.;self._soft_energy=0.;self.prescribed_work=0.
        self.links=[]
        stiffness=np.zeros(n)
        for link in payload['links']:
            a=self.index[link['a']];b=self.index[link['b']]
            point=np.asarray(link['point_m'])
            self.links.append((a,b,point-self.x0[a],point-self.x0[b],link['stiffness_n_m'],link['damping_ns_m']))
            stiffness[a]+=link['stiffness_n_m'];stiffness[b]+=link['stiffness_n_m']
        # CFL safety factor is a numerical choice, not a biological parameter.
        self.max_substep=.002  # linearized implicit support solve, refinement tested
        self.muscles=[]
        for m in payload['muscles']:
            anchors=[(self.index[a['entity_id']],np.asarray(a['point_m'])-self.x0[self.index[a['entity_id']]]) for a in m['anchors']]
            self.muscles.append((m,anchors))
        self._build_implicit_operators()
        self._link_a=np.array([z[0] for z in self.links]);self._link_b=np.array([z[1] for z in self.links])
        self._link_oa=np.array([z[2] for z in self.links]);self._link_ob=np.array([z[3] for z in self.links])
        self._link_k=np.array([z[4] for z in self.links]);self._link_c=np.array([z[5] for z in self.links])
        segments=[]
        for j,(m,anchors) in enumerate(self.muscles):
            for (a,oa),(b,ob) in zip(anchors,anchors[1:]):segments.append((a,b,oa,ob,j))
        self._seg_a=np.array([z[0] for z in segments]);self._seg_b=np.array([z[1] for z in segments])
        self._seg_oa=np.array([z[2] for z in segments]);self._seg_ob=np.array([z[3] for z in segments])
        self._seg_m=np.array([z[4] for z in segments])
        self._rest=np.array([m['rest_path_length_m'] for m,_ in self.muscles])
        self._fiber=np.maximum(1e-6,[m['optimal_fiber_length_m'] for m,_ in self.muscles])
        self._fmax=np.array([m['max_isometric_force_n']*math.cos(m['pennation_angle_rad']) for m,_ in self.muscles])
        self._passive=np.array([m['passive_stiffness_n_m'] for m,_ in self.muscles])
        self._soft_previous=set()

    def _build_implicit_operators(self):
        # This reduced skeleton constrains reference orientations. Translational
        # point-attachment stiffness is implicit; holding moments are audited.
        rows=[];cols=[];data=[];ks=[];cs=[];row=0
        pairs=[(a,b,k,c) for a,b,oa,ob,k,c in self.links]
        for m,anchors in self.muscles:
            for (a,oa),(b,ob) in zip(anchors,anchors[1:]):
                # Cauchy-Schwarz bounds the total polyline length Hessian by
                # segment count times the independent-segment bound.
                pairs.append((a,b,m['passive_stiffness_n_m']*(len(anchors)-1),0.))
        for a,b,k,c in pairs:
            rows.extend([row,row]);cols.extend([a,b]);data.extend([-1.,1.]);ks.append(k);cs.append(c);row+=1
        j=sparse.coo_matrix((data,(rows,cols)),shape=(row,len(self.ids))).tocsr()
        self._k=(j.T@sparse.diags(ks)@j).tocsc();self._c=(j.T@sparse.diags(cs)@j).tocsc()
        self._solvers={}

    def _advance(self,h,forces,torques,prescribed=None):
        prescribed={} if prescribed is None else prescribed
        fixed=tuple(sorted(prescribed));key=(round(h,12),fixed)
        if key not in self._solvers:
            matrix=(sparse.diags(self.mass)+h*self._c+h*h*self._k).tocsc()
            free=np.array([i for i in range(len(self.ids)) if i not in prescribed],int)
            solve=factorized(matrix[free][:,free].tocsc()) if len(free) else None
            coupling=matrix[free][:,list(fixed)]
            self._solvers[key]=(matrix,free,coupling,solve)
        # f already includes -Cv. Solve a velocity increment so damping is
        # counted once. Linearized backward Euler: (M+hC+h²K)dv=h(f-hKv).
        matrix,free,coupling,solve=self._solvers[key]
        rhs=h*(forces-h*(self._k@self.v));dv=np.zeros_like(self.v)
        if fixed:
            ix=list(fixed)
            targets=np.array([prescribed[i] for i in fixed])
            dv[ix]=(targets-self.x[ix])/h-self.v[ix]
            if len(free):dv[free]=solve(rhs[free]-coupling@dv[ix])
        else:dv=solve(rhs)
        self.prescribed_reactions=np.zeros_like(self.x)
        if fixed:self.prescribed_reactions[ix]=(matrix@dv-rhs)[ix]/h
        self.v+=dv;self.x+=h*self.v
        self.orientation_reactions=-torques

    @classmethod
    def from_dict(cls,payload):return cls(payload)

    def _point(self,i,offset):
        # Soft anchors follow affine deformation; rigid anchors have F=I.
        arm=self.r[i]@self.deformation[i]@offset
        return self.x[i]+arm,arm

    def _soft_solve(self,drivers,activations):
        energy=0.;pressure_reactions={};material_residual=0.
        volume=drivers.get('volume_ratios',{});pressure=drivers.get('pressure_pa',{});gradients=drivers.get('deformation_gradients',{})
        changed=set(volume)|set(pressure)|set(gradients)|{id for id,a in activations.items() if a}|self._soft_previous
        self._soft_previous=set(volume)|set(pressure)|set(gradients)|{id for id,a in activations.items() if a}
        for id in changed:
            i=self.index[id];e=self.specs[i]
            if e['constitutive']=='rigid':raise ValueError('Soft boundary requires a deformable entity')
            mu=_value(e['material'],'shear_modulus');lam=_value(e['material'],'lame_lambda')
            ratio=float(volume.get(e['id'],1.));p=float(pressure.get(e['id'],0.))
            if not math.isfinite(ratio) or not .25<=ratio<=4 or not math.isfinite(p):
                raise ValueError('Volume ratio must be .25..4 and pressure finite')
            if e['id'] in volume and e['id'] in pressure:raise ValueError('Choose volume or pressure boundary per tissue')
            if e['id'] not in volume and p:
                # Hydrostatic equilibrium Cauchy stress = transmural pressure.
                # This NH law loses pressure stability at large stretch: reject
                # outside the explicitly bounded stable small-strain branch.
                def stress(s):return (mu*(s*s-1)+3*lam*math.log(s))/(s**3)
                lo,hi=.8,1.2
                if not stress(lo)<=p<=stress(hi):raise ValueError('Pressure exceeds reduced model stable strain range')
                for _ in range(50):
                    mid=(lo+hi)/2
                    if stress(mid)<p:lo=mid
                    else:hi=mid
                ratio=((lo+hi)/2)**3
            # Volume-preserving fiber shape equilibrates passive NH stress with
            # an active axial stress a*specific_tension (quasistatic Hill prior).
            activation=activations.get(e['id'],0.)
            axial_stress=activation*_value(e['material'],'active_specific_tension') if 'active_specific_tension' in e['material'] else 0.
            q=0.
            if axial_stress:
                # For stretches [exp(q),exp(-q/2),exp(-q/2)], derivative
                # of W is mu*(exp(2q)-exp(-q)); active stress adds +sigma*q.
                lo,hi=-.35,0.
                for _ in range(45):
                    mid=(lo+hi)/2
                    if mu*ratio**(2/3)*(math.exp(2*mid)-math.exp(-mid))+axial_stress*(1+mid/.35)<0:lo=mid
                    else:hi=mid
                q=(lo+hi)/2
            axis=np.asarray(e['fiber_axis']);axial=np.outer(axis,axis)
            f=ratio**(1/3)*(math.exp(q)*axial+math.exp(-q/2)*(np.eye(3)-axial))
            if e['id'] in gradients:
                f=np.asarray(gradients[e['id']],float);ratio=float(np.linalg.det(f))
            self.deformation[i]=f
            w,piola=neo_hookean(f,mu,lam);energy+=w*e.get('material_volume_m3',e['volume_m3'])
            if e['id'] in volume or e['id'] in pressure or e['id'] in gradients:
                sigma=piola@f.T/ratio
                pressure_reactions[e['id']]=float(np.trace(sigma)/3)
                if p:material_residual=max(material_residual,abs(float(np.trace(sigma)/3)-p))
        return energy,pressure_reactions,material_residual

    def step(self,dt,drivers=None):
        if isinstance(dt,bool):raise ValueError('Boolean dt is not a time interval')
        dt=float(dt);drivers={} if drivers is None else drivers
        if not isinstance(drivers,dict):raise ValueError('Drivers must be an object')
        if not math.isfinite(dt) or not 0<dt<=1:raise ValueError('dt must be in (0,1] seconds')
        allowed={'activation','volume_ratios','pressure_pa','external_forces_n','prescribed_translations_m','deformation_gradients'}
        if set(drivers)-allowed:raise ValueError('Unknown mechanics driver')
        for key,value in drivers.items():
            if not isinstance(value,dict):raise ValueError('Driver maps must be objects')
        for key in ['volume_ratios','pressure_pa','external_forces_n','prescribed_translations_m','deformation_gradients']:
            if set(drivers.get(key,{}))-set(self.ids):raise ValueError('Unknown canonical entity')
        # Validate all boundary maps before any state, energy or time mutation.
        volumes=drivers.get('volume_ratios',{});pressures=drivers.get('pressure_pa',{});gradients=drivers.get('deformation_gradients',{})
        if set(volumes)&set(pressures):raise ValueError('Choose volume or pressure boundary per tissue')
        if set(gradients)&(set(volumes)|set(pressures)):raise ValueError('Choose a single affine, volume or pressure boundary per tissue')
        for id,value in gradients.items():
            f=np.asarray(value,float)
            if self.specs[self.index[id]]['constitutive']=='rigid':raise ValueError('Cannot prescribe rigid body deformation')
            if f.shape!=(3,3) or not np.isfinite(f).all() or not .25<=float(np.linalg.det(f))<=4:raise ValueError('Finite positive deformation with determinant .25..4 required')
        for value in drivers.get('prescribed_translations_m',{}).values():
            vector=np.asarray(value,float)
            if vector.shape!=(3,) or not np.isfinite(vector).all():raise ValueError('Expected finite prescribed translation vector')
        for id,value in volumes.items():
            if isinstance(value,bool) or not math.isfinite(float(value)) or not .25<=float(value)<=4:raise ValueError('Invalid finite volume ratio')
            if self.specs[self.index[id]]['constitutive']=='rigid':raise ValueError('Cannot deform rigid body')
        for id,value in pressures.items():
            if isinstance(value,bool) or not math.isfinite(float(value)):raise ValueError('Invalid pressure')
            e=self.specs[self.index[id]]
            if e['constitutive']=='rigid':raise ValueError('Cannot deform rigid body')
            mu=_value(e['material'],'shear_modulus');lam=_value(e['material'],'lame_lambda')
            stress=lambda s:(mu*(s*s-1)+3*lam*math.log(s))/(s**3)
            if not stress(.8)<=float(value)<=stress(1.2):raise ValueError('Pressure exceeds reduced model stable strain range')
        for id,value in drivers.get('external_forces_n',{}).items():
            vector=np.asarray(value,float)
            if vector.shape!=(3,) or not np.isfinite(vector).all():raise ValueError('Expected finite force vector')
        raw=drivers.get('activation',{});activations={};muscle_activation={}
        known={m['canonical_entity_id'] for m,_ in self.muscles}|{m['id'] for m,_ in self.muscles}|{m.get('source_name','') for m,_ in self.muscles}
        if set(raw)-known:raise ValueError('Unknown muscle activation key')
        for m,_ in self.muscles:
            value=raw.get(m['id'],raw.get(m.get('source_name'),raw.get(m['canonical_entity_id'],0.)))
            if isinstance(value,bool):raise ValueError('Boolean activation is not a scalar activation')
            a=float(value)
            if not math.isfinite(a) or not 0<=a<=1:raise ValueError('Activation must be 0..1')
            muscle_activation[m['id']]=a
            activations[m['canonical_entity_id']]=max(activations.get(m['canonical_entity_id'],0),a)
        soft_energy,pressures,pressure_residual=self._soft_solve(drivers,activations)
        self.affine_work+=soft_energy-self._soft_energy;self._soft_energy=soft_energy
        nstep=max(1,math.ceil(dt/self.max_substep));h=dt/nstep
        prescribed={self.index[id]:self.x0[self.index[id]]+np.asarray(vec,float) for id,vec in drivers.get('prescribed_translations_m',{}).items()}
        prescribed_start={i:self.x[i].copy() for i in prescribed}
        ext=np.zeros_like(self.x)
        for id,vec in drivers.get('external_forces_n',{}).items():
            vec=np.asarray(vec,float)
            if vec.shape!=(3,) or not np.isfinite(vec).all():raise ValueError('Expected finite force vector')
            ext[self.index[id]]=vec
        residual=torque_residual=0.;link_energy=0.;forces_by_muscle={}
        for substep in range(nstep):
            forces=np.zeros_like(self.x);torques=np.zeros_like(self.x);link_energy=0.;active_power=0.;damping_power=0.
            # Vectorized point attachments retain the same force-gradient law.
            def arms(index,offset):
                return np.einsum('nij,nj->ni',self.r[index],np.einsum('nij,nj->ni',self.deformation[index],offset))
            a=self._link_a;b=self._link_b;ra=arms(a,self._link_oa);rb=arms(b,self._link_ob)
            delta=self.x[b]+rb-self.x[a]-ra
            velocity=self.v[b]+np.cross(self.omega[b],rb)-self.v[a]-np.cross(self.omega[a],ra)
            force=self._link_k[:,None]*delta+self._link_c[:,None]*velocity
            np.add.at(forces,a,force);np.add.at(forces,b,-force)
            couple=.5*np.cross(delta,force)
            np.add.at(torques,a,np.cross(ra,force)+couple);np.add.at(torques,b,np.cross(rb,-force)+couple)
            link_energy=.5*float(np.sum(self._link_k*np.sum(delta*delta,axis=1)))
            damping_power=float(np.sum(self._link_c*np.sum(velocity*velocity,axis=1)))
            a=self._seg_a;b=self._seg_b;ra=arms(a,self._seg_oa);rb=arms(b,self._seg_ob)
            segment=self.x[b]+rb-self.x[a]-ra;lengths=np.linalg.norm(segment,axis=1)
            direction=segment/np.maximum(lengths,1e-12)[:,None]
            length=np.bincount(self._seg_m,weights=lengths,minlength=len(self.muscles))
            activation=np.array([muscle_activation[m['id']] for m,_ in self.muscles])
            active=activation*self._fmax*np.exp(-(((length-self._rest)/self._fiber)/.45)**2)
            stretch=np.maximum(0.,length-self._rest);passive=self._passive*stretch;tension=active+passive
            force=tension[self._seg_m,None]*direction
            np.add.at(forces,a,force);np.add.at(forces,b,-force)
            np.add.at(torques,a,np.cross(ra,force));np.add.at(torques,b,np.cross(rb,-force))
            link_energy+=.5*float(np.sum(self._passive*stretch*stretch))
            velocity=self.v[a]+np.cross(self.omega[a],ra)-self.v[b]-np.cross(self.omega[b],rb)
            active_power=float(np.sum(active[self._seg_m,None]*direction*velocity))
            forces_by_muscle={m['id']:float(tension[j]) for j,(m,_) in enumerate(self.muscles)}
            residual=max(residual,float(np.linalg.norm(forces.sum(axis=0))))
            torque_residual=max(torque_residual,float(np.linalg.norm((np.cross(self.x,forces)+torques).sum(axis=0))))
            # Reference orientations are constrained; holding moments are output.
            boundary={i:prescribed_start[i]+(target-prescribed_start[i])*(substep+1)/nstep for i,target in prescribed.items()}
            self._advance(h,forces+ext,torques,boundary)
            velocity=self.v[self._seg_a]-self.v[self._seg_b]
            active_power=float(np.sum(active[self._seg_m,None]*direction*velocity))
            boundary_work=h*float(np.sum(self.prescribed_reactions*self.v))
            self.prescribed_work+=boundary_work
            self.work+=h*(active_power+float(np.sum(ext*self.v)))+boundary_work;self.dissipation+=h*damping_power
        self.time+=dt
        # Evaluate stored energy at returned positions, not the previous substep.
        delta=self.x[self._link_b]+arms(self._link_b,self._link_ob)-self.x[self._link_a]-arms(self._link_a,self._link_oa)
        link_energy=.5*float(np.sum(self._link_k*np.sum(delta*delta,axis=1)))
        segment=self.x[self._seg_b]+arms(self._seg_b,self._seg_ob)-self.x[self._seg_a]-arms(self._seg_a,self._seg_oa)
        length=np.bincount(self._seg_m,weights=np.linalg.norm(segment,axis=1),minlength=len(self.muscles))
        link_energy+=.5*float(np.sum(self._passive*np.maximum(0.,length-self._rest)**2))
        if not np.isfinite(self.x).all():raise FloatingPointError('Mechanics diverged')
        kinetic=float(.5*np.sum(self.mass[:,None]*self.v*self.v)+.5*np.sum(self.inertia*self.omega*self.omega))
        return {'schema_version':1,'model_id':self.payload['model_id'],'time_s':self.time,
                'entities':{id:{'translation_m':(self.x[i]-self.x0[i]).tolist(),'centroid_m':self.x[i].tolist(),'rotation_matrix':self.r[i].tolist(),'deformation_gradient':self.deformation[i].tolist()} for i,id in enumerate(self.ids)},
                'muscle_forces_n':forces_by_muscle,'pressure_reactions_pa':pressures,
                'prescribed_reactions_n':{self.ids[i]:self.prescribed_reactions[i].tolist() for i in prescribed},
                'orientation_reaction_torques_nm':{id:self.orientation_reactions[i].tolist() for i,id in enumerate(self.ids)},
                'audit':{'internal_force_residual_n':residual,'internal_torque_residual_nm':torque_residual,
                         'elastic_energy_j':soft_energy+link_energy,'kinetic_energy_j':kinetic,
                         'accumulated_active_external_work_j':self.work,'accumulated_dissipation_j':self.dissipation,
                         'accumulated_prescribed_boundary_work_j':self.prescribed_work,
                         'dirichlet_position_residual_m':max((float(np.linalg.norm(self.x[i]-target)) for i,target in prescribed.items()),default=0.),
                         'prescribed_affine_boundary_work_j':self.affine_work,'energy_balance_residual_j':soft_energy+link_energy+kinetic+self.dissipation-self.work-self.affine_work,
                         'pressure_equilibrium_residual_pa':pressure_residual,'substeps':nstep,
                         'numerical_scope':'force/torque cancellation; affine boundary work equals quasistatic elastic energy change; dynamic work and dissipation use discrete quadrature; energy residual includes implicit numerical damping',
                         'biological_validation':False},'boundary_drivers':drivers}
