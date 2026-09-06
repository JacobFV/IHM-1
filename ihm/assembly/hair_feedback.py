"""Explicit opt-in hair/body port; never attached implicitly to the live app.

The caller supplies source vertices, topology and fixed registered owners.
Only single-owner rigid follicle triangles are supported. Native loads use the
same full-wrench contract as GarmentFeedback. No new body inertia is created.
Hair inertia is additional to the supplied native model: a source-bound factory
must reconcile existing anatomical hair proxy mass before live enablement.
"""
import copy
import numpy as np
from .garment_surface_contact import MovingSurfaceContact
from .garment_feedback import wrench_point_loads


class HairFeedback:
    root_tolerance_m=1e-9
    def __init__(self,hair,skin_positions_m,skin_triangles,skin_body_indices,registration,follicle_triangles,barycentric,*,identity,friction_static,friction_kinetic):
        self.hair=hair;self.reference=np.asarray(skin_positions_m,float).copy();tri=np.asarray(skin_triangles);owners=np.asarray(skin_body_indices)
        follicles=np.asarray(follicle_triangles);beta=np.asarray(barycentric,float)
        if self.reference.ndim!=2 or self.reference.shape[1]!=3 or not np.isfinite(self.reference).all():raise ValueError('Finite registered source vertices required')
        if tri.ndim!=2 or tri.shape[1]!=3 or not len(tri) or tri.dtype.kind not in 'iu' or tri.min()<0 or tri.max()>=len(self.reference):raise ValueError('Valid registered source faces required')
        if len(self.reference)>65536 or len(tri)>131072:raise ValueError('Contact patch exceeds explicit source geometry budget')
        if owners.dtype.kind not in 'iu' or owners.shape!=(len(self.reference),) or owners.min()<0 or owners.max()>=len(registration.bodies):raise ValueError('Fixed registered owner for every source vertex required')
        if follicles.dtype.kind not in 'iu' or follicles.shape!=(len(hair.rods),) or follicles.min()<0 or follicles.max()>=len(tri):raise ValueError('Exact source follicle face indices required')
        if beta.shape!=(len(hair.rods),3) or not np.isfinite(beta).all() or np.any(beta<0) or not np.allclose(beta.sum(1),1,rtol=0,atol=1e-12):raise ValueError('Convex follicle material weights required')
        if not isinstance(identity,dict) or not identity.get('source_id'):raise ValueError('Explicit source identity required')
        follicle_nodes=tri[follicles];root_owners=owners[follicle_nodes]
        if np.any(root_owners!=root_owners[:,0,None]):raise ValueError('Cross-owner follicle triangles have no rigid clamp derivative')
        roots=np.sum(beta[:,:,None]*self.reference[follicle_nodes],axis=1)
        if not np.allclose(roots,hair.root_m,rtol=0,atol=1e-10):raise ValueError('Hair roots do not reconstruct registered source material points')
        self.triangles=tri.copy();self.owners=owners.copy();self.registration=registration;self.follicle_nodes=follicle_nodes;self.barycentric=beta.copy();self.root_owners=root_owners[:,0];self.reference_tangents=hair.tangent.copy()
        self.owner_nodes=[np.flatnonzero(self.owners==i) for i in range(len(registration.bodies))]
        self.friction_static=friction_static;self.friction_kinetic=friction_kinetic;self.identity=copy.deepcopy(identity);self.last=None
        # Reuse contact validation without constructing any native body.
        MovingSurfaceContact(self.reference,self.reference,self.triangles,.001,friction_static=friction_static,friction_kinetic=friction_kinetic)

    def checkpoint(self):return {'hair':self.hair.checkpoint(),'last':copy.deepcopy(self.last)}
    def restore(self,state):self.hair.restore(state['hair']);self.last=copy.deepcopy(state['last'])

    def surface(self,native):
        result=np.empty_like(self.reference);transforms=self.registration.transforms(native)
        for body,ids in zip(self.registration.bodies,self.owner_nodes):
            transform=transforms[body];result[ids]=self.reference[ids]@transform[:3,:3].T+transform[:3,3]
        return result

    def roots(self,native):
        surface=self.surface(native);transforms=self.registration.transforms(native)
        roots=np.sum(self.barycentric[:,:,None]*surface[self.follicle_nodes],axis=1)
        directions=np.array([transforms[self.registration.bodies[owner]][:3,:3]@direction for owner,direction in zip(self.root_owners,self.reference_tangents)])
        return roots,directions

    def clamp_residual(self,native):
        roots,directions=self.roots(native);spacing=np.array([r['h'] for r in self.hair.rods])
        expected=np.stack((roots,roots+spacing[:,None]*directions),axis=1).reshape(-1,3)
        return float(np.linalg.norm(self.hair.position_m[self.hair.fixed]-expected,axis=1).max())

    def _integrate(self,start,end,dt):
        a=self.surface(start);b=self.surface(end);roots,directions=self.roots(end)
        initial_roots,_=self.roots(start)
        if float(np.linalg.norm(self.hair.root_m-initial_roots,axis=1).max())>self.root_tolerance_m:raise ValueError('Hair state is not synchronized to registered body roots')
        if self.clamp_residual(start)>self.root_tolerance_m:raise ValueError('Hair clamp direction is not synchronized to registered body')
        contact=MovingSurfaceContact(a,b,self.triangles,dt,friction_static=self.friction_static,friction_kinetic=self.friction_kinetic)
        report=self.hair.step(dt,roots_m=roots,root_tangents=directions,gravity_m_s2=self.registration.basis@start['gravity_m_s2'],contact=contact)
        impulses=report.pop('body_reaction_impulses_ns');angular=report.pop('body_reaction_angular_impulses_nms');follicle=report.pop('follicle_body_impulses_ns');follicle_angular=report.pop('follicle_body_angular_impulses_nms')
        owner_impulse=np.zeros((len(self.registration.bodies),3));owner_angular=np.zeros_like(owner_impulse)
        for i,ids in enumerate(self.owner_nodes):owner_impulse[i]+=impulses[ids].sum(0);owner_angular[i]+=angular[ids].sum(0)
        for s,owner in enumerate(self.root_owners):owner_impulse[owner]+=follicle[s];owner_angular[owner]+=follicle_angular[s]
        loads=[];canonical=[];basis=self.registration.basis;offset=self.registration.global_map[:3,3]
        for i,body in enumerate(self.registration.bodies):
            origin=np.asarray(start['bodies'][body]['transform_ground'])[:3,3];canonical_origin=basis@origin+offset
            force=owner_impulse[i]/dt;moment=owner_angular[i]/dt-np.cross(canonical_origin,force)
            loads.extend(wrench_point_loads(body,origin,basis.T@force,basis.T@moment))
            canonical.append({'body':body,'point_m':canonical_origin.tolist(),'force_n':force.tolist(),'moment_nm':moment.tolist(),'kind':'interval_mean_hair_reaction_wrench'})
        report['canonical_reaction_wrenches']=canonical;report['body_impulse_ns']=owner_impulse.sum(0);report['body_angular_impulse_nms']=owner_angular.sum(0)
        return loads,report

    def advance(self,native,dt,external_loads,actuation,*,tolerance_m=2e-7,max_iterations=5):
        if not np.isfinite(dt) or not 0<dt<=.02 or not np.isfinite(tolerance_m) or tolerance_m<=0 or isinstance(max_iterations,bool) or not isinstance(max_iterations,int) or max_iterations<1:raise ValueError('Invalid hair/native coupling interval')
        start=native.snapshot();saved=self.checkpoint()
        if abs(self.hair.time_s-start['time_s'])>1e-9:raise ValueError('Native and authoritative hair clocks differ')
        token=native.checkpoint()
        try:
            prediction=native.advance(dt,external_loads,actuation)
            for iteration in range(max_iterations):
                self.restore(saved);loads,report=self._integrate(start,prediction,dt)
                native.restore(token);result=native.advance(dt,[*external_loads,*loads],actuation)
                residual=float(np.linalg.norm(self.surface(result)-self.surface(prediction),axis=1).max())
                result_roots,_=self.roots(result);root_residual=float(np.linalg.norm(self.hair.root_m-result_roots,axis=1).max())
                clamp_residual=self.clamp_residual(result)
                if residual<=tolerance_m and max(root_residual,clamp_residual)<=self.root_tolerance_m:
                    native_work=0.
                    for load in loads:
                        old=start['bodies'][load['body']];new=result['bodies'][load['body']];t0=np.asarray(old['transform_ground']);t1=np.asarray(new['transform_ground']);p0=np.asarray(load['point_m']);p1=(t1@np.linalg.inv(t0)@np.r_[p0,1])[:3]
                        v0=np.asarray(old['origin_velocity_m_s'])+np.cross(old['angular_velocity_rad_s'],p0-t0[:3,3]);v1=np.asarray(new['origin_velocity_m_s'])+np.cross(new['angular_velocity_rad_s'],p1-t1[:3,3])
                        native_work+=dt*.5*float(np.dot(load['force_n'],v0+v1))
                    report.update({'native_hair_work_j':native_work,'hair_plus_native_interface_work_j':report['hair_energy_change_j']+native_work,
                                   'interface_work_quadrature_defect_j':native_work-report['interface_body_work_j'],'iterations':iteration+1,'boundary_position_residual_m':residual,'boundary_tolerance_m':tolerance_m,
                                   'root_position_residual_m':root_residual,'root_tolerance_m':self.root_tolerance_m,
                                   'clamp_position_residual_m':clamp_residual,
                                   'coverage':'Opt-in authoritative Python straight small-deflection beams; moving registered face-interior centerline contact and full single-owner follicle clamp reaction; zero radial clearance, no edges/CCD/self-contact, no viewer synchronization'})
                    self.last=copy.deepcopy(report);return result,copy.deepcopy(report)
                prediction=result
            raise RuntimeError('Hair/native interface failed convergence; interval rolled back')
        except BaseException:
            native.restore(token);self.restore(saved);raise
        finally:native.release(token)

    def frame(self):
        def serializable(value):
            if isinstance(value,np.ndarray):return value.tolist()
            if isinstance(value,np.generic):return value.item()
            if isinstance(value,dict):return {k:serializable(v) for k,v in value.items()}
            if isinstance(value,(list,tuple)):return [serializable(v) for v in value]
            return copy.deepcopy(value)
        return serializable({'hair':{'positions_m':self.hair.position_m,'velocities_m_s':self.hair.velocity_m_s,'time_s':self.hair.time_s},'audit':self.last,'identity':self.identity,'viewer_synchronized':False})
