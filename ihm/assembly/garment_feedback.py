"""Partitioned elastic whole-garment feedback into one articulated inertial owner.

The native trajectory and a subcycled cloth solve are iterated from the same
checkpoint. This is a numerical coupling, not a full containment certificate.
"""
from pathlib import Path
import copy,gzip,json,math
import numpy as np
from .garment_body_ports import BarycentricAttachment,attachment_forces,build_source_pairings
from .garment_surface_contact import MovingSurfaceContact
from .whole_garment import materialize_garments


def wrench_point_loads(body,origin_m,force_n,moment_nm,lever_m=.1):
    """Represent a complete ground-frame wrench by four point forces exactly."""
    origin=np.asarray(origin_m,float);force=np.asarray(force_n,float);moment=np.asarray(moment_nm,float)
    if any(a.shape!=(3,) or not np.isfinite(a).all() for a in (origin,force,moment)) or not np.isfinite(lever_m) or lever_m<=0:raise ValueError('Finite wrench and positive lever required')
    axes=np.eye(3);couple=np.cross(np.broadcast_to(moment,(3,3)),axes)/(2*lever_m)
    return [{'body':body,'point_m':origin.tolist(),'force_n':(force-couple.sum(0)).tolist()}]+[
        {'body':body,'point_m':(origin+lever_m*axis).tolist(),'force_n':f.tolist()} for axis,f in zip(axes,couple)]


class GarmentFeedback:
    def __init__(self,garments,skin_positions_m,skin_triangles,skin_body_indices,registration,attachments,*,friction_static,friction_kinetic,identity=None):
        self.garments=garments;self.reference=np.asarray(skin_positions_m,float).copy();self.triangles=np.asarray(skin_triangles,int).copy()
        self.owners=np.asarray(skin_body_indices,int).copy();self.registration=registration;self.attachments=attachments
        if self.owners.shape!=(len(self.reference),) or np.any(self.owners<0) or np.any(self.owners>=len(registration.bodies)):raise ValueError('Invalid material surface owner')
        if set(garments)!=set(attachments):raise ValueError('Missing garment support definition')
        self.friction_static=friction_static;self.friction_kinetic=friction_kinetic;self.identity=identity or {};self.last=None
        self.owner_nodes=[np.flatnonzero(self.owners==i) for i in range(len(registration.bodies))]
        self.hmax=min(b.max_explicit_dt_s for b in garments.values())
        # Account for added support stiffness, not just the cloth edge Laplacian.
        for name,b in garments.items():
            stiffness=np.zeros(len(b.mass_kg))
            for a in attachments[name]:stiffness[a.garment_node]+=a.stiffness_n_m
            active=stiffness>0
            if active.any():self.hmax=min(self.hmax,float(.1*np.sqrt(np.min(b.mass_kg[active]/stiffness[active]))))
    @classmethod
    def from_root(cls,root,registration,*,areal_density_kg_m2=.18,edge_stiffness_n_m=12.,support_stiffness_n_m=12.,friction_static=.4,friction_kinetic=.3):
        root=Path(root);garments,identity=materialize_garments(root,areal_density_kg_m2=areal_density_kg_m2,edge_stiffness_n_m=edge_stiffness_n_m)
        pairings=build_source_pairings(root);skin=json.loads(gzip.decompress((root/identity['skin']['path']).read_bytes()));x=np.asarray(skin['positions']).reshape(-1,3);tri=np.asarray(skin['indices']).reshape(-1,3)
        owners=np.empty(len(x),int)
        # Fixed reference assignment, bounded batch size. Never rebind a material
        # point to a different segment because its current position moved.
        lo=np.array([registration.groups[b]['bounds_min_m'] for b in registration.bodies]);hi=np.array([registration.groups[b]['bounds_max_m'] for b in registration.bodies])
        for start in range(0,len(x),2048):
            points=x[start:start+2048,None,:];distance=np.linalg.norm(np.maximum(np.maximum(lo-points,points-hi),0),axis=2);owners[start:start+len(points)]=np.argmin(distance,axis=1)
        attachments={name:[BarycentricAttachment(p['garment_node'],tuple(p['body_nodes']),tuple(p['barycentric']),p['distance_m'],support_stiffness_n_m) for p in rows['support_pairings']] for name,rows in pairings['garments'].items()}
        identity.update({'pairings':pairings,'support_stiffness_n_m':support_stiffness_n_m,'friction_static':friction_static,'friction_kinetic':friction_kinetic,
                         'binding_basis':'Fixed nearest reference named bone envelope per retained source skin vertex; inferred, no tissue inertia added',
                         'contact_scope':'Retained source face orientation; anatomical exterior/cavity semantics and intersegment skin continuity unresolved',
                         'material_parameters':'Explicit engineering cloth, tether and Coulomb priors; no calibrated whole-fabric or skin-region claim'})
        return cls(garments,x,tri,owners,registration,attachments,friction_static=friction_static,friction_kinetic=friction_kinetic,identity=identity)
    def checkpoint(self):
        return {'garments':{n:{'position_m':b.position_m.copy(),'velocity_m_s':b.velocity_m_s.copy(),'time_s':b.time_s} for n,b in self.garments.items()},'last':copy.deepcopy(self.last)}
    def restore(self,state):
        for n,v in state['garments'].items():
            b=self.garments[n];b.position_m=v['position_m'].copy();b.velocity_m_s=v['velocity_m_s'].copy();b.time_s=v['time_s']
        self.last=copy.deepcopy(state['last'])
    def surface(self,native):
        result=np.empty_like(self.reference);transforms=self.registration.transforms(native)
        for body,ids in zip(self.registration.bodies,self.owner_nodes):
            t=transforms[body];result[ids]=self.reference[ids]@t[:3,:3].T+t[:3,3]
        return result
    def _integrate(self,start,end,dt):
        a=self.surface(start);b=self.surface(end);count=math.ceil(dt/self.hmax);h=dt/count
        surface=MovingSurfaceContact(a,b,self.triangles,dt,friction_static=self.friction_static,friction_kinetic=self.friction_kinetic)
        impulse=np.zeros((len(self.owner_nodes),3));angular=np.zeros_like(impulse);gravity=self.registration.basis@start['gravity_m_s2']
        totals={'normal_impact_dissipation_j':0.,'friction_dissipation_j':0.,'cloth_numerical_energy_defect_j':0.,'attachment_numerical_energy_defect_j':0.,'interface_body_work_j':0.,'contact_count':0,'unresolved_edge_contacts':0,'maximum_position_correction_m':0.}
        for k in range(count):
            old_surface=a+(k/count)*(b-a);surface.fraction=(k+1)/count;current_surface=a+surface.fraction*(b-a)
            for name,cloth in self.garments.items():
                old_cloth=cloth.position_m.copy();spring=attachment_forces(old_cloth,old_surface,self.attachments[name]);result=cloth.step(h,gravity_m_s2=gravity,contact=surface,external_forces_n=spring['garment_force_n']);contact=result['contact']
                next_spring=attachment_forces(cloth.position_m,current_surface,self.attachments[name])
                spring_work=float(np.sum(spring['garment_force_n']*(cloth.position_m-old_cloth))+np.sum(spring['body_force_n']*(current_surface-old_surface)))
                totals['attachment_numerical_energy_defect_j']+=next_spring['energy_j']-spring['energy_j']+spring_work
                spring_impulse=h*spring['body_force_n'];contact_impulse=contact['body_reaction_impulses_ns'];j=spring_impulse+contact_impulse
                for i,ids in enumerate(self.owner_nodes):
                    impulse[i]+=j[ids].sum(0);angular[i]+=np.cross(old_surface[ids],spring_impulse[ids]).sum(0)+np.cross(current_surface[ids],contact_impulse[ids]).sum(0)
                totals['interface_body_work_j']+=float(np.sum(j*surface.velocity))
                totals['normal_impact_dissipation_j']+=result['normal_impact_dissipation_j'];totals['friction_dissipation_j']+=result['friction_dissipation_j']
                totals['cloth_numerical_energy_defect_j']+=result['numerical_energy_defect_j'];totals['contact_count']+=contact['contact_count'];totals['unresolved_edge_contacts']+=contact['unresolved_edge_contacts']
                totals['maximum_position_correction_m']=max(totals['maximum_position_correction_m'],result['max_position_correction_m'])
        loads=[];canonical=[];r=self.registration.basis;c=self.registration.global_map[:3,3]
        for i,body in enumerate(self.registration.bodies):
            origin=np.array(start['bodies'][body]['transform_ground'])[:3,3];canonical_origin=r@origin+c;force=impulse[i]/dt;moment=angular[i]/dt-np.cross(canonical_origin,force)
            loads.extend(wrench_point_loads(body,origin,r.T@force,r.T@moment))
            canonical.append({'body':body,'point_m':canonical_origin.tolist(),'force_n':force.tolist(),'moment_nm':moment.tolist(),'kind':'interval_mean_resultant_wrench'})
        totals.update({'substeps':count,'dt_s':dt,'canonical_reaction_wrenches':canonical,'body_impulse_ns':impulse.sum(0).tolist()})
        return loads,totals
    def advance(self,native,dt,external_loads,actuation,*,tolerance_m=2e-7,max_iterations=5):
        if not np.isfinite(dt) or not 0<dt<=.02 or not np.isfinite(tolerance_m) or tolerance_m<=0 or isinstance(max_iterations,bool) or not isinstance(max_iterations,int) or max_iterations<1:raise ValueError('Invalid coupling interval or iteration controls')
        start=native.snapshot();saved=self.checkpoint()
        if any(abs(b.time_s-start['time_s'])>1e-9 for b in self.garments.values()):raise ValueError('Native and cloth clocks differ')
        token=native.checkpoint();loads=[]
        try:
            prediction=native.advance(dt,external_loads,actuation)
            for iteration in range(max_iterations):
                self.restore(saved);loads,report=self._integrate(start,prediction,dt)
                native.restore(token);result=native.advance(dt,[*external_loads,*loads],actuation)
                residual=float(np.linalg.norm(self.surface(result)-self.surface(prediction),axis=1).max())
                if residual<=tolerance_m:
                    native_work=0.
                    for load in loads:
                        old=start['bodies'][load['body']];new=result['bodies'][load['body']];t0=np.array(old['transform_ground']);t1=np.array(new['transform_ground']);p0=np.array(load['point_m']);p1=(t1@np.linalg.inv(t0)@np.r_[p0,1])[:3]
                        v0=np.array(old['origin_velocity_m_s'])+np.cross(old['angular_velocity_rad_s'],p0-t0[:3,3]);v1=np.array(new['origin_velocity_m_s'])+np.cross(new['angular_velocity_rad_s'],p1-t1[:3,3])
                        native_work+=dt*.5*np.dot(load['force_n'],v0+v1)
                    report.update({'native_garment_work_j':float(native_work),'interface_work_quadrature_defect_j':float(native_work-report['interface_body_work_j']),'iterations':iteration+1,'boundary_position_residual_m':residual,'boundary_tolerance_m':tolerance_m,
                                   'coverage':'Partitioned endpoint trajectory, interval-mean work-conjugate body loads, explicit cloth springs and discrete face contact; CCD/self/edge and full containment remain unresolved'})
                    self.last=report;return result,copy.deepcopy(report)
                prediction=result
            raise RuntimeError('Garment/native interface failed convergence; interval rolled back')
        except BaseException:
            native.restore(token);self.restore(saved);raise
        finally:native.release(token)
    def frame(self):
        return {'garments':{n:{'positions_m':b.position_m.tolist(),'time_s':b.time_s} for n,b in self.garments.items()},'audit':copy.deepcopy(self.last),'identity':self.identity}
