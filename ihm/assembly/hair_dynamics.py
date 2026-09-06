"""Opt-in authoritative SI beam state for mechanical coupling experiments.

Initially straight, uniform, small-deflection rods; no JS state synchronization.
The first two nodes form a moving material clamp. Contact uses free centerline
nodes (zero clearance), not finite-radius shaft/edge/continuous contact.
"""
import copy,math
import numpy as np
from scipy.linalg import cho_factor,cho_solve


class ElasticHairState:
    def __init__(self,data,*,max_substep_s=1/240):
        x=np.asarray(data['centerlines_m'],float).reshape(-1,3)
        offsets=np.asarray(data['strand_offsets']);r=np.asarray(data['radius_m'],float)
        if x.ndim!=2 or not np.isfinite(x).all() or r.ndim!=1 or not 0<len(r)<=512 or not np.isfinite(r).all() or np.any(r<=0):raise ValueError('Bounded finite hair geometry required')
        if offsets.dtype.kind not in 'iu' or offsets.shape!=(len(r)+1,) or offsets[0]!=0 or offsets[-1]!=len(x) or np.any(np.diff(offsets)<3) or np.any(np.diff(offsets)>32):raise ValueError('Each guide needs 3..32 nodes and valid offsets')
        E,B,rho=[data[k] for k in ('tensile_modulus_pa','bending_modulus_pa','density_kg_m3')]
        if not np.isfinite([E,B,rho,max_substep_s]).all() or min(E,B,rho,max_substep_s)<=0 or max_substep_s>.02:raise ValueError('Positive finite material and bounded substep required')
        self.position_m=x.copy();self.velocity_m_s=np.zeros_like(x);self.mass_kg=np.zeros(len(x));self.offsets=offsets.copy();self.radius_m=r.copy();self.max_substep_s=max_substep_s;self.time_s=0.
        self.rods=[];self.root_m=[];self.tangent=[];self.fixed=[]
        for s,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
            n=b-a;delta=np.diff(x[a:b],axis=0);h=float(np.linalg.norm(delta[0]));d=delta[0]/h if h>0 else delta[0]
            if h<=0 or not np.allclose(delta,delta[0],rtol=0,atol=h*1e-6):raise ValueError('Initially straight uniform rods required')
            area=np.pi*r[s]**2;moment=np.pi*r[s]**4/4
            self.mass_kg[a:b]=rho*area*h;self.mass_kg[[a,b-1]]*=.5
            D1=np.diff(np.eye(n),axis=0);D2=np.diff(np.eye(n),n=2,axis=0)
            self.rods.append({'a':a,'b':b,'h':h,'axial':E*area/h*(D1.T@D1),'bend':B*moment/h**3*(D2.T@D2),'cache':{}})
            self.root_m.append(x[a].copy());self.tangent.append(d);self.fixed.extend([a,a+1])
        self.root_m=np.array(self.root_m);self.tangent=np.array(self.tangent);self.fixed=np.array(self.fixed)
        self.free=np.setdiff1d(np.arange(len(x)),self.fixed)

    def checkpoint(self):
        return {k:copy.deepcopy(getattr(self,k)) for k in ('position_m','velocity_m_s','root_m','tangent','time_s')}

    def restore(self,state):
        prepared={}
        for k in ('position_m','velocity_m_s','root_m','tangent'):
            value=np.asarray(state[k],float)
            if value.shape!=getattr(self,k).shape or not np.isfinite(value).all():raise ValueError('Invalid hair checkpoint')
            prepared[k]=value.copy()
        if not np.isfinite(state['time_s']) or state['time_s']<0:raise ValueError('Invalid hair checkpoint clock')
        for k,v in prepared.items():setattr(self,k,v)
        self.time_s=float(state['time_s'])

    def energy_j(self):
        energy=.5*float(np.sum(self.mass_kg[:,None]*self.velocity_m_s**2))
        for s,rod in enumerate(self.rods):
            a,b=rod['a'],rod['b'];d=self.tangent[s]
            y=self.position_m[a:b]-self.root_m[s]-np.arange(b-a)[:,None]*rod['h']*d
            along=y@d;transverse=y-along[:,None]*d
            energy+=.5*float(along@rod['axial']@along+np.sum(transverse*(rod['bend']@transverse)))
        return energy

    def step(self,dt_s,*,roots_m,root_tangents,gravity_m_s2,contact=None):
        roots=np.asarray(roots_m,float);directions=np.asarray(root_tangents,float);gravity=np.asarray(gravity_m_s2,float)
        if not np.isfinite(dt_s) or not 0<dt_s<=.02 or roots.shape!=self.root_m.shape or directions.shape!=self.tangent.shape or gravity.shape!=(3,) or not all(np.isfinite(v).all() for v in (roots,directions,gravity)):raise ValueError('Finite common-clock hair interval and boundary required')
        lengths=np.linalg.norm(directions,axis=1)
        if np.any(lengths<=0):raise ValueError('Nonzero root tangents required')
        directions=directions/lengths[:,None];start_roots=self.root_m.copy();start_directions=self.tangent.copy()
        if np.any(np.sum(start_directions*directions,axis=1)<=0):raise ValueError('Root rotation over 90 degrees needs smaller coupled intervals')
        saved=self.checkpoint();initial_energy=self.energy_j();initial_momentum=(self.mass_kg[:,None]*self.velocity_m_s).sum(0)
        initial_angular=np.cross(self.position_m,self.mass_kg[:,None]*self.velocity_m_s).sum(0);gravity_angular=np.zeros(3)
        steps=max(1,math.ceil(dt_s/self.max_substep_s-1e-10));h=dt_s/steps
        if steps>240:raise ValueError('Hair substep budget exceeded; reduce the parent interval')
        clamp=np.zeros((len(self.rods),3));clamp_angular=np.zeros_like(clamp)
        reaction=np.zeros_like(contact.start) if contact is not None else np.zeros((0,3))
        vertex_angular_reaction=np.zeros_like(reaction)
        angular_reaction=np.zeros(3);external_work=0.;clamp_work=0.;contact_work=0.;normal_loss=0.;friction_loss=0.;count=0;unresolved=0;maximum_correction=0.
        try:
            for step in range(1,steps+1):
                old_x=self.position_m.copy();old_v=self.velocity_m_s.copy();fraction=step/steps
                self.root_m=start_roots*(1-fraction)+roots*fraction
                self.tangent=start_directions*(1-fraction)+directions*fraction;self.tangent/=np.linalg.norm(self.tangent,axis=1)[:,None]
                for s,rod in enumerate(self.rods):
                    a,b=rod['a'],rod['b'];d=self.tangent[s];mass=self.mass_kg[a:b];n=b-a
                    rest=self.root_m[s]+np.arange(n)[:,None]*rod['h']*d
                    predicted=old_x[a:b]+h*old_v[a:b]+h*h*gravity-rest
                    axial=predicted@d;transverse=predicted-axial[:,None]*d
                    factors=rod['cache'].get(h)
                    if factors is None:
                        factors=[cho_factor(np.diag(mass[2:]/h**2)+rod[key][2:,2:],check_finite=False) for key in ('axial','bend')]
                        rod['cache'][h]=factors
                        if len(rod['cache'])>8:del rod['cache'][next(iter(rod['cache']))]
                    along=cho_solve(factors[0],mass[2:]/h**2*axial[2:],check_finite=False)
                    across=cho_solve(factors[1],mass[2:,None]/h**2*transverse[2:],check_finite=False)
                    self.position_m[a:b]=rest;self.position_m[a+2:b]+=along[:,None]*d+across
                    self.velocity_m_s[a:b]=(self.position_m[a:b]-old_x[a:b])/h
                    y=self.position_m[a:b]-rest;parallel=y@d
                    elastic=-(rod['axial']@parallel)[:,None]*d-rod['bend']@(y-parallel[:,None]*d)
                    j=mass[:2,None]*(self.velocity_m_s[a:a+2]-old_v[a:a+2])-h*(mass[:2,None]*gravity+elastic[:2])
                    clamp[s]+=j.sum(0);clamp_angular[s]+=np.cross(self.position_m[a:a+2],j).sum(0)
                    clamp_work+=float(np.sum(j*(self.position_m[a:a+2]-old_x[a:a+2])/h))
                external_work+=float(np.sum(self.mass_kg[:,None]*gravity*(self.position_m-old_x)))
                gravity_angular+=h*np.cross(self.position_m,self.mass_kg[:,None]*gravity).sum(0)
                if contact is not None:
                    contact.fraction=fraction
                    receipt=contact.resolve(self.position_m[self.free],self.velocity_m_s[self.free],self.mass_kg[self.free],h)
                    self.position_m[self.free]=receipt['positions_m'];self.velocity_m_s[self.free]=receipt['velocities_m_s']
                    reaction+=receipt['body_reaction_impulses_ns']
                    current_angular=np.cross(contact.start+fraction*contact.delta,receipt['body_reaction_impulses_ns'])
                    vertex_angular_reaction+=current_angular;angular_reaction+=current_angular.sum(0)
                    contact_work+=receipt['prescribed_surface_work_j'];normal_loss+=receipt['normal_impact_dissipation_j'];friction_loss+=receipt['friction_dissipation_j'];count+=receipt['contact_count'];unresolved+=receipt['unresolved_edge_contacts']
                    maximum_correction=max(maximum_correction,float(np.linalg.norm(receipt['position_correction_m'],axis=1).max(initial=0)))
                if not np.isfinite(self.position_m).all() or not np.isfinite(self.velocity_m_s).all():raise RuntimeError('Hair solve produced nonfinite state')
            self.time_s+=dt_s
            delta_energy=self.energy_j()-initial_energy
            gravity_impulse=dt_s*self.mass_kg.sum()*gravity
            delta_momentum=(self.mass_kg[:,None]*self.velocity_m_s).sum(0)-initial_momentum
            body_impulse=reaction.sum(0)-clamp.sum(0)
            angular_residual=np.cross(self.position_m,self.mass_kg[:,None]*self.velocity_m_s).sum(0)-initial_angular-gravity_angular+angular_reaction-clamp_angular.sum(0)
            slope=0.
            for s,rod in enumerate(self.rods):
                edges=np.diff(self.position_m[rod['a']:rod['b']],axis=0);d=self.tangent[s]
                slope=max(slope,float(np.linalg.norm(edges-(edges@d)[:,None]*d,axis=1).max()/rod['h']))
            return {'body_reaction_impulses_ns':reaction,'body_reaction_angular_impulses_nms':vertex_angular_reaction,'follicle_body_impulses_ns':-clamp,'follicle_body_angular_impulses_nms':-clamp_angular,
                    'surface_contact_impulse_ns':reaction.sum(0),'surface_contact_angular_impulse_nms':angular_reaction,
                    'hair_body_impulse_residual_ns':delta_momentum-gravity_impulse+body_impulse,
                    'hair_body_angular_impulse_residual_nms':angular_residual,'maximum_transverse_slope':slope,'within_small_deflection':slope<=.3,
                    'hair_energy_change_j':delta_energy,'interface_body_work_j':-clamp_work-contact_work,
                    'external_gravity_work_j':external_work,'normal_impact_dissipation_j':normal_loss,'friction_dissipation_j':friction_loss,
                    'hair_numerical_energy_defect_j':delta_energy-external_work-clamp_work-contact_work+normal_loss+friction_loss,
                    'contact_count':count,'unresolved_edge_contacts':unresolved,'maximum_position_correction_m':maximum_correction,'substeps':steps}
        except BaseException:
            self.restore(saved);raise
