"""Finite-mass elastic cloth and unilateral point contact, in SI units.

This independent solver is a mechanical building block, not yet the body's
two-way garment/genital contact model. Edge springs have no bending law and
node contact does not provide continuous triangle or self-collision detection.
All constructor coefficients are explicit engineering inputs, not calibration.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def _positive(value, name, *, zero=False):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value):
        raise ValueError(f'{name} must be a finite scalar')
    value=float(value)
    if not np.isfinite(value) or value<0 or (value==0 and not zero):
        raise ValueError(f'{name} must be finite and {"nonnegative" if zero else "positive"}')
    return value


def _vector(value,name):
    a=np.asarray(value,dtype=float)
    if a.shape!=(3,) or not np.isfinite(a).all():
        raise ValueError(f'{name} must have three finite components')
    return a


@dataclass
class PlaneContact:
    """Prescribed planar body surface: n·x >= offset + thickness.

    Resolve uses zero restitution and Coulomb impulses. Surface motion is an
    explicit velocity input; the caller must advance its geometric offset.
    Returned body impulses permit a future two-way owner to apply reactions.
    """
    normal: object
    offset_m: float
    friction_static: float
    friction_kinetic: float
    thickness_m: float=0.
    surface_velocity_m_s: object=(0.,0.,0.)

    def __post_init__(self):
        self.normal=_vector(self.normal,'normal')
        norm=np.linalg.norm(self.normal)
        if not np.isclose(norm,1.,rtol=0,atol=1e-12):
            raise ValueError('normal must have unit length')
        if not np.isfinite(self.offset_m):raise ValueError('offset_m must be finite')
        self.thickness_m=_positive(self.thickness_m,'thickness_m',zero=True)
        self.friction_static=_positive(self.friction_static,'friction_static',zero=True)
        self.friction_kinetic=_positive(self.friction_kinetic,'friction_kinetic',zero=True)
        if self.friction_kinetic>self.friction_static:
            raise ValueError('kinetic friction must not exceed static friction')
        self.surface_velocity_m_s=_vector(self.surface_velocity_m_s,'surface velocity')

    def resolve(self,positions_m,velocities_m_s,mass_kg,dt_s):
        dt_s=_positive(dt_s,'dt_s')
        x=np.array(positions_m,dtype=float,copy=True)
        v=np.array(velocities_m_s,dtype=float,copy=True)
        mass=np.asarray(mass_kg,dtype=float)
        if x.ndim!=2 or x.shape[1]!=3 or v.shape!=x.shape or mass.shape!=(len(x),):
            raise ValueError('contact arrays must be Nx3, Nx3 and N')
        if not (np.isfinite(x).all() and np.isfinite(v).all() and np.isfinite(mass).all() and (mass>0).all()):
            raise ValueError('contact inputs must be finite, with positive masses')
        gap=x@self.normal-self.offset_m-self.thickness_m
        active=gap<=1e-12
        correction=np.maximum(-gap,0)[:,None]*self.normal
        x+=correction
        relative=v-self.surface_velocity_m_s
        vn=relative@self.normal
        jn=np.where(active,mass*np.maximum(-vn,0),0.)
        tangential=relative-vn[:,None]*self.normal
        speed=np.linalg.norm(tangential,axis=1)
        required=mass*speed
        stick=active & (required<=self.friction_static*jn)
        magnitude=np.where(stick,required,np.minimum(required,self.friction_kinetic*jn))
        direction=np.divide(tangential,speed[:,None],out=np.zeros_like(tangential),where=speed[:,None]>0)
        jt=-magnitude[:,None]*direction
        impulse=jn[:,None]*self.normal+jt
        friction_loss=-np.sum(jt*(tangential+.5*jt/mass[:,None]))
        normal_loss=np.sum(.5*jn*jn/mass)
        v+=impulse/mass[:,None]
        return {
            'positions_m':x,'velocities_m_s':v,'impulses_ns':impulse,
            'body_reaction_impulses_ns':-impulse,
            'body_reaction_forces_n':-impulse/dt_s,
            'friction_dissipation_j':float(friction_loss),
            'normal_impact_dissipation_j':float(normal_loss),
            'prescribed_surface_work_j':float(np.sum(impulse*self.surface_velocity_m_s)),
            'position_correction_m':correction,
            'max_penetration_before_projection_m':float(np.max(np.maximum(-gap,0),initial=0)),
            'contact_count':int(np.count_nonzero(active)),
            'sticking_count':int(np.count_nonzero(stick & (jn>0))),
        }


class Cloth:
    """Triangulated, lumped-mass cloth with objective elastic edge springs.

    Energy = sum_e k_e (length_e - rest_e)^2 / 2. Symplectic Euler
    advances velocities then positions; dt must resolve the elastic frequencies.
    The caller supplies physical rest lengths (including elastic-band prestrain)
    by choosing the reference mesh, independently of initial deformed positions.
    """
    def __init__(self,reference_positions_m,triangles,*,areal_density_kg_m2,edge_stiffness_n_m):
        x=np.asarray(reference_positions_m,dtype=float)
        raw=np.asarray(triangles)
        if x.ndim!=2 or x.shape[1]!=3 or len(x)<3 or not np.isfinite(x).all():
            raise ValueError('reference positions must be finite Nx3')
        if raw.ndim!=2 or raw.shape[1]!=3 or not len(raw) or raw.dtype.kind not in 'iu':
            raise ValueError('triangles must be integer Mx3')
        if raw.min()<0 or raw.max()>=len(x):raise ValueError('triangle index outside cloth')
        tri=raw.astype(np.int64)
        area=.5*np.linalg.norm(np.cross(x[tri[:,1]]-x[tri[:,0]],x[tri[:,2]]-x[tri[:,0]]),axis=1)
        if (area<=1e-16).any():raise ValueError('degenerate cloth triangle')
        density=_positive(areal_density_kg_m2,'areal_density_kg_m2')
        self.stiffness_n_m=_positive(edge_stiffness_n_m,'edge_stiffness_n_m')
        self.mass_kg=np.zeros(len(x))
        for k in range(3):np.add.at(self.mass_kg,tri[:,k],area*density/3)
        if (self.mass_kg<=0).any():raise ValueError('all cloth vertices need positive mass')
        self.edges=np.unique(np.sort(np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]),axis=1),axis=0)
        self.reference_position_m=x.copy()
        self.rest_lengths_m=np.linalg.norm(x[self.edges[:,1]]-x[self.edges[:,0]],axis=1)
        self.position_m=x.copy();self.velocity_m_s=np.zeros_like(x);self.time_s=0.
        self.triangles=tri.copy()
        incident=np.bincount(self.edges.reshape(-1),minlength=len(x))
        # Conservative Gershgorin bound on the linearized spring Laplacian.
        self.max_explicit_dt_s=float(.2*np.sqrt(np.min(self.mass_kg/(incident*self.stiffness_n_m))))

    def elastic_forces(self,positions_m):
        x=np.asarray(positions_m,dtype=float)
        if x.shape!=self.position_m.shape or not np.isfinite(x).all():
            raise ValueError('elastic positions must match the finite cloth mesh')
        delta=x[self.edges[:,1]]-x[self.edges[:,0]]
        length=np.linalg.norm(delta,axis=1)
        if (length<1e-14).any():raise ValueError('collapsed cloth edge')
        extension=length-self.rest_lengths_m
        tension=(self.stiffness_n_m*extension/length)[:,None]*delta
        force=np.zeros_like(x)
        np.add.at(force,self.edges[:,0],tension);np.add.at(force,self.edges[:,1],-tension)
        return force,float(.5*self.stiffness_n_m*np.dot(extension,extension))

    def step(self,dt_s,*,gravity_m_s2=(0.,-9.81,0.),contact=None,external_forces_n=None):
        dt=_positive(dt_s,'dt_s');gravity=_vector(gravity_m_s2,'gravity')
        if dt>self.max_explicit_dt_s:
            raise ValueError(f'dt_s exceeds conservative elastic limit {self.max_explicit_dt_s:g}')
        external=np.zeros_like(self.position_m) if external_forces_n is None else np.asarray(external_forces_n,dtype=float)
        if external.shape!=self.position_m.shape or not np.isfinite(external).all():
            raise ValueError('external forces must match the finite cloth mesh')
        # Work on copies. Failed contact/geometry validation must not advance state.
        force,elastic_before=self.elastic_forces(self.position_m)
        total=force+external+self.mass_kg[:,None]*gravity
        old=self.position_m;oldv=self.velocity_m_s
        v=oldv+dt*total/self.mass_kg[:,None];x=old+dt*v
        result=contact.resolve(x,v,self.mass_kg,dt) if contact is not None else None
        if result is not None:x=result['positions_m'];v=result['velocities_m_s']
        _,elastic_after=self.elastic_forces(x)
        kinetic_before=float(.5*np.sum(self.mass_kg[:,None]*oldv**2))
        kinetic_after=float(.5*np.sum(self.mass_kg[:,None]*v**2))
        potential_change=float(-np.sum(self.mass_kg[:,None]*(x-old)*gravity))
        external_work=float(np.sum(external*(x-old)))
        impulses=np.zeros_like(x) if result is None else result['impulses_ns']
        dissipation=0. if result is None else result['friction_dissipation_j']+result['normal_impact_dissipation_j']
        surface_work=0. if result is None else result['prescribed_surface_work_j']
        momentum_residual=np.sum(self.mass_kg[:,None]*(v-oldv)-dt*total-impulses,axis=0)
        energy_defect=kinetic_after-kinetic_before+elastic_after-elastic_before+potential_change-external_work-surface_work+dissipation
        self.position_m=x.copy();self.velocity_m_s=v.copy();self.time_s+=dt
        return {
            'time_s':self.time_s,'kinetic_energy_j':kinetic_after,'elastic_energy_j':elastic_after,
            'body_reaction_impulse_ns':-impulses.sum(axis=0),
            'momentum_residual_ns':momentum_residual,
            'friction_dissipation_j':0. if result is None else result['friction_dissipation_j'],
            'normal_impact_dissipation_j':0. if result is None else result['normal_impact_dissipation_j'],
            'numerical_energy_defect_j':float(energy_defect),
            'max_position_correction_m':0. if result is None else float(np.linalg.norm(result['position_correction_m'],axis=1).max()),
            'contact':result,
        }
