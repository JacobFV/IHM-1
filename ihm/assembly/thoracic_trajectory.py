"""Bounded explicit midpoint free trajectories with SO(3) parent updates.

Second-order approximation, not an energy-exact or symplectic integrator.
This module has no native physiology, recoil, gravity, or external forcing.
"""
import numpy as np
from .thoracic_free_dynamics import pose_rates


def rotation_increment(rotation_vector):
    """Rodrigues exponential with continuous small-angle coefficients."""
    w=np.asarray(rotation_vector,float)
    if w.shape!=(3,) or not np.isfinite(w).all():raise ValueError('Finite rotation vector required')
    theta=np.linalg.norm(w);k=np.cross(w,np.eye(3)).T
    a=np.sinc(theta/np.pi);b=.5*np.sinc(theta/(2*np.pi))**2
    return np.eye(3)+a*k+b*(k@k)


class FreeThoraxTrajectory:
    def __init__(self,dynamics):
        self.dynamics=dynamics;self.model=dynamics.model

    def validate(self,state):
        q=self.model.coordinates(state['q']).copy();u=np.asarray(state['velocity'],float).copy()
        r=np.asarray(state['rotation'],float).copy();t=np.asarray(state['translation_m'],float).copy()
        pose_rates(r,u)
        if t.shape!=(3,) or not np.isfinite(t).all():raise ValueError('Finite world translation required')
        if any(abs(u[6+k])>1e-14 for k in self.model.locked):raise ValueError('Unresolved rib speeds are constrained')
        return {'q':q,'velocity':u,'rotation':r,'translation_m':t}

    def step(self,state,dt):
        """One free Lie midpoint step; never mutate or partially publish input.

Geometry is validated at midpoint and endpoint, including cavity facet/domain
checks. These do not prove absence of all material self-intersections.
        """
        if isinstance(dt,(bool,np.bool_)) or not np.isscalar(dt) or not np.isfinite(dt) or dt<=0:raise ValueError('Finite positive time step required')
        start=self.validate(state);q=start['q'];u=start['velocity'];r=start['rotation'];t=start['translation_m']
        initial=self.dynamics.evaluate(q,u)
        half_u=u+.5*dt*initial['velocity_derivative'];half_q=q+.5*dt*u[6:]
        half_r=r@rotation_increment(.5*dt*u[3:6])
        self.model.cavity(half_q)
        midpoint=self.dynamics.evaluate(half_q,half_u)
        end={'q':q+dt*half_u[6:],'velocity':u+dt*midpoint['velocity_derivative'],
             'translation_m':t+dt*(half_r@half_u[:3]),'rotation':r@rotation_increment(dt*half_u[3:6])}
        end=self.validate(end);self.model.cavity(end['q'])
        return end

    def invariants(self,state):
        """Energy and momenta about one fixed world origin for a free trajectory."""
        state=self.validate(state);u=state['velocity'];r=state['rotation'];t=state['translation_m']
        kinetic=self.model.kinetic(state['q'],u);momentum=kinetic['mass_matrix']@u
        linear=r@momentum[:3];angular=r@momentum[3:6]+np.cross(t,linear)
        return {'kinetic_energy_J':kinetic['kinetic_energy_J'],'mass_kg':kinetic['mass_kg'],
                'world_linear_momentum_kg_m_per_s':linear,'world_angular_momentum_kg_m2_per_s':angular}
