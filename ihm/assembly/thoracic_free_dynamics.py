"""Conservative free-parent thorax equations for the source-only material map.

u=(R.T tdot, omega_body, qdot_internal). Parent angular velocity is a
quasi-velocity, not an Euler-angle derivative. No native coupling, recoil,
activation, gravity, numerical integration, or additional mass is installed.
"""
import numpy as np
from .thoracic_mechanism import spatial_jacobian


def pose_rates(rotation, velocity):
    """Return world origin derivative and SO(3) matrix derivative R [omega]x."""
    r=np.asarray(rotation,float);u=np.asarray(velocity,float)
    if r.shape!=(3,3) or not np.isfinite(r).all() or not np.allclose(r.T@r,np.eye(3),atol=1e-10,rtol=0) or not np.isclose(np.linalg.det(r),1,atol=1e-10,rtol=0):
        raise ValueError('Parent rotation must be a finite proper orthogonal matrix')
    if u.shape!=(32,) or not np.isfinite(u).all():raise ValueError('Expected 32 finite generalized speeds')
    skew=np.cross(u[3:6],np.eye(3)).T
    return r@u[:3],r@skew


class ThoracicFreeDynamics:
    """Mass-preserving d'Alembert projection with fixed ambiguous rib locks."""
    def __init__(self,mechanism):
        self.model=mechanism

    def directional_curvature(self,jacobian,internal_speed):
        """Exact H(q)[s,s] for independent rotations and linear material modes.

Each retained anchor map has constant interpolation weights. Every nonlinear
term depends on exactly one rib angle; its derivative is axis cross J_j.
There are no mixed internal Hessians in this particular source recipe.
        """
        result=np.zeros(jacobian.shape[:2])
        for index,hinge in self.model.hinges.items():
            if index not in self.model.locked:
                result+=np.cross(hinge['axis'],jacobian[:,:,index])*internal_speed[index]**2
        return result

    def evaluate(self,q,velocity,generalized_force=None):
        """Evaluate M, inertial bias, and free or explicitly forced u derivative.

The material acceleration is R(A udot+c), with c=omega cross v +
omega cross (omega cross x)+2 omega cross (J s)+H[s,s]. Constant reference
triangle masses use the same positive degree-two quadrature as kinetic().
        """
        m=self.model;q=m.coordinates(q);u=np.asarray(velocity,float)
        kinetic=m.kinetic(q,u)  # Also validates finite speeds and exact locks.
        v=u[:3];omega=u[3:6];s=u[6:]
        core=m.anatomy['replace_reduced_torso_with'];center=np.asarray(core['center_m']);inertia=np.asarray(core['inertia_kg_m2'])
        a=spatial_jacobian(center[None,:],np.zeros((1,3,26)))[0]
        core_bias=np.cross(omega,v)+np.cross(omega,np.cross(omega,center))
        bias=core['mass_kg']*a.T@core_bias
        bias[3:6]+=np.cross(omega,inertia@omega)
        bary=np.full((3,3),1/6);np.fill_diagonal(bary,2/3)
        for ident,material in m.material.items():
            x,j=m.material_state(ident,q);h=self.directional_curvature(j,s)
            internal_velocity=np.einsum('ijk,k->ij',j,s)
            convective=np.cross(omega,v)+np.cross(omega,np.cross(omega,x))+2*np.cross(omega,internal_velocity)+h
            for start in range(0,len(material['faces']),1500):
                faces=material['faces'][start:start+1500]
                weights=np.repeat(material['triangle_mass'][start:start+1500]/3,3)
                p=np.einsum('ab,tbc->tac',bary,x[faces]).reshape(-1,3)
                pj=np.einsum('ab,tbcd->tacd',bary,j[faces]).reshape(-1,3,26)
                c=np.einsum('ab,tbc->tac',bary,convective[faces]).reshape(-1,3)
                jac=spatial_jacobian(p,pj)
                bias+=np.einsum('ijk,ij,i->k',jac,c,weights)
        result={**kinetic,'velocity':u.copy(),'internal_coordinate_derivative':s.copy(),'inertial_bias':bias}
        return self.solve(result,np.zeros(32) if generalized_force is None else generalized_force)

    def solve(self,evaluated,generalized_force):
        """Reuse an evaluated state for another explicit load, without integration.

Loads on locked coordinates are reported as ideal constraint reactions and do
no work. Returned parent acceleration is a body-speed derivative; world origin
acceleration is R(udot_linear+omega cross v), not simply R udot_linear.
        """
        force=np.asarray(generalized_force,float)
        if force.shape!=(32,) or not np.isfinite(force).all():raise ValueError('Expected 32 finite generalized forces')
        matrix=evaluated['mass_matrix'];bias=evaluated['inertial_bias'];active=self.model.active
        acceleration=np.zeros(32)
        acceleration[active]=np.linalg.solve(matrix[np.ix_(active,active)],(force-bias)[active])
        reaction=matrix@acceleration+bias-force
        return {**evaluated,'velocity_derivative':acceleration,'generalized_force':force.copy(),
                'constraint_generalized_reaction':reaction,'external_power_W':float(evaluated['velocity']@force),
                'equation_residual_active':reaction[active]}
