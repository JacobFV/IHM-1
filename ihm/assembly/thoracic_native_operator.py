"""Exact replacement operator for a native tree's torso inertia contribution.

This computes an enlarged mass/bias operator, never an additive force fed to
an unchanged native acceleration solve. Native constraint projection and
state extraction/integration remain the caller's explicit responsibilities.
"""
import numpy as np
from .thoracic_cervical_metric import skew


def rigid_parent_operator(moment_ledger,twist):
    m=float(moment_ledger['mass_kg']);h=np.asarray(moment_ledger['H_first_moment_kg_m']);q=np.asarray(moment_ledger['Q_second_moment_kg_m2']);center=h/m
    second_c=q-m*np.outer(center,center);inertia=np.trace(second_c)*np.eye(3)-second_c;a=np.c_[np.eye(3),-skew(center)]
    matrix=m*a.T@a;matrix[3:,3:]+=inertia;v=twist[:3];omega=twist[3:]
    bias=m*a.T@(np.cross(omega,v)+np.cross(omega,np.cross(omega,center)));bias[3:]+=np.cross(omega,inertia@omega)
    return matrix,bias


def native_frame_to_body_twist_map(rotation,world_frame_jacobian,world_frame_bias,native_velocity):
    """Convert Simbody angular-first world frame data to linear-first body data.

Native frame bias is physical spatial acceleration with native udot=0.
The body linear-speed bias subtracts omega cross v; angular bias needs only R^T.
    """
    r=np.asarray(rotation,float);j=np.asarray(world_frame_jacobian,float);a=np.asarray(world_frame_bias,float);u=np.asarray(native_velocity,float)
    if r.shape!=(3,3) or not np.allclose(r.T@r,np.eye(3),atol=1e-10,rtol=0) or np.linalg.det(r)<=0 or j.shape!=(6,len(u)) or a.shape!=(6,) or u.ndim!=1 or not all(np.isfinite(x).all() for x in [r,j,a,u]):raise ValueError('Invalid native frame snapshot')
    body_map=np.r_[r.T@j[3:],r.T@j[:3]];twist=body_map@u
    beta=np.r_[r.T@a[3:]-np.cross(twist[3:],twist[:3]),r.T@a[:3]]
    return body_map,beta


def compose_native_operator(metric,native_mass,native_bias,parent_map,parent_bias,native_velocity,internal_q,internal_speed,*,native_model_sha256,pressure_pa=0.):
    """Replace the original torso within same-state, unconstrained native operators.

Inputs must share one native state, realized at Velocity stage. The native
inverse-dynamics bias must use zero applied forces and zero udot. This function
does not silently discard native constraints or prescribe native accelerations.
    """
    if native_model_sha256!=metric.plan['target_native_model_identity']['sha256']:raise ValueError('Native identity differs from composition target')
    native_u=np.asarray(native_velocity,float);n=len(native_u);matrix=np.asarray(native_mass,float);bias=np.asarray(native_bias,float);j=np.asarray(parent_map,float);beta=np.asarray(parent_bias,float);speed=np.asarray(internal_speed,float)
    if native_u.ndim!=1 or matrix.shape!=(n,n) or bias.shape!=(n,) or j.shape!=(6,n) or beta.shape!=(6,) or speed.shape!=(32,) or not all(np.isfinite(x).all() for x in [native_u,matrix,bias,j,beta,speed]) or not np.allclose(matrix,matrix.T,atol=1e-10,rtol=0):raise ValueError('Invalid native mass/bias/parent snapshot')
    twist=j@native_u;coupled=metric.evaluate(internal_q,np.r_[twist,speed],pressure_pa);old_matrix,old_bias=rigid_parent_operator(metric.plan['current_native_torso'],twist)
    pullback=np.zeros((38,n+32));pullback[:6,:n]=j;pullback[6:,n:]=np.eye(32);frame_bias=np.r_[beta,np.zeros(32)]
    enlarged=np.zeros((n+32,n+32));enlarged[:n,:n]=matrix-j.T@old_matrix@j;enlarged+=pullback.T@coupled['mass_matrix']@pullback
    inertial=np.zeros(n+32);inertial[:n]=bias-j.T@(old_matrix@beta+old_bias);inertial+=pullback.T@(coupled['mass_matrix']@frame_bias+coupled['inertial_bias'])
    force=pullback.T@coupled['pressure_port']['generalized_force'];velocity=np.r_[native_u,speed]
    return {'mass_matrix':enlarged,'inertial_bias':inertial,'geometric_pressure_force':force,
            'pressure_power_W':float(force@velocity),'composite_velocity_pullback':pullback,
            'composite_state_evaluation':coupled,'eliminated_internal_indices':[n+k for k in sorted(metric.thorax.locked)],
            'native_constraints_projected':False,'native_acceleration_solved':False,
            'scope':'Full replacement operator; caller must compose native constraints/Ndot, other forces and global state integration. Never apply as a force to the original mass operator.'}
