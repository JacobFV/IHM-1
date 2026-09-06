"""Bounded stationary geodesic shooting on the pinned BIClong ellipsoid.

This source-only reference solves physical endpoint tangency and surface
geodesic equations. It never adds an artificial generalized force correction.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
AXES=np.array([.027559229367297551,.02204738349383804,.02204738349383804])

def solve_path(p1,p2,seed_contact,seed_arc,*,axes=AXES):
    p=np.asarray(p1,float);q=np.asarray(p2,float);a=np.asarray(axes,float);seed=np.asarray(seed_contact,float)
    if any(x.shape!=(3,) for x in (p,q,a,seed)) or not np.isfinite([p,q,a,seed]).all() or np.any(a<=0):raise ValueError('Finite ellipsoid geometry required')
    if min(np.sum((p/a)**2),np.sum((q/a)**2))<=1+1e-7:raise ValueError('Outside endpoint domain required')
    if not 1e-6<seed_arc<np.pi*a.min():raise ValueError('Short native branch required')
    y0=seed/a;y0/=np.linalg.norm(y0);axis=np.eye(3)[np.argmin(abs(y0))];e1=np.cross(y0,axis);e1/=np.linalg.norm(e1);e2=np.cross(y0,e1)
    d=1/(a*a);scale=a.min();evaluations=0
    def integrate(z,dense=False):
        nonlocal evaluations
        evaluations+=1
        y=y0+z[0]*e1+z[1]*e2;y/=np.linalg.norm(y);x=a*y
        normal=d*x;normal/=np.linalg.norm(normal);incoming=x-p;incoming/=np.linalg.norm(incoming)
        tangent=incoming-normal*(normal@incoming);tangent/=np.linalg.norm(tangent)
        def rhs(_,state):
            xx=state[:3];v=state[3:];n=d*xx
            return np.r_[v,-(v@(d*v))/(n@n)*n]
        arc=z[2]*scale
        result=solve_ivp(rhs,(0,arc),np.r_[x,tangent],method='DOP853',rtol=2e-11,atol=2e-13,max_step=arc/12,dense_output=dense)
        if not result.success:raise ValueError('Geodesic integration failed')
        end=result.y[:3,-1];v=result.y[3:,-1]
        residual=np.r_[normal@incoming,(end+z[3]*scale*v-q)/scale]
        return residual,result,x,tangent
    initial=np.array([0.,0.,seed_arc/scale,np.linalg.norm(q-seed)/scale])
    solution=least_squares(lambda z:integrate(z)[0],initial,bounds=([-.4,-.4,1e-6/scale,1e-6/scale],[.4,.4,np.pi*a.min()/scale,10.]),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=35,diff_step=1e-5)
    residual,result,x,tangent=integrate(solution.x,True);arc=solution.x[2]*scale;end=result.y[:3,-1];samples=result.sol(np.linspace(0,arc,65)).T
    if not solution.success or np.max(abs(residual))>1e-9:raise ValueError('Stationary shooting did not converge')
    if np.max(samples[:,1])>=-1e-8:raise ValueError('Native minus-y quadrant branch violated')
    surface=np.max(abs(np.sum((samples[:,:3]/a)**2,axis=1)-1));speed=np.max(abs(np.linalg.norm(samples[:,3:],axis=1)-1))
    if max(surface,speed)>1e-8:raise ValueError('Geodesic drift exceeds domain')
    length=np.linalg.norm(x-p)+arc+np.linalg.norm(q-end)
    gradients=np.array([(p-x)/np.linalg.norm(p-x),(q-end)/np.linalg.norm(q-end)])
    return {'length_m':float(length),'arc_m':float(arc),'r1':x,'r2':end,'endpoint_gradients':gradients,'surface_points':samples[:,:3],
            'shooting_residual_max':float(max(abs(residual))),'surface_residual_max':float(surface),'unit_speed_residual_max':float(speed),
            'native_seed_contact_shift_m':float(np.linalg.norm(x-seed)),'rhs_shooting_integrations':evaluations,'branch':'native-seeded short geodesic, entire sampled arc in minus-y quadrant',
            'scope':'Converged local stationary branch; global shortest-route uniqueness and full source-domain coverage not certified'}
