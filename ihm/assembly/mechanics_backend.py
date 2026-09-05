"""Regional compressible hyperelastic tetrahedra with frictionless tool contact.

Quasistatic reference patch solver, not a whole-body articulated/contact backend.
Material inputs are explicit SI parameters; no empirical coefficient defaults.
"""
import itertools
import math
import time
import numpy as np
from scipy.optimize import minimize


def tetra_box(size_m,divisions):
    size=np.asarray(size_m,float)
    if size.shape!=(3,) or not np.isfinite(size).all() or np.any(size<=0):raise ValueError('Positive box dimensions required')
    if len(divisions)!=3 or any(isinstance(n,bool) or not isinstance(n,(int,np.integer)) or not 1<=n<=64 for n in divisions):raise ValueError('Integer divisions in1..64 required')
    nx,ny,nz=divisions
    x=np.stack(np.meshgrid(*(np.linspace(0,s,n+1) for s,n in zip(size,divisions)),indexing='ij'),axis=-1).reshape(-1,3)
    index=lambda v:int(v[0]*(ny+1)*(nz+1)+v[1]*(nz+1)+v[2])
    t=[]
    for base in itertools.product(range(nx),range(ny),range(nz)):
        base=np.array(base)
        for order in itertools.permutations(range(3)):
            v=base.copy();vertices=[index(v)]
            for axis in order:v=v.copy();v[axis]+=1;vertices.append(index(v))
            t.append(vertices)
    t=np.asarray(t,int)
    negative=np.linalg.det(np.swapaxes(x[t[:,1:]]-x[t[:,0,None]],1,2))<0
    t[negative]=t[negative][:,[0,2,1,3]]
    return x,t


class DeformableRegion:
    def __init__(self,vertices_m,tetrahedra,*,mu_pa,lambda_pa,density_kg_m3):
        self.reference=np.asarray(vertices_m,float).copy();self.tets=np.asarray(tetrahedra)
        if self.reference.ndim!=2 or self.reference.shape[1]!=3 or not np.isfinite(self.reference).all():raise ValueError('Finite reference vertices required')
        if self.tets.ndim!=2 or self.tets.shape[1]!=4 or self.tets.dtype.kind not in 'iu' or not len(self.tets) or self.tets.min()<0 or self.tets.max()>=len(self.reference):raise ValueError('Valid tetrahedral connectivity required')
        self.mu=np.broadcast_to(np.asarray(mu_pa,float),(len(self.tets),)).copy()
        self.lam=np.broadcast_to(np.asarray(lambda_pa,float),(len(self.tets),)).copy()
        rho=np.broadcast_to(np.asarray(density_kg_m3,float),(len(self.tets),))
        if not all(np.isfinite(a).all() for a in [self.mu,self.lam,rho]) or np.any(self.mu<=0) or np.any(self.lam<0) or np.any(rho<=0):raise ValueError('Finite positive material parameters required')
        dm=np.swapaxes(self.reference[self.tets[:,1:]]-self.reference[self.tets[:,0,None]],1,2)
        determinant=np.linalg.det(dm)
        if np.any(determinant<=1e-20):raise ValueError('Positive oriented reference tetrahedra required')
        self.volumes=determinant/6;self.inverse=np.linalg.inv(dm)
        self.nodal_mass=np.zeros(len(self.reference))
        for i in range(4):np.add.at(self.nodal_mass,self.tets[:,i],rho*self.volumes/4)
        self.positions=self.reference.copy()

    def checkpoint(self):return self.positions.copy()

    def restore(self,checkpoint):
        values=np.asarray(checkpoint,float)
        if values.shape!=self.positions.shape or not np.isfinite(values).all():raise ValueError('Invalid mechanical checkpoint')
        self.energy_gradient(values)
        self.positions=values.copy()

    def deformation(self,positions):
        y=np.asarray(positions,float)
        if y.shape!=self.reference.shape or not np.isfinite(y).all():raise ValueError('Invalid deformed positions')
        return np.swapaxes(y[self.tets[:,1:]]-y[self.tets[:,0,None]],1,2)@self.inverse

    def energy_gradient(self,positions,*,trial=False):
        f=self.deformation(positions);j=np.linalg.det(f)
        if not trial and np.any(j<=0):raise ValueError('Inverted mechanical element')
        # Smooth trial-only continuation lets line searches recover from invalid
        # trial tetrahedra. Accepted results must satisfy the original log(J) law.
        floor=.2
        safe=np.maximum(j,floor) if trial else j
        logj=np.log(safe)
        volumetric=-self.mu*logj+self.lam/2*logj**2
        derivative=(-self.mu+self.lam*logj)/safe
        if trial:
            delta=np.minimum(j-floor,0.)
            second=(self.mu+self.lam*(1-np.log(floor)))/floor**2
            volumetric += derivative*delta+.5*second*delta**2
            derivative += second*delta
        cofactor=np.stack((np.cross(f[:,:,1],f[:,:,2]),np.cross(f[:,:,2],f[:,:,0]),np.cross(f[:,:,0],f[:,:,1])),axis=2)
        density=self.mu/2*(np.sum(f*f,axis=(1,2))-3)+volumetric
        piola=self.mu[:,None,None]*f+derivative[:,None,None]*cofactor
        h=self.volumes[:,None,None]*piola@np.swapaxes(self.inverse,1,2)
        local=np.concatenate((-h.sum(axis=2)[:,None,:],np.swapaxes(h,1,2)),axis=1)
        gradient=np.zeros_like(self.reference)
        for k in range(4):np.add.at(gradient,self.tets[:,k],local[:,k])
        return float(self.volumes@density),gradient

    def solve(self,*,fixed_nodes,indenter=None,gravity_m_s2=(0.,0.,0.),force_tolerance_n=1e-7):
        fixed=np.asarray(fixed_nodes)
        gravity=np.asarray(gravity_m_s2,float)
        if fixed.ndim!=1 or fixed.dtype.kind not in 'iu' or len(fixed)<3 or fixed.min()<0 or fixed.max()>=len(self.reference):raise ValueError('At least three valid fixed substrate nodes required')
        if gravity.shape!=(3,) or not np.isfinite(gravity).all():raise ValueError('Finite gravity vector required')
        if isinstance(force_tolerance_n,bool) or not math.isfinite(force_tolerance_n) or force_tolerance_n<=0:raise ValueError('Positive force tolerance required')
        lo=np.full_like(self.reference,-np.inf);hi=np.full_like(self.reference,np.inf)
        lo[fixed]=0.;hi[fixed]=0.;contact=np.array([],int)
        if indenter is not None:
            if not isinstance(indenter,dict) or set(indenter)!={'center_m','radius_m','height_m'}:raise ValueError('Specify circular planar indenter center, radius, height')
            center=np.asarray(indenter['center_m'],float);radius=indenter['radius_m'];height=indenter['height_m']
            if center.shape!=(2,) or not np.isfinite(center).all() or isinstance(radius,bool) or not np.isfinite(radius) or radius<=0 or not np.isfinite(height):raise ValueError('Invalid indenter')
            if height<=self.reference[:,2].min():raise ValueError('Indenter crosses fixed substrate')
            contact=np.flatnonzero((np.linalg.norm(self.reference[:,:2]-center,axis=1)<=radius)&np.isclose(self.reference[:,2],self.reference[:,2].max(),atol=1e-12,rtol=0))
            if not len(contact):raise ValueError('Indenter misses all top-surface nodes; refine mesh')
            hi[contact,2]=height-self.reference[contact,2]
            if np.any(hi<lo):raise ValueError('Conflicting contact and fixed constraints')
        external=self.nodal_mass[:,None]*gravity
        initial=np.clip(self.positions-self.reference,lo,hi)
        # Under zero external load the undeformed configuration is an exact
        # quasistatic minimum. Use it as the initial iterate on tool release;
        # this is not a dynamic relaxation or an imposed transient trajectory.
        if indenter is None and not np.any(external):initial=np.zeros_like(initial)
        self.energy_gradient(self.reference+initial)
        scale=float(np.min(np.ptp(self.reference,axis=0)))
        energy_scale=float(np.mean(self.mu))*scale**3
        def objective(z):
            u=z*scale
            y=self.reference+u.reshape(-1,3)
            energy,gradient=self.energy_gradient(y,trial=True)
            return (energy-float(np.sum(external*u.reshape(-1,3))))/energy_scale,(gradient-external).ravel()*scale/energy_scale
        start=time.perf_counter()
        result=minimize(objective,initial.ravel()/scale,jac=True,method='L-BFGS-B',bounds=list(zip(lo.ravel()/scale,hi.ravel()/scale)),options={'ftol':1e-15,'gtol':force_tolerance_n*.1*scale/energy_scale,'maxiter':3000,'maxls':50,'maxcor':20})
        displacement=result.x.reshape(-1,3)*scale
        y=self.reference+displacement
        if np.linalg.det(self.deformation(y)).min()<=.2:raise RuntimeError('Contact left the original constitutive domain J>0.2')
        energy,elastic_gradient=self.energy_gradient(y);gradient=elastic_gradient-external
        lower=displacement<=lo+1e-12;upper=displacement>=hi-1e-12
        projected=gradient.copy();projected[lower&(gradient>0)]=0;projected[upper&(gradient<0)]=0;projected[fixed]=0
        residual=float(np.max(np.abs(projected)))
        if residual>force_tolerance_n:raise RuntimeError(f'Regional contact failed force tolerance: {residual:g} N ({result.message})')
        self.positions=y.copy()
        force=-float(np.sum(gradient[contact,2][upper[contact,2]])) if len(contact) else 0.
        return {'positions_m':y.tolist(),'tetrahedra':self.tets.tolist(),'elastic_energy_j':energy,'indenter_reaction_n':force,
            'minimum_jacobian':float(np.linalg.det(self.deformation(y)).min()),'force_balance_residual_n':residual,
            'maximum_penetration_m':float(max(0,np.max(y[contact,2]-indenter['height_m']))) if len(contact) else 0.,
            'wall_seconds':time.perf_counter()-start,'iterations':int(result.nit),'solver':'quasistatic tetrahedral neo-Hookean energy minimization; nodal frictionless planar contact',
            'limitations':['Reference patch only; no all-body collision, articulated joints or dynamic inertia.','Nodal tool contact needs mesh refinement; no arbitrary tool surface contact.','Material inputs are explicit; numerical convergence is not empirical skin validation.']}
