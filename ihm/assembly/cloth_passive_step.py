"""Experimental sparse constrained backward-Euler cloth step.

Based on the implicit integration approach in Baraff & Witkin (1998),
https://www.cs.cmu.edu/~baraff/papers/sig98.pdf . This is a separate experimental
solver, not their continuum model or a claim to reproduce their implementation.

Contact geometry here is prescribed translating planes, NOT the production
sphere envelope. Friction and cloth self-collision are unresolved. Initial
positions must already satisfy contacts/material limits. Nonconvex compressed
springs prevent a blanket mathematical passivity guarantee: acceptance checks
actual energy and rejects failed steps atomically. No velocity clipping or
post-projection energy masking occurs. Callers may retry a smaller timestep.
"""
import numpy as np
from scipy import sparse
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, OptimizeResult
from scipy.sparse.linalg import spsolve, MatrixRankWarning
import warnings
from itertools import combinations


class PassiveStepRejected(RuntimeError):
    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics or {}


class ClothPassiveStepper:
    def __init__(self, mesh, max_extension_ratio=1.12, tolerance_m=1e-7,
                 energy_tolerance_j=1e-8, max_iterations=200, allow_generic_fallback=False,
                 bending='mesh'):
        if not np.isfinite(max_extension_ratio) or max_extension_ratio < 1:
            raise ValueError('Finite extension limit at least one required')
        if not np.isfinite(tolerance_m) or tolerance_m<=0 or not np.isfinite(energy_tolerance_j) or energy_tolerance_j<0:
            raise ValueError('Positive feasibility tolerance and nonnegative energy tolerance required')
        if type(max_iterations) is not int or max_iterations<1:
            raise ValueError('Positive integer maximum iterations required')
        self.n = len(mesh.x)
        self.a = np.asarray(mesh.a, int).copy(); self.b = np.asarray(mesh.b, int).copy()
        self.length = np.asarray(mesh.length, float).copy()
        self.stiffness = np.asarray(mesh.stiffness, float).copy()
        self.mass = np.broadcast_to(np.asarray(mesh.mass, float), (self.n,)).copy()
        if np.any(self.mass <= 0) or not np.isfinite(self.mass).all():
            raise ValueError('Finite positive masses required')
        if np.any(self.length <= 0) or not np.isfinite(self.length).all():
            raise ValueError('Finite positive rest lengths required')
        if np.any(self.stiffness<0) or not np.isfinite(self.stiffness).all():
            raise ValueError('Finite nonnegative spring stiffness required')
        self.material = np.asarray(mesh.edges)[:, 2] >= .99
        self.ca = self.a[self.material]; self.cb = self.b[self.material]
        self.limit = self.length[self.material]*max_extension_ratio
        self.ratio = float(max_extension_ratio)
        self.tolerance = float(tolerance_m); self.energy_tolerance = float(energy_tolerance_j)
        self.max_iterations = int(max_iterations)
        self.allow_generic_fallback = bool(allow_generic_fallback)
        # Distance springs leave an isometric crease free. The hinge term has a
        # constant positive semidefinite Hessian, so including it can only make
        # this implicit system better conditioned; it cannot create energy.
        self.bending = getattr(mesh, 'bending', None) if bending == 'mesh' else bending
        if self.bending is not None and getattr(self.bending, 'count', self.n) != self.n:
            raise ValueError('Bending model vertex count differs from the mesh')
        self.bending_hessian = None if self.bending is None else self.bending.hessian()
        self.dof_mass = np.repeat(self.mass, 3)
        self.incidence = sparse.coo_matrix((np.tile([-1.,1.],len(self.ca)),
            (np.repeat(np.arange(len(self.ca)),2),np.column_stack([self.ca,self.cb]).ravel())),
            shape=(len(self.ca),self.n)).tocsr()

    def _potential(self, q, gravity):
        d=q[self.b]-q[self.a]; length=np.linalg.norm(d,axis=1)
        normal=d/np.maximum(length[:,None],1e-12)
        extension=length-self.length
        energy=float(.5*np.sum(self.stiffness*extension**2)-np.sum(self.mass[:,None]*q*gravity))
        force=self.stiffness[:,None]*extension[:,None]*normal
        gradient=-self.mass[:,None]*np.broadcast_to(gravity,q.shape).copy()
        np.add.at(gradient,self.a,-force);np.add.at(gradient,self.b,force)
        if self.bending is not None:
            energy+=self.bending.energy(q);gradient-=self.bending.forces(q)
        return energy,gradient,d,length,normal

    def _hessian(self,q,gravity):
        _,_,d,length,normal=self._potential(q,gravity)
        outer=normal[:,:,None]*normal[:,None,:]
        blocks=self.stiffness[:,None,None]*(outer+(1-self.length/np.maximum(length,1e-12))[:,None,None]*(np.eye(3)-outer))
        rows=[];cols=[];data=[]
        for a,b,sign in ((self.a,self.a,1),(self.b,self.b,1),(self.a,self.b,-1),(self.b,self.a,-1)):
            rows.append(np.broadcast_to((3*a[:,None]+np.arange(3))[...,None],blocks.shape).ravel())
            cols.append(np.broadcast_to((3*b[:,None]+np.arange(3))[:,None,:],blocks.shape).ravel())
            data.append((sign*blocks).ravel())
        assembled=sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(3*self.n,3*self.n)).tocsr()
        return assembled if self.bending_hessian is None else assembled+self.bending_hessian

    def _kkt_solve(self,H,J,grad,c):
        """Solve one KKT system that may have linearly dependent active rows.

        Overlapping contact constraints make the active Jacobian rank
        deficient: a vertex held by several nearly tangent sampled spheres, or
        a material edge at its limit whose two endpoints are already fully
        determined by anchors or contacts. The KKT matrix is then exactly
        singular and the sparse factorization fails, which previously rejected
        the whole step with the misleading advice to retry a smaller timestep.
        A smaller timestep cannot repair a rank deficiency.

        Redundant rows are not dropped by inspection, because deciding which of
        several physically real contacts is 'the' duplicate is arbitrary and
        because inconsistent rows must stay visible. Instead the system is
        solved with dual regularization and refined against the ORIGINAL rows,
        which converges whenever the redundant set is consistent and leaves a
        measurable residual when it is not. Regularization only selects one
        member of a non-unique multiplier set; the total reaction, and hence
        the reported impulse and momentum residual, is unchanged. The caller
        still validates stationarity, feasibility, energy and momentum against
        the unregularized constraints, so an inconsistent set is still rejected.
        """
        size=3*self.n
        rhs=np.concatenate([-grad,-c]) if len(c) else -grad
        system=sparse.bmat([[H,J.T],[J,None]],format='csc') if len(c) else H.tocsc()
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            try:
                answer=spsolve(system,rhs)
                if np.isfinite(answer).all():return answer,'exact'
                failure='Nonfinite KKT solution'
            except (MatrixRankWarning,RuntimeError) as error:failure=str(error)
        if not len(c):return None,failure
        magnitude=max(1.,float(np.max(np.abs(rhs),initial=0.)))
        scale=max(1.,float(abs(H).max()) if H.nnz else 1.)
        best=None
        for exponent in (12,10,8):
            regularized=sparse.bmat([[H,J.T],[J,-(10.**-exponent)*scale*sparse.eye(len(c))]],format='csc')
            trial=np.zeros(size+len(c))
            for _ in range(16):
                with warnings.catch_warnings():
                    warnings.simplefilter('error',MatrixRankWarning)
                    try:correction=spsolve(regularized,rhs-system@trial)
                    except (MatrixRankWarning,RuntimeError):break
                if not np.isfinite(correction).all():break
                trial=trial+correction
                residual=float(np.max(np.abs(rhs-system@trial),initial=0.))
                if best is None or residual<best[0]:best=(residual,trial.copy())
                if residual<=1e-13*magnitude:break
            if best is not None and best[0]<=1e-11*magnitude:
                return best[1],'dual regularization resolved redundant active rows'
        if best is None:return None,failure
        return None,'Redundant active constraints are inconsistent; residual %.3e'%best[0]

    def step(self,mesh,dt_s,gravity=(0,0,-9.81),planes=(),forces_n=None):
        """Advance atomically; return force/work diagnostics on accepted solve.

        Each plane dict has node, unit normal, offset_m at step START,
        velocity_m_s (constant), and optional owner. The allowed halfspace is
        normal dot x >= offset + dt*normal dot velocity. Returned plane impulses
        act on prescribed supports (opposite impulses imparted to cloth).
        """
        h=float(dt_s);g=np.asarray(gravity,float)
        if not np.isfinite(h) or h<=0 or g.shape!=(3,) or not np.isfinite(g).all():
            raise ValueError('Positive finite dt and finite gravity required')
        x=np.asarray(mesh.x,float).copy();v=np.asarray(mesh.v,float).copy()
        if x.shape!=(self.n,3) or v.shape!=x.shape or not np.isfinite(x).all() or not np.isfinite(v).all():
            raise ValueError('Finite mesh state must match cached topology')
        if np.any(np.abs(v[mesh.fixed])>self.tolerance):
            raise ValueError('Fixed anchors require zero velocity')
        if np.any(np.linalg.norm(x[self.b]-x[self.a],axis=1)<1e-10):
            raise PassiveStepRejected('Degenerate initial spring direction; no state changed')
        initial_lengths=np.linalg.norm(x[self.cb]-x[self.ca],axis=1)
        if np.any(initial_lengths>self.limit+self.tolerance):
            raise PassiveStepRejected('Initial material constraints are infeasible; no state changed')
        external=np.zeros_like(x) if forces_n is None else np.asarray(forces_n,float)
        if external.shape!=x.shape or not np.isfinite(external).all():
            raise ValueError('External forces require finite N by 3 values')
        pred=x+h*v;constraints=[];metadata=[];row_data=[];row_indices=[];col_indices=[];lower=[];upper=[]
        for plane in planes:
            node=plane['node'];normal=np.asarray(plane['normal'],float);velocity=np.asarray(plane.get('velocity_m_s',[0,0,0]),float)
            offset=float(plane['offset_m'])
            if type(node) is not int or not 0<=node<self.n or normal.shape!=(3,) or velocity.shape!=(3,) or not np.isfinite(normal).all() or not np.isfinite(velocity).all() or not np.isfinite(offset) or abs(np.linalg.norm(normal)-1)>1e-8:
                raise ValueError('Plane requires valid node, unit normal and finite offset/velocity')
            if np.dot(normal,x[node])<offset-self.tolerance:
                raise PassiveStepRejected('Initial contact geometry is infeasible; no state changed')
            row=len(lower);row_indices.extend([row]*3);col_indices.extend((3*node+np.arange(3)).tolist());row_data.extend(normal.tolist())
            lower.append(offset+h*np.dot(normal,velocity));upper.append(np.inf)
            metadata.append({'node':node,'normal':normal,'velocity':velocity,'owner':plane.get('owner',str(row)),'kind':'plane'})
        for node in np.flatnonzero(mesh.fixed):
            for axis in range(3):
                normal=np.eye(3)[axis];row=len(lower);row_indices.append(row);col_indices.append(3*int(node)+axis);row_data.append(1.)
                lower.append(x[node,axis]);upper.append(x[node,axis]);metadata.append({'node':int(node),'normal':normal,'velocity':np.zeros(3),'owner':'fixed:'+str(node),'kind':'anchor'})
        matrix=sparse.coo_matrix((row_data,(row_indices,col_indices)),shape=(len(lower),3*self.n)).tocsr()
        if len(lower):constraints.append(LinearConstraint(matrix,lower,upper))
        if len(self.ca):
            def lengths_squared(flat):
                q=flat.reshape(-1,3);return np.sum((q[self.cb]-q[self.ca])**2,axis=1)
            def jacobian(flat):
                q=flat.reshape(-1,3);d=2*(q[self.cb]-q[self.ca]);rows=np.repeat(np.arange(len(d)),6)
                cols=np.concatenate((3*self.ca[:,None]+np.arange(3),3*self.cb[:,None]+np.arange(3)),axis=1).ravel()
                return sparse.coo_matrix((np.concatenate((-d,d),axis=1).ravel(),(rows,cols)),shape=(len(d),3*self.n)).tocsr()
            def constraint_hessian(flat,multipliers):
                return sparse.kron(self.incidence.T@sparse.diags(2*multipliers)@self.incidence,sparse.eye(3),format='csr')
            constraints.append(NonlinearConstraint(lengths_squared,-np.inf,self.limit**2,jac=jacobian,hess=constraint_hessian,keep_feasible=False))
        def objective(flat):
            q=flat.reshape(-1,3);potential,gradient,*_=self._potential(q,g)
            potential-=float(np.sum(external*q));gradient-=external
            delta=q-pred
            return float(.5*np.sum(self.mass[:,None]*delta**2)+h*h*potential),(self.mass[:,None]*delta+h*h*gradient).ravel()
        def hessian(flat):return sparse.diags(self.dof_mass)+h*h*self._hessian(flat.reshape(-1,3),g)
        potential0=self._potential(x,g)[0];energy0=potential0+.5*np.sum(self.mass[:,None]*v**2)
        active_failure={};redundant_solves=[]
        plane_groups={}
        for index,meta in enumerate(metadata):plane_groups.setdefault(meta['node'],[]).append(index)
        def independent_contacts(active, reference, all_active=False):
            selected=active.copy()
            for node,indices in plane_groups.items():
                current=[i for i in indices if selected[i]]
                if not current:continue
                normals=np.array([metadata[i]['normal'] for i in current])
                if not all_active and len(current)<=3 and np.linalg.matrix_rank(normals,tol=1e-10)==len(current):continue
                # Closest-point QP onto ALL halfspaces for this vertex. In 3D
                # its optimum has at most3 independent active normals. Keep all
                # original inequalities in the global feasibility validation.
                normals=np.array([metadata[i]['normal'] for i in indices]);offsets=np.asarray(lower)[indices]
                equality_local=np.array([metadata[i]['kind']=='anchor' for i in indices])
                target=reference.reshape(-1,3)[node];best=None
                for count in range(min(3,len(indices))+1):
                    for subset in combinations(range(len(indices)),count):
                        if not set(np.flatnonzero(equality_local)).issubset(subset):continue
                        if count:
                            N=normals[list(subset)];gram=N@N.T
                            if np.linalg.matrix_rank(gram,tol=1e-11)<count:continue
                            dual=np.linalg.solve(gram,offsets[list(subset)]-N@target)
                            if np.any(dual[~equality_local[list(subset)]] < -1e-10):continue
                            trial=target+N.T@dual
                        else:trial=target
                        if np.any(normals@trial<offsets-1e-10):continue
                        cost=float(np.sum((trial-target)**2))
                        if best is None or cost<best[0]:best=(cost,subset)
                if best is None:
                    active_failure.update(reason='Per-vertex contact halfspaces infeasible or degenerate',node=node)
                    return None
                selected[indices]=False
                selected[[indices[i] for i in best[1]]]=True
            return selected
        def active_newton():
            # Sparse active-set Newton avoids interior-barrier forces on a cloth
            # already resting against its support. This is a local nonlinear
            # solve; difficult contact configurations can still fail atomically.
            qflat=x.ravel().copy();active_linear=np.zeros(len(lower),bool)
            equality=np.isfinite(upper) if len(lower) else np.zeros(0,bool)
            if len(lower):active_linear=equality | (matrix@pred.ravel()<=np.asarray(lower)+self.tolerance)
            active_linear=independent_contacts(active_linear,pred.ravel(),all_active=True)
            if active_linear is None:return None
            active_edges=np.linalg.norm(pred[self.cb]-pred[self.ca],axis=1)>=self.limit
            lmult=np.zeros(len(lower));emult=np.zeros(len(self.ca));optimality=np.inf;gap_error=np.inf
            for iteration in range(min(self.max_iterations,60)):
                active_linear=independent_contacts(active_linear,qflat-objective(qflat)[1]/self.dof_mass)
                if active_linear is None:return None
                li=np.flatnonzero(active_linear);ei=np.flatnonzero(active_edges)
                blocks=[];residual=[]
                if len(li):blocks.append(matrix[li]);residual.append((matrix@qflat-np.asarray(lower))[li])
                if len(ei):blocks.append(jacobian(qflat)[ei]);residual.append((lengths_squared(qflat)-self.limit**2)[ei])
                J=sparse.vstack(blocks,format='csr') if blocks else sparse.csr_matrix((0,3*self.n))
                c=np.concatenate(residual) if residual else np.zeros(0)
                emult[~active_edges]=0
                grad=objective(qflat)[1];H=hessian(qflat)
                if len(ei):H=H+constraint_hessian(qflat,emult)
                answer,note=self._kkt_solve(H,J,grad,c)
                if note!='exact':redundant_solves.append({'iteration':iteration,'note':note})
                if answer is None:
                    active_failure.update(reason=note,iteration=iteration,active_planes=len(li),active_edges=len(ei),planes=[dict(index=int(i),node=metadata[i]["node"],normal=metadata[i]["normal"].tolist(),lower=float(lower[i])) for i in li]);return None
                delta=answer[:3*self.n];dual=answer[3*self.n:]
                lmult[:]=0;emult[:]=0
                if len(li):lmult[li]=dual[:len(li)]
                if len(ei):emult[ei]=dual[len(li):]
                bad_linear=active_linear & (lmult>1e-12) & ~equality
                bad_edges=active_edges & (emult < -1e-12)
                if np.any(bad_linear) or np.any(bad_edges):
                    active_linear[bad_linear]=False;active_edges[bad_edges]=False
                    continue
                qflat+=delta
                # Wrong-sign inequality reactions leave the active set. Add
                # newly violated constraints before any convergence decision.
                old_l=active_linear.copy();old_e=active_edges.copy()
                if len(lower):
                    gap=matrix@qflat-np.asarray(lower)
                    active_linear=(active_linear & ((lmult<=1e-12)|equality)) | (gap < -self.tolerance) | equality
                if len(self.ca):
                    excess=lengths_squared(qflat)-self.limit**2
                    active_edges=(active_edges & (emult>=-1e-12)) | (excess>self.tolerance*np.maximum(self.limit,1e-12))
                stationarity=objective(qflat)[1]
                if len(lower):stationarity=stationarity+matrix.T@lmult
                if len(self.ca):stationarity=stationarity+jacobian(qflat).T@emult
                optimality=float(np.max(np.abs(stationarity),initial=0.))
                gap_error=float(np.max(np.abs(c),initial=0.))
                if np.array_equal(old_l,active_linear) and np.array_equal(old_e,active_edges) and optimality<1e-11 and gap_error<1e-10 and np.max(np.abs(delta),initial=0.)<1e-8:
                    multipliers=[]
                    if len(lower):multipliers.append(lmult.copy())
                    if len(self.ca):multipliers.append(emult.copy())
                    return OptimizeResult(x=qflat,success=True,message='Sparse active-set Newton converged',nit=iteration+1,v=multipliers,optimality=optimality)
            active_failure.update(reason="Active Newton iteration limit",iteration=iteration,active_planes=len(li),active_edges=len(ei),optimality=optimality,step_max=float(np.max(np.abs(delta))),constraint_residual=gap_error)
            return None
        result=active_newton()
        if result is None and not self.allow_generic_fallback:
            raise PassiveStepRejected("Active contact solve did not converge; retry a smaller timestep",{
                "accepted":False,"solver_success":False,"active_set_failure":active_failure,"redundant_constraint_solves":list(redundant_solves)})
        if result is None:
            result=minimize(objective,x.ravel(),jac=True,hess=hessian,method='trust-constr',constraints=constraints,
            options={'maxiter':self.max_iterations,'gtol':1e-15,'xtol':1e-12,'barrier_tol':1e-12,
                     'initial_tr_radius':.01,'initial_barrier_parameter':1e-16,'initial_barrier_tolerance':1e-12,'sparse_jacobian':True})
        q=result.x.reshape(-1,3);next_v=(q-x)/h
        reaction=np.zeros((len(metadata),3));work=0.
        if len(metadata):
            for i,(meta,multiplier) in enumerate(zip(metadata,result.v[0])):
                reaction[i]=float(multiplier)/h*meta['normal']
                work-=float(reaction[i]@meta['velocity'])
        potential1=self._potential(q,g)[0];energy1=potential1+.5*np.sum(self.mass[:,None]*next_v**2)
        lengths=np.linalg.norm(q[self.cb]-q[self.ca],axis=1)
        violation=float(np.maximum(0,lengths-self.limit).max(initial=0.))
        if len(lower):violation=max(violation,float(np.maximum(0,np.asarray(lower)-matrix@result.x).max(initial=0.)),float(np.maximum(0,matrix@result.x-np.asarray(upper)).max(initial=0.)))
        external_work=float(np.sum(external*(q-x)))
        momentum=np.sum(self.mass[:,None]*(next_v-v-h*g)-h*external,axis=0)+reaction.sum(axis=0)
        diagnostics={'accepted':False,'iterations':int(result.nit),'solver_success':bool(result.success),'solver_message':str(result.message),
            'active_set_failure':active_failure,'redundant_constraint_solves':list(redundant_solves),'constraint_violation_m':violation,'stationarity_residual':float(result.optimality),
            'energy_before_j':float(energy0),'energy_after_j':float(energy1),'prescribed_boundary_work_j':work,'external_force_work_j':external_work,
            'energy_excess_j':float(energy1-energy0-work-external_work),'momentum_residual_ns':momentum.tolist(),
            'max_extension_ratio':float((lengths/self.length[self.material]).max(initial=0.)),
            'plane_reactions':[{'owner':m['owner'],'node':m['node'],'kind':m['kind'],'impulse_ns':r.tolist()} for m,r in zip(metadata,reaction)],
            'bending':None if self.bending is None else self.bending.state(q,getattr(mesh,'areal_density_kg_m2',None)),
            'scope':'Experimental constrained implicit Euler, actual spring plus hinge bending energy, prescribed translating planes only; no friction/self-collision/nonlinear cover geometry. Failed feasibility, convergence or passivity leaves state unchanged.'}
        if not result.success or not np.isfinite(result.optimality) or result.optimality>1e-9 or not np.isfinite(momentum).all() or np.linalg.norm(momentum)>1e-6 or not np.isfinite(q).all() or not np.isfinite(diagnostics['energy_excess_j']) or violation>self.tolerance or diagnostics['energy_excess_j']>self.energy_tolerance:
            raise PassiveStepRejected('Implicit step rejected: convergence, feasibility or passivity failed',diagnostics)
        mesh.x[:]=q;mesh.v[:]=next_v;diagnostics['accepted']=True
        return diagnostics
