"""Bounded-memory face-interior cloth contact against a moving source surface.

Prescribed boundary within a partitioned interval. Reaction impulses and work
are returned for the finite native body owner; no target mass is invented.
Closest edge/vertex cases and continuous/self contact remain unresolved.
"""
import numpy as np
from scipy.spatial import cKDTree

class MovingSurfaceContact:
    def __init__(self,start_m,end_m,triangles,duration_s,*,friction_static,friction_kinetic,search_distance_m=.004,max_candidate_pairs=100000):
        a=np.asarray(start_m,float);b=np.asarray(end_m,float);t=np.asarray(triangles)
        if a.ndim!=2 or a.shape[1]!=3 or b.shape!=a.shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('Finite compatible surface endpoints required')
        if t.ndim!=2 or t.shape[1]!=3 or not len(t) or t.dtype.kind not in 'iu' or t.min()<0 or t.max()>=len(a):raise ValueError('Invalid source triangles')
        for value in (duration_s,search_distance_m):
            if isinstance(value,bool) or not np.isfinite(value) or value<=0:raise ValueError('Positive interval/search distance required')
        if any(isinstance(v,(bool,np.bool_)) for v in (friction_kinetic,friction_static)) or not 0<=friction_kinetic<=friction_static or not np.isfinite([friction_kinetic,friction_static]).all():raise ValueError('Finite Coulomb coefficients required')
        if isinstance(max_candidate_pairs,(bool,np.bool_)) or not isinstance(max_candidate_pairs,(int,np.integer)) or max_candidate_pairs<1:raise ValueError('Positive integer candidate budget required')
        self.start=a.copy();self.delta=b-a;self.triangles=t.copy();self.duration=duration_s;self.velocity=self.delta/duration_s
        self.mu_s=friction_static;self.mu_k=friction_kinetic;self.distance=search_distance_m;self.max_pairs=max_candidate_pairs;self.fraction=0.
        vertices=a[t];centers=vertices.mean(1);self.radius=float(np.linalg.norm(vertices-centers[:,None],axis=2).max())
        self.motion_bound=float(np.linalg.norm(self.delta,axis=1).max());self.tree=cKDTree(centers)
    def resolve(self,positions_m,velocities_m_s,mass_kg,dt_s):
        x=np.array(positions_m,float,copy=True);v=np.array(velocities_m_s,float,copy=True);m=np.asarray(mass_kg,float)
        if x.ndim!=2 or x.shape[1]!=3 or v.shape!=x.shape or m.shape!=(len(x),) or not np.isfinite(x).all() or not np.isfinite(v).all() or not np.isfinite(m).all() or np.any(m<=0):raise ValueError('Invalid finite cloth state')
        if isinstance(dt_s,(bool,np.bool_)) or not np.isfinite(dt_s) or dt_s<=0 or dt_s>self.duration:raise ValueError('Contact substep outside positive parent interval')
        if not 0<=self.fraction<=1:raise ValueError('Surface interpolation outside interval')
        impulses=np.zeros_like(x);reaction=np.zeros_like(self.start);correction=np.zeros_like(x);normal_loss=0.;friction_loss=0.;work=0.;count=0;unresolved=0;max_pairs_seen=0
        for lower in range(0,len(x),128):
            upper=min(lower+128,len(x));radius=self.radius+self.distance+self.motion_bound
            lengths=self.tree.query_ball_point(x[lower:upper],radius,return_length=True);total=int(lengths.sum());max_pairs_seen=max(max_pairs_seen,total)
            if total>self.max_pairs:raise MemoryError('Contact candidate budget exceeded; reduce geometry-query batch before retrying')
            if total==0:continue
            queries=self.tree.query_ball_point(x[lower:upper],radius,return_sorted=True)
            fi=np.concatenate([q for q in queries if len(q)]);ni=np.repeat(np.arange(lower,upper),lengths)
            nodes=self.triangles[fi];tri=self.start[nodes]+self.fraction*self.delta[nodes]
            e0=tri[:,1]-tri[:,0];e1=tri[:,2]-tri[:,0];normal=np.cross(e0,e1);norm=np.linalg.norm(normal,axis=1)
            if np.any(norm<1e-16):raise ValueError('Collapsed moving source contact triangle')
            normal/=norm[:,None];r=x[ni]-tri[:,0];gap=np.sum(r*normal,axis=1)
            aa=np.sum(e0*e0,axis=1);bb=np.sum(e0*e1,axis=1);cc=np.sum(e1*e1,axis=1);dd=np.sum(r*e0,axis=1);ee=np.sum(r*e1,axis=1);den=aa*cc-bb*bb
            w1=(cc*dd-bb*ee)/den;w2=(aa*ee-bb*dd)/den;weights=np.column_stack((1-w1-w2,w1,w2));inside=np.all(weights>=-1e-10,axis=1)
            d2=np.where(inside,gap*gap,np.inf)
            for a,b in ((0,1),(1,2),(2,0)):
                edge=tri[:,b]-tri[:,a];u=np.clip(np.sum((x[ni]-tri[:,a])*edge,axis=1)/np.sum(edge*edge,axis=1),0,1)
                q=tri[:,a]+u[:,None]*edge;d2=np.minimum(d2,np.sum((x[ni]-q)**2,axis=1))
            starts=np.r_[0,np.cumsum(lengths[lengths>0])[:-1]];mins=np.minimum.reduceat(d2,starts)
            matching=np.flatnonzero(d2==np.repeat(mins,lengths[lengths>0]));_,first=np.unique(ni[matching],return_index=True);k=matching[first]
            active=(d2[k]<=self.distance**2)&(gap[k]<=1e-12);unresolved+=int(np.count_nonzero(active&~inside[k]));k=k[active&inside[k]]
            if not len(k):continue
            ids=ni[k];beta=np.maximum(weights[k],0);beta/=beta.sum(1)[:,None];target_nodes=nodes[k];n=normal[k]
            body_v=np.sum(beta[:,:,None]*self.velocity[target_nodes],axis=1);relative=v[ids]-body_v;vn=np.sum(relative*n,axis=1)
            jn=m[ids]*np.maximum(-vn,0);tangent=relative-vn[:,None]*n;speed=np.linalg.norm(tangent,axis=1);needed=m[ids]*speed
            magnitude=np.where(needed<=self.mu_s*jn,needed,np.minimum(needed,self.mu_k*jn));direction=np.divide(tangent,speed[:,None],out=np.zeros_like(tangent),where=speed[:,None]>0)
            jt=-magnitude[:,None]*direction;j=jn[:,None]*n+jt;delta=np.maximum(-gap[k],0)[:,None]*n
            x[ids]+=delta;v[ids]+=j/m[ids,None];impulses[ids]+=j;correction[ids]+=delta
            np.add.at(reaction,target_nodes.reshape(-1),(-beta[:,:,None]*j[:,None,:]).reshape(-1,3))
            normal_loss+=float(.5*np.sum(jn*jn/m[ids]));friction_loss-=float(np.sum(jt*(tangent+.5*jt/m[ids,None])));work+=float(np.sum(j*body_v));count+=len(k)
        return {'positions_m':x,'velocities_m_s':v,'impulses_ns':impulses,'body_reaction_impulses_ns':reaction,
                'position_correction_m':correction,'normal_impact_dissipation_j':normal_loss,'friction_dissipation_j':friction_loss,
                'prescribed_surface_work_j':work,'contact_count':count,'unresolved_edge_contacts':unresolved,'maximum_candidate_pairs_in_batch':max_pairs_seen,
                'paired_impulse_residual_ns':impulses.sum(0)+reaction.sum(0),
                'paired_contact_point_angular_impulse_residual_nms':np.cross(x,impulses).sum(0)+np.cross(self.start+self.fraction*self.delta,reaction).sum(0),
                'coverage':'Discrete nearest face-interior nodal contact; recurrent contacts re-evaluated without persistent friction history; edge/vertex/CCD/self-contact not certified'}
