"""Opt-in floating-point swept edge impact, with explicit unresolved geometry.

No contact thickness, self-contact schedule or containment guarantee. Positions
at a verified event share a contact point; finite-mass impulses preserve force
and angular ledgers without an invented finite-clearance friction couple.
"""
import numpy as np

def _edge(x,name):
    a=np.asarray(x,float)
    if a.shape!=(2,3) or not np.isfinite(a).all():raise ValueError('Finite2x3 '+name+' required')
    return a

def _positive(x,name):
    if isinstance(x,(bool,np.bool_)) or not np.isscalar(x) or not np.isfinite(x) or x<=0:raise ValueError('Positive finite '+name+' required')
    return float(x)

def swept_edge_event(position_a_m,velocity_a_m_s,position_b_m,velocity_b_m_s,duration_s):
    """First nondegenerate linear-trajectory edge crossing during (0,duration].

    Coplanarity is a cubic in normalized time; roots are checked against actual
    segment parameters and point coincidence. This is floating-point CCD, not
    certified exact predicates. Persistent/initial contact is a separate problem.
    """
    a=_edge(position_a_m,'positions');b=_edge(position_b_m,'positions')
    va=_edge(velocity_a_m_s,'velocities');vb=_edge(velocity_b_m_s,'velocities');dt=_positive(duration_s,'duration')
    ea=a[1]-a[0];eb=b[1]-b[0];da=(va[1]-va[0])*dt;db=(vb[1]-vb[0])*dt
    q=b[0]-a[0];dq=(vb[0]-va[0])*dt
    scale=max(np.linalg.norm(ea),np.linalg.norm(eb),np.linalg.norm(da),np.linalg.norm(db),np.linalg.norm(q),np.linalg.norm(dq),1e-12)
    tolerance=256*np.finfo(float).eps*scale
    for edge,change in ((ea,da),(eb,db)):
        at=float(np.clip(-edge@change/(change@change),0,1)) if change@change else 0.
        if np.linalg.norm(edge+at*change)<=tolerance:return {'status':'unresolved','reason':'collapsed_edge_during_interval'}
    aa=np.r_[a,a+dt*va];bb=np.r_[b,b+dt*vb]
    if np.any(aa.max(axis=0)<bb.min(axis=0)-tolerance) or np.any(bb.max(axis=0)<aa.min(axis=0)-tolerance):return {'status':'no_impact'}
    cross0=np.cross(ea,eb);cross1=np.cross(da,eb)+np.cross(ea,db);cross2=np.cross(da,db)
    coefficients=np.array([q@cross0,dq@cross0+q@cross1,dq@cross1+q@cross2,dq@cross2])
    magnitude=float(np.max(np.abs(coefficients)));precision=256*np.finfo(float).eps*scale**3
    if magnitude<=precision:return {'status':'unresolved','reason':'coplanar_or_unresolved_polynomial','coplanarity_coefficients':coefficients.tolist()}
    degree=3
    while degree and abs(coefficients[degree])<=precision:degree-=1
    roots=np.roots((coefficients[:degree+1]/magnitude)[::-1]) if degree else np.array([])
    candidates=sorted(float(r.real) for r in roots if abs(r.imag)<1e-9 and -1e-10<=r.real<=1+1e-10)
    for u in candidates:
        u=float(np.clip(u,0,1));pa=a+u*dt*va;pb=b+u*dt*vb
        e=pa[1]-pa[0];f=pb[1]-pb[0];normal=np.cross(e,f);normal_length=np.linalg.norm(normal)
        if normal_length<=256*np.finfo(float).eps*np.linalg.norm(e)*np.linalg.norm(f):
            return {'status':'unresolved','reason':'parallel_or_collinear_event','time_s':u*dt}
        normal/=normal_length
        matrix=np.stack((e,-f),axis=1);parameters=np.linalg.lstsq(matrix,pb[0]-pa[0],rcond=None)[0]
        if np.any(parameters<-1e-10) or np.any(parameters>1+1e-10):continue
        sa,sb=np.clip(parameters,0,1);wa=np.array([1-sa,sa]);wb=np.array([1-sb,sb]);qa=wa@pa;qb=wb@pb
        separation=float(np.linalg.norm(qa-qb))
        if separation>max(1e-12,tolerance*8):return {'status':'unresolved','reason':'root_did_not_produce_coincident_points','time_s':u*dt,'point_separation_m':separation}
        if u<=1e-10:return {'status':'unresolved','reason':'initial_contact_requires_history','time_s':0.}
        relative=wa@va-wb@vb;vn=float(relative@normal)
        if abs(vn)<=256*np.finfo(float).eps*max(np.linalg.norm(relative),1e-12):
            return {'status':'unresolved','reason':'grazing_event','time_s':u*dt}
        if vn>0:normal=-normal
        return {'status':'impact','time_s':u*dt,'weights_a':wa.tolist(),'weights_b':wb.tolist(),
                'normal':normal.tolist(),'contact_point_m':((qa+qb)/2).tolist(),'point_separation_m':separation,
                'coplanarity_coefficients':coefficients.tolist(),
                'method':'normalized-time cubic coplanarity, segment-parameter and coincidence checks; floating-point predicates'}
    return {'status':'no_impact'}

def resolve_edge_impact(position_a_m,velocity_a_m_s,mass_a_kg,position_b_m,velocity_b_m_s,mass_b_kg,event,
                        *,friction_static,friction_kinetic,mobile_a=None,mobile_b=None):
    """Pure zero-restitution impulse at an already-verified common event point."""
    a=_edge(position_a_m,'positions');b=_edge(position_b_m,'positions');va=_edge(velocity_a_m_s,'velocities');vb=_edge(velocity_b_m_s,'velocities')
    ma=np.asarray(mass_a_kg,float);mb=np.asarray(mass_b_kg,float)
    if ma.shape!=(2,) or mb.shape!=(2,) or not np.isfinite(ma).all() or not np.isfinite(mb).all() or np.any(ma<=0) or np.any(mb<=0):raise ValueError('Positive endpoint masses required')
    if event.get('status')!='impact':raise ValueError('A verified impact event is required')
    alpha=np.asarray(event['weights_a'],float);beta=np.asarray(event['weights_b'],float);normal=np.asarray(event['normal'],float)
    if any(w.shape!=(2,) or not np.isfinite(w).all() or np.any(w<0) or not np.isclose(w.sum(),1,atol=1e-12,rtol=0) for w in (alpha,beta)):
        raise ValueError('Convex endpoint weights required')
    if normal.shape!=(3,) or not np.isfinite(normal).all() or not np.isclose(normal@normal,1,rtol=0,atol=1e-10):raise ValueError('Unit contact normal required')
    point=np.asarray(event['contact_point_m'],float)
    if point.shape!=(3,) or not np.isfinite(point).all() or max(np.linalg.norm(alpha@a-point),np.linalg.norm(beta@b-point))>1e-10:
        raise ValueError('Impact positions do not share the verified contact point')
    for coefficient in (friction_static,friction_kinetic):
        if isinstance(coefficient,(bool,np.bool_)) or not np.isscalar(coefficient) or not np.isfinite(coefficient) or coefficient<0:raise ValueError('Finite nonnegative friction required')
    if friction_kinetic>friction_static:raise ValueError('Kinetic friction exceeds static friction')
    movable=[]
    for values in (mobile_a,mobile_b):
        raw=np.ones(2,dtype=bool) if values is None else np.asarray(values)
        if raw.shape!=(2,) or raw.dtype.kind!='b':raise ValueError('Boolean endpoint mobility masks required')
        movable.append(raw)
    am,bm=movable;ia=am/ma;ib=bm/mb;w=float(alpha**2@ia+beta**2@ib)
    if w<=0:raise ValueError('Impact has no mobile endpoint')
    relative=alpha@va-beta@vb;vn=float(relative@normal)
    if vn>1e-12:raise ValueError('Impact normal is separating, not approaching')
    jn=max(-vn,0)/w;tangent=relative-vn*normal;speed=float(np.linalg.norm(tangent));needed=speed/w
    stick=needed<=friction_static*jn
    magnitude=needed if stick else min(needed,friction_kinetic*jn)
    jt=-magnitude*tangent/speed if speed>0 else np.zeros(3)
    impulse=jn*normal+jt;ja=alpha[:,None]*impulse;jb=-beta[:,None]*impulse
    new_a=va+ia[:,None]*ja;new_b=vb+ib[:,None]*jb
    loss=.5*jn**2*w-float(jt@(tangent+.5*w*jt))
    before=float(.5*np.sum(ma[:,None]*va**2)+.5*np.sum(mb[:,None]*vb**2))
    after=float(.5*np.sum(ma[:,None]*new_a**2)+.5*np.sum(mb[:,None]*new_b**2))
    support=-(ja[~am].sum(axis=0)+jb[~bm].sum(axis=0))
    support_work=-float(np.sum(ja[~am]*va[~am])+np.sum(jb[~bm]*vb[~bm]))
    support_angular=-(np.cross(a[~am],ja[~am]).sum(axis=0)+np.cross(b[~bm],jb[~bm]).sum(axis=0))
    delta_a=ma[:,None]*(new_a-va);delta_b=mb[:,None]*(new_b-vb)
    return {'velocity_a_m_s':new_a,'velocity_b_m_s':new_b,'impulse_a_ns':ja,'impulse_b_ns':jb,
            'normal_impulse_ns':jn,'tangential_impulse_ns':jt,'friction_regime':'sticking' if stick else 'sliding',
            'effective_inverse_mass_kg_inverse':w,'dissipation_j':float(loss),
            'kinetic_energy_before_j':before,'kinetic_energy_after_j':after,'support_work_j':support_work,
            'support_impulse_ns':support,'support_angular_impulse_nms':support_angular,
            'paired_impulse_residual_ns':ja.sum(axis=0)+jb.sum(axis=0),
            'momentum_residual_ns':delta_a.sum(axis=0)+delta_b.sum(axis=0)-support,
            'angular_impulse_residual_nms':np.cross(a,delta_a).sum(axis=0)+np.cross(b,delta_b).sum(axis=0)-support_angular,
            'energy_residual_j':after-before-support_work+loss}
