"""Atomic first interior node-face impact for explicitly owned finite bodies.

One elastic/external kick, drift to event, existing inelastic Coulomb impulse,
then a validated ballistic remainder. No native skin masses are synthesized.
"""
import numpy as np
from .garment_swept_vertex_face import swept_vertex_face
from .contact_dynamics import resolve_node_triangle_contact


def _indices(values,count,name):
    ids=np.asarray(values)
    if ids.ndim!=1 or not len(ids) or ids.dtype.kind not in 'iu' or ids.min()<0 or ids.max()>=count or len(np.unique(ids))!=len(ids):raise ValueError('Unique valid '+name+' required')
    return ids


def _bounds_overlap(a,b,t,u):
    return not (np.any(np.maximum(a,b)<np.minimum(t.min(0),u.min(0))-1e-12) or np.any(np.minimum(a,b)>np.maximum(t.max(0),u.max(0))+1e-12))


def _remaining_contact(point,velocity,tri,velocities,duration):
    """Supported postimpact case: translating planar face, interior linear path."""
    if not np.allclose(velocities,velocities[0],rtol=0,atol=1e-12):raise ValueError('Deforming postimpact triangle requires contact schedule')
    normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);length=np.linalg.norm(normal)
    if length<=1e-14:raise ValueError('Degenerate postimpact face')
    normal/=length
    if (velocity-velocities[0])@normal < -1e-12:raise ValueError('Postimpact approach remains unresolved')
    relative_end=point+duration*(velocity-velocities[0])
    matrix=np.stack((tri[1]-tri[0],tri[2]-tri[0]),axis=1)
    # A line inside a convex triangle at both endpoints remains inside. Positive
    # separation is allowed, tangential drift through its edge is not certified.
    for p in (point,relative_end):
        beta=np.linalg.lstsq(matrix,p-tri[0],rcond=None)[0];beta=np.r_[1-beta.sum(),beta]
        if beta.min()<=1e-9:raise ValueError('Postimpact edge exit requires contact schedule')


def step_vertex_face_impact(bodies,dt_s,*,cloth_owner,surface_owner,node_ids,face_ids,
                           friction_static,friction_kinetic,max_energy_defect_j,external_forces_n=None,
                           gravity_m_s2=(0.,0.,0.),max_candidate_pairs=4096):
    """Commit at most one supported impact for all-mobile finite owner inputs.

    `bodies` expose Cloth-compatible position/velocity/mass, elastic_forces,
    max_explicit_dt_s and common time_s. Contact target exposes triangles.
    Caller supplies the exact candidate node/face inventory; no completeness or
    self-contact guarantee is inferred. All candidates are reconsidered after
    the impact. Additional impact/ambiguous geometry rejects the whole step.
    """
    if not bodies or len({id(b) for b in bodies.values()})!=len(bodies) or cloth_owner==surface_owner or cloth_owner not in bodies or surface_owner not in bodies:raise ValueError('Distinct finite owning bodies required')
    if isinstance(dt_s,bool) or not np.isfinite(dt_s) or dt_s<=0:raise ValueError('Positive finite interval required')
    if isinstance(max_energy_defect_j,bool) or not np.isfinite(max_energy_defect_j) or max_energy_defect_j<0:raise ValueError('Explicit finite nonnegative energy defect budget required')
    if isinstance(max_candidate_pairs,bool) or not isinstance(max_candidate_pairs,int) or max_candidate_pairs<1:raise ValueError('Positive finite candidate budget required')
    if not np.isfinite([friction_static,friction_kinetic]).all() or not 0<=friction_kinetic<=friction_static:raise ValueError('Finite ordered Coulomb coefficients required')
    owners=[o for b in bodies.values() for o in getattr(b,'material_owner_ids',())]
    if len(set(owners))!=len(owners):raise ValueError('Duplicate material owner identity')
    gravity=np.asarray(gravity_m_s2,float)
    if gravity.shape!=(3,) or not np.isfinite(gravity).all():raise ValueError('Finite gravity required')
    ext={} if external_forces_n is None else external_forces_n
    if set(ext)-set(bodies):raise ValueError('Unknown external force owner')
    states={};before=0.;external_work=0.;gravity_work=0.;external_impulse=np.zeros(3);external_angular=np.zeros(3)
    for name,body in bodies.items():
        if len(getattr(body,'fixed_nodes',())):raise ValueError('Prescribed support requires a separate work port')
        if not np.isfinite(body.time_s) or body.time_s<0 or not np.isfinite(body.max_explicit_dt_s) or dt_s>body.max_explicit_dt_s:raise ValueError('Invalid clock or elastic timestep')
        x=np.array(body.position_m,float,copy=True);v=np.array(body.velocity_m_s,float,copy=True);m=np.asarray(body.mass_kg,float)
        if x.ndim!=2 or x.shape[1]!=3 or v.shape!=x.shape or m.shape!=(len(x),) or not all(np.isfinite(z).all() for z in (x,v,m)) or np.any(m<=0):raise ValueError('Finite positive owner arrays required')
        force,energy=body.elastic_forces(x);force=np.asarray(force,float);applied=np.asarray(ext.get(name,np.zeros_like(x)),float)
        if force.shape!=x.shape or applied.shape!=x.shape or not np.isfinite(force).all() or not np.isfinite(applied).all() or not np.isfinite(energy):raise ValueError('Invalid elastic or applied force')
        kick=dt_s*(force+applied+m[:,None]*gravity);newv=v+kick/m[:,None]
        # Work during the discrete kick uses midpoint kick velocity, not the
        # end-of-step contact-altered displacement. Elastic defect is separate.
        average=.5*(v+newv);external_work+=float(dt_s*np.sum(applied*average));gravity_work+=float(dt_s*np.sum(m[:,None]*gravity*average))
        imposed=dt_s*(applied+m[:,None]*gravity);external_impulse+=imposed.sum(0);external_angular+=np.cross(x,imposed).sum(0)
        before+=float(energy+.5*np.sum(m[:,None]*v*v))
        states[name]=dict(x=x,v=newv,m=m,oldx=x.copy(),oldv=v,time=body.time_s)
    times=[s['time'] for s in states.values()]
    if max(times)-min(times)>1e-12:raise ValueError('Owner clocks differ')
    a=states[cloth_owner];b=states[surface_owner];tri=np.asarray(bodies[surface_owner].triangles)
    if tri.ndim!=2 or tri.shape[1]!=3 or tri.dtype.kind not in 'iu' or not len(tri) or tri.min()<0 or tri.max()>=len(b['x']):raise ValueError('Valid owning surface faces required')
    nodes=_indices(node_ids,len(a['x']),'contact nodes');faces=_indices(face_ids,len(tri),'contact faces')
    if len(nodes)*len(faces)>max_candidate_pairs:raise ValueError('Contact candidate budget exceeded')
    def search(duration,previous=None):
        events=[]
        for node in nodes:
            for face in faces:
                ids=tri[face];ta=b['x'][ids];tb=ta+duration*b['v'][ids];p=a['x'][node];q=p+duration*a['v'][node]
                if previous==(int(node),int(face)):
                    _remaining_contact(p,a['v'][node],ta,b['v'][ids],duration);continue
                if not _bounds_overlap(p,q,ta,tb):continue
                event=swept_vertex_face(p,q,ta,tb)
                if event['status']=='unresolved':raise ValueError('Unresolved contact: '+event['reason'])
                if event['status']=='crossing':
                    if event['direction']!='inward':raise ValueError('Outward crossing indicates unresolved initial containment')
                    events.append((event['fraction']*duration,int(node),int(face),event))
        return sorted(events,key=lambda r:r[0])
    events=search(dt_s);record=None;loss=0.;normal_loss=0.;pair_residual=np.zeros(3)
    if len(events)>1:raise ValueError('Multiple candidate impacts require contact schedule')
    if events:
        at,node,face,event=events[0]
        if at<=1e-12:raise ValueError('Initial impact requires persistent contact history')
        for s in states.values():s['x']+=at*s['v']
        ids=tri[face]
        response=resolve_node_triangle_contact(a['x'][[node]],a['v'][[node]],a['m'][[node]],b['x'][ids],b['v'][ids],b['m'][ids],np.array([[0,1,2]]),node_ids=np.array([0]),friction_static=friction_static,friction_kinetic=friction_kinetic,search_distance_m=1e-9)
        if response['contact_count']!=1 or response['unresolved_edge_contacts'] or max(np.linalg.norm(response['position_correction_a_m']),np.linalg.norm(response['position_correction_b_m']))>1e-11:raise ValueError('Impact root did not produce coincident interior contact')
        # Preserve event geometry: tiny projection artifacts in the discrete
        # kernel are not a separate correction or source of impulse work.
        a['v'][node]=response['velocity_a_m_s'][0];b['v'][ids]=response['velocity_b_m_s']
        loss=response['dissipation_j'];pair_residual=response['paired_impulse_residual_ns']
        beta=np.asarray(event['barycentric']);jn=float(response['impulse_a_ns'][0]@np.asarray(event['normal']))
        normal_loss=.5*jn*jn*(1/a['m'][node]+np.sum(beta*beta/b['m'][ids]))
        remaining=dt_s-at
        if remaining>0 and search(remaining,(node,face)):raise ValueError('Additional postimpact event requires contact schedule')
        for s in states.values():s['x']+=remaining*s['v']
        record=dict(time_s=times[0]+at,elapsed_time_s=at,node=node,face=face,**event,normal_restitution=0.,impulse_cloth_ns=response['impulse_a_ns'][0].tolist(),impulse_surface_nodes_ns=response['impulse_b_ns'].tolist())
    else:
        for s in states.values():s['x']+=dt_s*s['v']
    after=0.;momentum=np.zeros(3);angular=np.zeros(3)
    for name,s in states.items():
        force,energy=bodies[name].elastic_forces(s['x'])
        if not all(np.isfinite(z).all() for z in (s['x'],s['v'],force)) or not np.isfinite(energy):raise ValueError('Invalid proposed finite state')
        after+=float(energy+.5*np.sum(s['m'][:,None]*s['v']**2))
        momentum+=np.sum(s['m'][:,None]*(s['v']-s['oldv']),axis=0)
        angular+=(np.cross(s['x'],s['m'][:,None]*s['v'])-np.cross(s['oldx'],s['m'][:,None]*s['oldv'])).sum(0)
    defect=after-before-external_work-gravity_work+loss
    if not np.isfinite(defect) or abs(defect)>max_energy_defect_j:raise ValueError('Numerical energy defect exceeds explicit acceptance budget')
    momentum_error=momentum-external_impulse;angular_error=angular-external_angular
    if np.linalg.norm(momentum_error)>1e-10 or np.linalg.norm(angular_error)>1e-10:raise ValueError('Momentum ledger exceeds absolute numerical tolerance')
    for name,s in states.items():bodies[name].position_m=s['x'];bodies[name].velocity_m_s=s['v'];bodies[name].time_s=times[0]+dt_s
    return dict(time_s=times[0]+dt_s,event=record,normal_restitution=0.,contact_dissipation_j=float(loss),normal_impact_dissipation_j=float(normal_loss),friction_dissipation_j=float(loss-normal_loss),
        max_energy_defect_j=max_energy_defect_j,momentum_tolerance_ns=1e-10,angular_tolerance_nms=1e-10,
        paired_contact_impulse_ns=pair_residual,momentum_residual_ns=momentum-external_impulse,
        angular_impulse_residual_nms=angular-external_angular,numerical_energy_defect_j=float(defect),
        external_work_j=external_work,gravity_work_j=gravity_work,
        scope='One supported finite-owner interior event plus validated linear remainder; no self/edge contact or native garment coupling')
