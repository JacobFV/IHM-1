"""Transactional opt-in elastic-body step containing one swept edge impact.

This first-event fixture path is not a sustained-contact or self-contact solver.
The post-event remainder is ballistic and its integration defect is reported.
"""
from dataclasses import dataclass
import numpy as np
from .swept_edge_contact import swept_edge_event,resolve_edge_impact

@dataclass(frozen=True)
class SweptEdgePair:
    a_id:str
    edge_a:tuple
    b_id:str
    edge_b:tuple
    friction_static:float=0.
    friction_kinetic:float=0.

class UnresolvedSweptContact(ValueError):
    def __init__(self,records):
        self.records=records
        super().__init__('Unsupported swept contact; no body advanced: '+str(records))

def step_first_impact(bodies,dt_s,*,contacts=(),external_forces_n=None,gravity_m_s2=(0.,0.,0.)):
    """Kick/drift with earliest single impact and atomic body-state commit.

    Fixed nodes retain their given velocities, including prescribed moving
    support. No implicit contact state is carried into subsequent calls.
    """
    if not bodies or isinstance(dt_s,bool) or not np.isfinite(dt_s) or dt_s<=0:raise ValueError('Bodies and positive finite step required')
    gravity=np.asarray(gravity_m_s2,float)
    if gravity.shape!=(3,) or not np.isfinite(gravity).all():raise ValueError('Finite gravity required')
    if len({id(b) for b in bodies.values()})!=len(bodies):raise ValueError('Duplicate body instance')
    owners=[o for b in bodies.values() for o in getattr(b,'material_owner_ids',())]
    if len(owners)!=len(set(owners)):raise ValueError('Duplicate material owner')
    times=[b.time_s for b in bodies.values()]
    if not np.isfinite(times).all() or min(times)<0 or max(times)-min(times)>1e-12:raise ValueError('Body clocks disagree')
    if any(dt_s>b.max_explicit_dt_s for b in bodies.values()):raise ValueError('Step exceeds elastic restriction')
    ext={} if external_forces_n is None else external_forces_n
    if set(ext)-set(bodies):raise ValueError('Unknown external force owner')
    states={};before=0.;support=np.zeros(3);support_angular=np.zeros(3);support_work=0.;external_impulse=np.zeros(3);external_angular=np.zeros(3)
    for name,b in bodies.items():
        x=np.array(b.position_m,float,copy=True);v=np.array(b.velocity_m_s,float,copy=True);m=np.asarray(b.mass_kg,float)
        if x.ndim!=2 or x.shape[1]!=3 or v.shape!=x.shape or m.shape!=(len(x),) or not np.isfinite(x).all() or not np.isfinite(v).all() or not np.isfinite(m).all() or np.any(m<=0):raise ValueError('Invalid finite-mass body state')
        force,energy=b.elastic_forces(x);force=np.asarray(force,float);external=np.asarray(ext.get(name,np.zeros_like(x)),float)
        if force.shape!=x.shape or not np.isfinite(force).all() or not np.isscalar(energy) or not np.isfinite(energy):raise ValueError('Invalid elastic force or energy')
        if external.shape!=x.shape or not np.isfinite(external).all():raise ValueError('Invalid external force')
        fixed=np.asarray(getattr(b,'fixed_nodes',()));mobile=np.ones(len(x),bool)
        if fixed.ndim!=1 or (fixed.size and fixed.dtype.kind not in 'iu') or np.any(fixed<0) or np.any(fixed>=len(x)):raise ValueError('Invalid fixed support indices')
        fixed=fixed.astype(np.int64)
        mobile[fixed]=False;total=force+external+m[:,None]*gravity
        kick=dt_s*total;reaction=-kick[fixed]
        support+=reaction.sum(axis=0);support_angular+=np.cross(x[fixed],reaction).sum(axis=0);support_work+=float(np.sum(reaction*v[fixed]))
        trial_v=v+kick/m[:,None];trial_v[fixed]=v[fixed]
        states[name]={'oldx':x.copy(),'oldv':v,'x':x,'v':trial_v,'mass':m,'mobile':mobile,'external':external}
        before+=energy+.5*np.sum(m[:,None]*v*v)-np.sum(m[:,None]*x*gravity)
        impulse=dt_s*(external+m[:,None]*gravity);external_impulse+=impulse.sum(axis=0);external_angular+=np.cross(x,impulse).sum(axis=0)
    events=[];unsupported=[]
    for index,pair in enumerate(contacts):
        if pair.a_id==pair.b_id or pair.a_id not in states or pair.b_id not in states:raise ValueError('Two distinct known bodies required')
        sa,sb=states[pair.a_id],states[pair.b_id];ia=np.asarray(pair.edge_a);ib=np.asarray(pair.edge_b)
        for ids,s in ((ia,sa),(ib,sb)):
            if ids.shape!=(2,) or ids.dtype.kind not in 'iu' or ids[0]==ids[1] or np.any(ids<0) or np.any(ids>=len(s['x'])):raise ValueError('Two distinct valid endpoint indices required')
        e=swept_edge_event(sa['x'][ia],sa['v'][ia],sb['x'][ib],sb['v'][ib],dt_s)
        if e['status']=='unresolved':unsupported.append({'pair_index':index,**e})
        elif e['status']=='impact':events.append((e['time_s'],index,pair,ia,ib,e))
    if unsupported:raise UnresolvedSweptContact(unsupported)
    if len(events)>1:raise UnresolvedSweptContact([{'reason':'multiple_events_require_contact_schedule','event_times_s':[e[0] for e in events]}])
    impact=None;event_record=None;loss=0.
    if events:
        at,index,pair,ia,ib,event=events[0]
        for s in states.values():s['x']+=at*s['v']
        a,b=states[pair.a_id],states[pair.b_id]
        impact=resolve_edge_impact(a['x'][ia],a['v'][ia],a['mass'][ia],b['x'][ib],b['v'][ib],b['mass'][ib],event,
                                  mobile_a=a['mobile'][ia],mobile_b=b['mobile'][ib],friction_static=pair.friction_static,friction_kinetic=pair.friction_kinetic)
        a['v'][ia]=impact['velocity_a_m_s'];b['v'][ib]=impact['velocity_b_m_s']
        support+=impact['support_impulse_ns'];support_angular+=impact['support_angular_impulse_nms'];support_work+=impact['support_work_j'];loss=impact['dissipation_j']
        for s in states.values():s['x']+=(dt_s-at)*s['v']
        event_record={'pair_index':index,'body_time_s':times[0]+at,**event}
    else:
        for s in states.values():s['x']+=dt_s*s['v']
    after=0.;work=0.;momentum=np.zeros(3);angular=np.zeros(3);elastic_total=0.;kinetic_total=0.
    for name,b in bodies.items():
        s=states[name];x,v,m=s['x'],s['v'],s['mass'];force,energy=b.elastic_forces(x)
        if not np.isfinite(x).all() or not np.isfinite(v).all() or np.asarray(force).shape!=x.shape or not np.isfinite(force).all() or not np.isscalar(energy) or not np.isfinite(energy):raise ValueError('Invalid proposed elastic state')
        kinetic=float(.5*np.sum(m[:,None]*v*v));elastic_total+=energy;kinetic_total+=kinetic
        after+=energy+kinetic-np.sum(m[:,None]*x*gravity);work+=np.sum(s['external']*(x-s['oldx']))
        momentum+=np.sum(m[:,None]*(v-s['oldv']),axis=0)
        angular+=(np.cross(x,m[:,None]*v)-np.cross(s['oldx'],m[:,None]*s['oldv'])).sum(axis=0)
    # Commit only after every proposed geometry and every contact was validated.
    for name,b in bodies.items():b.position_m=states[name]['x'];b.velocity_m_s=states[name]['v'];b.time_s=times[0]+dt_s
    return {'time_s':times[0]+dt_s,'event':event_record,'impact':impact,'elastic_energy_j':float(elastic_total),
            'kinetic_energy_j':kinetic_total,'total_energy_j':float(after),'external_work_j':float(work),'support_work_j':support_work,
            'support_impulse_ns':support,'support_angular_impulse_nms':support_angular,'contact_dissipation_j':float(loss),
            'momentum_residual_ns':momentum-external_impulse-support,'angular_impulse_residual_nms':angular-external_angular-support_angular,
            'numerical_energy_defect_j':float(after-before-work-support_work+loss),
            'post_event_model':'ballistic remainder of this single-impact step; sustained contact and additional post-impact collision events are not certified'}
