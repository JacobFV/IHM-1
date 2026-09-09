"""Read-only native balance diagnostics using source contact-sphere proxies.

Support polygons use loaded source foot contact centers, not measured foot meshes.
A positive static COM margin is not proof of dynamic balance or feasibility.
"""
import math
import numpy as np
from .moment_arm_control import native_center_of_mass

FEET={'calcn_r':'r','toes_r':'r','calcn_l':'l','toes_l':'l'}
def _finite(value):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value) or not math.isfinite(value):raise ValueError('Finite scalar diagnostic input required')
    return float(value)
def convex_hull_xz(points):
    """Counterclockwise unique monotone-chain hull; collinear endpoints retained."""
    pts=np.asarray(points,dtype=float)
    if pts.size==0:return []
    if pts.ndim!=2 or pts.shape[1]!=2 or not np.isfinite(pts).all():raise ValueError('Finite x/z point pairs required')
    values=sorted(set(map(tuple,pts.tolist())))
    if len(values)<=1:return [list(p) for p in values]
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower=[];upper=[]
    for p in values:
        while len(lower)>=2 and cross(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    for p in reversed(values):
        while len(upper)>=2 and cross(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    return [list(p) for p in lower[:-1]+upper[:-1]]
def signed_polygon_margin(point,hull):
    """Signed Euclidean distance to hull boundary, positive inside, None degenerate."""
    p=np.asarray(point,dtype=float);h=np.asarray(hull,dtype=float)
    if p.shape!=(2,) or not np.isfinite(p).all():raise ValueError('Finite x/z COM point required')
    if len(hull)<3:return None
    if h.ndim!=2 or h.shape[1]!=2 or not np.isfinite(h).all():raise ValueError('Finite planar hull required')
    area=float(np.sum(h[:,0]*np.roll(h[:,1],-1)-h[:,1]*np.roll(h[:,0],-1))*.5)
    if abs(area)<=1e-14:return None
    inside=True;distance=float('inf');orientation=1 if area>0 else -1
    for a,b in zip(h,np.roll(h,-1,axis=0)):
        edge=b-a;norm2=float(edge@edge)
        if norm2==0:continue
        rel=p-a;cross=float(edge[0]*rel[1]-edge[1]*rel[0])
        if orientation*cross < -1e-12:inside=False
        t=float(np.clip((rel@edge)/norm2,0,1));distance=min(distance,float(np.linalg.norm(p-(a+t*edge))))
    if distance<1e-12:return 0.
    return distance if inside else -distance

def observe_balance(state):
    """Return JSON-safe COM/support/actuation metrics without changing native state."""
    com,velocity=native_center_of_mass(state)
    points=[];loads={'r':0.,'l':0.};clearances=[];active=[]
    for contact in state.get('contacts',[]):
        point=np.asarray(contact['center_m'],dtype=float);force=np.asarray(contact['force_n'],dtype=float)
        if point.shape!=(3,) or force.shape!=(3,) or not np.isfinite(point).all() or not np.isfinite(force).all():raise ValueError('Invalid native contact observation')
        if state.get('environment')=='upright' and 'radius_m' in contact:
            radius=_finite(contact['radius_m'])
            if radius<0:raise ValueError('Nonnegative contact proxy radius required')
            clearances.append(float(point[1]-radius))
        side=FEET.get(contact.get('body_frame'))
        if side is not None:
            loads[side]+=max(0.,float(force[1]))
            if force[1]>1.:
                points.append([float(point[0]),float(point[2])]);active.append(contact.get('name'))
    hull=convex_hull_xz(points)
    muscles=state['muscles'];lags={};tendon={};low=high=0
    for name,m in muscles.items():
        activation=_finite(m['activation']);excitation=_finite(m['excitation'])
        if not 0<=activation<=1 or not 0<=excitation<=1:raise ValueError('Native activation/excitation outside [0,1]')
        lags[name]=excitation-activation;tendon[name]=_finite(m['tendon_force_n'])
        low+=excitation<=1e-8;high+=excitation>=1-1e-8
    max_name=max(lags,key=lambda n:abs(lags[n])) if lags else None
    return {'schema':'ihm.native-balance-observation.v1','time_s':_finite(state['time_s']),
        'com_position_ground_m':com.tolist(),'com_velocity_ground_m_s':velocity.tolist(),
        'support_hull_xz_m':hull,'com_support_margin_m':signed_polygon_margin(com[[0,2]],hull),
        'active_foot_contact_names':active,'active_foot_contact_count':len(active),
        'per_foot_normal_force_n':loads,'min_floor_clearance_m':min(clearances) if clearances else None,
        'muscle_activation_excitation_lag_max_abs':abs(lags[max_name]) if max_name is not None else 0.,
        'muscle_activation_excitation_lag_max_muscle':max_name,
        'muscle_activation_excitation_lag_at_max':lags[max_name] if max_name is not None else 0.,
        'tendon_force_max_n':max(tendon.values()) if tendon else None,
        'excitation_lower_bound_count':int(low),'excitation_upper_bound_count':int(high),'muscle_count':len(muscles),
        'scope':{'support':'Convex hull of calcn/toes source sphere centers projected to xz, only contacts with upward force >1N; not measured foot mesh geometry',
            'margin':'Signed static COM projection distance to polygon boundary; positive inside, None for fewer than three noncollinear support points; not dynamic stability or actuator-feasibility proof',
            'floor_clearance':'Minimum source contact-sphere center_y minus radius for upright ground y=0; includes contact penetration; not anatomical skin clearance; None for other environments',
            'excitation_bounds':'Counts commands at numeric bounds, not evidence that clipping occurred; lag is excitation minus activation'}}
