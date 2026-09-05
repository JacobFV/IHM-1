"""Delayed spinal stretch-feedback primitive on explicitly observed muscle paths.

The controller owns delay/activation state only. BodyMechanics owns position,
force and work. The fiber-length sensor is a declared stiff-tendon proxy, not a
native OpenSim fiber state or a calibrated spindle firing-rate model.
"""
from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import numpy as np

GEYER_HERR_URL='https://dam-prod.media.mit.edu/x/files/wp-content/uploads/sites/3/2013/04/A-Muscle-Reflex-Model-that-Encodes-Princlples-of-Legged-Mechanics-Produces-Human-Walking-Dynamics-and-Muscle-Activities.pdf'


def finite(value,lo,hi,name):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not lo<=value<=hi:
        raise ValueError(f'{name} must be finite in [{lo},{hi}]')
    return float(value)


@dataclass(frozen=True)
class ReflexParameters:
    length_gain: float=1.1
    length_offset: float=.71
    prestimulation: float=.01
    minimum_stimulation: float=.01
    loop_delay_s: float=.02
    activation_tau_s: float=.01

    def __post_init__(self):
        for key in ('length_gain','length_offset'):
            finite(getattr(self,key),0,100,key)
        for key in ('prestimulation','minimum_stimulation'):
            finite(getattr(self,key),0,1,key)
        finite(self.loop_delay_s,0,1,'loop_delay_s')
        finite(self.activation_tau_s,1e-6,1,'activation_tau_s')


class BodyReflex:
    """One TA-like length-feedback primitive, with a single lumped neural delay.

    At time t, consume a computed mechanical observation at t. Its sampled
    sensory feature arrives at t+loop_delay_s. Advance activation exactly across
    arrivals through t+dt. Caller applies returned activation on its next
    mechanical interval (one explicit exchange interval of additional latency).
    Descending inputs are external motoneuron drive/gain controls, not inferred
    from a cortical regional label. Nerve blocks clear in-flight sensory state.
    """
    def __init__(self,muscle,reference_centroids,parameters=None,nerve_id='explicit-test-nerve'):
        self.muscle=deepcopy(muscle)
        self.reference_centroids=deepcopy(reference_centroids)
        self.parameters=parameters or ReflexParameters()
        if not isinstance(self.parameters,ReflexParameters):raise ValueError('Expected ReflexParameters')
        if not isinstance(nerve_id,str) or not nerve_id:raise ValueError('Explicit nerve identity required')
        self.nerve_id=nerve_id
        for key in ('max_isometric_force_n','optimal_fiber_length_m','rest_path_length_m'):
            finite(muscle[key],1e-9,1e9,key)
        finite(muscle['pennation_angle_rad'],0,math.pi/2-1e-6,'pennation')
        if len(muscle['anchors'])<2:raise ValueError('A muscle requires a path')
        for anchor in muscle['anchors']:
            for value in (anchor['point_m'],reference_centroids[anchor['entity_id']]):
                a=np.asarray(value,float)
                if a.shape!=(3,) or not np.isfinite(a).all():raise ValueError('Invalid reference anchor')
        identity=dict(muscle=self.muscle,centroids=self.reference_centroids,parameters=asdict(self.parameters),nerve_id=nerve_id)
        self.model_sha256=hashlib.sha256(json.dumps(identity,sort_keys=True,allow_nan=False).encode()).hexdigest()
        self.time_s=0.;self.activation=0.;self.arrived_length=self.parameters.length_offset
        self.events=[];self.serial=0

    def observe(self,mechanical_state):
        if not isinstance(mechanical_state,dict):raise ValueError('Mechanical observation required')
        time=finite(mechanical_state.get('time_s'),0,1e9,'mechanical time')
        if abs(time-self.time_s)>1e-9:raise ValueError('Mechanics/controller clocks differ')
        points=[]
        for anchor in self.muscle['anchors']:
            key=anchor['entity_id']
            try:
                state=mechanical_state['entities'][key]
                x=np.asarray(state['translation_m'],float)
                r=np.asarray(state['rotation_matrix'],float)
                f=np.asarray(state['deformation_gradient'],float)
            except (KeyError,TypeError,ValueError):raise ValueError('Missing mechanical attachment state') from None
            if x.shape!=(3,) or r.shape!=(3,3) or f.shape!=(3,3) or not all(np.isfinite(a).all() for a in (x,r,f)):
                raise ValueError('Invalid attachment transform')
            if not np.allclose(r.T@r,np.eye(3),atol=1e-8) or np.linalg.det(r)<=0 or np.linalg.det(f)<=0:
                raise ValueError('Invalid mechanical rotation/deformation')
            c=np.asarray(self.reference_centroids[key])
            points.append(c+x+r@f@(np.asarray(anchor['point_m'])-c))
        length=sum(float(np.linalg.norm(b-a)) for a,b in zip(points,points[1:]))
        if length<=0:raise ValueError('Collapsed muscle path')
        # Calibrate zero path change to l_opt; tendon compliance is NOT solved here.
        fiber=self.muscle['optimal_fiber_length_m']+(length-self.muscle['rest_path_length_m'])/math.cos(self.muscle['pennation_angle_rad'])
        if not math.isfinite(fiber) or fiber<=0:raise ValueError('Fiber proxy left positive-length domain')
        try:force=mechanical_state['muscle_forces_n'][self.muscle['id']]
        except (KeyError,TypeError):raise ValueError('Missing source muscle force') from None
        finite(force,0,1e9,'muscle force')
        return dict(path_length_m=length,fiber_length_proxy_m=fiber,normalized_fiber_length=fiber/self.muscle['optimal_fiber_length_m'],
                    tendon_force_n=force,normalized_tendon_force=force/self.muscle['max_isometric_force_n'],
                    length_sensor_basis='reference-calibrated stiff-tendon path proxy; not native CE or measured spindle firing')

    def step(self,dt_s,mechanical_state,*,descending_drive=0.,descending_gain=1.,
             sensory_block=False,motor_block=False,controller_enabled=True):
        dt=finite(dt_s,1e-9,1,'dt_s')
        drive=finite(descending_drive,0,1,'descending_drive');gain=finite(descending_gain,0,2,'descending_gain')
        if any(type(x) is not bool for x in (sensory_block,motor_block,controller_enabled)):
            raise ValueError('Block/controller switches must be booleans')
        sensor=self.observe(mechanical_state)
        p=self.parameters;end=self.time_s+dt
        if end<=self.time_s:raise ValueError('Clock must advance')
        # All inputs validated above before queue, clock, or activation mutation.
        if sensory_block or motor_block:
            self.events=[];self.arrived_length=p.length_offset
        else:
            self.serial+=1
            self.events.append((self.time_s+p.loop_delay_s,self.serial,sensor['normalized_fiber_length']))
        def stimulation():
            reflex=gain*p.length_gain*(self.arrived_length-p.length_offset) if controller_enabled else 0.
            return 0. if motor_block else min(1.,max(p.minimum_stimulation,p.prestimulation+drive+reflex))
        cursor=self.time_s
        def advance(h):
            target=stimulation()
            self.activation=target+(self.activation-target)*math.exp(-h/p.activation_tau_s)
        while self.events and self.events[0][0]<=end+1e-12:
            arrival,_,length=self.events.pop(0)
            advance(max(0.,arrival-cursor));cursor=max(cursor,arrival)
            self.arrived_length=length
        advance(max(0.,end-cursor));self.time_s=end
        return dict(time_s=end,muscle_id=self.muscle['id'],nerve_id=self.nerve_id,sensor=sensor,
            delayed_normalized_length=self.arrived_length,stimulation=stimulation(),activation=self.activation,
            motor_activations={self.muscle['id']:self.activation},pending_samples=len(self.events),
            sensory_block=sensory_block,motor_block=motor_block,controller_enabled=controller_enabled,
            descending_drive=drive,descending_gain=gain,loop_delay_s=p.loop_delay_s,
            source='Geyer and Herr 2010 TA length-feedback primitive, transferred onto declared mechanical reduction',
            experimental_validation=False,cortical_controller=False)

    def checkpoint(self):
        return dict(schema='ihm.body-reflex.v1',model_sha256=self.model_sha256,time_s=self.time_s,
            activation=self.activation,arrived_length=self.arrived_length,events=[list(e) for e in self.events],serial=self.serial)

    def restore(self,checkpoint):
        c=deepcopy(checkpoint)
        if not isinstance(c,dict) or c.get('schema')!='ihm.body-reflex.v1' or c.get('model_sha256')!=self.model_sha256:
            raise ValueError('Checkpoint model identity mismatch')
        time=finite(c.get('time_s'),0,1e9,'checkpoint time');a=finite(c.get('activation'),0,1,'activation')
        length=finite(c.get('arrived_length'),0,1e6,'arrived length')
        serial=c.get('serial')
        if type(serial) is not int or serial<0 or not isinstance(c.get('events'),list):raise ValueError('Invalid event state')
        events=[];previous=(time,-1)
        for event in c['events']:
            if not isinstance(event,list) or len(event)!=3:raise ValueError('Invalid delayed sample')
            arrival=finite(event[0],time,1e9,'sample arrival');n=event[1]
            if type(n) is not int or not 0<n<=serial or (arrival,n)<=previous:raise ValueError('Invalid sample order')
            sample=finite(event[2],1e-12,1e6,'sample length');previous=(arrival,n);events.append((arrival,n,sample))
        self.time_s=time;self.activation=a;self.arrived_length=length;self.serial=serial;self.events=events
