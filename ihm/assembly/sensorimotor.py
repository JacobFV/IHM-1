"""Source ankle reflexes plus an explicitly engineered pinned-IBM descending bridge.

Excitation is dimensionless. Native mechanics owns activation, fiber dynamics,
force and work. No cortical label is treated as an identified motor policy.
"""
from ihm.brain.active_source import DEFAULT_SOURCE,resolve_source
from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import heapq
import json
from pathlib import Path
import numpy as np
from .brain import BodyBrain
from .reflexes import finite
from .sensorimotor_catalog import native_muscle_catalog

MUSCLES=tuple(f'{name}_{side}' for side in ('r','l') for name in ('tibant','soleus','gasmed','gaslat'))
SOURCE_URL='https://www.cs.cmu.edu/~cga/tmp-public/song.pdf'
SOURCE_SHA256='8a60efec9ceba5120c1679610f676ab297e8dce887ed4ae6f0890f7e17e13d94'

@dataclass(frozen=True)
class SensorimotorParameters:
    afferent_delay_s: float=.01
    efferent_delay_s: float=.01
    stance_threshold_n: float=5.
    cortical_afferent_hz_per_strain: float=200.
    cortical_afferent_hz_per_normalized_force: float=40.
    cortical_command_hz: float=100.
    cortical_gain_per_hz: float=.005
    cortical_drive_per_hz: float=.002

    def __post_init__(self):
        for key,value in asdict(self).items():finite(value,0,1000,key)
        for key in ('afferent_delay_s','efferent_delay_s'):
            finite(getattr(self,key),0,1,key)

class SensorimotorController:
    """Source catalog effectors, eight ankle reflexes, transactional neural state.

    At t observe the plant, deliver sensory samples after afferent delay, execute
    the preserved IBM ODE over this exchange interval, derive spinal excitations
    and queue them after efferent delay. Outputs at t+dt apply on the NEXT plant
    interval. The exchange discretization adds latency; it is reported, not
    subtracted from measured time. Blocks discard in-flight signals immediately.
    ``descending`` contains per-effector requested drive fractions: these enter
    the corresponding precentral population and gate an explicit rate decoder.
    Additional sensory inputs are already-delayed external receptor endpoints
    from the previous exchange; they share this brain integration and saturation.
    """
    @classmethod
    def from_root(cls,root,*,source_pin=DEFAULT_SOURCE,**kwargs):
        root=Path(root)
        source_pin=resolve_source(root,source_pin)
        paper=root/'data/raw/sensorimotor/geyer_herr_2010.pdf'
        if hashlib.sha256(paper.read_bytes()).hexdigest()!=SOURCE_SHA256:
            raise ValueError('Reflex primary source hash mismatch')
        kwargs.setdefault('muscle_catalog',native_muscle_catalog(root))
        return cls(BodyBrain(json.loads((root/'data/derived/canonical/brain.json').read_text()),root=root,source_pin=source_pin),**kwargs)

    def __init__(self,brain,parameters=None,*,muscle_catalog=None):
        if not isinstance(brain,BodyBrain) or brain.time_s!=0:raise ValueError('Fresh pinned BodyBrain required')
        self.brain=brain;self.parameters=parameters or SensorimotorParameters()
        if not isinstance(self.parameters,SensorimotorParameters):raise ValueError('Invalid controller parameters')
        for side in ('lh','rh'):
            for region in ('postcentral','precentral'):
                if f'brain-{side}-{region}' not in brain.ids:raise ValueError('Required sensorimotor population missing')
        self.catalog=deepcopy(muscle_catalog or [{'id':k,'side':k[-1],'body_group':'ankle_foot',
            'motor_region':'brain-'+('lh' if k[-1]=='r' else 'rh')+'-precentral',
            'sensory_region':'brain-'+('lh' if k[-1]=='r' else 'rh')+'-postcentral',
            'assignment_basis':'explicit ankle-only engineering fixture'} for k in MUSCLES])
        self.muscles=tuple(row['id'] for row in self.catalog)
        if len(set(self.muscles))!=len(self.muscles) or not set(MUSCLES)<=set(self.muscles):raise ValueError('Catalog must contain unique ankle source effectors')
        self.bindings={row['id']:row for row in self.catalog}
        for row in self.catalog:
            if row['side'] not in ('r','l') or row['motor_region'] not in brain.ids or row['sensory_region'] not in brain.ids:raise ValueError('Invalid sensorimotor assignment')
        self.reference_rates=dict(zip(brain.ids,brain.state[:,1].tolist()))
        identity={'brain':brain.data,'brain_source_identity':brain.source_identity,'brain_step_s':brain.max_step_s,'parameters':asdict(self.parameters),
                  'catalog':self.catalog,'reflex_source_sha256':SOURCE_SHA256,'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'brain_implementation_sha256':hashlib.sha256(Path(__file__).with_name('brain.py').read_bytes()).hexdigest()}
        self.model_sha256=hashlib.sha256(json.dumps(identity,sort_keys=True,allow_nan=False).encode()).hexdigest()
        self.time_s=0.;self.events=[];self.serial=0
        self.arrived={};self.excitations={k:0. for k in self.muscles}

    def _blocks(self,value,label):
        if not isinstance(value,(list,tuple,set)) or any(not isinstance(k,str) for k in value) or set(value)-set(self.muscles):
            raise ValueError(f'Unknown {label} muscle')
        return set(value)

    def _validate_observation(self,observation):
        if not isinstance(observation,dict):raise ValueError('Mechanical observation required')
        time=finite(observation.get('time_s'),0,1e9,'mechanical time')
        if abs(time-self.time_s)>1e-8:raise ValueError('Mechanics and neural clocks differ')
        muscles=observation.get('muscles');contact=observation.get('foot_contact_force_n')
        if not isinstance(muscles,dict) or not isinstance(contact,dict):raise ValueError('Muscles and actual foot contact required')
        loads={s:finite(contact.get(s),0,1e9,'foot normal force') for s in ('r','l')}
        sensors={}
        for key in self.muscles:
            m=muscles.get(key)
            if not isinstance(m,dict) or not isinstance(m.get('sensor_basis'),str) or not m['sensor_basis']:
                raise ValueError(f'Missing sensor provenance: {key}')
            opt=finite(m.get('optimal_fiber_length_m'),1e-9,100,'optimal fiber length')
            fmax=finite(m.get('max_isometric_force_n'),1e-9,1e9,'maximum force')
            # A proxy must remain explicitly named. Never reconstruct CE from path silently.
            field='fiber_length_m' if 'fiber_length_m' in m else 'fiber_length_proxy_m'
            length=finite(m.get(field),1e-9,100,'fiber length')
            force=finite(m.get('tendon_force_n'),0,1e9,'tendon force')
            sensors[key]={'length':length/opt,'force':force/fmax,'force_n':force,'fmax_n':fmax,
                          'stance':loads[self.bindings[key]['side']]>self.parameters.stance_threshold_n,
                          'sensor_basis':m['sensor_basis'],'length_field':field}
        return sensors

    def _queue(self,arrival,kind,key,value):
        self.serial+=1;heapq.heappush(self.events,(arrival,self.serial,kind,key,deepcopy(value)))

    def step(self,dt_s,mechanical_observation,*,descending=None,sensory_blocks=(),motor_blocks=(),physiology=None,additional_sensory_inputs_hz=None):
        dt=finite(dt_s,1e-9,.1,'dt_s');sensors=self._validate_observation(mechanical_observation)
        sb=self._blocks(sensory_blocks,'sensory block');mb=self._blocks(motor_blocks,'motor block')
        descending={} if descending is None else descending
        if not isinstance(descending,dict) or set(descending)-set(self.muscles):raise ValueError('Unknown descending effector')
        descending={k:finite(v,0,1,'descending drive fraction') for k,v in descending.items()}
        additional={} if additional_sensory_inputs_hz is None else additional_sensory_inputs_hz
        if not isinstance(additional,dict) or set(additional)-set(self.brain.ids):raise ValueError('Unknown additional sensory population')
        additional={k:finite(v,0,1000,'additional sensory rate Hz') for k,v in additional.items()}
        saved=self.checkpoint()
        try:return self._step(dt,sensors,descending,sb,mb,physiology,additional)
        except Exception:
            self.restore(saved)
            raise

    def _step(self,dt,sensors,descending,sb,mb,physiology,additional):
        p=self.parameters;end=self.time_s+dt
        if end<=self.time_s:raise ValueError('Clock must advance')
        self.events=[e for e in self.events if not(e[2]=='sensor' and e[3] in sb or e[2]=='motor' and e[3] in mb)]
        heapq.heapify(self.events)
        for k in sb:self.arrived.pop(k,None)
        for k in mb:self.excitations[k]=0.
        for k,v in sensors.items():
            if k not in sb:self._queue(self.time_s+p.afferent_delay_s,'sensor',k,v)
        # Brain integrates the previously arrived afferents (explicit causal exchange).
        inputs={}
        for key,s in self.arrived.items():
            region=self.bindings[key]['sensory_region']
            rate=p.cortical_afferent_hz_per_strain*max(0,s['length']-1)+p.cortical_afferent_hz_per_normalized_force*s['force']
            inputs[region]=min(1000.,inputs.get(region,0)+rate)
        for key,value in descending.items():
            region=self.bindings[key]['motor_region']
            inputs[region]=min(1000.,inputs.get(region,0)+p.cortical_command_hz*value)
        # Caller supplies already-delayed receptor output from the previous exchange.
        # Pool here so the existing shared brain advances exactly once.
        for region,rate in additional.items():
            inputs[region]=min(1000.,inputs.get(region,0)+rate)
        neural=self.brain.step(dt,physiology=physiology,sensory_inputs_hz=inputs)
        while self.events and self.events[0][0]<=end+1e-12:
            _,_,kind,key,value=heapq.heappop(self.events)
            if kind=='sensor':self.arrived[key]=value
            else:self.excitations[key]=value
        rates=dict(zip(neural['regional_state']['node_ids'],neural['regional_state']['activity_hz']))
        targets={};gains={};drives={}
        for key in self.muscles:
            side=self.bindings[key]['side'];region=self.bindings[key]['motor_region']
            gains[key]=float(np.clip(1+p.cortical_gain_per_hz*(rates[region]-self.reference_rates[region]),0,2))
            drives[key]=descending.get(key,0)*float(np.clip(p.cortical_drive_per_hz*rates[region],0,1))
            sensor=self.arrived.get(key);reflex=0.
            if sensor and key in MUSCLES:
                if key.startswith('tibant'):
                    reflex=1.1*(sensor['length']-.71)
                    soleus=self.arrived.get('soleus_'+side)
                    if soleus and sensor['stance']:reflex-=.3*soleus['force']
                elif sensor['stance']:
                    if key.startswith('soleus'):reflex=1.2*sensor['force']
                    else:
                        group=[self.arrived.get(n+'_'+side) for n in ('gasmed','gaslat')]
                        # Both heads required for the explicitly lumped GAS transfer.
                        if all(group):reflex=1.1*sum(s['force_n'] for s in group)/sum(s['fmax_n'] for s in group)
            baseline=.01 if key in MUSCLES else 0.
            target=float(np.clip(baseline+gains[key]*reflex+drives[key],baseline,1))
            targets[key]=0. if key in mb else target
            if key not in mb:self._queue(end+p.efferent_delay_s,'motor',key,target)
        # Zero-delay motor commands are available at this endpoint only.
        while self.events and self.events[0][0]<=end+1e-12:
            _,_,kind,key,value=heapq.heappop(self.events)
            if kind=='motor':self.excitations[key]=value
            else:self.arrived[key]=value
        self.time_s=end
        return {'schema':'ihm.sensorimotor.v1','time_s':end,'motor_excitations':dict(self.excitations),
            'requested_excitations':targets,'sensors':deepcopy(sensors),'delayed_sensors':deepcopy(self.arrived),
            'brain':neural,'additional_sensory_inputs_hz':dict(additional),'brain_sensory_inputs_hz':dict(inputs),'descending_gain':gains,'descending_drive':drives,'pending_events':len(self.events),
            'sensory_blocks':sorted(sb),'motor_blocks':sorted(mb),'exchange_interval_s':dt,
            'neural_delay_s':p.afferent_delay_s+p.efferent_delay_s,'activation_owner':'mechanical_plant',
            'biological_validation':False,'decoder_basis':'engineered regional-rate gain and effector-gated drive; not identified motor recruitment',
            'scope':'all catalog effectors have engineered sensory/descending ports; eight ankle effectors have source reflex primitives; no autonomous walking policy',
            'spinal_reflex_effectors':list(MUSCLES),'descending_effector_count':len(self.muscles),
            'model_sha256':self.model_sha256}

    def checkpoint(self):
        return {'schema':'ihm.sensorimotor-state.v1','model_sha256':self.model_sha256,'time_s':self.time_s,
                'brain_state':self.brain.state.tolist(),'arrived':deepcopy(self.arrived),
                'excitations':dict(self.excitations),'events':deepcopy([list(e) for e in sorted(self.events)]),'serial':self.serial}

    def restore(self,checkpoint):
        c=deepcopy(checkpoint)
        if not isinstance(c,dict) or c.get('schema')!='ihm.sensorimotor-state.v1' or c.get('model_sha256')!=self.model_sha256:
            raise ValueError('Checkpoint identity mismatch')
        time=finite(c.get('time_s'),0,1e9,'checkpoint time');brain=np.asarray(c.get('brain_state'),float)
        if brain.shape!=self.brain.state.shape or not np.isfinite(brain).all():raise ValueError('Invalid brain checkpoint')
        def sensor(s):
            if not isinstance(s,dict) or type(s.get('stance')) is not bool:raise ValueError('Invalid sensor checkpoint')
            for k in ('length','fmax_n'):finite(s.get(k),1e-12,1e12,k)
            for k in ('force','force_n'):finite(s.get(k),0,1e12,k)
            if not isinstance(s.get('sensor_basis'),str) or s.get('length_field') not in ('fiber_length_m','fiber_length_proxy_m'):
                raise ValueError('Invalid sensor provenance')
        arrived=c.get('arrived');exc=c.get('excitations')
        if not isinstance(arrived,dict) or set(arrived)-set(self.muscles) or not isinstance(exc,dict) or set(exc)!=set(self.muscles):raise ValueError('Invalid effector state')
        for s in arrived.values():sensor(s)
        for v in exc.values():finite(v,0,1,'excitation')
        serial=c.get('serial');events=c.get('events')
        if type(serial) is not int or serial<0 or not isinstance(events,list):raise ValueError('Invalid event state')
        previous=(time,-1);seen=set();validated=[]
        for e in events:
            if not isinstance(e,list) or len(e)!=5:raise ValueError('Invalid event')
            t,n,kind,key,value=e;finite(t,time,1e9,'arrival')
            if type(n) is not int or not 0<n<=serial or n in seen or (t,n)<=previous or key not in self.muscles or kind not in ('sensor','motor'):raise ValueError('Invalid event order/owner')
            seen.add(n);previous=(t,n)
            if kind=='sensor':sensor(value)
            else:finite(value,0,1,'queued excitation')
            validated.append(tuple(e))
        self.time_s=time;self.brain.time_s=time;self.brain.state=brain.copy()
        self.arrived=arrived;self.excitations=exc;self.serial=serial;self.events=validated
