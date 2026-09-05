"""Causal excess-rate somatic peripheral reduction, with named inferred routes.

No baseline firing or tonic motor command is invented. Inputs are held over each
step; finite transport delays queue timestamped rate changes. Output at t+dt is
suitable for the NEXT brain/mechanics interval (explicit co-simulation).
"""
from __future__ import annotations
import heapq
import math
import numpy as np


def scalar(value,lo,hi,label):
    if isinstance(value,bool):raise ValueError(f'{label} must be numeric')
    try:value=float(value)
    except (TypeError,ValueError):raise ValueError(f'{label} must be numeric') from None
    if not math.isfinite(value) or not lo<=value<=hi:raise ValueError(f'{label} outside [{lo},{hi}]')
    return value


class BodyPeripheral:
    @classmethod
    def from_dict(cls,data):return cls(data)

    def __init__(self,data):
        self.data=data;self.parameters=data['parameters'];self.time_s=0.
        self.patches={p['id']:p for p in data['receptor_patches']}
        self.bindings={m['muscle_id']:m for m in data['muscle_bindings']}
        self.receptors={};self.proprioceptors={};self.arrived={};self.motor_targets={}
        self.activations={};self.sent={};self.events=[];self.serial=0

    def _send(self,kind,key,value,arrival,meta):
        identity=(kind,key)
        if abs(value-self.sent.get(identity,0.))<=1e-9:return
        self.sent[identity]=value;self.serial+=1
        heapq.heappush(self.events,(arrival,self.serial,kind,key,value,meta))

    def _advance_activation(self,dt):
        decay=math.exp(-dt/self.parameters['activation_tau_s'])
        for key in self.motor_targets.keys()|self.activations.keys():
            target=self.motor_targets.get(key,0.)
            self.activations[key]=target+(self.activations.get(key,0.)-target)*decay

    def step(self,dt_s,stimuli=None,mechanical_state=None,brain_state=None):
        dt=scalar(dt_s,.0001,1.,'dt_s');stimuli={} if stimuli is None else stimuli
        mechanical_state={} if mechanical_state is None else mechanical_state
        brain_state={} if brain_state is None else brain_state
        if not all(isinstance(x,dict) for x in [stimuli,mechanical_state,brain_state]):raise ValueError('Inputs must be objects')
        if set(stimuli)-self.patches.keys():raise ValueError('Unknown receptor patch')
        validated={}
        for key,value in stimuli.items():
            if not isinstance(value,dict) or set(value)-{'pressure_pa','temperature_C','stretch_fraction'}:raise ValueError('Unknown stimulus modality')
            validated[key]={k:scalar(v,*{'pressure_pa':(0.,1e6),'temperature_C':(-20.,80.),'stretch_fraction':(0.,1.)}[k],k) for k,v in value.items()}
        commands=brain_state.get('motor_commands',{})
        regional=brain_state.get('regional_motor_drive_hz',{})
        if not isinstance(commands,dict) or set(commands)-self.bindings.keys():raise ValueError('Unknown or unsupported motor muscle')
        commands={k:scalar(v,0.,1.,'motor command') for k,v in commands.items()}
        known_regions={b['brain_motor_id'] for b in self.bindings.values()}
        if not isinstance(regional,dict) or set(regional)-known_regions:raise ValueError('Unknown motor brain region')
        regional={k:scalar(v,0.,1000.,'regional motor drive') for k,v in regional.items()}
        # Compute and validate feedback before committing any runtime mutation.
        feedback={};entities=mechanical_state.get('entities',{});forces=mechanical_state.get('muscle_forces_n',{})
        if not isinstance(entities,dict) or not isinstance(forces,dict):raise ValueError('Malformed mechanics feedback')
        for key,b in self.bindings.items():
            e=entities.get(b['canonical_entity_id'],{});f=np.asarray(e.get('deformation_gradient',np.eye(3)),float)
            if f.shape!=(3,3) or not np.isfinite(f).all():raise ValueError('Invalid mechanical deformation')
            axis=np.asarray(b['fiber_axis']);stretch=max(0.,float(np.linalg.norm(f@axis))-1.)
            # Actual mechanics path length includes bone motion and soft anchors.
            points=[]
            for a in b['anchors']:
                ent=entities.get(a['entity_id'],{})
                shift=np.asarray(ent.get('translation_m',[0.,0.,0.]),float)
                if shift.shape!=(3,) or not np.isfinite(shift).all():raise ValueError('Invalid mechanical translation')
                points.append(np.asarray(a['point_m'])+shift)
            length=sum(float(np.linalg.norm(y-x)) for x,y in zip(points,points[1:]))
            stretch=max(stretch,max(0.,length/b['rest_path_length_m']-1.))
            force=scalar(forces.get(key,0.),0.,1e9,'muscle force')
            feedback[key]=min(200.,200.*stretch+40.*force/max(1.,b['max_isometric_force_n']))
        p=self.parameters;end=self.time_s+dt;decay=math.exp(-dt/p['receptor_tau_s'])
        # Four channels distinguish warm C and cold slow afferent delay.
        for key,patch in self.patches.items():
            v=validated.get(key,{});temp=v.get('temperature_C',p['baseline_skin_temperature_C'])
            targets={'pressure':v.get('pressure_pa',0.)*p['pressure_gain_hz_pa'],
                     'stretch':v.get('stretch_fraction',0.)*p['stretch_gain_hz'],
                     'warm':max(0.,temp-p['baseline_skin_temperature_C'])*p['temperature_gain_hz_C'],
                     'cold':max(0.,p['baseline_skin_temperature_C']-temp)*p['temperature_gain_hz_C']}
            for modality,target in targets.items():
                channel=f'{key}:{modality}';target=min(p['max_receptor_rate_hz'],target)
                rate=target+(self.receptors.get(channel,0.)-target)*decay;self.receptors[channel]=rate
                speed=p.get(f'{modality}_velocity_m_s',p['tactile_velocity_m_s'])
                delay=patch['path_length_m']/speed+p['central_afferent_delay_s']
                self._send('afferent',channel,rate,end+delay,patch)
        cutaneous_gates={(meta['side'],meta.get('sensorimotor_region'))
                         for channel,(value,meta) in self.arrived.items()
                         if not channel.startswith('proprio:') and value>1e-8}
        for key,b in self.bindings.items():
            rate=feedback[key]+(self.proprioceptors.get(key,0.)-feedback[key])*decay
            self.proprioceptors[key]=rate
            self._send('afferent','proprio:'+key,rate,end+b['afferent_delay_s'],b)
            gated=(b['side'],b.get('sensorimotor_region')) in cutaneous_gates
            command=commands.get(key,min(1.,regional.get(b['brain_motor_id'],0.)*b['somatic_readout_gain']) if gated else 0.)
            self._send('motor',key,command,self.time_s+b['motor_delay_s'],b)
        cursor=self.time_s
        while self.events and self.events[0][0]<=end+1e-12:
            arrival,_,kind,key,value,meta=heapq.heappop(self.events)
            self._advance_activation(max(0.,arrival-cursor));cursor=max(cursor,arrival)
            if kind=='motor':self.motor_targets[key]=value
            else:self.arrived[key]=(value,meta)
        self._advance_activation(max(0.,end-cursor));self.time_s=end
        brain={};relay={};nerve={}
        for value,meta in self.arrived.values():
            for output,k in [(brain,'brain_target_id'),(relay,'relay_id'),(nerve,'nerve_id')]:
                output[meta[k]]=output.get(meta[k],0.)+value
        for key,value in self.activations.items():
            b=self.bindings[key];nerve[b['nerve_id']]=nerve.get(b['nerve_id'],0.)+value*100.
        return {'schema_version':1,'model_id':self.data['id'],'time_s':end,
                'brain_inputs_hz':{k:min(1000.,v) for k,v in brain.items()},
                'motor_activations':dict(self.activations),'receptor_rates_hz':dict(self.receptors),
                'proprioceptor_rates_hz':dict(self.proprioceptors),'relay_activity_hz':relay,'nerve_activity_hz':nerve,
                'pending_events':len(self.events),'stimuli':validated,
                'scope':'Excess evoked-rate reduction; delayed somatic motor and sensory priors; no autonomic controller',
                'biological_validation':False,'coupling_application':'Caller feeds brain_inputs_hz and motor_activations to subsequent brain/mechanics intervals'}
