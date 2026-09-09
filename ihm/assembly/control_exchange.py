"""Explicit current-boundary feedback inside a slower physiology interval."""
from copy import deepcopy
import math
import numpy as np
from .world_exchange import advance_body_world


def advance_feedback_exchange(plant,controller,world,dt,state,forces,inputs,object_forces=(),respiratory_projector=None):
    h=controller.control_interval_s
    if isinstance(h,bool) or not isinstance(h,(int,float)) or not math.isfinite(h) or h<=0 or h>dt:
        raise ValueError('Invalid feedback control interval')
    count=round(dt/h)
    if abs(count*h-dt)>1e-10:raise ValueError('Feedback interval must divide physiology exchange')
    anchors=[]
    for force in forces:
        entity=state['entities'].get(force['id']);local=None
        if entity is not None:
            local=np.asarray(entity.get('rotation_matrix',np.eye(3))).T@(np.asarray(force['point_m'])-entity['centroid_m'])
        anchors.append((deepcopy(force),local))
    loads=[];samples=[];positive_work=0.;substeps=0;control_times=[]
    observation_ids=set(getattr(controller,'observation_entity_ids',()))
    observation_ids.update(force['id'] for force in forces)
    for index in range(count):
        start=state['time_s'];control_times.append(start)
        neural=controller.step(h,state,**inputs)
        if abs(neural['time_s']-start-h)>1e-8:raise RuntimeError('Feedback control clock diverged')
        actuation={k:v for k,v in neural['motor_excitations'].items() if k not in inputs.get('motor_blocks',())}
        applied=[]
        for force,local in anchors:
            port=deepcopy(force)
            if local is not None:
                entity=state['entities'][force['id']]
                port['point_m']=(np.asarray(entity['centroid_m'])+np.asarray(entity.get('rotation_matrix',np.eye(3)))@local).tolist()
            applied.append(port)
        if world is None:
            if respiratory_projector is not None:
                sample=respiratory_projector(applied,state['entities'])
                for port in sample['force_ports']:port['sample_time_s']=start
                samples.append(sample)
            if hasattr(plant,'advance_observation') and index<count-1:
                state=plant.advance_observation(h,forces=applied,actuation=actuation,entity_ids=observation_ids)
            else:state=plant.advance(h,forces=applied,actuation=actuation)
            ports=applied;substeps+=1
        else:
            state,ports,audit=advance_body_world(plant,world,h,state,applied,actuation,object_forces,
                respiratory_projector=respiratory_projector,final_full=index==count-1,
                observation_entity_ids=observation_ids)
            substeps+=audit['substeps']
            if 'respiratory_load' in audit:samples.append(audit['respiratory_load'])
        if abs(state['time_s']-start-h)>1e-8:raise RuntimeError('Feedback mechanical clock diverged')
        positive_work+=state['positive_muscle_work_j']
        loads.extend({**port,'force_n':(np.asarray(port['force_n'])/count).tolist()} for port in ports)
    state['positive_muscle_work_j']=positive_work
    audit={'interval_s':dt/substeps,'substeps':substeps,'control_interval_s':h,
        'control_substeps':count,'control_observation_times_s':control_times,
        'motor_exchange_latency_s':0.,'actuation_basis':'Feedback computed at current control boundary and held through that interval',
        'biological_validation':False}
    if samples:
        load=deepcopy(samples[0]);load.update(external_pressure_pa=0.,force_ports=[],ignored_nonrespiratory_ids=[])
        for sample in samples:
            load['external_pressure_pa']+=sample['external_pressure_pa']/count
            for port in sample['force_ports']:
                weighted=deepcopy(port)
                for key in ('force_n','moment_nm'):
                    if key in weighted:weighted[key]=(np.asarray(weighted[key])/count).tolist()
                weighted['generalized_force_pa']/=count
                weighted['quadrature_weight']=weighted.get('quadrature_weight',1.)/count
                load['force_ports'].append(weighted)
            load['ignored_nonrespiratory_ids'].extend(sample['ignored_nonrespiratory_ids'])
        load['ignored_nonrespiratory_ids']=list(dict.fromkeys(load['ignored_nonrespiratory_ids']))
        load['quadrature']={'substeps':substeps,'interval_s':dt/substeps,
            'basis':'Mean per-pose generalized load across control and mechanical substeps; native lung volume held over physiology exchange'}
        audit['respiratory_load']=load
    return neural,state,loads,audit
