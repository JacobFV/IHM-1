"""Actual trained IBM E/I cortex owns an ankle motor primitive in native runtime.

The target and coordinate feedback are engineering ports, not anatomical sensory
reconstruction. Cortical E/I integration is persistent and uses the trained
copied association embedding; the original fused checkpoint is never modified.
"""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import torch
from ihm.assembly.reflexes import finite
from .cortical_motor import load_cortical_policy


class IBMCorticalAnkleController:
    @classmethod
    def from_root(cls, root, *, muscle_catalog, artifact_path=None, target_rad=.12, **kwargs):
        from ihm.assembly.ibm_controller import IBMImplicitController
        path=Path(artifact_path or Path(root)/'data/models/ibm_cortical_ankle_v1/cortical_motor.pt')
        target_rad=finite(target_rad,-.25,.25,'ankle target')
        kwargs=dict(kwargs,sites=128)
        self=cls()
        self.inner=IBMImplicitController.from_root(root,muscle_catalog=muscle_catalog,**kwargs)
        self.policy,self.artifact=load_cortical_policy(path)
        self.artifact_bytes=path.read_bytes()
        if hashlib.sha256(self.artifact_bytes).hexdigest()!=self.artifact['artifact_sha256']:
            raise ValueError('Trained artifact changed while loading')
        self.frozen_sources={name:Path(__file__).with_name(name).read_bytes() for name in
                             ('cortical_motor.py','cortical_motor_controller.py')}
        self.frozen_sources['pretrain_video_loop.py']=(path.parent/'pretrain_video_loop.py').read_bytes()
        self.target_rad=target_rad
        self.identity={'inner':self.inner.identity,'artifact_sha256':self.artifact['artifact_sha256'],
            'target_rad':target_rad,'sources':{k:hashlib.sha256(v).hexdigest() for k,v in self.frozen_sources.items()}}
        self.model_sha256=hashlib.sha256(json.dumps(self.identity,sort_keys=True).encode()).hexdigest()
        self.controller_metadata=dict(self.inner.controller_metadata,kind='implicit_cortical_ankle',
            trained_motor_policy=True,cortical_motor_output_active=True,
            motor_owner='trained128site-IBM-EI-cortex',model_sha256=self.model_sha256,
            artifact_sha256=self.artifact['artifact_sha256'],target_rad=target_rad,
            training='Synthetic ankle PD imitation; copied association embedding and precentral decoder',
            sensory_basis='Privileged native ankle coordinate/speed; fixed disjoint sensory/motor cortical ports',
            scope='One ankle target primitive; not walking or general motor competence')
        self.inner.cortical_state=self.policy.state()
        self.excitations=dict(self.inner.excitations)
        with torch.no_grad(): self.weights=self.policy.dyn.edge_weights().detach()
        return self

    def __getattr__(self,name): return getattr(self.inner,name)

    def retain_sources(self,output):
        self.inner.retain_sources(output)
        target=Path(output)/'trained-cortical-motor';target.mkdir(parents=True,exist_ok=True)
        (target/'cortical_motor.pt').write_bytes(self.artifact_bytes)
        (target/'manifest.json').write_text(json.dumps(self.identity,indent=2)+'\n')
        for name,raw in self.frozen_sources.items(): (target/name).write_bytes(raw)
        return str(target)

    def step(self,dt_s,mechanical_observation,*,descending=None,sensory_blocks=(),motor_blocks=(),
             physiology=None,additional_sensory_inputs_hz=None):
        dt=finite(dt_s,.001,.1,'dt_s');ticks=round(dt/.001)
        if abs(dt-ticks*.001)>1e-10: raise ValueError('Integral1ms ticks required')
        c=self.inner
        sensors=c._validate_observation(mechanical_observation)
        sb=c._blocks(sensory_blocks,'sensory block');mb=c._blocks(motor_blocks,'motor block')
        descending={} if descending is None else descending
        if not isinstance(descending,dict) or set(descending)-{'tibant_r','soleus_r'}:
            raise ValueError('Learned cortical primitive accepts only right-ankle descending requests')
        requests={k:finite(v,0,1,'descending request') for k,v in descending.items()}
        target=float(np.clip(self.target_rad+.2*(requests.get('tibant_r',0)-requests.get('soleus_r',0)),-.5,.5))
        joints=mechanical_observation.get('joints',mechanical_observation.get('coordinates',{}))
        q=joints.get('ankle_angle_r')
        if not isinstance(q,dict): raise ValueError('Native right ankle coordinate required')
        angle=finite(q.get('value'),-10,10,'ankle angle');speed=finite(q.get('speed'),-100,100,'ankle speed')
        additional={} if additional_sensory_inputs_hz is None else additional_sensory_inputs_hz
        if not isinstance(additional,dict) or set(additional)-set(c.channels): raise ValueError('Unknown sensory population')
        rates={k:finite(v,0,1000,'additional sensory rate') for k,v in additional.items()}
        from ihm.assembly.brain import INPUT_BASELINES
        supplied={} if physiology is None else physiology
        if not isinstance(supplied,dict) or set(supplied)-set(INPUT_BASELINES):raise ValueError('Unknown physiology input')
        values=dict(INPUT_BASELINES,**supplied)
        pressure=finite(values['mean_arterial_pressure_mmHg'],0,300,'MAP')
        oxygen=finite(values['oxygen_saturation'],0,1,'oxygen saturation')
        temperature=finite(values['core_temperature_C'],20,45,'temperature')
        priors=c.brain.data['parameters']['body_transfer_priors']
        availability=min(1.,pressure/priors['map_reference_mmHg'])*min(1.,oxygen/.98)
        temperature_factor=priors['temperature_Q10']**((temperature-37)/10)
        sensory_blocked=bool(sb & {'tibant_r','soleus_r'})
        x=torch.tensor([[0. if sensory_blocked else target-angle,0. if sensory_blocked else speed]],dtype=torch.float32)
        stretch=np.array([0. if m in sb else min(1.,max(0.,sensors[m]['length']-1)*2) for m in c.muscles],np.float32)
        force=np.array([0. if m in sb else min(1.,sensors[m]['force']) for m in c.muscles],np.float32)
        saved=self.checkpoint()
        try:
            arcs={k:0. for k in ('stretch','reciprocal','autogenic','renshaw')}
            for name in ('stretch','autogenic','reciprocal'):
                for sample in c.cord._delay_buf.get(name,[]):
                    for i,m in enumerate(c.muscles):
                        if name!='reciprocal' and m in sb or name=='reciprocal' and c.reciprocal_pairs.get(m) in sb: sample[i]=0.
            with torch.no_grad():
                for _ in range(ticks):
                    output,c.cortical_state=self.policy.advance(x,c.cortical_state,ticks=1,sever=c.sever,
                        availability=availability,temperature_factor=temperature_factor,
                        supplemental_drive=min(5.,.002*sum(rates.values())),weights=self.weights)
                    if not all(torch.isfinite(s).all() for s in c.cortical_state):raise ValueError('Nonfinite trained cortex')
                    cortical=np.zeros(len(c.muscles),np.float32)
                    cortical[c.muscles.index('tibant_r')]=float(output[0,0])
                    cortical[c.muscles.index('soleus_r')]=float(output[0,1])
                    if c.no_cord:alpha=cortical.copy()
                    else:
                        response=c.cord.step(cortical,stretch=stretch,force=force,antagonist=c.antagonist)
                        alpha=response['alpha'].copy()
                        for name in arcs:arcs[name]=max(arcs[name],float(np.abs(response[name]).max()))
                    for i,m in enumerate(c.muscles):
                        if m in mb:alpha[i]=0.
                self.excitations=dict(zip(c.muscles,map(float,alpha)))
                c.excitations=dict(self.excitations)
                c.time_s+=dt;c.brain.time_s=c.time_s
                region={'node_ids':['implicit-postcentral','implicit-precentral']}
                for field,n in (('potential_mV',0),('activity_hz',1),('adaptation_mV',2)):
                    region[field]=[float(c.cortical_state[n][0,ix].mean()) for ix in (self.policy.sensory_sites,self.policy.motor_sites)]
                return {'schema':'ihm.sensorimotor.v1','time_s':c.time_s,'motor_excitations':dict(self.excitations),
                    'requested_excitations':dict(self.excitations),'cortical_commands':dict(zip(c.muscles,map(float,cortical))),
                    'controller':deepcopy(self.controller_metadata),'arc_max':arcs,'model_sha256':self.model_sha256,
                    'brain':{'time_s':c.time_s,'model_id':'ibm-trained-cortical-ankle','source_identity':self.identity,
                        'regional_state':region,'physiology_inputs':values,'physiology_coupling_applied':True,
                        'oxygen_perfusion_availability':availability,'temperature_factor':temperature_factor},
                    'sensors':sensors,'brain_sensory_inputs_hz':rates,'activation_owner':'mechanical_plant',
                    'motor_primitive':{'target_rad':target,'angle_rad':angle,'sensory_blocked':sensory_blocked},
                    'exchange_interval_s':dt,'neural_delay_s':.030 if not c.no_cord else 0.,
                    'biological_validation':False,'scope':self.controller_metadata['scope']}
        except Exception:self.restore(saved);raise

    def checkpoint(self):
        return {'schema':'ihm.cortical-ankle-state.v1','model_sha256':self.model_sha256,
                'inner':self.inner.checkpoint(),'excitations':dict(self.excitations)}

    def restore(self,checkpoint):
        # Validate the complete two-layer state before changing either owner.
        # Excitations are duplicated for the public controller and inner neural
        # checkpoint, so accepting disagreeing copies would create a false replay.
        checkpoint=deepcopy(checkpoint)
        if (not isinstance(checkpoint,dict) or
                set(checkpoint)!={'schema','model_sha256','inner','excitations'} or
                checkpoint.get('schema')!='ihm.cortical-ankle-state.v1' or
                checkpoint.get('model_sha256')!=self.model_sha256):
            raise ValueError('Trained cortical checkpoint identity mismatch')
        commands=checkpoint['excitations'];inner=checkpoint['inner']
        if not isinstance(commands,dict) or set(commands)!=set(self.inner.muscles):
            raise ValueError('Invalid cortical motor checkpoint muscles')
        for value in commands.values():finite(value,0,1,'checkpoint excitation')
        if not isinstance(inner,dict) or set(inner)!=set(self.inner.checkpoint()):
            raise ValueError('Incomplete inner cortical checkpoint')
        if inner.get('excitations')!=commands:
            raise ValueError('Public and inner checkpoint commands differ')
        time=finite(inner.get('time_s'),0,1e9,'checkpoint time')
        if abs(time-round(time/.001)*.001)>1e-9:
            raise ValueError('Checkpoint time must lie on1ms neural grid')
        try:self.inner.restore(inner)
        except (TypeError,KeyError,RuntimeError) as exc:
            raise ValueError('Malformed inner cortical checkpoint') from exc
        self.excitations=dict(commands)
