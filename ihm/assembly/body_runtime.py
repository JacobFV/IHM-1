"""Transactional reduced-body execution against immutable native telemetry.

Native gas/blood/heat remain externally owned. Prescribed thoracic work is
reported as boundary work, not silently fed back or counted as conserved native
energy. This runtime does not pretend that replay is a two-way native solve.
"""
import copy
import math
from .mechanics import BodyMechanics
from .brain import BodyBrain
from .respiration import BodyRespiration
from .peripheral import BodyPeripheral


class BodyRuntime:
    _mutable={
        'mechanics':('x','v','omega','r','deformation','time','work','dissipation','affine_work','_soft_energy','prescribed_work','_soft_previous','prescribed_reactions','orientation_reactions'),
        'brain':('state','time_s'),
        'respiration':('q','v','time','work','dissipation','numerical_dissipation'),
        'peripheral':('time_s','receptors','proprioceptors','arrived','motor_targets','activations','sent','events','serial'),
    }

    def __init__(self,root,assets,bindings,reference):
        from .body_states import body_state_registry
        from .body import digest
        from pathlib import Path
        self.registry=body_state_registry(assets,digest(Path(root)/'data/derived/canonical/anatomy.json'))
        self.mechanics=BodyMechanics.from_dict(assets['mechanics'])
        self.brain=BodyBrain.from_dict(assets['brain'],root=root)
        self.respiration=BodyRespiration.from_dict(assets['respiration'])
        self.peripheral=BodyPeripheral.from_dict(assets['peripheral'])
        self.bindings=bindings;self.reference=copy.deepcopy(reference)
        self.time_s=0.;self.output=None

    def checkpoint(self):
        return {'time_s':self.time_s,'output':copy.deepcopy(self.output),
            'engines':{name:{key:copy.deepcopy(getattr(getattr(self,name),key)) for key in keys if hasattr(getattr(self,name),key)} for name,keys in self._mutable.items()}}

    def restore(self,checkpoint):
        self.time_s=checkpoint['time_s'];self.output=copy.deepcopy(checkpoint['output'])
        for name,keys in self._mutable.items():
            engine=getattr(self,name);values=checkpoint['engines'][name]
            for key in keys:
                if key in values:setattr(engine,key,copy.deepcopy(values[key]))
                elif hasattr(engine,key):delattr(engine,key)

    def step(self,dt,physiology,compartments,inputs):
        from .body import project_volumes
        if isinstance(dt,bool) or not math.isfinite(dt) or not 0<dt<=1:raise ValueError('Body interval must be in (0,1] s')
        snapshot=self.checkpoint()
        try:
            ratios=project_volumes(self.bindings,compartments,self.reference)
            lungs={key:value for key,value in ratios.items() if key in self.respiration.lung_ids}
            delta=sum(compartments[side+'LungPulmonaryGasVolume(mL)']-self.reference[side+'LungPulmonaryGasVolume(mL)'] for side in ['Left','Right'])*1e-6
            respiratory=self.respiration.step(dt,{'volume_change_m3':delta,'lung_volume_ratios':lungs})
            prior=self.output['peripheral'] if self.output else {}
            drivers=dict(respiratory['mechanics_drivers'])
            drivers['volume_ratios']={key:value for key,value in ratios.items() if key not in lungs}
            drivers['activation']=prior.get('motor_activations',{})
            mechanical=self.mechanics.step(dt,drivers)
            neural=self.brain.step(dt,{'mean_arterial_pressure_mmHg':physiology['MeanArterialPressure(mmHg)'],'oxygen_saturation':physiology['OxygenSaturation'],'core_temperature_C':physiology['CoreTemperature(degC)']},sensory_inputs_hz=prior.get('brain_inputs_hz',{}))
            peripheral=self.peripheral.step(dt,inputs['stimuli'],mechanical_state=mechanical,brain_state={'motor_commands':inputs['motor_commands']},blocked_nerves=inputs['blocked_nerves'])
            end=self.time_s+dt
            if any(abs(state['time_s']-end)>1e-8 for state in [mechanical,neural,respiratory,peripheral]):raise RuntimeError('Body component clocks diverged')
            self.output={'time_s':end,'mechanics':mechanical,'brain':neural,'respiration':respiratory,'peripheral':peripheral,'volume_ratios':ratios}
            self.time_s=end
            return self.output
        except Exception:
            self.restore(snapshot)
            raise
