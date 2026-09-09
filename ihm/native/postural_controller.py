"""Unvalidated engineering stance candidate on the unified native protocol.

Regional brain activity is observed separately. Muscle allocation owns motor
output; no trained brain, walking, or stable balance claim follows from this mode.
"""
from copy import deepcopy
import hashlib,io,json,math
from pathlib import Path
SOURCE_BYTES=Path(__file__).read_bytes()

class EngineeringStanceController:
    @classmethod
    def from_root(cls,root,*,muscle_catalog,policy_config=None,**kwargs):
        from ihm.assembly.sensorimotor import SensorimotorController
        root=Path(root);self=cls()
        config={} if policy_config is None else deepcopy(policy_config)
        allowed={'kp','kd','baseline','pelvis_gain','pelvis_damping','regularization','targets','feedforward_torques','com_position_gain','com_velocity_gain','com_target_x_m','com_lateral_position_gain','com_lateral_velocity_gain','com_target_z_m','equilibrium_excitations'}
        if not isinstance(config,dict) or set(config)-allowed:raise ValueError('Unknown engineering stance policy configuration')
        self.inner=SensorimotorController.from_root(root,muscle_catalog=muscle_catalog,**kwargs)
        self.policy_config=config
        source=root/'ihm/native/moment_arm_control.py';raw=source.read_bytes()
        namespace={'__name__':'_ihm_frozen_stance_policy','__file__':str(source)}
        exec(compile(raw,str(source),'exec'),namespace);self.policy_class=namespace['JointPosturalController']
        xml=root/'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled_FunctionBasedPathSet.xml'
        self.path_bytes=xml.read_bytes()
        self.source_bytes={'postural_controller.py':SOURCE_BYTES,'moment_arm_control.py':raw,'muscle_paths.xml':self.path_bytes}
        self.identity={'schema':'ihm.engineering-stance-controller.v1','inner_model_sha256':self.inner.model_sha256,
            'policy_config':config,'source_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in self.source_bytes.items()},
            'muscles':list(self.inner.muscles),'native_moment_arm_refresh_steps':5}
        self.model_sha256=hashlib.sha256(json.dumps(self.identity,sort_keys=True,allow_nan=False).encode()).hexdigest()
        self.controller_metadata={'kind':'engineering_stance','motor_owner':'engineering_joint_torque_allocator',
            'trained_motor_policy':False,'cortical_motor_output_active':False,'stable_balance_validated':False,
            'walking_demonstrated':False,'biological_validation':False,'model_sha256':self.model_sha256,
            'sensory_basis':'privileged native coordinates and COM; any sensory block disables all stance output',
            'arc_availability':dict.fromkeys(('stretch','reciprocal','autogenic','renshaw'),False)}
        self.native=None;self.reference=None;self.policy=None;self.cached_arms={};self.steps=0
        self.excitations={m:0. for m in self.inner.muscles}
        return self
    def __getattr__(self,name):return getattr(self.inner,name)
    def bind_native(self,native):
        if self.native is not None and self.native is not native:raise ValueError('Cannot replace bound native plant')
        state=native.snapshot()
        if set(state['muscles'])!=set(self.inner.muscles):raise ValueError('Stance/native catalog mismatch')
        lumbar=[m for m in state['muscles'] if m.startswith('gait2392_')]
        if lumbar and len(lumbar)!=6:raise ValueError('Complete six-muscle donor lumbar set required')
        self.native=native;self.lumbar=lumbar
        return self
    def _new_policy(self):
        return self.policy_class(self.reference,path=io.BytesIO(self.path_bytes),native_moment_arms=self.cached_arms,**self.policy_config)
    def retain_sources(self,output):
        folder=Path(output)/'engineering-stance';folder.mkdir(parents=True,exist_ok=True)
        for name,raw in self.source_bytes.items():(folder/name).write_bytes(raw)
        (folder/'manifest.json').write_text(json.dumps(self.identity,indent=2)+'\n')
        if hasattr(self.inner,'retain_sources'):self.inner.retain_sources(output)
        return str(folder)
    def step(self,dt_s,mechanical_observation,*,descending=None,sensory_blocks=(),motor_blocks=(),physiology=None,additional_sensory_inputs_hz=None):
        if self.native is None:raise ValueError('Engineering stance requires bind_native before stepping')
        state=self.native.snapshot()
        if state['time_s']!=mechanical_observation['time_s'] or state['coordinates']!=mechanical_observation['joints'] or set(state['muscles'])!=set(mechanical_observation['muscles']):
            raise ValueError('Stance observation differs from bound current native plant')
        saved=self.checkpoint()
        try:
            result=self.inner.step(dt_s,mechanical_observation,descending=descending,sensory_blocks=sensory_blocks,
                motor_blocks=motor_blocks,physiology=physiology,additional_sensory_inputs_hz=additional_sensory_inputs_hz)
            if self.steps%5==0 and self.lumbar:
                self.cached_arms=self.native.moment_arms(muscles=self.lumbar,coordinates=['lumbar_extension','lumbar_bending','lumbar_rotation'])['moment_arms_m']
            if self.reference is None:self.reference=deepcopy(state);self.policy=self._new_policy()
            elif self.lumbar:self.policy.update_native_moment_arms(self.cached_arms)
            blocked=bool(sensory_blocks)
            commands=({m:0. for m in self.inner.muscles} if blocked else self.policy.commands(state))
            for name in motor_blocks:commands[name]=0.
            self.excitations=commands;self.steps+=1
            result['inactive_neural_commands']=result.get('motor_excitations',{})
            result['inactive_neural_arc_max']=result.get('arc_max',{})
            result.update(motor_excitations=dict(commands),requested_excitations=dict(commands),
                arc_max=dict.fromkeys(('stretch','reciprocal','autogenic','renshaw'),0.),
                controller=deepcopy(self.controller_metadata),model_sha256=self.model_sha256,
                engineering_stance={'allocation':{} if blocked else deepcopy(self.policy.last_allocation),'sensory_blocked':blocked,
                    'moment_arm_refresh_steps':5,'native_moment_arm_muscles':list(self.cached_arms)},
                neural_delay_s=0.,scope='Unvalidated engineering stance candidate; privileged native feedback; regional brain motor output inactive')
            return result
        except Exception:self.restore(saved);raise
    def checkpoint(self):
        return {'schema':'ihm.engineering-stance-state.v1','model_sha256':self.model_sha256,
            'native_identity':None if self.native is None else self.native.identity,
            'inner':self.inner.checkpoint(),'reference':deepcopy(self.reference),'cached_arms':deepcopy(self.cached_arms),
            'steps':self.steps,'excitations':dict(self.excitations),'last_allocation':{} if self.policy is None else deepcopy(self.policy.last_allocation)}
    def restore(self,checkpoint):
        c=deepcopy(checkpoint)
        if c.get('schema')!='ihm.engineering-stance-state.v1' or c.get('model_sha256')!=self.model_sha256 or c.get('native_identity')!=(None if self.native is None else self.native.identity):raise ValueError('Stance checkpoint identity mismatch')
        if type(c['steps']) is not int or c['steps']<0 or set(c['excitations'])!=set(self.inner.muscles) or any(isinstance(v,bool) or not math.isfinite(v) or not 0<=v<=1 for v in c['excitations'].values()):raise ValueError('Invalid stance checkpoint state')
        old_reference,old_arms=self.reference,self.cached_arms
        try:
            self.reference=c['reference'];self.cached_arms=c['cached_arms']
            policy=None if self.reference is None else self._new_policy()
            self.inner.restore(c['inner'])
        except Exception:self.reference,self.cached_arms=old_reference,old_arms;raise
        self.policy=policy
        if policy is not None:policy.last_allocation=c['last_allocation']
        self.steps=c['steps'];self.excitations=c['excitations']
