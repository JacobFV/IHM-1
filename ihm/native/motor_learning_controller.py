"""Explicit learned ankle primitive in the unified native brain/body/world loop.

The persistent IBM E/I cortex still advances, but its untrained motor readout is
inactive in this mode. The separate reduced trained kernel owns motor output.
"""
from copy import deepcopy
import hashlib,json,math
from pathlib import Path
from .motor_learning import load_motor_policy_bytes, SOURCE_BYTES as POLICY_SOURCE_BYTES
SOURCE_BYTES=Path(__file__).read_bytes()

class IBMAnklePrimitiveController:
    @classmethod
    def from_root(cls,root,*,muscle_catalog,artifact_path=None,target_rad=.12,**kwargs):
        from ihm.assembly.ibm_controller import IBMImplicitController
        path=Path(artifact_path or Path(root)/'data/models/ibm_ankle_primitive_v1/motor_kernel.pt')
        if not math.isfinite(target_rad) or abs(target_rad)>.5:raise ValueError('Bounded finite ankle target required')
        self=cls();self.inner=IBMImplicitController.from_root(root,muscle_catalog=muscle_catalog,**kwargs)
        self.artifact_bytes=path.read_bytes()
        self.policy,self.artifact=load_motor_policy_bytes(self.artifact_bytes)
        self.source_bytes={'motor_learning.py':POLICY_SOURCE_BYTES,'motor_learning_controller.py':SOURCE_BYTES}
        self.target_rad=float(target_rad)
        self.artifact_sha256=hashlib.sha256(self.artifact_bytes).hexdigest()
        self.identity={'inner':self.inner.identity,'artifact_sha256':self.artifact_sha256,
            'target_rad':self.target_rad,'source_provenance':self.artifact['provenance'],
            'adapter_sha256':hashlib.sha256(SOURCE_BYTES).hexdigest(),
            'policy_sha256':hashlib.sha256(POLICY_SOURCE_BYTES).hexdigest()}
        self.model_sha256=hashlib.sha256(json.dumps(self.identity,sort_keys=True).encode()).hexdigest()
        self.controller_metadata=dict(self.inner.controller_metadata,kind='implicit_ankle_primitive',
            trained_motor_policy=True,cortical_motor_output_active=False,motor_owner='reduced-eight-site-ankle-kernel',
            model_sha256=self.model_sha256,artifact_sha256=self.artifact_sha256,target_rad=self.target_rad,
            training='synthetic PD imitation; copied embedding only',cord_applied_to_motor_output=False,
            sensory_basis='privileged native ankle coordinate and speed; no anatomical encoder')
        self.controller_metadata['inactive_cord_arc_availability']=self.controller_metadata.get('arc_availability',{})
        self.controller_metadata['arc_availability']={name:False for name in ('stretch','reciprocal','autogenic','renshaw')}
        self.excitations=dict(self.inner.excitations)
        return self
    def __getattr__(self,name):return getattr(self.inner,name)
    def retain_sources(self,output):
        self.inner.retain_sources(output);folder=Path(output)/'ankle-motor-primitive';folder.mkdir(exist_ok=True)
        (folder/'motor_kernel.pt').write_bytes(self.artifact_bytes)
        (folder/'manifest.json').write_text(json.dumps(self.identity,indent=2)+'\n')
        for filename,raw in self.source_bytes.items():
            (folder/filename).write_bytes(raw)
        return str(folder)
    def step(self,dt_s,mechanical_observation,*,descending=None,sensory_blocks=(),motor_blocks=(),physiology=None,additional_sensory_inputs_hz=None):
        # Reject before changing cortical state if privileged port is unavailable.
        if 'joints' not in mechanical_observation:raise ValueError('Projected native joints required')
        joints=mechanical_observation['joints']
        if 'coordinates' in mechanical_observation and mechanical_observation['coordinates']!=joints:
            raise ValueError('Ambiguous native joint/coordinate observation')
        policy_observation={'coordinates':joints,'muscles':mechanical_observation['muscles']}
        q=joints['ankle_angle_r']
        if not all(math.isfinite(q[k]) for k in ('value','speed')):raise ValueError('Nonfinite ankle observation')
        saved=self.checkpoint()
        try:
            result=self.inner.step(dt_s,mechanical_observation,descending=descending,sensory_blocks=sensory_blocks,
                motor_blocks=motor_blocks,physiology=physiology,additional_sensory_inputs_hz=additional_sensory_inputs_hz)
            # Blocking either antagonist sensory channel invalidates this joint's
            # privileged feedback entirely; do not leak coordinates around block.
            sever=self.inner.sever or bool(set(sensory_blocks)&{'tibant_r','soleus_r'})
            commands=self.policy.commands(policy_observation,self.target_rad,sever=sever)
            availability=result['brain'].get('oxygen_perfusion_availability',1.)
            for name in commands:commands[name]=0. if name in motor_blocks else commands[name]*availability
            self.excitations=commands
            result.update(motor_excitations=dict(commands),requested_excitations=dict(commands),
                controller=deepcopy(self.controller_metadata),model_sha256=self.model_sha256,
                motor_primitive={'target_rad':self.target_rad,'angle_rad':q['value'],'artifact_sha256':self.artifact_sha256,
                    'sensory_blocked':sever and not self.inner.sever},
                scope='Persistent IBM E/I observed; reduced trained kernel owns ankle motor output; no walking',
                neural_delay_s=0.)
            # Retain inactive cortex/cord signals by explicit names, never pass
            # them off as the command actually delivered to native muscles.
            result['inactive_cortical_commands']=result.pop('cortical_commands')
            result['inactive_cord_arc_max']=result.pop('arc_max')
            result['arc_max']={key:0. for key in result['inactive_cord_arc_max']}
            return result
        except Exception:self.restore(saved);raise
    def checkpoint(self):
        return {'schema':'ihm.ibm-ankle-primitive-state.v1','model_sha256':self.model_sha256,
                'inner':self.inner.checkpoint(),'excitations':dict(self.excitations)}
    def restore(self,checkpoint):
        if checkpoint.get('schema')!='ihm.ibm-ankle-primitive-state.v1' or checkpoint.get('model_sha256')!=self.model_sha256:raise ValueError('Primitive checkpoint identity mismatch')
        commands=checkpoint['excitations']
        if set(commands)!=set(self.inner.muscles) or any(not math.isfinite(v) or not 0<=v<=1 for v in commands.values()):raise ValueError('Invalid primitive checkpoint commands')
        self.inner.restore(checkpoint['inner']);self.excitations=dict(commands)
