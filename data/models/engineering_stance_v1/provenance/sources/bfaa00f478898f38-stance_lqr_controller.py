"""Strict 10 ms engineered native LQR stance adapter; no cortical motor claim."""
from copy import deepcopy
import hashlib,json,math,tempfile
from pathlib import Path
SOURCE_BYTES=Path(__file__).read_bytes()

class EngineeredLQRStanceController:
    control_interval_s=.01
    current_interval_actuation=True
    @classmethod
    def from_root(cls,root,*,muscle_catalog,artifact_path=None,registration_path=None,**kwargs):
        from ihm.assembly.sensorimotor import SensorimotorController
        root=Path(root);self=cls()
        artifact=Path(artifact_path or root/'data/models/engineering_stance_v1/linearization.npz')
        registration=Path(registration_path or root/'data/models/engineering_stance_v1/registration.json')
        self.artifact_bytes=artifact.read_bytes();self.registration_bytes=registration.read_bytes()
        self.registration=json.loads(self.registration_bytes)
        source=root/'ihm/native/stance_lqr.py';raw=source.read_bytes()
        namespace={'__name__':'_ihm_frozen_stance_lqr','__file__':str(source)}
        exec(compile(raw,str(source),'exec'),namespace)
        # NativeStanceLQR accepts a path; the temporary contains exactly the
        # already frozen artifact bytes and is removed immediately after load.
        with tempfile.TemporaryDirectory(prefix='ihm-lqr-load-') as temporary:
            retained=Path(temporary)/'linearization.npz';retained.write_bytes(self.artifact_bytes)
            self.policy=namespace['NativeStanceLQR'](retained,model_sha256=self.registration['model_sha256'],dt_s=.01)
        if self.policy.target_mass_kg is None:raise ValueError('LQR artifact lacks explicit native identification mass')
        if set(self.policy.muscle_names)!={row['id'] for row in muscle_catalog}:raise ValueError('LQR controller catalog differs from registered native model')
        self.inner=SensorimotorController.from_root(root,muscle_catalog=muscle_catalog,**kwargs)
        self.source_bytes={'stance_lqr_controller.py':SOURCE_BYTES,'stance_lqr.py':raw,
            'linearization.npz':self.artifact_bytes,'registration.json':self.registration_bytes}
        self.identity={'schema':'ihm.engineered-lqr-stance-controller.v1','inner_model_sha256':self.inner.model_sha256,
            'source_sha256':{name:hashlib.sha256(data).hexdigest() for name,data in self.source_bytes.items()},
            'native_model_sha256':self.policy.model_sha256,'sampling_interval_s':.01,'native_target_mass_kg':self.policy.target_mass_kg,
            'muscles':list(self.policy.muscle_names),'state_names':list(self.policy.state_names)}
        self.model_sha256=hashlib.sha256(json.dumps(self.identity,sort_keys=True).encode()).hexdigest()
        self.controller_metadata={'kind':'engineering_stance','motor_owner':'native_local_discrete_lqr',
            'cortical_motor_output_active':False,'trained_motor_policy':False,'walking_demonstrated':False,
            'biological_validation':False,'model_sha256':self.model_sha256,'sampling_interval_s':.01,
            'native_target_mass_kg':self.policy.target_mass_kg,
            'native_model_sha256':self.policy.model_sha256,'artifact_sha256':hashlib.sha256(self.artifact_bytes).hexdigest(),
            'scope':'Bounded engineering stance feedback; nonlinear validation applies only to separately recorded runs',
            'sensory_basis':'Privileged native joint/muscle state; any sensory block disables all stance output',
            'arc_availability':dict.fromkeys(('stretch','reciprocal','autogenic','renshaw'),False)}
        self.native=None;self.steps=0;self.excitations=dict.fromkeys(self.policy.muscle_names,0.);self.last_diagnostics={}
        return self
    def __getattr__(self,name):return getattr(self.inner,name)
    def bind_native(self,native):
        if self.native is not None and self.native is not native:raise ValueError('Cannot replace bound LQR plant')
        execution=json.loads((native.output/'execution.json').read_bytes())
        digest=execution['source_sha256']['subject_walk_scaled.osim']
        actual=hashlib.sha256((native.output/'inputs/subject_walk_scaled.osim').read_bytes()).hexdigest()
        registration=(native.output/'inputs/augmentation_registration.json').read_bytes()
        if digest!=self.policy.model_sha256 or actual!=digest or registration!=self.registration_bytes:raise ValueError('LQR/native source model or registration identity mismatch')
        if abs(execution['target_mass_kg']-self.policy.target_mass_kg)>1e-9:raise ValueError('LQR native target mass differs from artifact identification mass')
        state=native.snapshot()
        if state.get('environment')!='upright':raise ValueError('LQR stance requires upright native environment')
        self.policy.state_vector(state)
        self.native=native;return self
    def retain_sources(self,output):
        folder=Path(output)/'engineering-lqr-stance';folder.mkdir(parents=True,exist_ok=True)
        for name,data in self.source_bytes.items():(folder/name).write_bytes(data)
        (folder/'manifest.json').write_text(json.dumps(self.identity,indent=2)+'\n')
        if hasattr(self.inner,'retain_sources'):self.inner.retain_sources(output)
        return str(folder)
    def step(self,dt_s,mechanical_observation,*,descending=None,sensory_blocks=(),motor_blocks=(),physiology=None,additional_sensory_inputs_hz=None):
        if isinstance(dt_s,bool) or not math.isfinite(dt_s) or abs(dt_s-.01)>1e-12:raise ValueError('LQR requires exact10ms sampling')
        if self.native is None:raise ValueError('LQR requires bind_native before stepping')
        state=self.native.snapshot()
        if state['time_s']!=mechanical_observation['time_s'] or state['coordinates']!=mechanical_observation['joints'] or state['muscles']!=mechanical_observation['muscles']:raise ValueError('LQR observation differs from bound native state')
        # Mass transfer changes the identified plant and is not covered by K.
        if abs(float(state['mass_kg'])-self.policy.target_mass_kg)>1e-6:raise ValueError('LQR native mass changed from artifact identification basis')
        saved=self.checkpoint()
        try:
            result=self.inner.step(dt_s,mechanical_observation,descending=descending,sensory_blocks=sensory_blocks,
                motor_blocks=motor_blocks,physiology=physiology,additional_sensory_inputs_hz=additional_sensory_inputs_hz)
            blocked=bool(sensory_blocks)
            if blocked:commands=dict.fromkeys(self.policy.muscle_names,0.);diagnostics={'sensory_blocked':True}
            else:commands,diagnostics=self.policy.commands(state)
            for name in motor_blocks:commands[name]=0.
            self.excitations=commands;self.last_diagnostics=diagnostics;self.steps+=1
            result['inactive_neural_commands']=result.get('motor_excitations',{})
            result.update(motor_excitations=dict(commands),requested_excitations=dict(commands),controller=deepcopy(self.controller_metadata),
                model_sha256=self.model_sha256,arc_max=dict.fromkeys(('stretch','reciprocal','autogenic','renshaw'),0.),
                lqr_stance=deepcopy(diagnostics),neural_delay_s=0.,scope=self.controller_metadata['scope'])
            return result
        except Exception:self.restore(saved);raise
    def checkpoint(self):
        return {'schema':'ihm.engineered-lqr-stance-state.v1','model_sha256':self.model_sha256,
            'native_identity':None if self.native is None else self.native.identity,'inner':self.inner.checkpoint(),
            'steps':self.steps,'excitations':dict(self.excitations),'diagnostics':deepcopy(self.last_diagnostics)}
    def restore(self,checkpoint):
        c=deepcopy(checkpoint)
        if c.get('schema')!='ihm.engineered-lqr-stance-state.v1' or c.get('model_sha256')!=self.model_sha256 or c.get('native_identity')!=(None if self.native is None else self.native.identity):raise ValueError('LQR checkpoint identity mismatch')
        if type(c['steps']) is not int or c['steps']<0 or set(c['excitations'])!=set(self.policy.muscle_names) or any(isinstance(v,bool) or not math.isfinite(v) or not 0<=v<=1 for v in c['excitations'].values()):raise ValueError('Invalid LQR checkpoint state')
        self.inner.restore(c['inner']);self.steps=c['steps'];self.excitations=c['excitations'];self.last_diagnostics=c['diagnostics']
