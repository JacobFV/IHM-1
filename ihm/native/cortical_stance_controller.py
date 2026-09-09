"""Persistent IBM E/I motor owner for a bounded native standing task.

The native state encoder and zero-input neural reference are engineering ports.
Offline LQR supplied decoder targets; no LQR executes on the runtime motor path.
"""
from copy import deepcopy
import hashlib,json,math,tempfile
from pathlib import Path
import torch
from ihm.assembly.reflexes import finite
SOURCE_BYTES=Path(__file__).read_bytes()

# One adapter, several retained stance bundles. Each entry names the trained
# artifact and the IBM association kernel it was initialised from; the artifact's
# own provenance is what actually enforces that pairing at load time.
BUNDLES={
    'implicit_cortical_stance':{'artifact':'data/models/ibm_cortical_stance_v1/cortical_stance.pt',
        'kernel':None,
        'kernel_identity':'IBM-1 ckpt/ibm1_implicit.pt, the fused 34-source implicit kernel'},
    'implicit_curriculum16_stance':{'artifact':'data/models/ibm_curriculum16_stance_v1/cortical_stance.pt',
        'kernel':'data/models/ibm_curriculum16_kernel_v1/kernel.pt',
        'kernel_identity':'IBM-1 ckpt/ibm1_curriculum16.pt, the 16-objective consolidated kernel, retained locally'},
}

class IBMCorticalStanceController:
    control_interval_s=.01
    current_interval_actuation=True

    @classmethod
    def from_root(cls,root,*,muscle_catalog,kind='implicit_cortical_stance',artifact_path=None,registration_path=None,sever=False,no_cord=True,**kwargs):
        from ihm.assembly.ibm_controller import IBMImplicitController
        if type(sever) is not bool or type(no_cord) is not bool:raise ValueError('Controller ablations must be boolean')
        if kind not in BUNDLES:raise ValueError('Unknown cortical stance bundle selection')
        bundle=BUNDLES[kind]
        root=Path(root);self=cls()
        path=Path(artifact_path or root/bundle['artifact'])
        if bundle['kernel'] is not None:
            kernel=root/bundle['kernel']
            if not kernel.is_file():raise ValueError('Retained IBM kernel for this stance bundle is unavailable')
            kwargs=dict(kwargs,checkpoint_path=kernel.resolve())
        registration=Path(registration_path or root/'data/models/engineering_stance_v1/registration.json')
        artifact_bytes=path.read_bytes();equations=(path.parent/'pretrain_video_loop.py').read_bytes()
        self.registration_bytes=registration.read_bytes();self.registration=json.loads(self.registration_bytes)
        raw=(root/'ihm/native/cortical_stance.py').read_bytes()
        namespace={'__name__':'ihm.native._frozen_stance_policy','__package__':'ihm.native','__file__':str(root/'ihm/native/cortical_stance.py')}
        exec(compile(raw,namespace['__file__'],'exec'),namespace)
        with tempfile.TemporaryDirectory(prefix='ihm-cortical-stance-load-') as directory:
            temporary=Path(directory);(temporary/'cortical_stance.pt').write_bytes(artifact_bytes)
            (temporary/'pretrain_video_loop.py').write_bytes(equations)
            self.policy,self.artifact=namespace['load_cortical_stance'](temporary/'cortical_stance.pt',model_sha256=self.registration['model_sha256'],dt_s=.01)
        if not self.policy.reference_normalization:raise ValueError('Cortical stance runtime requires paired neural reference')
        if self.artifact['provenance'].get('target_mass_kg') is None:raise ValueError('Cortical stance artifact must bind identified native mass')
        if set(self.policy.muscle_names)!={row['id'] for row in muscle_catalog}:raise ValueError('Cortical stance muscle catalog differs from native model')
        # Reuse the peripheral catalog, donor identity and physiology parameters;
        # the raw implicit dynamics/cord never execute or own motor commands.
        kwargs=dict(kwargs,sites=self.policy.dyn.n)
        self.inner=IBMImplicitController.from_root(root,muscle_catalog=muscle_catalog,sever=sever,no_cord=True,**kwargs)
        if self.inner.identity['checkpoint_sha256']!=self.artifact['provenance']['brain_checkpoint_sha256']:raise ValueError('Shared IBM checkpoint differs from trained lineage')
        self.source_bytes={'cortical_stance_controller.py':SOURCE_BYTES,'cortical_stance.py':raw,
            'cortical_stance.pt':artifact_bytes,'pretrain_video_loop.py':equations,'registration.json':self.registration_bytes}
        self.identity={'schema':'ihm.ibm-cortical-stance-controller.v1','kind':kind,
            'bundle':bundle['artifact'],'kernel_bundle':bundle['kernel'],
            'brain_checkpoint_sha256':self.artifact['provenance']['brain_checkpoint_sha256'],
            'kernel_identity':bundle['kernel_identity'],'catalog_identity':self.inner.model_sha256,
            'source_sha256':{name:hashlib.sha256(value).hexdigest() for name,value in self.source_bytes.items()},
            'native_model_sha256':self.registration['model_sha256'],'target_mass_kg':self.artifact['provenance']['target_mass_kg'],
            'sampling_interval_s':.01,'sever':bool(sever),'muscles':list(self.policy.muscle_names),'state_names':list(self.policy.state_names)}
        self.model_sha256=hashlib.sha256(json.dumps(self.identity,sort_keys=True).encode()).hexdigest()
        self.controller_metadata=dict(self.inner.controller_metadata,kind=kind,
            model_sha256=self.model_sha256,artifact_sha256=self.artifact['artifact_sha256'],
            bundle=bundle['artifact'],kernel_bundle=bundle['kernel'],kernel_identity=bundle['kernel_identity'],
            brain_checkpoint_sha256=self.artifact['provenance']['brain_checkpoint_sha256'],
            motor_owner='trained1024site-persistent-IBM-EI-cortex',cortical_motor_output_active=True,
            trained_motor_policy=True,walking_demonstrated=False,biological_validation=False,
            sampling_interval_s=.01,native_model_sha256=self.registration['model_sha256'],
            native_target_mass_kg=self.identity['target_mass_kg'],sites=self.policy.dyn.n,
            state_width=len(self.policy.state_names),muscle_count=len(self.policy.muscle_names),
            sensory_sites=len(self.policy.sensory_sites),motor_sites=len(self.policy.motor_sites),
            computation_dtype=str(self.policy.dyn.embed.dtype),sever=bool(sever),no_cord=True,
            arc_availability=dict.fromkeys(('stretch','reciprocal','autogenic','renshaw'),False),
            sensory_basis='Privileged full native state through disjoint signed cortical ports; any sensory block flushes correction to tonic baseline',
            additional_sensory_motor_coupling=False,descending_motor_coupling=False,
            baseline_owner='Explicit fixed native equilibrium tonic excitation, retained in severed arm',
            reference_basis='Same-kernel persistent zero-input counterfactual subtracts intrinsic neural baseline',
            physiology_basis='Engineering MAP/O2 drive availability and temperature Q10 derivative modulation; recovery receipts cover baseline physiology only',
            scope='Bounded native stance and recorded small pelvis pushes; no walking or general brain competence claim')
        self.native=None;self.steps=0;self.time_s=0.;self.sever=bool(sever)
        self.cortical_state=self.policy.state();self.excitations=dict(zip(self.policy.muscle_names,map(float,self.policy.u0)))
        self.last_diagnostics={}
        with torch.no_grad():self.weights=self.policy.dyn.edge_weights().detach()
        return self

    def __getattr__(self,name):return getattr(self.inner,name)

    def bind_native(self,native):
        if self.native is not None and self.native is not native:raise ValueError('Cannot replace bound cortical plant')
        execution=json.loads((native.output/'execution.json').read_bytes())
        model=(native.output/'inputs/subject_walk_scaled.osim').read_bytes()
        if hashlib.sha256(model).hexdigest()!=self.identity['native_model_sha256'] or execution['source_sha256']['subject_walk_scaled.osim']!=self.identity['native_model_sha256']:raise ValueError('Cortical/native source model mismatch')
        if (native.output/'inputs/augmentation_registration.json').read_bytes()!=self.registration_bytes:raise ValueError('Cortical/native registration differs')
        actual_mass=finite(execution['target_mass_kg'],.001,1000,'native mass')
        if abs(actual_mass-self.identity['target_mass_kg'])>1e-9:raise ValueError('Cortical/native mass differs')
        snapshot=native.snapshot()
        if snapshot.get('environment')!='upright':raise ValueError('Cortical stance requires upright mechanics')
        self.policy.snapshot_vector(snapshot);self.native=native
        return self

    def retain_sources(self,output):
        folder=Path(output)/'trained-cortical-stance';folder.mkdir(parents=True,exist_ok=True)
        for name,raw in self.source_bytes.items():(folder/name).write_bytes(raw)
        (folder/'manifest.json').write_text(json.dumps(self.identity,indent=2)+'\n')
        self.inner.retain_sources(output)
        return str(folder)

    def step(self,dt_s,mechanical_observation,*,descending=None,sensory_blocks=(),motor_blocks=(),physiology=None,additional_sensory_inputs_hz=None):
        dt=finite(dt_s,.01,.01,'cortical stance interval')
        if self.native is None:raise ValueError('Cortical stance requires bound native plant')
        state=self.native.snapshot()
        if state['time_s']!=mechanical_observation['time_s'] or state['coordinates']!=mechanical_observation['joints'] or state['muscles']!=mechanical_observation['muscles']:raise ValueError('Cortical observation differs from bound native state')
        if abs(finite(state['mass_kg'],.001,1000,'native mass')-self.identity['target_mass_kg'])>1e-6:raise ValueError('Native mass differs from cortical training basis')
        if descending is not None and not isinstance(descending,dict):raise ValueError('Descending requests must be a dictionary')
        if descending:raise ValueError('Standing policy has no trained descending target port')
        c=self.inner;sensors=c._validate_observation(mechanical_observation)
        sb=c._blocks(sensory_blocks,'sensory block');mb=c._blocks(motor_blocks,'motor block')
        additional={} if additional_sensory_inputs_hz is None else additional_sensory_inputs_hz
        if not isinstance(additional,dict) or set(additional)-set(c.channels):raise ValueError('Unknown additional sensory population')
        rates={k:finite(v,0,1000,'additional sensory rate') for k,v in additional.items()}
        from ihm.assembly.brain import INPUT_BASELINES
        supplied={} if physiology is None else physiology
        if not isinstance(supplied,dict) or set(supplied)-set(INPUT_BASELINES):raise ValueError('Unknown physiology input')
        values={k:finite(v,-1e9,1e9,'physiology '+k) for k,v in dict(INPUT_BASELINES,**supplied).items()}
        pressure=finite(values['mean_arterial_pressure_mmHg'],0,300,'MAP');oxygen=finite(values['oxygen_saturation'],0,1,'oxygen saturation')
        temperature=finite(values['core_temperature_C'],20,45,'temperature');priors=c.brain.data['parameters']['body_transfer_priors']
        availability=min(1.,pressure/priors['map_reference_mmHg'])*min(1.,oxygen/.98)
        temperature_factor=priors['temperature_Q10']**((temperature-37.)/10.)
        vector=self.policy.snapshot_vector(state);saved=self.checkpoint()
        try:
            if sb:
                # Invalidate in-flight afference as well as the current vector.
                reference=self.cortical_state[4:]
                self.cortical_state=tuple(x.clone() for x in (*reference,*reference));vector=self.policy.x0[None]
            with torch.no_grad():
                output,next_state=self.policy.advance(vector,self.cortical_state,ticks=10,sever=self.sever,weights=self.weights,
                    availability=availability,temperature_factor=temperature_factor)
            if not torch.isfinite(output).all() or any(not torch.isfinite(x).all() for x in next_state):raise ValueError('Nonfinite cortical stance state')
            commands=dict(zip(self.policy.muscle_names,(max(.01,min(1.,float(v))) for v in output[0])))
            for name in mb:commands[name]=0.
            self.cortical_state=next_state;self.excitations=commands;self.steps+=1;self.time_s=self.steps*.01
            c.time_s=c.brain.time_s=self.time_s
            region={'node_ids':['implicit-engineering-sensory','implicit-engineering-motor']}
            for field,index in (('potential_mV',0),('activity_hz',1),('adaptation_mV',2)):
                region[field]=[float(next_state[index][0,ports].mean()) for ports in (self.policy.sensory_sites,self.policy.motor_sites)]
            diagnostics={'sensory_blocked':bool(sb),'normalized_state_peak':float(((vector-self.policy.x0)/self.policy.state_scale).abs().max()),
                'max_cortical_correction':float((output[0]-self.policy.u0).abs().max()),'sever':self.sever,
                'additional_sensory_motor_coupling':False,'motor_blocks':sorted(mb)}
            self.last_diagnostics=diagnostics
            return {'schema':'ihm.sensorimotor.v1','time_s':self.time_s,'motor_excitations':dict(commands),'requested_excitations':dict(commands),
                'cortical_commands':dict(commands),'controller':deepcopy(self.controller_metadata),'model_sha256':self.model_sha256,
                'arc_max':dict.fromkeys(('stretch','reciprocal','autogenic','renshaw'),0.),
                'brain':{'time_s':self.time_s,'model_id':'ibm-trained-cortical-stance','source_identity':self.identity,'regional_state':region,
                    'physiology_inputs':values,'physiology_coupling_applied':True,'oxygen_perfusion_availability':availability,'temperature_factor':temperature_factor},
                'sensors':sensors,'brain_sensory_inputs_hz':rates,'activation_owner':'mechanical_plant','cortical_stance':diagnostics,
                'exchange_interval_s':dt,'neural_delay_s':0.,'biological_validation':False,'scope':self.controller_metadata['scope']}
        except Exception:self.restore(saved);raise

    def checkpoint(self):
        return {'schema':'ihm.cortical-stance-state.v1','model_sha256':self.model_sha256,'native_identity':None if self.native is None else self.native.identity,
            'steps':self.steps,'time_s':self.time_s,'cortical_state':[x.tolist() for x in self.cortical_state],
            'excitations':dict(self.excitations),'diagnostics':deepcopy(self.last_diagnostics)}

    def restore(self,checkpoint):
        c=deepcopy(checkpoint)
        if not isinstance(c,dict) or set(c)!=set(self.checkpoint()) or c.get('schema')!='ihm.cortical-stance-state.v1' or c.get('model_sha256')!=self.model_sha256 or c.get('native_identity')!=(None if self.native is None else self.native.identity):raise ValueError('Cortical stance checkpoint identity mismatch')
        time=finite(c['time_s'],0,1e9,'checkpoint time')
        if type(c['steps']) is not int or c['steps']<0 or abs(time-c['steps']*.01)>1e-10:raise ValueError('Cortical checkpoint time differs from integral steps')
        try:states=[torch.tensor(x,dtype=self.policy.dyn.embed.dtype) for x in c['cortical_state']]
        except (TypeError,ValueError,RuntimeError) as exc:raise ValueError('Malformed cortical checkpoint tensors') from exc
        if len(states)!=8 or any(x.shape!=(1,self.policy.dyn.n) or not torch.isfinite(x).all() for x in states):raise ValueError('Malformed cortical state checkpoint')
        if not isinstance(c['excitations'],dict) or set(c['excitations'])!=set(self.policy.muscle_names):raise ValueError('Wrong checkpoint motor width')
        for value in c['excitations'].values():finite(value,0,1,'checkpoint excitation')
        if not isinstance(c['diagnostics'],dict):raise ValueError('Malformed controller diagnostics')
        try:json.dumps(c['diagnostics'],allow_nan=False)
        except (ValueError,TypeError) as exc:raise ValueError('Nonfinite or nonserializable controller diagnostics') from exc
        self.cortical_state=tuple(states);self.steps=c['steps'];self.time_s=c['time_s'];self.excitations=dict(c['excitations']);self.last_diagnostics=c['diagnostics']
        self.inner.time_s=self.inner.brain.time_s=self.time_s
