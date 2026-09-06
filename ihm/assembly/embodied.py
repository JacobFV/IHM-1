"""One-clock articulated, sensorimotor and native physiological execution.

Mechanical/neural checkpoints are exact for their adapters. Native serializer
exactness is not established: an uncertain native command aborts this runtime
instead of claiming a whole-body rollback. The caller retains its native journal.
"""
from copy import deepcopy
import math


def finite(value,label,low=None,high=None):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('Invalid '+label)
    if low is not None and value<low or high is not None and value>high:raise ValueError('Out-of-domain '+label)
    return float(value)


def native_field_metadata(values):
    suffixes=[('compliance_ml_per_mmhg','mL/mmHg'),('concentration_mg_per_dl','mg/dL'),
        ('molarity_mmol_per_l','mmol/L'),('concentration_g_per_l','g/L'),('ml_per_min','mL/min'),
        ('ml_per_s','mL/s'),('l_per_min','L/min'),('l_per_s','L/s'),('pmol_per_min','pmol/min'),
        ('per_min','1/min'),('mmhg','mmHg'),('cmh2o','cmH2O'),('_pa','Pa'),('_ml','mL'),
        ('_mg','mg'),('_g','g'),('_w','W'),('_j','J'),('_c','degC'),('_mv','mV')]
    result={}
    for name in values:
        unit=next((u for suffix,u in suffixes if name.endswith(suffix)),None)
        if unit is None and (name.endswith(('_fraction','_scale','_saturation')) or name.startswith('nervous_') or name in ('arterial_ph','respiratory_exchange_ratio')):unit='1'
        result[name]={'label':name.replace('_',' ').replace('.',' · '),'unit':unit,
            'owner':'external mechanical boundary' if name.startswith('coupling.') else 'BioGears',
            'evidence':'current source-model observation; no implied empirical calibration'}
    return result


class CleanupOwners:
    """Retain partial startup ownership until every owner confirms cleanup."""
    def __init__(self,native,plant):self.native=native;self.plant=plant
    def close(self):
        errors=[]
        for name in ('native','plant'):
            owner=getattr(self,name)
            if owner is None:continue
            try:
                if name=='native':owner.close(graceful=False)
                else:owner.close()
            except BaseException as error:errors.append(f'{name}: {error}')
            else:setattr(self,name,None)
        if errors:raise RuntimeError('Cleanup remains unconfirmed: '+'; '.join(errors))


class EmbodiedRuntime:
    @classmethod
    def from_workspace(cls,root,output,*,environment='supine',state_path=None):
        from pathlib import Path
        import hashlib,json,sys
        from ihm.native.session import SessionConfig
        from ihm.native.coupled_session import CoupledNativeSession
        from .articulated import ArticulatedBodyPlant
        from .sensorimotor import SensorimotorController
        from .body_exchange import NativeTissueExchange
        from .embodied_respiration import EmbodiedRespiration
        from .interactive_scene import _loaded_source
        root=Path(root).resolve();output=Path(output).resolve()
        if not output.is_relative_to(root) or output.exists():raise ValueError('Fresh retained embodied output required')
        state=Path(state_path) if state_path else root/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'
        names=(__name__,'ihm.assembly.articulated','ihm.native.mechanical_stream','ihm.native.coupled_session',
            'ihm.native.session','ihm.assembly.sensorimotor','ihm.assembly.sensorimotor_catalog',
            'ihm.assembly.brain','ihm.assembly.body_exchange','ihm.assembly.body_microstructure',
            'ihm.assembly.respiratory_feedback','ihm.assembly.embodied_respiration','ihm.app.embodied')
        receipts=[_loaded_source(sys.modules[name]) for name in names if name in sys.modules]
        frozen={r['path']:r['bytes'] for r in receipts}
        for name in ('respiration','brain','anatomy','mechanics','microvascular','profile'):
            p=root/f'data/derived/canonical/{name}.json';frozen[p]=p.read_bytes()
        output.mkdir(parents=True)
        hashes={}
        for source,raw in frozen.items():
            relative=source.relative_to(root);destination=output/'inputs'/relative
            destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(raw)
            hashes[str(relative)]=hashlib.sha256(raw).hexdigest()
        native=plant=None
        try:
            native=CoupledNativeSession(SessionConfig(state_path=state,engine_variant='whole_body_integrity_depletion',horizon_s=120),output/'physiology')
            manifest=json.loads((output/'physiology/manifest.json').read_text())
            weight=manifest['patient_identity']['Weight']
            if weight['unit']!='kg':raise ValueError('Expected explicit native initial mass in kg')
            mass=finite(float(weight['value']),'native initial body mass',1,500)
            plant=ArticulatedBodyPlant(root,output/'mechanics',environment=environment,target_mass_kg=mass,
                augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
            neural=SensorimotorController.from_root(root,muscle_catalog=plant.muscle_catalog)
            reference=native.snapshot()
            identity={key:manifest[key] for key in ('library_sha256','executable_sha256','state_sha256')}
            identity['manifest_sha256']=hashlib.sha256((output/'physiology/manifest.json').read_bytes()).hexdigest()
            exchange=NativeTissueExchange.from_workspace(root,reference,identity)
            respiratory_path=root/'data/derived/canonical/respiration.json'
            respiratory=EmbodiedRespiration(json.loads(frozen[respiratory_path]),reference['values']['lung_volume_ml'])
            body=cls(plant,neural,native,exchange,respiratory)
            if any(p.read_bytes()!=raw for p,raw in frozen.items()):raise ValueError('Embodied source changed during initialization; reopen with a stable revision')
            (output/'manifest.json').write_text(json.dumps({'schema':'ihm.embodied-runtime.v1','sources':hashes,
                'loaded_code':{str(r['path'].relative_to(root)):r['loaded_code_sha256'] for r in receipts},
                'source_receipts':{str(r['path'].relative_to(root)):{k:v for k,v in r.items() if k not in ('path','bytes')} for r in receipts},
                'environment':environment,'native_identity':identity,'effective_mechanical_mass_kg':mass,
                'mass_mapping':'Initial native patient mass, including native initial GI contents, uniformly scales source segment inertia; local mass distribution is an engineering prior',
                'native_checkpoint_exact':False,'exchange_dt_s':.02},indent=2)+'\n')
            return body
        except BaseException as error:
            cleanup=CleanupOwners(native,plant)
            try:cleanup.close()
            except BaseException as cleanup_error:
                error.cleanup_owner=cleanup
                error.add_note(str(cleanup_error))
            raise

    def __init__(self,plant,neural,native,exchange,respiratory_load):
        self.plant,self.neural,self.native,self.exchange,self.respiratory_load=plant,neural,native,exchange,respiratory_load
        self.time_s=0.;self.sequence=0;self.failed=False;self.closed=False;self.next_excitation={};self.frame=None
        self.native_state=native.snapshot();self.mechanical_state=plant.snapshot()
        self.cleanup_owner=CleanupOwners(native,plant)
        self.reference_metabolic_w=finite(self.mechanical_state['total_muscle_metabolic_w'],'initial native muscle metabolic power',0)
        if abs(self.native_state['elapsed_s'])>1e-9 or abs(self.mechanical_state['time_s'])>1e-9:
            raise ValueError('Fresh common native and mechanical clocks required')

    def _abort(self):
        self.failed=True;self.cleanup_failures=[]
        for label,close in [('physiology',lambda:self.native.close(graceful=False)),('mechanics',self.plant.close)]:
            try:close()
            except BaseException as error:self.cleanup_failures.append({'owner':label,'error':str(error)})

    def _validate(self,data):
        if not isinstance(data,dict) or set(data)-{'seconds','forces','descending','sensory_blocks','motor_blocks','skin_compression_pa'}:
            raise ValueError('Unknown embodied input')
        dt=finite(data.get('seconds',.02),'embodied interval',.02,.02)
        forces=data.get('forces',[])
        if not isinstance(forces,list) or len(forces)>32:raise ValueError('Expected at most32 force ports')
        for f in forces:
            if not isinstance(f,dict) or set(f)!={'id','force_n','point_m'} or not isinstance(f['id'],str):raise ValueError('Invalid force port')
            for key in ('force_n','point_m'):
                if not isinstance(f[key],(list,tuple)) or len(f[key])!=3:raise ValueError('Expected spatial force/point vector')
                for value in f[key]:finite(value,key,-1000 if key=='force_n' else -5,1000 if key=='force_n' else 5)
        if 'skin_compression_pa' in data:finite(data['skin_compression_pa'],'whole-skin pressure',0,5000)
        horizon=getattr(getattr(self.native,'config',None),'horizon_s',None)
        if horizon is not None and self.time_s+dt>horizon+1e-9:raise ValueError('Native horizon reached; no owners advanced')
        return dt,deepcopy(forces)

    def step(self,data):
        if self.failed or self.closed:raise RuntimeError('Embodied runtime is no longer advancing')
        dt,forces=self._validate(data)
        p_checkpoint=n_checkpoint=None;native_touched=False
        try:
            p_checkpoint=self.plant.checkpoint();n_checkpoint=self.neural.checkpoint()
            v=self.native_state['values']
            physiology={'mean_arterial_pressure_mmHg':v['mean_arterial_pressure_mmhg'],
                'oxygen_saturation':v['oxygen_saturation'],'core_temperature_C':v['core_temperature_c']}
            neural=self.neural.step(dt,self.mechanical_state,descending=data.get('descending',{}),
                sensory_blocks=data.get('sensory_blocks',()),motor_blocks=data.get('motor_blocks',()),physiology=physiology)
            # Neural output belongs to the next exchange interval. Current
            # motor blocks still suppress already-delivered excitations now.
            actuation={k:v for k,v in self.next_excitation.items() if k not in data.get('motor_blocks',())}
            mechanical=self.plant.advance(dt,forces=forces,actuation=actuation)
            end=self.time_s+dt
            if abs(mechanical['time_s']-end)>1e-8:raise RuntimeError('Mechanical exchange clock diverged')
            load=self.respiratory_load.project_load(forces,mechanical['entities'],v['lung_volume_ml'])
            # Segment resultants cannot determine skin/organ strain work. Keep
            # these visible until a resolved traction mapping is available.
            load['unresolved_contact_wrenches']=deepcopy(mechanical.get('body_environment',{}).get('canonical_wrenches',[]))
            pressure=finite(load['external_pressure_pa'],'respiratory load',-5000,5000)
            work=finite(mechanical['positive_muscle_work_j'],'positive muscle work',0)
            energy=finite(mechanical['muscle_metabolic_energy_j'],'native muscle metabolic energy')
            previous=finite(self.mechanical_state['muscle_metabolic_energy_j'],'previous native muscle metabolic energy')
            metabolic_w=(energy-previous)/dt
            incremental_w=metabolic_w-self.reference_metabolic_w
            if incremental_w<0:raise ValueError('Muscle metabolic demand fell below initial reference; native exercise port cannot represent signed decrement')
            max_work=finite(v['maximum_work_rate_w'],'native maximum work rate',1e-12)
            intensity=incremental_w/max_work
            if intensity>.5:raise ValueError('Computed metabolic demand exceeds the supported native exercise range')
            native_touched=True
            self.native.respiratory_load(pressure)
            self.native.exercise(intensity)
            if 'skin_compression_pa' in data:self.native.skin_compression(data['skin_compression_pa'])
            native=self.native.step(dt)
            native['signal_metadata']=native_field_metadata(native['values'])
            if abs(native['elapsed_s']-end)>1e-8:raise RuntimeError('Native exchange clock diverged')
            tissue=self.exchange.observe(native)
            geometry=self.respiratory_load.geometry(native['values']['lung_volume_ml'],mechanical['entities'],end)
            self.next_excitation=deepcopy(neural['motor_excitations'])
            self.native_state=native;self.mechanical_state=mechanical;self.time_s=end;self.sequence+=1
            self.frame={'schema':'ihm.embodied-frame.v1','time_s':end,'sequence':self.sequence,
                'entities':geometry['entities'],'skin_field':geometry['skin_field'],'respiration':geometry,'mechanics':mechanical,'neural':neural,'physiology':native,
                'tissue_exchange':tissue,'respiratory_load':load,
                'coupling':{'exchange_interval_s':dt,'motor_exchange_latency_s':dt,
                    'positive_muscle_work_j':work,'native_extra_metabolic_demand_w':incremental_w,
                    'muscle_metabolic_reference_w':self.reference_metabolic_w,'interval_muscle_metabolic_w':metabolic_w,
                    'native_exercise_intensity':intensity,
                    'metabolic_law':'Native Umberger muscle energy increment / interval minus fixed initial muscle reference. BioGears owns basal metabolism and ramps its exercise setpoint; this is not instantaneous energy equality. Negative increments rejected.',
                    'storage_owners':{'articulation_muscle':'native mechanical plant','neural':'pinned IBM plus declared decoder/reflexes',
                        'blood_gas_nutrients_heat':'BioGears','tissue_views':'native-owned compartments, no duplicate storage'},
                    'rollback':'An uncertain native commit terminates this runtime; no serializer-exactness claim'}}
            result=deepcopy(self.frame)
            if hasattr(self.plant,'release'):self.plant.release(p_checkpoint)
            return result
        except BaseException:
            if native_touched:
                self._abort()
            else:
                try:
                    if p_checkpoint is None or n_checkpoint is None:raise RuntimeError('Checkpoint acquisition failed')
                    self.plant.restore(p_checkpoint);self.neural.restore(n_checkpoint)
                    if hasattr(self.plant,'release'):self.plant.release(p_checkpoint)
                except BaseException:self._abort()
            raise

    def snapshot(self):
        if self.frame:return deepcopy(self.frame)
        native=deepcopy(self.native_state);native['signal_metadata']=native_field_metadata(native['values'])
        geometry=self.respiratory_load.geometry(native['values']['lung_volume_ml'],self.mechanical_state['entities'],0.)
        return {'schema':'ihm.embodied-frame.v1','time_s':0.,'sequence':0,
            'entities':geometry['entities'],'skin_field':geometry['skin_field'],'respiration':geometry,
            'mechanics':deepcopy(self.mechanical_state),'physiology':native,
            'tissue_exchange':self.exchange.observe(self.native_state)}

    def close(self):
        if self.closed:return
        self.cleanup_owner.close()
        self.closed=True
