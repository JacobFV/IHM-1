"""One-clock articulated, sensorimotor and native physiological execution.

Mechanical/neural checkpoints are exact for their adapters. Native serializer
exactness is not established: an uncertain native command aborts this runtime
instead of claiming a whole-body rollback. The caller retains its native journal.
"""
from ihm.brain.active_source import DEFAULT_SOURCE,resolve_source
from copy import deepcopy
from .snapshot_data import clone_snapshot_data
import math
from .intake_schedule import IntakeSchedule,IntakeEvent
from ihm.native.session import Meal


def finite(value,label,low=None,high=None):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('Invalid '+label)
    if low is not None and value<low or high is not None and value>high:raise ValueError('Out-of-domain '+label)
    return float(value)


# The mechanical plant reports muscle metabolic power as a 50 Hz finite
# difference of a cumulative Umberger/Uchida ledger. That signal fluctuates by
# tens of watts between adjacent 20 ms intervals on a body holding a posture,
# while the native muscle compartment's discretionary energy above its
# obligatory amino-acid draw is roughly 8.7 W. Exchanging the instantaneous
# difference therefore rejects a physically ordinary interval. This first-order
# lag carries the same increments with a substrate-delivery time constant.
# It is an explicit engineering choice, not an identified substrate kinetic:
# it delays energy, never creates or destroys it, and its mean is the mean of
# the raw increments.
METABOLIC_EXCHANGE_TAU_S=2.


def measure_resting_metabolic_reference(plant,*,settling_s=.2,window_s=.6,dt_s=.02):
    """Measure a fixed zero-excitation counterfactual, then restore the plant.

    Native default activation is an initial condition, not the resting metabolic
    budget already represented by physiology. No task/control outputs enter this
    reference and it is never adapted during a running task. The short window is
    an engineering baseline, not a claim of equilibrated human resting metabolism.
    """
    owner=getattr(plant,'native',plant)
    initial=owner.snapshot()
    if abs(initial['time_s'])>1e-9:raise ValueError('Fresh mechanical reference required')
    for value,label in ((settling_s,'settling duration'),(window_s,'reference window'),(dt_s,'reference step')):
        finite(value,label,.001,10.)
        if abs(value/dt_s-round(value/dt_s))>1e-8:raise ValueError('Reference durations must be whole steps')
    commands={key:0. for key in initial['muscles']}
    checkpoint=owner.checkpoint()
    try:
        start=initial
        for _ in range(round(settling_s/dt_s)):start=owner.advance(dt_s,actuation=commands)
        end=start
        for _ in range(round(window_s/dt_s)):end=owner.advance(dt_s,actuation=commands)
        m=(end['muscle_metabolic_energy_j']-start['muscle_metabolic_energy_j'])/window_s
        w=(end['signed_active_fiber_work_j']-start['signed_active_fiber_work_j'])/window_s
        h=(end['muscle_heat_energy_j']-start['muscle_heat_energy_j'])/window_s
        for value in (m,w,h):finite(value,'resting reference power')
        if abs(m-w-h)>1e-9*(1+abs(m)+abs(w)+abs(h)):raise ValueError('Resting reference energy ledger mismatch')
        return {'M0_w':m,'W0_w':w,'H0_w':m-w,
            'basis':'Fixed zero-excitation native mechanical counterfactual; integrated after declared settling, exact initial checkpoint restored; no task-policy adaptation or budget clamp. Body may move under gravity: this is not an equilibrated resting or standing measurement.',
            'environment':initial.get('environment',initial.get('body_environment',{}).get('kind')),
            'settling_s':settling_s,'window_s':window_s,'dt_s':dt_s,
            'zero_excitation_muscle_ids':sorted(commands),
            'measurement_energy_j':{'chemical':end['muscle_metabolic_energy_j']-start['muscle_metabolic_energy_j'],
                'signed_work':end['signed_active_fiber_work_j']-start['signed_active_fiber_work_j'],
                'heat':end['muscle_heat_energy_j']-start['muscle_heat_energy_j']},
            'initial_instantaneous_reference':deepcopy(initial['metabolic_reference'])}
    finally:
        owner.restore(checkpoint)
        if hasattr(owner,'release'):owner.release(checkpoint)


def native_field_metadata(values):
    suffixes=[('compliance_ml_per_mmhg','mL/mmHg'),('concentration_mg_per_dl','mg/dL'),
        ('molarity_mmol_per_l','mmol/L'),('concentration_g_per_l','g/L'),('ml_per_min','mL/min'),
        ('ml_per_s','mL/s'),('l_per_min','L/min'),('l_per_s','L/s'),('pmol_per_min','pmol/min'),
        ('per_min','1/min'),('mmhg','mmHg'),('cmh2o','cmH2O'),('_pa','Pa'),('_ml','mL'),
        ('_kcal','kcal'),('_mg','mg'),('_g','g'),('_w','W'),('_j','J'),('_c','degC'),('_mv','mV'),('_v','V')]
    result={}
    for name in values:
        unit=next((u for suffix,u in suffixes if name.endswith(suffix)),None)
        if unit is None and (name.endswith(('_fraction','_scale','_saturation')) or name.startswith('nervous_') or name in ('arterial_ph','respiratory_exchange_ratio')):unit='1'
        result[name]={'label':name.replace('_',' ').replace('.',' · '),'unit':unit,
            'owner':'external mechanical boundary' if name.startswith('coupling.') else 'BioGears',
            'evidence':'current source-model observation; no implied empirical calibration'}
    return result


def bind_cutaneous(root,contacts,configuration,*,source_pin=DEFAULT_SOURCE):
    """Bind explicit cortical recruitment priors to exact native material sites."""
    from .cutaneous_feedback import CutaneousFeedback
    if not isinstance(configuration,dict) or set(configuration)!={'regions','recruitment_hz_per_response','reference_temperature_C'}:
        raise ValueError('Explicit cutaneous regions, recruitment and thermal reference required')
    regions=configuration['regions']
    if not isinstance(regions,dict) or not 1<=len(regions)<=64:raise ValueError('Expected 1 to 64 skin sensory regions')
    if not isinstance(contacts,list) or len(contacts)!=len(regions) or {p['id'] for p in contacts}!=set(regions):
        raise ValueError('Native skin sensor identities differ from requested cortical mapping')
    sites=[]
    for point in contacts:
        sites.append({'id':point['id'],'position_m':point['point_m'],'normal':point['normal'],
            'contact_area_m2':point['contact_area_m2'],'mechanical_input':'native_indentation',
            'material_identity':point['material_identity'],'indentation_basis':point['indentation_basis'],
            'area_basis':point['area_basis'],'sensory_region':regions[point['id']],
            'reference_temperature_C':configuration['reference_temperature_C'],
            'support_basis':'Retained native skin quadrature, rigid canonical registration; cortical mapping and recruitment are explicit engineering priors'})
    return CutaneousFeedback(root,sites=sites,recruitment_hz_per_response=configuration['recruitment_hz_per_response'],source_pin=source_pin)


def bind_intake_mass(plant,initial_intake,canonical_raw):
    """Bind the retained neutral canonical stomach centroid to its native torso."""
    import hashlib,json
    import numpy as np
    from .intake_mass import IntakeMassBridge
    if plant.snapshot()['time_s']!=0:raise ValueError('Fresh neutral mechanical reference required for intake binding')
    if (plant.output/'canonical_mechanics.json').read_bytes()!=canonical_raw:
        raise ValueError('Mechanical canonical reference differs from frozen intake source')
    canonical=json.loads(canonical_raw)
    candidates=[e for e in canonical['entities'] if e.get('name')=='stomach']
    if len(candidates)!=1:raise ValueError('Exactly one retained canonical stomach required')
    stomach=candidates[0];registration=plant.registration
    if registration.rows[stomach['id']]['body']!='torso':raise ValueError('Canonical stomach must be bound to native torso')
    transform=np.array(plant.native.snapshot()['bodies']['torso']['transform_ground'])
    station=(np.linalg.inv(transform)@np.linalg.inv(registration.global_map)@np.r_[stomach['centroid_m'],1.])[:3].tolist()
    registration_raw=(plant.output/'registration.json').read_bytes()
    registration_hash=hashlib.sha256(registration_raw).hexdigest()
    binding={'canonical_entity_id':stomach['id'],'body':'torso','station_source_m':station,
        'canonical_sha256':hashlib.sha256(canonical_raw).hexdigest(),'registration_sha256':registration_hash,
        'registration_basis':deepcopy(registration.rows[stomach['id']]),
        'incoming_velocity_basis':'co_moving_at_ingestion_assumption',
        'scope':'Neutral canonical stomach centroid on inferred rigid torso support; no internal organ deformation or measured swallowing momentum'}
    bridge=IntakeMassBridge(plant,initial_intake,body='torso',station_m=station,
        registration_identity=registration_hash,incoming_velocity_basis=binding['incoming_velocity_basis'])
    return bridge,binding


def _prepare_brain_candidate(root,source_pin):
    """Verify/load one candidate before native startup and freeze its exact bytes."""
    if source_pin is None:return {},{}
    import hashlib,json,sys
    from ihm.brain.candidate import SourcePin,verify_pin
    from ihm.brain.ibm_backend import IBMBackend
    if not isinstance(source_pin,SourcePin):raise ValueError('Explicit SourcePin required')
    if not source_pin.artifact_dir.is_relative_to(root):
        raise ValueError('Body source candidate must be inside the configured workspace for retained relative receipts')
    manifest,snapshots=verify_pin(source_pin)
    pin_path=source_pin.artifact_dir/'source_pin.json'
    pin_raw=pin_path.read_bytes()
    if json.loads(pin_raw)!=source_pin.to_dict():raise ValueError('Retained source pin file differs from selected pin')
    manifest_path=source_pin.artifact_dir/'manifest.json'
    manifest_raw=manifest_path.read_bytes()
    if hashlib.sha256(manifest_raw).hexdigest()!=source_pin.manifest_sha256:
        raise ValueError('Candidate manifest changed during preparation')
    # SnapshotLoader rejects an already-loaded different package here, before
    # output creation, source receipt capture or either native owner is opened.
    IBMBackend(source_pin=source_pin)
    frozen={manifest_path:manifest_raw,pin_path:pin_raw}
    frozen.update({source_pin.artifact_dir/'source'/name:raw for name,raw in snapshots.items()})
    loaded={}
    for name,module in tuple(sys.modules.items()):
        if name=='ibm' or name.startswith('ibm.'):
            if getattr(module,'__ihm_source_identity__',None)!=source_pin.package_sha256:
                raise RuntimeError('Loaded IBM module differs from selected source pin')
            from pathlib import Path
            path=Path(module.__file__).resolve()
            if path not in frozen:raise ValueError('Loaded IBM module lacks a captured source receipt')
            loaded[name]={'path':str(path.relative_to(root)),
                          'source_sha256':hashlib.sha256(frozen[path]).hexdigest(),
                          'package_sha256':source_pin.package_sha256}
    return frozen,loaded


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


def _prepare_mechanical_registration(root,relative):
    """Bounded source-only opt-in preflight before either native owner starts."""
    import xml.etree.ElementTree as ET
    from pathlib import Path
    import math,hashlib,json
    root=Path(root).resolve();frozen={};total=0
    def read(name,digest=None):
        nonlocal total
        if not isinstance(name,str) or Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError('Mechanical registration paths must be workspace-relative')
        path=root/name
        if not path.is_relative_to(root) or any(p.is_symlink() for p in (path,*path.parents)):
            raise ValueError('Mechanical registration may not use symlink redirects')
        if not path.is_file() or path.stat().st_size>32*1024*1024:
            raise ValueError('Mechanical registration source must be a bounded ordinary file')
        raw=path.read_bytes()
        if digest is not None and hashlib.sha256(raw).hexdigest()!=digest:raise ValueError('Mechanical registration source hash differs')
        if path not in frozen:total+=len(raw)
        if total>64*1024*1024:raise ValueError('Mechanical registration exceeds 64 MiB source budget')
        frozen[path]=raw;return raw
    manifest=json.loads(read(relative))
    if manifest.get('schema') not in ('ihm.upperbody-registration.v1','ihm.lumbar-muscle-variant.v1'):
        raise ValueError('Unsupported mechanical registration schema')
    model=ET.fromstring(read(manifest['model_path'],manifest['model_sha256']))
    rows=json.loads(read(manifest['catalog_path'],manifest['catalog_sha256']))
    for path,digest in manifest['sources'].items():read(path,digest)
    if 'insert_path' in manifest:read(manifest['insert_path'],manifest['insert_sha256'])
    muscles=[m for m in model.findall('.//ForceSet/objects/*') if 'Muscle' in m.tag]
    if (not isinstance(rows,list) or len(rows)!=manifest['muscle_count']
        or len({r['id'] for r in rows})!=len(rows) or {m.get('name') for m in muscles}!={r['id'] for r in rows}):
        raise ValueError('Mechanical catalog and model muscle identities differ')
    by_id={m.get('name'):m for m in muscles}
    for row in rows:
        if manifest['sources'].get(row['source_path'])!=row['source_sha256']:
            raise ValueError('Mechanical catalog source owner differs')
        for field,tag in (('max_isometric_force_n','max_isometric_force'),('optimal_fiber_length_m','optimal_fiber_length')):
            value=row[field]
            if isinstance(value,bool) or not math.isfinite(value) or value<=0 or not math.isclose(value,float(by_id[row['id']].findtext(tag)),rel_tol=1e-12):
                raise ValueError('Mechanical catalog normalization differs from model')
    return frozen,manifest,rows


class EmbodiedRuntime:
    @classmethod
    def from_workspace(cls,root,output,*,environment='supine',state_path=None,surface_contact_manifest=None,cutaneous_configuration=None,bed_material=None,regional_skin=False,source_pin=DEFAULT_SOURCE,intake_mass=False,augmented_registration=None,native_afferent_allocation=None,environment_selection=None,controller=None,initial_pose=None):
        from .controller_selection import resolve_controller
        controller=resolve_controller(controller)
        if controller['kind']!='regional' and native_afferent_allocation is not None:
            raise ValueError('Selected controller does not yet support native afferent allocation')
        if type(regional_skin) is not bool:raise ValueError('regional_skin must be a bool')
        if type(intake_mass) is not bool:raise ValueError('intake_mass must be a bool')
        from .environment_dynamics import EnvironmentDynamics, resolve_selection
        from . import rigid_contact, cloth_contact, cloth_stretch, cloth_cover, world_frame, world_exchange, control_exchange
        environment_selection, environment_options = resolve_selection(root,environment,environment_selection)
        if environment_options:
            surface_contact_manifest=environment_options.get('surface_contact_manifest',surface_contact_manifest)
            bed_material=environment_options.get('bed_material',bed_material)
        from pathlib import Path
        import hashlib,json,sys
        from ihm.native.session import SessionConfig
        from ihm.native.coupled_session import SignedCoupledNativeSession
        native_session_type=SignedCoupledNativeSession
        if regional_skin:
            from ihm.native.regional_session import RegionalSignedNativeSession
            native_session_type=RegionalSignedNativeSession
        if native_afferent_allocation is not None:
            from ihm.native.afferent_session import AfferentSignedNativeSession,AfferentRegionalNativeSession
            native_session_type=AfferentRegionalNativeSession if regional_skin else AfferentSignedNativeSession
        from .selective_projection import SelectiveProjectionPlant
        if intake_mass:
            from .intake_mass import IntakeMassBridge
        from .sensorimotor import SensorimotorController
        from .body_exchange import NativeTissueExchange
        from .embodied_respiration import EmbodiedRespiration
        from .interactive_scene import _loaded_source
        from .cutaneous_feedback import CutaneousFeedback
        root=Path(root).resolve();output=Path(output).resolve()
        if not output.is_relative_to(root) or output.exists():raise ValueError('Fresh retained embodied output required')
        from .controller_selection import STANCE_KINDS,CORTICAL_STANCE_KINDS
        if controller['kind'] in STANCE_KINDS:
            if environment!='upright' or intake_mass or surface_contact_manifest is not None:
                raise ValueError('Stance controllers require upright fixed-mass native floor mechanics')
            stance_path='data/models/engineering_stance_v1/registration.json'
            if augmented_registration is not None or initial_pose is not None:
                raise ValueError('Stance controller owns its identified model and initial pose')
            if not (root/stance_path).is_file():raise ValueError('Validated engineered stance bundle is unavailable')
            augmented_registration=stance_path
            initial_pose=json.loads((root/stance_path).read_text())['initial_pose']
        source_pin=None if controller['kind'] not in ('regional','engineering_stance') else resolve_source(root,source_pin)
        # Resolve optional controller code before freezing source receipts. A
        # long-lived workbench must reject stale imported code after an edit,
        # rather than retain new bytes for an older loaded implementation.
        controller_modules=()
        if controller['kind']!='regional':
            import importlib
            controller_modules=(('ihm.native.stance_lqr','ihm.native.stance_lqr_controller') if controller['kind']=='engineering_stance'
                                else ('ihm.assembly.ibm_controller',))
            if controller['kind'] in CORTICAL_STANCE_KINDS:
                controller_modules+=('ihm.native.cortical_stance','ihm.native.cortical_stance_controller')
            elif controller['kind']=='implicit_cortical_ankle':
                controller_modules+=('ihm.native.cortical_motor','ihm.native.cortical_motor_controller')
            elif controller['kind']=='implicit_ankle_primitive':
                controller_modules+=('ihm.native.motor_learning','ihm.native.motor_learning_controller')
            for name in controller_modules:importlib.import_module(name)
        mechanical_frozen,mechanical_manifest,mechanical_catalog=({},None,None) if augmented_registration is None else _prepare_mechanical_registration(root,augmented_registration)
        candidate_frozen,candidate_loaded=_prepare_brain_candidate(root,source_pin)
        if native_afferent_allocation is not None:
            from .native_afferents import NativeAfferentBridge
            from types import SimpleNamespace
            # Static allocation validation occurs before either native owner starts.
            NativeAfferentBridge(SimpleNamespace(data=json.loads((root/'data/derived/canonical/brain.json').read_bytes()),source_identity={}),
                native_identity={k:'0'*64 for k in ('library_sha256','executable_sha256','state_sha256','manifest_sha256')},
                source_sha256='0'*64,origin_s=0.,allocation=native_afferent_allocation)
        sensor_indices=[]
        if cutaneous_configuration is not None:
            if surface_contact_manifest is None:raise ValueError('Cutaneous binding requires native surface contact')
            regions=cutaneous_configuration.get('regions') if isinstance(cutaneous_configuration,dict) else None
            if not isinstance(regions,dict) or not 1<=len(regions)<=64:raise ValueError('Expected explicit bounded skin sensor mapping')
            for key in regions:
                if not isinstance(key,str) or not key.startswith('skin-contact-') or not key[13:].isdigit():raise ValueError('Invalid native skin sensor ID')
                index=int(key[13:])
                if key!='skin-contact-'+str(index):raise ValueError('Noncanonical native skin sensor ID')
                sensor_indices.append(index)
        reference_path=root/'data/derived/systemic/exertion_v3/exercise/native/manifest.json'
        reference_raw=reference_path.read_bytes();reference_manifest=json.loads(reference_raw)
        base_variant=reference_manifest['configuration']['engine_variant']
        if base_variant!='whole_body_integrity_evaporation_humidity':raise ValueError('Expected retained final thermal-corrected research variant')
        engine_variant='whole_body_integrity_regional_skin_graph_v2' if regional_skin else 'whole_body_integrity_gi_absorption'
        state=Path(state_path) if state_path else Path(reference_manifest['configuration']['state_path'])
        if state_path is None and hashlib.sha256(state.read_bytes()).hexdigest()!=reference_manifest['state_sha256']:
            raise ValueError('Paired native initial state changed')
        names=(__name__,'ihm.assembly.environment_dynamics','ihm.assembly.rigid_contact','ihm.assembly.cloth_contact','ihm.assembly.cloth_stretch','ihm.assembly.cloth_cover','ihm.assembly.world_frame','ihm.assembly.controller_selection','ihm.assembly.world_exchange','ihm.assembly.articulated','ihm.native.mechanical_stream','ihm.native.coupled_session',
            'ihm.native.session','ihm.assembly.sensorimotor','ihm.assembly.sensorimotor_catalog',
            'ihm.assembly.brain','ihm.assembly.body_exchange','ihm.assembly.regional_exchange','ihm.assembly.body_microstructure',
            'ihm.assembly.cutaneous_feedback','ihm.brain.causal','ihm.brain.ibm_backend','ihm.brain.port_mapping',
            'ihm.brain.active_source','ihm.assembly.respiratory_feedback','ihm.assembly.embodied_respiration','ihm.assembly.intake_schedule','ihm.app.embodied')
        if regional_skin:names+=('ihm.native.regional_session',)
        if intake_mass:names+=('ihm.assembly.intake_mass','ihm.native.instance_mass',)
        if native_afferent_allocation is not None:names+=('ihm.assembly.native_afferents','ihm.native.afferent_session',)
        if source_pin is not None:names+=('ihm.brain.candidate','ihm.brain.source_loader',)
        names+=controller_modules+('ihm.assembly.snapshot_data','ihm.assembly.control_exchange','ihm.assembly.selective_projection','ihm.assembly.surface_binding','ihm.assembly.continuous_surface_binding')
        receipts=[_loaded_source(sys.modules[name]) for name in names if name in sys.modules]
        frozen={r['path']:r['bytes'] for r in receipts}
        frozen.update(candidate_frozen)
        frozen.update(mechanical_frozen)
        frozen[reference_path]=reference_raw
        variant_dir=root/'data/runtime/physiology/variants';current=engine_variant;seen=set()
        while True:
            if current in seen:raise ValueError('Cyclic native variant lineage')
            seen.add(current);path=variant_dir/current/'manifest.json';raw=path.read_bytes();entry=json.loads(raw);frozen[path]=raw
            if hashlib.sha256((path.parent/'libbiogears.so.8.0.0').read_bytes()).hexdigest()!=entry['library_sha256']:raise ValueError('Native lineage library changed')
            if current==base_variant:
                if entry['library_sha256']!=reference_manifest['library_sha256']:raise ValueError('Paired native base library changed')
                break
            parent=entry['parent_variant'];parent_path=variant_dir/parent/'manifest.json';parent_raw=parent_path.read_bytes()
            if hashlib.sha256(parent_raw).hexdigest()!=entry['parent_manifest_sha256'] or json.loads(parent_raw)['library_sha256']!=entry['parent_library_sha256']:raise ValueError('Native lineage manifest changed')
            current=parent
        for name in ('respiration','brain','anatomy','mechanics','microvascular','profile'):
            p=root/f'data/derived/canonical/{name}.json';frozen[p]=p.read_bytes()
        surface_asset=root/'data/derived/canonical/continuous_surface_binding.json.gz'
        frozen[surface_asset]=surface_asset.read_bytes()
        output.mkdir(parents=True)
        hashes={}
        for source,raw in frozen.items():
            relative=source.relative_to(root);destination=output/'inputs'/relative
            destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(raw)
            hashes[str(relative)]=hashlib.sha256(raw).hexdigest()
        native=plant=None
        try:
            native=native_session_type(SessionConfig(state_path=state,engine_variant=engine_variant,horizon_s=120),output/'physiology')
            manifest=json.loads((output/'physiology/manifest.json').read_text())
            if manifest['library_sha256']!=json.loads(frozen[variant_dir/engine_variant/'manifest.json'])['library_sha256']:raise ValueError('Signed native library changed')
            weight=manifest['patient_identity']['Weight']
            if weight['unit']!='kg':raise ValueError('Expected explicit native initial mass in kg')
            mass=finite(float(weight['value']),'native initial body mass',1,500)
            plant=SelectiveProjectionPlant(root,output/'mechanics',environment=environment,target_mass_kg=mass,
                augmented_registration=augmented_registration or 'data/derived/mechanics/whole_body_arm26_v2/registration.json',
                surface_contact_manifest=surface_contact_manifest,surface_sensor_indices=sensor_indices,bed_material=bed_material,
                instance_mass_variant='data/runtime/opensim/variants/instance_mass_v1' if intake_mass else None,
                **({'initial_pose':initial_pose} if initial_pose is not None else {}))
            if mechanical_catalog is not None and plant.muscle_catalog!=mechanical_catalog:
                raise ValueError('Native plant catalog differs from preflight mechanical identity')
            metabolic_reference=measure_resting_metabolic_reference(plant)
            (output/'metabolic_reference.json').write_text(json.dumps(metabolic_reference,indent=2,allow_nan=False)+'\n')
            if controller['kind']!='regional':
                if controller['kind']=='engineering_stance':
                    from ihm.native.stance_lqr_controller import EngineeredLQRStanceController
                    neural=EngineeredLQRStanceController.from_root(root,muscle_catalog=plant.muscle_catalog,source_pin=source_pin)
                    neural.bind_native(plant.native)
                elif controller['kind'] in CORTICAL_STANCE_KINDS:
                    from ihm.native.cortical_stance_controller import IBMCorticalStanceController
                    neural=IBMCorticalStanceController.from_root(root,muscle_catalog=plant.muscle_catalog,
                        kind=controller['kind'],sever=controller['sever'])
                    neural.bind_native(plant.native)
                elif controller['kind']=='implicit_cortical_ankle':
                    from ihm.native.cortical_motor_controller import IBMCorticalAnkleController
                    neural=IBMCorticalAnkleController.from_root(root,muscle_catalog=plant.muscle_catalog,
                        target_rad=controller['target_rad'],sever=controller['sever'],no_cord=controller['no_cord'])
                elif controller['kind']=='implicit_ankle_primitive':
                    from ihm.native.motor_learning_controller import IBMAnklePrimitiveController
                    neural=IBMAnklePrimitiveController.from_root(root,muscle_catalog=plant.muscle_catalog,
                        target_rad=controller['target_rad'],sever=controller['sever'])
                else:
                    from .ibm_controller import IBMImplicitController
                    from .controller_selection import RAW_IMPLICIT_KERNELS
                    kernel=RAW_IMPLICIT_KERNELS[controller['kind']]
                    if kernel is not None and not (root/kernel).is_file():
                        raise ValueError('Retained IBM kernel for this controller kind is unavailable')
                    neural=IBMImplicitController.from_root(root,muscle_catalog=plant.muscle_catalog,
                        sever=controller['sever'],no_cord=controller['no_cord'],
                        kind=controller['kind'],
                        checkpoint_path=None if kernel is None else (root/kernel).resolve())
                neural.retain_sources(output)
            else:
                neural=SensorimotorController.from_root(root,muscle_catalog=plant.muscle_catalog,source_pin=source_pin)
            reference=native.snapshot()
            intake_bridge=intake_binding=None
            if intake_mass:
                intake_bridge,intake_binding=bind_intake_mass(plant,reference['intake'],frozen[root/'data/derived/canonical/mechanics.json'])
                (output/'intake_mass_binding.json').write_text(json.dumps(intake_binding,indent=2,allow_nan=False)+'\n')
            identity={key:manifest[key] for key in ('library_sha256','executable_sha256','state_sha256')}
            identity['manifest_sha256']=hashlib.sha256((output/'physiology/manifest.json').read_bytes()).hexdigest()
            afferents=None
            if native_afferent_allocation is not None:
                observer_sha=manifest['adapter_sha256']['native_nervous_afferents.h']
                afferents=NativeAfferentBridge(neural.brain,native_identity=identity,source_sha256=observer_sha,
                    origin_s=reference['origin_s'],allocation=native_afferent_allocation)
            exchange=NativeTissueExchange.from_workspace(root,reference,identity)
            respiratory_path=root/'data/derived/canonical/respiration.json'
            respiratory=EmbodiedRespiration(json.loads(frozen[respiratory_path]),reference['values']['lung_volume_ml'])
            cutaneous=None if cutaneous_configuration is None else bind_cutaneous(root,plant.snapshot().get('cutaneous_contacts'),cutaneous_configuration,source_pin=source_pin)
            environment_owner=None if environment_selection is None else EnvironmentDynamics(root,environment,environment_selection,plant.registration)
            if environment_owner is not None:
                for relative,digest in environment_owner.sources.items():
                    raw=(root/relative).read_bytes()
                    if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError('Environment asset changed during startup')
                    destination=output/'inputs'/relative
                    destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(raw)
                    hashes[relative]=digest
            body=cls(plant,neural,native,exchange,respiratory,cutaneous=cutaneous,afferents=afferents,intake_mass_bridge=intake_bridge,intake_mass_binding=intake_binding,reference_identity=hashlib.sha256((output/'mechanics/native/execution.json').read_bytes()).hexdigest(),environment_dynamics=environment_owner,metabolic_reference=metabolic_reference)
            if any(p.read_bytes()!=raw for p,raw in frozen.items()):raise ValueError('Embodied source changed during initialization; reopen with a stable revision')
            (output/'manifest.json').write_text(json.dumps({'schema':'ihm.embodied-runtime.v1','sources':hashes,
                'loaded_code':{str(r['path'].relative_to(root)):r['loaded_code_sha256'] for r in receipts},
                'source_receipts':{str(r['path'].relative_to(root)):{k:v for k,v in r.items() if k not in ('path','bytes')} for r in receipts},
                'mechanical_registration_override':None if mechanical_manifest is None else {'path':augmented_registration,'sha256':hashlib.sha256(mechanical_frozen[root/augmented_registration]).hexdigest(),'model_sha256':mechanical_manifest['model_sha256'],'catalog_sha256':mechanical_manifest['catalog_sha256']},
                'controller_selection':controller,'initial_pose':initial_pose,
                'brain_source_pin':None if source_pin is None else source_pin.to_dict(),'brain_candidate_loaded_modules':candidate_loaded,
                'native_afferent_allocation':native_afferent_allocation,'native_afferent_model_sha256':None if afferents is None else afferents.model_sha256,
                'environment':environment,'environment_selection':environment_selection,'environment_sources':None if environment_owner is None else environment_owner.sources,'regional_skin':regional_skin,'intake_mass':intake_mass,'intake_mass_binding':intake_binding,
                'intake_mass_initial_bridge':None if intake_bridge is None else intake_bridge.snapshot(),
                'intake_mass_variant':None if intake_bridge is None else plant.native.instance_mass_variant,
                'intake_mass_binding_sha256':None if intake_binding is None else hashlib.sha256((output/'intake_mass_binding.json').read_bytes()).hexdigest(),'cutaneous_materialization':None if cutaneous is None else cutaneous.audit,'native_identity':identity,'effective_mechanical_mass_kg':mass,
                'physiology_scope':'Paired retained thermal-corrected research initial state/library; known long-run glucose and acid-base failures remain unresolved',
                'metabolic_reference':metabolic_reference,
                'metabolic_reference_sha256':hashlib.sha256((output/'metabolic_reference.json').read_bytes()).hexdigest(),
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

    def __init__(self,plant,neural,native,exchange,respiratory_load,*,reference_identity=None,cutaneous=None,intake_mass_bridge=None,intake_mass_binding=None,afferents=None,environment_dynamics=None,metabolic_reference=None):
        self.plant,self.neural,self.native,self.exchange,self.respiratory_load=plant,neural,native,exchange,respiratory_load
        self.environment_dynamics=environment_dynamics
        self.time_s=0.;self.sequence=0;self.failed=False;self.closed=False;self.next_excitation={};self.frame=None
        self.native_state=native.snapshot();self.mechanical_state=plant.snapshot()
        self.intake_mass_bridge=intake_mass_bridge;self.intake_mass_binding=deepcopy(intake_mass_binding)
        if intake_mass_bridge is not None:
            if intake_mass_bridge.plant is not plant or intake_mass_bridge.previous!=self.native_state.get('intake'):
                raise ValueError('Intake bridge must bind these exact fresh owners')
            if intake_mass_bridge.failed or intake_mass_bridge.sequence!=0 or intake_mass_bridge.applied_mass!=0 or self.native_state['intake'].get('consumed_count')!=0:
                raise ValueError('Fresh unused intake bridge required')
            if not isinstance(intake_mass_binding,dict) or intake_mass_binding.get('registration_sha256')!=intake_mass_bridge.registration_identity:
                raise ValueError('Intake binding provenance differs from bridge')
        elif intake_mass_binding is not None:raise ValueError('Intake binding requires enabled bridge')
        self.afferents=afferents;self.afferent_input=None;self.afferent_receipt=None
        if afferents is not None:
            if afferents.time_s!=0 or afferents.brain is not neural.brain:raise ValueError('Fresh shared-brain native afferent bridge required')
            from .native_afferents import endpoint_receipt
            self.afferent_receipt=endpoint_receipt(self.native_state,afferents.native_identity,afferents.source_sha256)
        self.cutaneous=cutaneous;self.cutaneous_state=None
        if cutaneous is not None and abs(cutaneous.time_s)>1e-9:raise ValueError('Fresh cutaneous clock required')
        self.cleanup_owner=CleanupOwners(native,plant)
        self.intakes=IntakeSchedule(horizon_s=getattr(getattr(native,'config',None),'horizon_s',86400))
        import hashlib,json
        self.metabolic_reference=deepcopy(self.mechanical_state['metabolic_reference'] if metabolic_reference is None else metabolic_reference)
        for key in ('M0_w','H0_w','W0_w'):finite(self.metabolic_reference[key],'native reference '+key)
        if abs(self.metabolic_reference['M0_w']-self.metabolic_reference['H0_w']-self.metabolic_reference['W0_w'])>1e-9*(1+sum(abs(self.metabolic_reference[k]) for k in ('M0_w','H0_w','W0_w'))):
            raise ValueError('Metabolic reference chemical/heat/work ledger mismatch')
        self.reference_metabolic_w=self.metabolic_reference['M0_w']
        # Zero is the exact initial condition: at t=0 the plant is the reference.
        self.metabolic_exchange_tau_s=METABOLIC_EXCHANGE_TAU_S
        self.metabolic_filter={'m_w':0.,'h_w':0.}
        self.metabolic_pending_energy_j={'m_j':0.,'h_j':0.,'w_j':0.}
        self.metabolic_reference_id=hashlib.sha256(json.dumps({'native_execution_sha256':reference_identity,'reference':self.metabolic_reference},sort_keys=True).encode()).hexdigest()
        if abs(self.native_state['elapsed_s'])>1e-9 or abs(self.mechanical_state['time_s'])>1e-9:
            raise ValueError('Fresh common native and mechanical clocks required')

    def schedule_intakes(self,data):
        if self.failed or self.closed:raise RuntimeError('Embodied runtime is no longer accepting inputs')
        if not isinstance(data,dict) or set(data)!={'events'} or not isinstance(data['events'],list) or not 1<=len(data['events'])<=256:
            raise ValueError('Expected 1–256 intake events')
        events=[]
        for event in data['events']:
            if not isinstance(event,dict) or set(event)!={'event_id','time_s','meal'} or not isinstance(event['meal'],dict):raise ValueError('Invalid intake event')
            events.append(IntakeEvent(event['event_id'],event['time_s'],Meal(**event['meal'])))
        if self.intake_mass_bridge is not None:
            from .intake_mass import planned_nutrition_mass_kg
            planned=math.fsum(planned_nutrition_mass_kg(e.meal) for e in (*self.intakes.events,*events))
            if planned>self.intake_mass_bridge.snapshot()['capacity_kg']:
                raise ValueError('Intake schedule exceeds validated mechanical payload capacity; no events added')
        self.intakes.add_events(events,round(self.time_s*50))
        self.sequence+=1
        try:return self.snapshot()
        except BaseException as error:
            self.failed=True
            raise RuntimeError('Intake schedule changed but frame publication failed') from error

    def _abort(self):
        self.failed=True;self.cleanup_failures=[]
        for label,close in [('physiology',lambda:self.native.close(graceful=False)),('mechanics',self.plant.close)]:
            try:close()
            except BaseException as error:self.cleanup_failures.append({'owner':label,'error':str(error)})

    def _validate(self,data):
        if not isinstance(data,dict) or set(data)-{'seconds','forces','descending','sensory_blocks','motor_blocks','skin_compression_pa','skin_sensory_blocks','regional_skin_pressures','native_sensory_blocks'}:
            raise ValueError('Unknown embodied input')
        dt=finite(data.get('seconds',.02),'embodied interval',.02,.02)
        forces=data.get('forces',[])
        if not isinstance(forces,list) or len(forces)>32:raise ValueError('Expected at most32 force ports')
        for f in forces:
            if not isinstance(f,dict) or set(f)!={'id','force_n','point_m'} or not isinstance(f['id'],str):raise ValueError('Invalid force port')
            for key in ('force_n','point_m'):
                if not isinstance(f[key],(list,tuple)) or len(f[key])!=3:raise ValueError('Expected spatial force/point vector')
                for value in f[key]:finite(value,key,-1000 if key=='force_n' else -5,1000 if key=='force_n' else 5)
        if 'skin_compression_pa' in data:
            finite(data['skin_compression_pa'],'whole-skin pressure',0,5000)
            if not getattr(self.native,'supports_whole_skin_compression',True):raise ValueError('Whole-Skin pressure topology is unavailable')
        if 'regional_skin_pressures' in data:
            pressures=data['regional_skin_pressures']
            if not hasattr(self.native,'regional_skin_pressure') or not isinstance(pressures,dict) or set(pressures)-set(self.native.regions):raise ValueError('Unavailable or unknown regional skin pressure boundary')
            if 'skin_compression_pa' in data:raise ValueError('Whole-Skin and regional pressure topologies cannot be combined')
            for pressure in pressures.values():finite(pressure,'regional skin pressure',0,5000)
        blocks=data.get('skin_sensory_blocks',[])
        if not isinstance(blocks,(list,tuple)) or any(not isinstance(k,str) or self.cutaneous is None or k not in self.cutaneous.sites for k in blocks):raise ValueError('Unknown skin sensory block')
        if 'native_sensory_blocks' in data and self.afferents is None:raise ValueError('Native afferent input is unavailable')
        horizon=getattr(getattr(self.native,'config',None),'horizon_s',None)
        if horizon is not None and self.time_s+dt>horizon+1e-9:raise ValueError('Native horizon reached; no owners advanced')
        return dt,deepcopy(forces)

    def step(self,data):
        if self.failed or self.closed:raise RuntimeError('Embodied runtime is no longer advancing')
        dt,forces=self._validate(data)
        object_forces=[]
        if self.environment_dynamics is not None and hasattr(self.environment_dynamics,'owns'):
            raw_objects=[force for force in forces if self.environment_dynamics.owns(force['id'])]
            object_forces=self.environment_dynamics.validate_object_forces(raw_objects)
            forces=[force for force in forces if not self.environment_dynamics.owns(force['id'])]
        p_checkpoint=n_checkpoint=c_checkpoint=a_checkpoint=e_checkpoint=None;native_touched=False
        try:
            p_checkpoint=self.plant.checkpoint();n_checkpoint=self.neural.checkpoint()
            if self.environment_dynamics is not None:
                e_checkpoint=self.environment_dynamics.checkpoint()
            v=self.native_state['values']
            cutaneous=None;additional={}
            if self.cutaneous is not None:
                c_checkpoint=self.cutaneous.checkpoint()
                contacts=self.mechanical_state.get('cutaneous_contacts')
                if not isinstance(contacts,list):raise ValueError('Explicit mechanical cutaneous contacts required')
                required={key for key,site in self.cutaneous.sites.items() if site.get('mechanical_input')=='native_indentation'}
                observed=[row.get('id') for row in contacts if isinstance(row,dict)]
                if any(observed.count(key)!=1 for key in required):raise ValueError('Missing or duplicate native skin sensor observation')
                blocks=data.get('skin_sensory_blocks',())
                # Previously accepted receptor endpoints feed this exchange.
                # Immediate blocks also remove already queued site contributions.
                for row in (self.cutaneous_state or {}).get('sites',[]):
                    if row['id'] not in blocks:
                        region=row['sensory_region']
                        additional[region]=min(1000.,additional.get(region,0.)+row['sensory_input_hz'])
                cutaneous=self.cutaneous.step(dt,{'time_s':self.time_s,'contacts':contacts,
                    'skin_temperature_C':v.get('skin_temperature_c')},sensory_blocks=blocks)

            afferent_input=None
            if self.afferents is not None:
                a_checkpoint=self.afferents.checkpoint()
                afferent_input=self.afferents.step(dt,self.afferent_receipt,sensory_blocks=data.get('native_sensory_blocks',()))
                for region,rate in afferent_input['sensory_inputs_hz'].items():
                    additional[region]=finite(additional.get(region,0.)+rate,'combined afferent Hz',0,1000)
            physiology={'mean_arterial_pressure_mmHg' :v['mean_arterial_pressure_mmhg'],
                'oxygen_saturation':v['oxygen_saturation'],'core_temperature_C':v['core_temperature_c']}
            neural_inputs=dict(descending=data.get('descending',{}),sensory_blocks=data.get('sensory_blocks',()),
                motor_blocks=data.get('motor_blocks',()),physiology=physiology,additional_sensory_inputs_hz=additional)
            current_feedback=getattr(self.neural,'current_interval_actuation',False) is True
            if current_feedback:
                from .control_exchange import advance_feedback_exchange
                neural,mechanical,forces,mechanical_exchange=advance_feedback_exchange(
                    self.plant,self.neural,self.environment_dynamics,dt,self.mechanical_state,forces,neural_inputs,object_forces,
                    respiratory_projector=lambda ports,entities:self.respiratory_load.project_load(ports,entities,v['lung_volume_ml']))
            else:
                neural=self.neural.step(dt,self.mechanical_state,**neural_inputs)
                # Legacy output belongs to the next exchange; current blocks
                # still suppress already-delivered excitations now.
                actuation={k:v for k,v in self.next_excitation.items() if k not in data.get('motor_blocks',())}
                mechanical_exchange={'interval_s':dt,'substeps':1}
                if self.environment_dynamics is None:
                    mechanical=self.plant.advance(dt,forces=forces,actuation=actuation)
                else:
                    from .world_exchange import advance_body_world
                    mechanical,forces,mechanical_exchange=advance_body_world(self.plant,self.environment_dynamics,
                        dt,self.mechanical_state,forces,actuation,object_forces,
                        respiratory_projector=lambda ports,entities:self.respiratory_load.project_load(ports,entities,v['lung_volume_ml']))
            load=mechanical_exchange.pop('respiratory_load',None)
            mechanical['world_exchange']=mechanical_exchange
            end=self.time_s+dt
            if abs(mechanical['time_s']-end)>1e-8:raise RuntimeError('Mechanical exchange clock diverged')
            if load is None:load=self.respiratory_load.project_load(forces,mechanical['entities'],v['lung_volume_ml'])
            # Segment resultants cannot determine skin/organ strain work. Keep
            # these visible until a resolved traction mapping is available.
            load['unresolved_contact_wrenches']=deepcopy(mechanical.get('body_environment',{}).get('canonical_wrenches',[]))
            pressure=finite(load['external_pressure_pa'],'respiratory load',-5000,5000)
            work=finite(mechanical['positive_muscle_work_j'],'positive muscle work',0)
            energy=finite(mechanical['muscle_metabolic_energy_j'],'native muscle metabolic energy')
            previous=finite(self.mechanical_state['muscle_metabolic_energy_j'],'previous native muscle metabolic energy')
            metabolic_w=(energy-previous)/dt
            incremental_w=metabolic_w-self.reference_metabolic_w
            signed_work=finite(mechanical['signed_active_fiber_work_j'],'signed active fiber work')-finite(self.mechanical_state['signed_active_fiber_work_j'],'previous signed work')
            heat=finite(mechanical['muscle_heat_energy_j'],'native muscle heat')-finite(self.mechanical_state['muscle_heat_energy_j'],'previous muscle heat')
            if abs((energy-previous)-signed_work-heat)>1e-10*(1+abs(energy-previous)+abs(signed_work)+abs(heat)):raise ValueError('Mechanical chemical/heat/work ledger mismatch')
            delta_w=signed_work/dt-self.metabolic_reference['W0_w']
            delta_h=heat/dt-self.metabolic_reference['H0_w']
            # One shared first-order lag on the chemical and heat increments.
            # The filter is linear and both channels use the same coefficient,
            # so deriving work as their difference keeps the exchanged triple on
            # the same chemical = heat + work identity the raw increments hold.
            # Committed only once the whole interval succeeds, so a rejected
            # step leaves the lag exactly where the last accepted interval did.
            blend=1.-math.exp(-dt/self.metabolic_exchange_tau_s)
            exchanged_m=self.metabolic_filter['m_w']+blend*(incremental_w-self.metabolic_filter['m_w'])
            exchanged_h=self.metabolic_filter['h_w']+blend*(delta_h-self.metabolic_filter['h_w'])
            exchanged_w=exchanged_m-exchanged_h
            native_touched=True
            self.native.respiratory_load(pressure)
            if 'skin_compression_pa' in data:self.native.skin_compression(data['skin_compression_pa'])
            for region,pressure in data.get('regional_skin_pressures',{}).items():self.native.regional_skin_pressure(region,pressure)
            if self.native_state.get('pending_meal',False) is False:
                for event in self.intakes.due(round(self.time_s*50)):
                    try:self.intakes.record_accepted(event.event_id,self.native.meal(event.meal))
                    except BaseException as error:
                        self.intakes.record_uncertain(event.event_id,str(error)[:1024] or type(error).__name__)
                        raise
            native=self.native.signed_step(self.metabolic_reference_id,exchanged_m,exchanged_h,exchanged_w)
            unmet=finite(native['values']['coupling.muscle_unmet_kcal'],'unmet native muscle energy')
            if unmet>1e-12:raise RuntimeError('Native muscle energy demand is unmet; mechanical supply feedback is not yet supported')
            native['signal_metadata']=native_field_metadata(native['values'])
            if abs(native['elapsed_s']-end)>1e-8:raise RuntimeError('Native exchange clock diverged')
            if self.intake_mass_bridge is not None:
                # The plant owns the final mechanical substep; ``mechanical``
                # carries whole-exchange quadrature. Compare like intervals
                # when auditing an endpoint mass update, then retain the sum.
                endpoint_work=finite(self.plant.snapshot()['positive_muscle_work_j'],'endpoint positive muscle work',0)
                self.intake_mass_bridge.apply(native['intake'])
                refreshed=self.plant.snapshot()
                for key in ('time_s','muscle_metabolic_energy_j','signed_active_fiber_work_j','muscle_heat_energy_j','metabolic_reference'):
                    if refreshed[key]!=mechanical[key]:raise RuntimeError('Intake endpoint changed completed mechanical interval '+key)
                if refreshed['positive_muscle_work_j']!=endpoint_work:
                    raise RuntimeError('Intake endpoint changed completed mechanical interval positive_muscle_work_j')
                refreshed['positive_muscle_work_j']=work
                refreshed['world_exchange']=deepcopy(mechanical_exchange)
                mechanical=refreshed
            tissue=self.exchange.observe(native)
            geometry=self.respiratory_load.geometry(native['values']['lung_volume_ml'],mechanical['entities'],end)
            self.next_excitation={} if current_feedback else deepcopy(neural['motor_excitations'])
            next_afferent_receipt=None
            if self.afferents is not None:
                from .native_afferents import endpoint_receipt
                next_afferent_receipt=endpoint_receipt(native,self.afferents.native_identity,self.afferents.source_sha256)
            self.afferent_receipt=next_afferent_receipt;self.afferent_input=afferent_input
            self.metabolic_filter={'m_w':exchanged_m,'h_w':exchanged_h}
            self.metabolic_pending_energy_j={
                'm_j':self.metabolic_pending_energy_j['m_j']+(incremental_w-exchanged_m)*dt,
                'h_j':self.metabolic_pending_energy_j['h_j']+(delta_h-exchanged_h)*dt,
                'w_j':self.metabolic_pending_energy_j['w_j']+(delta_w-exchanged_w)*dt}
            self.native_state=native;self.mechanical_state=mechanical;self.time_s=end;self.sequence+=1
            self.frame={'schema':'ihm.embodied-frame.v1','time_s':end,'sequence':self.sequence,'input_capabilities':self._input_capabilities(),
                'environment_state':None if self.environment_dynamics is None else self.environment_dynamics.frame(),
                'entities':geometry['entities'],'skin_field':geometry['skin_field'],'respiration':geometry,'mechanics':mechanical,'neural':neural,'physiology':native,
                'native_afferents':self._afferent_state(),'cutaneous':cutaneous,'intake_mass':self._intake_mass_audit(),'tissue_exchange':tissue,'respiratory_load':load,'intake_schedule':self.intakes.snapshot(),
                'coupling':{'exchange_interval_s':dt,'motor_exchange_latency_s':0. if current_feedback else dt,
                    'control_interval_s':self.neural.control_interval_s if current_feedback else dt,
                    'current_interval_actuation':current_feedback,
                    'positive_muscle_work_j':work,'native_extra_metabolic_demand_w':incremental_w,
                    'muscle_metabolic_reference_w':self.reference_metabolic_w,'interval_muscle_metabolic_w':metabolic_w,
                    'signed_work_increment_w':delta_w,'muscle_heat_increment_w':delta_h,'metabolic_reference_id':self.metabolic_reference_id,
                    'exchanged_metabolic_increment_w':exchanged_m,'exchanged_heat_increment_w':exchanged_h,'exchanged_work_increment_w':exchanged_w,
                    'metabolic_exchange_tau_s':self.metabolic_exchange_tau_s,
                    'metabolic_reference':deepcopy(self.metabolic_reference),
                    'metabolic_pending_energy_j':deepcopy(self.metabolic_pending_energy_j),
                    'metabolic_exchange_basis':'The raw increments above are the plant\'s instantaneous interval differences. A first-order lag of the stated time constant carries them to the native substrate; it delays energy without creating or destroying it, and is an uncalibrated engineering choice, not an identified substrate kinetic.',
                    'metabolic_law':'Native Umberger signed chemical/heat/work increments relative to fixed reference, exchanged through a declared first-order substrate lag; native muscle-only substrate budget and thermal source. Rejects excessive decrement and unmet supply; no basal/stress overwrite.',
                    'storage_owners':{'articulation_muscle':'native mechanical plant','neural':'pinned IBM plus declared decoder/reflexes',
                        'blood_gas_nutrients_heat':'BioGears','tissue_views':'native-owned compartments, no duplicate storage'},
                    'rollback':'An uncertain native commit terminates this runtime; no serializer-exactness claim'}}
            self.cutaneous_state=cutaneous
            result=clone_snapshot_data(self.frame)
            if hasattr(self.plant,'release'):self.plant.release(p_checkpoint)
            return result
        except BaseException:
            if native_touched:
                self._abort()
            else:
                try:
                    if p_checkpoint is None or n_checkpoint is None:raise RuntimeError('Checkpoint acquisition failed')
                    self.plant.restore(p_checkpoint);self.neural.restore(n_checkpoint)
                    if e_checkpoint is not None:self.environment_dynamics.restore(e_checkpoint)
                    if c_checkpoint is not None:self.cutaneous.restore(c_checkpoint)
                    if a_checkpoint is not None:self.afferents.restore(a_checkpoint)
                    if hasattr(self.plant,'release'):self.plant.release(p_checkpoint)
                except BaseException:self._abort()
            raise

    def _intake_mass_audit(self):
        if self.intake_mass_bridge is None:return {'enabled':False}
        return {'enabled':True,'binding':deepcopy(self.intake_mass_binding),'bridge':self.intake_mass_bridge.snapshot()}

    def _afferent_state(self):
        return {'available':self.afferents is not None,'endpoint_receipt':deepcopy(self.afferent_receipt),
            'applied_input':deepcopy(self.afferent_input),'efferent_control':False}

    def _input_capabilities(self):
        regional=getattr(self.native,'regions',())
        regional=list(regional) if (callable(getattr(self.native,'regional_skin_pressure',None))
            and isinstance(regional,(tuple,list)) and all(isinstance(k,str) and k for k in regional)
            and len(set(regional))==len(regional)) else []
        whole=(callable(getattr(self.native,'skin_compression',None))
            and getattr(self.native,'supports_whole_skin_compression',True) is True and not regional)
        return {'skin_pressure':{'whole_skin':whole,'regional_ids':regional}}

    def snapshot(self):
        if self.frame:
            frame=clone_snapshot_data(self.frame);frame.update(controller=deepcopy(getattr(self.neural,'controller_metadata',{'kind':'regional'})),sequence=self.sequence,intake_schedule=self.intakes.snapshot(),intake_mass=self._intake_mass_audit(),input_capabilities=self._input_capabilities())
            return frame
        native=clone_snapshot_data(self.native_state);native['signal_metadata']=native_field_metadata(native['values'])
        geometry=self.respiratory_load.geometry(native['values']['lung_volume_ml'],self.mechanical_state['entities'],0.)
        return {'schema':'ihm.embodied-frame.v1','controller':deepcopy(getattr(self.neural,'controller_metadata',{'kind':'regional'})),'time_s':0.,'sequence':self.sequence,'input_capabilities':self._input_capabilities(),
            'environment_state':None if self.environment_dynamics is None else self.environment_dynamics.frame(),
            'entities':geometry['entities'],'skin_field':geometry['skin_field'],'respiration':geometry,
            'native_afferents':self._afferent_state(),'mechanics':clone_snapshot_data(self.mechanical_state),'physiology':native,
            'tissue_exchange':self.exchange.observe(self.native_state),'intake_schedule':self.intakes.snapshot(),'intake_mass':self._intake_mass_audit()}

    def close(self):
        if self.closed:return
        self.cleanup_owner.close()
        self.closed=True
