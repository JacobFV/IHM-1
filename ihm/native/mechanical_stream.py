"""Continuing source-native articulated dynamics with opaque State checkpoints."""
from pathlib import Path
import copy,gzip,hashlib,json,os,re,selectors,shutil,subprocess,threading,uuid
import numpy as np
from .instance_mass import variant_identity,pointer_directory
from ..assembly.snapshot_data import clone_snapshot_data

# Wall-clock deadline for one native response. A single ``advance`` is an
# error-controlled Simbody integration whose cost is not bounded by dt: when the
# plant is stiff the controller collapses the step and the engine can spend
# minutes inside one 10 ms interval. The deadline is therefore a diagnostic
# knob, not a correctness parameter -- raising it converts a fast failure into a
# slow one and nothing else.
RESPONSE_TIMEOUT_S=float(os.environ.get('IHM_NATIVE_RESPONSE_TIMEOUT_S','120'))

SOURCE_FILES=('subject_walk_scaled.osim','subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml','subject_walk_scaled_FunctionBasedPathSet.xml','subject_walk_scaled_ContactForceSet.xml','subject_walk_scaled_ContactGeometrySet.xml')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def finite(value):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value) or not np.isfinite(value):raise ValueError('Finite scalar required')
    return float(value)
def vec(value):
    x=np.asarray(value,float)
    if x.shape!=(3,) or not np.isfinite(x).all():raise ValueError('Three finite components required')
    return x

class NativeCommandRejected(ValueError):
    """Native command explicitly rejected after transactional rollback."""

class NativeMechanicalStream:
    """The crude 22-segment plant. THIS IS THE SCAFFOLD, NOT THE BODY.

    22 rigid bodies, 80 muscles, 33 coordinates -- what the engine actually
    integrates, and the only thing that touches the world.  The real body of
    3,816 anatomical entities is posed KINEMATICALLY from these segments via
    `data/derived/anatomy-segment-binding/` and feels nothing itself.

    The crude body has always been a scaffold with a planned disposal; its jobs
    move to the real body in stages, contact before actuation.  The end state is
    the brain driving individual muscles whose contractions move the skeleton.
    Make this good enough to generate honest load and afference; do not optimise
    its fidelity as an end, and never report what it did as what the body did.
    See docs/ACTUATION_STAGES.md.
    """
    def __init__(self,root,output,*,environment='supine',target_mass_kg,augmented_registration=None,surface_contact_manifest=None,surface_sensor_indices=(),bed_material=None,instance_mass_variant=None,initial_pose=None,coordinate_limits=None,scene_objects=None,scene_contact_material=None):
        self.root=Path(root).resolve();self.output=Path(output).resolve()
        if self.output.exists() or not self.output.is_relative_to(self.root):raise ValueError('Fresh owned native output directory required')
        if environment not in ('free','supine','upright') or finite(target_mass_kg)<=0:raise ValueError('Invalid native environment or mass')
        self.identity=uuid.uuid4().hex
        self.instance_mass_variant=None if instance_mass_variant is None else variant_identity(self.root,instance_mass_variant)
        pointer_base=self.root/'data/runtime/mechanical-stream' if self.instance_mass_variant is None else pointer_directory(self.root,self.instance_mass_variant)
        pointer=json.loads((pointer_base/'latest.json').read_text());build=self.root/pointer['build'];manifest=json.loads((build/'manifest.json').read_text())
        if manifest.get('instance_mass_variant')!=self.instance_mass_variant:raise ValueError('Native adapter/instance mass variant binding mismatch')
        for path,digest in manifest['files'].items():
            if sha(self.root/path)!=digest:raise ValueError('Stale native mechanical build: '+path)
        augmentation=None;augmentation_bytes=None;self.muscle_catalog=None;self.source_overrides={}
        if augmented_registration is not None:
            record_path=(self.root/augmented_registration).resolve()
            if not record_path.is_relative_to(self.root):raise ValueError('Augmented registration must be owned')
            augmentation_bytes=record_path.read_bytes();augmentation=json.loads(augmentation_bytes)
            for key in ('model','catalog'):
                relative=Path(augmentation[key+'_path'])
                if relative.is_absolute() or '..' in relative.parts or sha(self.root/relative)!=augmentation[key+'_sha256']:raise ValueError('Augmented source identity mismatch')
            for path,digest in augmentation['sources'].items():
                relative=Path(path)
                if relative.is_absolute() or '..' in relative.parts or sha(self.root/relative)!=digest:raise ValueError('Augmented donor identity mismatch')
            catalog_bytes=(self.root/augmentation['catalog_path']).read_bytes()
            if hashlib.sha256(catalog_bytes).hexdigest()!=augmentation['catalog_sha256']:raise ValueError('Augmented catalog changed while copying')
            self.muscle_catalog=json.loads(catalog_bytes)
            # A registration may replace source files other than the model. The
            # only reason this exists is geometric scaling, and the guard below
            # is the whole point of it: the muscle paths this engine runs on are
            # NOT the model's own GeometryPath point sets. They are 80 fitted
            # polynomials in subject_walk_scaled_FunctionBasedPathSet.xml, which
            # replacePathsWithFunctionBasedPaths substitutes in after the model
            # is loaded. A scaled model with the source path set is a body whose
            # skeleton was resized and whose muscles were not, and nothing in the
            # native engine would complain -- it would integrate, and every
            # musculotendon length in the lower limb would be wrong.
            overrides=augmentation.get('source_overrides') or {}
            if not isinstance(overrides,dict) or set(overrides)-set(SOURCE_FILES[1:]):raise ValueError('Source overrides must name non-model native source files')
            for name,record in overrides.items():
                relative=Path(record['path'])
                if relative.is_absolute() or '..' in relative.parts or sha(self.root/relative)!=record['sha256']:raise ValueError('Overridden source identity mismatch: '+name)
            geometric_scale=augmentation.get('geometric_scale')
            if geometric_scale is not None and finite(geometric_scale)!=1.0 and 'subject_walk_scaled_FunctionBasedPathSet.xml' not in overrides:
                raise ValueError('A geometrically scaled model must supply a scaled FunctionBasedPathSet')
            self.source_overrides=overrides
        self.surface_sensor_identity={}
        surface_manifest=None;surface_bytes=None;surface_input=None
        if surface_contact_manifest is not None:
            if environment!='supine':raise ValueError('Surface foundation requires supine environment')
            surface_path=(self.root/surface_contact_manifest).resolve()
            if not surface_path.is_relative_to(self.root):raise ValueError('Owned surface contact manifest required')
            surface_bytes=surface_path.read_bytes();surface_manifest=json.loads(surface_bytes)
            if surface_manifest.get('schema')!='ihm.supine-skin-foundation.v1':raise ValueError('Unknown surface foundation schema')
            for path,digest in surface_manifest['source_files'].items():
                relative=Path(path)
                if relative.is_absolute() or '..' in relative.parts or sha(self.root/relative)!=digest:raise ValueError('Surface foundation source identity mismatch: '+path)
            for key in ('native_input','arrays'):
                relative=Path(surface_manifest[key+'_path'])
                if relative.is_absolute() or '..' in relative.parts or sha(self.root/relative)!=surface_manifest[key+'_sha256']:raise ValueError('Surface foundation artifact identity mismatch')
            surface_input=(self.root/surface_manifest['native_input_path']).read_bytes()
            selected=list(surface_sensor_indices)
            if len(selected)>128 or len(set(selected))!=len(selected) or any(type(i)!=int or not 0<=i<surface_manifest['points'] for i in selected):raise ValueError('At most 128 unique valid skin sensor indices required')
            with np.load(self.root/surface_manifest['arrays_path']) as arrays:
                for index in selected:self.surface_sensor_identity[index]={'manifest_sha256':hashlib.sha256(surface_bytes).hexdigest(),'quadrature_index':index,'triangle_index':int(arrays['face_indices'][index])}
        elif surface_sensor_indices:raise ValueError('Skin sensor selection requires surface contact manifest')
        bed=None
        if bed_material is not None:
            if surface_manifest is None:raise ValueError('Measured bed requires explicit surface foundation')
            from ihm.assembly.bed_compression import load_bed
            bed=load_bed(self.root,bed_material)
            bed['implementation_sha256']=sha(self.root/'ihm/assembly/bed_compression.py')
        self.output.mkdir(parents=True);source=self.output/'inputs';source.mkdir()
        stops=None
        if coordinate_limits is not None:
            stops=[]
            for item in coordinate_limits:
                if set(item)!={'coordinate','lower_rad','upper_rad','stiffness_nm_per_rad','damping_nm_s_per_rad','transition_rad'}:raise ValueError('Coordinate limit requires coordinate/lower_rad/upper_rad/stiffness_nm_per_rad/damping_nm_s_per_rad/transition_rad')
                if not isinstance(item['coordinate'],str) or re.fullmatch(r'[A-Za-z0-9_]+',item['coordinate']) is None:raise ValueError('Invalid coordinate limit name')
                row={k:(item[k] if k=='coordinate' else finite(item[k])) for k in ('coordinate','lower_rad','upper_rad','stiffness_nm_per_rad','damping_nm_s_per_rad','transition_rad')}
                if row['upper_rad']<=row['lower_rad'] or row['stiffness_nm_per_rad']<=0 or row['damping_nm_s_per_rad']<0 or row['transition_rad']<=0:raise ValueError('Invalid coordinate limit values')
                stops.append(row)
            if not stops or len({r['coordinate'] for r in stops})!=len(stops):raise ValueError('Nonempty unique coordinate limits required')
            (source/'coordinate_limits.txt').write_text('IHM_COORDINATE_LIMITS_V1 '+str(len(stops))+'\n'+''.join(
                ' '.join([r['coordinate'],*(str(r[k]) for k in ('lower_rad','upper_rad','stiffness_nm_per_rad','damping_nm_s_per_rad','transition_rad'))])+'\n' for r in stops))
        objects=None
        if scene_objects is not None:
            material=dict(stiffness_pa=1000000.,dissipation_s_m=2.,static_friction=.8,dynamic_friction=.8,viscous_friction=.5,transition_velocity_m_s=.2)
            if scene_contact_material is not None:
                if set(scene_contact_material)-set(material):raise ValueError('Unknown scene contact material field')
                material.update({k:finite(v) for k,v in scene_contact_material.items()})
            if any(material[k]<0 for k in material) or material['stiffness_pa']<=0 or material['transition_velocity_m_s']<=0:raise ValueError('Invalid scene contact material')
            objects=[]
            for item in scene_objects:
                if set(item)!={'id','radius_m','mass_kg','position_m'}:raise ValueError('Scene object requires id/radius_m/mass_kg/position_m')
                if not isinstance(item['id'],str) or re.fullmatch(r'[A-Za-z0-9_]+',item['id']) is None:raise ValueError('Invalid scene object id')
                row=dict(id=item['id'],radius_m=finite(item['radius_m']),mass_kg=finite(item['mass_kg']),position_m=list(vec(item['position_m'])))
                if row['radius_m']<=0 or row['mass_kg']<=0:raise ValueError('Scene object needs positive radius and mass')
                objects.append(row)
            if not objects or len(objects)>8 or len({r['id'] for r in objects})!=len(objects):raise ValueError('One to eight uniquely named scene objects required')
            self.scene_contact_material=material
            (source/'scene_objects.txt').write_text(
                ' '.join(map(str,['IHM_SCENE_OBJECTS_V1',len(objects),material['stiffness_pa'],material['dissipation_s_m'],
                                  material['static_friction'],material['dynamic_friction'],material['viscous_friction'],
                                  material['transition_velocity_m_s']]))+'\n'
                +''.join(' '.join(map(str,[r['id'],r['radius_m'],r['mass_kg'],*r['position_m']]))+'\n' for r in objects))
        pose=None
        if initial_pose is not None:
            if not isinstance(initial_pose,dict) or not initial_pose or any(not isinstance(name,str) or re.fullmatch(r'[A-Za-z0-9_]+',name) is None for name in initial_pose):raise ValueError('Initial pose requires a nonempty coordinate mapping')
            pose={name:finite(value) for name,value in initial_pose.items()}
            (source/'initial_pose.txt').write_text('IHM_INITIAL_POSE_V1 '+str(len(pose))+'\n'+''.join(name+' '+str(value)+'\n' for name,value in pose.items()))
        original=self.root/'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking';inputs={}
        if pose is not None:inputs['initial_pose.txt']=sha(source/'initial_pose.txt')
        if stops is not None:inputs['coordinate_limits.txt']=sha(source/'coordinate_limits.txt')
        if objects is not None:inputs['scene_objects.txt']=sha(source/'scene_objects.txt')
        overrides=self.source_overrides
        for name in SOURCE_FILES:
            if augmentation is not None and name=='subject_walk_scaled.osim':origin=self.root/augmentation['model_path']
            elif name in overrides:origin=self.root/overrides[name]['path']
            else:origin=original/name
            data=origin.read_bytes();(source/name).write_bytes(data);inputs[name]=hashlib.sha256(data).hexdigest()
        if augmentation is not None:
            (source/'augmentation_registration.json').write_bytes(augmentation_bytes)
            (source/'augmentation_catalog.json').write_bytes(catalog_bytes)
            if inputs['subject_walk_scaled.osim']!=augmentation['model_sha256']:raise ValueError('Augmented model changed while copying')
        if surface_manifest is not None:
            (source/'supine_surface_foundation.txt').write_bytes(surface_input)
            (source/'surface_contact_manifest.json').write_bytes(surface_bytes)
            inputs['supine_surface_foundation.txt']=hashlib.sha256(surface_input).hexdigest()
            sensor_bytes=(' '.join(map(str,[len(self.surface_sensor_identity),*self.surface_sensor_identity]))+'\n').encode()
            (source/'surface_sensor_indices.txt').write_bytes(sensor_bytes);inputs['surface_sensor_indices.txt']=hashlib.sha256(sensor_bytes).hexdigest()
        if bed is not None:
            lines=[' '.join(map(str,['IHM_BED_COMPRESSION_V1',bed['material'],bed['thickness_m'],len(bed['strain'])]))]
            lines.extend(' '.join(map(str,pair)) for pair in zip(bed['strain'],bed['pressure_pa']))
            bed_bytes=('\n'.join(lines)+'\n').encode();(source/'bed_compression.txt').write_bytes(bed_bytes)
            (source/'bed_material_identity.json').write_text(json.dumps(bed,indent=2)+'\n');inputs['bed_compression.txt']=hashlib.sha256(bed_bytes).hexdigest()
        runtime=self.root/'data/runtime/opensim';libdirs=[runtime/'install/opensim/lib',runtime/'install/simbody/lib']+[runtime/'sysroot/usr/lib/aarch64-linux-gnu'/p for p in ('lapack','blas','')]
        if self.instance_mass_variant is not None:libdirs.insert(0,self.root/self.instance_mass_variant['path'])
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LD_LIBRARY_PATH=':'.join(map(str,libdirs)),IHM_INSTANCE_MASS_MODE='1' if self.instance_mass_variant is not None else '0',IHM_MASS_REFERENCE_ID=self.identity if self.instance_mass_variant is not None else '')
        limiter=shutil.which('prlimit')
        if not limiter:raise RuntimeError('prlimit is required to enforce the native mechanical memory ceiling')
        engine_command=[str(build/'native_mechanical_stream'),str(source),str(self.output),environment,str(float(target_mass_kg))]
        command=[limiter,'--as=4294967296','--','nice','-n','10',*engine_command]
        self.lock=threading.RLock();self.closed=False;self.tokens=set();self.log=(self.output/'engine.log').open('w')
        execution={'schema':'ihm.native-mechanical-stream.v1','command':command,'engine_command':engine_command,'address_space_limit_bytes':4294967296,'source_sha256':inputs,'build':manifest,'target_mass_kg':target_mass_kg,
                   'scene_objects':objects,'scene_contact_material':None if objects is None else self.scene_contact_material,
                   'scene_object_contact_basis':'Free rigid spheres with SimTK sphere-sphere collision through OpenSim HuntCrossleyForce, one force per object/body-element pair so no object is tested against another object or against the body\'s own overlapping proxies. Contact material transferred from the source foot contact set; the body elements touched are engineering proxies, not an anatomical skin surface.',
                   'coordinate_limits':stops,'coordinate_limits_basis':'Joint stops at the coordinate ranges the source model already declares; nothing else in this plant enforces them. The LIMIT is the model\'s own -- the stiffness, damping and transition width are explicit engineering constants stated by the caller, not measured ligament properties.',
                   'initial_pose':pose,'initial_pose_basis':'Explicit source coordinate initialization after contact reference construction, before muscle equilibrium and energy reference; no ongoing pose constraint or equilibrium claim',
                   'instance_mass_variant':self.instance_mass_variant,'mass_reference_id':self.identity if self.instance_mass_variant is not None else None,'bed_material':bed,'surface_contact_manifest':surface_manifest,'augmented_registration':augmentation,'checkpoint_scope':'Complete in-process SimTK State including effective mass/inertia, excitation/load commands, work accumulators and local mass owner inventory/sequence/receipt; native process must remain alive. Call release(checkpoint) after accepted intervals.',
                   'support_scope':('Retained posterior skin quadrature with prior-based confined neo-Hookean layers in series with the explicitly identified measured conservative mattress compression curve; no calibrated damping or hysteresis, and an infinite native support-plane footprint.' if bed is not None else 'Retained posterior skin quadrature with prior-based confined neo-Hookean layers against an infinite rigid support plane; no mattress compliance.' if surface_manifest is not None else 'Upright source foot contacts plus unilateral non-foot COM spheres inscribed in inertia-derived ellipsoids; transferred source foot material, no balance support or anatomical skin fidelity.' if environment=='upright' else 'No contact support in free environment.' if environment=='free' else 'Supine unilateral engineering posterior spheres from source COM and inertia ellipsoid approximation; source foot contact parameters transferred, not calibrated mattress.'),
                   'external_work_scope':'Endpoint trapezoidal point-force power; not exact integration or metabolic energy'}
        (self.output/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')
        self.process=subprocess.Popen(command,cwd=self.output,env=env,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,bufsize=0)
        self._buffer=b''
        self._selector=selectors.DefaultSelector();self._selector.register(self.process.stdout,selectors.EVENT_READ)
        try:
            self.state=self._read()
            if self.state.get('mass_transfer',{}).get('enabled')!=(self.instance_mass_variant is not None):raise ValueError('Native mass mode differs from selected adapter')
            if self.instance_mass_variant is not None:
                selected=(self.root/self.instance_mass_variant['library_path']).resolve()
                loaded=[Path(line.split(maxsplit=5)[5]).resolve() for line in Path(f'/proc/{self.process.pid}/maps').read_text().splitlines() if len(line.split(maxsplit=5))==6 and line.split(maxsplit=5)[5].startswith('/')]
                if selected not in loaded or sha(selected)!=self.instance_mass_variant['library_sha256']:raise ValueError('Selected instance mass library not loaded')
                if self.state['mass_transfer']['reference_id']!=self.identity:raise ValueError('Native mass reference mismatch')
                execution['loaded_instance_mass_library_sha256']=sha(selected);execution['loaded_instance_mass_library_path']=str(selected)
                (self.output/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')
            if self.muscle_catalog is not None and set(self.state['muscles'])!={m['id'] for m in self.muscle_catalog}:raise ValueError('Native muscle coverage differs from selected catalog')
        except BaseException:self.close();raise
    def _read(self):
        while True:
            while b'\n' not in self._buffer:
                if not self._selector.select(RESPONSE_TIMEOUT_S):raise TimeoutError('Native mechanical stream response timed out after '+str(RESPONSE_TIMEOUT_S)+'s')
                block=os.read(self.process.stdout.fileno(),65536)
                if not block:raise RuntimeError('Native mechanical stream terminated; see '+str(self.output/'engine.log'))
                self._buffer+=block
            raw,self._buffer=self._buffer.split(b'\n',1);line=raw.decode()
            if not line.startswith('@IHM '):self.log.write(line);self.log.flush();continue
            data=json.loads(line[5:],parse_constant=lambda value:(_ for _ in ()).throw(ValueError('Nonfinite native response')))
            if 'error' in data:raise NativeCommandRejected(data['error'])
            self._bind_surface_sensors(data)
            return data
    def _bind_surface_sensors(self,data):
        points=data.get('surface_foundation',{}).get('sensor_points',[])
        if not isinstance(points,list) or any(not isinstance(point,dict) or type(point.get('quadrature_index'))!=int for point in points):
            raise ValueError('Invalid native skin sensor records')
        indices=[point['quadrature_index'] for point in points]
        if len(indices)!=len(set(indices)) or any(index not in self.surface_sensor_identity for index in indices):
            raise ValueError('Duplicate or unrequested native skin sensor')
        if data.get('kind') in ('initialized','observed','restored','advanced','mass_transferred') and set(indices)!=set(self.surface_sensor_identity):
            raise ValueError('Missing selected native skin sensor observation')
        for point in points:
            index=point['quadrature_index']
            point.update(id='skin-contact-'+str(index),material_identity=copy.deepcopy(self.surface_sensor_identity[index]),
                indentation_basis='native nonlinear confined-layer modeled compression',area_basis='projected reference quadrature area',
                point_basis='rigid reference material attachment in native source world')

    def _request(self,command):
        with self.lock:
            if self.closed:raise ValueError('Native mechanical stream closed')
            try:
                payload=(command+'\n').encode()
                if self.process.stdin.write(payload)!=len(payload):raise OSError('Short native command write')
                self.process.stdin.flush();result=self._read()
                expected={'advance':'advanced','observe':'observed','checkpoint':'checkpointed','restore':'restored','drop':'dropped','body_point':'body_point','mass_transfer':'mass_transferred','moment_arms':'moment_arms'}.get(command.split()[0])
                if expected is not None and result.get('kind')!=expected:raise ValueError('Unexpected native command response kind')
                return result
            except NativeCommandRejected:raise
            except BaseException:
                # Dispatch may have succeeded. Never reuse a stream with an uncertain reply.
                self.close();raise
    def snapshot(self):
        # The native state is plain parsed JSON that is only ever rebound, never
        # mutated in place, so a pickle round trip is an identical deep copy and
        # about four times cheaper than tree-walking deepcopy on this shape.
        with self.lock:return clone_snapshot_data(self.state)
    def advance(self,dt_s,forces=(),actuation=None,coordinate_actuation=None):
        """Integrate one interval under muscle excitation and torque-port commands.

        ``coordinate_actuation`` commands the source model's declared
        CoordinateActuators (lumbar, and both arms) in [-1,1] of optimal force.
        They are torque ports, not muscles, and are reported separately in
        ``state['coordinate_actuators']`` so nothing can count them as muscle.
        Omitting the argument leaves the previous commands in force; they start
        at zero, which is exactly what the model did before the ports existed.
        """
        dt=finite(dt_s)
        if not 0<dt<=.02:raise ValueError('Native mechanical step must be in (0,.02]s')
        commands=[]
        for item in forces:
            if set(item)!={'body','point_m','force_n'} or item['body'] not in self.state['bodies']:raise ValueError('Unknown force body/fields')
            commands += [item['body'],*map(str,vec(item['point_m'])),*map(str,vec(item['force_n']))]
        activation={} if actuation is None else actuation
        for name,value in activation.items():
            if name not in self.state['muscles'] or not 0<=finite(value)<=1:raise ValueError('Unknown muscle excitation or out-of-range value')
        line=['advance',str(dt),str(len(forces)),*commands,str(len(activation))]
        for name,value in activation.items():line += [name,str(float(value))]
        if coordinate_actuation is not None:
            ports=self.state.get('coordinate_actuators',{})
            for name,value in coordinate_actuation.items():
                if name not in ports or not -1<=finite(value)<=1:raise ValueError('Unknown coordinate actuator or out-of-range command')
            line += [str(len(coordinate_actuation))]
            for name,value in coordinate_actuation.items():line += [name,str(float(value))]
        with self.lock:self.state=self._request(' '.join(line));return self.snapshot()
    def body_point(self,*,body,station_m):
        """Read current station position/velocity in native source ground frame."""
        if body not in self.state['bodies']:raise ValueError('Unknown native body')
        station=list(vec(station_m))
        with self.lock:
            result=self._request(' '.join(['body_point',body,*map(str,station)]))
            try:
                if result.get('kind')!='body_point' or result.get('body')!=body or result.get('station_m')!=station or result.get('time_s')!=self.state['time_s']:raise ValueError('Native body-point receipt mismatch')
                vec(result['point_source_m']);vec(result['velocity_source_m_s'])
            except BaseException:self.close();raise
            return result
    def moment_arms(self,*,muscles,coordinates):
        """Query actual rotational source path moment arms, without advancing."""
        muscles=list(muscles);coordinates=list(coordinates)
        if not muscles or len(muscles)!=len(set(muscles)) or any(m not in self.state['muscles'] for m in muscles):raise ValueError('Unique known moment arm muscles required')
        if not coordinates or len(coordinates)!=len(set(coordinates)) or any(c not in self.state['coordinates'] for c in coordinates):raise ValueError('Unique known moment arm coordinates required')
        if any(self.state['coordinates'][c]['unit']!='rad' for c in coordinates):raise ValueError('Moment arms in metres require rotational coordinates')
        with self.lock:
            result=self._request(' '.join(['moment_arms',str(len(muscles)),*muscles,str(len(coordinates)),*coordinates]))
            try:
                if result.get('time_s')!=self.state['time_s'] or set(result.get('moment_arms_m',{}))!=set(muscles):raise ValueError('Moment arm receipt mismatch')
                for row in result['moment_arms_m'].values():
                    if set(row)!=set(coordinates):raise ValueError('Moment arm coordinate receipt mismatch')
                    for value in row.values():finite(value)
            except BaseException:self.close();raise
            return result
    def evaluate_static_pose(self,coordinates,*,activations=None):
        """Read-only native pose/activation residual on a copied equilibrium State.

        Activations modify candidate muscle states, never continuing excitations.
        Native muscle laws may clamp activation floors; receipt reports actuals.
        """
        if not isinstance(coordinates,dict) or not coordinates or any(k not in self.state['coordinates'] for k in coordinates):raise ValueError('Nonempty known static coordinates required')
        line=['evaluate_static_pose',str(len(coordinates))]
        for name,value in coordinates.items():line += [name,str(finite(value))]
        if activations is not None:
            if not isinstance(activations,dict) or any(k not in self.state['muscles'] for k in activations):raise ValueError('Known static muscle activation names required')
            line += [str(len(activations))]
            for name,value in activations.items():
                value=finite(value)
                if not 0<=value<=1:raise ValueError('Static activation must be in[0,1]')
                line += [name,str(value)]
        with self.lock:
            result=self._request(' '.join(line))
            try:
                if result.get('kind')!='static_pose_evaluated' or result.get('time_s')!=self.state['time_s'] or result.get('continuing_state_unchanged') is not True or result.get('physical_time_advanced_s')!=0:raise ValueError('Static candidate receipt mismatch')
                if activations is not None and set(result.get('activation_overrides',{}))!=set(activations):raise ValueError('Static activation receipt mismatch')
            except BaseException:self.close();raise
            return result
    def transfer_mass(self,*,sequence,owner,delta_mass_kg,station_m,velocity_source_m_s,body='torso'):
        """Explicit endpoint payload transaction; no physiology mass inference."""
        if self.instance_mass_variant is None:raise ValueError('Instance mass mode was not explicitly selected')
        if type(sequence) is not int or not 0<sequence<2**64:raise ValueError('Positive uint64 mass sequence required')
        if not isinstance(owner,str) or re.fullmatch(r'[A-Za-z0-9_.:-]{1,128}',owner) is None or body!='torso':raise ValueError('Bound owner and supported torso target required')
        dm=finite(delta_mass_kg)
        if abs(dm)>.5:raise ValueError('Local mass transfer must be within +/-0.5 kg per transaction')
        line=['mass_transfer',self.identity,str(sequence),owner,body,str(dm),*map(str,vec(station_m)),*map(str,vec(velocity_source_m_s))]
        with self.lock:
            result=self._request(' '.join(line))
            try:
                ledger=result['mass_transfer'];receipt=ledger['last_receipt']
                if result.get('kind')!='mass_transferred' or ledger['enabled'] is not True or ledger['reference_id']!=self.identity or ledger['last_sequence']!=sequence or receipt['sequence']!=sequence or receipt['zero'] is not (dm==0):raise ValueError('Native mass transaction receipt mismatch')
                if dm!=0:
                    prior=self.state['mass_transfer']['owners'].get(owner,{}).get('mass_kg',0)
                    if receipt['owner']!=owner or receipt['body']!=body or receipt['delta_mass_kg']!=dm or ledger['owners'][owner]['station_m']!=list(vec(station_m)) or ledger['owners'][owner]['mass_kg']!=prior+dm:raise ValueError('Native mass payload receipt mismatch')
                if abs(result['mass_kg']-self.state['mass_kg']-dm)>1e-9:raise ValueError('Native effective mass receipt mismatch')
            except BaseException:self.close();raise
            self.state=result;return self.snapshot()
    def checkpoint(self):
        with self.lock:
            token=uuid.uuid4().hex;self._request('checkpoint '+token);self.tokens.add(token)
            # Rollback is owned by the native process against the token. A Python
            # copy of the state was never read back by restore or by any caller,
            # so it is not taken; time_s is kept as a cheap ownership witness.
            return {'session':self.identity,'token':token,'time_s':self.state['time_s']}
    def restore(self,checkpoint):
        with self.lock:
            if checkpoint.get('session')!=self.identity or checkpoint.get('token') not in self.tokens:raise ValueError('Unknown mechanical checkpoint ownership')
            self.state=self._request('restore '+checkpoint['token']);return self.snapshot()
    def release(self,checkpoint):
        with self.lock:
            if checkpoint.get('session')!=self.identity or checkpoint.get('token') not in self.tokens:raise ValueError('Unknown mechanical checkpoint ownership')
            self._request('drop '+checkpoint['token']);self.tokens.remove(checkpoint['token'])
    def close(self):
        with self.lock:
            if self.closed:return
            self.closed=True
            try:
                if self.process.poll() is None:self.process.stdin.write(b'close\n');self.process.stdin.flush();self.process.wait(timeout=5)
            except (OSError,subprocess.TimeoutExpired):self.process.kill();self.process.wait()
            finally:
                self._selector.close();self.process.stdin.close();self.process.stdout.close();self.log.close()
