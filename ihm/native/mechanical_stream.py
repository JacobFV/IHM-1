"""Continuing source-native articulated dynamics with opaque State checkpoints."""
from pathlib import Path
import copy,gzip,hashlib,json,os,re,selectors,shutil,subprocess,threading,uuid
import numpy as np
from .instance_mass import variant_identity,pointer_directory

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
    def __init__(self,root,output,*,environment='supine',target_mass_kg,augmented_registration=None,surface_contact_manifest=None,surface_sensor_indices=(),bed_material=None,instance_mass_variant=None):
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
        augmentation=None;augmentation_bytes=None;self.muscle_catalog=None
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
        original=self.root/'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking';inputs={}
        for name in SOURCE_FILES:
            data=((self.root/augmentation['model_path']) if augmentation is not None and name=='subject_walk_scaled.osim' else original/name).read_bytes();(source/name).write_bytes(data);inputs[name]=hashlib.sha256(data).hexdigest()
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
                   'instance_mass_variant':self.instance_mass_variant,'mass_reference_id':self.identity if self.instance_mass_variant is not None else None,'bed_material':bed,'surface_contact_manifest':surface_manifest,'augmented_registration':augmentation,'checkpoint_scope':'Complete in-process SimTK State including effective mass/inertia, excitation/load commands, work accumulators and local mass owner inventory/sequence/receipt; native process must remain alive. Call release(checkpoint) after accepted intervals.',
                   'support_scope':('Retained posterior skin quadrature with prior-based confined neo-Hookean layers; rigid plane, uncalibrated mattress.' if surface_manifest is not None else 'Supine unilateral engineering posterior spheres from source COM and inertia ellipsoid approximation; source foot contact parameters transferred, not calibrated mattress.'),
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
                if not self._selector.select(120):raise TimeoutError('Native mechanical stream response timed out')
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
                expected={'advance':'advanced','observe':'observed','checkpoint':'checkpointed','restore':'restored','drop':'dropped','body_point':'body_point','mass_transfer':'mass_transferred'}.get(command.split()[0])
                if expected is not None and result.get('kind')!=expected:raise ValueError('Unexpected native command response kind')
                return result
            except NativeCommandRejected:raise
            except BaseException:
                # Dispatch may have succeeded. Never reuse a stream with an uncertain reply.
                self.close();raise
    def snapshot(self):
        with self.lock:return copy.deepcopy(self.state)
    def advance(self,dt_s,forces=(),actuation=None):
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
            return {'session':self.identity,'token':token,'snapshot':self.snapshot()}
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
