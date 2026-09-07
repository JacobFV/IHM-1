"""Continuing, recorded force experiments on the canonical reduced mechanics.

An environment supplies explicit ideal supports, not an unvalidated whole-body
collision replacement. Free spheres have finite mass, rotational inertia and
Coulomb ground contact. The native blood/gas solver is not duplicated here.
"""
from pathlib import Path
import copy
import gzip
import hashlib
import json
import math
import os
import re
import sys
import threading
import time
import types
import uuid

import numpy as np

from .mechanics import BodyMechanics
from . import mechanics as _mechanics_module


def _code_hash(code):
    # marshal encodes string interning/reference sharing, which can differ
    # between a pyc load and fresh compilation of identical implementation.
    def normalized(value):
        if isinstance(value,types.CodeType):
            return {key:normalized(getattr(value,key)) for key in (
                'co_argcount','co_posonlyargcount','co_kwonlyargcount','co_nlocals',
                'co_stacksize','co_flags','co_code','co_consts','co_names',
                'co_varnames','co_freevars','co_cellvars','co_qualname',
                'co_firstlineno','co_linetable','co_exceptiontable')}
        if isinstance(value,tuple):return [normalized(x) for x in value]
        if isinstance(value,frozenset):return sorted((normalized(x) for x in value),key=repr)
        return {'type':type(value).__name__,'value':repr(value)}
    return hashlib.sha256(json.dumps(normalized(code),sort_keys=True).encode()).hexdigest()


def _dataclass_generated(module,raw):
    """Rebuild a restricted inert dataclass declaration, never a module body.

    Dynamic defaults, inheritance and custom decorators fail closed. This is
    sufficient for the configuration records in the native body factory.
    """
    import ast,dataclasses,reprlib
    generated={};hashes={}
    def mismatch():raise ValueError('Loaded code differs from source; unsupported or changed dataclass')
    for node in ast.parse(raw).body:
        cls=vars(module).get(getattr(node,'name',None))
        if not isinstance(node,ast.ClassDef) or not isinstance(cls,type) or not dataclasses.is_dataclass(cls):continue
        if node.bases or node.keywords or len(node.decorator_list)!=1:mismatch()
        decorator=node.decorator_list[0];options={}
        if isinstance(decorator,ast.Call):
            if decorator.args:mismatch()
            for option in decorator.keywords:
                if option.arg is None or not isinstance(option.value,ast.Constant) or type(option.value.value) is not bool:mismatch()
                options[option.arg]=option.value.value
            decorator=decorator.func
        if not isinstance(decorator,ast.Name) or vars(module).get(decorator.id) is not dataclasses.dataclass:mismatch()
        namespace={'__module__':module.__name__,'__qualname__':node.name,'__annotations__':{}}
        user_methods=set()
        for statement in node.body:
            if isinstance(statement,ast.AnnAssign) and isinstance(statement.target,ast.Name):
                name=statement.target.id
                namespace['__annotations__'][name]=object
                if statement.value is not None:
                    try:namespace[name]=ast.literal_eval(statement.value)
                    except (ValueError,TypeError):mismatch()
            elif isinstance(statement,(ast.FunctionDef,ast.AsyncFunctionDef)):
                user_methods.add(statement.name);namespace[statement.name]=object()
            elif isinstance(statement,ast.Pass):pass
            elif isinstance(statement,ast.Expr) and isinstance(statement.value,ast.Constant) and isinstance(statement.value.value,str):pass
            else:mismatch()
        if list(namespace['__annotations__'])!=[field.name for field in dataclasses.fields(cls)]:mismatch()
        try:expected=dataclasses.dataclass(type(node.name,(),namespace),**options)
        except (ValueError,TypeError):mismatch()
        def compare(actual,wanted,label):
            if actual is cls and wanted is expected:return
            if isinstance(actual,types.FunctionType) and isinstance(wanted,types.FunctionType):
                digest=_code_hash(actual.__code__)
                if digest!=_code_hash(wanted.__code__):mismatch()
                hashes[label]=digest
                compare(actual.__defaults__,wanted.__defaults__,label+'.defaults')
                compare(actual.__kwdefaults__,wanted.__kwdefaults__,label+'.kwdefaults')
                if actual.__code__.co_freevars!=wanted.__code__.co_freevars:mismatch()
                for name,a,b in zip(actual.__code__.co_freevars,actual.__closure__ or (),wanted.__closure__ or ()):
                    compare(a.cell_contents,b.cell_contents,label+'.closure.'+name)
                return
            if type(actual) is not type(wanted):mismatch()
            if isinstance(actual,(tuple,list)):
                if len(actual)!=len(wanted):mismatch()
                for index,(a,b) in enumerate(zip(actual,wanted)):compare(a,b,label+'.'+str(index))
            elif isinstance(actual,dict):
                if actual.keys()!=wanted.keys():mismatch()
                for key in actual:compare(actual[key],wanted[key],label+'.'+str(key))
            elif actual is wanted:return
            elif isinstance(actual,(str,int,float,bool,bytes,type(None),set)):
                if actual!=wanted:mismatch()
            else:mismatch()
        for name,method in vars(expected).items():
            if name in user_methods or not isinstance(method,types.FunctionType):continue
            actual=vars(cls).get(name)
            compare(actual,method,node.name+'.'+name)
            generated[(node.name,name)]=actual
    generators={}
    if generated:
        generators={'python_version':sys.version,'implementation':sys.implementation.name,'modules':{}}
        for generator in (dataclasses,reprlib):
            path=Path(generator.__file__).resolve()
            generators['modules'][generator.__name__]={'path':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    return generated,hashes,generators


def _loaded_source(module):
    """Bind retained source to loaded function bytecode, including nested code.

    This records Python method/function implementation identity, not a complete
    environment reconstruction or an identity for NumPy/SciPy native libraries.
    """
    path=Path(module.__file__).resolve();raw=path.read_bytes();compiled={}
    def visit(code):
        compiled[code.co_qualname]=code
        for value in code.co_consts:
            if isinstance(value,types.CodeType):visit(value)
    visit(compile(raw,str(path),'exec',optimize=sys.flags.optimize))
    generated,generated_hashes,generators=_dataclass_generated(module,raw)
    loaded={}
    for name,value in vars(module).items():
        if getattr(value,'__module__',None)!=module.__name__:continue
        if isinstance(value,types.FunctionType):loaded[value.__qualname__]=value.__code__
        elif isinstance(value,type):
            for method_name,method in vars(value).items():
                if (name,method_name) in generated:continue
                if isinstance(method,(classmethod,staticmethod)):method=method.__func__
                if isinstance(method,types.FunctionType):loaded[method.__qualname__]=method.__code__
    hashes={name:_code_hash(code) for name,code in loaded.items()}
    if any(name not in compiled or _code_hash(compiled[name])!=digest for name,digest in hashes.items()):
        raise ValueError('Loaded code differs from source; restart the server before creating a scene')
    return {'path':path,'bytes':raw,'source_sha256':hashlib.sha256(raw).hexdigest(),'loaded_code_sha256':hashes,
            'generated_code_sha256':generated_hashes,'code_generators':generators}

ENVIRONMENTS={
    'studio':dict(label='Studio',gravity=[0.,0.,0.],axis=1,plane=-.96,
                  supports=[],description='Free translation; reference orientations constrained'),
    'floor':dict(label='Floor · supported stance',gravity=[0.,-9.81,0.],axis=1,plane=-.96,
                 supports=['body-bp3d-FJ3256','body-bp3d-FJ3360'],
                 description='Incremental body motion about a supported reference; ground contact for free objects'),
    'bed':dict(label='Bed · supported supine',gravity=[0.,0.,-9.81],axis=2,plane=-.24,
               supports=['body-bp3d-FJ3309','body-bp3d-FJ3393','body-bp3d-FJ3256','body-bp3d-FJ3360'],
               description='Incremental body motion about a supported reference; mattress mechanics unresolved'),
}


def vector(value,name):
    if not isinstance(value,(list,tuple,np.ndarray)) or len(value)!=3 or any(isinstance(v,(bool,np.bool_)) for v in value):
        raise ValueError(name+' requires three finite numbers')
    result=np.asarray(value,dtype=float)
    if result.shape!=(3,) or not np.isfinite(result).all():raise ValueError(name+' requires three finite numbers')
    return result


class Sphere:
    def __init__(self,ident,position,radius=.065,mass=.4):
        self.id=ident;self.position=vector(position,'position');self.velocity=np.zeros(3)
        self.omega=np.zeros(3);self.rotation=np.eye(3)
        if not 0<radius<=1 or not 0<mass<=100:raise ValueError('Positive bounded sphere radius and mass required')
        self.radius=radius;self.mass=mass;self.inertia=.4*mass*radius**2
        self.work=0.;self.contact_loss=0.;self.projection_work=0.;self.support_impulse=np.zeros(3)

    def kinetic(self):
        return .5*self.mass*float(self.velocity@self.velocity)+.5*self.inertia*float(self.omega@self.omega)

    def step(self,dt,force,point,gravity,*,axis,plane,contact=True):
        from .mechanics import _rotation
        before=self.position.copy();omega_before=self.omega.copy()
        torque=np.cross(point-before,force)
        acceleration=gravity+force/self.mass
        self.position+=self.velocity*dt+.5*acceleration*dt**2
        self.velocity+=acceleration*dt;self.omega+=torque/self.inertia*dt
        self.rotation=_rotation(.5*(omega_before+self.omega)*dt)@self.rotation
        self.work+=float(force@(self.position-before)+torque@(.5*(omega_before+self.omega)*dt))
        normal=np.eye(3)[axis]
        penetration=plane+self.radius-self.position[axis]
        if contact and penetration>0:
            self.position+=penetration*normal
            self.projection_work-=self.mass*float(gravity@normal)*penetration
            lever=-self.radius*normal
            relative=self.velocity+np.cross(self.omega,lever)
            vn=float(relative@normal)
            jn=max(0.,-vn*self.mass)  # zero restitution, no manufactured bounce
            tangent=relative-vn*normal;speed=float(np.linalg.norm(tangent))
            inverse_tangent_mass=1/self.mass+self.radius**2/self.inertia
            needed=speed/inverse_tangent_mass
            magnitude=needed if needed<=.5*jn else min(needed,.35*jn)
            jt=-magnitude*tangent/speed if speed else np.zeros(3)
            impulse=jn*normal+jt
            energy=self.kinetic()
            self.velocity+=impulse/self.mass;self.omega+=np.cross(lever,impulse)/self.inertia
            self.contact_loss+=energy-self.kinetic();self.support_impulse-=impulse

    def snapshot(self):
        return dict(id=self.id,kind='sphere',position_m=self.position.tolist(),velocity_m_s=self.velocity.tolist(),
                    angular_velocity_rad_s=self.omega.tolist(),rotation_matrix=self.rotation.tolist(),
                    radius_m=self.radius,mass_kg=self.mass,kinetic_energy_j=self.kinetic(),
                    applied_work_j=self.work,contact_dissipation_j=self.contact_loss,
                    contact_projection_potential_change_j=self.projection_work,
                    ground_support_impulse_ns=self.support_impulse.tolist())


SCENE_CATALOGUE='data/derived/environment-catalogue-v1/catalogue.json'
MAX_SCENE_OBJECTS=16


def scene_catalogue(root):
    path=Path(root)/SCENE_CATALOGUE
    if not path.is_file():raise ValueError('Scene catalogue unavailable; run scripts/build_environment_catalogue.py')
    raw=path.read_bytes();return json.loads(raw),hashlib.sha256(raw).hexdigest()


def insertable_object(catalogue,ident):
    entry=next((o for o in catalogue['objects'] if o['id']==ident),None)
    if entry is None:raise ValueError('Unknown scene object')
    if entry.get('collider')!='sphere':raise ValueError('Object is display only; the engine has no collider for it')
    return entry


def scene_objects(root,scene,environment):
    """Objects a composed scene adds. Only sphere colliders exist, so only they are instantiated."""
    catalogue,digest=scene_catalogue(root)
    record=next((s for s in catalogue['scenes'] if s['id']==scene),None)
    if record is None:raise ValueError('Unknown scene')
    if record['base_environment']!=environment:raise ValueError('Scene requires environment '+record['base_environment'])
    objects=[]
    for placement in record['placements']:
        entry=next((o for o in catalogue['objects'] if o['id']==placement['object']),None)
        if entry is None:raise ValueError('Scene names an object the catalogue does not carry')
        if entry.get('collider')!='sphere':continue  # display only: no collider exists for that shape
        engine_id=placement['instance']
        if engine_id=='scene-ball':continue  # already instantiated as the scene default
        if not re.fullmatch(r'[a-z0-9-]{1,64}',engine_id):raise ValueError('Invalid scene object id')
        objects.append(Sphere(engine_id,placement['start_m'],radius=entry['radius_m'],mass=entry['mass_kg']))
    summary={k:record[k] for k in ('id','label','base_environment','objects','simulated_objects','display_only_objects')}
    summary['catalogue']=SCENE_CATALOGUE
    return objects,summary,digest


class InteractiveScene:
    def __init__(self,root,output,environment='studio',scene=None):
        self.root=Path(root).resolve();self.output=Path(output).resolve()
        if environment not in ENVIRONMENTS:raise ValueError('Unknown scene environment')
        if scene is not None and (not isinstance(scene,str) or not re.fullmatch(r'[a-z0-9-]{1,64}',scene)):raise ValueError('Invalid scene id')
        if not self.output.is_relative_to(self.root) or self.output.exists():raise ValueError('Fresh retained scene directory required')
        self.scene_id=scene;self.scene=None;composed=[];catalogue_digest=None
        if scene is not None:composed,self.scene,catalogue_digest=scene_objects(self.root,scene,environment)
        self.environment=copy.deepcopy(ENVIRONMENTS[environment]);self.environment_id=environment
        path=self.root/'data/derived/canonical/mechanics.json';mechanics_bytes=path.read_bytes();payload=json.loads(mechanics_bytes)
        for receipt in _IMPORT_SOURCES:
            if receipt['path'].read_bytes()!=receipt['bytes']:raise ValueError('Scene source changed after import; restart the server')
        self.body=BodyMechanics(payload);self.lock=threading.Lock();self.sequence=0;self.closed=False
        if set(self.environment['supports'])-set(self.body.ids):raise ValueError('Environment support anatomy missing')
        self.output.mkdir(parents=True);(self.output/'inputs').mkdir();(self.output/'events').mkdir()
        self._event_index=0;self._last_event_sha256=None
        sources={}
        for source,raw in [(path,mechanics_bytes)]+[(r['path'],r['bytes']) for r in _IMPORT_SOURCES]:
            target=self.output/'inputs'/source.name;target.write_bytes(raw)
            sources[str(source.relative_to(self.root))]=hashlib.sha256(target.read_bytes()).hexdigest()
        if scene is not None:
            (self.output/'inputs'/'scene_catalogue.json').write_bytes((self.root/SCENE_CATALOGUE).read_bytes())
            sources[SCENE_CATALOGUE]=catalogue_digest
        self.objects={s.id:s for s in [Sphere('scene-ball', [.42,.1,.3])]+composed}
        self.initial_object_energy={s.id:s.kinetic()-s.mass*np.dot(self.environment['gravity'],s.position) for s in self.objects.values()}
        self.frame=dict(model_id='ihm-body',time_s=0.,entities={ident:dict(
            translation_m=[0.,0.,0.],centroid_m=self.body.x0[i].tolist(),
            rotation_matrix=np.eye(3).tolist(),deformation_gradient=np.eye(3).tolist())
            for i,ident in enumerate(self.body.ids)},audit={})
        self.scope=dict(body_mechanics='Canonical linked rigid translations and affine soft regions',
            body_rotations='Reference orientations constrained; applied moments carried by constraints',
            body_environment='Explicit named ideal supports; no body-surface mattress or floor contact solve',
            body_gravity='Incremental motion assumes a balanced reference preload; gravitational prestress and support distribution are not solved',
            objects='Finite-mass spheres with rotational inertia and Coulomb ground contact',
            scene_objects='A composed scene adds only its sphere-collider objects; its constructed furniture is display geometry the engine never instantiates',
            body_object_contact=False,clothing_contact=False,physiology_feedback=False,
            force_location='Force acts on the selected entity translation; off-centroid moment is an explicit constraint reaction',
            sphere_force_location='A point inside/on the sphere, captured as a material offset and followed through substeps; explicit distant wrench ports unsupported',
            calibration='Existing source/engineering parameters; no new human calibration')
        (self.output/'manifest.json').write_text(json.dumps(dict(schema='ihm.interactive-scene.v2',
            environment_id=environment,environment=self.environment,sources=sources,scope=self.scope,
            scene_id=scene,scene=self.scene,
            simulated_objects=[s.id for s in self.objects.values()],
            loaded_code={r['path'].name:r['loaded_code_sha256'] for r in _IMPORT_SOURCES},
            journal='Immutable individually compressed events/*.json.gz; SHA256 predecessor chain; atomic publication',
            created_unix=time.time(),timestep_s=.002),indent=2)+'\n')
        self._record({'kind':'initialization','frame':self.snapshot()})

    def _record(self,value):
        record={'schema':'ihm.scene-event.v1','event_index':self._event_index,
                'previous_event_sha256':self._last_event_sha256,**value}
        raw=gzip.compress((json.dumps(record,allow_nan=False)+'\n').encode(),compresslevel=3,mtime=0)
        digest=hashlib.sha256(raw).hexdigest();directory=self.output/'events'
        destination=directory/f'{self._event_index:08d}-{value["kind"]}.json.gz'
        temporary=directory/f'.pending-{uuid.uuid4().hex}'
        if destination.exists():raise FileExistsError('Immutable scene event already exists')
        try:
            with temporary.open('xb') as stream:stream.write(raw);stream.flush();os.fsync(stream.fileno())
            # Atomic no-overwrite publication. The temporary hard link is then
            # removed; an accepted event is never opened for writing again.
            os.link(temporary,destination)
        finally:
            try:temporary.unlink(missing_ok=True)
            except OSError:pass  # A pending orphan cannot invalidate publication.
        self._event_index+=1;self._last_event_sha256=digest

    def snapshot(self):
        return {**self.frame,'objects':[s.snapshot() for s in self.objects.values()],
                'environment_id':self.environment_id,'environment':self.environment,
                'scene_id':self.scene_id,'scene':self.scene,
                'scope':self.scope,'sequence':self.sequence,'closed':self.closed}

    def insert(self,data):
        """Additive object insert. Each call is a new instance with its own id; nothing is replaced."""
        if not isinstance(data,dict) or set(data)-{'object','offset_m'}:raise ValueError('Insert accepts object and offset_m')
        catalogue,digest=scene_catalogue(self.root)
        entry=insertable_object(catalogue,data.get('object'))
        offset=vector(data.get('offset_m',[0.,0.,0.]),'offset')
        if np.max(np.abs(offset))>2:raise ValueError('Insert offset exceeds the scene domain')
        with self.lock:
            if self.closed:raise ValueError('Scene is closed')
            if len(self.objects)>=MAX_SCENE_OBJECTS:raise ValueError('Scene object limit reached; close or start a new scene')
            base=entry.get('engine_object_id',entry['id'])
            index=1
            while (base if index==1 else base+'-'+str(index)) in self.objects:index+=1
            ident=base if index==1 else base+'-'+str(index)
            sphere=Sphere(ident,np.asarray(entry['start_m'],float)+offset,radius=entry['radius_m'],mass=entry['mass_kg'])
            self.objects[ident]=sphere
            self.initial_object_energy[ident]=sphere.kinetic()-sphere.mass*np.dot(self.environment['gravity'],sphere.position)
            self.sequence+=1
            result=self.snapshot()
            self._record({'kind':'insert','command':data,'object_id':ident,'catalogue_sha256':digest,'frame':result})
            return result

    def step(self,data):
        if not isinstance(data,dict) or set(data)-{'seconds','forces','sequence'}:raise ValueError('Unknown scene command')
        if type(data.get('sequence')) is not int or data['sequence']<0:raise ValueError('Scene sequence must be a nonnegative integer')
        seconds=data.get('seconds',.02)
        if isinstance(seconds,bool) or not isinstance(seconds,(int,float)) or not math.isfinite(seconds) or not 0<seconds<=.1:
            raise ValueError('Scene steps must be in (0, 0.1] seconds')
        ticks=round(seconds/.002)
        if abs(seconds-ticks*.002)>1e-10:raise ValueError('Scene steps require whole 2 ms intervals')
        commands=data.get('forces',[])
        if not isinstance(commands,list) or len(commands)>8:raise ValueError('At most eight force ports per step')
        validated=[];seen=set()
        for command in commands:
            if not isinstance(command,dict) or set(command)!={'id','force_n','point_m'}:raise ValueError('Force requires id, force_n, point_m')
            ident=command['id']
            if not isinstance(ident,str) or ident in seen or ident not in set(self.body.ids)|set(self.objects):raise ValueError('Unknown or duplicate force owner')
            seen.add(ident);force=vector(command['force_n'],'force');point=vector(command['point_m'],'point')
            if np.linalg.norm(force)>100 or np.max(np.abs(point))>5:raise ValueError('Force or application point exceeds scene domain')
            validated.append((ident,force,point))
        with self.lock:
            if self.closed:raise ValueError('Scene is closed')
            if data.get('sequence')!=self.sequence:raise ValueError('Scene sequence changed; refresh before applying force')
            sphere_ports={}
            for ident,force,point in validated:
                if ident in self.objects:
                    sphere=self.objects[ident];offset=point-sphere.position
                    if np.linalg.norm(offset)>sphere.radius+1e-7:raise ValueError('Force application point is outside its sphere owner; distant wrench unsupported')
                    sphere_ports[ident]=(force,sphere.rotation.T@offset)
            if self.body.time+seconds>120:raise ValueError('Start a new scene after the 120 s research horizon')
            # Capture before mutation so validation/numerical failure cannot
            # leave an unacknowledged partially advanced model.
            state_keys=('x','v','omega','r','deformation','time','work','dissipation','affine_work',
                        '_soft_energy','_soft_previous','prescribed_work','prescribed_reactions','orientation_reactions')
            saved={k:copy.deepcopy(getattr(self.body,k)) for k in state_keys if hasattr(self.body,k)}
            objects=copy.deepcopy(self.objects);oldframe=self.frame;oldsequence=self.sequence
            gravity=np.asarray(self.environment['gravity'])
            # This solver describes increments about the existing supported
            # reference. Its gravity prestress has never been solved. Retain the
            # assumed balancing preload explicitly rather than treating an
            # unsupported gravitational collapse as standing/crawling dynamics.
            reference_gravity=self.body.mass[:,None]*gravity
            reference_support=-reference_gravity
            external=reference_gravity+reference_support
            constraint_moments={}
            for ident,force,point in validated:
                if ident in self.body.index:
                    i=self.body.index[ident];external[i]+=force
                    constraint_moments[ident]=(-np.cross(point-self.body.x[i],force)).tolist()
            drivers={'external_forces_n':{ident:external[i].tolist() for i,ident in enumerate(self.body.ids) if np.any(external[i])},
                     'prescribed_translations_m':{ident:[0.,0.,0.] for ident in self.environment['supports']}}
            try:
                frame=self.body.step(seconds,drivers)
                if np.max(np.linalg.norm(self.body.x-self.body.x0,axis=1))>.15:
                    raise ValueError('Reduced-body displacement domain exceeded; reset or use smaller force')
                for _ in range(ticks):
                    for ident,sphere in self.objects.items():
                        force,offset=sphere_ports.get(ident,(np.zeros(3),np.zeros(3)))
                        applied=(force,sphere.position+sphere.rotation@offset)
                        sphere.step(.002,*applied,gravity,axis=self.environment['axis'],plane=self.environment['plane'],
                                    contact=self.environment_id!='studio')
                frame['external_orientation_constraint_reactions_nm']=constraint_moments
                frame['audit']['assumed_reference_support_force_n']=reference_support.sum(axis=0).tolist()
                frame['audit']['reference_gravity_force_n']=reference_gravity.sum(axis=0).tolist()
                frame['audit']['object_energy_balance_residual_j']={ident:s.kinetic()-s.mass*np.dot(gravity,s.position)
                    -self.initial_object_energy[ident]-s.work+s.contact_loss-s.projection_work for ident,s in self.objects.items()}
                self.frame=frame;self.sequence+=1
                result=self.snapshot();self._record({'kind':'advance','command':data,'frame':result})
                return result
            except BaseException:
                for key in state_keys:
                    if key in saved:setattr(self.body,key,saved[key])
                    elif hasattr(self.body,key):delattr(self.body,key)
                self.objects=objects;self.frame=oldframe;self.sequence=oldsequence
                raise


class SceneSessions:
    def __init__(self,root):self.root=Path(root);self.sessions={};self.lock=threading.Lock()
    def catalog(self):return dict(environments=[dict(id=k,**v) for k,v in ENVIRONMENTS.items()],force_limit_n=100.,step_s=.02)
    def current(self,ident):
        with self.lock:scene=self.sessions.get(ident)
        if scene is None:raise ValueError('Unknown scene session')
        with scene.lock:return copy.deepcopy(dict(id=ident,**scene.snapshot()))
    def create(self,data):
        if not isinstance(data,dict) or set(data)-{'environment','scene'}:raise ValueError('Unknown scene configuration')
        with self.lock:
            if len(self.sessions)>=4:raise ValueError('Close a scene before opening another')
            ident=uuid.uuid4().hex
            scene=InteractiveScene(self.root,self.root/'data/derived/interactive-scenes'/ident,data.get('environment','studio'),data.get('scene'))
            self.sessions[ident]=scene
        return dict(id=ident,**scene.snapshot())
    def command(self,ident,action,data):
        if action=='close' and (not isinstance(data,dict) or data):raise ValueError('Close accepts an empty object')
        with self.lock:scene=self.sessions.get(ident)
        if scene is None and action=='close':return dict(id=ident,closed=True,already_absent=True)
        if scene is None:raise ValueError('Unknown scene session')
        if action=='step':return dict(id=ident,**scene.step(data))
        if action=='insert':return dict(id=ident,**scene.insert(data))
        if action=='close':
            if not isinstance(data,dict) or data:raise ValueError('Close accepts an empty object')
            with scene.lock:
                if not scene.closed:
                    scene._record({'kind':'close','sequence':scene.sequence,'time_s':scene.body.time})
                    scene.closed=True
            with self.lock:self.sessions.pop(ident,None)
            return dict(id=ident,closed=True)
        raise ValueError('Unknown scene action')


_IMPORT_SOURCES=(_loaded_source(_mechanics_module),_loaded_source(sys.modules[__name__]))
