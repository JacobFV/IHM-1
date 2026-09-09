"""Demand-stepped native body actors: one controller thread per live body."""
from concurrent.futures import Future
from pathlib import Path
import gzip,hashlib,json,os,queue,threading,time,uuid


from ihm.brain.active_source import resolve_ibm_candidate,resolve_source,ACTIVE_COMMIT
from ihm.app.surface_assets import SurfaceAssetRegistry,register_body_surface_assets


class BodyActor:
    def __init__(self,factory,output):
        self.output=Path(output);self.events=self.output/'events';self.events.mkdir(parents=True)
        self.pending=queue.Queue(maxsize=8);self.ready=Future();self.closed=False;self.close_requested=False;self.lifecycle=threading.Lock()
        self.index=0;self.previous=None;self.body=None;self.cleanup_error=None;self.error=None
        self.thread=threading.Thread(target=self._run,args=(factory,),name='ihm-embodied-owner',daemon=True)
        self.thread.start()

    def _publish(self,event):
        payload={'schema':'ihm.embodied-event.v1','index':self.index,'previous_sha256':self.previous,**event}
        raw=gzip.compress(json.dumps(payload,allow_nan=False,separators=(',',':')).encode(),compresslevel=1,mtime=0)
        target=self.events/f'{self.index:08d}.json.gz';temporary=self.events/f'.pending-{uuid.uuid4().hex}'
        try:
            with temporary.open('xb') as f:f.write(raw);f.flush();os.fsync(f.fileno())
            os.link(temporary,target)
        finally:temporary.unlink(missing_ok=True)
        self.previous=hashlib.sha256(raw).hexdigest();self.index+=1

    def _cleanup(self):
        """Only successful owner cleanup releases the machine resource slot."""
        with self.lifecycle:self.close_requested=True
        try:
            if self.body is not None:self.body.close()
        except BaseException as error:
            with self.lifecycle:self.cleanup_error=str(error)
            try:self._publish({'kind':'cleanup_error','error':str(error)})
            except BaseException:pass
            return error
        try:self._publish({'kind':'close'})
        except BaseException as error:self.error=str(error)
        with self.lifecycle:
            self.cleanup_error=None;self.closed=True
            while not self.pending.empty():
                action,_,future=self.pending.get_nowait()
                if action=='close':future.set_result({'closed':True})
                else:future.set_exception(RuntimeError('Body actor closed'))
        return None

    def _run(self,factory):
        try:
            self.body=factory();self._publish({'kind':'initialization','frame':self.body.snapshot()});self.ready.set_result(True)
        except BaseException as error:
            self.error=str(error)
            # A failed constructor may hand back its partially initialized owner.
            if self.body is None:self.body=getattr(error,'cleanup_owner',None)
            self.ready.set_exception(error)
            with self.lifecycle:self.close_requested=True
        if self.close_requested:
            if self._cleanup() is None:return
        while True:
            action,data,future=self.pending.get()
            if action=='close':
                error=self._cleanup()
                if error is None:future.set_result({'closed':True});return
                future.set_exception(error)
                continue
            try:
                if self.close_requested:raise RuntimeError('Body actor is closing')
                if action=='snapshot':result=self.body.snapshot()
                elif action in ('step','intakes'):
                    if not isinstance(data,dict) or type(data.get('sequence')) is not int:raise ValueError('Integer body sequence required')
                    if data['sequence']!=self.body.snapshot()['sequence']:raise ValueError('Body sequence changed; read current state')
                    command={k:v for k,v in data.items() if k!='sequence'}
                    self._publish({'kind':'command','operation':action,'command':data})
                    result=self.body.step(command) if action=='step' else self.body.schedule_intakes(command)
                    try:self._publish({'kind':'advance' if action=='step' else 'intake_schedule','frame':result})
                    except BaseException:
                        self.body.failed=True
                        raise
                else:raise ValueError('Unknown embodied operation')
                future.set_result(result)
            except BaseException as error:
                try:self._publish({'kind':'error','operation':action,'error':str(error),'body_failed':getattr(self.body,'failed',False)})
                except BaseException:self.body.failed=True
                if getattr(self.body,'failed',False):
                    self.error=str(error)
                    future.set_exception(RuntimeError('Body operation outcome uncertain: '+str(error)))
                else:future.set_exception(error)
            if self.cleanup_error is None and (self.close_requested or getattr(self.body,'failed',False)):
                if self._cleanup() is None:return

    def status(self):
        with self.lifecycle:
            error=self.cleanup_error or self.error
            if self.closed:return {'status':'closed','closed':True,'error':error}
            if self.cleanup_error:return {'status':'failed','closed':False,'error':error}
            if self.close_requested:return {'status':'closing','closed':False,'error':error}
            if not self.ready.done():return {'status':'initializing','closed':False}
            if self.ready.exception():return {'status':'failed','closed':False,'error':error}
            return {'status':'ready','closed':False}

    def request_close(self,wait=True):
        with self.lifecycle:
            self.close_requested=True
            initializing=not self.ready.done()
            closed=self.closed
        if initializing or closed:return self.status()
        return self.call('close',wait=wait)

    def call(self,action,data=None,wait=True):
        if action!='close':self.ready.result(timeout=180)
        future=Future()
        with self.lifecycle:
            if self.closed:
                if action=='close':return {'closed':True,'already_absent':True}
                raise RuntimeError('Body actor closed; retained records remain available')
            if self.close_requested and action!='close':raise RuntimeError('Body actor is closing')
            try:self.pending.put_nowait((action,data,future))
            except queue.Full:raise ValueError('Body command queue is full')
        if not wait:return self.status()
        result=future.result(timeout=180)
        if action=='close':self.thread.join(timeout=5)
        return result


class EmbodiedSessions:
    def __init__(self,root):
        self.root=Path(root);self.lock=threading.Lock();self.actors={};self.creating=False;self.shutting_down=False
        self.creation_done=threading.Condition(self.lock)
        self.surface_assets=SurfaceAssetRegistry()
    def create(self,data):
        if not isinstance(data,dict) or set(data)-{'environment','regional_skin','intake_mass','ibm_candidate','environment_selection','controller'}:raise ValueError('Unknown embodied configuration')
        from ihm.assembly.controller_selection import resolve_controller
        controller=resolve_controller(data.get('controller'))
        if controller['kind']!='regional' and 'ibm_candidate' in data:
            raise ValueError('Implicit checkpoint controller cannot select a regional IBM candidate')
        intake_mass=data.get('intake_mass',False)
        if type(intake_mass) is not bool:raise ValueError('intake_mass must be boolean')
        regional_skin=data.get('regional_skin',False)
        if type(regional_skin) is not bool:raise ValueError('regional_skin must be boolean')
        environment=data.get('environment','supine')
        if environment not in ('free','supine','upright'):raise ValueError('Unknown articulated environment')
        from ihm.assembly.controller_selection import STANCE_KINDS
        if controller['kind'] in STANCE_KINDS and (environment!='upright' or intake_mass):
            raise ValueError('Stance controllers require upright fixed-mass mechanics')
        from ihm.assembly.environment_dynamics import resolve_selection
        environment_selection,environment_options=resolve_selection(self.root,environment,data.get('environment_selection'))
        explicit='ibm_candidate' in data
        source_pin=None if controller['kind'] not in ('regional','engineering_stance') else (resolve_ibm_candidate(self.root,data['ibm_candidate']) if explicit else resolve_source(self.root))
        with self.lock:
            if self.shutting_down:raise RuntimeError('Embodied service is shutting down')
            if self.creating or any(not a.closed for a in self.actors.values()):raise ValueError('One native body at a time while sharing machine resources')
            self.creating=True
        try:
            from ihm.assembly.embodied import EmbodiedRuntime
            ident=uuid.uuid4().hex;output=self.root/'data/derived/embodied-sessions'/ident
            options={'controller':controller,'environment':environment,'regional_skin':regional_skin,'intake_mass':intake_mass,'source_pin':source_pin,**({'environment_selection':environment_selection,**environment_options} if environment_selection is not None else {})}
            # The active default is a candidate like any other: re-resolved on the
            # owner thread and disclosed, never an unreported implicit selection.
            selected=None if source_pin is None else (dict(data['ibm_candidate']) if explicit else {'commit':ACTIVE_COMMIT,'manifest_sha256':source_pin.manifest_sha256})
            def factory():
                if selected is not None and resolve_ibm_candidate(self.root,selected)!=source_pin:
                    raise ValueError('Candidate source pin changed before initialization')
                body=EmbodiedRuntime.from_workspace(self.root,output/'runtime',**options)
                try:register_body_surface_assets(body,self.surface_assets)
                except BaseException as error:
                    error.cleanup_owner=body
                    raise
                return body
            actor=BodyActor(factory,output)
            actor.controller_selection=dict(controller)
            actor.brain_source_selection=({'mode':'implicit_checkpoint','kind':controller['kind']} if source_pin is None else {'mode':'immutable_candidate' if explicit else 'active_default',
                'commit':selected['commit'],'manifest_sha256':source_pin.manifest_sha256,
                'package_sha256':source_pin.package_sha256,'neural_source_sha256':source_pin.neural_source_sha256})
            if controller['kind']=='engineering_stance':
                actor.brain_source_selection.update(mode='engineered_control',kind='engineering_stance',cortical_motor_output_active=False)

            # Timeout must not orphan initialization or free its resource slot.
            with self.lock:
                self.actors[ident]=actor
                shutting_down=self.shutting_down
            if shutting_down:actor.request_close(wait=False)
            return {'id':ident,**actor.status(),**self._identity(actor)}
        finally:
            with self.creation_done:self.creating=False;self.creation_done.notify_all()
    @staticmethod
    def _identity(actor):
        value=getattr(actor,'brain_source_selection',None)
        result={'brain_source_selection':dict(value)} if value is not None else {}
        controller=getattr(actor,'controller_selection',None)
        if controller is not None:result['controller_selection']=dict(controller)
        return result
    def command(self,ident,action,data=None):
        with self.lock:actor=self.actors.get(ident)
        if action=='close' and data not in ({},None):raise ValueError('Close requires an empty object')
        if actor is None:
            if action=='close':return {'id':ident,'closed':True,'already_absent':True}
            raise ValueError('Unknown embodied session')
        status=actor.status()
        if action=='snapshot' and status['status']!='ready':return {'id':ident,**status,**self._identity(actor)}
        result=actor.request_close() if action=='close' else actor.call(action,data)
        if action=='close' and result.get('closed'):
            with self.lock:self.actors.pop(ident,None)
        return {'id':ident,**result,**self._identity(actor)}
    def list(self):
        with self.lock:return {'sessions':[{'id':ident,**actor.status(),**self._identity(actor)} for ident,actor in self.actors.items()]}
    def close(self,timeout=30):
        deadline=time.monotonic()+timeout
        with self.creation_done:
            self.shutting_down=True
            self.creation_done.wait_for(lambda:not self.creating,timeout=max(0,deadline-time.monotonic()))
            creating=self.creating
            actors=tuple(self.actors.items())
        for _,actor in actors:
            try:actor.request_close(wait=False)
            except ValueError:pass
        for _,actor in actors:actor.thread.join(timeout=max(0,deadline-time.monotonic()))
        unresolved=[ident for ident,actor in actors if not actor.closed or actor.thread.is_alive()]
        if creating:unresolved.append('pending creation')
        if unresolved:raise RuntimeError('Body termination unconfirmed: '+', '.join(unresolved))
        with self.lock:
            for ident,actor in actors:
                if self.actors.get(ident) is actor:self.actors.pop(ident)
