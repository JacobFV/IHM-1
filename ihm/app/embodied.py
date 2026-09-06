"""Demand-stepped native body actors: one controller thread per live body."""
from concurrent.futures import Future
from pathlib import Path
import gzip,hashlib,json,os,queue,threading,time,uuid


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
    def create(self,data):
        if not isinstance(data,dict) or set(data)-{'environment'}:raise ValueError('Unknown embodied configuration')
        environment=data.get('environment','supine')
        if environment not in ('free','supine','upright'):raise ValueError('Unknown articulated environment')
        with self.lock:
            if self.shutting_down:raise RuntimeError('Embodied service is shutting down')
            if self.creating or any(not a.closed for a in self.actors.values()):raise ValueError('One native body at a time while sharing machine resources')
            self.creating=True
        try:
            from ihm.assembly.embodied import EmbodiedRuntime
            ident=uuid.uuid4().hex;output=self.root/'data/derived/embodied-sessions'/ident
            actor=BodyActor(lambda:EmbodiedRuntime.from_workspace(self.root,output/'runtime',environment=environment),output)
            # Timeout must not orphan initialization or free its resource slot.
            with self.lock:
                self.actors[ident]=actor
                shutting_down=self.shutting_down
            if shutting_down:actor.request_close(wait=False)
            return {'id':ident,**actor.status()}
        finally:
            with self.creation_done:self.creating=False;self.creation_done.notify_all()
    def command(self,ident,action,data=None):
        with self.lock:actor=self.actors.get(ident)
        if action=='close' and data not in ({},None):raise ValueError('Close requires an empty object')
        if actor is None:
            if action=='close':return {'id':ident,'closed':True,'already_absent':True}
            raise ValueError('Unknown embodied session')
        status=actor.status()
        if action=='snapshot' and status['status']!='ready':return {'id':ident,**status}
        result=actor.request_close() if action=='close' else actor.call(action,data)
        if action=='close' and result.get('closed'):
            with self.lock:self.actors.pop(ident,None)
        return {'id':ident,**result}
    def list(self):
        with self.lock:return {'sessions':[{'id':ident,**actor.status()} for ident,actor in self.actors.items()]}
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
