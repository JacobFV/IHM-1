"""Concurrent HTTP callers share one native-controller thread and command clock."""
from pathlib import Path
from concurrent.futures import Future
from unittest.mock import patch
import tempfile,threading,unittest
from ihm.app.embodied import BodyActor,EmbodiedSessions

class FakeBody:
    def __init__(self):self.owner=threading.get_ident();self.sequence=0;self.failed=False
    def snapshot(self):return {'sequence':self.sequence,'time_s':self.sequence*.02,'entities':{}}
    def step(self,data):
        assert threading.get_ident()==self.owner;self.sequence+=1;return self.snapshot()
    def schedule_intakes(self,data):
        assert threading.get_ident()==self.owner;self.sequence+=1;return self.snapshot()
    def close(self):assert threading.get_ident()==self.owner

class Tests(unittest.TestCase):
    def test_intake_mass_configuration_is_explicit_bool_and_forwarded(self):
        with tempfile.TemporaryDirectory() as p:
            sessions=EmbodiedSessions(p)
            with patch('ihm.app.embodied.BodyActor') as actor:
                for value in (1,None,'true',[]):
                    with self.subTest(value=value),self.assertRaisesRegex(ValueError,'intake_mass must be boolean'):
                        sessions.create({'intake_mass':value})
                actor.assert_not_called()
            with patch('ihm.assembly.embodied.EmbodiedRuntime.from_workspace',side_effect=lambda *args,**kwargs:FakeBody()) as factory:
                try:
                    result=sessions.create({'intake_mass':True,'regional_skin':True})
                    sessions.actors[result['id']].ready.result(2)
                    self.assertIs(factory.call_args.kwargs['intake_mass'],True)
                    self.assertIs(factory.call_args.kwargs['regional_skin'],True)
                finally:sessions.close(timeout=2)

    def test_regional_configuration_validates_before_starting_owner(self):
        with tempfile.TemporaryDirectory() as p:
            sessions=EmbodiedSessions(p)
            with patch('ihm.app.embodied.BodyActor') as actor:
                for value in (1,0,None,'true',{},[]):
                    with self.subTest(value=value),self.assertRaisesRegex(ValueError,'regional_skin must be boolean'):
                        sessions.create({'regional_skin':value})
                actor.assert_not_called()
            self.assertEqual(sessions.list()['sessions'],[])

    def test_regional_configuration_reaches_single_owner_factory(self):
        for config,expected in (({},False),({'regional_skin':True},True)):
            with self.subTest(config=config),tempfile.TemporaryDirectory() as p:
                sessions=EmbodiedSessions(p)
                with patch('ihm.assembly.embodied.EmbodiedRuntime.from_workspace',side_effect=lambda *args,**kwargs:FakeBody()) as factory:
                    try:
                        result=sessions.create(config)
                        sessions.actors[result['id']].ready.result(2)
                        self.assertEqual(factory.call_args.kwargs['regional_skin'],expected)
                        self.assertEqual(factory.call_args.kwargs['environment'],'supine')
                        with self.assertRaisesRegex(ValueError,'One native body'):sessions.create(config)
                    finally:sessions.close(timeout=2)

    def test_intake_schedule_uses_owner_sequence_and_distinct_journal_event(self):
        import gzip,json
        with tempfile.TemporaryDirectory() as p:
            actor=BodyActor(FakeBody,Path(p))
            try:
                result=actor.call('intakes',{'sequence':0,'events':[]})
                self.assertEqual(result['sequence'],1)
                with self.assertRaises(ValueError):actor.call('intakes',{'sequence':0,'events':[]})
                kinds=[json.loads(gzip.decompress(f.read_bytes()))['kind'] for f in sorted(actor.events.glob('*.gz'))]
                self.assertIn('intake_schedule',kinds)
            finally:actor.request_close()
    def test_shutdown_tracks_creation_before_actor_registration(self):
        with tempfile.TemporaryDirectory() as p:
            entered=threading.Event();release=threading.Event();result=Future()
            sessions=EmbodiedSessions(p)
            def delayed_actor(factory,output):
                entered.set();release.wait(2);return BodyActor(FakeBody,output)
            def create():
                try:result.set_result(sessions.create({}))
                except BaseException as error:result.set_exception(error)
            with patch('ihm.app.embodied.BodyActor',side_effect=delayed_actor):
                creator=threading.Thread(target=create);creator.start()
                try:
                    self.assertTrue(entered.wait(2))
                    with self.assertRaisesRegex(RuntimeError,'pending creation'):sessions.close(timeout=.01)
                finally:
                    release.set();creator.join(2);result.result(2);sessions.close(timeout=2)
            self.assertEqual(sessions.list()['sessions'],[])

    def test_shutdown_prevents_new_native_startup(self):
        with tempfile.TemporaryDirectory() as p:
            sessions=EmbodiedSessions(p);sessions.close()
            with self.assertRaisesRegex(RuntimeError,'shutting down'):sessions.create({})

    def test_cancellation_is_not_lost_when_command_queue_is_full(self):
        with tempfile.TemporaryDirectory() as p:
            entered=threading.Event();release=threading.Event()
            class SlowBody(FakeBody):
                def step(self,data):entered.set();release.wait(2);return super().step(data)
            actor=BodyActor(SlowBody,Path(p)/'actor');actor.ready.result(2)
            actor.call('step',{'sequence':0},wait=False)
            self.assertTrue(entered.wait(2))
            for _ in range(8):actor.pending.put_nowait(('snapshot',None,Future()))
            with self.assertRaisesRegex(ValueError,'queue is full'):actor.request_close(wait=False)
            release.set();actor.thread.join(.1)
            try:self.assertTrue(actor.closed)
            finally:
                actor.request_close();actor.thread.join(2)

    def test_partial_factory_owner_can_retry_failed_cleanup(self):
        with tempfile.TemporaryDirectory() as p:
            attempted=threading.Event()
            class PartialOwner(FakeBody):
                allow_close=False
                def close(self):
                    super().close();attempted.set()
                    if not self.allow_close:raise RuntimeError('partial cleanup failed')
            def factory():
                owner=PartialOwner()
                error=RuntimeError('startup failed');error.cleanup_owner=owner
                raise error
            actor=BodyActor(factory,Path(p)/'actor')
            with self.assertRaisesRegex(RuntimeError,'startup failed'):actor.ready.result(2)
            self.assertTrue(attempted.wait(2))
            try:
                self.assertFalse(actor.closed)
                with self.assertRaisesRegex(RuntimeError,'partial cleanup failed'):actor.request_close()
            finally:
                actor.body.allow_close=True;actor.request_close();actor.thread.join(2)
            self.assertTrue(actor.closed)
            self.assertEqual(actor.status()['error'],'startup failed')

    def test_pending_cancellation_cleanup_failure_keeps_owner(self):
        with tempfile.TemporaryDirectory() as p:
            gate=threading.Event();attempted=threading.Event()
            class RetryBody(FakeBody):
                allow_close=False
                def close(self):
                    super().close();attempted.set()
                    if not self.allow_close:raise RuntimeError('cancel cleanup failed')
            def factory():gate.wait(2);return RetryBody()
            actor=BodyActor(factory,Path(p)/'actor');actor.request_close();gate.set()
            self.assertTrue(attempted.wait(2))
            try:
                self.assertFalse(actor.closed)
                with self.assertRaisesRegex(RuntimeError,'cancel cleanup failed'):actor.request_close()
            finally:
                actor.body.allow_close=True;actor.request_close();actor.thread.join(2)
            self.assertFalse(actor.thread.is_alive())

    def test_shutdown_waits_for_initialization_and_reports_timeout(self):
        with tempfile.TemporaryDirectory() as p:
            gate=threading.Event()
            def factory():gate.wait(2);return FakeBody()
            actor=BodyActor(factory,Path(p)/'actor')
            sessions=EmbodiedSessions(p);sessions.actors['a']=actor
            try:
                with self.assertRaisesRegex(RuntimeError,'termination unconfirmed'):sessions.close(timeout=.01)
                self.assertFalse(actor.closed)
            finally:
                gate.set();sessions.close(timeout=2)
            self.assertFalse(actor.thread.is_alive())

    def test_failed_cleanup_retains_slot_and_retries_on_owner(self):
        with tempfile.TemporaryDirectory() as p:
            class RetryBody(FakeBody):
                allow_close=False
                def close(self):
                    super().close()
                    if not self.allow_close:raise RuntimeError('termination unconfirmed')
            actor=BodyActor(RetryBody,Path(p)/'actor');actor.ready.result(2)
            sessions=EmbodiedSessions(p);sessions.actors['a']=actor
            try:
                with self.assertRaisesRegex(RuntimeError,'termination unconfirmed'):sessions.command('a','close',{})
                actor.thread.join(.02)
                self.assertFalse(actor.closed)
                self.assertTrue(actor.thread.is_alive())
                self.assertIn('termination unconfirmed',actor.status()['error'])
                with self.assertRaisesRegex(ValueError,'One native body'):sessions.create({})
            finally:
                actor.body.allow_close=True
                sessions.command('a','close',{})
            self.assertFalse(actor.thread.is_alive())
            self.assertTrue(sessions.command('a','close',{})['already_absent'])

    def test_pending_initialization_can_be_closed_without_losing_owner(self):
        with tempfile.TemporaryDirectory() as p:
            gate=threading.Event()
            def factory():gate.wait(timeout=2);return FakeBody()
            actor=BodyActor(factory,Path(p))
            self.assertEqual(actor.status()['status'],'initializing')
            self.assertEqual(actor.request_close()['status'],'closing')
            self.assertFalse(actor.closed)
            gate.set();actor.thread.join(timeout=2)
            self.assertTrue(actor.closed)
            self.assertFalse(actor.thread.is_alive())
    def test_sequence_read_recovery_and_same_owner_thread(self):
        with tempfile.TemporaryDirectory() as p:
            actor=BodyActor(lambda:FakeBody(),Path(p))
            try:
                first=actor.call('step',{'sequence':0})
                self.assertEqual(first['sequence'],1)
                with self.assertRaises(ValueError):actor.call('step',{'sequence':0})
                self.assertEqual(actor.call('snapshot')['sequence'],1)
                with self.assertRaises(ValueError):actor.call('step',{'sequence':True})
            finally:actor.call('close')
            self.assertFalse(actor.thread.is_alive())

if __name__=='__main__':unittest.main()
