"""Source-only API candidate validation; no native startup or donor imports."""
from pathlib import Path
import json,shutil,tempfile,unittest
from unittest.mock import patch
from ihm.app.embodied import EmbodiedSessions
from scripts.verify_embodied_actor import FakeBody

ROOT=Path(__file__).resolve().parents[1]
COMMIT='398375cc993ac17907d4ddd9f6d53eed2fab31de'

class Tests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'data/derived/ibm-candidates'/COMMIT
        shutil.copytree(ROOT/'data/derived/ibm-candidates'/COMMIT,self.path)
        raw=json.loads((self.path/'source_pin.json').read_text());raw['artifact_dir']=str(self.path)
        (self.path/'source_pin.json').write_text(json.dumps(raw))
        self.selector={'commit':COMMIT,'manifest_sha256':raw['manifest_sha256']}
        self.sessions=EmbodiedSessions(self.root)

    def rejected(self,selector=None):
        with patch('ihm.app.embodied.BodyActor') as actor:
            with self.assertRaises(ValueError):self.sessions.create({'ibm_candidate':self.selector if selector is None else selector})
            actor.assert_not_called()
        self.assertFalse((self.root/'data/derived/embodied-sessions').exists())

    def test_candidate_forwarded_and_identity_retained(self):
        with patch('ihm.assembly.embodied.EmbodiedRuntime.from_workspace',side_effect=lambda *a,**kw:FakeBody()) as factory:
            try:
                result=self.sessions.create({'ibm_candidate':self.selector,'regional_skin':True})
                self.sessions.actors[result['id']].ready.result(2)
                pin=factory.call_args.kwargs['source_pin']
                self.assertEqual(pin.artifact_dir,self.path)
                self.assertEqual(result['brain_source_selection']['manifest_sha256'],self.selector['manifest_sha256'])
                for response in (self.sessions.list()['sessions'][0],self.sessions.command(result['id'],'snapshot')):
                    self.assertEqual(response['brain_source_selection'],result['brain_source_selection'])
                self.assertEqual(result['brain_source_selection']['package_sha256'],pin.package_sha256)
            finally:self.sessions.close(timeout=2)

    def test_strict_selector(self):
        for value in (None,{},'path',{'commit':'../outside','manifest_sha256':'a'*64},dict(self.selector,path='/tmp'),dict(self.selector,commit=COMMIT.upper()),dict(self.selector,manifest_sha256='0'*64)):
            with self.subTest(value=value):self.rejected(value) if value is not None else self.rejected_explicit_null()

    def rejected_explicit_null(self):
        with patch('ihm.app.embodied.BodyActor') as actor:
            with self.assertRaises(ValueError):self.sessions.create({'ibm_candidate':None})
            actor.assert_not_called()

    def test_extra_stale_and_symlink_files(self):
        extra=self.path/'unexpected';extra.write_text('extra');self.rejected();extra.unlink()
        extra=self.path/'source/empty';extra.mkdir();self.rejected();extra.rmdir()
        source=self.path/'source/ibm/__init__.py';original=source.read_bytes()
        source.write_bytes(original+b'\n# changed');self.rejected();source.write_bytes(original)
        outside=self.root/'outside.py';outside.write_bytes(original);source.unlink();source.symlink_to(outside);self.rejected()

    def test_symlink_candidate_and_pin_redirect(self):
        raw=json.loads((self.path/'source_pin.json').read_text());raw['artifact_dir']=str(self.root/'outside')
        (self.path/'source_pin.json').write_text(json.dumps(raw));self.rejected()
        outside=self.root/'moved';self.path.rename(outside);self.path.symlink_to(outside,target_is_directory=True);self.rejected()

    def test_changed_candidate_rejected_again_on_owner_thread(self):
        from ihm.app.embodied import BodyActor
        def delayed(factory,output):
            (self.path/'extra').write_text('arrived after request validation')
            return BodyActor(factory,output)
        with patch('ihm.app.embodied.BodyActor',side_effect=delayed),patch('ihm.assembly.embodied.EmbodiedRuntime.from_workspace') as native:
            try:
                result=self.sessions.create({'ibm_candidate':self.selector})
                with self.assertRaises(ValueError):self.sessions.actors[result['id']].ready.result(2)
                native.assert_not_called()
            finally:self.sessions.close(timeout=2)

    def test_default_omits_candidate_pin(self):
        with patch('ihm.assembly.embodied.EmbodiedRuntime.from_workspace',side_effect=lambda *a,**kw:FakeBody()) as factory:
            try:
                result=self.sessions.create({});self.sessions.actors[result['id']].ready.result(2)
                self.assertNotIn('source_pin',factory.call_args.kwargs)
                self.assertNotIn('brain_source_selection',result)
            finally:self.sessions.close(timeout=2)

if __name__=='__main__':unittest.main()
