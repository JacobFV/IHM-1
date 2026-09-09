"""Controller choices validate before owner startup and preserve model selection."""
import tempfile
import unittest
from unittest.mock import patch

from ihm.app.embodied import EmbodiedSessions
from ihm.assembly.controller_selection import resolve_controller
from scripts.verify_embodied_actor import FakeBody


class Tests(unittest.TestCase):
    def test_strict_choices(self):
        for bad in ([], 'implicit', {'kind': 'unknown'}, {'path': '/tmp/model'},
                    {'sever': 1}, {'no_cord': 'false'}, {'kind': 'regional', 'sever': True}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                resolve_controller(bad)
        self.assertEqual(resolve_controller(), {'kind': 'regional', 'sever': False, 'no_cord': False})
        self.assertEqual(resolve_controller({'kind':'engineering_stance'}),{'kind':'engineering_stance','sever':False,'no_cord':False})
        for kind in ('implicit_cortical_stance','implicit_curriculum16_stance'):
            self.assertTrue(resolve_controller({'kind':kind,'sever':True})['sever'])
            for flag in ('no_cord','target_rad'):
                with self.subTest(kind=kind,flag=flag),self.assertRaises(ValueError):
                    resolve_controller({'kind':kind,flag:True})
        for flag in ('sever','no_cord','target_rad'):
            with self.assertRaises(ValueError):resolve_controller({'kind':'engineering_stance',flag:True})
        # The retained-kernel raw path keeps both kernel and cord ablations.
        self.assertEqual(resolve_controller({'kind':'implicit_curriculum16','sever':True,'no_cord':True}),
                         {'kind':'implicit_curriculum16','sever':True,'no_cord':True})
        with self.assertRaises(ValueError):resolve_controller({'kind':'implicit_curriculum16','target_rad':.1})

    def test_stance_rejects_unidentified_environment_and_mass_before_owner(self):
        with tempfile.TemporaryDirectory() as root:
            for kind,config in ((kind,config) for kind in ('engineering_stance','implicit_cortical_stance','implicit_curriculum16_stance')
                    for config in ({},{'environment':'supine'},{'environment':'upright','intake_mass':True})):
                with patch('ihm.app.embodied.BodyActor') as actor:
                    with self.assertRaisesRegex(ValueError,'upright fixed-mass'):
                        EmbodiedSessions(root).create({'controller':{'kind':kind},**config})
                    actor.assert_not_called()

    def test_implicit_does_not_resolve_regional_candidate(self):
        with tempfile.TemporaryDirectory() as root:
            sessions = EmbodiedSessions(root)
            selection = {'kind': 'implicit', 'sever': True, 'no_cord': False}
            with patch('ihm.assembly.embodied.EmbodiedRuntime.from_workspace', side_effect=lambda *a, **kw: FakeBody()) as factory:
                try:
                    result = sessions.create({'controller': selection})
                    sessions.actors[result['id']].ready.result(2)
                    self.assertEqual(factory.call_args.kwargs['controller'], selection)
                    self.assertIsNone(factory.call_args.kwargs['source_pin'])
                    self.assertEqual(result['controller_selection'], selection)
                    self.assertEqual(result['brain_source_selection']['mode'], 'implicit_checkpoint')
                finally:
                    sessions.close(timeout=2)

    def test_conflicting_brains_rejected_before_owner(self):
        with tempfile.TemporaryDirectory() as root:
            with patch('ihm.app.embodied.BodyActor') as actor:
                with self.assertRaisesRegex(ValueError, 'cannot select'):
                    EmbodiedSessions(root).create({'controller': {'kind': 'implicit'}, 'ibm_candidate': {}})
                actor.assert_not_called()


if __name__ == '__main__':
    unittest.main()
