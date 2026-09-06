"""Regional boundary protocol tests without native process execution."""
import unittest
from types import SimpleNamespace

class Tests(unittest.TestCase):
    def session(self):
        from ihm.native.regional_session import RegionalSignedNativeSession
        s=RegionalSignedNativeSession.__new__(RegionalSignedNativeSession)
        s.config=SimpleNamespace(horizon_s=.1);s._elapsed_ticks=0;s.calls=[]
        def command(*args,**kwargs):s.calls.append((args,kwargs));return {'status':'ok','elapsed_s':0}
        s._command=command
        return s
    def test_regional_pressure_is_validated_without_advancing(self):
        s=self.session();s.regional_skin_pressure('region_a',133.322)
        self.assertEqual(s.calls,[(('regional_skin_pressure','region_a',133.322),{})])
        for region,pressure in [('Skin',1),('region_a',-1),('region_b',5001),('residual',float('nan')),('region_a',True)]:
            with self.assertRaises(ValueError):s.regional_skin_pressure(region,pressure)
        self.assertEqual(len(s.calls),1)
    def test_signed_guard_and_topology_restrictions(self):
        s=self.session();s.signed_step('a'*64,.01,.02,-.01)
        for action in (lambda:s.step(.02),lambda:s.skin_compression(0),s.save_state):
            with self.assertRaises(ValueError):action()
        self.assertEqual(len(s.calls),1)
    def test_variant_and_adapter_are_separate(self):
        from unittest.mock import patch
        from ihm.native.regional_session import RegionalSignedNativeSession
        from ihm.native.session import VARIANTS
        self.assertIn('whole_body_integrity_regional_skin_graph_v2',VARIANTS)
        with patch('ihm.native.session.NativeSession.__init__',return_value=None) as initialize:
            RegionalSignedNativeSession(SimpleNamespace(engine_variant='whole_body_integrity_regional_skin_graph_v2',state_path='fixture.xml'),'unused')
            self.assertEqual(initialize.call_count,1)
        with self.assertRaises(ValueError):RegionalSignedNativeSession(SimpleNamespace(engine_variant='whole_body_integrity_gi_absorption'),'unused')
        with self.assertRaises(ValueError):RegionalSignedNativeSession(SimpleNamespace(engine_variant='whole_body_integrity_regional_skin_graph_v2',state_path=None),'unused')
        self.assertEqual(RegionalSignedNativeSession.executable_name,'native_biogears_regional_signed')

if __name__=='__main__':unittest.main()
