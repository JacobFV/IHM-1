"""Signed command boundary contracts without starting a native process."""
import unittest
from types import SimpleNamespace
from ihm.native.coupled_session import SignedCoupledNativeSession

class Tests(unittest.TestCase):
    def session(self):
        session=SignedCoupledNativeSession.__new__(SignedCoupledNativeSession)
        session.config=SimpleNamespace(horizon_s=.04);session._elapsed_ticks=0;session.calls=[]
        def command(*args,**kwargs):session.calls.append((args,kwargs));return {'status':'ok'}
        session._command=command
        return session
    def test_signed_incidence_and_fixed_reference(self):
        s=self.session();s.signed_step('a'*64,-2,3,-5)
        self.assertEqual(s.calls,[(('signed_step','a'*64,-2,3,-5),{'advance_ticks':1})])
        with self.assertRaises(ValueError):s.signed_step('b'*64,-2,3,-5)
        with self.assertRaises(ValueError):s.save_state()
    def test_invalid_values_and_horizon_do_not_send(self):
        s=self.session()
        for args in [('bad',0,0,0),('a'*64,1,0,0),('a'*64,float('nan'),0,0),('a'*64,6000,6000,0)]:
            with self.assertRaises(ValueError):s.signed_step(*args)
        s._elapsed_ticks=2
        with self.assertRaises(ValueError):s.signed_step('a'*64,0,0,0)
        self.assertEqual(s.calls,[])
if __name__=='__main__':unittest.main()
