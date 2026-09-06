"""Signed command boundary contracts without starting a native process."""
import unittest
from types import SimpleNamespace
from ihm.native.coupled_session import SignedCoupledNativeSession

class Tests(unittest.TestCase):
    def test_verified_gi_descendant_is_an_explicit_signed_variant(self):
        from unittest.mock import patch
        from ihm.native.session import VARIANTS
        variant='whole_body_integrity_gi_absorption'
        self.assertIn(variant,VARIANTS)
        with patch('ihm.native.session.NativeSession.__init__',return_value=None) as initialize:
            SignedCoupledNativeSession(SimpleNamespace(engine_variant=variant),'unused')
            self.assertEqual(initialize.call_count,1)

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
    def test_unsigned_step_allowed_only_before_reference_binding(self):
        s=self.session()
        s.step(.02)
        self.assertEqual(s.calls,[(('step',1),{'advance_ticks':1})])
        s.signed_step('a'*64,0,0,0)
        calls=list(s.calls)
        with self.assertRaisesRegex(ValueError,'signed_step'):
            s.step(.02)
        self.assertEqual(s.calls,calls)
        s.signed_step('a'*64,0,0,0)
        self.assertEqual(len(s.calls),len(calls)+1)
    def test_invalid_values_and_horizon_do_not_send(self):
        s=self.session()
        for args in [('bad',0,0,0),('a'*64,1,0,0),('a'*64,float('nan'),0,0),('a'*64,6000,6000,0)]:
            with self.assertRaises(ValueError):s.signed_step(*args)
        s._elapsed_ticks=2
        with self.assertRaises(ValueError):s.signed_step('a'*64,0,0,0)
        self.assertEqual(s.calls,[])
def verify_native():
    """Bypass the Python guard to exercise the native command rejection itself."""
    import hashlib,json,math,tempfile,time
    from pathlib import Path
    from ihm.native.session import SessionConfig
    root=Path(__file__).resolve().parents[1]
    out=Path(tempfile.mkdtemp(prefix='signed-command-guard-',dir=root/'data/derived/audits'))
    source=json.loads((root/'data/derived/systemic/exertion_v3/exercise/native/manifest.json').read_text())
    config=SessionConfig(state_path=source['configuration']['state_path'],engine_variant='whole_body_integrity_signed_muscle_v2',horizon_s=.08)
    ref=hashlib.sha256(b'protocol guard fixture; frozen identity without reference calibration').hexdigest()
    session=None;started=time.monotonic()
    try:
        session=SignedCoupledNativeSession(config,out/'native')
        ordinary=session.step(.02)
        signed=session.signed_step(ref,.01,.01,0)
        before_ticks=session._elapsed_ticks
        with unittest.TestCase().assertRaisesRegex(ValueError,'signed_step'):
            session.step(.02)
        assert session._elapsed_ticks==before_ticks
        # Exercise the wire command directly; Python's public guard cannot hide a native defect.
        rejected=session._command('step',1,status='rejected')
        assert rejected['time_s']==signed['time_s'] and rejected['elapsed_s']==signed['elapsed_s']
        assert rejected['values']==signed['values']
        assert session._elapsed_ticks==before_ticks
        continued=session.signed_step(ref,-.01,-.01,0)
        assert math.isclose(continued['elapsed_s'],.06,abs_tol=1e-12)
        observations={}
        for label,frame,previous,delta in [('positive',signed,ordinary,.01),('negative',continued,signed,-.01)]:
            values=frame['values'];readers={}
            for reader in ('cardiovascular','endocrine','nervous_metabolic_fraction','nervous_tmr'):
                key='coupling.effective_reader_'+reader
                count=values[key+'_count'];assert count>=0 and count==int(count)
                assert (key+'_w' in values)==(count>0)
                readers[reader]={'count':count,'watts':values.get(key+'_w')}
            assert readers['endocrine']['count']>0
            # Native MetabolicToneResponse is conditional on drug-driven resistance
            # change (Cardiovascular.cpp:2080-2094); an idle consumer must stay absent.
            if readers['cardiovascular']['count']:
                assert math.isclose(readers['cardiovascular']['watts'],previous['values']['metabolic_rate_w']+delta,abs_tol=1e-10)
            assert math.isclose(readers['endocrine']['watts'],values['coupling.native_heat_w']+delta,abs_tol=1e-10)
            for reader in ('nervous_metabolic_fraction','nervous_tmr'):
                if readers[reader]['count']:
                    assert math.isclose(readers[reader]['watts'],values['coupling.native_heat_w']+delta,abs_tol=1e-10)
            observations[label]=readers
        session.close();session=None
        report={'passed':True,'native_ordinary_step_rejected_after_binding':True,'rejected_step_state_exact':True,
                'signed_continuation_elapsed_s':continued['elapsed_s'],'effective_readers':observations,
                'wall_s':time.monotonic()-started,'scope':'Command rejection and native effective-reader phase observations, no rollback or general stress validation claim'}
        (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'passed':True,'output':str(out)},indent=2))
    except BaseException as error:
        (out/'failure.json').write_text(json.dumps({'error':str(error),'wall_s':time.monotonic()-started},indent=2)+'\n')
        raise
    finally:
        if session:session.close(graceful=False)

if __name__=='__main__':
    import sys
    if sys.argv[1:]==['--native']:verify_native()
    else:unittest.main()
