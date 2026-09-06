"""Bounded native receptor receipt/IBM causal bridge fixtures; no native process."""
from copy import deepcopy
from pathlib import Path
import unittest
from ihm.assembly.native_afferents import CHANNELS,NativeAfferentBridge,endpoint_receipt
from ihm.assembly.sensorimotor import SensorimotorController
from verify_sensorimotor import observation
ROOT=Path(__file__).resolve().parents[1]
IDENTITY={k:'a'*64 for k in ('library_sha256','executable_sha256','state_sha256','manifest_sha256')}
SOURCE='b'*64

def controller():return SensorimotorController.from_root(ROOT)
def allocation():return {k:{'brain-bp3d-FJ1769':.5,'brain-bp3d-FJ1831':.5} for k in CHANNELS}
def bridge(c):return NativeAfferentBridge(c.brain,native_identity=IDENTITY,source_sha256=SOURCE,origin_s=100.,allocation=allocation())
def receipt(tick=0,sequence=None):
    return endpoint_receipt({'elapsed_s':tick*.02,'origin_s':100.,'time_s':100+tick*.02,'sequence':tick if sequence is None else sequence,
        'values':{'nervous.afferent.'+k+'_hz':10+i for i,k in enumerate(CHANNELS)}},IDENTITY,SOURCE)

class Tests(unittest.TestCase):
    def test_actual_shared_brain_causality_and_blocks_without_extra_step(self):
        a,b=controller(),controller();ports,blocked=bridge(a),bridge(b)
        for tick in range(4):
            live=ports.step(.02,receipt(tick));off=blocked.step(.02,receipt(tick),sensory_blocks=CHANNELS)
            a.step(.02,observation(a.time_s),additional_sensory_inputs_hz=live['sensory_inputs_hz'])
            b.step(.02,observation(b.time_s),additional_sensory_inputs_hz=off['sensory_inputs_hz'])
            self.assertAlmostEqual(a.time_s,a.brain.time_s);self.assertAlmostEqual(a.time_s,ports.time_s)
        self.assertNotEqual(a.brain.state.tolist(),b.brain.state.tolist())
        release=blocked.step(.02,receipt(4));self.assertGreater(sum(release['sensory_inputs_hz'].values()),0)

    def test_missing_stale_tampered_wrong_owner_are_atomic(self):
        p=bridge(controller());saved=p.checkpoint()
        bad=receipt();bad['rates_hz']['chemoreceptor']+=1
        wrong=receipt();wrong['identity']['library_sha256']='c'*64
        for value in (None,bad,wrong,receipt(1)):
            with self.assertRaises(ValueError):p.step(.02,value)
            self.assertEqual(p.checkpoint(),saved)
        p.step(.02,receipt());accepted=p.checkpoint()
        with self.assertRaises(ValueError):p.step(.02,receipt())
        self.assertEqual(p.checkpoint(),accepted)
        with self.assertRaises(ValueError):endpoint_receipt({'elapsed_s':0,'time_s':100,'origin_s':100,'sequence':0,'values':{}},IDENTITY,SOURCE)

    def test_exact_replay_and_brain_source_pin_identity(self):
        c=controller();p=bridge(c);p.step(.02,receipt());saved=p.checkpoint()
        expected=p.step(.02,receipt(1));p.restore(saved);self.assertEqual(expected,p.step(.02,receipt(1)))
        c.brain.source_identity={'selection':'different source pin'}
        with self.assertRaises(ValueError):bridge(c).restore(saved)

    def test_runtime_uses_prior_endpoint_and_restores_before_native_failure(self):
        from verify_embodied_runtime import Plant,Native,Load,Exchange
        from ihm.assembly.embodied import EmbodiedRuntime
        c=controller();port=bridge(c)
        class BodyPlant(Plant):
            def snapshot(self):return {**super().snapshot(),**observation(self.t)}
        class BodyNative(Native):
            missing=False
            def snapshot(self):
                frame=super().snapshot();frame.update(sequence=round(self.t*50)+1,origin_s=100.)
                if not self.missing:frame['values'].update({'nervous.afferent.'+k+'_hz':10 for k in CHANNELS})
                return frame
        native=BodyNative();body=EmbodiedRuntime(BodyPlant(),c,native,Exchange(),Load(),afferents=port)
        frame=body.step({})
        self.assertEqual(frame['native_afferents']['applied_input']['source_endpoint_s'],0.)
        self.assertEqual(frame['native_afferents']['endpoint_receipt']['elapsed_s'],.02)
        self.assertAlmostEqual(c.time_s,c.brain.time_s)
        saved=port.checkpoint();original=body.plant.advance
        body.plant.advance=lambda *a,**k:(_ for _ in ()).throw(ValueError('pre-native plant rejection'))
        with self.assertRaises(ValueError):body.step({})
        self.assertEqual(port.checkpoint(),saved);self.assertEqual(native.t,.02)
        body.plant.advance=original
        frame=body.step({'native_sensory_blocks':list(CHANNELS)})
        self.assertEqual(sum(frame['native_afferents']['applied_input']['sensory_inputs_hz'].values()),0)
        native.missing=True
        with self.assertRaises(ValueError):body.step({})
        self.assertTrue(body.failed);self.assertTrue(native.closed)

    def test_explicit_allocation_and_blocks_reject_unknown_ids(self):
        c=controller();bad=allocation();bad['chemoreceptor']={'invented-medulla':1}
        with self.assertRaises(ValueError):NativeAfferentBridge(c.brain,native_identity=IDENTITY,source_sha256=SOURCE,origin_s=100,allocation=bad)
        with self.assertRaises(ValueError):bridge(c).step(.02,receipt(),sensory_blocks=['invented'])

if __name__=='__main__':unittest.main()
