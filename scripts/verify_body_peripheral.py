"""Causal somatic conduction, canonical binding and executable coupling checks."""
from pathlib import Path
import importlib.util
import json
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

class PeripheralVerification(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('ihm.assembly.peripheral'), 'Missing executable peripheral model')
        from ihm.assembly.peripheral import BodyPeripheral
        self.data=json.loads((ROOT/'data/derived/canonical/peripheral.json').read_text())
        self.make=lambda:BodyPeripheral.from_dict(self.data)
        self.patch=next(p for p in self.data['receptor_patches'] if p['id']=='peripheral-skin-left-palm')

    def test_quiet_and_delayed_lateral_sensation(self):
        quiet=self.make(); p=self.make()
        self.assertFalse(any(quiet.step(.02)['brain_inputs_hz'].values()))
        stimulus={self.patch['id']:{'pressure_pa':20000.}}
        first=p.step(.001,stimulus)
        self.assertGreater(sum(first['receptor_rates_hz'].values()),0)
        self.assertFalse(any(first['brain_inputs_hz'].values()))
        for _ in range(100):out=p.step(.005,stimulus)
        self.assertGreater(out['brain_inputs_hz']['brain-rh-postcentral'],0)
        self.assertEqual(out['brain_inputs_hz'].get('brain-lh-postcentral',0),0)
        self.assertFalse(any(out['motor_activations'].values()))
        for _ in range(300):out=p.step(.01)
        self.assertLess(sum(out['brain_inputs_hz'].values()),1e-8)

    def test_sensory_input_changes_real_brain_ode(self):
        from ihm.assembly.brain import BodyBrain
        d=json.loads((ROOT/'data/derived/canonical/brain.json').read_text())
        a,b=BodyBrain(d),BodyBrain(d)
        aa=a.step(.1);bb=b.step(.1,sensory_inputs_hz={'brain-rh-postcentral':100.})
        i=a.ids.index('brain-rh-postcentral')
        self.assertGreater(bb['regional_state']['activity_hz'][i],aa['regional_state']['activity_hz'][i])
        with self.assertRaises(ValueError):b.step(.02,sensory_inputs_hz={'not-a-region':2.})

    def test_motor_drives_mechanics_after_delay_and_feedback(self):
        from ihm.assembly.mechanics import BodyMechanics
        d=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
        binding=next(m for m in self.data['muscle_bindings'] if m['muscle_id']=='body-muscle-opensim-tibant_l')
        p=self.make();command={'motor_commands':{binding['muscle_id']:.5}}
        self.assertFalse(any(p.step(.001,brain_state=command)['motor_activations'].values()))
        for _ in range(30):out=p.step(.005,brain_state=command)
        self.assertGreater(out['motor_activations'][binding['muscle_id']],.1)
        self.assertEqual(out['motor_activations'].get('body-muscle-opensim-tibant_r',0),0)
        mechanics=BodyMechanics(d)
        state=mechanics.step(.01,{'activation':out['motor_activations']})
        self.assertGreater(state['muscle_forces_n'][binding['muscle_id']],1.)
        for _ in range(40):feedback=p.step(.01,mechanical_state=state)
        self.assertGreater(sum(feedback['proprioceptor_rates_hz'].values()),0)

    def test_touch_does_not_invent_a_motor_policy(self):
        from ihm.assembly.brain import BodyBrain
        brain=BodyBrain(json.loads((ROOT/'data/derived/canonical/brain.json').read_text()))
        p=self.make();neural={}
        for _ in range(100):
            out=p.step(.01,{self.patch['id']:{'pressure_pa':20000.}},brain_state=neural)
            neural=brain.step(.01,sensory_inputs_hz=out['brain_inputs_hz'])
        active={k:v for k,v in out['motor_activations'].items() if v>1e-6}
        self.assertFalse(active, 'Sensory brain activity is not a descending motor command')
        self.assertEqual(out['motor_activations'].get('body-muscle-opensim-tibant_l',0),0)
        self.assertTrue(any(neural['regional_motor_drive_hz'].values()))

    def test_nerve_block_selectively_stops_arrivals_and_motor_drive(self):
        p=self.make()
        binding=next(m for m in self.data['muscle_bindings'] if m['muscle_id']=='body-muscle-opensim-tibant_l')
        other=next(m for m in self.data['muscle_bindings'] if m['muscle_id']=='body-muscle-opensim-tibant_r')
        commands={'motor_commands':{binding['muscle_id']:.5,other['muscle_id']:.5}}
        for _ in range(50):out=p.step(.01,brain_state=commands,blocked_nerves=[binding['nerve_id']])
        self.assertEqual(out['motor_activations'].get(binding['muscle_id'],0),0)
        self.assertGreater(out['motor_activations'][other['muscle_id']],.4)
        with self.assertRaises(ValueError):p.step(.01,blocked_nerves=['unknown'])

    def test_thermal_channels_and_invalid_input_transaction(self):
        p=self.make();s={self.patch['id']:{'temperature_C':42.}}
        for _ in range(400):out=p.step(.01,s)
        self.assertGreater(out['brain_inputs_hz']['brain-rh-postcentral'],0)
        t=p.time_s
        for bad in [{'unknown':{'pressure_pa':1}}, {self.patch['id']:{'temperature_C':float('nan')}}, {self.patch['id']:{'pressure_pa':-1}}]:
            with self.assertRaises(ValueError):p.step(.01,bad)
            self.assertEqual(t,p.time_s)
        with self.assertRaises(ValueError):p.step(2.)
        with self.assertRaises(ValueError):p.step(.01,brain_state={'motor_commands':{'missing':.1}})

    def test_artifact_identity_and_prior_labels(self):
        anatomy=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
        mechanics=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
        ids={e['id'] for e in anatomy['entities']};muscles={m['id'] for m in mechanics['muscles']}
        self.assertGreater(len(self.data['muscle_bindings']),100)
        self.assertGreaterEqual(len(self.data['receptor_patches']),16)
        for b in self.data['muscle_bindings']:
            self.assertIn(b['muscle_id'],muscles);self.assertIn(b['canonical_entity_id'],ids)
            self.assertGreater(b['motor_delay_s'],0);self.assertEqual(b['evidence_kind'],'named_anatomical_prior')
        for p in self.data['receptor_patches']:self.assertIn(p['body_entity_id'],ids)
        self.assertFalse(self.data['biological_validation'])

if __name__=='__main__':unittest.main()
