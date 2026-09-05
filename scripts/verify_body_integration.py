"""Exercise real reduced components on a short retained native fixture."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

class IntegrationTests(unittest.TestCase):
    def make(self):
        from ihm.assembly.body import CanonicalBody, read_native
        directory=ROOT/'data/derived/canonical'
        payload=json.loads((directory/'body.json').read_text())
        assets={name:json.loads((directory/(name+'.json')).read_text()) for name in ['anatomy','profile','mechanics','brain','respiration','peripheral']}
        body=CanonicalBody(ROOT,payload,assets)
        series=read_native(directory/'native_baseline_v1')
        for key in ['time_s','physiology','compartments']:series[key]=series[key][:100]
        return body,series

    def test_chest_skin_and_peripheral_share_the_body_clock(self):
        body,series=self.make()
        with tempfile.TemporaryDirectory() as tmp, patch('ihm.assembly.body.read_native',return_value=series):
            result=body.simulate('fixture',Path(tmp)/'run.json')
        for frame in result['frames']:
            self.assertAlmostEqual(frame['respiration']['time_s'],frame['time_s'])
            self.assertAlmostEqual(frame['peripheral']['time_s'],frame['time_s'])
            self.assertLess(abs(frame['respiration']['audit']['energy_balance_residual_j']),1e-6)
            self.assertFalse(any(frame['peripheral']['motor_activations'].values()))
        self.assertGreater(np.ptp([f['respiration']['diaphragm_descent_m'] for f in result['frames']]),.001)
        moving=set().union(*(set(f['entities']) for f in result['frames']))
        self.assertGreater(len(moving),7)
        self.assertEqual(result['clock']['coupling_mode'],'native_replay')

    def test_explicit_motor_command_changes_force_after_nerve_delay(self):
        body,series=self.make()
        muscle='body-muscle-opensim-tibant_l'
        protocol=[{'start_s':.1,'end_s':.4,'motor_commands':{muscle:.5}}]
        with tempfile.TemporaryDirectory() as tmp, patch('ihm.assembly.body.read_native',return_value=series):
            result=body.simulate('fixture',Path(tmp)/'run.json',body_interventions=protocol)
        self.assertEqual(result['body_interventions'],protocol)
        self.assertTrue(all(f['peripheral']['motor_activations'].get(muscle,0)==0 for f in result['frames'] if f['time_s']<=.1))
        self.assertGreater(max(f['muscle_forces_n'].get(muscle,0) for f in result['frames']),1.)

    def test_failed_body_step_restores_every_clock_and_signal_queue(self):
        from ihm.assembly.body_runtime import BodyRuntime
        body,series=self.make()
        runtime=BodyRuntime(ROOT,body.assets,body.payload['volume_bindings'],series['initial_compartments'])
        before=runtime.checkpoint()
        with self.assertRaises(ValueError):
            runtime.step(.02,series['physiology'][0],series['compartments'][0],{'stimuli':{'missing':{'pressure_pa':100}},'motor_commands':{},'blocked_nerves':[]})
        after=runtime.checkpoint()
        self.assertFalse(hasattr(runtime.mechanics,'orientation_reactions'))
        self.assertEqual(before['time_s'],after['time_s'])
        for name,state in before['engines'].items():
            for key,value in state.items():
                if isinstance(value,np.ndarray):np.testing.assert_array_equal(value,after['engines'][name][key])
                else:self.assertEqual(value,after['engines'][name][key])

    def test_submillisecond_intervention_boundary_is_executable(self):
        body,series=self.make()
        protocol=[{'start_s':.10001,'end_s':.10009,'stimuli':{'peripheral-skin-left-palm':{'pressure_pa':1000.}}}]
        with tempfile.TemporaryDirectory() as tmp, patch('ihm.assembly.body.read_native',return_value=series):
            result=body.simulate('fixture',Path(tmp)/'run.json',body_interventions=protocol)
        self.assertEqual(result['body_interventions'],protocol)

    def test_invalid_protocol_does_not_write_output(self):
        body,series=self.make()
        with tempfile.TemporaryDirectory() as tmp, patch('ihm.assembly.body.read_native',return_value=series):
            path=Path(tmp)/'bad.json'
            with self.assertRaises(ValueError):body.simulate('fixture',path,body_interventions=[{'start_s':0,'end_s':1,'motor_commands':{'unknown':1}}])
            self.assertFalse(path.exists())

if __name__=='__main__':unittest.main()
