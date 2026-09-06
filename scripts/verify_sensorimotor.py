"""Tiny deterministic sensory/motor/IBM ownership tests; no optimizer or browser."""
from copy import deepcopy
import json
from pathlib import Path
import unittest
from ihm.assembly.sensorimotor import SensorimotorController
ROOT=Path(__file__).resolve().parents[1]
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog
CATALOG=native_muscle_catalog(ROOT)

def observation(time=0., extension=1.0, load=100.):
    muscles={}
    for side in ('r','l'):
        for name in ('tibant','soleus','gasmed','gaslat'):
            muscles[f'{name}_{side}']={'fiber_length_m':.1*(extension if name=='tibant' and side=='r' else 1.),
                'optimal_fiber_length_m':.1,'tendon_force_n':load if name=='soleus' and side=='r' else 0.,
                'max_isometric_force_n':1000.,'sensor_basis':'synthetic unit fixture CE'}
    for row in CATALOG:
        muscles.setdefault(row['id'],{'fiber_length_m':.1,'optimal_fiber_length_m':.1,'tendon_force_n':0.,'max_isometric_force_n':1000.,'sensor_basis':'synthetic unit fixture CE'})
    return {'time_s':time,'muscles':muscles,'foot_contact_force_n':{'r':100.,'l':0.}}

class SensorimotorTests(unittest.TestCase):
    def test_delays_specificity_and_blocks(self):
        c=SensorimotorController.from_root(ROOT);inp=observation(extension=1.1);before=deepcopy(inp)
        a=c.step(.01,inp);self.assertEqual(inp,before)
        self.assertTrue(all(v==0 for v in a['motor_excitations'].values()))
        for _ in range(5):a=c.step(.01,observation(c.time_s,1.1))
        self.assertGreater(a['motor_excitations']['soleus_r'],a['motor_excitations']['soleus_l'])
        self.assertGreater(a['motor_excitations']['tibant_r'],a['motor_excitations']['tibant_l'])
        a=c.step(.01,observation(c.time_s,1.1),motor_blocks=['tibant_r'])
        self.assertEqual(a['motor_excitations']['tibant_r'],0)
        self.assertGreater(a['motor_excitations']['soleus_r'],0)
        for _ in range(5):a=c.step(.01,observation(c.time_s,1.1),sensory_blocks=['soleus_r'])
        self.assertAlmostEqual(a['motor_excitations']['soleus_r'],.01)
        self.assertGreater(a['motor_excitations']['tibant_r'],.3)
    def test_catalog_wide_targeting_and_held_source_laws(self):
        from ihm.assembly.sensorimotor import SensorimotorParameters
        c=SensorimotorController.from_root(ROOT,parameters=SensorimotorParameters(cortical_gain_per_hz=0,cortical_drive_per_hz=0))
        self.assertEqual(len(c.catalog),80)
        self.assertTrue(all(row['attachment_bodies'] and row['source_sha256'] for row in c.catalog))
        for _ in range(5):x=c.step(.02,observation(c.time_s,1.1,100.))
        self.assertAlmostEqual(x['motor_excitations']['soleus_r'],.01+1.2*.1)
        self.assertAlmostEqual(x['motor_excitations']['tibant_r'],.01+1.1*(1.1-.71)-.3*.1)
        self.assertEqual(x['motor_excitations']['addbrev_r'],0.)
        c=SensorimotorController.from_root(ROOT)
        for _ in range(5):x=c.step(.02,observation(c.time_s),descending={'addbrev_r':1.})
        self.assertGreater(x['motor_excitations']['addbrev_r'],0.)
        self.assertEqual(x['motor_excitations']['addbrev_l'],0.)
        self.assertEqual(x['motor_excitations']['recfem_r'],0.)
        before=c.checkpoint()
        with self.assertRaises(ValueError):c.step(.02,observation(c.time_s),physiology={'oxygen_saturation':float('nan')})
        self.assertEqual(c.checkpoint(),before)

    def test_checkpoint_and_failure_atomicity(self):
        c=SensorimotorController.from_root(ROOT)
        c.step(.02,observation());snapshot=c.checkpoint()
        a=c.step(.02,observation(c.time_s),descending={'tibant_r':.5})
        c.restore(snapshot);b=c.step(.02,observation(c.time_s),descending={'tibant_r':.5})
        self.assertEqual(a,b)
        before=c.checkpoint();bad=observation(c.time_s);bad['muscles']['soleus_r']['tendon_force_n']=float('nan')
        with self.assertRaises(ValueError):c.step(.02,bad)
        self.assertEqual(before,c.checkpoint())
        bad=deepcopy(before);bad['brain_state'][0][0]=float('nan')
        with self.assertRaises(ValueError):c.restore(bad)
        self.assertEqual(before,c.checkpoint())
        bad=deepcopy(before);bad['model_sha256']='wrong'
        with self.assertRaises(ValueError):c.restore(bad)
    def test_brain_is_causal_and_decoder_is_effectorspecific(self):
        a=SensorimotorController.from_root(ROOT);b=SensorimotorController.from_root(ROOT)
        for _ in range(8):
            x=a.step(.02,observation(a.time_s),descending={'tibant_r':1.})
            y=b.step(.02,observation(b.time_s))
        self.assertNotEqual(x['brain']['regional_state']['activity_hz'],y['brain']['regional_state']['activity_hz'])
        self.assertGreater(x['motor_excitations']['tibant_r'],y['motor_excitations']['tibant_r'])
        self.assertEqual(x['descending_drive']['tibant_l'],0.)
        self.assertFalse(x['biological_validation'])

    def test_additional_sensory_inputs_validate_and_preserve_single_clock(self):
        c=SensorimotorController.from_root(ROOT)
        for bad in ({'missing-region':1.}, {'brain-rh-postcentral':float('nan')},
                    {'brain-rh-postcentral':-1.}, {'brain-rh-postcentral':1001.},
                    {'brain-rh-postcentral':True}, []):
            before=c.checkpoint()
            with self.assertRaises(ValueError):
                c.step(.01,observation(c.time_s),additional_sensory_inputs_hz=bad)
            self.assertEqual(before,c.checkpoint())
        result=c.step(.01,observation(),additional_sensory_inputs_hz={'brain-rh-postcentral':1000.})
        self.assertEqual(c.brain.time_s,.01)
        self.assertEqual(result['additional_sensory_inputs_hz'],{'brain-rh-postcentral':1000.})
        self.assertLessEqual(result['brain_sensory_inputs_hz']['brain-rh-postcentral'],1000.)

    def test_skin_force_shared_brain_motor_path_against_blocked_control(self):
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        site={'id':'left-skin-fixture','position_m':[0.,0.,0.],'normal':[0.,0.,1.],
            'contact_area_m2':.001,'stiffness_pa_per_m':1e7,
            'sensory_region':'brain-rh-postcentral','reference_temperature_C':33.,
            'support_basis':'explicit synthetic causal fixture'}
        skin=[CutaneousFeedback(ROOT,sites=[site],recruitment_hz_per_response=1.,delay_s=.01)
              for _ in range(2)]
        controllers=[SensorimotorController.from_root(ROOT) for _ in range(2)]
        pending=[{},{}];outputs=[]
        for step in range(40):
            outputs=[]
            for i,(receptors,controller) in enumerate(zip(skin,controllers)):
                # Only previous receptor endpoint can drive this brain interval.
                output=controller.step(.01,observation(controller.time_s),
                    descending={'addbrev_l':.5},additional_sensory_inputs_hz=pending[i])
                contact=receptors.step(.01,{'time_s':receptors.time_s,
                    'contacts':[{'id':site['id'],'force_n':[0.,0.,-1.]}],
                    'skin_temperature_C':33.},sensory_blocks=[site['id']] if i else [])
                pending[i]=contact['sensory_inputs_hz'];outputs.append(output)
                self.assertAlmostEqual(controller.brain.time_s,(step+1)*.01)
                self.assertAlmostEqual(controller.brain.time_s,receptors.time_s)
            if step==0:
                self.assertEqual(outputs[0]['brain'],outputs[1]['brain'])
        self.assertNotEqual(outputs[0]['brain']['regional_state']['activity_hz'],
                            outputs[1]['brain']['regional_state']['activity_hz'])
        self.assertNotEqual(outputs[0]['motor_excitations']['addbrev_l'],
                            outputs[1]['motor_excitations']['addbrev_l'])
        self.assertEqual(outputs[1]['additional_sensory_inputs_hz'],{})

if __name__=='__main__':unittest.main()
