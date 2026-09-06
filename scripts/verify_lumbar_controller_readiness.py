"""Bounded retained 98-muscle controller/replay checks; no native launch."""
from copy import deepcopy
from pathlib import Path
import json,unittest,tempfile
from unittest.mock import patch
import numpy as np
from ihm.assembly.sensorimotor import SensorimotorController
ROOT=Path(__file__).resolve().parents[1]
REGISTRATION='data/derived/lumbar-muscle-native-lb45uirs/variant/registration.json'
BASE=json.loads((ROOT/'data/derived/mechanics/whole_body_arm26_v2/catalog.json').read_text())
ROWS=json.loads((ROOT/'data/derived/lumbar-muscle-native-lb45uirs/variant/catalog.json').read_text())
ADDED=ROWS[92:]

def observation(controller,stimulus=()):
    return {'time_s':controller.time_s,'foot_contact_force_n':{'r':0.,'l':0.},'muscles':{
        r['id']:{'fiber_length_m':r['optimal_fiber_length_m']*(1.1 if r['id'] in stimulus else 1.),
            'optimal_fiber_length_m':r['optimal_fiber_length_m'],'max_isometric_force_n':r['max_isometric_force_n'],
            'tendon_force_n':r['max_isometric_force_n']*(.2 if r['id'] in stimulus else 0.),
            'sensor_basis':'Synthetic CE fixture using retained native catalog normalization; no anatomical measurement'} for r in controller.catalog}}

def controller(rows=ROWS):return SensorimotorController.from_root(ROOT,muscle_catalog=rows)

class Tests(unittest.TestCase):
    def test_exact_catalog_and_inferred_regions(self):
        from ihm.assembly.embodied import _prepare_mechanical_registration
        frozen,manifest,rows=_prepare_mechanical_registration(ROOT,REGISTRATION)
        self.assertEqual(rows,ROWS);self.assertEqual(ROWS[:92],BASE);self.assertEqual(len(rows),98)
        self.assertIn(ROOT/REGISTRATION,frozen)
        self.assertEqual({r['id'] for r in ADDED},{'gait2392_'+n+'_'+s for n in ('ercspn','intobl','extobl') for s in ('l','r')})
        for row in ADDED:
            hemi='lh' if row['side']=='r' else 'rh'
            self.assertEqual(row['sensory_region'],f'brain-{hemi}-postcentral')
            self.assertEqual(row['motor_region'],f'brain-{hemi}-precentral')
            self.assertIn('engineering prior',row['assignment_basis'])

    def test_factory_override_receipt_and_pre_native_rejection(self):
        from ihm.assembly.embodied import EmbodiedRuntime,_prepare_mechanical_registration
        from scripts.verify_regional_embodied_factory import FactoryTests
        frozen,manifest,rows=_prepare_mechanical_registration(ROOT,REGISTRATION)
        fixture=FactoryTests()
        with tempfile.TemporaryDirectory() as directory:
            root,manifests=fixture.fixture(directory)
            for source,raw in frozen.items():
                target=root/source.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
            calls=[]
            with fixture.patches(root,manifests,calls):
                from ihm.assembly.articulated import ArticulatedBodyPlant
                ArticulatedBodyPlant.muscle_catalog=rows
                body=EmbodiedRuntime.from_workspace(root,root/'out',augmented_registration=REGISTRATION)
                self.assertEqual(body.plant.arguments['augmented_registration'],REGISTRATION)
                receipt=json.loads((root/'out/manifest.json').read_text())
                self.assertEqual(receipt['mechanical_registration_override']['model_sha256'],manifest['model_sha256'])
                for source,raw in frozen.items():self.assertEqual((root/'out/inputs'/source.relative_to(ROOT)).read_bytes(),raw)
                body.close()
            model=root/manifest['model_path'];raw=model.read_bytes();model.write_bytes(raw+b'changed')
            with patch('ihm.native.coupled_session.SignedCoupledNativeSession') as native:
                with self.assertRaisesRegex(ValueError,'hash differs'):
                    EmbodiedRuntime.from_workspace(root,root/'bad',augmented_registration=REGISTRATION)
                native.assert_not_called();self.assertFalse((root/'bad').exists())
            model.write_bytes(raw)
            outside=root/'outside';outside.write_bytes(raw);model.unlink();model.symlink_to(outside)
            with self.assertRaisesRegex(ValueError,'symlink'):_prepare_mechanical_registration(root,REGISTRATION)
            with self.assertRaises(ValueError):_prepare_mechanical_registration(root,'../outside')

    def test_each_added_sensory_descending_block_and_release(self):
        for row in ADDED:
            with self.subTest(muscle=row['id']):
                c=controller();key=row['id'];region=row['sensory_region']
                for _ in range(8):out=c.step(.01,observation(c,[key]),descending={key:.3})
                self.assertGreater(out['brain_sensory_inputs_hz'][region],0)
                self.assertAlmostEqual(out['sensors'][key]['length'],1.1)
                self.assertAlmostEqual(out['sensors'][key]['force'],.2)
                self.assertGreater(out['motor_excitations'][key],0)
                self.assertTrue(all(out['motor_excitations'][r['id']]==0 for r in ADDED if r['id']!=key))
                out=c.step(.01,observation(c,[key]),descending={key:.3},sensory_blocks=[key],motor_blocks=[key])
                self.assertEqual(out['motor_excitations'][key],0);self.assertNotIn(key,out['delayed_sensors'])
                self.assertEqual(out['brain_sensory_inputs_hz'].get(region,0),0)
                self.assertFalse(any(e[3]==key for e in c.events))
                out=c.step(.01,observation(c,[key]),descending={key:.3})
                self.assertEqual(out['motor_excitations'][key],0)
                for _ in range(8):out=c.step(.01,observation(c,[key]),descending={key:.3})
                self.assertGreater(out['motor_excitations'][key],0)
                self.assertAlmostEqual(c.time_s,c.brain.time_s)

    def test_existing_92_exact_behavior_when_added_ports_are_silent(self):
        old,new=controller(BASE),controller()
        for _ in range(8):
            a=old.step(.01,observation(old,['soleus_r']),descending={BASE[-1]['id']:.2})
            b=new.step(.01,observation(new,['soleus_r']),descending={BASE[-1]['id']:.2})
            self.assertEqual(a['motor_excitations'],{k:b['motor_excitations'][k] for k in a['motor_excitations']})
            np.testing.assert_array_equal(old.brain.state,new.brain.state)

    def test_order_independent_outputs_and_exact_replay_identity(self):
        a,b=controller(),controller(list(reversed(ROWS)))
        keys=[r['id'] for r in ADDED];command={k:.23 for k in keys}
        for _ in range(8):
            x=a.step(.01,observation(a,keys),descending=command);y=b.step(.01,observation(b,keys),descending=command)
            for key in a.muscles:self.assertAlmostEqual(x['motor_excitations'][key],y['motor_excitations'][key],places=13)
            np.testing.assert_allclose(a.brain.state,b.brain.state,rtol=0,atol=1e-12)
        saved=a.checkpoint();expected=a.step(.01,observation(a,keys),descending=command)
        a.restore(saved);self.assertEqual(expected,a.step(.01,observation(a,keys),descending=command))
        for other in (controller(BASE),b):
            with self.assertRaises(ValueError):other.restore(saved)
        changed=deepcopy(ROWS);changed[-1]['assignment_basis']+=' revised prior'
        with self.assertRaises(ValueError):controller(changed).restore(saved)

if __name__=='__main__':unittest.main()
