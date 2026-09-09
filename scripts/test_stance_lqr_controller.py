"""Protocol/source-binding tests; synthetic fixture metadata is not native proof."""
import sys,json,tempfile,unittest,hashlib
from pathlib import Path
from copy import deepcopy
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.stance_lqr_controller import EngineeredLQRStanceController
ROOT=Path(__file__).resolve().parents[1]
class Inner:
    model_sha256='test-inner'
    def __init__(self,muscles):self.muscles=muscles;self.time_s=0.
    def checkpoint(self):return {'time_s':self.time_s}
    def restore(self,c):self.time_s=c['time_s']
    def step(self,dt,obs,**kw):
        self.time_s+=dt;return {'motor_excitations':dict.fromkeys(self.muscles,.5)}
class Native:
    identity='fixture-native'
    def __init__(self,output,state):self.output=output;self.state=state
    def snapshot(self):return deepcopy(self.state)
class Tests(unittest.TestCase):
    def test_source_sampling_blocks_restore(self):
        registration=ROOT/'data/derived/mechanics/resting_stance98/registration.json';reg=json.loads(registration.read_bytes())
        catalog=json.loads((ROOT/reg['catalog_path']).read_bytes())
        source=ROOT/'data/research/locomotion_control/linearization_eyxmbi_z/discrete_margin/linearization.npz'
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);artifact=out/'fixture.npz'
            with np.load(source) as d:values={k:d[k].copy() for k in d.files}
            values['target_mass_kg']=np.array(70.) # explicit synthetic test fixture
            np.savez(artifact,**values)
            with patch('ihm.assembly.sensorimotor.SensorimotorController.from_root',return_value=Inner([m['id'] for m in catalog])):
                p=EngineeredLQRStanceController.from_root(ROOT,muscle_catalog=catalog,artifact_path=artifact,registration_path=registration)
            original=artifact.read_bytes();artifact.write_bytes(b'changed')
            output=out/'native';(output/'inputs').mkdir(parents=True)
            (output/'inputs/subject_walk_scaled.osim').write_bytes((ROOT/reg['model_path']).read_bytes())
            (output/'inputs/augmentation_registration.json').write_bytes(registration.read_bytes())
            (output/'execution.json').write_text(json.dumps({'source_sha256':{'subject_walk_scaled.osim':reg['model_sha256']},'target_mass_kg':70.}))
            state=json.loads((ROOT/'data/derived/mechanics/resting_stance98/initial_snapshot.json').read_bytes())
            native=Native(output,state);p.bind_native(native)
            obs={'time_s':state['time_s'],'joints':state['coordinates'],'muscles':state['muscles']}
            checkpoint=p.checkpoint()
            with self.assertRaises(ValueError):p.step(.02,obs)
            result=p.step(.01,obs,sensory_blocks=[p.policy.muscle_names[0]])
            self.assertTrue(all(v==0 for v in result['motor_excitations'].values()))
            p.restore(checkpoint);self.assertEqual(p.steps,0)
            result=p.step(.01,obs,motor_blocks=[p.policy.muscle_names[0]])
            self.assertEqual(result['motor_excitations'][p.policy.muscle_names[0]],0.)
            p.restore(checkpoint);p.retain_sources(out)
            retained=out/'engineering-lqr-stance'
            self.assertEqual((retained/'linearization.npz').read_bytes(),original)
            for name,digest in p.identity['source_sha256'].items():self.assertEqual(hashlib.sha256((retained/name).read_bytes()).hexdigest(),digest)
            native.state['mass_kg']=77.
            with self.assertRaises(ValueError):p.step(.01,dict(obs,muscles=native.state['muscles']))
            self.assertEqual(p.steps,0)
if __name__=='__main__':unittest.main()
