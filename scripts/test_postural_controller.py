"""Boundary tests for unvalidated engineering stance controller adapter."""
import sys,unittest
from pathlib import Path
from copy import deepcopy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.postural_controller import EngineeringStanceController
class Policy:
    def __init__(self,reference,**kwargs):self.muscles=reference['muscles'];self.last_allocation={};self.arms=kwargs['native_moment_arms']
    def update_native_moment_arms(self,arms):self.arms=deepcopy(arms)
    def commands(self,state):self.last_allocation={'test':True};return dict.fromkeys(self.muscles,.2)
class Native:
    identity='test-native'
    def __init__(self):
        self.calls=0;self.state={'time_s':0.,'coordinates':{},'muscles':{f'gait2392_{i}':{} for i in range(6)}}
    def snapshot(self):return deepcopy(self.state)
    def moment_arms(self,*,muscles,coordinates):
        self.calls+=1;return {'moment_arms_m':{m:{q:.1 for q in coordinates} for m in muscles}}
class Inner:
    def __init__(self,muscles):self.muscles=tuple(muscles);self.time_s=0.
    def checkpoint(self):return {'time_s':self.time_s}
    def restore(self,c):self.time_s=c['time_s']
    def step(self,dt,obs,**kwargs):
        self.time_s+=dt;return {'motor_excitations':dict.fromkeys(self.muscles,.5),'arc_max':{'stretch':.1}}
def controller():
    native=Native();p=EngineeringStanceController();p.inner=Inner(native.state['muscles'])
    p.model_sha256='test';p.controller_metadata={};p.native=None;p.reference=None;p.policy=None;p.cached_arms={};p.steps=0
    p.excitations=dict.fromkeys(p.inner.muscles,0.);p.policy_class=Policy;p.policy_config={};p.path_bytes=b'<xml/>'
    p.bind_native(native)
    return p,native,{'time_s':0.,'joints':{},'muscles':native.state['muscles']}
class Tests(unittest.TestCase):
    def test_cache_rollback_blocks(self):
        p,native,obs=controller();start=p.checkpoint()
        for _ in range(5):p.step(.01,obs)
        self.assertEqual(native.calls,1)
        saved=p.checkpoint();result=p.step(.01,obs,motor_blocks=['gait2392_0'])
        self.assertEqual(native.calls,2);self.assertEqual(result['motor_excitations']['gait2392_0'],0.)
        p.restore(saved);self.assertEqual(p.steps,5);self.assertEqual(p.inner.time_s,.05)
        result=p.step(.01,obs,sensory_blocks=['gait2392_1'])
        self.assertTrue(all(x==0 for x in result['motor_excitations'].values()))
        p.restore(start);self.assertIsNone(p.reference);self.assertEqual(p.steps,0)
    def test_reject_stale_and_cross_plant(self):
        p,native,obs=controller()
        with self.assertRaises(ValueError):p.step(.01,dict(obs,time_s=1.))
        with self.assertRaises(ValueError):p.bind_native(Native())
        saved=p.checkpoint();saved['native_identity']='another'
        with self.assertRaises(ValueError):p.restore(saved)
if __name__=='__main__':unittest.main()
