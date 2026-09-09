"""Current-boundary control scheduling and whole-exchange rollback fixtures."""
from copy import deepcopy
import unittest
import numpy as np
from scripts.verify_embodied_runtime import Plant,Native,Neural,Load,Exchange
from ihm.assembly.embodied import EmbodiedRuntime


class MovingPlant(Plant):
    def __init__(self):super().__init__();self.ports=[]
    def snapshot(self):
        return {**super().snapshot(),'positive_muscle_work_j':.01,
            'entities':{'chest':{'centroid_m':[self.t,0,0],'rotation_matrix':np.eye(3).tolist()}}}
    def advance(self,dt_s,forces=(),actuation=None):
        self.ports.append(deepcopy(forces));super().advance(dt_s,forces,actuation)
        return self.snapshot()


class Feedback(Neural):
    control_interval_s=.01
    current_interval_actuation=True
    def __init__(self):super().__init__();self.seen=[];self.fail_second=False
    def checkpoint(self):return deepcopy(self.__dict__)
    def restore(self,state):self.__dict__=deepcopy(state)
    def step(self,dt,observation,**inputs):
        assert dt==.01
        assert abs(observation['time_s']-self.t)<1e-10
        self.seen.append(observation['entities']['chest']['centroid_m'][0])
        if self.fail_second and len(self.seen)==2:raise ValueError('second feedback failed')
        command=.2+10*self.t;self.t+=dt
        return {'time_s':self.t,'motor_excitations':{'muscle':command}}


class World:
    def __init__(self):self.t=0.;self.seen=[]
    def checkpoint(self):return deepcopy(self.__dict__)
    def restore(self,state):self.__dict__=deepcopy(state)
    def advance(self,dt,entities,**kwargs):self.seen.append(entities['chest']['centroid_m'][0]);self.t+=dt;return []
    def frame(self):return {'time_s':self.t}


class Projection(Load):
    def project_load(self,forces,entities,volume):
        ports=[{**f,'generalized_force_pa':f['point_m'][0]**2} for f in forces]
        return {'external_pressure_pa':sum(p['generalized_force_pa'] for p in ports),
            'force_ports':ports,'ignored_nonrespiratory_ids':[]}


class Tests(unittest.TestCase):
    force={'id':'chest','point_m':[1.,0,0],'force_n':[2.,0,0]}
    def make(self,world=True):
        body=EmbodiedRuntime(MovingPlant(),Feedback(),Native(),Exchange(),Projection(),environment_dynamics=World() if world else None)
        body.next_excitation={'muscle':.99}
        return body
    def test_feedback_sees_fresh_pose_and_applies_at_each_boundary(self):
        body=self.make();frame=body.step({'forces':[self.force]})
        np.testing.assert_allclose(body.neural.seen,[0,.01])
        np.testing.assert_allclose([c['muscle'] for c in body.plant.commands],[.2,.2,.3,.3])
        np.testing.assert_allclose(body.environment_dynamics.seen,[0,.005,.01,.015])
        np.testing.assert_allclose([p[0]['point_m'][0] for p in body.plant.ports],[1,1.005,1.01,1.015])
        self.assertEqual(len(body.native.demands),1)
        self.assertAlmostEqual(body.native.t,.02)
        self.assertAlmostEqual(frame['coupling']['positive_muscle_work_j'],.04)
        self.assertEqual(frame['mechanics']['world_exchange']['substeps'],4)
        self.assertEqual(frame['coupling']['control_interval_s'],.01)
        self.assertEqual(frame['coupling']['motor_exchange_latency_s'],0.)
        self.assertEqual(body.next_excitation,{})
        load=frame['respiratory_load']
        self.assertAlmostEqual(load['external_pressure_pa'],np.mean(np.array([1,1.005,1.01,1.015])**2))
        self.assertAlmostEqual(sum(p['quadrature_weight'] for p in load['force_ports']),1.)
        np.testing.assert_allclose(np.sum([p['force_n'] for p in load['force_ports']],axis=0),[2,0,0])
    def test_no_world_keeps_two_current_control_intervals_and_blocks(self):
        body=self.make(False);frame=body.step({'motor_blocks':['muscle']})
        self.assertEqual(body.plant.commands,[{},{}])
        np.testing.assert_allclose(body.neural.seen,[0,.01])
        self.assertAlmostEqual(frame['coupling']['positive_muscle_work_j'],.02)
        self.assertEqual(len(body.native.demands),1)
    def test_second_feedback_failure_rolls_back_all_pre_native_owners(self):
        body=self.make();body.neural.fail_second=True
        before=(body.plant.checkpoint(),body.neural.checkpoint(),body.environment_dynamics.checkpoint())
        with self.assertRaisesRegex(ValueError,'second feedback'):body.step({'forces':[self.force]})
        self.assertEqual(body.plant.checkpoint(),before[0]);self.assertEqual(body.neural.checkpoint(),before[1])
        self.assertEqual(body.environment_dynamics.checkpoint(),before[2])
        self.assertEqual(body.native.demands,[]);self.assertEqual(body.native.loads,[])
        self.assertEqual(body.time_s,0);self.assertFalse(body.failed)
        self.assertEqual(body.next_excitation,{'muscle':.99})
        body.neural.fail_second=False;body.step({'forces':[self.force]})
        self.assertEqual(len(body.native.demands),1)
    def test_legacy_controller_still_uses_next_twenty_ms_interval(self):
        body=self.make(False);body.neural=Neural()
        first=body.step({});self.assertEqual(body.plant.commands,[{'muscle':.99}])
        self.assertEqual(first['coupling']['motor_exchange_latency_s'],.02)
        self.assertFalse(first['coupling']['current_interval_actuation'])
        body.step({});self.assertEqual(body.plant.commands[-1],{'muscle':.5})
    def test_nondividing_control_interval_rejects_before_owner_advances(self):
        body=self.make();body.neural.control_interval_s=.007
        with self.assertRaisesRegex(ValueError,'divide physiology'):body.step({})
        self.assertEqual(body.plant.commands,[]);self.assertEqual(body.neural.seen,[])
        self.assertEqual(body.environment_dynamics.t,0.);self.assertEqual(body.native.demands,[])
        self.assertFalse(body.failed)


if __name__=='__main__':unittest.main()
