"""Small causal exchange/failed-native checks; no full model or subprocess."""
from copy import deepcopy
import unittest
from ihm.assembly.embodied import EmbodiedRuntime

class Plant:
    def __init__(self):self.t=0.;self.force=0.;self.commands=[];self.rate=105.
    def snapshot(self):return {'time_s':self.t,'entities':{},'muscles':{},'foot_contact_force_n':{'r':0,'l':0},'total_muscle_metabolic_w':100.,'muscle_metabolic_energy_j':self.rate*self.t,'signed_active_fiber_work_j':2*self.t,'muscle_heat_energy_j':(self.rate-2)*self.t,'metabolic_reference':{'M0_w':100.,'W0_w':1.,'H0_w':99.}}
    def checkpoint(self):return deepcopy(self.__dict__)
    def restore(self,state):self.__dict__=deepcopy(state)
    def advance(self,dt_s,forces=(),actuation=None):
        self.commands.append(dict(actuation or {}));self.t+=dt_s
        return {**self.snapshot(),'positive_muscle_work_j':.1}
    def close(self):pass

class Neural:
    def __init__(self):self.t=0;self.inputs=[]
    def checkpoint(self):return self.t
    def restore(self,state):self.t=state
    def step(self,dt,observation,**inputs):
        assert observation['time_s']==self.t;self.t+=dt;self.inputs.append(deepcopy(inputs))
        return {'time_s':self.t,'motor_excitations':{'muscle':.5}}

class Native:
    def __init__(self):self.t=0;self.loads=[];self.demands=[];self.closed=False;self.fail=False;self.meals=[];self.meal_fail=False
    def snapshot(self):return {'elapsed_s':self.t,'time_s':100+self.t,'values':{
        'mean_arterial_pressure_mmhg':90,'oxygen_saturation':.98,'core_temperature_c':37,
        'maximum_work_rate_w':100,'lung_volume_ml':3000,'coupling.muscle_unmet_kcal':0}}
    def respiratory_load(self,p):self.loads.append(p)
    def exercise(self,f):raise AssertionError('Generic exercise must not own signed muscle demand')
    def signed_step(self,reference,m,h,w):
        assert abs(m-h-w)<1e-9;self.demands.append((m,h,w));return self.step(.02)
    def meal(self,meal):
        self.meals.append(meal)
        if self.meal_fail:raise RuntimeError('lost meal acknowledgment')
        return {'status':'ok','sequence':len(self.meals),'elapsed_s':self.t,'pending_meal':True}
    def step(self,dt):
        self.t+=dt
        if self.fail:raise RuntimeError('Native interrupted after advance')
        return self.snapshot()
    def close(self,graceful=True):self.closed=True

class Load:
    bindings={'chest':None}
    def project_load(self,forces,entities,volume):return {'external_pressure_pa':sum(f['force_n'][0] for f in forces)}

    def geometry(self,volume,entities,time):return {'entities':deepcopy(entities),'skin_field':{},'time_s':time}

class Exchange:
    def observe(self,snapshot):return {'time_s':snapshot['time_s'],'owner':'native'}

class Tests(unittest.TestCase):
    def body(self):return EmbodiedRuntime(Plant(),Neural(),Native(),Exchange(),Load())
    def test_cutaneous_receptor_endpoint_feeds_next_brain_exchange_and_blocks(self):
        from pathlib import Path
        from ihm.assembly.cutaneous_feedback import CutaneousFeedback
        site={'id':'fixture','position_m':[0,0,0],'normal':[0,0,1],
              'contact_area_m2':.001,'stiffness_pa_per_m':1e7,
              'sensory_region':'brain-rh-postcentral','reference_temperature_C':33,
              'support_basis':'synthetic runtime timing fixture'}
        receptor=CutaneousFeedback(Path(__file__).resolve().parents[1],sites=[site],recruitment_hz_per_response=.1)
        body=EmbodiedRuntime(Plant(),Neural(),Native(),Exchange(),Load(),cutaneous=receptor)
        body.mechanical_state['cutaneous_contacts']=[{'id':'fixture','force_n':[0,0,-1]}]
        first=body.step({})
        self.assertEqual(body.neural.inputs[0]['additional_sensory_inputs_hz'],{})
        self.assertGreater(first['cutaneous']['sensory_inputs_hz']['brain-rh-postcentral'],0)
        body.mechanical_state['cutaneous_contacts']=[]
        body.step({})
        self.assertGreater(body.neural.inputs[1]['additional_sensory_inputs_hz']['brain-rh-postcentral'],0)
        body.mechanical_state['cutaneous_contacts']=[]
        body.step({'skin_sensory_blocks':['fixture']})
        self.assertEqual(body.neural.inputs[2]['additional_sensory_inputs_hz'],{})
        self.assertAlmostEqual(receptor.time_s,body.time_s)
        before=receptor.checkpoint()
        body.mechanical_state['cutaneous_contacts']=[{'id':'fixture','force_n':[0,0,-2]}]
        project=body.respiratory_load.project_load
        def fail(*args):raise ValueError('load projection rejected')
        body.respiratory_load.project_load=fail
        with self.assertRaisesRegex(ValueError,'load projection'):body.step({})
        self.assertEqual(before,receptor.checkpoint())
        self.assertAlmostEqual(body.plant.t,body.time_s)
        self.assertAlmostEqual(body.neural.t,body.time_s)
        body.respiratory_load.project_load=project
        body.mechanical_state.pop('cutaneous_contacts')
        # Missing physical observation is unknown, never silently a release.
        with self.assertRaisesRegex(ValueError,'cutaneous contacts'):body.step({})
        self.assertEqual(before,receptor.checkpoint())
        self.assertFalse(body.failed)

    def test_missing_native_sensor_is_unknown_not_release(self):
        from pathlib import Path
        from ihm.assembly.embodied import bind_cutaneous
        contact={'id':'skin-contact-3','point_m':[0,0,0],'normal':[0,0,-1],
            'force_n':[0,0,1],'contact_area_m2':.001,'indentation_m':.0001,
            'material_identity':{'manifest_sha256':'a'*64,'quadrature_index':3,'triangle_index':8},
            'indentation_basis':'native modeled compression','area_basis':'reference quadrature'}
        receptor=bind_cutaneous(Path(__file__).resolve().parents[1],[contact],
            {'regions':{'skin-contact-3':'brain-rh-postcentral'},'recruitment_hz_per_response':.1,'reference_temperature_C':33})
        body=EmbodiedRuntime(Plant(),Neural(),Native(),Exchange(),Load(),cutaneous=receptor)
        body.mechanical_state['cutaneous_contacts']=[]
        with self.assertRaisesRegex(ValueError,'native skin sensor'):body.step({})
        self.assertEqual(body.native.loads,[])
        self.assertEqual(body.time_s,0)
        self.assertEqual(receptor.time_s,0)
        body.mechanical_state['cutaneous_contacts']=[{**contact,'force_n':[0,0,0],'indentation_m':0}]
        frame=body.step({})
        self.assertEqual(frame['cutaneous']['sites'][0]['indentation_um'],0)

    def test_delayed_actuation_native_load_and_work(self):
        body=self.body();first=body.step({'forces':[{'id':'chest','force_n':[2,0,0],'point_m':[0,0,0]}]})
        self.assertEqual(body.plant.commands,[{}]);self.assertEqual(body.native.loads,[2])
        self.assertAlmostEqual(body.native.demands[0][0],5.);self.assertEqual(first['time_s'],.02)
        body.step({});self.assertEqual(body.plant.commands[-1],{'muscle':.5})
        self.assertEqual(body.native.loads[-1],0)
    def test_scheduled_intake_changes_sequence_without_advancing_then_delivers_once(self):
        body=self.body()
        frame=body.schedule_intakes({'events':[{'event_id':'water_1','time_s':.02,'meal':{'water_ml':50}}]})
        self.assertEqual(frame['sequence'],1);self.assertEqual(frame['time_s'],0)
        self.assertEqual(frame['intake_schedule']['events'][0]['state'],'queued')
        body.step({});self.assertEqual(body.native.meals,[])
        frame=body.step({});self.assertEqual(len(body.native.meals),1)
        self.assertEqual(frame['intake_schedule']['events'][0]['state'],'accepted')
        body.step({});self.assertEqual(len(body.native.meals),1)
        with self.assertRaises(ValueError):body.schedule_intakes({'events':[{'event_id':'water_1','time_s':.06,'meal':{'water_ml':50}}]})
    def test_lost_meal_ack_is_retained_uncertain_and_aborts(self):
        body=self.body();body.native.meal_fail=True
        body.schedule_intakes({'events':[{'event_id':'meal','time_s':0,'meal':{'carbohydrate_g':10}}]})
        with self.assertRaisesRegex(RuntimeError,'lost meal'):body.step({})
        self.assertTrue(body.failed);self.assertEqual(len(body.native.meals),1)
        self.assertEqual(body.intakes.snapshot()['events'][0]['state'],'uncertain')
    def test_signed_negative_increment_is_preserved_without_generic_exercise(self):
        body=self.body();body.plant.rate=95
        frame=body.step({})
        self.assertEqual(body.native.demands,[(-5.,-6.,1.)])
        self.assertEqual(frame['coupling']['native_extra_metabolic_demand_w'],-5)
        self.assertFalse(body.failed)
    def test_horizon_preflight_has_no_side_effects(self):
        from types import SimpleNamespace
        body=self.body();body.native.config=SimpleNamespace(horizon_s=0)
        with self.assertRaisesRegex(ValueError,'horizon'):body.step({})
        self.assertEqual(body.plant.t,0);self.assertEqual(body.native.loads,[])
    def test_close_failure_retains_owner_for_retry(self):
        body=self.body();original=body.plant.close
        def fail():raise RuntimeError('cleanup failed')
        body.plant.close=fail
        with self.assertRaisesRegex(RuntimeError,'unconfirmed'):body.close()
        self.assertFalse(body.closed);self.assertTrue(body.native.closed)
        body.plant.close=original;body.close();self.assertTrue(body.closed)
    def test_uncertain_native_commit_aborts_instead_of_fake_rollback(self):
        body=self.body();body.native.fail=True
        with self.assertRaises(RuntimeError):body.step({})
        self.assertTrue(body.failed);self.assertTrue(body.native.closed)
        with self.assertRaises(RuntimeError):body.step({})
    def test_invalid_input_does_not_advance_any_owner(self):
        body=self.body()
        for data in [{'seconds':.03},{'bad':2},{'forces':[{'id':'chest','force_n':[float('nan'),0,0],'point_m':[0,0,0]}]}]:
            with self.assertRaises(ValueError):body.step(data)
        self.assertEqual(body.plant.t,0);self.assertEqual(body.neural.t,0);self.assertEqual(body.native.t,0)

if __name__=='__main__':unittest.main()
