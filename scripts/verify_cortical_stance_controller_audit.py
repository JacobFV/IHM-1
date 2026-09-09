"""Independent actual-weight stance adapter audit; native ports are fixtures."""
from copy import deepcopy
from pathlib import Path
import hashlib,json,tempfile,unittest
from unittest.mock import patch
import torch
from ihm.native.cortical_stance_controller import IBMCorticalStanceController
from ihm.native.cortical_stance import load_cortical_stance
from ihm.native.stance_lqr import NativeStanceLQR
from ihm.assembly.embodied import EmbodiedRuntime
from scripts.verify_embodied_runtime import Native,Load,Exchange
from scripts.verify_control_exchange import World

ROOT=Path(__file__).resolve().parents[1]


class NativePorts:
    def __init__(self,controller,folder):
        self.output=Path(folder);self.identity='independent-native-port-fixture'
        (self.output/'inputs').mkdir()
        model=(ROOT/controller.registration['model_path']).read_bytes()
        (self.output/'inputs/subject_walk_scaled.osim').write_bytes(model)
        (self.output/'inputs/augmentation_registration.json').write_bytes(controller.registration_bytes)
        self.execution={'source_sha256':{'subject_walk_scaled.osim':hashlib.sha256(model).hexdigest()},'target_mass_kg':controller.identity['target_mass_kg']}
        self.write_execution()
        catalog={row['id']:row for row in controller.catalog}
        muscles={name:{'activation':.01,'fiber_length_m':catalog[name]['optimal_fiber_length_m'],
            'optimal_fiber_length_m':catalog[name]['optimal_fiber_length_m'],'max_isometric_force_n':catalog[name]['max_isometric_force_n'],
            'tendon_force_n':0.,'sensor_basis':'independent actual-weight native-port fixture'} for name in controller.policy.muscle_names}
        coordinates={}
        for path,value in zip(controller.policy.state_names,controller.policy.x0.tolist()):
            parts=path.strip('/').split('/');name,variable=parts[-2:]
            if parts[0]=='jointset':coordinates.setdefault(name,{})[variable]=value
            else:muscles[name]['fiber_length_m' if variable=='fiber_length' else variable]=value
        self.state={'time_s':0.,'environment':'upright','mass_kg':controller.identity['target_mass_kg'],
            'muscles':muscles,'coordinates':coordinates}
    def write_execution(self):(self.output/'execution.json').write_text(json.dumps(self.execution))
    def snapshot(self):return deepcopy(self.state)
    def observation(self):
        return {'time_s':self.state['time_s'],'joints':deepcopy(self.state['coordinates']),
            'muscles':deepcopy(self.state['muscles']),'entities':{'chest':{'centroid_m':[0.,0.,0.]}},
            'foot_contact_force_n':{'r':0.,'l':0.}}


class PlantPorts:
    def __init__(self,native):self.native=native;self.m=0.;self.commands=[];self.interval_work=0.
    def snapshot(self):return {**self.native.observation(),'positive_muscle_work_j':self.interval_work,
        'muscle_metabolic_energy_j':self.m,'signed_active_fiber_work_j':0.,'muscle_heat_energy_j':self.m,
        'metabolic_reference':{'M0_w':0.,'H0_w':0.,'W0_w':0.}}
    def checkpoint(self):return deepcopy((self.native.state,self.m,self.commands,self.interval_work))
    def restore(self,saved):self.native.state,self.m,self.commands,self.interval_work=deepcopy(saved)
    def advance(self,dt,forces=(),actuation=None):
        assert set(actuation)==set(self.native.state['muscles'])
        self.commands.append(dict(actuation));self.m+=sum(actuation.values())*dt
        self.native.state['time_s']+=dt
        for name,value in actuation.items():self.native.state['muscles'][name]['activation']=value
        self.interval_work=0.
        return self.snapshot()
    def close(self):pass


class Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)
        cls.artifact=ROOT/'data/models/ibm_cortical_stance_v1/cortical_stance.pt'
        catalog=json.loads((ROOT/'data/models/engineering_stance_v1/catalog.json').read_text())
        cls.controller=IBMCorticalStanceController.from_root(ROOT,muscle_catalog=catalog)
        cls.folder=tempfile.TemporaryDirectory(prefix='independent-cortical-audit-')
        cls.native=NativePorts(cls.controller,cls.folder.name);cls.controller.bind_native(cls.native)
        cls.initial=cls.controller.checkpoint();cls.initial_native=cls.native.snapshot()
    @classmethod
    def tearDownClass(cls):cls.folder.cleanup()
    def setUp(self):
        self.native.state=deepcopy(self.initial_native);self.controller.restore(self.initial)
    def perturb(self):
        name=next(k for k in self.native.state['coordinates'] if 'pelvis_tilt' in k)
        self.native.state['coordinates'][name]['value']+=.003
    def test_exact_eight_state_replay_and_no_engineering_motor_fallback(self):
        c=self.controller;self.perturb()
        with (patch.object(c.inner,'step',side_effect=AssertionError('raw implicit motor fallback')),
             patch.object(c.inner.cord,'step',side_effect=AssertionError('cord fallback')),
             patch.object(NativeStanceLQR,'commands',side_effect=AssertionError('teacher fallback'))):
            first=c.step(.01,self.native.observation())
        endpoint=c.checkpoint();self.assertEqual(len(endpoint['cortical_state']),8)
        self.assertTrue(any(a!=b for a,b in zip(self.initial['cortical_state'],endpoint['cortical_state'])))
        c.restore(self.initial);self.assertEqual(c.step(.01,self.native.observation()),first)
        self.assertEqual(c.checkpoint(),endpoint);self.assertEqual(len(first['motor_excitations']),98)
    def test_runtime_ten_ms_control_and_all_98_energy_channels(self):
        c=self.controller;plant=PlantPorts(self.native);native=Native()
        body=EmbodiedRuntime(plant,c,native,Exchange(),Load(),environment_dynamics=World())
        frame=body.step({});self.assertEqual(c.steps,2);self.assertEqual(len(plant.commands),4)
        self.assertTrue(all(len(row)==98 for row in plant.commands))
        self.assertEqual(plant.commands[0],plant.commands[1]);self.assertEqual(plant.commands[2],plant.commands[3])
        expected=sum(sum(row.values())*.005 for row in plant.commands)
        self.assertAlmostEqual(frame['mechanics']['muscle_metabolic_energy_j'],expected)
        self.assertAlmostEqual(frame['coupling']['native_extra_metabolic_demand_w']*.02,expected)
        self.assertEqual(len(native.demands),1);self.assertEqual(frame['coupling']['control_interval_s'],.01)
        self.assertEqual(frame['coupling']['motor_exchange_latency_s'],0.)
    def test_severed_edges_and_sensory_flush_return_only_tonic_baseline(self):
        c=self.controller;self.perturb();self.assertFalse(c.sever)
        c.sever=True
        try:severed=c.step(.01,self.native.observation())
        finally:c.sever=False
        expected=dict(zip(c.policy.muscle_names,map(float,c.policy.u0)))
        self.assertEqual(severed['motor_excitations'],expected)
        c.restore(self.initial);c.step(.01,self.native.observation());self.native.state['time_s']=.01
        flushed=c.step(.01,self.native.observation(),sensory_blocks=[c.muscles[0]])
        self.assertEqual(flushed['motor_excitations'],expected)
    def test_all_eight_states_restore_atomically_after_late_failure(self):
        c=self.controller;self.perturb();advance=c.policy.advance
        def fail(vector,state,**kwargs):
            advance(vector,state,**kwargs)
            state[7].add_(1.)
            raise RuntimeError('late cortical failure')
        with patch.object(c.policy,'advance',side_effect=fail):
            with self.assertRaisesRegex(RuntimeError,'late cortical'):c.step(.01,self.native.observation())
        self.assertEqual(c.checkpoint(),self.initial)
        for index in range(8):
            bad=deepcopy(self.initial);bad['cortical_state'][index][0][0]=float('nan')
            with self.assertRaises(ValueError):c.restore(bad)
            self.assertEqual(c.checkpoint(),self.initial)
    def test_second_control_failure_rolls_back_cortex_world_and_native_ports(self):
        c=self.controller;plant=PlantPorts(self.native);native=Native();world=World()
        body=EmbodiedRuntime(plant,c,native,Exchange(),Load(),environment_dynamics=world)
        before=plant.checkpoint();advance=c.policy.advance;calls=[0]
        def fail(*args,**kwargs):
            result=advance(*args,**kwargs);calls[0]+=1
            if calls[0]==2:raise RuntimeError('second real cortical update failed')
            return result
        with patch.object(c.policy,'advance',side_effect=fail):
            with self.assertRaisesRegex(RuntimeError,'second real'):body.step({})
        self.assertEqual(c.checkpoint(),self.initial);self.assertEqual(plant.checkpoint(),before)
        self.assertEqual(world.t,0.);self.assertEqual(world.seen,[])
        self.assertEqual(native.demands,[]);self.assertEqual(native.loads,[]);self.assertFalse(body.failed)
    def test_validation_and_physiology_causality(self):
        c=self.controller;self.perturb()
        for dt,kwargs in ((.02,{}),(.01,{'additional_sensory_inputs_hz':{'wrong':1.}}),
            (.01,{'motor_blocks':['wrong']}),(.01,{'descending':{'muscle':1.}}),
            (.01,{'physiology':{'oxygen_saturation':float('nan')}}),(.01,{'descending':[]})):
            with self.assertRaises(ValueError):c.step(dt,self.native.observation(),**kwargs)
            self.assertEqual(c.checkpoint(),self.initial)
        baseline=c.step(.01,self.native.observation());c.restore(self.initial)
        hypoxic=c.step(.01,self.native.observation(),physiology={'oxygen_saturation':.2})
        self.assertNotEqual(baseline['motor_excitations'],hypoxic['motor_excitations']);c.restore(self.initial)
        warm=c.step(.01,self.native.observation(),physiology={'core_temperature_C':40.})
        self.assertNotEqual(baseline['motor_excitations'],warm['motor_excitations'])
    def test_artifact_and_native_binding_rejections(self):
        for kwargs in ({'dt_s':.02},{'model_sha256':'wrong'}):
            with self.assertRaises(ValueError):load_cortical_stance(self.artifact,**kwargs)
        saved=deepcopy(self.native.execution)
        for value in (76.,float('nan'),True):
            self.native.execution['target_mass_kg']=value;self.native.write_execution()
            with self.assertRaises(ValueError):self.controller.bind_native(self.native)
        self.native.execution=saved;self.native.write_execution()
        path=self.native.output/'inputs/augmentation_registration.json';original=path.read_bytes()
        path.write_bytes(original+b' ')
        try:
            with self.assertRaises(ValueError):self.controller.bind_native(self.native)
        finally:path.write_bytes(original)
    def test_adversarial_artifact_numeric_bindings_fail_closed(self):
        original=torch.load(self.artifact,map_location='cpu',weights_only=True)
        with tempfile.TemporaryDirectory(prefix='cortical-artifact-audit-') as directory:
            path=Path(directory)/'cortical_stance.pt'
            (path.parent/'pretrain_video_loop.py').write_bytes((self.artifact.parent/'pretrain_video_loop.py').read_bytes())
            for key,value in (('dt_s',float('nan')),('dt_s',True),('target_mass_kg',float('nan')),
                ('target_mass_kg',True),('reference_normalization','true')):
                bad=deepcopy(original);bad['provenance'][key]=value;torch.save(bad,path)
                with self.assertRaises(ValueError):load_cortical_stance(path,dt_s=.01)
            bad=deepcopy(original);bad['provenance']['reference_normalization']=False;torch.save(bad,path)
            with self.assertRaisesRegex(ValueError,'paired neural reference'):
                IBMCorticalStanceController.from_root(ROOT,muscle_catalog=self.controller.catalog,artifact_path=path)

    def test_each_bundle_runs_only_against_the_kernel_it_was_trained_from(self):
        from ihm.native.cortical_stance_controller import BUNDLES
        catalog=self.controller.catalog
        with self.assertRaisesRegex(ValueError,'Unknown cortical stance bundle'):
            IBMCorticalStanceController.from_root(ROOT,muscle_catalog=catalog,kind='not-a-bundle')
        for kind,bundle in BUNDLES.items():
            artifact=ROOT/bundle['artifact']
            if not artifact.is_file():continue
            controller=IBMCorticalStanceController.from_root(ROOT,muscle_catalog=catalog,kind=kind)
            metadata=controller.controller_metadata
            self.assertEqual(metadata['kind'],kind)
            self.assertEqual(metadata['bundle'],bundle['artifact'])
            self.assertEqual(metadata['kernel_bundle'],bundle['kernel'])
            # The kernel the cortex actually loaded is the one the artifact declares.
            self.assertEqual(controller.inner.identity['checkpoint_sha256'],
                             metadata['brain_checkpoint_sha256'])
            self.assertTrue(metadata['kernel_identity'])
            # Every other bundle's artifact must be refused under this kind, because
            # a bundle is only valid against the kernel it was trained from.
            for other,alternative in BUNDLES.items():
                path=ROOT/alternative['artifact']
                if other==kind or not path.is_file():continue
                with self.assertRaisesRegex(ValueError,'Shared IBM checkpoint differs'):
                    IBMCorticalStanceController.from_root(ROOT,muscle_catalog=catalog,
                                                          kind=kind,artifact_path=path)


if __name__=='__main__':unittest.main()
