"""Embodied endpoint intake mass integration; fake owners, no native jobs."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import hashlib,json,tempfile,unittest
from unittest.mock import patch
import numpy as np
from ihm.assembly.embodied import EmbodiedRuntime
from ihm.assembly.intake_mass import IntakeMassBridge
from verify_embodied_runtime import Plant,Neural,Native,Load,Exchange
from verify_intake_mass import initial,consumed

class MassPlant(Plant):
    def __init__(self):
        super().__init__();self.closed=False;self.transfers=[];self.restore_calls=0;self.transfer_fail=False
        self.mass={'enabled':True,'reference_id':'mass-fixture','last_sequence':0,'owned_payload_mass_kg':0,'owners':{}}
    def snapshot(self):return {**super().snapshot(),'positive_muscle_work_j':.1,'mass_transfer':deepcopy(self.mass),'effective_native_body_mass_kg':70+self.mass['owned_payload_mass_kg']}
    def body_point(self,*,body,station_m):return {'kind':'body_point','body':body,'station_m':station_m,'time_s':self.t,'velocity_source_m_s':[.1,0,0]}
    def transfer_mass(self,**payload):
        self.transfers.append(payload);self.mass['last_sequence']=payload['sequence'];self.mass['owned_payload_mass_kg']+=payload['delta_mass_kg']
        self.mass['owners'][payload['owner']]={'body':payload['body'],'station_m':payload['station_m'],'mass_kg':self.mass['owned_payload_mass_kg']}
        self.mass['last_receipt']={'sequence':payload['sequence'],'owner':payload['owner'],'body':payload['body'],'delta_mass_kg':payload['delta_mass_kg']}
        if self.transfer_fail:raise TimeoutError('uncertain endpoint transfer')
        return self.snapshot()
    def restore(self,state):super().restore(state);self.restore_calls+=1
    def close(self):self.closed=True
class IntakeNative(Native):
    def __init__(self):super().__init__();self.intake=initial();self.bad_clock=False
    def snapshot(self):return {**super().snapshot(),'intake':deepcopy(self.intake)}
    def step(self,dt):
        self.t+=dt
        if self.intake['consumed_count']==0:self.intake=consumed()
        result=self.snapshot()
        if self.bad_clock:result['elapsed_s']+=.02
        return result
class RuntimeTests(unittest.TestCase):
    def make(self):
        plant=MassPlant();native=IntakeNative();bridge=IntakeMassBridge(plant,native.intake,body='torso',station_m=[.1,.2,.3],registration_identity='fixture-registration',incoming_velocity_basis='co_moving_at_ingestion_assumption')
        body=EmbodiedRuntime(plant,Neural(),native,Exchange(),Load(),intake_mass_bridge=bridge,intake_mass_binding={'registration_sha256':'fixture-registration'})
        return body,plant,native,bridge
    def test_endpoint_latched_mass_refresh_preserves_signed_interval(self):
        body,plant,native,bridge=self.make();self.assertEqual(body.snapshot()['intake_mass']['bridge']['applied_mass_kg'],0)
        body.schedule_intakes({'events':[{'event_id':'later','time_s':1.,'meal':{'water_ml':10}}]});self.assertEqual(plant.transfers,[])
        frame=body.step({});self.assertEqual(len(plant.transfers),1);self.assertEqual(native.demands,[(5.,4.,1.)])
        self.assertEqual(frame['mechanics']['effective_native_body_mass_kg'],70.013);self.assertEqual(frame['intake_mass']['bridge']['applied_mass_kg'],.013)
        self.assertEqual(frame['mechanics']['muscle_metabolic_energy_j'],2.1);self.assertEqual(body.time_s,.02)
        body.snapshot();body.snapshot();body.step({});self.assertEqual(len(plant.transfers),1);self.assertFalse(body.failed)
    def test_schedule_capacity_rejects_before_native_consumption_or_queue_mutation(self):
        body,plant,native,bridge=self.make()
        body.schedule_intakes({'events':[{'event_id':'water','time_s':0.,'meal':{'water_ml':500.}}]})
        before=body.snapshot();sequence=body.sequence
        with self.assertRaisesRegex(ValueError,'mechanical payload capacity'):
            body.schedule_intakes({'events':[{'event_id':'calcium','time_s':.02,'meal':{'calcium_mg':1.}}]})
        self.assertEqual(body.sequence,sequence);self.assertEqual(body.snapshot(),before)
        self.assertFalse(body.failed);self.assertEqual(plant.transfers,[]);self.assertEqual(native.intake['consumed_count'],0)

    def test_uncertain_post_consumption_transfer_aborts_both_without_rollback(self):
        body,plant,native,bridge=self.make();plant.transfer_fail=True
        with self.assertRaisesRegex(TimeoutError,'uncertain endpoint'):body.step({})
        self.assertTrue(body.failed);self.assertTrue(native.closed);self.assertTrue(plant.closed);self.assertTrue(bridge.failed)
        self.assertEqual(plant.restore_calls,0);self.assertEqual(plant.mass['owned_payload_mass_kg'],.013)
        self.assertEqual(body.time_s,0);self.assertIsNone(body.frame)
        with self.assertRaises(RuntimeError):body.step({})
    def test_bad_native_clock_aborts_before_mass_credit(self):
        body,plant,native,bridge=self.make();native.bad_clock=True
        with self.assertRaisesRegex(RuntimeError,'clock diverged'):body.step({})
        self.assertTrue(body.failed);self.assertTrue(native.closed);self.assertTrue(plant.closed);self.assertEqual(plant.transfers,[])
    def test_constructor_rejects_used_bridge_and_mismatched_provenance(self):
        body,plant,native,bridge=self.make()
        with self.assertRaisesRegex(ValueError,'provenance'):EmbodiedRuntime(plant,Neural(),native,Exchange(),Load(),intake_mass_bridge=bridge,intake_mass_binding={'registration_sha256':'wrong'})
        bridge.failed=True
        with self.assertRaisesRegex(ValueError,'unused intake bridge'):EmbodiedRuntime(plant,Neural(),native,Exchange(),Load(),intake_mass_bridge=bridge,intake_mass_binding={'registration_sha256':'fixture-registration'})
    def test_disabled_default_keeps_existing_contract(self):
        body=EmbodiedRuntime(Plant(),Neural(),Native(),Exchange(),Load());self.assertEqual(body.snapshot()['intake_mass'],{'enabled':False})
        self.assertEqual(body.step({})['intake_mass'],{'enabled':False})
    def test_flag_rejects_nonbool_before_paths_or_owners(self):
        for value in (None,0,1,'true',[],{}):
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'intake_mass.*bool'):EmbodiedRuntime.from_workspace('/nonexistent','/nonexistent/out',intake_mass=value)
    def test_binding_exact_neutral_station_and_source_identity(self):
        from ihm.assembly.embodied import bind_intake_mass
        with tempfile.TemporaryDirectory() as d:
            output=Path(d);raw=json.dumps({'entities':[{'id':'stomach-id','name':'stomach','centroid_m':[4.,6.,8.]}]}).encode();(output/'canonical_mechanics.json').write_bytes(raw);(output/'registration.json').write_bytes(b'fixture registration')
            plant=MassPlant();plant.output=output;plant.registration=SimpleNamespace(rows={'stomach-id':{'body':'torso','prior':'rigid'}},global_map=np.eye(4));transform=np.eye(4);transform[:3,3]=[1,2,3]
            plant.native=SimpleNamespace(snapshot=lambda:{'bodies':{'torso':{'transform_ground':transform.tolist()}}})
            bridge,binding=bind_intake_mass(plant,initial(),raw)
            self.assertEqual(bridge.station,[3.,4.,5.]);self.assertEqual(binding['canonical_sha256'],hashlib.sha256(raw).hexdigest())
            self.assertEqual(binding['registration_sha256'],hashlib.sha256(b'fixture registration').hexdigest())
            with self.assertRaisesRegex(ValueError,'canonical'):bind_intake_mass(plant,initial(),raw+b' ')

class FactoryTests(unittest.TestCase):
    def test_opt_in_selects_native_variant_retains_binding_and_regional_source(self):
        from verify_regional_embodied_factory import FactoryTests as Fixtures
        fixtures=Fixtures()
        for regional in (False,True):
            with self.subTest(regional=regional),tempfile.TemporaryDirectory() as d:
                root,manifests=fixtures.fixture(d);calls=[]
                canonical={'entities':[{'id':'stomach-id','name':'stomach','centroid_m':[.1,.2,.3]}]}
                (root/'data/derived/canonical/mechanics.json').write_text(json.dumps(canonical))
                class FactoryPlant(MassPlant):
                    muscle_catalog={}
                    def __init__(self,root,output,**options):
                        super().__init__();self.output=output;self.options=options;(output/'native').mkdir(parents=True)
                        (output/'native/execution.json').write_text('{}');(output/'canonical_mechanics.json').write_bytes((root/'data/derived/canonical/mechanics.json').read_bytes());(output/'registration.json').write_bytes(b'exact neutral registration')
                        self.registration=SimpleNamespace(rows={'stomach-id':{'body':'torso'}},global_map=np.eye(4))
                        self.native=SimpleNamespace(snapshot=lambda:{'bodies':{'torso':{'transform_ground':np.eye(4).tolist()}}},instance_mass_variant={'path':options['instance_mass_variant']})
                class FactoryNative(IntakeNative):
                    def __init__(self,config,output):
                        super().__init__();self.config=config;calls.append(config.engine_variant);output.mkdir(parents=True)
                        (output/'manifest.json').write_text(json.dumps({'library_sha256':manifests[config.engine_variant]['library_sha256'],'executable_sha256':'fixture','state_sha256':'fixture','patient_identity':{'Weight':{'unit':'kg','value':70}}}))
                with fixtures.patches(root,manifests,[]),patch('ihm.assembly.articulated.ArticulatedBodyPlant',FactoryPlant),patch('ihm.native.coupled_session.SignedCoupledNativeSession',FactoryNative),patch('ihm.native.regional_session.RegionalSignedNativeSession',FactoryNative):
                    body=EmbodiedRuntime.from_workspace(root,root/'out',intake_mass=True,regional_skin=regional)
                receipt=json.loads((root/'out/manifest.json').read_text());binding_raw=(root/'out/intake_mass_binding.json').read_bytes()
                self.assertTrue(receipt['intake_mass']);self.assertEqual(receipt['regional_skin'],regional)
                self.assertEqual(body.plant.options['instance_mass_variant'],'data/runtime/opensim/variants/instance_mass_v1')
                self.assertEqual(receipt['intake_mass_binding_sha256'],hashlib.sha256(binding_raw).hexdigest())
                for name in ('ihm/assembly/intake_mass.py','ihm/native/instance_mass.py'):
                    self.assertIn(name,receipt['sources']);self.assertIn(name,receipt['loaded_code'])
                self.assertEqual(receipt['intake_mass_binding']['station_source_m'],[.1,.2,.3])
                self.assertEqual(receipt['intake_mass_initial_bridge']['applied_mass_kg'],0)
                self.assertEqual(calls,['whole_body_integrity_regional_skin_graph_v2' if regional else 'whole_body_integrity_gi_absorption'])
                body.close()

def verify_native():
    """Coordinated combined regional/intake factory initialization, zero advances."""
    import time,resource
    root=Path(__file__).resolve().parents[1];output=Path(tempfile.mkdtemp(prefix='embodied-intake-mass-factory-',dir=root/'data/derived/audits'))
    started=time.monotonic();body=None
    try:
        body=EmbodiedRuntime.from_workspace(root,output/'body',regional_skin=True,intake_mass=True)
        frame=body.snapshot();audit=frame['intake_mass'];assert audit['enabled'] is True
        assert frame['time_s']==frame['mechanics']['time_s']==frame['physiology']['elapsed_s']==0
        assert audit['bridge']['applied_mass_kg']==0;assert audit['bridge']['mechanical_sequence']==0
        assert frame['mechanics']['mass_transfer']['enabled'] is True
        assert frame['physiology']['intake']['consumed_count']==0
        queued=body.schedule_intakes({'events':[{'event_id':'fixture_water','time_s':0,'meal':{'water_ml':10}}]})
        assert queued['intake_mass']['bridge']==audit['bridge'];assert queued['physiology']['intake']['pending'] is None
        manifest=json.loads((output/'body/manifest.json').read_text());assert manifest['intake_mass'] is True and manifest['regional_skin'] is True
        for source in ('ihm/assembly/intake_mass.py','ihm/native/instance_mass.py','ihm/native/regional_session.py'):
            assert source in manifest['sources'] and source in manifest['loaded_code']
        assert manifest['intake_mass_binding_sha256']==hashlib.sha256((output/'body/intake_mass_binding.json').read_bytes()).hexdigest()
        native_process=body.native.process;mechanical_process=body.plant.native.process;body.close();body=None
        assert native_process.poll() is not None and mechanical_process.poll() is not None
        report={'passed':True,'advances':0,'regional_skin':True,'intake_mass':True,'queued_events':1,'applied_mass_kg':0,'both_processes_reaped':True,
            'binding':manifest['intake_mass_binding'],'variant':manifest['intake_mass_variant'],'factory_manifest_sha256':hashlib.sha256((output/'body/manifest.json').read_bytes()).hexdigest(),
            'wall_s':time.monotonic()-started,'parent_max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'Actual combined factory initialization, initial/queued snapshots and cleanup only; no supported equilibrium or whole-body signed interval acceptance'}
        (output/'verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({'passed':True,'output':str(output),'wall_s':report['wall_s'],'advances':0},indent=2))
    except BaseException as error:
        (output/'failure.json').write_text(json.dumps({'error':str(error),'wall_s':time.monotonic()-started},indent=2)+'\n');raise
    finally:
        if body is not None:body.close()

if __name__=='__main__':
    import sys
    if '--run-native' in sys.argv:verify_native()
    else:unittest.main()
