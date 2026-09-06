"""Focused factory routing/provenance tests; all native/mechanical owners are fakes."""
from pathlib import Path
from contextlib import ExitStack
from hashlib import sha256
import json,tempfile,unittest
from unittest.mock import patch
from ihm.assembly.embodied import EmbodiedRuntime

class FactoryTests(unittest.TestCase):
    def fixture(self,base):
        root=Path(base);variants=root/'data/runtime/physiology/variants'
        manifests={};parent=None
        for name in ('whole_body_integrity_evaporation_humidity','whole_body_integrity_gi_absorption','whole_body_integrity_regional_skin_graph_v2'):
            directory=variants/name;directory.mkdir(parents=True);library=name.encode();(directory/'libbiogears.so.8.0.0').write_bytes(library)
            entry={'library_sha256':sha256(library).hexdigest()}
            if parent:
                raw=(variants/parent/'manifest.json').read_bytes()
                entry.update(parent_variant=parent,parent_manifest_sha256=sha256(raw).hexdigest(),parent_library_sha256=manifests[parent]['library_sha256'])
            (directory/'manifest.json').write_text(json.dumps(entry));manifests[name]=entry;parent=name
        state=root/'initial.xml';state.write_bytes(b'explicit fake native state')
        ref=root/'data/derived/systemic/exertion_v3/exercise/native/manifest.json';ref.parent.mkdir(parents=True)
        ref.write_text(json.dumps({'configuration':{'engine_variant':'whole_body_integrity_evaporation_humidity','state_path':str(state)},
            'state_sha256':sha256(state.read_bytes()).hexdigest(),'library_sha256':manifests['whole_body_integrity_evaporation_humidity']['library_sha256']}))
        for name in ('respiration','brain','anatomy','mechanics','microvascular','profile'):
            path=root/f'data/derived/canonical/{name}.json';path.parent.mkdir(parents=True,exist_ok=True);path.write_text('{}')
        return root,manifests

    def patches(self,root,manifests,calls):
        class Native:
            def __init__(self,config,output,kind):
                calls.append((kind,config));self.config=config;self.kind=kind;self.closed=False
                output.mkdir(parents=True)
                (output/'manifest.json').write_text(json.dumps({'library_sha256':manifests[config.engine_variant]['library_sha256'],
                    'executable_sha256':'fake-executable','state_sha256':'fake-state','patient_identity':{'Weight':{'unit':'kg','value':70}}}))
            def snapshot(self):return {'elapsed_s':0,'time_s':0,'values':{'lung_volume_ml':3000}}
            def close(self,graceful=True):self.closed=True
        class Plant:
            muscle_catalog={}
            def __init__(self,root,output,**kwargs):
                self.arguments=kwargs;(output/'native').mkdir(parents=True);(output/'native/execution.json').write_text('{}')
            def snapshot(self):return {'time_s':0,'entities':{},'metabolic_reference':{'M0_w':1,'H0_w':1,'W0_w':0}}
            def close(self):pass
        # Factory routing fixtures intentionally replace source attestation with named
        # file receipts. Actual loaded-code attestation remains covered by native runs.
        def source_receipt(module):
            path=root/(module.__name__.replace('.','/')+'.py');path.parent.mkdir(parents=True,exist_ok=True)
            raw=('# fake receipt for '+module.__name__).encode();path.write_bytes(raw)
            return {'path':path,'bytes':raw,'loaded_code_sha256':sha256(raw).hexdigest()}
        stack=ExitStack()
        stack.enter_context(patch('ihm.native.coupled_session.SignedCoupledNativeSession',side_effect=lambda c,o:Native(c,o,'default')))
        stack.enter_context(patch('ihm.native.regional_session.RegionalSignedNativeSession',side_effect=lambda c,o:Native(c,o,'regional')))
        stack.enter_context(patch('ihm.assembly.articulated.ArticulatedBodyPlant',Plant))
        stack.enter_context(patch('ihm.assembly.sensorimotor.SensorimotorController.from_root',return_value=object()))
        stack.enter_context(patch('ihm.assembly.body_exchange.NativeTissueExchange.from_workspace',return_value=object()))
        stack.enter_context(patch('ihm.assembly.embodied_respiration.EmbodiedRespiration',return_value=object()))
        stack.enter_context(patch('ihm.assembly.interactive_scene._loaded_source',side_effect=source_receipt))
        return stack

    def test_default_and_opt_in_select_distinct_native_owners_and_capture_source(self):
        for enabled,kind,variant in [(False,'default','whole_body_integrity_gi_absorption'),(True,'regional','whole_body_integrity_regional_skin_graph_v2')]:
            with self.subTest(regional_skin=enabled),tempfile.TemporaryDirectory() as temporary:
                root,manifests=self.fixture(temporary);calls=[]
                with self.patches(root,manifests,calls):
                    options={'regional_skin':True} if enabled else {}
                    body=EmbodiedRuntime.from_workspace(root,root/'output',source_pin=None,**options)
                self.assertEqual([(k,c.engine_variant) for k,c in calls],[(kind,variant)])
                self.assertIsNone(body.plant.arguments['surface_contact_manifest'])
                receipt=json.loads((root/'output/manifest.json').read_text())
                self.assertIs(receipt['regional_skin'],enabled)
                if enabled:
                    self.assertIn('ihm/native/regional_session.py',receipt['sources'])
                    self.assertIn('ihm/native/regional_session.py',receipt['loaded_code'])
                self.assertIn('ihm/assembly/regional_exchange.py',receipt['sources'])
                body.close()

    def test_regional_flag_requires_bool_before_files_or_native_owners(self):
        for value in (None,0,1,'true',[],{}):
            with self.subTest(value=value),tempfile.TemporaryDirectory() as temporary:
                root=Path(temporary)
                with self.assertRaisesRegex(ValueError,'regional_skin.*bool'):
                    EmbodiedRuntime.from_workspace(root,root/'output',regional_skin=value)
                self.assertFalse((root/'output').exists())

    def test_opt_in_preserves_lineage_failure_before_native_creation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root,manifests=self.fixture(temporary);calls=[]
            path=root/'data/runtime/physiology/variants/whole_body_integrity_regional_skin_graph_v2/manifest.json'
            bad=json.loads(path.read_text());bad['parent_manifest_sha256']='changed';path.write_text(json.dumps(bad))
            with self.patches(root,manifests,calls),self.assertRaisesRegex(ValueError,'lineage manifest changed'):
                EmbodiedRuntime.from_workspace(root,root/'output',source_pin=None,regional_skin=True)
            self.assertEqual(calls,[])
            self.assertFalse((root/'output').exists())

def verify_native():
    """Actual opt-in factory init/snapshot/cleanup; deliberately zero advances."""
    import time,resource
    root=Path(__file__).resolve().parents[1]
    parent=Path(tempfile.mkdtemp(prefix='regional-embodied-factory-',dir=root/'data/derived/audits'))
    started=time.monotonic();body=None
    try:
        body=EmbodiedRuntime.from_workspace(root,parent/'body',regional_skin=True)
        initial=body.snapshot()
        assert initial['time_s']==initial['physiology']['elapsed_s']==initial['mechanics']['time_s']==0
        regional=initial['tissue_exchange']['regional_skin'];assert regional['available'] is True
        values=initial['physiology']['values']
        assert values['tissue.regional_skin.completed_native_steps']==0
        for region,fraction in [('region_a',.2),('region_b',.3),('residual',.5)]:
            assert values['tissue.regional_skin.'+region+'.fraction']==fraction
            assert values['tissue.regional_skin.'+region+'.requested_pressure_pa']==0
            assert values['tissue.regional_skin.'+region+'.applied_pressure_pa']==0
        manifest=json.loads((parent/'body/manifest.json').read_text())
        assert manifest['regional_skin'] is True
        assert 'ihm/native/regional_session.py' in manifest['sources']
        assert 'ihm/native/regional_session.py' in manifest['loaded_code']
        assert 'ihm/assembly/regional_exchange.py' in manifest['sources']
        native=body.native;mechanical=body.plant.native
        native_pid=native.process.pid;mechanical_pid=mechanical.process.pid
        (parent/'initial-native-snapshot.json').write_text(json.dumps({'time_s':initial['time_s'],'physiology':initial['physiology'],
            'tissue_exchange':initial['tissue_exchange'],'mechanical_time_s':initial['mechanics']['time_s']},indent=2)+'\n')
        body.close();assert native.process.poll() is not None and mechanical.process.poll() is not None
        body=None
        report={'passed':True,'regional_skin':True,'advances':0,'native_pid':native_pid,'mechanical_pid':mechanical_pid,
            'native_exit_code':native.process.returncode,'mechanical_exit_code':mechanical.process.returncode,
            'both_processes_reaped':True,'source_receipts_verified':True,'factory_manifest_sha256':sha256((parent/'body/manifest.json').read_bytes()).hexdigest(),
            'wall_s':time.monotonic()-started,'parent_peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'Opt-in native regional factory initialization, initial snapshot and cleanup only; no equilibrium or coupled step acceptance'}
        (parent/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'passed':True,'output':str(parent),'wall_s':report['wall_s']},indent=2))
    except BaseException as error:
        (parent/'failure.json').write_text(json.dumps({'error':str(error),'wall_s':time.monotonic()-started},indent=2)+'\n')
        raise
    finally:
        if body:body.close()

if __name__=='__main__':
    import sys
    if sys.argv[1:]==['--native']:verify_native()
    else:unittest.main()
