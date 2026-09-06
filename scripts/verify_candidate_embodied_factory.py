"""Candidate body factory tests with fake native owners and fresh IBM processes."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]

class Tests(unittest.TestCase):
    def test_factory_exposes_source_pin(self):
        import inspect
        from ihm.assembly.embodied import EmbodiedRuntime
        self.assertIn('source_pin',inspect.signature(EmbodiedRuntime.from_workspace).parameters)

    def test_fresh_process_candidate_factory_receipts(self):
        result=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--probe'],cwd=ROOT,
            capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_fresh_process_hot_swap_fails_before_native(self):
        result=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--hot-swap'],cwd=ROOT,
            capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)


def probe(hot_swap=False):
    from dataclasses import replace
    from unittest.mock import patch
    from types import SimpleNamespace
    from scripts.verify_regional_embodied_factory import FactoryTests
    from ihm.brain.candidate import capture_candidate
    from ihm.brain.ibm_backend import IBMBackend
    from ihm.assembly.embodied import EmbodiedRuntime,bind_cutaneous
    fixture=FactoryTests()
    with tempfile.TemporaryDirectory() as directory:
        root,manifests=fixture.fixture(directory)
        pin=capture_candidate(ROOT.parent/'IBM-1',root/'candidate')
        calls=[]
        if hot_swap:
            IBMBackend()
            with fixture.patches(root,manifests,calls):
                try:EmbodiedRuntime.from_workspace(root,root/'out',source_pin=pin)
                except RuntimeError as error:assert 'fresh process' in str(error)
                else:raise AssertionError('Hot-swap accepted')
            assert calls==[] and not (root/'out').exists()
            return
        for bad in (replace(pin,package_sha256='0'*64),object()):
            with fixture.patches(root,manifests,calls):
                try:EmbodiedRuntime.from_workspace(root,root/'bad',source_pin=bad)
                except ValueError:pass
                else:raise AssertionError('Invalid source pin accepted')
            assert calls==[] and not (root/'bad').exists()
        # The candidate must lie inside the configured workspace for replay receipts.
        with fixture.patches(root,manifests,calls):
            try:EmbodiedRuntime.from_workspace(root/'other',root/'other/out',source_pin=pin)
            except ValueError as error:assert 'workspace' in str(error)
            else:raise AssertionError('Outside-workspace candidate accepted')
        assert calls==[]
        config={'regions':{'skin-contact-3':'brain-rh-postcentral'},
                'recruitment_hz_per_response':.1,'reference_temperature_C':33.}
        with fixture.patches(root,manifests,calls),patch('ihm.assembly.embodied.bind_cutaneous',return_value=SimpleNamespace(time_s=0.,audit={})) as bind:
            from ihm.assembly.sensorimotor import SensorimotorController
            body=EmbodiedRuntime.from_workspace(root,root/'out',source_pin=pin,
                surface_contact_manifest='explicit fake surface',cutaneous_configuration=config,regional_skin=True)
            assert SensorimotorController.from_root.call_args.kwargs['source_pin'] is pin
            assert bind.call_args.kwargs['source_pin'] is pin
        receipt=json.loads((root/'out/manifest.json').read_text())
        assert receipt['brain_source_pin']==pin.to_dict()
        assert receipt['regional_skin'] is True
        manifest=json.loads((pin.artifact_dir/'manifest.json').read_text())
        for relative in ('manifest.json','source_pin.json',*('source/'+p for p in manifest['files'])):
            source=pin.artifact_dir/relative
            archived=root/'out/inputs'/source.relative_to(root)
            assert archived.read_bytes()==source.read_bytes()
            assert str(source.relative_to(root)) in receipt['sources']
        assert receipt['brain_candidate_loaded_modules']['ibm.processes.neural']['package_sha256']==pin.package_sha256
        assert 'ihm/brain/candidate.py' in receipt['source_receipts']
        body.close()
        with patch('ihm.assembly.cutaneous_feedback.CutaneousFeedback',return_value=object()) as constructor:
            point={'id':'skin-contact-3','point_m':[0,0,0],'normal':[0,0,1],
                'contact_area_m2':.001,'material_identity':{},'indentation_basis':'fixture','area_basis':'fixture'}
            bind_cutaneous(root,[point],config,source_pin=pin)
            assert constructor.call_args.kwargs['source_pin'] is pin

def native_probe():
    """Actual candidate/regional/one-material-site init, snapshot and cleanup only."""
    import hashlib,resource,signal,time
    from ihm.brain.candidate import SourcePin
    from ihm.assembly.embodied import EmbodiedRuntime
    pin=SourcePin.load(ROOT/'data/derived/ibm-candidates/398375cc993ac17907d4ddd9f6d53eed2fab31de/source_pin.json')
    reference=ROOT/'data/derived/audits/cutaneous-factory-ylro2d66/body/manifest.json'
    retained=json.loads(reference.read_bytes())['cutaneous_materialization']
    selected=retained['sites'][0]
    configuration={'regions':{selected['id']:selected['sensory_region']},
        'recruitment_hz_per_response':retained['recruitment_hz_per_response'],
        'reference_temperature_C':selected['reference_temperature_C']}
    surface=ROOT/'data/derived/supine-surface-contact-exmzq9pq/manifest.json'
    assert hashlib.sha256(surface.read_bytes()).hexdigest()==selected['material_identity']['manifest_sha256']
    output=Path(tempfile.mkdtemp(prefix='candidate-regional-factory-',dir=ROOT/'data/derived/audits'))
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,4*1024**3))
    def deadline(signum,frame):raise TimeoutError('30 second candidate factory acceptance deadline')
    old=signal.signal(signal.SIGALRM,deadline);signal.alarm(30)
    started=time.monotonic();body=None
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    try:
        body=EmbodiedRuntime.from_workspace(ROOT,output/'body',source_pin=pin,
            regional_skin=True,surface_contact_manifest=str(surface.relative_to(ROOT)),
            cutaneous_configuration=configuration)
        initial=body.snapshot();receipt=json.loads((output/'body/manifest.json').read_bytes())
        assert initial['time_s']==initial['physiology']['elapsed_s']==initial['mechanics']['time_s']==0.
        assert body.neural.time_s==body.neural.brain.time_s==body.cutaneous.time_s==0.
        assert receipt['brain_source_pin']==pin.to_dict()
        assert body.neural.brain.source_identity['package_sha256']==body.cutaneous.audit['package_sha256']==pin.package_sha256
        contacts=initial['mechanics']['cutaneous_contacts']
        assert len(contacts)==1 and contacts[0]['material_identity']==selected['material_identity']
        assert contacts[0]['coordinate_frame']=='canonical_current_world'
        assert body.cutaneous.audit['sites'][0]['material_identity']==contacts[0]['material_identity']
        assert initial['tissue_exchange']['regional_skin']['available'] is True
        assert initial['physiology']['values']['tissue.regional_skin.completed_native_steps']==0
        candidate=json.loads((pin.artifact_dir/'manifest.json').read_bytes())
        for relative in ('manifest.json','source_pin.json',*('source/'+path for path in candidate['files'])):
            source=pin.artifact_dir/relative;archived=output/'body/inputs'/source.relative_to(ROOT)
            assert archived.read_bytes()==source.read_bytes()
            assert receipt['sources'][str(source.relative_to(ROOT))]==hashlib.sha256(source.read_bytes()).hexdigest()
        assert receipt['brain_candidate_loaded_modules']['ibm.processes.neural']['package_sha256']==pin.package_sha256
        write('initial.json',{'time_s':initial['time_s'],'physiology':initial['physiology'],
            'tissue_exchange':initial['tissue_exchange'],'cutaneous_contacts':contacts,
            'brain_source_identity':body.neural.brain.source_identity,'cutaneous_materialization':body.cutaneous.audit})
        physiology=body.native.process;mechanics=body.plant.native.process
        body.close();body=None
        assert physiology.poll() is not None and mechanics.poll() is not None
        report={'passed':True,'advances':0,'regional_skin':True,'selected_cutaneous_site':selected['id'],
            'prior_configuration':configuration,'prior_reference_manifest_sha256':hashlib.sha256(reference.read_bytes()).hexdigest(),
            'package_sha256':pin.package_sha256,'neural_source_sha256':pin.neural_source_sha256,
            'candidate_files_verified':len(candidate['files'])+2,'loaded_modules':len(receipt['brain_candidate_loaded_modules']),
            'factory_manifest_sha256':hashlib.sha256((output/'body/manifest.json').read_bytes()).hexdigest(),
            'native_pid':physiology.pid,'mechanical_pid':mechanics.pid,
            'native_exit_code':physiology.returncode,'mechanical_exit_code':mechanics.returncode,
            'both_processes_reaped':True,'wall_s':time.monotonic()-started,
            'parent_peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'maximum_child_peak_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            'parent_soft_address_space_bytes':2*1024**3,'hard_address_space_bytes':4*1024**3,
            'scope':'Actual candidate/regional/selected-skin factory init, snapshot, source retention and cleanup; no exchange step, receptor response, equilibrium or default promotion'}
        write('verification.json',report)
        print(json.dumps({'passed':True,'output':str(output),'wall_s':report['wall_s'],'both_processes_reaped':True}))
    except BaseException as error:
        write('failure.json',{'error':str(error),'wall_s':time.monotonic()-started})
        raise
    finally:
        if body:body.close()
        signal.alarm(0);signal.signal(signal.SIGALRM,old)

if __name__=='__main__':
    if sys.argv[1:]==['--native']:native_probe()
    elif sys.argv[1:]==['--probe']:probe()
    elif sys.argv[1:]==['--hot-swap']:probe(True)
    else:unittest.main()
