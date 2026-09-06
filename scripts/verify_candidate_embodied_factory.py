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

if __name__=='__main__':
    if sys.argv[1:]==['--probe']:probe()
    elif sys.argv[1:]==['--hot-swap']:probe(True)
    else:unittest.main()
