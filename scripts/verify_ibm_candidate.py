"""Source-only immutable capture and opt-in compatibility checks; no graph builds."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
DONOR=ROOT.parent/'IBM-1'

class CandidateTests(unittest.TestCase):
    def test_candidate_helper_exists(self):
        self.assertTrue((ROOT/'ihm/brain/candidate.py').is_file(),'Immutable candidate helper missing')

    def test_capture_is_immutable_and_exact_commit_not_dirty_worktree(self):
        from ihm.brain.candidate import capture_candidate, verify_pin
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);donor=root/'donor';donor.mkdir()
            def git(*args):return subprocess.run(['git','-C',str(donor),*args],check=True,capture_output=True,text=True).stdout.strip()
            git('init');git('config','user.name','Fixture');git('config','user.email','fixture@example.invalid')
            (donor/'ibm/processes').mkdir(parents=True);(donor/'ibm/__init__.py').write_text('')
            source=donor/'ibm/processes/neural.py'
            source.write_text('def _sigmoid(x): return x\ndef wilson_cowan_excitatory(x, theta): return x\ndef shunting_inhibition_rate(x, theta): return x\n')
            git('add','.');git('commit','-m','fixture');committed=source.read_bytes()
            source.write_text('uncommitted training development must not be captured\n')
            dest=root/'candidate';pin=capture_candidate(donor,dest,baseline_neural=ROOT/'data/derived/canonical/brain-sources/ibm-neural.py')
            identity, snapshots=verify_pin(pin)
            self.assertEqual(snapshots['ibm/processes/neural.py'],committed)
            self.assertEqual(identity['donor_commit'],git('rev-parse','HEAD'))
            self.assertTrue(identity['donor_status'])
            from dataclasses import replace
            with self.assertRaises(ValueError):verify_pin(replace(pin,package_sha256='0'*64))
            self.assertEqual(identity['evidence_policy'],'compiled_source_priors_only')
            self.assertFalse(identity['training_or_datasets_included'])
            extra=dest/'source/ibm/extra.py';extra.write_text('unexpected=1')
            with self.assertRaises(ValueError):verify_pin(pin)
            extra.unlink()
            with self.assertRaises(ValueError):capture_candidate(donor,dest)
            (dest/'source/ibm/processes/neural.py').write_text('changed')
            with self.assertRaises(ValueError):verify_pin(pin)

    def test_shared_pin_in_fresh_process(self):
        from ihm.brain.candidate import capture_candidate
        with tempfile.TemporaryDirectory() as directory:
            pin=capture_candidate(DONOR,Path(directory)/'candidate',baseline_neural=ROOT/'data/derived/canonical/brain-sources/ibm-neural.py')
            command=[sys.executable,str(Path(__file__).resolve()),'--probe',str(pin.artifact_dir/'source_pin.json')]
            result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            receipt=json.loads(result.stdout.strip().splitlines()[-1])
            self.assertEqual(receipt['package_sha256'],pin.package_sha256)
            self.assertTrue(receipt['shared_identity'])


def probe(path):
    import numpy as np
    from ihm.brain.candidate import SourcePin
    from ihm.brain.ibm_backend import IBMBackend
    from ihm.assembly.brain import BodyBrain
    from ihm.assembly.cutaneous_feedback import CutaneousFeedback
    from ihm.brain.causal import CausalIBM
    pin=SourcePin.load(path);backend=IBMBackend(source_pin=pin)
    from ibm.topologies import microcircuit_prior, vascular_prior
    assert not any(item.from_evidence for item in (microcircuit_prior.MEASURED,
        microcircuit_prior.HUMAN,microcircuit_prior.KINETICS,vascular_prior.MEASURED))
    from ihm.assembly.sensorimotor import SensorimotorController
    controller=SensorimotorController.from_root(ROOT,source_pin=pin)
    legacy=SensorimotorController.from_root(ROOT)
    assert controller.model_sha256!=legacy.model_sha256
    before=controller.checkpoint()
    try:controller.restore(legacy.checkpoint())
    except ValueError:pass
    else:raise AssertionError('Cross-law checkpoint replay accepted')
    assert controller.checkpoint()==before
    assert controller.brain.source_identity['package_sha256']==pin.package_sha256
    brain=BodyBrain(json.loads((ROOT/'data/derived/canonical/brain.json').read_text()),root=ROOT,source_pin=pin)
    site={'id':'fixture','position_m':[0.,0.,0.],'normal':[0.,0.,1.],
          'contact_area_m2':.001,'stiffness_pa_per_m':1e7,'sensory_region':'brain-rh-postcentral',
          'reference_temperature_C':33.,'support_basis':'synthetic source parity fixture'}
    skin=CutaneousFeedback(ROOT,sites=[site],recruitment_hz_per_response=.1,source_pin=pin)
    for kind in ('rapid','slow'):
        model=CausalIBM(backend,sites_m=[[0,0,0]],kind=kind)
        assert model.audit['donor_relative_error']<1e-9
        assert model.advance([1.],.01)[0]>0
    pulse=np.zeros((1,64));pulse[:,8:16]=1.
    for kind in ('rapid','slow'):
        window=backend.run_window(pulse,kind=kind,sites_m=[[0,0,0]])
        basis=backend.basis(64,.002)
        expected=basis.synthesize(basis.analyze(pulse)*backend.transfer(kind,basis)/(1+1j*basis.omega))
        np.testing.assert_allclose(window.values,expected,atol=1e-10,rtol=1e-8)
    sample={'neural.exc.potential':np.array([-60.]),'neural.inh.gaba_a':np.array([.2])}
    theta={'tau_membrane_s':.015,'g_leak':1.,'e_gaba_a_mv':-70.}
    old=float(legacy.brain.source_inhibition_law(sample,theta)['neural.exc.potential'][0])
    new=float(brain.source_inhibition_law(sample,theta)['neural.exc.potential'][0])
    result=skin.step(.01,{'time_s':0.,'contacts':[{'id':'fixture','force_n':[0,0,-1]}],'skin_temperature_C':33.})
    output=brain.step(.01,sensory_inputs_hz=result['sensory_inputs_hz'])
    assert np.isfinite(output['regional_state']['activity_hz']).all()
    try:IBMBackend()
    except RuntimeError:pass
    else:raise AssertionError('Source identity hot-swap accepted')
    print(json.dumps({'package_sha256':pin.package_sha256,'inhibition_probe_old_mV_s':old,'inhibition_probe_candidate_mV_s':new,'shared_identity':brain.source_identity['package_sha256']==skin.audit['package_sha256']==backend.identity['package_sha256']}))

if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--probe':probe(sys.argv[2])
    else:unittest.main()
