"""Fresh-process active default and explicit legacy source selection; no native runs."""
from pathlib import Path
import subprocess,sys,unittest
ROOT=Path(__file__).resolve().parents[1]
class Tests(unittest.TestCase):
    def test_default_one_identity_and_explicit_legacy_hotswap_rejection(self):
        code='''
from ihm.brain.active_source import resolve_source,ACTIVE_PACKAGE_SHA256
from ihm.brain.ibm_backend import IBMBackend
from ihm.assembly.sensorimotor import SensorimotorController
from ihm.assembly.cutaneous_feedback import CutaneousFeedback
from scripts.verify_cutaneous_feedback import sites
from pathlib import Path
root=Path.cwd();pin=resolve_source(root)
brain=SensorimotorController.from_root(root).brain
backend=IBMBackend()
skin=CutaneousFeedback(root,sites=sites(),recruitment_hz_per_response=.1)
assert brain.source_identity==pin.to_dict()==backend.source_pin.to_dict()==skin.audit['source_pin']
assert pin.package_sha256==ACTIVE_PACKAGE_SHA256==skin.audit['package_sha256']
try:IBMBackend(source_pin=None)
except RuntimeError:pass
else:raise AssertionError('Cross-source hot swap accepted')
'''
        subprocess.run([sys.executable,'-c',code],cwd=ROOT,check=True,capture_output=True,text=True)
    def test_explicit_legacy_remains_available_in_fresh_process(self):
        code='''
from pathlib import Path
from ihm.brain.ibm_backend import IBMBackend,PINNED_PACKAGE_SHA256
from ihm.assembly.sensorimotor import SensorimotorController
backend=IBMBackend(source_pin=None)
c=SensorimotorController.from_root(Path.cwd(),source_pin=None)
assert backend.identity['package_sha256']==PINNED_PACKAGE_SHA256
assert c.brain.source_identity['selection']=='legacy independently preserved regional law'
'''
        subprocess.run([sys.executable,'-c',code],cwd=ROOT,check=True,capture_output=True,text=True)
    def test_missing_active_registry_fails_closed_legacy_is_explicit(self):
        from ihm.brain.active_source import resolve_source
        import tempfile
        with tempfile.TemporaryDirectory() as p:
            self.assertIsNone(resolve_source(p,None))
            with self.assertRaises(ValueError):resolve_source(p)
if __name__=='__main__':unittest.main()
