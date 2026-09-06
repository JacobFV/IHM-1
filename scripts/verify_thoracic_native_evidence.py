"""Retained real native extraction receipts and narrowly stated proof scope."""
import gzip,hashlib,json,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

class NativeEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'data/research/thoracic_mechanism/native_operator_proof_v1'
        cls.manifest=json.loads((cls.path/'manifest.json').read_bytes())
    def test_retained_actual_bytes_and_failure(self):
        for row in self.manifest['artifacts']:
            raw=gzip.decompress((self.path/row['retained_copy']).read_bytes())
            self.assertEqual(hashlib.sha256(raw).hexdigest(),row['sha256']);self.assertEqual(len(raw),row['bytes'])
        error=gzip.decompress((self.path/'v1-run.log.gz').read_bytes()).decode()
        self.assertIn("object type 'ExcitationPorts' is not a registered Object",error)
        build=json.loads(gzip.decompress((self.path/'v2-manifest.json.gz').read_bytes()))
        self.assertEqual(build['status'],'complete');self.assertTrue(build['native_snapshot_run'])
        self.assertTrue(any('libSimTKsimbody.so' in r['path'] for r in build['inputs']))
    def test_actual_snapshot_scope_and_replacement(self):
        raw=gzip.decompress((self.path/'v2-snapshot.json.gz').read_bytes());snapshot=json.loads(raw);proof=self.manifest['proof']
        self.assertEqual(hashlib.sha256(raw).hexdigest(),proof['snapshot_sha256'])
        for key in ['native_integrated','native_constraints_projected','private_force_control_state_restored','force_or_coupled_state_equivalence_claimed']:self.assertIs(snapshot[key],False)
        self.assertIs(snapshot['private_callbacks_throw_if_requested'],True);self.assertIs(snapshot['zero_applied_force_inverse_dynamics'],True)
        self.assertEqual(snapshot['body_count'],22);self.assertEqual(snapshot['source_muscle_count'],92)
        self.assertEqual(snapshot['declared_native_constraint_components'],2)
        self.assertEqual(len(snapshot['u']),33);self.assertEqual(np.asarray(snapshot['mass_matrix']).shape,(33,33))
        self.assertEqual(len(proof['active_eigenvalues']),63);self.assertGreater(min(proof['active_eigenvalues']),0)
        self.assertLess(proof['reference_native_mass_block_max_error'],2e-10);self.assertLess(proof['reference_native_bias_max_error'],2e-10)
        self.assertGreater(proof['native_to_new_internal_cross_mass_norm'],1.)

if __name__=='__main__':unittest.main()
