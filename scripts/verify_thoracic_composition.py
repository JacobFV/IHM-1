"""Independent retained cervical/thoracic reference moment ledger checks."""
import hashlib,json,unittest
from pathlib import Path
import numpy as np
from scripts.build_thoracic_composition_plan import build
ROOT=Path(__file__).resolve().parents[1]

class CompositionTests(unittest.TestCase):
    def test_receipts_and_builder_reproducibility(self):
        plan=json.loads((ROOT/'data/research/thoracic_mechanism/native_composition_v1/plan.json').read_bytes())
        self.assertEqual(plan,build())
        for receipt in plan['source_receipts'].values():
            raw=(ROOT/receipt['path']).read_bytes();self.assertEqual(hashlib.sha256(raw).hexdigest(),receipt['sha256'])
        self.assertFalse(plan['native_activation_allowed'])

    def test_reference_inertia_closure_in_one_parent_frame(self):
        p=build();neck=p['cervical_aggregate'];thorax=p['thoracic_subsystem_aggregate']
        self.assertEqual(len(p['cervical_bodies']),9);self.assertEqual(len(p['thoracic_material_partitions']),48)
        mass=neck['mass_kg']+thorax['mass_kg'];first=np.array(neck['H_first_moment_kg_m'])+thorax['H_first_moment_kg_m']
        second=np.array(neck['Q_second_moment_kg_m2'])+thorax['Q_second_moment_kg_m2'];center=first/mass
        central=second-mass*np.outer(center,center);inertia=np.trace(central)*np.eye(3)-central
        self.assertAlmostEqual(mass,27.654676965260336,places=12)
        np.testing.assert_allclose(center,[-.03,.32,0],atol=1e-14)
        np.testing.assert_allclose(inertia,np.diag([1.520014507406313,.7788205903211117,1.4755841071015214]),atol=2e-12)
        self.assertGreater(np.linalg.eigvalsh(central)[0],0)
        self.assertAlmostEqual(neck['mass_kg'],7.418919568221101,places=12)
        self.assertAlmostEqual(thorax['mass_kg'],20.235757397039237,places=12)

if __name__=='__main__':unittest.main()
