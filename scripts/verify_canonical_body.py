"""Canonical coupling checks: identity, clocks, projection and provenance."""
from pathlib import Path
import sys, json, tempfile, unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

class BodyTests(unittest.TestCase):
    def test_native_clock_and_projection(self):
        from ihm.assembly.body import read_native, build_bindings, project_volumes
        a=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
        if not (ROOT/'data/derived/canonical/native_baseline_v1/summary.json').exists():self.skipTest('Native baseline has not been executed')
        series=read_native(ROOT/'data/derived/canonical/native_baseline_v1')
        self.assertEqual(len(series['time_s']),1500)
        np.testing.assert_allclose(np.diff(series['time_s']),.02,atol=1e-9)
        bindings=build_bindings(a)
        self.assertTrue(any(b['entity_id']=='body-bp3d-FJ2422' for b in bindings))
        for row in [series['compartments'][0],series['compartments'][-1]]:
            ratios=project_volumes(bindings,row,series['initial_compartments'])
            for b in bindings:
                self.assertAlmostEqual(ratios[b['entity_id']],row[b['column']]/series['initial_compartments'][b['column']])
        self.assertEqual(set(project_volumes(bindings,series['initial_compartments'],series['initial_compartments']).values()),{1.})
        self.assertTrue(all(b['evidence_kind']=='synthesized_transfer' for b in bindings))

    def test_clock_mismatch_rejected(self):
        from ihm.assembly.body import read_native
        import shutil
        src=ROOT/'data/derived/canonical/native_baseline_v1'
        if not (src/'summary.json').exists():self.skipTest('Native baseline has not been executed')
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)
            for name in ['native_multisystem.csv','body_compartments.csv','summary.json']:shutil.copyfile(src/name,p/name)
            target=p/'body_compartments.csv';text=target.read_text();target.write_text(text.replace('0.02,','0.021,',1))
            from ihm.assembly.body import digest
            summary=json.loads((p/'summary.json').read_text());summary['compartment_telemetry']['csv_sha256']=digest(target);(p/'summary.json').write_text(json.dumps(summary))
            with self.assertRaisesRegex(ValueError,'different clocks'):read_native(p)

    def test_changed_native_receipt_rejected(self):
        from ihm.assembly.body import read_native
        import shutil
        src=ROOT/'data/derived/canonical/native_baseline_v1'
        if not (src/'summary.json').exists():self.skipTest('Native baseline not executed')
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)
            for name in ['native_multisystem.csv','body_compartments.csv','summary.json']:shutil.copyfile(src/name,p/name)
            target=p/'native_multisystem.csv';target.write_text(target.read_text().replace('0.980295','0.200000',1))
            with self.assertRaisesRegex(ValueError,'receipt'):read_native(p)

    def test_materialization_integrity(self):
        from ihm.assembly.body import CanonicalBody
        body=CanonicalBody.from_workspace(ROOT)
        self.assertEqual(body.describe()['entity_count'],len(body.entities))
        # The entity count is not a literal here. It was 2408, then 2403 after five
        # duplicate-authored BodyParts3D surfaces collapsed, and 4000 once the display
        # promotion landed; a literal goes stale on the next legitimate change and then
        # gets edited to whatever the build now says, which checks nothing. What does not
        # go stale is the arithmetic the count has to satisfy: every entity has exactly one
        # evidence kind, the acquired rows are the full BP scaffold minus the recorded
        # collapse, and the promoted rows are exactly the promote decisions the recorded
        # ledger carries. A silent drop, a silent addition or a divergence from the ledger
        # still fails all three.
        from ihm.assembly.anatomy import BP_SOURCE_SURFACE_COUNT
        anatomy=json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
        counts=anatomy['counts'];collapse=anatomy['duplicate_surface_collapse']['dropped']
        promotion=anatomy['display_structure_promotion']
        self.assertEqual(len(body.entities),len(anatomy['entities']))
        self.assertEqual(len(body.entities),counts['entities'])
        self.assertEqual(sum(counts['evidence_kinds'].values()),len(body.entities))
        self.assertEqual(sum(counts['roles'].values()),len(body.entities))
        self.assertEqual(sum(counts['systems'].values()),len(body.entities))
        self.assertEqual(counts['evidence_kinds']['source_geometry'],
                         BP_SOURCE_SURFACE_COUNT-len(collapse))
        self.assertEqual(len(collapse),counts['collapsed_duplicate_surfaces'])
        self.assertEqual(counts['promoted_display_structures'],
                         promotion['decision_counts']['promote'])
        self.assertEqual(len(promotion['ids']),counts['promoted_display_structures'])
        self.assertEqual(len(promotion['promoted']),counts['promoted_display_structures'])
        self.assertIsNone(body.describe()['certainty']['calibrated_confidence_score'])
        self.assertFalse(body.describe()['validated_digital_twin'])
        self.assertTrue(all(b['entity_id'] in body.entities for b in body.payload['volume_bindings']))
        certainty=body.certainty('body-skin-dermis')
        self.assertEqual(certainty['evidence_kind'],'synthesized_layer')
        self.assertIn('gpu_coordinate_rounding_bound_m',certainty['display_numerics'])
        with self.assertRaises(ValueError):body.certainty('../not-an-entity')

    def test_recorded_execution(self):
        p=ROOT/'data/derived/canonical/trajectory.json'
        if not p.exists():self.skipTest('Full canonical trajectory not materialized')
        data=json.loads(p.read_text());self.assertEqual(data['clock']['end_s'],30.)
        self.assertEqual(data['clock']['native_samples'],1500)
        self.assertGreaterEqual(len(data['frames']),300)
        self.assertLess(data['audit']['maximum_internal_force_residual_n'],1e-6)
        self.assertLess(data['audit']['maximum_internal_torque_residual_nm'],1e-6)
        for frame in data['frames']:
            self.assertAlmostEqual(frame['time_s'],frame['brain']['time_s'])
            self.assertAlmostEqual(frame['time_s'],frame['physiology']['#Time(s)'])
            self.assertFalse(frame['brain']['autonomic_commands']['applied_to_body'])
            for identity,ratio in frame['volume_ratios'].items():
                self.assertAlmostEqual(float(np.linalg.det(frame['entities'][identity]['deformation_gradient'])),ratio,places=9)
        self.assertGreater(np.ptp([f['compartments']['LeftHeartVolume(mL)'] for f in data['frames']]),10.)
        self.assertGreater(np.ptp([f['compartments']['LeftLungPulmonaryGasVolume(mL)'] for f in data['frames']]),100.)
        spectra=json.loads(p.with_name('trajectory-spectra.json').read_text())
        self.assertEqual(spectra['model_id'],'ihm-body')
        self.assertEqual(len(spectra['variables']),11)
        self.assertAlmostEqual(spectra['sample_interval_s'],.02)

if __name__=='__main__':unittest.main()
