"""Pinned source/control-path regressions; no native compilation or run."""
import math,tempfile,unittest
from pathlib import Path
from patch_signed_cardiovascular_reader import SOURCE,patch,prepare,PRIOR,READER
class Tests(unittest.TestCase):
    def test_reader_added_only_to_nondrug_else_with_no_physical_writes(self):
        before=SOURCE.read_text();after=patch(before)
        self.assertEqual(after.count('MetabolicToneResponse();'),before.count('MetabolicToneResponse();'))
        extra=READER.split('} else if',1)[1];self.assertNotIn('SetValue',extra)
        self.assertIn('ihm_signed::active() && m_data.GetEnergy().HasTotalMetabolicRate()',extra)
        self.assertIn('Convert(1.0, PowerUnit::W, PowerUnit::kcal_Per_day), 0)',extra)
    def test_guard_rejects_changed_source(self):
        with self.assertRaises(ValueError):patch(SOURCE.read_text()+'\n')
    def test_original_tone_law_and_energy_owners_untouched(self):
        before=SOURCE.read_text();after=patch(before,vascular=True)
        a='void Cardiovascular::MetabolicToneResponse()';b='void Cardiovascular::TuneCircuit()'
        self.assertEqual(before[before.index(a):before.index(b)],after[after.index(a):after.index(b)])
        self.assertNotIn('GetExerciseMeanArterialPressureDelta',PRIOR);self.assertNotIn('GetNextHeatSource',PRIOR)
        self.assertIn('signedMuscle->delta_m_W != 0',PRIOR);self.assertIn('requires native regional baseline resets',PRIOR)
    def test_prepare_keeps_distinct_unbuilt_receipts(self):
        with tempfile.TemporaryDirectory() as d:
            out=prepare(Path(d)/'overlay');self.assertTrue((out/'reader_only/Cardiovascular.cpp').exists());self.assertTrue((out/'vascular_prior/native_signed_vascular_prior.h').exists())
    def test_transferred_prior_algebra_baseline_positive_negative(self):
        # Algebraic audit only. Exact C++ implementation needs separately granted probe.
        for background,basal,delta in [(80,80,0),(90,80,5),(90,80,-5),(90,80,-90)]:
            q=lambda w:(1.5*w/basal+5.5)/7
            normalized=q(background+delta)/q(background)
            ratio=1+1.5*delta/(1.5*background+5.5*basal)
            self.assertAlmostEqual(normalized,ratio);self.assertGreater(ratio,0)
            if delta==0:self.assertEqual(ratio,1)
            elif delta>0:self.assertLess(1/ratio**2,1)
            else:self.assertGreater(1/ratio**2,1)
if __name__=='__main__':unittest.main()
