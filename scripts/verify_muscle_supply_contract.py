"""Light local contract tests; no native process, no fabricated force model."""
from dataclasses import replace
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.coupling.muscle_supply import TrialBinding,MuscleDemandTrial,SupplyPreview,assess_supply

class SupplyContractTests(unittest.TestCase):
    def setUp(self):
        self.binding=TrialBinding('phys','mech','control','reference','build',1,0.,.02)
        self.trial=MuscleDemandTrial(self.binding,.2,.3,-.1)
        self.preview=SupplyPreview(self.binding,.2,1.,1.2,1.2,0.)
    def test_no_native_preview_is_unavailable(self):
        self.assertEqual(assess_supply(self.trial,None).reason,'native_supply_preview_unavailable')
    def test_supplied_eccentric_candidate(self):
        self.assertTrue(assess_supply(self.trial,self.preview).admissible)
    def test_unmet_demands_retry_without_changing_candidate(self):
        preview=replace(self.preview,provided_muscle_j=1.,unmet_muscle_j=.2)
        self.assertFalse(assess_supply(self.trial,preview).admissible)
        self.assertEqual(self.trial.delta_signed_work_j,-.1)
    def test_signed_decrease_is_not_supply_credit(self):
        trial=replace(self.trial,delta_metabolic_j=-.2,delta_heat_j=-.1)
        preview=replace(self.preview,delta_metabolic_j=-.2,requested_muscle_j=.8,provided_muscle_j=.8)
        self.assertTrue(assess_supply(trial,preview).admissible)
        with self.assertRaises(ValueError):assess_supply(trial,replace(preview,provided_muscle_j=1.))
    def test_every_binding_dimension_invalidates_preview(self):
        for name,value in [('physiology_state_id','other'),('mechanics_state_id','other'),('control_id','other'),('reference_id','other'),('native_build_id','other'),('sequence',2),('start_s',.001),('end_s',.03)]:
            with self.subTest(name=name),self.assertRaises(ValueError):
                assess_supply(self.trial,replace(self.preview,binding=replace(self.binding,**{name:value})))
    def test_invalid_ledgers_rejected(self):
        for change in ({'requested_muscle_j':1.3},{'provided_muscle_j':-1.},{'unmet_muscle_j':float('nan')},{'delta_metabolic_j':.1},{'origin':'previous-step-unmet'}):
            with self.subTest(change=change),self.assertRaises(ValueError):assess_supply(self.trial,replace(self.preview,**change))
    def test_inconsistent_heat_work_rejected(self):
        with self.assertRaises(ValueError):assess_supply(replace(self.trial,delta_heat_j=.4),self.preview)
    def test_invalid_binding_and_interval_rejected(self):
        for change in ({'sequence':True},{'reference_id':''},{'end_s':0.},{'start_s':float('nan')}):
            with self.subTest(change=change),self.assertRaises(ValueError):assess_supply(replace(self.trial,binding=replace(self.binding,**change)),None)
if __name__=='__main__':unittest.main()
