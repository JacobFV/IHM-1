#!/usr/bin/env python3
"""Bounded native Energy donor/waste patch and engine lifecycle source checks."""
import unittest
from build_biogears_regional_sweat_variant import corrected_source,ENERGY_SOURCE_PIN
from verify_native_regional_skin import RUNTIME,sha
from verify_native_regional_skin_engine import step_source,ENGINE


class Checks(unittest.TestCase):
    def test_paired_sweat_source_only(self):
        source=RUNTIME/'variants/whole_body_integrity_signed_muscle_v2/Energy.cpp'
        self.assertEqual(sha(source),ENERGY_SOURCE_PIN);before=source.read_text();after=corrected_source(before)
        restored=after.removeprefix('#define IHM_REGIONAL_SPECIES_IMPLEMENTATION\n#include "native_regional_species.h"\n')
        for index,(species,var,z) in enumerate([('Sodium','sodium',1),('Potassium','potassium',1),('Chloride','chloride',-1)]):
            old=f'  m_Skin{species}->GetMass().IncrementValue(-{var}Lost_mg, MassUnit::mg);\n  Get{species}LostToSweat().IncrementValue({var}Lost_mg, MassUnit::mg);'
            new=f'  ihm_regional::withdraw_skin_ion(*m_data.GetCompartments().GetLiquidCompartment(BGE::ExtravascularCompartment::SkinExtracellular), *m_Skin{species}, {var}Lost_mg, Get{species}LostToSweat(), {index}, {z});'
            self.assertEqual(restored.count(new),1);restored=restored.replace(new,old)
        self.assertEqual(restored,before)
        with self.assertRaises(ValueError):corrected_source(before+before)

    def test_exact_native_lifecycle_plus_two_hooks(self):
        before=ENGINE.read_text();after=step_source(before)
        signature='bool BioGearsEngine::AdvanceModelTime(bool appendDataTrack)'
        start=before.index(signature);end=before.index('\n//-------------------------------------------------------------------------------',start)
        restored=after.rstrip('\n').replace('bool AdvanceModelTime(bool appendDataTrack=false) override',signature)
        restored=restored.replace('\n  if(regional)regional->after_preprocess();','').replace('\n  if(regional)regional->after_postprocess();','')
        self.assertEqual(restored,before[start:end])


if __name__=='__main__':unittest.main()
