#!/usr/bin/env python3
"""Small source-transformation checks; never compiles or loads a native patient."""
import unittest
from build_biogears_regional_skin_variant import regional_source, LOOKUPS
from build_biogears_shared_donor_variant import corrected_source
from verify_native_regional_skin import SOURCE


class Checks(unittest.TestCase):
    def test_only_native_lookup_owner_changes(self):
        donor=SOURCE/'projects/biogears/libBiogears/src/engine/Systems/Diffusion.cpp'
        before=corrected_source(donor.read_text())
        after=regional_source(before)
        restored=after
        for name in LOOKUPS:
            new='m_data.GetCircuits().GetFluidPath('+name+')'
            old='m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath('+name+')'
            self.assertEqual(after.count(new),1)
            restored=restored.replace(new,old)
        self.assertEqual(restored,before)
        self.assertIn('combinedWithdrawal_ug',after)
        self.assertIn('DistributeMassbyMassWeighted(source, sub, -massToMove_ug',after)

    def test_missing_or_duplicate_anchors_fail(self):
        for text in ['', '\n'.join('m_data.GetCircuits().GetActiveCardiovascularCircuit().GetPath('+name+')' for name in (*LOOKUPS,LOOKUPS[0]))]:
            with self.assertRaises(ValueError):regional_source(text)


if __name__=='__main__':unittest.main()
