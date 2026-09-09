"""Physical support filtering preserves first-hit provenance and missing cells."""
import unittest
import numpy as np
from scripts.build_supine_surface_contact import restrict_posterior_envelope


class Tests(unittest.TestCase):
    def test_excluded_first_hit_is_omitted_without_substitution(self):
        original={'points_source_m':np.array([[0.,0,0],[1.,0,1]]),
                  'face_indices':np.array([165251,7]),'area_m2':np.array([.000025,.000025])}
        kept,missing=restrict_posterior_envelope(original,[7,8])
        np.testing.assert_array_equal(kept['face_indices'],[7])
        np.testing.assert_array_equal(kept['original_quadrature_indices'],[1])
        np.testing.assert_array_equal(kept['points_source_m'],original['points_source_m'][[1]])
        self.assertEqual(missing['source_face_indices'],[165251])
        self.assertAlmostEqual(missing['projected_area_m2'],.000025)
        self.assertNotIn(8,kept['face_indices'])

    def test_no_exclusions_preserves_rows_exactly(self):
        original={'points_source_m':np.array([[.1,.2,.3],[.4,.5,.6]]),
                  'face_indices':np.array([2,4]),'area_m2':np.array([.1,.2])}
        kept,missing=restrict_posterior_envelope(original,[4,2])
        for key in original:np.testing.assert_array_equal(kept[key],original[key])
        self.assertEqual(missing['original_quadrature_indices'],[])
        self.assertEqual(missing['projected_area_m2'],0)


if __name__=='__main__':unittest.main()
