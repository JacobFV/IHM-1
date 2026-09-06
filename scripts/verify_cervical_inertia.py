"""Analytic geometric moments and exclusive-inertia accounting; no native run."""
import itertools,unittest
import numpy as np
from ihm.assembly.cervical_inertia import convex_geometry_prior,segment_prior,partition_body,combine_bodies,transform_body


class Tests(unittest.TestCase):
    def box(self):
        return np.array(list(itertools.product((-1.,1.),(-2.,2.),(-3.,3.))))

    def test_convex_volume_moments_match_analytic_box(self):
        prior=convex_geometry_prior(self.box(),np.array([[0,1,2]]))
        self.assertAlmostEqual(prior['volume_m3'],48.)
        np.testing.assert_allclose(prior['centroid_m'],[0,0,0],atol=1e-14)
        np.testing.assert_allclose(prior['covariance_m2'],np.diag([1/3,4/3,3]),atol=1e-14)
        self.assertFalse(prior['source_topology']['closed_oriented_edge_manifold'])
        self.assertTrue(prior['hull_topology']['closed_oriented_edge_manifold'])

    def test_moments_are_rigid_transform_covariant(self):
        v=self.box();rotation=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        original=convex_geometry_prior(v,[])
        shifted=convex_geometry_prior(v@rotation.T+[7,-9,11],[])
        np.testing.assert_allclose(shifted['centroid_m'],[7,-9,11],atol=1e-13)
        np.testing.assert_allclose(shifted['covariance_m2'],rotation@original['covariance_m2']@rotation.T,atol=1e-13)

    def test_soft_shell_mass_and_second_moments_are_conserved(self):
        geometry=convex_geometry_prior(self.box(),[])
        body=segment_prior(geometry,192.,bone_density_kg_m3=2.,soft_density_kg_m3=1.)
        scale=3**(1/3)
        self.assertAlmostEqual(body['components']['bone_proxy_mass_kg'],96.)
        self.assertAlmostEqual(body['components']['soft_proxy_mass_kg'],96.)
        self.assertAlmostEqual(body['components']['soft_shell_scale'],scale)
        expected_second=(2*48+48*(scale**5-1))*np.diag([1/3,4/3,3])
        np.testing.assert_allclose(body['inertia_kg_m2'],np.trace(expected_second)*np.eye(3)-expected_second,rtol=1e-13)
        with self.assertRaisesRegex(ValueError,'bone proxy exceeds'):segment_prior(geometry,1.)

    def test_exclusive_partition_reconstructs_mass_first_and_second_moments(self):
        children=[{'mass_kg':2.,'center_m':[-.2,.3,.1],'inertia_kg_m2':np.diag([.1,.12,.14]).tolist()},
                  {'mass_kg':1.,'center_m':[.1,-.1,.2],'inertia_kg_m2':np.diag([.02,.03,.04]).tolist()}]
        remainder={'mass_kg':7.,'center_m':[0,.05,0],'inertia_kg_m2':np.diag([.3,.4,.5]).tolist()}
        parent=combine_bodies(children+[remainder]);residual=partition_body(parent,children)
        self.assertAlmostEqual(residual['mass_kg'],remainder['mass_kg'])
        np.testing.assert_allclose(residual['center_m'],remainder['center_m'],atol=1e-14)
        np.testing.assert_allclose(residual['inertia_kg_m2'],remainder['inertia_kg_m2'],atol=1e-14)
        with self.assertRaisesRegex(ValueError,'physical inertia'):partition_body(remainder,[children[0]])

    def test_positive_definite_but_impossible_inertia_is_rejected(self):
        bad={'mass_kg':1.,'center_m':[0,0,0],'inertia_kg_m2':np.diag([.02,.08,.02]).tolist()}
        with self.assertRaisesRegex(ValueError,'physical inertia'):combine_bodies([bad])

    def test_reflection_is_not_a_registration(self):
        body={'mass_kg':1.,'center_m':[0,0,0],'inertia_kg_m2':np.eye(3).tolist()}
        reflection=np.eye(4);reflection[0,0]=-1
        with self.assertRaisesRegex(ValueError,'proper rigid'):transform_body(body,reflection)


if __name__=='__main__':unittest.main()
