"""Bounded anatomy recipe and geometric work checks; no native engine."""
import unittest
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_anatomy import surface_mass_prior,hinge_motion,nearest_node_pair,descent_weights,diaphragm_motion,validate_manifest

class ThoracicTests(unittest.TestCase):
    def test_exact_lamina_moments(self):
        v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0.]])
        p=surface_mass_prior(v,[[0,1,2],[0,2,3]],2.)
        np.testing.assert_allclose(p['center_m'],[.5,.5,0])
        np.testing.assert_allclose(p['inertia_kg_m2'],np.diag([1/6,1/6,1/3]),atol=1e-14)
        self.assertEqual(p['source_area_m2'],1)
    def test_hinge_jacobian_and_virtual_work(self):
        v=np.array([[1,0,0],[0,2,1.]])
        origin=[.1,.2,.3];axis=[0,1,0];q=.2
        x,j=hinge_motion(v,origin,axis,q)
        h=1e-7
        fd=(hinge_motion(v,origin,axis,q+h)[0]-hinge_motion(v,origin,axis,q-h)[0])/(2*h)
        np.testing.assert_allclose(j,fd,atol=2e-9)
        force=np.array([[1,2,3],[-1,3,2.]])
        self.assertAlmostEqual(float(np.sum(force*j))*.3,float(np.sum(force*(j*.3))))
    def test_exact_node_binding_and_anchor_zero(self):
        a=np.array([[0,0,0],[2,0,0],[0,2,0.]])
        b=a+[0,0,1]
        pair=nearest_node_pair(a,b)
        self.assertEqual(pair['gap_m'],1)
        self.assertEqual(pair['source_vertex'],0)
        w=descent_weights(a,[0,1])
        np.testing.assert_array_equal(w,[0,0,1])
    def test_diaphragm_shape_derivative(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,1.]])
        w=[0,.5,1]
        x,j=diaphragm_motion(v,w,[0,-1,0],.003)
        np.testing.assert_allclose(x-v,.003*j,atol=1e-15)
        np.testing.assert_array_equal(x[0],v[0])

    def test_retained_bindings_and_mass_ledger(self):
        path=Path(__file__).resolve().parents[1]/'data/research/thoracic_anatomy/v3/manifest.json'
        result=validate_manifest(path)
        self.assertEqual(result['source_entities'],63)
        self.assertEqual(result['verified_attachment_pairs'],202)

    def test_invalid_geometry_and_axis(self):
        with self.assertRaises(ValueError):surface_mass_prior([[0,0,0]]*3,[[0,1,2]],1)
        with self.assertRaises(ValueError):hinge_motion([[0,0,0]],[0,0,0],[0,2,0],0)

if __name__=='__main__':unittest.main()
