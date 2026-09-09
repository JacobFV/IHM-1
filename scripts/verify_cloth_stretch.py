"""Material strain-constraint checks, including incompatible physical boundaries."""
import unittest
from types import SimpleNamespace
import numpy as np
from ihm.assembly.cloth_stretch import ClothStretchConstraint
from ihm.assembly.environment_dynamics import SpringMesh


def chain(points, fixed=None):
    x=np.asarray(points,float)
    a=np.arange(len(x)-1);b=a+1
    return SimpleNamespace(x=x,rest=x.copy(),v=np.zeros_like(x),mass=np.ones(len(x)),
        fixed=np.zeros(len(x),bool) if fixed is None else np.asarray(fixed,bool),
        a=a,b=b,edges=np.column_stack([a,b,np.ones(len(a))]),length=np.ones(len(a)),stiffness=np.ones(len(a))*45)


class Tests(unittest.TestCase):
    def test_extreme_extension_preserves_momentum_and_com(self):
        mesh=chain([[0,0,0],[5,0,0]])
        mesh.mass=np.array([2.,3.]);mesh.v[:]=[[-1,2,0],[2,3,0]]
        com=np.sum(mesh.mass[:,None]*mesh.x,axis=0)
        momentum=np.sum(mesh.mass[:,None]*mesh.v,axis=0)
        result=ClothStretchConstraint(mesh).resolve(mesh)
        self.assertTrue(result['converged']);self.assertAlmostEqual(result['max_extension_ratio'],1.12)
        np.testing.assert_allclose(np.sum(mesh.mass[:,None]*mesh.x,axis=0),com)
        np.testing.assert_allclose(np.sum(mesh.mass[:,None]*mesh.v,axis=0),momentum)
        self.assertLess(result['kinetic_energy_change_j'],0)
        self.assertLess(result['spring_energy_change_j'],0)

    def test_projection_does_not_create_motion_or_change_rest(self):
        mesh=chain([[0,0,0],[5,0,0]])
        rest=mesh.length.copy()
        ClothStretchConstraint(mesh).resolve(mesh,gravity=[0,0,-9.81])
        np.testing.assert_array_equal(mesh.v,0)
        np.testing.assert_array_equal(mesh.length,rest)

    def test_compression_and_bending_remain_free(self):
        mesh=chain([[0,0,0],[.3,0,0],[.3,.2,0]])
        before=mesh.x.copy()
        result=ClothStretchConstraint(mesh).resolve(mesh)
        self.assertEqual(result['projected_edges'],0)
        np.testing.assert_array_equal(mesh.x,before)

    def test_fixed_anchor_reaction(self):
        mesh=chain([[0,0,0],[2,0,0]],fixed=[True,False]);mesh.v[1,0]=2
        momentum=mesh.v.sum(axis=0)
        result=ClothStretchConstraint(mesh).resolve(mesh)
        np.testing.assert_array_equal(mesh.x[0],[0,0,0])
        np.testing.assert_allclose(mesh.v.sum(axis=0)+result['support_impulse_ns'],momentum)

    def test_incompatible_fixed_anchors_report_failure_without_rest_rewrite(self):
        mesh=chain([[0,0,0],[2.5,0,0],[5,0,0]],fixed=[True,False,True])
        result=ClothStretchConstraint(mesh,iterations=30).resolve(mesh)
        self.assertFalse(result['converged']);self.assertGreater(result['violating_edges'],0)
        np.testing.assert_array_equal(mesh.x[[0,2]],[[0,0,0],[5,0,0]])
        np.testing.assert_array_equal(mesh.length,1)

    def test_actual_mesh_local_fivefold_strain_reduces_and_replays(self):
        mesh=SpringMesh('cloth','cloth',[-.42,-.9,.164],[.42,-.04,.164],1.5)
        mesh.x[-1,2]+=.2
        solver=ClothStretchConstraint(mesh,iterations=120)
        before=np.max(np.linalg.norm(mesh.x[solver.b]-mesh.x[solver.a],axis=1)/solver.rest)
        x=mesh.x.copy();v=mesh.v.copy()
        result=solver.resolve(mesh)
        self.assertGreater(before,5)
        self.assertTrue(result['converged'],result)
        expected=mesh.x.copy()
        mesh.x[:]=x;mesh.v[:]=v
        self.assertEqual(solver.resolve(mesh),result)
        np.testing.assert_array_equal(mesh.x,expected)


if __name__=='__main__':unittest.main()
