"""Finite-owner first-impact acceptance, no native engine or body geometry."""
from pathlib import Path
import sys, unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.garment_vertex_face_step import step_vertex_face_impact

class Body:
    def __init__(self,x,v,m):
        self.position_m=np.array(x,float);self.velocity_m_s=np.array(v,float);self.mass_kg=np.array(m,float)
        self.time_s=0.;self.max_explicit_dt_s=2.;self.triangles=np.array([[0,1,2]],int)
    def elastic_forces(self,x):return np.zeros_like(x),0.

def fixture(tangent=0):
    # Centroid target weights yield equal triangle nodal velocities after impact.
    tri=[[-1,-1,0],[1,-1,0],[0,2,0]]
    return dict(cloth=Body([[0,0,.1]],[[tangent,0,-1]],[1]),skin=Body(tri,[[0,0,0]]*3,[1,1,1]))
def call(b,**kw):
    kw.setdefault('max_energy_defect_j',1e-12)
    return step_vertex_face_impact(b,.2,cloth_owner='cloth',surface_owner='skin',node_ids=[0],face_ids=[0],friction_static=.4,friction_kinetic=.3,**kw)

class Checks(unittest.TestCase):
    def test_first_event_remaining_drift_and_ledger(self):
        b=fixture();r=call(b)
        self.assertAlmostEqual(r['event']['time_s'],.1)
        self.assertAlmostEqual(b['cloth'].position_m[0,2],-.025)
        np.testing.assert_allclose(b['skin'].position_m[:,2],-.025,atol=1e-14)
        self.assertAlmostEqual(r['contact_dissipation_j'],.375)
        self.assertLess(abs(r['numerical_energy_defect_j']),1e-13)
        self.assertLess(np.linalg.norm(r['momentum_residual_ns']),1e-13)
        self.assertLess(np.linalg.norm(r['angular_impulse_residual_nms']),1e-13)
        self.assertEqual(b['cloth'].time_s,.2)
    def test_sliding_friction_existing_kernel(self):
        b=fixture();b['cloth'].position_m[0,0]=-.1;b['cloth'].velocity_m_s[0,0]=1
        r=call(b)
        self.assertGreater(r['contact_dissipation_j'],.375)
        self.assertLess(abs(r['numerical_energy_defect_j']),1e-13)
        np.testing.assert_allclose(r['paired_contact_impulse_ns'],0,atol=1e-14)
    def test_unsupported_or_multiple_event_atomicity(self):
        b=fixture();b['cloth'].position_m[0,:2]=[0,-1]
        before={k:v.position_m.copy() for k,v in b.items()}
        with self.assertRaises(ValueError):call(b)
        for k,v in b.items():np.testing.assert_array_equal(v.position_m,before[k]);self.assertEqual(v.time_s,0)
    def test_event_uses_common_owner_clock(self):
        b=fixture()
        for owner in b.values():owner.time_s=10.
        r=call(b)
        self.assertAlmostEqual(r['event']['time_s'],10.1)
        self.assertAlmostEqual(r['event']['elapsed_time_s'],.1)
        self.assertAlmostEqual(r['time_s'],10.2)
    def test_rigid_and_galilean_invariance(self):
        a=fixture();b=fixture();rotation=np.array([[0,0,1],[1,0,0],[0,1,0]]);boost=np.array([.7,-.3,.2]);shift=np.array([2.,-4,1])
        for owner in b.values():
            owner.position_m=owner.position_m@rotation+shift
            owner.velocity_m_s=owner.velocity_m_s@rotation+boost
        ra=call(a);rb=call(b)
        self.assertAlmostEqual(ra['contact_dissipation_j'],rb['contact_dissipation_j'])
        self.assertAlmostEqual(ra['event']['time_s'],rb['event']['time_s'])
        for k in a:
            np.testing.assert_allclose(b[k].position_m,a[k].position_m@rotation+shift+.2*boost,atol=1e-13)
            np.testing.assert_allclose(b[k].velocity_m_s,a[k].velocity_m_s@rotation+boost,atol=1e-13)
        self.assertLess(abs(rb['numerical_energy_defect_j']),1e-12)
    def test_multiple_events_reject(self):
        b=fixture();b['cloth']=Body([[0,0,.1],[0,0,.15]],[[0,0,-1]]*2,[1,1])
        with self.assertRaisesRegex(ValueError,'Multiple'):
            step_vertex_face_impact(b,.2,cloth_owner='cloth',surface_owner='skin',node_ids=[0,1],face_ids=[0],friction_static=.4,friction_kinetic=.3,max_energy_defect_j=1e-12)
        self.assertEqual(b['cloth'].time_s,0)
    def test_unsupported_remainder_rolls_back(self):
        b=fixture();b['cloth'].position_m[0,0]=.1
        with self.assertRaisesRegex(ValueError,'Deforming postimpact'):
            call(b)
        self.assertEqual(b['cloth'].time_s,0);np.testing.assert_array_equal(b['skin'].velocity_m_s,0)
        b=fixture();b['cloth']=Body([[0,0,.1],[0,0,-.01]],[[0,0,-1],[0,0,0]],[1,1])
        with self.assertRaises(ValueError):
            step_vertex_face_impact(b,.2,cloth_owner='cloth',surface_owner='skin',node_ids=[0,1],face_ids=[0],friction_static=.4,friction_kinetic=.3,max_energy_defect_j=1e-12)
        self.assertEqual(b['cloth'].time_s,0)
    def test_candidate_budget_and_no_mass_duplication(self):
        b=fixture()
        with self.assertRaises(ValueError):call(b,max_candidate_pairs=0)
        b['skin']=b['cloth']
        with self.assertRaises(ValueError):call(b)
    def test_energy_budget_rejection_is_atomic(self):
        b=fixture();b['cloth'].position_m[0,2]=1;b['cloth'].velocity_m_s[:]=0
        b['cloth'].elastic_forces=lambda x:(np.zeros_like(x),float(x[0,0]))
        with self.assertRaisesRegex(ValueError,'energy defect'):
            call(b,external_forces_n={'cloth':np.array([[1,0,0]])})
        self.assertEqual(b['cloth'].time_s,0);self.assertEqual(b['cloth'].position_m[0,0],0)
    def test_external_kick_accounting_without_contact(self):
        b=fixture();b['cloth'].position_m[0,2]=1;b['cloth'].velocity_m_s[:]=0
        r=call(b,external_forces_n={'cloth':np.array([[1,0,0]])})
        self.assertIsNone(r['event']);self.assertAlmostEqual(r['external_work_j'],.02)
        self.assertLess(abs(r['numerical_energy_defect_j']),1e-14)
if __name__=='__main__':unittest.main()
