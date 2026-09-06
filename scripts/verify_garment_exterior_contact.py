"""Tiny source-exterior and swept-contact rejection fixtures; no native launch."""
from pathlib import Path
import sys, unittest, tempfile, gzip, json, hashlib, copy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.garment_swept_vertex_face import swept_vertex_face
from ihm.assembly.garment_exterior_contact import StrictGarmentSurfaceContact, exterior_source

T=np.array([[-1.,-1,0],[1,-1,0],[0,1,0]])
class Checks(unittest.TestCase):
    def test_source_filter_and_new_tether_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);raw=gzip.compress(json.dumps(dict(positions=np.vstack((T,T+[0,0,1])).tolist(),indices=[[0,1,2],[3,4,5]])).encode())
            path='skin.json.gz';(root/path).write_bytes(raw);digest=hashlib.sha256(raw).hexdigest()
            evidence=dict(schema='engineered_skin_territories_v1',source_receipts=[dict(path=path,sha256=digest)],contact_eligible_triangle_ids=[0],surface_diagnostic=dict(face_component_ids=[0,1],selected_component_id=0,selection_evidence='inferred_largest_component_positive_signed_volume_integral',selected_area_m2=2.,physical_surface_exclusivity_validated=False,self_intersections_tested=False,signed_integral_is_closed_volume=False))
            out=root/'data/research/engineered_skin_territories/materialization.json';out.parent.mkdir(parents=True);out.write_text(json.dumps(evidence))
            previous=dict(garment_node=0,triangle_index=1,body_nodes=[3,4,5],barycentric=[.25,.25,.5],distance_m=.1)
            pairs=dict(garments=dict(shirt=dict(support_pairings=[previous],contact_samples=[])))
            x,tri,ids,new,support=exterior_source(root,dict(path=path,geometry_sha256=digest),pairs,dict(shirt=[[0,0,1.1]]))
            self.assertEqual(ids.tolist(),[0]);self.assertEqual(tri.tolist(),[[0,1,2]])
            row=new['garments']['shirt']['support_pairings'][0]
            self.assertEqual(row['previous_raw_source_pairing'],previous);self.assertEqual(row['triangle_index'],0)
            self.assertNotEqual(row['body_nodes'],previous['body_nodes']);self.assertIn('new_exterior',row['mapping_status'])
            self.assertEqual(pairs['garments']['shirt']['support_pairings'][0],previous)
            from unittest.mock import patch
            from types import SimpleNamespace
            from ihm.assembly.garment_feedback import GarmentFeedback
            from ihm.assembly.clothing import Cloth
            cloth=Cloth(T+[0,0,1.1],np.array([[0,1,2]]),areal_density_kg_m2=1,edge_stiffness_n_m=1)
            registration=SimpleNamespace(bodies=['body'],groups={'body':dict(bounds_min_m=[-2,-2,-2],bounds_max_m=[2,2,2])})
            with patch('ihm.assembly.garment_feedback.materialize_garments',return_value=({'shirt':cloth},{'skin':dict(path=path,geometry_sha256=digest)})),patch('ihm.assembly.garment_feedback.build_source_pairings',return_value=pairs):
                g=GarmentFeedback.from_root(root,registration)
            self.assertIs(g.surface_contact_class,StrictGarmentSurfaceContact)
            self.assertEqual(g.triangles.tolist(),[[0,1,2]])
            self.assertEqual(g.identity['contact_source_face_indices'],[0])
            self.assertEqual(g.attachments['shirt'][0].body_nodes,(0,1,2))

            with self.assertRaises(ValueError):exterior_source(root,dict(path=path,geometry_sha256='bad'),pairs,dict(shirt=[[0,0,1.1]]))
    def test_feedback_rejects_and_rolls_back(self):
        from verify_garment_feedback import fixture
        native,g=fixture();g.surface_contact_class=StrictGarmentSurfaceContact
        g.garments['patch'].velocity_m_s[:,2]=-60
        initial=native.snapshot();saved=g.checkpoint()
        with self.assertRaises(ValueError):g.advance(native,.001,[],{})
        self.assertEqual(native.snapshot(),initial);self.assertEqual(native.tokens,0)
        np.testing.assert_array_equal(g.garments['patch'].position_m,saved['garments']['patch']['position_m'])
        np.testing.assert_array_equal(g.garments['patch'].velocity_m_s,saved['garments']['patch']['velocity_m_s'])
    def test_endpoint_tunneling(self):
        e=swept_vertex_face([0,0,1],[0,0,-1],T,T)
        self.assertEqual(e['status'],'crossing');self.assertAlmostEqual(e['fraction'],.5)
        np.testing.assert_allclose(np.array(e['barycentric'])@T,e['point_m'],atol=1e-14)
    def test_moving_face_and_rigid_invariance(self):
        e=swept_vertex_face([0,0,0],[0,0,0],T+[0,0,-1],T+[0,0,1])
        self.assertEqual(e['status'],'crossing')
        q=np.array([[0,0,1],[1,0,0],[0,1,0]])
        r=swept_vertex_face(np.array([0,0,1])@q+3,np.array([0,0,-1])@q+3,T@q+3,T@q+3)
        self.assertAlmostEqual(r['fraction'],.5)
    def test_interior_tangency_is_unresolved(self):
        t=np.array([[0.,0,0],[1,0,-.5],[0,1,0]])
        u=np.array([[0.,0,0],[1,0,.5],[0,1,0]])
        e=swept_vertex_face([-.25,.25,-.125],[.75,.25,.125],t,u)
        self.assertEqual(e['status'],'unresolved')
        self.assertIn(e['reason'],['tangent_or_ill_conditioned_root','near_multiple_root'])
    def test_outside_edge_coplanar_and_degenerate(self):
        self.assertEqual(swept_vertex_face([3,0,1],[3,0,-1],T,T)['status'],'clear')
        self.assertEqual(swept_vertex_face([0,-1,1],[0,-1,-1],T,T)['status'],'unresolved')
        self.assertEqual(swept_vertex_face([0,0,0],[.1,0,0],T,T)['status'],'unresolved')
        self.assertEqual(swept_vertex_face([0,0,1],[0,0,-1],T,T[[0,2,1]])['status'],'unresolved')
    def test_strict_adapter_rejects_deep_crossing(self):
        c=StrictGarmentSurfaceContact(T,T,np.array([[0,1,2]]),1,friction_static=.4,friction_kinetic=.3);c.fraction=1
        with self.assertRaises(ValueError):c.resolve([[0,0,-.1]],[[0,0,-.2]],[1],1)
    def test_nearest_edge_cannot_be_accepted(self):
        c=StrictGarmentSurfaceContact(T,T,np.array([[0,1,2]]),1,friction_static=.4,friction_kinetic=.3);c.fraction=1
        with self.assertRaisesRegex(ValueError,'nearest-edge'):
            c.resolve([[0,-1.001,-.001]],[[0,0,0]],[1],1)
    def test_strict_adapter_clear_and_budget(self):
        c=StrictGarmentSurfaceContact(T,T,np.array([[0,1,2]]),1,friction_static=.4,friction_kinetic=.3);c.fraction=1
        r=c.resolve([[0,0,.1]],[[0,0,0]],[1],1);self.assertEqual(r['contact_count'],0)
        c.max_pairs=1
        with self.assertRaises(MemoryError):c.resolve([[0,0,.001],[.1,0,.001]],[[0,0,0]]*2,[1,1],1)
if __name__=='__main__':unittest.main()
