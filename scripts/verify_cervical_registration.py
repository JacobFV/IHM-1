"""Bounded source-only registration regressions."""
import unittest
import copy,gzip,hashlib,json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from ihm.assembly.cervical_registration import proper_fit, scalar_function, reference_joint,validate_recipe

class RegistrationTests(unittest.TestCase):
    def test_proper_fit_recovers_rigid_frame(self):
        x=np.array([[0,0,0],[1,0,0],[0,2,0],[0,0,3.]])
        r=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        t,stats=proper_fit(x,x@r.T+[4,5,6])
        np.testing.assert_allclose(t[:3,:3],r,atol=1e-14)
        np.testing.assert_allclose(t[:3,3],[4,5,6],atol=1e-14)
        self.assertLess(stats['rms_m'],1e-14)
    def test_collinear_or_nonfinite_correspondence_rejected(self):
        x=np.array([[0,0,0],[1,0,0],[2,0,0.]])
        with self.assertRaises(ValueError):proper_fit(x,x)
        x[0,0]=np.nan
        with self.assertRaises(ValueError):proper_fit(x,x)
    def test_reflection_never_returned(self):
        x=np.array([[0,0,0],[1,0,0],[0,2,0],[0,0,3.]])
        t,s=proper_fit(x,x*[-1,1,1])
        self.assertAlmostEqual(np.linalg.det(t[:3,:3]),1)
        self.assertGreater(s['rms_m'],.1)
    def test_function_scope_and_neutral_frame(self):
        f=ET.fromstring('<SimmSpline><x>-2 2</x><y>-1 1</y></SimmSpline>')
        self.assertEqual(scalar_function(f,0),0)
        with self.assertRaises(ValueError):scalar_function(ET.fromstring('<SimmSpline><x>0 1 2</x><y>0 1 2</y></SimmSpline>'),0)
        j=ET.fromstring('<WeldJoint><parent_body>spine</parent_body><location_in_parent>1 2 3</location_in_parent><location>0 1 0</location><orientation_in_parent>0 0 0</orientation_in_parent><orientation>0 0 0</orientation></WeldJoint>')
        np.testing.assert_allclose(reference_joint(j,{})[:3,3],[1,1,3])
        j.tag='PinJoint'
        with self.assertRaises(ValueError):reference_joint(j,{})


    def test_retained_recipe_identity_ledger_and_attachment_mapping(self):
        root=Path(__file__).resolve().parents[1]
        directory=root/'data/research/cervical_registration/v2'
        recipe=json.loads((directory/'recipe.json').read_bytes())
        for entry in recipe['inputs']+[r['source'] for r in recipe['donor_geometry'].values()]:
            raw=gzip.decompress((directory/entry['retained_copy']).read_bytes())
            self.assertEqual(hashlib.sha256(raw).hexdigest(),entry['sha256'])
            self.assertEqual(len(raw),entry['bytes'])
        entry=recipe['target_model_identity']
        target=gzip.decompress((root/'data/research/cervical_inertia/v2'/entry['retained_copy']).read_bytes())
        self.assertEqual(validate_recipe(recipe,target)['attachments'],198)
        with self.assertRaises(ValueError):validate_recipe(recipe,target+b' ')
        corrupt=copy.deepcopy(recipe)
        corrupt['bodies']['cerv1']['mass_properties_in_new_body_frame']['mass_kg']+=.1
        with self.assertRaises(ValueError):validate_recipe(corrupt,target)
        corrupt=copy.deepcopy(recipe);corrupt['muscles'][0]['path_points'][0]['target_body']='nonexistent'
        with self.assertRaises(ValueError):validate_recipe(corrupt,target)
        self.assertEqual(recipe['unsupported_geometry'][0]['missing_exact_geometry'],'rotatedcerv7.vtp')
        for name,body in recipe['bodies'].items():
            transform=np.array(body['body_to_target_torso_reference'])
            for muscle in recipe['muscles']:
                for point in muscle['path_points']:
                    if point['source_body']==name:
                        np.testing.assert_array_equal(point['target_location_m'],point['source_location_m'])
            np.testing.assert_allclose(transform[:3,:3].T@transform[:3,:3],np.eye(3),atol=1e-12)

if __name__=='__main__':unittest.main()
