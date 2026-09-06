"""Native material observations retain identity and work under rigid frame mapping."""
import unittest
import numpy as np
from ihm.assembly.articulated import CanonicalRegistration

class Tests(unittest.TestCase):
    def test_contact_frame_covariance_and_identity(self):
        r=CanonicalRegistration.__new__(CanonicalRegistration)
        r.basis=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
        r.global_map=np.eye(4);r.global_map[:3,:3]=r.basis;r.global_map[:3,3]=[1,2,3]
        point={'id':'skin-contact-3','point_source_m':[.1,.2,.3],'normal_source':[-1,0,0],
               'force_n':[2,3,4],'indentation_m':.001,'contact_area_m2':.0001,
               'material_identity':{'manifest_sha256':'a'*64,'quadrature_index':3,'triangle_index':8},
               'indentation_basis':'native nonlinear confined-layer modeled compression'}
        native={'surface_foundation':{'sensor_points':[point]}}
        out=r.cutaneous_contacts(native)[0]
        np.testing.assert_allclose(out['point_m'],[.8,2.1,3.3])
        np.testing.assert_allclose(out['force_n'],[-3,2,4])
        np.testing.assert_allclose(out['normal'],[0,-1,0])
        self.assertEqual(out['material_identity'],point['material_identity'])
        self.assertAlmostEqual(-np.dot(out['normal'],out['force_n']),2)
        velocity=np.array([.3,.2,.1]);self.assertAlmostEqual(np.dot(out['force_n'],r.basis@velocity),np.dot(point['force_n'],velocity))
        out['material_identity']['triangle_index']=9
        self.assertEqual(point['material_identity']['triangle_index'],8)
        self.assertIsNone(r.cutaneous_contacts({}))
        self.assertEqual(r.cutaneous_contacts({'surface_foundation':{'sensor_points':[]}}),[])

    def test_factory_binding_requires_exact_registered_materials(self):
        from pathlib import Path
        from ihm.assembly.embodied import bind_cutaneous
        contact={'id':'skin-contact-3','point_m':[0,0,0],'normal':[0,0,-1],
            'force_n':[0,0,1],'contact_area_m2':.001,'indentation_m':.0001,
            'material_identity':{'manifest_sha256':'a'*64,'quadrature_index':3,'triangle_index':8},
            'indentation_basis':'native nonlinear confined-layer modeled compression',
            'area_basis':'projected reference quadrature area'}
        config={'regions':{'skin-contact-3':'brain-rh-postcentral'},
                'recruitment_hz_per_response':.1,'reference_temperature_C':33}
        root=Path(__file__).resolve().parents[1]
        model=bind_cutaneous(root,[contact],config)
        frame=model.step(.02,{'time_s':0,'contacts':[contact]})
        self.assertEqual(frame['sites'][0]['indentation_um'],100.)
        with self.assertRaises(ValueError):bind_cutaneous(root,[],config)
        with self.assertRaises(ValueError):bind_cutaneous(root,[contact,contact],config)

if __name__=='__main__':unittest.main()
