"""Selected material sensor coverage is mandatory on physical state responses."""
from copy import deepcopy
import unittest
from ihm.native.mechanical_stream import NativeMechanicalStream

class Tests(unittest.TestCase):
    def stream(self):
        stream=NativeMechanicalStream.__new__(NativeMechanicalStream)
        stream.surface_sensor_identity={3:{'manifest_sha256':'a'*64,'quadrature_index':3,'triangle_index':8}}
        return stream
    def test_missing_duplicate_and_unrequested_rows_reject(self):
        stream=self.stream()
        for data in ({'kind':'advanced'}, {'kind':'observed','surface_foundation':{'sensor_points':[]}},
                     {'kind':'restored','surface_foundation':{'sensor_points':[{'quadrature_index':3}]*2}},
                     {'kind':'initialized','surface_foundation':{'sensor_points':[{'quadrature_index':4}]}}):
            with self.assertRaises(ValueError):stream._bind_surface_sensors(deepcopy(data))
    def test_zero_observation_is_present_and_receipt_is_independent(self):
        stream=self.stream()
        data={'kind':'advanced','surface_foundation':{'sensor_points':[{'quadrature_index':3,'indentation_m':0}]}}
        stream._bind_surface_sensors(data)
        row=data['surface_foundation']['sensor_points'][0]
        self.assertEqual(row['id'],'skin-contact-3')
        row['material_identity']['triangle_index']=99
        self.assertEqual(stream.surface_sensor_identity[3]['triangle_index'],8)
    def test_control_ack_and_static_diagnostics_do_not_claim_sensor_samples(self):
        stream=self.stream()
        for kind in ('checkpointed','dropped','static_pose_evaluated'):
            data={'kind':kind};stream._bind_surface_sensors(data)
            self.assertEqual(data,{'kind':kind})

if __name__=='__main__':unittest.main()
