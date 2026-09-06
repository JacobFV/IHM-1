"""Actual-source atomic rollback and geometric pressure work checks."""
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from ihm.assembly.thoracic_mechanism import ThoracicMechanism
from ihm.assembly.thoracic_cervical_metric import CoupledThoracicMetric
from ihm.assembly.thoracic_composite_owner import ThoracicCompositeOwner
ROOT=Path(__file__).resolve().parents[1]

class OwnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.metric=CoupledThoracicMetric(ThoracicMechanism(ROOT/'data/research/thoracic_mechanism/v2/manifest.json'),ROOT/'data/research/thoracic_mechanism/native_composition_v1/plan.json')
    def assert_state_equal(self,a,b):
        for key in a:np.testing.assert_array_equal(a[key],b[key])
    def test_actual_pressure_step_and_exact_rollback(self):
        owner=ThoracicCompositeOwner(self.metric);start=owner.state;start['velocity'][31]=.005
        owner=ThoracicCompositeOwner(self.metric,start);checkpoint=owner.checkpoint();v0=owner.pressure_port(2.)['geometric_volume_m3']
        result=owner.advance(.001,2.);v1=owner.pressure_port(2.)['geometric_volume_m3']
        self.assertAlmostEqual(result['mechanical_pressure_work_J'],2*(v1-v0),places=15)
        self.assertGreater(result['state']['time_s'],0);self.assertFalse(result['native_state_advanced'])
        self.assert_state_equal(owner.restore(checkpoint),start)
        self.assertGreater(owner.revision,result['revision']);owner.release_checkpoint(checkpoint)
        with self.assertRaises(ValueError):owner.restore(checkpoint)
    def test_state_copy_failed_step_and_foreign_checkpoint(self):
        owner=ThoracicCompositeOwner(self.metric);before=owner.state;external=owner.state;external['q'][0]=.001;self.assert_state_equal(owner.state,before)
        other=ThoracicCompositeOwner(self.metric)
        with self.assertRaises(ValueError):owner.restore(other.checkpoint())
        with patch.object(self.metric,'evaluate',side_effect=ValueError('forced mechanical failure')):
            with self.assertRaises(ValueError):owner.advance(.001,2.)
        first={'velocity_derivative':np.zeros(38),'pressure_port':{'geometric_volume_m3':.006}}
        with patch.object(self.metric,'evaluate',side_effect=[first,ValueError('midpoint failure')]):
            with self.assertRaises(ValueError):owner.advance(.001,2.)
        self.assert_state_equal(owner.state,before);self.assertEqual(owner.revision,0)
        tokens=[owner.checkpoint() for _ in range(8)]
        with self.assertRaises(ValueError):owner.checkpoint()
        for token in tokens:owner.release_checkpoint(token)

if __name__=='__main__':unittest.main()
