"""Work-conjugate mechanical/native respiratory boundary; no native process."""
import unittest
import numpy as np
from ihm.assembly.respiratory_feedback import RespiratoryLoadPort


class Tests(unittest.TestCase):
    def port(self):
        return RespiratoryLoadPort({'volume_jacobian_m2':[.1,.2,.3],
            'stiffness_n_m':np.diag([100.,200.,300.]).tolist(),
            'bindings':[{'entity_id':'sternum','translation_basis':np.eye(3).tolist()},
                        {'entity_id':'diaphragm','translation_basis':np.diag([0,0,-1]).tolist()}]})

    def test_external_work_and_pressure_are_conjugate(self):
        port=self.port();force=np.array([1.,-2.,3.]);dvolume=.00012
        load=port.evaluate([{'id':'sternum','force_n':force.tolist()}])
        actual_work=force@(np.array(load['entity_displacement_per_volume_m_per_m3']['sternum'])*dvolume)
        self.assertAlmostEqual(actual_work,-load['external_pressure_pa']*dvolume,places=14)
        self.assertAlmostEqual(np.dot([.1,.2,.3],load['mode_displacement_per_volume_m_per_m3']),1.)

    def test_released_force_and_repeated_owners(self):
        port=self.port();self.assertEqual(port.evaluate([])['external_pressure_pa'],0)
        whole=port.evaluate([{'id':'sternum','force_n':[0,0,2]}])
        pieces=port.evaluate([{'id':'sternum','force_n':[0,0,1]},{'id':'sternum','force_n':[0,0,1]}])
        self.assertEqual(whole['external_pressure_pa'],pieces['external_pressure_pa'])

    def test_invalid_inputs_do_not_create_pressure(self):
        port=self.port()
        for forces in [[{'id':'unknown','force_n':[1,0,0]}],[{'id':'sternum','force_n':[float('nan'),0,0]}]]:
            with self.assertRaises(ValueError):port.evaluate(forces)
        self.assertEqual(port.evaluate([])['external_pressure_pa'],0)


if __name__=='__main__':unittest.main()
