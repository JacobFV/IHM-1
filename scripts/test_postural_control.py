"""Run with python -m unittest discover -s scripts -p test_postural_control.py."""
import unittest
from ihm.native.postural_control import PosturalController,PosturalConfig

def state(length=.1,velocity=0):
    return {'muscles':{'soleus_r':{'fiber_length_m':length,'optimal_fiber_length_m':.1,'fiber_velocity_m_s':velocity}}}
class PosturalControlTests(unittest.TestCase):
    def test_native_stretch_and_velocity_increase_command(self):
        p=PosturalController(state(),PosturalConfig(.1,2,.1))
        self.assertAlmostEqual(p.commands(state())['soleus_r'],.1)
        self.assertAlmostEqual(p.commands(state(.11))['soleus_r'],.3)
        self.assertAlmostEqual(p.commands(state(.1,.1))['soleus_r'],.2)
        self.assertEqual(p.commands(state(.2,1))['soleus_r'],1)
        self.assertEqual(p.commands(state(.01,-1))['soleus_r'],0)
    def test_catalog_mismatch_is_not_silent_afference_loss(self):
        p=PosturalController(state())
        with self.assertRaisesRegex(ValueError,'catalog'):p.commands({'muscles':{}})
        with self.assertRaisesRegex(ValueError,'Unknown'):PosturalController(state(),muscle_baselines={'wrong':.1})
        with self.assertRaises(ValueError):PosturalConfig(length_gain=float('nan'))
if __name__=='__main__':unittest.main()
