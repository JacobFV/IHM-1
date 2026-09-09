import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.balance_observation import convex_hull_xz,signed_polygon_margin,observe_balance
class Tests(unittest.TestCase):
    def test_hull_margin(self):
        hull=convex_hull_xz([[1,1],[0,0],[0,1],[1,0],[.5,.5],[0,0]])
        self.assertEqual(len(hull),4);self.assertAlmostEqual(signed_polygon_margin([.5,.5],hull),.5)
        self.assertAlmostEqual(signed_polygon_margin([2,2],hull),-2**.5)
        self.assertEqual(signed_polygon_margin([0,.5],hull),0.)
        self.assertIsNone(signed_polygon_margin([0,0],convex_hull_xz([[0,0],[1,0],[2,0]])))
    def test_real_native_snapshot(self):
        root=Path(__file__).resolve().parents[1]
        state=json.loads((root/'data/derived/mechanics/resting_stance86/initial_snapshot.json').read_text())
        result=observe_balance(state);json.dumps(result,allow_nan=False)
        self.assertEqual(result['muscle_count'],86);self.assertGreater(result['active_foot_contact_count'],2)
        self.assertIsNotNone(result['com_support_margin_m']);self.assertLess(result['muscle_activation_excitation_lag_max_abs'],1e-8)
        for contact in state['contacts']:contact['force_n']=[0.,1.,0.]
        result=observe_balance(state);self.assertIsNone(result['com_support_margin_m']);self.assertEqual(result['active_foot_contact_count'],0)
if __name__=='__main__':unittest.main()
