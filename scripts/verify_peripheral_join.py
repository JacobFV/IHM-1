"""Join completeness, geometry scope and fibre timing regression checks."""
import ast
from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.assembly.peripheral import BodyPeripheral
from ihm.assembly.reflexes import ReflexParameters

class JoinVerification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/'data/derived/canonical/peripheral.json').read_text())

    def test_complete_bilateral_join_and_no_silent_unknowns(self):
        tree=ast.parse((ROOT.parent/'IBM-1/ibm/topologies/nerve.py').read_text())
        expected=next({k.value for k in n.value.keys} for n in tree.body
                      if isinstance(n,ast.AnnAssign) and n.target.id=='TRUNK_COMPOSITION')
        aliases=self.data['route_contract']['aliases']
        for side in ('left','right'):
            names={n['id'].removeprefix(f'peripheral-nerve-{side}-') for n in self.data['nerves'] if n['side']==side}
            joined={aliases.get(n,n) for n in names}
            self.assertEqual(joined-expected,{'sciatic_fibular'})
            self.assertFalse(expected-joined)
        self.assertEqual(len({n['id'] for n in self.data['nerves']}),144)

    def test_lengths_and_relays(self):
        relays={r['id']:r for r in self.data['relays']}
        brain={n['id'] for n in json.loads((ROOT/'data/derived/canonical/brain.json').read_text())['nodes']}
        for n in self.data['nerves']:
            self.assertIn(n['relay_id'],relays)
            self.assertTrue(math.isfinite(n['path_length_m']) and n['path_length_m']>0)
            self.assertFalse(n['measured_axon_geometry'])
            self.assertEqual(n['geometry_kind'],'schematic_route')
            self.assertEqual(n['path_length_scope'],'representative_endpoint_to_relay')
            self.assertTrue(n['limitations'])
            if 'points_m' in n:
                self.assertEqual(n['points_m'][-1],relays[n['relay_id']]['position_m'])
                self.assertAlmostEqual(n['path_length_m'],sum(math.dist(a,b) for a,b in zip(n['points_m'],n['points_m'][1:])))
            for target in n.get('downstream_brain_target_ids',[]):self.assertIn(target,brain)
        for name in ('optic','cochlear','vestibular','olfactory'):
            n=next(n for n in self.data['nerves'] if n['id']==f'peripheral-nerve-left-{name}')
            self.assertNotIn('cranial',n['relay_id'])
            self.assertTrue(all('postcentral' not in t for t in n['downstream_brain_target_ids']))

    def test_runtime_uses_class_delays_not_legacy_scalar(self):
        data=deepcopy(self.data)
        b=next(b for b in data['muscle_bindings'] if b['muscle_id']=='body-muscle-opensim-tibant_l')
        b['motor_delay_s']=b['afferent_delay_s']=999.
        body=BodyPeripheral(data)
        body.step(.001,brain_state={'motor_commands':{b['muscle_id']:.5}},
                  mechanical_state={'muscle_forces_n':{b['muscle_id']:100.}})
        motor=next(e for e in body.events if e[2]=='motor' and e[3]==b['muscle_id'])
        self.assertAlmostEqual(motor[0],b['delays_s']['alpha']+b['central_delay_s'])
        afferent=next(e for e in body.events if e[2]=='afferent' and e[3]=='proprio:'+b['muscle_id'])
        self.assertAlmostEqual(afferent[0],.001+b['delays_s']['ia']+b['central_delay_s'])
        self.assertGreater(b['delays_s']['gamma'],b['delays_s']['alpha'])
        p=ReflexParameters.from_binding(b)
        self.assertAlmostEqual(p.loop_delay_s,b['delays_s']['ia']+b['delays_s']['alpha']+.001)
        for c,v in b['delays_s'].items():
            self.assertAlmostEqual(v,b['path_length_m']/data['fibre_velocity_m_s'][c]['typical'])

if __name__=='__main__':unittest.main()
