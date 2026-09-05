"""Numerical/software contract checks; not biological validation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unittest
import importlib.util


class Contracts(unittest.TestCase):
    def test_frame_execution_semantics(self):
        from ihm.assembly import contracts
        self.assertTrue(hasattr(contracts, 'Frame'))
        frame=contracts.Frame(0.,'receipt',{}, {}, {}, {}, {}, {}, 'recorded_source_replay')
        self.assertEqual(frame.playback_kind,'recorded_source_replay')
        with self.assertRaises(ValueError): contracts.Frame(0.,'receipt',{}, {}, {}, {}, {}, {}, 'live')

    def test_material_identity_and_observation_receipt(self):
        from dataclasses import asdict, FrozenInstanceError
        import json
        from ihm.assembly.contracts import MaterialPoint, ObservationRequest, MaterializationReceipt
        point = MaterialPoint('v1', 'a', 'tet-1', [.1, .2, .3])
        self.assertEqual(hash(point), hash(MaterialPoint('v1', 'a', 'tet-1', (.1,.2,.3))))
        with self.assertRaises(FrozenInstanceError): point.element_id = 'display-vertex-4'
        with self.assertRaises(ValueError): MaterialPoint('v1', '', 'tet-1', (.1,))
        with self.assertRaises(ValueError): ObservationRequest('v1','body',('mass',),1.,float('nan'),.1,{'seconds':1.})
        request = ObservationRequest('v1','body',('mass',),.01,1.,.1,{'seconds':1.})
        receipt = MaterializationReceipt('r', {'data':'abc'}, {'code':'def'}, {'fixture':'1'}, asdict(request), 42, {'water':'a'}, ('flow',), {'temporal':None}, {'status':'unknown'}, 'analytic fixture')
        self.assertEqual(json.loads(json.dumps(asdict(receipt)))['uncertainty']['status'], 'unknown')

    def test_contract_module_exists(self):
        self.assertIsNotNone(importlib.util.find_spec('ihm.assembly.contracts'))

    def test_identity_dimension_ownership_and_scenario_gate(self):
        from ihm.assembly.contracts import MaterialPoint, StateSpec, Interaction
        from ihm.assembly.interfaces import ContractRegistry, Port
        r = ContractRegistry('v1', {'a', 'b'})
        p = MaterialPoint('v1', 'a', 'tet-1', (.2, .3, .1))
        r.validate_point(p)
        with self.assertRaises(ValueError): r.validate_point(MaterialPoint('v2', 'a', 'tet-1', (0.,)))
        state = StateSpec('water', 'engine-a', 'a', 'kg', 'mass', 'lumped', 'time_domain', 'reference')
        r.register_state(state)
        with self.assertRaises(ValueError): r.register_state(StateSpec('water', 'engine-b', 'a', 'kg', 'mass', 'lumped', 'time_domain', 'reference'))
        with self.assertRaises(ValueError): StateSpec('bad', 'engine-a', 'a', 's', 'mass', 'lumped', 'time_domain', 'reference')
        r.register_state(StateSpec('water-b', 'engine-b', 'b', 'kg', 'mass', 'lumped', 'time_domain', 'reference'))
        r.register_interaction(Interaction('flow', 'a', 'b', 'flow', 'a_to_b', 'law', 'params', 'reference', ('source',)))
        r.register_port(Port('pa', 'flow', 'engine-a', 'water', 'mass', 'kg'))
        r.register_port(Port('pb', 'flow', 'engine-b', 'water-b', 'mass', 'kg'))
        r.validate_scenario({'engine-a', 'engine-b'}, {'law'}, {'params'}, {'source'})
        with self.assertRaises(ValueError): r.validate_scenario({'engine-a'}, {'law'}, {'params'}, {'source'})
        with self.assertRaises(ValueError): r.validate_scenario({'engine-a', 'engine-b'}, set(), {'params'}, {'source'})
        with self.assertRaises(ValueError): r.require_ready({'engine-a','engine-b'})
        with self.assertRaises(ValueError): StateSpec('temperature','engine-a','a','K','temperature','lumped','time','reference')

    def test_both_interface_endpoints_need_ports(self):
        from ihm.assembly.contracts import StateSpec,Interaction
        from ihm.assembly.interfaces import ContractRegistry,Port
        r=ContractRegistry('v1',{'a','b'})
        for owner in ['x','y']:
            r.register_state(StateSpec(owner,owner,'a','kg','mass','lumped','time','fixture'))
        r.register_interaction(Interaction('flow','a','b','flow','a_to_b','law','p','fixture',('e',)))
        for owner in ['x','y']:r.register_port(Port(owner,'flow',owner,owner,'mass','kg'))
        with self.assertRaises(ValueError):r.validate_scenario({'x','y'},{'law'},{'p'},{'e'})

    def test_unknown_fit_data_is_parameter_ancestry(self):
        from ihm.assembly.evidence import EvidenceGraph,EvidenceNode,Uncertainty,ParameterCard
        g=EvidenceGraph()
        g.add(EvidenceNode('source','source',(),Uncertainty('quantified',1.,'kg')))
        g.add(EvidenceNode('fitdata','source',(),Uncertainty('unknown',reason='unmeasured error')))
        card=ParameterCard('p',1.,'kg','fixture',('source',),'unknown',('fitdata',),(),Uncertainty('quantified',0.,'kg'))
        with self.assertRaises(ValueError):g.validate_parameter(card)

    def test_exchange_uses_opposite_endpoints(self):
        from ihm.assembly.contracts import StateSpec,Interaction,Exchange
        from ihm.assembly.interfaces import ContractRegistry,Port
        r=ContractRegistry('v1',{'a','b'})
        r.register_interaction(Interaction('f','a','b','flow','a_to_b','law','p','fixture',('e',)))
        for owner,support in [('x','a'),('y','a'),('z','b')]:
            r.register_state(StateSpec(owner,owner,support,'kg','mass','lumped','time','fixture'))
            r.register_port(Port(owner,'f',owner,owner,'mass','kg'))
        r.validate_scenario({'x','y','z'},{'law'},{'p'},{'e'})
        with self.assertRaises(ValueError):r.validate_exchange(Exchange('f',0.,1.,'mass',1.,'x','y'))
        r.validate_exchange(Exchange('f',0.,1.,'mass',1.,'x','z'))

    def test_unknown_ancestry_and_cycle(self):
        from ihm.assembly.evidence import EvidenceGraph, EvidenceNode, Uncertainty, ParameterCard
        g = EvidenceGraph()
        g.add(EvidenceNode('raw', 'source', (), Uncertainty('unknown', reason='not measured')))
        g.add(EvidenceNode('fit', 'fit', ('raw',), Uncertainty('quantified', variance=1., unit='kg')))
        self.assertEqual(g.uncertainty('fit').status, 'unknown')
        with self.assertRaises(ValueError): g.add(EvidenceNode('loop', 'repair', ('loop',), Uncertainty('unknown', reason='unknown')))
        with self.assertRaises(ValueError): Uncertainty('unknown', variance=.3)
        card = ParameterCard('p', 1., 'kg', 'reference', ('fit',), 'unknown', (), (), Uncertainty('unknown', reason='no holdout'))
        g.validate_parameter(card)
        self.assertEqual(g.ancestors('fit'), {'raw'})
        with self.assertRaises(ValueError): g.add(EvidenceNode('missing', 'fit', ('absent',), Uncertainty('unknown', reason='missing')))
        self.assertNotIn('missing', g.nodes)


if __name__ == '__main__': unittest.main()
