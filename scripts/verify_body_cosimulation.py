"""Analytic two-reservoir fixture for software coupling verification only."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unittest
import importlib.util
import copy
import math


class Reservoir:
    def __init__(self, owner, value): self.owner, self.value, self.events = owner, value, []
    def checkpoint(self): return copy.deepcopy((self.value, self.events))
    def restore(self, state): self.value, self.events = copy.deepcopy(state)
    def advance(self, start, end, incoming, event_queue):
        from ihm.assembly.cosimulation import EngineStep
        from ihm.assembly.contracts import Exchange
        self.events.append(end)
        event_queue.append((self.owner, end))
        if self.owner == 'a':
            # Backward Euler leak, with its integrated amount sent to b.
            amount = self.value * (end-start) / (1 + end-start)
            self.value -= amount
            return EngineStep((Exchange('flow', start, end, 'mass', amount, 'a', 'b', 2*amount),), {'value': self.value})
        self.value += sum(x.amount_SI for x in incoming)
        return EngineStep((), {'value': self.value})


class Coupling(unittest.TestCase):
    def test_adapter_failure_restores_other_engine_and_queues(self):
        c,e=self.make()
        def fail(*args):
            e['b'].value=99.; e['b'].events.append('bad'); c.event_queue.append('bad')
            raise RuntimeError('adapter failure')
        e['b'].advance=fail
        with self.assertRaisesRegex(RuntimeError, 'adapter failure'): c.step(.2)
        self.assertEqual((e['a'].value,e['b'].value),(1.,0.))
        self.assertEqual(c.event_queue,[])
        self.assertEqual(e['a'].events + e['b'].events, [])

    def test_rate_and_wrong_interval_rejected_before_commit(self):
        from ihm.assembly.contracts import RateSample, Exchange
        from ihm.assembly.cosimulation import EngineStep
        for transfer in (RateSample('flow',0.,'mass',1.), Exchange('flow',0.,.1,'mass',1.,'a','b')):
            c,e=self.make()
            e['a'].advance=lambda *args: EngineStep((transfer,), {})
            with self.assertRaises(ValueError): c.step(.2)
            self.assertEqual(c.time_s,0.); self.assertEqual(c.exchanges,[])

    def test_registration_change_requires_revalidation(self):
        from ihm.assembly.contracts import StateSpec
        c,e=self.make()
        c.registry.register_state(StateSpec('extra','a','ea','kg','mass','lumped','time_domain','fixture'))
        with self.assertRaises(ValueError): c.step(.2)

    def test_nonlinear_exchange_iteration_and_adaptive_rejection(self):
        from ihm.assembly.cosimulation import EngineStep, StepRejected
        from ihm.assembly.contracts import Exchange
        c,e=self.make()
        # Each source's outgoing amount depends on the opposite incoming
        # amount. A large step deliberately rejects; smaller steps converge.
        def coupled(owner):
            def run(start,end,incoming,queue):
                if end-start > .11:
                    e[owner].events.append('rejected'); queue.append('rejected')
                    raise StepRejected('fixture step limit')
                amount=(end-start)*(.1 + .2*sum(x.amount_SI for x in incoming)**2)
                e[owner].value += sum(x.amount_SI for x in incoming)-amount
                e[owner].events.append(end); queue.append((owner,end))
                return EngineStep((Exchange('flow',start,end,'mass',amount,owner,'b' if owner=='a' else 'a'),),{'value':e[owner].value})
            return run
        for owner in e: e[owner].advance=coupled(owner)
        audits=c.advance(.2,.2)
        self.assertEqual(len(audits),2)
        self.assertTrue(all(a.iterations > 2 for a in audits))
        self.assertAlmostEqual(e['a'].value+e['b'].value,1.,places=9)
        self.assertEqual(len(c.event_queue),4)
        self.assertEqual(len(e['a'].events),2)

    def test_module_exists(self): self.assertIsNotNone(importlib.util.find_spec('ihm.assembly.cosimulation'))

    def make(self):
        from ihm.assembly.cosimulation import Coordinator
        from ihm.assembly.interfaces import ContractRegistry, Port
        from ihm.assembly.contracts import StateSpec, Interaction
        r = ContractRegistry('v1', {'ea', 'eb'})
        for owner in ('a', 'b'):
            r.register_state(StateSpec(owner, owner, 'e'+owner, 'kg', 'mass', 'lumped', 'time_domain', 'fixture'))
        r.register_interaction(Interaction('flow', 'ea', 'eb', 'flow', 'a_to_b', 'leak', 'p', 'fixture', ('analytic',)))
        for owner in ('a', 'b'): r.register_port(Port(owner, 'flow', owner, owner, 'mass', 'kg'))
        r.validate_scenario({'a', 'b'}, {'leak'}, {'p'}, {'analytic'})
        engines = {'a': Reservoir('a', 1.), 'b': Reservoir('b', 0.)}
        return Coordinator(engines, ('a', 'b'), r), engines

    def test_storage_and_work_cancel(self):
        c, e = self.make(); audit = c.step(.2)
        self.assertAlmostEqual(e['a'].value + e['b'].value, 1.)
        self.assertEqual(audit.balance_residuals['mass'], 0.)
        self.assertEqual(audit.work_residual_J, 0.)
        self.assertEqual(len(e['a'].events), 1)
        self.assertEqual(len(c.event_queue), 2)

    def test_rejection_restores_all_engines_and_queues(self):
        from ihm.assembly.cosimulation import StepRejected
        c, e = self.make(); c.event_queue.append(('initial', 0.))
        before = {k: v.checkpoint() for k,v in e.items()}
        with self.assertRaises(StepRejected): c.step(.2, max_iterations=1)
        self.assertEqual({k:v.checkpoint() for k,v in e.items()}, before)
        self.assertEqual(c.event_queue, [('initial', 0.)]); self.assertEqual(c.time_s, 0.)
        self.assertEqual(c.exchanges, [])

    def test_halving_converges_to_analytic_solution(self):
        errors=[]
        for dt in (.2, .1, .05):
            c,e=self.make(); c.advance(1., dt)
            errors.append(abs(e['a'].value - math.exp(-1)))
        self.assertLess(errors[1], .6*errors[0]); self.assertLess(errors[2], .6*errors[1])


if __name__ == '__main__': unittest.main()
