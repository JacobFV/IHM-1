"""Transactional partitioned fixed-point coupling of checkpointable adapters.

Each adapter emits its outgoing integrated transfers exactly once and accounts
for its own debit. Incoming transfers are applied by recipients. Iterations
restore the interval-start state, including all adapter-owned event queues.
This coordinator verifies ledger consistency and iterative residuals; physical
storage and constitutive validity require independent adapter audits.
"""
from dataclasses import dataclass, field
from typing import Mapping
from copy import deepcopy
from .contracts import Exchange, finite


@dataclass(frozen=True)
class EngineStep:
    exchanges: tuple[Exchange, ...] = ()
    observations: Mapping[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class StepAudit:
    interval_start_s: float
    interval_end_s: float
    iterations: int
    scaled_iteration_residual: float
    balance_residuals: Mapping[str, float]
    work_residual_J: float
    error_estimates: Mapping[str, object]
    playback_kind: str = 'cosimulation'


class StepRejected(RuntimeError):
    pass


class Coordinator:
    def __init__(self, engines, order, registry, *, time_s=0.):
        finite(time_s)
        self.engines, self.order, self.registry = dict(engines), tuple(order), registry
        if len(self.order) != len(set(self.order)) or set(self.order) != set(self.engines):
            raise ValueError('Explicit coupling order must contain each engine once')
        registry.require_ready(self.engines)
        self.time_s, self.event_queue, self.exchanges, self.audits = time_s, [], [], []

    def checkpoint(self):
        return deepcopy((self.time_s, {k:e.checkpoint() for k,e in self.engines.items()}, self.event_queue, self.exchanges, self.audits))

    def restore(self, checkpoint):
        time_s, states, queue, exchanges, audits = deepcopy(checkpoint)
        # Attempt every restore even if one adapter fails; failure is fatal.
        failures = []
        for owner in self.order:
            try: self.engines[owner].restore(states[owner])
            except Exception as exc: failures.append((owner, exc))
        self.time_s, self.event_queue, self.exchanges, self.audits = time_s, queue, exchanges, audits
        if failures: raise RuntimeError('Adapter checkpoint restoration failed: ' + ', '.join(k for k,_ in failures)) from failures[0][1]

    @staticmethod
    def _values(exchanges, observations):
        values = dict(observations)
        for x in exchanges:
            key = ('exchange', x.interface_id, x.quantity, x.source_owner, x.target_owner)
            if key in values: raise ValueError('Aggregate each directed interface quantity into one interval exchange')
            values[key] = x.amount_SI
            values[key + ('work',)] = x.work_J_or_none or 0.
        return values

    def step(self, end_time_s, *, atol=1e-10, rtol=1e-8, max_iterations=30):
        finite(end_time_s); finite(atol, positive=True); finite(rtol)
        if end_time_s <= self.time_s or rtol < 0: raise ValueError('Invalid coupling interval or tolerance')
        if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1: raise ValueError('Positive iteration count required')
        self.registry.require_ready(self.engines)
        initial = self.checkpoint()
        start = self.time_s
        previous, previous_values = (), None
        try:
            for iteration in range(1, max_iterations + 1):
                self.restore(initial)
                current, observations = [], {}
                for owner in self.order:
                    incoming = tuple(x for x in previous if x.target_owner == owner)
                    result = self.engines[owner].advance(start, end_time_s, incoming, self.event_queue)
                    if not isinstance(result, EngineStep): raise ValueError('Adapter must return EngineStep')
                    for x in result.exchanges:
                        self.registry.validate_exchange(x)
                        if x.source_owner != owner: raise ValueError('Only source owner can emit a transfer')
                        if x.interval_start_s != start or x.interval_end_s != end_time_s: raise ValueError('Transfer interval mismatch')
                        current.append(x)
                    for key, value in result.observations.items():
                        finite(value); observations[('observation', owner, key)] = value
                values = self._values(current, observations)
                residual = float('inf')
                if previous_values is not None and values.keys() == previous_values.keys():
                    residual = max((abs(value-previous_values[key])/(atol+rtol*max(abs(value), abs(previous_values[key]))) for key,value in values.items()), default=0.)
                if residual <= 1.:
                    balances, work = {}, 0.
                    transfer_mismatch = {}
                    for x in previous:
                        transfer_mismatch[x.quantity] = transfer_mismatch.get(x.quantity, 0.) + x.amount_SI
                    for x in current:
                        transfer_mismatch[x.quantity] = transfer_mismatch.get(x.quantity, 0.) - x.amount_SI
                        for _, amount, signed_work in x.signed_entries():
                            balances[x.quantity] = balances.get(x.quantity, 0.) + amount
                            work += signed_work
                    audit = StepAudit(start, end_time_s, iteration, residual, balances, work,
                        {'coupling_iteration_scaled_residual': residual,
                         'accepted_incoming_minus_outgoing_SI': transfer_mismatch,
                         'time_discretization': None, 'biological_prediction': None,
                         'note': 'Ledger cancellation is not an independent physical storage audit; unknown errors remain unquantified.'})
                    self.time_s = end_time_s
                    self.exchanges.extend(current); self.audits.append(audit)
                    return audit
                previous, previous_values = tuple(current), values
            raise StepRejected('Coupling did not converge within iteration limit')
        except Exception:
            self.restore(initial)
            raise

    def advance(self, end_time_s, step_s, *, min_step_s=1e-8, **step_options):
        """Advance accepted intervals; halve rejected steps, retaining prior accepts."""
        finite(end_time_s); finite(step_s, positive=True); finite(min_step_s, positive=True)
        if end_time_s < self.time_s: raise ValueError('Cannot integrate backwards')
        audits = []
        dt = step_s
        while self.time_s < end_time_s:
            target = min(end_time_s, self.time_s + dt)
            if target <= self.time_s: raise StepRejected('Step below clock precision')
            try: audits.append(self.step(target, **step_options))
            except StepRejected:
                dt *= .5
                if dt < min_step_s: raise
        return audits
