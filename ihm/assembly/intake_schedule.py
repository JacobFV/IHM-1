"""Deterministic Meal command scheduling; owns no physiological quantities.

Single-controller API. Claim before sending, then record the native outcome.
A claimed command is never returned again, including after uncertain delivery.
"""
from dataclasses import asdict, dataclass
import re

from ihm.native import number
from ihm.native.session import Meal

STEP_S = .02
MAX_HORIZON_S = 86400
MAX_EVENTS = 10000


def _ticks(seconds, maximum, name):
    number(seconds, 0, maximum, name)
    ticks = round(seconds * 50)
    if abs(ticks - seconds * 50) > 1e-7:
        raise ValueError(f'{name} must align to the native 0.02 s step')
    return ticks


@dataclass(frozen=True)
class IntakeEvent:
    event_id: str
    time_s: float
    meal: Meal

    def __post_init__(self):
        if not isinstance(self.event_id, str) or not re.fullmatch(r'[A-Za-z0-9_-]{1,64}', self.event_id):
            raise ValueError('event_id must be 1–64 letters, digits, underscores or hyphens')
        _ticks(self.time_s, MAX_HORIZON_S, 'time_s')
        if not isinstance(self.meal, Meal):
            raise ValueError('Expected native Meal')
        # Revalidate even if an externally constructed instance bypassed __init__.
        Meal(**asdict(self.meal))

    @property
    def due_tick(self):
        return round(self.time_s * 50)


@dataclass(frozen=True)
class IntakeReceipt:
    event_id: str
    scheduled_tick: int
    issued_tick: int
    outcome: str
    native_sequence: int | None = None
    reason: str | None = None


class IntakeSchedule:
    """One-shot reservations with caller-recorded accepted/uncertain outcomes."""

    def __init__(self, events=(), *, horizon_s=MAX_HORIZON_S):
        self.horizon_ticks = _ticks(horizon_s, MAX_HORIZON_S, 'horizon_s')
        if self.horizon_ticks < 1:
            raise ValueError('horizon_s must be at least 0.02')
        validated = []
        ids = set()
        for event in events:
            if len(validated) >= MAX_EVENTS:
                raise ValueError('Too many intake events')
            if not isinstance(event, IntakeEvent):
                raise ValueError('Expected IntakeEvent')
            IntakeEvent(event.event_id, event.time_s, event.meal)
            if event.event_id in ids:
                raise ValueError('Duplicate intake event identity')
            if event.due_tick > self.horizon_ticks:
                raise ValueError('Intake event exceeds schedule horizon')
            ids.add(event.event_id)
            validated.append(event)
        self.events = tuple(sorted(validated, key=lambda e: (e.due_tick, e.event_id)))
        self._cursor = 0
        self._last_tick = 0
        self._pending = None
        self._receipts = []
        self._uncertain = False

    def add_events(self, events, current_tick):
        """Atomically append future/current-clock events while retaining all IDs."""
        if type(current_tick) is not int or not self._last_tick <= current_tick <= self.horizon_ticks:
            raise ValueError('current_tick must be a bounded monotonic integer tick')
        if self._pending is not None or self._uncertain:
            raise RuntimeError('Cannot modify unresolved or uncertain intake schedule')
        candidate = IntakeSchedule(events, horizon_s=self.horizon_ticks * STEP_S)
        existing = {e.event_id for e in self.events}
        if len(self.events) + len(candidate.events) > MAX_EVENTS:
            raise ValueError('Too many intake events')
        if any(e.event_id in existing for e in candidate.events):
            raise ValueError('Duplicate intake event identity')
        if any(e.due_tick < current_tick for e in candidate.events):
            raise ValueError('Cannot add intake events in the past')
        self.events = self.events[:self._cursor] + tuple(sorted(
            self.events[self._cursor:] + candidate.events, key=lambda e: (e.due_tick, e.event_id)))
        self._last_tick = current_tick

    def snapshot(self):
        """JSON-compatible observation, never a resumable native checkpoint."""
        receipts = {r.event_id: r for r in self._receipts}
        result = []
        for event in sorted(self.events, key=lambda e: (e.due_tick, e.event_id)):
            entry = asdict(event)
            entry['state'] = 'queued'
            if event.event_id in receipts:
                receipt = asdict(receipts[event.event_id])
                entry.update(receipt)
                entry['state'] = receipt['outcome']
            elif self._pending is not None and self._pending[0].event_id == event.event_id:
                entry.update(state='issued', issued_tick=self._pending[1], scheduled_tick=event.due_tick)
            result.append(entry)
        return {'events': result}

    @property
    def receipts(self):
        return tuple(self._receipts)

    def due(self, body_tick):
        """Reserve at most one due event before caller sends its native command.

        Caller must first establish the native process has no pending meal.
        Unreceipted or uncertain commands block all subsequent reservations.
        """
        if type(body_tick) is not int or not 0 <= body_tick <= self.horizon_ticks:
            raise ValueError('body_tick must be a bounded integer native tick')
        if body_tick < self._last_tick:
            raise ValueError('Body clock cannot go backwards')
        if self._pending is not None or self._uncertain:
            raise RuntimeError('Unresolved or uncertain intake command; no replay permitted')
        self._last_tick = body_tick
        if self._cursor == len(self.events) or self.events[self._cursor].due_tick > body_tick:
            return ()
        event = self.events[self._cursor]
        self._cursor += 1
        self._pending = (event, body_tick)
        return (event,)

    def _claimed(self, event_id):
        if self._pending is None or self._pending[0].event_id != event_id:
            raise ValueError('Receipt must match the outstanding event')
        return self._pending

    def record_accepted(self, event_id, native_receipt):
        """Record a successful NativeSession.meal return at the issuing body tick."""
        event, tick = self._claimed(event_id)
        if not isinstance(native_receipt, dict):
            raise ValueError('Expected native acknowledgment mapping')
        sequence = native_receipt.get('sequence')
        if (native_receipt.get('status') != 'ok' or type(sequence) is not int or sequence < 1
                or native_receipt.get('pending_meal') is not True):
            raise ValueError('Expected successful native meal acknowledgment')
        elapsed = native_receipt.get('elapsed_s')
        number(elapsed, 0, self.horizon_ticks * STEP_S, 'native elapsed_s')
        if abs(elapsed - tick * STEP_S) > 1e-8:
            raise ValueError('Native acknowledgment differs from issuing body tick')
        prior = [r.native_sequence for r in self._receipts if r.native_sequence is not None]
        if prior and sequence <= prior[-1]:
            raise ValueError('Native acknowledgment sequence must increase')
        receipt = IntakeReceipt(event_id, event.due_tick, tick, 'accepted', sequence)
        self._receipts.append(receipt)
        self._pending = None
        return receipt

    def record_uncertain(self, event_id, reason):
        """Terminal outcome: neither retry this event nor issue further commands."""
        event, tick = self._claimed(event_id)
        if not isinstance(reason, str) or not reason.strip() or len(reason) > 1024:
            raise ValueError('Uncertainty reason must contain 1–1024 characters')
        receipt = IntakeReceipt(event_id, event.due_tick, tick, 'uncertain', reason=reason)
        self._receipts.append(receipt)
        self._pending = None
        self._uncertain = True
        return receipt

    def checkpoint(self):
        """Configuration-only checkpoint, permitted only before any issuance."""
        if self._cursor:
            raise RuntimeError('Cannot checkpoint issued intake state or native physiology')
        return dict(schema='ihm.intake-schedule.unissued.v1', horizon_s=self.horizon_ticks * STEP_S,
                    body_tick=self._last_tick, events=[asdict(e) for e in self.events])

    @classmethod
    def from_checkpoint(cls, checkpoint):
        if (not isinstance(checkpoint, dict)
                or set(checkpoint) != {'schema', 'horizon_s', 'body_tick', 'events'}
                or checkpoint['schema'] != 'ihm.intake-schedule.unissued.v1'
                or not isinstance(checkpoint['events'], list)
                or len(checkpoint['events']) > MAX_EVENTS):
            raise ValueError('Expected unissued schedule checkpoint')
        events = []
        try:
            for entry in checkpoint['events']:
                if not isinstance(entry, dict) or set(entry) != {'event_id', 'time_s', 'meal'}:
                    raise ValueError('Malformed checkpoint event')
                events.append(IntakeEvent(entry['event_id'], entry['time_s'], Meal(**entry['meal'])))
        except TypeError as exc:
            raise ValueError('Malformed checkpoint Meal') from exc
        result = cls(events, horizon_s=checkpoint['horizon_s'])
        tick = checkpoint['body_tick']
        if type(tick) is not int or not 0 <= tick <= result.horizon_ticks:
            raise ValueError('Invalid checkpoint body_tick')
        result._last_tick = tick
        return result
