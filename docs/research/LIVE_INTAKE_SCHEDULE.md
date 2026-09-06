# Live intake schedule contract

`ihm.assembly.intake_schedule` schedules existing `ihm.native.session.Meal`
actions. It does not compute digestion, nutrient uptake, energy, fluid stores,
or any replacement physiology. Native BioGears remains the owner of those
quantities. This module starts no processes and changes no body runtime state.

## Controller integration

Create `IntakeEvent(event_id, time_s, Meal(...))` entries and pass them to
`IntakeSchedule(events, horizon_s=...)`. Meal quantities and names use the
existing native `Meal` validation without alternative units: carbohydrate,
protein, fat and sodium in grams, calcium in milligrams, water in milliliters.
Every amount is finite and between zero and 10,000; at least one is positive.
A water-only event uses `Meal(name='drink', water_ml=250)`.

`time_s` is elapsed body time from the start of the continuing native session,
not its physiological `time_s` origin, wall time, or a fresh origin after every
control batch. Times must align to 20 ms. The horizon is 0.02–86,400 seconds;
there are at most 10,000 events. Event IDs are unique, case-sensitive strings of
1–64 letters, digits, underscores or hyphens. Ordering is `(due_tick, event_id)`.

The single controller must:

1. Split native steps at scheduled boundaries when exact delivery time matters.
   Events become eligible when `body_tick >= event.due_tick`; late delivery is
   explicit in the receipt's `issued_tick`. Body ticks are monotonic integers,
   derived from acknowledged elapsed solver time, with one tick equal to 20 ms.
2. Establish from native acknowledgment that `pending_meal` is false before
   calling `schedule.due(body_tick)`. Native `meal()` rejects a second meal
   pending digestion. Concurrent controllers are unsupported.
3. For the returned zero-or-one-element tuple, call `session.meal(event.meal)`
   exactly once. `due()` already reserves the event, so it cannot be returned
   again. Do not call `due()` while a reservation lacks its outcome.
4. On successful native return, call
   `schedule.record_accepted(event.event_id, acknowledgment)`. This checks the
   `ok` status, positive increasing sequence, `pending_meal=True`, and elapsed
   time matching the issuing tick. Keep the full native command journal as
   evidence; the schedule receipt is a compact association, not a substitute.
5. If command completion cannot be established, call
   `schedule.record_uncertain(event.event_id, reason)` and stop further intake
   issuance. This includes a failure between reserving and sending, a native
   timeout, or an invalid acknowledgment. There is deliberately no retry,
   rollback, or uncertainty-to-accepted transition. Preserve native logs and
   investigate; recreating a schedule must not be used to replay commands.

Accepted means the native command was acknowledged, not that the meal has been
absorbed. Step the continuing native session to consume its pending action
before delivering another event, including simultaneous events. The scheduler
does not merge meals to evade the native pending-meal contract. Record the
actual issue ticks when simultaneous events therefore arrive on later ticks.

## Persistence boundary

`checkpoint()` and `IntakeSchedule.from_checkpoint(mapping)` round-trip only
unissued scheduling configuration and its observed clock. Checkpointing is
refused after any reservation, including accepted or uncertain ones. A saved
schedule is not a native checkpoint, restart recipe, transaction log, or an
exactly-once guarantee across process crashes. Never restore an older schedule
against a native process that may already have consumed its commands. Native
state saving separately rejects pending meals and does not establish exact
controller restoration. Persistent receipts or post-issuance resumption need
an external coordinated controller design and are outside this module.

The module is bounded by validated event counts, times, identifiers, amounts,
and uncertainty messages. `receipts` returns immutable compact receipts;
configuration events and Meal objects are frozen. The public Python object is
intended for a trusted single controller, not hostile mutation or forgery.

## Verification

Run `.venv/bin/python scripts/verify_intake_schedule.py`. The pure fixture
covers food at tick 2, water at tick 4, deterministic simultaneous ordering,
duplicate IDs, fractional/nonfinite/out-of-bounds values, backward clock
rejection, explicit acknowledgment validation, due-once reservations, terminal
uncertainty, and the unissued checkpoint boundary. It launches no native jobs
and makes no claims about physiological fidelity or nutrient dynamics.
