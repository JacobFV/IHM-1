# Persistent native physiological state

2026-09-05. `NativeSession` owns a continuing BioGears process and advances its native 0.02 s solver on integer ticks. Commands change that process, and snapshots observe the same blood, gas, nutrient, endocrine, renal and thermal state. This is an executable interface for online coupling; it does not itself supply a whole-body mechanical feedback law.

The bridge exposes 180 explicitly named quantities, including chemical respiratory requests, applied respiratory pressure, airflow, circulating nutrients and hormones, stomach/chyme contents, glycogen, exercise demand, basal metabolic power and maximum work rate. Optional absent/native nonfinite quantities are JSON null. Required storage quantities must be present in systemic experiment acceptance. Overlapping compartment names are not independent mass reservoirs and must not be summed.

`SessionConfig` permits up to 24 simulated hours per process; output cadence is independent of the native timestep. `Meal` specifies carbohydrate/protein/fat/sodium in g, calcium in mg and water in mL. Apnea severity uses the native action. Exercise intensity requests native lumped physiological demand; it is not calculated mechanical work or an IBM cortical command.

## Reproducible boundary

The build script writes a source/executable manifest. Startup verifies those identities, the selected variant library, and actual resolved dependency paths and hashes. The saved state supplies patient identity when loading a state; an unused initialization default cannot relabel it. Command receipts record attempts before writing to stdin and acknowledgments afterward. Protocol sequences and native elapsed time are checked. A reader thread drains native output, preserves logs and propagates engine failures. Invalid Python inputs do not issue commands or advance the native clock.

Before spawning the process, startup archives the resolved native environment,
materializes a detached tree of its runtime resources and selects that tree
through session-owned links. A saved input state is copied and hash-checked;
the actual command reads that copy. `execution-inputs.json` binds these selected
paths to the native and environment manifests. Dynamic libraries still load
from their recorded resolved paths; they are archived for reproduction, not
redirected to a new loader environment. Earlier post-start archives retain
their original capture phase and cannot acquire a retrospective startup claim.

Only one controller thread may command a session. A timeout or failed acknowledgment terminates that owned process. This is fail-closed execution, not rollback of already advanced native state. Checkpoints retain unique copies of native XML and hashes; exact continuation of all upstream controller latches remains unestablished.

Native meal actions are deferred until the next physiological PreProcess. The upstream action collection has one pending nutrition slot. Two consecutive meal commands previously replaced the first; both adapter layers now reject a second pending meal. Native serialization of pending nutrition produces XML that its own loader rejects, so save is also rejected until the action is consumed. Closing with an unconsumed meal terminates the owned process, preserves its receipts and raises an explicit error. Raw bridge `quit` rejects a pending meal, and stdin EOF with a pending meal exits nonzero. No additional solver step is silently inserted to make an action look consumed.

## Native evidence

`scripts/verify_native_session.py` compares resting stream and legacy batch outputs from the same saved state. At one second the observed differences are 4.51e-7 mL total lung volume, 1.42e-14 /min heart rate and 8.20e-8 mmHg MAP, within the batch CSV precision. A 10 g carbohydrate meal enters once after one native step. Invalid clock requests leave state time unchanged; repeated 0.1 s steps can reach a 0.3 s horizon without a floating upper-bound rejection.

A separate 75.04 s action experiment suppresses respiratory volume excursion to zero in its settled apneic window and observes 666.225 mL excursion after restoration. These are source-model execution checks, not clinical validation or population confidence intervals.

The first one-hour exercise contrast exposed an adapter defect: changing the
generic intensity scalar on a default `SEExercise` leaves its action mode at
`NONE`. Native acknowledgment therefore did not establish execution. The stream
now constructs the action from `SEExercise::SEGeneric`, as the batch adapter
already did. A matched 120 s native regression at intensity 0.15 observes
exercise-minus-rest differences of 3.02431 W in reported metabolic rate and
207.80081 mL/min in oxygen consumption. Both downstream effects are required by
the regression. The original no-effect runs remain in
`data/derived/systemic/exertion_v1/` with failed contrast checks; they are not
exercise validation. This comparison tests causal execution, not physiological
agreement between those two energy-related quantities.

The independent source audit then identified a more consequential native
bookkeeping defect: thermal control repeatedly reset total metabolic rate while
exercise demand accumulated, and stopping the action left demand active.
`whole_body_integrity_energy` separates non-exercise and exercise requested
power and clears stopped demand. The updated persistent test observes
148.68468 W and 217.24680 mL/min exercise-minus-rest differences after 120 s,
requires demand to reach the requested bounded target, and verifies exact
zero demand on the first step after stop. The native branch/full-body audit is
documented in [ENERGY_DEMAND_INTEGRITY.md](ENERGY_DEMAND_INTEGRITY.md).

Longer execution also exposed mixed-unit water-depletion roundoff: a complete
debit from an L-stored scalar could leave a negative remainder of about
5.3e-20 mL. The next preflight correctly rejected that value.
`whole_body_integrity_depletion` sets a known completely emptied water pool to
exact zero; it does not normalize arbitrary negative input. Earlier failed
one-hour runs are retained under `data/derived/systemic/exertion_v2/`.

The `whole_body_integrity` and `whole_body_integrity_renal` variants preserve prior source libraries and record isolated calcium and renal mass-transfer corrections. Their local native branch evidence is documented separately in [GI_MASS_INTEGRITY.md](GI_MASS_INTEGRITY.md) and [RENAL_MASS_INTEGRITY.md](RENAL_MASS_INTEGRITY.md).

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_native_adapter.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_session.py
```
