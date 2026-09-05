# Isolated exercise-demand integrity variant

`whole_body_integrity_energy` layers a single rebuilt `Energy.cpp` object over `whole_body_integrity_gi_water`. Every inherited object is hash-checked; original donor, GI-water library and other variants remain unchanged. It addresses the accumulated exercise demand and stale demand after stopping documented in `EXERCISE_ENERGY_AUDIT.md`.

## Explicit power partition

Let E be the previous exercise-demand power and T the previous total requested metabolic power. Before the existing thermoregulation calculation, non-exercise power is T − E. The existing basal/shivering/hypothermic equations operate on that component. With an active action, exercise demand updates independently:

`E_next = E + dt_s * (intensity * max_work_rate_W − E)`

This preserves the source gain (`MetabolicRateGain = m_dT_s`, nominal one-second relaxation) and its native 0.02 s integration step. The new total is the updated non-exercise component plus E_next. The thermal heat source is assigned after composition, using that same final requested total. No gain, work capacity, temperature threshold, metabolic substrate equation, tissue allocation fraction or gas conversion is fitted or replaced.

When the action is absent, E_next is zero. This removes requested exercise demand immediately on the next native update; it is not a newly invented post-exercise oxygen recovery model. Native circulatory, substrate and gas states continue from their current values. Existing tissue mechanics and anaerobic effects remain.

A nonfinite/negative exercise component or negative T − E throws an explicit inconsistent-partition exception. Historical states saved during the original runaway exercise are incompatible; start from a clean baseline or a state produced by the corrected variant. The variant does not silently reinterpret or repair those states.

## Native verification

`scripts/build_biogears_energy_variant.py` compiles the isolated object and relinks the complete parent object list. Its manifest retains the exact patch, compiler/link commands, source and object hashes, parent receipt and library hash. `scripts/verify_native_energy_integrity.py` compiles an independent probe using the donor's existing friend-test access. It runs actual native functions and whole-body advancement with explicit 0.02 s ticks, under an explicitly selected library and checked `ldd` identity.

Frozen-temperature branch tests hold the thermal node at 36.5°C and 37°C. They check onset, ten-second approach to the source target, exact total/non-exercise partition, final heat assignment and zero demand after stopping. The parent library reproduces the failures; the new variant satisfies these invariants. The fixture temporarily unlocks only its test node's temperature scalar; it does not change donor source or parameters.

Matched full-body rest, continuous 0.15 intensity and stop-at-30-second trials begin with the same saved state. They check finite outputs, bounded exercise demand, unchanged initial/pre-stop state, nonnegative stores, downstream oxygen response and reduced glycogen use after stopping. The retained probe `thermal_next_source_w` is inspected during PreProcess in branch tests; after a complete body tick the circuit has consumed/reset that next-step field, so it is not a cumulative thermal-energy observation.

## Limits

The source uses requested total metabolism as thermal input. This correction makes its bookkeeping internally consistent but does not equate thermal input with achieved tissue oxidation or subtract exported mechanical work. Native oxygen averages, substrate switching, biological debt and anaerobic demand still require their own validation. No body-wide energy conservation, clinical exercise calibration, long-duration endurance or mechanical-work-to-metabolism mapping is claimed.

## Retained acceptance result

The accepted native report is `data/derived/audits/energy-integrity-tfawudcr/report.json`.
All 26 checks pass, including unchanged rest trajectories, cold/neutral parent failure reproduction, corrected branch invariants, incompatible-state rejection, matched onset/stop histories and finite/nonnegative downstream stores. The new library SHA-256 is `23e7b7d9fdf3d3cb1b29b15f0cb393f0b290c48f180edad21bac84f9a6e95030`. All 360 linked objects are recorded. Inherited calcium transfer passes 6/6 cases and renal transfer passes 20/20.

| At 120 s | Exercise demand W | Total requested metabolism W | O₂ mL/min | Muscle glycogen g | Fatigue |
|---|---:|---:|---:|---:|---:|
| Rest | 0 | 92.900 | 262.938 | 592.159000 | 0 |
| Continuous intensity 0.15 | 150.120 | 241.585 | 480.184 | 592.158700 | 0 |
| Stop at 30 s | 0 | 92.511 | 261.931 | 592.158984 | 0 |

Rest is exactly equal to the GI-water parent at every sampled output. The post-exercise rest difference is a continuing-state effect; no baseline reset was performed. Oxygen consumption remains a realized, averaged tissue output, distinct from requested metabolic heat. These numbers are regression evidence, not a newly validated physiological exercise dataset.
