# Native exercise demand audit

The 120-second oxygen/metabolic-rate mismatch is primarily a reproduced demand bookkeeping defect, not evidence of a valid exercise response. The stream constructor correction makes the action active, but does not correct the donor energy lifecycle.

Source: retained BioGears revision `3f16a5fa1dade9c511b88d923606fa51cc35e95d`. Paths below are relative to `projects/biogears/libBiogears/src` in the acquired donor.

## Causal chain

1. `engine/Systems/Energy.cpp:200` calls metabolic heat calculation before exercise. At the saved baseline core temperature 36.575°C, the shivering branch (`:570–577`) resets total metabolic rate each step, even with an active exercise action.
2. `Energy.cpp:341–347` increments both total metabolic rate and a separate `ExerciseEnergyDemand` power scalar by `dt * (BMR + desired_work - current_total_rate)`. The next shivering update discards the exercise contribution to total rate while the separate demand continues accumulating.
3. `engine/Systems/Tissue.cpp:896` reads that accumulated power and multiplies by the timestep. `:1152–1159` allocates exercise demand 80% to muscle, 20% to fat. The demand is read regardless of whether an exercise action remains active.
4. `cdm/patient/actions/SEExercise.cpp:83` marks zero generic intensity inactive. `cdm/scenario/SEPatientActionCollection.cpp:675–685` successfully removes it and returns true. However, the no-action branch in `Energy.cpp:331–337` returns without clearing exercise demand.

This is not a missing factor of 60 in the demand conversion. The exercise scalar is power, and Tissue's conversion to interval energy is dimensionally appropriate. Its lifecycle is inconsistent with the separately overwritten total metabolic rate.

## Matched native evidence

`scripts/verify_native_exercise_energy.py` loads the same canonical stabilized state and selected library into three continuing sessions: rest, intensity 0.15, and intensity 0.15 stopped at 30 s. It samples native telemetry and saved XML, with hashes of state, donor files and session manifests. The requested work is 0.15 × 1000.8 = 150.12 W.

| Trial/time | Exercise demand W | Total metabolic rate W | O₂ mL/min | Muscle glycogen g |
|---|---:|---:|---:|---:|
| Rest, 120 s | 0 | 92.900 | 262.938 | 592.159 |
| Continuous, 1 s | 144.170 | 90.928 | 251.157 | 592.200 |
| Continuous, 30 s | 4310.744 | 91.915 | 471.417 | 551.093 |
| Continuous, 120 s | 17017.858 | 95.925 | 470.738 | 0 |
| Stop at 30 s, observed 120 s | 4310.744 | 93.013 | 463.891 | 199.022 |

The stop acknowledgment is successful and saved `ExerciseData` disappears. Demand nevertheless remains exactly 4310.744 W through 120 s. The two exercising branches have identical states before stopping. Continuous exercise reaches fatigue 0.99048, versus zero at rest; the small 3.024 W total-rate contrast therefore masks a severe unphysical accumulated demand.

## Legitimate distinctions and remaining limits

Oxygen consumption is calculated from actual substrate/O₂-limited reactions, includes hepatic flux, and is averaged over approximately 50 native steps (`Tissue.cpp:1483–1509`). It is not an instantaneous mirror of requested metabolic power. Substrate pathways have different stoichiometry; anaerobic glucose/glycogen contributes energy without proportional oxygen consumption (`:1383–1450`). These mechanisms explain delayed/limited oxygen response to the runaway request, not why a 150.12 W request becomes 17 kW.

`heatGenerated_kcal` is local accounting, with its output probe commented out (`:1568`); the thermal circuit heat source is instead set from requested total metabolic rate (`Energy.cpp:590`). Consequently, correcting demand bookkeeping alone does not establish an achieved-energy or external mechanical-work balance. No new caloric-equivalent equation or clinical calibration is justified by this audit.

Required correction is an explicit non-exercise/exercise power partition, independent bounded exercise-demand relaxation, and removal of the requested exercise component when the action stops. Preserve the source shivering and substrate equations. Test cold and neutral branches, onset/steady/stop, actual downstream stores and inherited fixes in a separately named library. Stream zero-intensity action acceptance itself requires no correction; add regression coverage for downstream demand removal rather than accepting an acknowledgment as physiological recovery.

The audit script intentionally exits successfully when the known defects are reproduced. Its `audit_reproduced` field is not a physiological acceptance result. The original reproduction is retained at `data/derived/audits/exercise-energy-5fc8vhal/report.json`; all donor sources remain unchanged.
