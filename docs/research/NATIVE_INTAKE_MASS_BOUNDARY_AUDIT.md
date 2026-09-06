# Native intake mass boundary audit

Read-only source audit, 2026-09-05. No native process, build, or physiological
acceptance run was performed. Paths below are repository-relative. Source root
`BG` denotes `data/raw/physiology/biogears/projects/biogears/libBiogears/src`.

## Established native boundary

| Stage | Source location | What actually changes |
| --- | --- | --- |
| Meal command | `scripts/native_biogears_stream.cpp:57`, signed adapter `:121`, regional signed adapter `:129` | Builds `SEConsumeNutrients` in native units and calls `ProcessAction`; does not advance time. |
| Action storage | `BG/cdm/scenario/SEPatientActionCollection.cpp:624` | Copies the action into the pending ConsumeNutrients slot. Acknowledgment is acceptance into that slot. |
| Consumption | `BG/engine/Systems/Gastrointestinal.cpp:276` | Active PreProcess, lines 289–291: increment stomach contents, increment patient weight by native nutrition weight, remove the pending action. This is the external intake boundary. |
| Payload weight | `BG/cdm/patient/SENutrition.cpp:281` | Sums carbohydrate, fat, protein, calcium and sodium masses; converts water mL as grams. This is the donor's exact convention, not a new density law. |
| Payload transfer | `BG/cdm/patient/SENutrition.cpp:63` | Adds each native nutrient scalar and water to stomach contents. |

Consequently the native mass addition in kg is the native `GetWeight(kg)` of
the consumed action. For the present Meal interface, its source-equivalent
expression is `(carbohydrate_g + protein_g + fat_g + sodium_g + calcium_mg/1000
+ water_ml)/1000`. Configuration alone does **not** prove that this addition
occurred. Prefer querying the actual held native action, then establishing its
consumption; never increment mechanical mass on the scheduler's accepted state.

The retained patched GI file
`data/runtime/physiology/variants/whole_body_integrity_gi_water/Gastrointestinal.cpp:277`
contains the same boundary at lines 290–292. This inspection does not certify
that an arbitrary loaded library matches that file; existing source/build and
library manifests must bind the selected runtime implementation.

## Internal redistribution and other exits

The same PreProcess immediately calls `DigestNutrient`, then `ChymeSecretion`
and `AbsorbNutrients` (`Gastrointestinal.cpp:294`, `:302`). Stomach nutrient
losses become chyme masses (`:373`, `:387`, `:400`); stomach water moves into
chyme volume (`:450`). Absorption transfers chyme nutrients to vascular owners
(`:515` onward), and water uses the GI-to-CV circuit (`:550` onward). These are
internal transfers. They must not create additional external mechanical mass.
Raw stomach differences mix consumption and digestion within one solver step.
Chyme, blood, tissue stores and native patient weight are overlapping accounting
views, not independent masses to sum.

Patient weight alone is not a complete intake receipt or validated all-boundary
balance. For example native sweat decrements it in
`BG/engine/Systems/Energy.cpp:655`, while upstream urination directly resets
bladder volume and rebalances substances in `BG/engine/Systems/Renal.cpp:1565`.
Urine production and bladder storage are not themselves external outflow.
The audit did not establish a complete gas, sweat, urine or fecal mass ledger.

## Current observability gap

`scripts/native_body_ports.h:56` exports stomach quantities, absorption rate,
selected compartment volume/substance masses and urine observations. It exports
neither native patient weight nor a cumulative consumed-intake mass/count or
consumed command identity. `native_tissue_ports.h` likewise provides no intake
boundary counter. The inspected ConsumeNutrients consumer contains no cumulative
intake counter. This is a scoped finding about the held action path and current
adapter exports, not a claim that every donor subsystem lacks counters.

Stream receipts contain command sequence, simulation time, elapsed time, origin
and `pending_meal` (`native_biogears_stream.cpp:28`). Signed and regional adapters
emit the corresponding envelope at `:35`. A successful meal receipt has unchanged
elapsed time and a pending action. A subsequent acknowledged solver step may
show pending false, but the current envelope lacks an explicit association to
the consumed meal's command sequence. A failed/uncertain advance may have already
consumed it and must never be retried to obtain a receipt.

## Bounded source design for a later adapter-only implementation

Preserve native libraries and their physics. Add one shared observational helper
used by signed and regional adapters at **every** `AdvanceModelTime` call:
`scripts/native_biogears_signed.cpp:87`, `:103` and
`scripts/native_biogears_regional_signed.cpp:88`, `:104`.

Before advancing, copy the actual held action's native nutrient values and
`GetWeight(kg)`, paired with its accepted meal command sequence and process owner
identity. After a successful advance, establish that the previously held action
was removed by the audited active GI consumer. Commit the receipt only after
all signed-step validation also succeeds. Retain a consumption interval
`[start_tick,end_tick]`, rather than claiming an unobserved substep instant.

Export an integer consumed count, cumulative consumed mass kg, cumulative native
nutrient quantities/water volume, and last-consumed record containing meal
command sequence, advance sequence, interval ticks and native mass. The scheduler
can associate its event ID with the meal command sequence. Count/sequence must
remain exact integers, not floating observational-port values. At most one meal
is pending; later snapshots must retain cumulative values so missed responses
do not require replay. Keep counters process-scoped; loaded state initializes
an explicitly new counter epoch and does not retroactively credit stomach mass.

Also export native patient weight as a diagnostic comparison, not the intake
counter. Missing identity, nonfinite payload, unexpected pending transition,
failed advancement or signed validation failure must close/poison ownership
without returning a successful mass credit. Do not claim rollback. Persist any
attempt evidence before native advancement under the existing journal contract.
Adapter inference must remain explicitly tied to the audited sole consumer;
if consumer/removal paths change, require a native consumer-level export instead.

Mechanical integration subsequently needs its own once-only application keyed
by `(native owner, consumed count/meal sequence)`, explicit recipient location,
boundary momentum/velocity and energy handling. This audit establishes no
anatomical mass distribution or automatic whole-body conservation. Those are
separate from proving that the native external intake action was consumed.
