# Exertion v3: independent trajectory and initialization audit

**Recommendation: retain these as executing research experiments; do not promote
them as healthy generic-body predictions.** The thermal corrections improve core
temperature and the numerical/causal checks pass, but observed exercise glucose
and acid–base behavior are physiological readiness blockers. Meal thermogenesis
is structurally absent from requested metabolic power. No numerical code,
default, model input, or historical trajectory was changed by this audit.

The independent reader checked every value in all five 721-frame trajectories
(0–3600 seconds, 5-second samples) against the append-only `frames.jsonl`
journals. It did not use the cached acceptance predicates. Complete field
minima, maxima, their times, selected observations, source hashes and frozen
input copies are in
`data/derived/audits/exertion-v3-physiology-audit-j33kwere/audit.json`.

## Observed trajectories

All runs use the final humidity variant, library SHA-256
`ab485812968b17e16915b04a6e5144bb1e1ec7f5bf78213936c3334d087c0a5a`,
and the same saved initial state `cba7ffb5…68b3`. Exercise requests intensity
0.15 from 1800 to 2400 seconds; its additional requested power reaches
150.12 W. This is a lumped physiological request, not measured skeletal work.

| Protocol | Final core / skin °C | Aortic glucose min → final mg/dL | Arterial pH min–max → final | Peak HR / MAP |
|---|---|---|---|---|
| Rest | 37.070 / 29.724 | 73.119 → 75.769 | 7.410–7.447 → 7.445 | 73.84/min / 98.77 mmHg |
| Hydration | 37.067 / 29.779 | 73.128 → 75.553 | 7.410–7.472 → 7.471 | 73.84/min / 101.19 mmHg |
| Meal | 37.059 / 29.802 | 81.571 → 83.166 | 7.409–7.497 → 7.485 | 73.72/min / 102.18 mmHg |
| Exercise | 37.378 / 29.666 | 32.896 → 48.683 | 7.299–7.498 → 7.323 | 94.78/min / 107.58 mmHg |
| Meal + exercise | 37.366 / 29.752 | 44.127 → 58.091 | 7.342–7.546 → 7.354 | 92.79/min / 108.06 mmHg |

Core peaks are 37.412 °C and 37.410 °C in the exercise cases. All gases are
finite, but this is weaker than physiological agreement: exercise arterial O₂
reaches 79.16 mmHg and CO₂ spans 33.22–42.86 mmHg. Combined exercise spans
79.63–104.95 mmHg O₂ and 33.78–43.36 mmHg CO₂. Peak cardiac output is about
6.74 L/min and respiratory rate about 22/min. These feedback paths are active;
their gains, transients and recovery are not validated by merely changing in
the expected direction.

Exercise glucose is below 54 mg/dL from the 2035-second sample through the end,
with nadir at 2460 seconds. Combined exercise is below 54 from 2175 through
3365 seconds, with nadir at 2620. The simulated muscle vascular glucose becomes
nearly exhausted: minima 0.0164 and 0.0234 mg/dL. The
[ADA/EASD joint position statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC6518070/)
identifies glucose below 54 mg/dL as clinically important hypoglycemia. This is
a screening comparison: a source-model aortic concentration is not a validated
clinical plasma measurement, and no diagnosis of a real person is implied.

Both exercise runs also cross outside the usual arterial pH interval of
7.35–7.45 ([Clinical Methods](https://www.ncbi.nlm.nih.gov/books/NBK371/)).
Hydration and meal develop alkalemic drift without exercise. At 3600 seconds,
20 minutes after exercise stops, aortic lactate is still rising and reaches
65.14 mg/dL (exercise) or 86.06 mg/dL (combined), compared with 3.33 mg/dL at
rest. Native substrate availability, glucose production/transport, lactate
clearance, buffering and gas-feedback timing need a separate source audit.
This review establishes the anomalies, not their unique cause.

The skin output is the actual source area-weighted mean of six thermal nodes,
not torso temperature alone. It falls from 31.837 to about 29.7–29.8 °C during
the hour, so the initial state is not thermally stationary. The held matched
JOS-3 comparison has mean skin 32.861 °C standing or 33.128 °C lying after an
hour. Different physiological states/controllers prevent treating that
difference as a calibrated correction factor. Native posture remains standing
for heat transfer, with no mattress or wet-garment contact boundary.

## Meal, stores and requested energy

The nominal meal adds 60 g carbohydrate, 20 g protein, 20 g fat, 1 g sodium,
300 mg calcium and 500 mL water. At an hour, stomach carbohydrate is still
35.977 g and fat 8.069 g, whereas protein is almost depleted (0.00785 g).
The trajectory has not completed end-to-end digestion. Relative to matched
hydration, meal glucose differs by up to 14.462 mg/dL and insulin synthesis by
1.114 pmol/min. These are responses of the held equations, not human fits.

Meal liver glycogen stays exactly 117 g. Source `Tissue::Initialize` and
`Hepatic::SetUp` initialize and cap it at 0.065 × 1800 g. Hepatic glycogenesis
requires stores below that cap. Hydration loses about 1.357 g, so the contrast
is **sparing an initially full store**, not observed refilling of depleted
glycogen. This saturation should constrain claims about meal storage responses.

Requested metabolic power in rest, hydration and meal remains exactly
82.1822014 W. Source `Energy::CalculateMetabolicHeatGeneration` uses temperature
and basal metabolic rate, and `Exercise` adds exercise demand. There is no meal
or digestive-energy term in this producer. Thus absent meal thermogenesis is
a missing interaction, not a fitted zero response. The source's metabolic
consumption code does process nutrients; nutrient processing alone does not
create a meal-to-metabolic-heat feedback path.

The requested-power partition BMR + exercise demand closes numerically. This
is not a whole-body chemical-energy ledger. `Energy::PreProcess` comments out
`ManageEnergyDeficit()`, so the observed all-zero `energy_deficit_w` cannot
demonstrate adequate substrate delivery. Tissue anaerobic branches can consume
glucose/glycogen and produce lactate while that field remains zero. Also,
stored-fat changes include release into other pools; they cannot be equated
with fat oxidation from this single scalar.

All observed local masses, concentrations and volumes remain nonnegative and
finite. Only `respiratory_cycle_fraction` is missing, at the initial observation
in each run. Nonnegativity does not establish complete mass conservation:
the 180 ports do not enumerate every intracellular or overlapping vascular
store. Urine and filtration outputs show substantial ongoing dynamics; for
example resting GFR spans 79.9–167.1 mL/min. Their calibration remains separate.

## Initial-state meaning and default promotion

The common initial stomach already contains **500 mL water, 1 g sodium and
500 mg calcium**, even in the rest and exercise controls. Hydration and meal
add their loads on top. Rest blood volume consequently rises from 5534.14 to
5763.15 mL. Calling this a fasting or empty-gut resting control would be wrong.
A fasting experiment needs an explicitly generated and accepted physiological
fixture, not silent deletion of these stores from a checkpoint.

Current code paths differ:

| Entry point | Current initialization | Required treatment for any later promotion |
|---|---|---|
| Web batch scenario (`Jobs.submit`) | Rejects client `state_path`; `NativeConfig.state_path=None` causes fresh initialization under selected library. Canonical body forces `IHMGenericMale`. | Adding a research selector does not automatically load the legacy baseline. Validate each freshly initialized profile and retain its state identity; distinguish this path from checkpoint replay. |
| `Jobs.list` / scenario selector | Prefers `whole_body_integrity_depletion`; final thermal variants are not listed here. | A later selector/default change must be tied to accepted variant/profile evidence, not just presence of a library file. |
| Low-level `NativeConfig` and `SessionConfig` | Default `upstream`, `StandardMale`, no saved state. | Keep source-reference behavior explicit, or introduce a separately named canonical initialization policy. Do not conflate StandardMale with the generic canonical body. |
| `SystemicConfig`, including implicit `body-systemic` materialization | Defaults to depletion library plus legacy `native_baseline_v1`. | Change variant and initialization policy together, or require both explicitly. Swapping only the variant silently carries an incompatible old thermal state forward. |
| Explicit saved-state `NativeSession` | Copies the supplied state and records its actual patient identity; variant and state are independently selectable. | Bind an accepted initialization ID to state SHA, library SHA, patient and boundary conditions. Reject an unintended mismatch; keep deliberate historical controls explicitly labeled. |

The final research fixture is indexed by
`data/derived/audits/thermal-final-initial-state-ctkr8bap/initial-state.json`.
It is fresh source stabilization under the final library, not thermal or
all-system equilibrium. Native checkpoint continuation is not established as
exact for every controller latch. A live session should remain continuous for
coupled experiments; loading the same fixture makes paired comparisons fair,
but does not prove exact equivalence to uninterrupted post-initialization state.

Do not mutate any default on the evidence reviewed here. Keep thermal
corrections available as research variants, retain these failing physiological
responses, and require corrected/validated glucose and acid–base behavior plus
final six-hour integrated acceptance before ordinary healthy-body publication.
Meal thermogenesis should remain explicitly absent until a source-attached
mechanism and measurement comparison are implemented.

Reproduce the independent read-only audit:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_exertion_trajectories.py
```
