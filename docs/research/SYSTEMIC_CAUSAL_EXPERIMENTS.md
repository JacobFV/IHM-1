# Shared systemic experiments and canonical projection

2026-09-05. `human.materialize('body-systemic', ...)` constructs an evidence-bound native experiment. Rest, matched hydration, mixed meal, apnea/restoration, exercise/recovery and mixed-meal/exercise protocols use one continuing physiological state per condition. The default corrected source variant is `whole_body_integrity_depletion`. Samples retain source identities, native units, exact action times, required-store checks, and state checkpoints. Exercise and long-horizon thermal acceptance require separate source audits; a protocol being executable does not establish physiological fidelity.

Source hashes are captured before execution. A causal comparison requires identical initial state, native library, executable, resolved dependencies, timestep and executing-source identities, equal initial observations, matching sample clocks, and equal observations before the first differing intervention. Empty collections and unpaired interventions cannot report a successful causal contrast. Required stomach and glycogen/protein/fat stores cannot pass by being null; negative stores remain explicit failures. These checks do not establish global elemental or energy conservation.

Input bytes are copied into each experiment's `inputs/` tree and bound by
`frozen-sources.json`. Historical results therefore retain their actual source
code, executable, variant library and initial state when current development
changes those files. Altered archive bytes are rejected; an unarchived stale
input is still rejected. The original executing-source identities remain in the
record. Display projection has its own current code/anatomy identities. This
separates reproducible historical execution from the projection that displays it.

Default display publication also requires a passed matched-contrast report
bound to the exact input result hashes. It recomputes those contrast checks,
writes immutable display records, and updates the index only after all requested
projections succeed. Failed runs remain available on disk as research evidence.

Publication recomputes required-store validity, uniform sample clocks and the
named protocol's actual intervention schedule from observations. A cached
`passed` value cannot override negative/missing stores. The embedded native
manifest must equal the retained manifest. Sampled frames must agree with the
retained native command acknowledgments; commands, clocks, intervention values,
graceful termination and both checkpoint hashes are checked and included in
the display's transitive source identities. The independent regression changes
derived samples, labels, manifests and checkpoints in temporary copies and
requires rejection. Historical input hashes remain distinct from the current
code performing these acceptance checks.

The mixed meal contains 60 g carbohydrate, 20 g protein, 20 g fat, 1 g sodium, 300 mg calcium and 500 mL water. Its hydration control receives the same water and electrolytes. This separates a macronutrient perturbation from that added fluid/electrolyte input. The shared initial state already contains source stomach water/electrolytes and is not described as an empty gut. Six-hour runs preserve the native 0.02 s timestep and sample every 30 s. Glucose, insulin synthesis and glycogen contrasts require changes after intervention; their magnitudes are execution evidence, not fitted human response curves.

The three-minute respiratory pair samples at 10 Hz. Apnea is applied at 30 s and removed at 90 s. Native output shows a maximum between-condition lung-volume difference of 737.112 mL and arterial CO2 pressure difference of 6.17523 mmHg. Required local stores remain present and nonnegative. Raw traces, actions, states, manifests and contrasts are under `data/derived/systemic/respiratory_v3/`.

The fresh `respiratory_v4` pair repeats this experiment on the depletion-corrected
library with180 ports and prestart detached input selection. The two causal
contrasts reproduce the same magnitudes. This pair now supplies the default
rest/apnea projections; all v3 artifacts remain available and verifiable.

## One-body inspection

`build_systemic_display.py` creates views in the same canonical anatomical frame. Dense respiratory records drive the existing three-mode thoracic solver from actual native total-volume change. All five lobes receive a common fractional expansion: this is an explicit uncalibrated projection assumption, not individually observed lobe gas state. The projection includes diaphragm, ribs and skin, with a discrete mechanical energy audit. It is one-way native-to-mechanical projection; its constraint reaction is not returned as native pleural pressure or metabolic work.

Coarse nutrient records retain reference geometry. No synthetic heartbeats or respiratory cycles are inserted between 30 s observations. Missing transforms mean identity, rather than stale displacement from a previous experiment. This distinction remains visible in the app.

Actual sampled variables also produce finite-window Laplace transforms and Hann Welch spectra using the existing temporal library. Time is seconds, frequency is Hz, and `s = sigma + 2 pi i f`. Laplace integration subtracts the full-record mean and uses the trapezoidal rule. Maximum frequency respects output Nyquist: 5 Hz for 0.1 s samples, 1/60 Hz for 30 s samples. These are finite-record descriptors, not physiological poles, causal transfer identification, or evidence of coupling on their own.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_systemic_experiments.py \
  --protocols rest apnea --seconds 180 --sample-interval .1 \
  --output data/derived/systemic/respiratory_v3
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_systemic_experiments.py \
  --protocols hydration meal --seconds 21600 --sample-interval 30 \
  --output data/derived/systemic/six_hour_v3 --workers 2
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_systemic_display.py \
  data/derived/systemic/respiratory_v3/rest data/derived/systemic/respiratory_v3/apnea
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_systemic_projection.py
```

These commands require fresh output directories for experiments. Completed canonical display records may be rebuilt from unchanged source experiments. An earlier six-hour attempt was explicitly interrupted for provenance correction, and a first respiratory run exposed a floating remaining-horizon comparison. Their partial data and failures are retained; neither is reported as completed evidence.

The completed `six_hour_v2` pair failed because depleted stomach sodium became
an invalid scalar. Its native macronutrient perturbation produced maximum
between-condition differences of 9.73913 mg/dL aortic glucose, 1.75324 pmol/min
insulin synthesis and 18.1011 g liver glycogen, but the failed local-store check
prevents publication as an accepted default experiment. The dry-gut correction
preserves verified depletion as zero and unknown inputs as unknown; it also fixes
zero-water division and carrier-availability defects, with native branch tests.

The corrected `six_hour_v3` pair completes all 721 samples per condition with
present, nonnegative required stores and the same three macronutrient contrasts.
It remains excluded from ordinary meal-response display: core temperature falls
to 24.163 °C with the meal and 24.194 °C with hydration. Passing local mass and
causal checks is insufficient evidence of physiological homeostasis. Native
thermal boundary algebra and initial-state compatibility are being audited
before this run can support a normothermic generic-body prediction.

The first energy-corrected one-hour group `exertion_v2` records completed
hydration, meal and meal/exercise trajectories, but its rest and exercise
conditions stop near 1199 s on the mixed-unit water-depletion error. Those
conditions have no completed `systemic.json`; the whole group is unsuccessful.
New group runners retain per-protocol errors while collecting other workers'
outcomes, reject changed sources before a queued worker starts, and cannot
declare a partial group successful. Exercise contrast acceptance additionally
requires finite bounded demand, a nonnegative non-exercise power partition,
explicit active/post-stop samples and cleared demand after stop.

## Fresh thermal initialization and integrated replay

The runner accepts `--state` to bind every condition to the same explicitly
retained initial state. Its group receipt now includes that state, the actual
variant library and stream executable before any workers start, in addition to
the implementation files. A queued condition rejects changed input bytes.
The four successive thermal research variants are available to native batch
and session configuration; registering them does not change the app default.

The final variant, `whole_body_integrity_evaporation_humidity`, includes clothing
and film resistance, regional skin perfusion, sweat-area and evaporative
capacity, and local humidity corrections. See
[the thermal audit](LONG_HORIZON_THERMAL_AUDIT.md) and
[the evaporation audit](NATIVE_EVAPORATION_AUDIT.md). Its fresh initialization is
retained under `data/derived/audits/thermal-final-initial-state-ctkr8bap/`, with
state SHA-256 `cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3`.
It starts at 37.0481°C and is not claimed to be a thermal equilibrium.

`exertion_v3` completes all five one-hour conditions from that same state, with
721 observations each. All tracked required local stores remain nonnegative
and seven named intervention contrasts pass. Core temperatures stay between
37.0481 and 37.4119°C. Exercise raises requested metabolic power by 150.12 W;
the demand clears after its ten-minute interval. These are numerical execution
results, not independently fitted human responses. In particular, reported
total metabolic power does not rise after the meal alone; this record does not
establish diet-induced thermogenesis. Independent trajectory review additionally
finds exercise aortic glucose as low as32.896mg/dL and arterial pH7.2987; the
combined meal/exercise condition reaches44.127mg/dL and pH7.5456. Successful
execution and finite stores do not make these ordinary healthy responses.
Reported zero energy deficit is also insufficient: the source call to
`ManageEnergyDeficit()` is commented out. The six-hour `six_hour_v4` replay
also completes721 observations per condition with nonnegative tracked stores
and three causal contrasts. Meal/hydration final core temperatures are36.7791
and36.8491°C, respectively, but arterial pH reaches7.5088/7.5108. It remains
excluded by the independent homeostasis screen. The source counterion audit
finds nutrient cotransport moving existing gastric/chyme sodium into blood
without an implemented chloride counterpart; this needs a separate explicit
electrolyte/charge mechanism, not a pH reset.

The default-view gate independently requires finite observed aortic glucose
at least70mg/dL and, for ordinary rest/meal/exercise protocols, arterial
pH7.35–7.45. These values correspond to retained native screening thresholds;
the glucose screen uses the observed aorta while the native event uses vena
cava, and the pH screen does not diagnose a metabolic disorder without
bicarbonate and cause attribution. The deliberately imposed apnea protocol
permits respiratory pH perturbation and still requires finite observations.
The glucose boundary is consistent with the
[NIDDK screening description](https://www.niddk.nih.gov/health-information/diabetes/overview/preventing-problems/low-blood-glucose-hypoglycemia).
These are rejection screens, not population calibration or a complete set of
physiological validity tests. No part of the failed five-condition group is
published by extracting only favorable conditions. Complete inputs and
trajectories remain available for diagnosis.

Unresolved biology remains substantial: swallowing, peristalsis and enteric control, enzymatic compartments and bile, intestinal lipid lymph transport, incretins, colon/microbiome/stool, independently calibrated macro digestion kinetics, exact elemental budgets, individual phrenic anatomy, and force/work exchange with whole-body mechanics. Native physiology implements coupled reduced mechanisms; it is not yet a complete calibrated human.
