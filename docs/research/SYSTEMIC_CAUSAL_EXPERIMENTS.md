# Shared systemic experiments and canonical projection

2026-09-05. `human.materialize('body-systemic', ...)` constructs an evidence-bound native experiment. Rest, matched hydration, mixed meal, apnea/restoration, exercise/recovery and mixed-meal/exercise protocols use one continuing physiological state per condition. The default corrected source variant is `whole_body_integrity_renal`. Samples retain source identities, native units, exact action times, required-store checks, and state checkpoints.

Source hashes are captured before execution. A causal comparison requires identical initial state, native library, executable, resolved dependencies, timestep and executing-source identities, equal initial observations, matching sample clocks, and equal observations before the first differing intervention. Empty collections and unpaired interventions cannot report a successful causal contrast. Required stomach and glycogen/protein/fat stores cannot pass by being null; negative stores remain explicit failures. These checks do not establish global elemental or energy conservation.

The mixed meal contains 60 g carbohydrate, 20 g protein, 20 g fat, 1 g sodium, 300 mg calcium and 500 mL water. Its hydration control receives the same water and electrolytes. This separates a macronutrient perturbation from that added fluid/electrolyte input. The shared initial state already contains source stomach water/electrolytes and is not described as an empty gut. Six-hour runs preserve the native 0.02 s timestep and sample every 30 s. Glucose, insulin synthesis and glycogen contrasts require changes after intervention; their magnitudes are execution evidence, not fitted human response curves.

The three-minute respiratory pair samples at 10 Hz. Apnea is applied at 30 s and removed at 90 s. Native output shows a maximum between-condition lung-volume difference of 737.112 mL and arterial CO2 pressure difference of 6.17523 mmHg. Required local stores remain present and nonnegative. Raw traces, actions, states, manifests and contrasts are under `data/derived/systemic/respiratory_v2/`.

## One-body inspection

`build_systemic_display.py` creates views in the same canonical anatomical frame. Dense respiratory records drive the existing three-mode thoracic solver from actual native total-volume change. All five lobes receive a common fractional expansion: this is an explicit uncalibrated projection assumption, not individually observed lobe gas state. The projection includes diaphragm, ribs and skin, with a discrete mechanical energy audit. It is one-way native-to-mechanical projection; its constraint reaction is not returned as native pleural pressure or metabolic work.

Coarse nutrient records retain reference geometry. No synthetic heartbeats or respiratory cycles are inserted between 30 s observations. Missing transforms mean identity, rather than stale displacement from a previous experiment. This distinction remains visible in the app.

Actual sampled variables also produce finite-window Laplace transforms and Hann Welch spectra using the existing temporal library. Time is seconds, frequency is Hz, and `s = sigma + 2 pi i f`. Laplace integration subtracts the full-record mean and uses the trapezoidal rule. Maximum frequency respects output Nyquist: 5 Hz for 0.1 s samples, 1/60 Hz for 30 s samples. These are finite-record descriptors, not physiological poles, causal transfer identification, or evidence of coupling on their own.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_systemic_experiments.py \
  --protocols rest apnea --seconds 180 --sample-interval .1 \
  --output data/derived/systemic/respiratory_v2
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_systemic_experiments.py \
  --protocols hydration meal --seconds 21600 --sample-interval 30 \
  --output data/derived/systemic/six_hour_v2 --workers 2
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_systemic_display.py \
  data/derived/systemic/respiratory_v2/rest data/derived/systemic/respiratory_v2/apnea
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_systemic_projection.py
```

These commands require fresh output directories for experiments. Completed canonical display records may be rebuilt from unchanged source experiments. An earlier six-hour attempt was explicitly interrupted for provenance correction, and a first respiratory run exposed a floating remaining-horizon comparison. Their partial data and failures are retained; neither is reported as completed evidence.

Unresolved biology remains substantial: swallowing, peristalsis and enteric control, enzymatic compartments and bile, intestinal lipid lymph transport, incretins, colon/microbiome/stool, independently calibrated macro digestion kinetics, exact elemental budgets, individual phrenic anatomy, and force/work exchange with whole-body mechanics. Native physiology implements coupled reduced mechanisms; it is not yet a complete calibrated human.
