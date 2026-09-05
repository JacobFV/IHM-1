# Native whole-body backend execution

The collected **BioGears C++ engine was built natively for ARM64** and exercised through its official `PhysiologyEngine` API. The thin adapter in `scripts/native_biogears_rest.cpp` initializes the upstream `StandardMale.xml` patient, requests physiological outputs and advances the original integrated engine. It supplies no replacement equations or fitted constants.

Source authority: `BioGearsEngine/core` commit `3f16a5fa1dade9c511b88d923606fa51cc35e95d`. The compiled version function reports `8.0.0+Source-Unknown` because this shallow checkout lacks the tag/version metadata expected by the upstream CMake script. Do not identify this executable as the downloaded 8.2.0 release.

## Reproduce

From the repository root, with the official source collection present:

```sh
.venv/bin/python scripts/build_native_biogears.py
.venv/bin/python scripts/run_native_biogears.py --seconds 60
.venv/bin/python scripts/run_native_biogears.py --summarize-only
```

The build script downloads Ubuntu dependency `.deb` files and extracts them into `data/runtime/physiology/sysroot`; it does not install system packages. It currently targets the Ubuntu/Debian package environment used here (Ubuntu 24.04, GCC 13, CMake 3.28, ARM64), and uses 16 compile jobs. The required packages are Eigen3, Xerces-C runtime/headers, CodeSynthesis XSD (`xsdcxx`) and log4cpp runtime/headers. Native artifacts, package hashes, build logs and dependency files remain in `data/runtime/physiology`. The adapter has an absolute RPATH into this workspace's build library directory, so rebuild after relocating the workspace.

The run script creates symlinks to staged upstream resources in its output folder, then saves the **original BioGears tracker CSV**, engine log, stdout, stabilized/final engine XML states, execution status and numeric summary under `data/derived/physiology/biogears_native_run/`. The CSV is not resampled or normalized. The runner rejects nonnumeric/nonfinite values, non-increasing time, ragged rows and fatal/unknown-request log messages, even when BioGears itself returns zero.

## Build findings and resolved failures

The commands and corresponding log files are retained under `data/runtime/physiology`:

1. `docker info --format '{{.Architecture}} {{.DockerRootDir}}'` failed: permission denied connecting to `/var/run/docker.sock`. Docker was not used.
2. The official `Biogears-8.2.0-Linux.tar.gz` release was downloaded and inspected. `file .../bin/bg-cli` identified x86-64, incompatible with this ARM64 host without emulation. It was not executed.
3. Initial `cmake -S data/raw/physiology/biogears -B data/runtime/physiology/biogears-build -DCMAKE_BUILD_TYPE=Release -DBiogears_BUILD_DOCUMENTATION=OFF` failed: `Unable to find Eigen3 which is marked REQUIRED` (`cmake-configure.log`). Local `.deb` extraction and `CMAKE_PREFIX_PATH`/`CodeSynthesis_EXECUTABLE` resolved dependencies (`cmake-configure-deps.log`).
4. `cmake --build ... --target bg-cli -j 16` initially failed with XSD `unable to open in write mode`, because generated schema subdirectories were absent (`build.log`). Creating `generated/Release/biogears/schema/{cdm,biogears}` resolved it.
5. Compilation then failed to locate `biogears/schema/cdm/Properties.hxx` (`build-retry.log`). The source build mixes configuration-specific and configuration-independent generated include paths. A **build-directory-only** symlink `generated/biogears -> Release/biogears` resolves this.
6. GCC 13 rejected `DataTrack.cpp` because `std::istream_iterator` was used without `<iterator>` (`build-retry2.log`). The generated Makefile compilation recipe for this single file was given `-include iterator`. This changes an include dependency, not physiological behavior. The archived source remains unchanged.
7. The native engine libraries built, but `bg-cli` failed on the missing optional embedded-IO `biogears/io/directories/config.h` with `Biogears_BUILD_IO_LIBRARY=OFF` (`build-retry3.log`). The direct API adapter uses the complete engine libraries and filesystem resources; the CLI is not required.
8. The consumer adapter needed Eigen's explicit include path and `<cassert>` before upstream template headers. Both were supplied in the adapter/build command. The final adapter compile succeeded (`adapter-build-final.log`).

The reproducible build script incorporates the three build-directory/include fixes and builds `libbiogears` plus runtime staging directly. No physiological source file was edited. The downloaded model asset hashes can still be checked with `scripts/collect_physiology.py --verify`.

## Execution meaning and limitations

This is an **integrated whole-body engine run**, rather than independent submodels with copied defaults. Its state and original solver retain cross-system cardiovascular, respiratory, renal, endocrine, metabolism/heat and blood-chemistry behavior.

The patient is resting with no insult or intervention. No explicit posture action was configured, so this is not independently established as a supine-posture simulation. Initialization includes BioGears' own circuit tuning and dynamic stabilization; simulation time is reset to zero once active. The log reports primary/secondary convergence durations and the final homeostasis interval; these are distinct from the subsequent 60-second recorded trajectory.

The corrected adapter requests 50 samples per second, matching the 0.02-second native step, and includes instantaneous `ArterialPressure` and `TotalLungVolume` to resolve heartbeat and breathing waveforms. Its other requested channels include heart rate, mean/systolic/diastolic arterial pressure, cardiac output, blood volume, respiration rate, tidal volume, oxygen saturation, arterial gases and pH, renal filtration/flow/urine production, insulin synthesis, core/skin temperature, metabolic power and sweat. The engine adds a `CTSresistance` diagnostic column; this is preserved separately from the 22 requested physiological channels.

Execution and finite outputs establish that the acquired original engine runs locally. They do **not** establish independent experimental validation, subject-specific calibration, adequate behavior under interventions or correctness of every internal coupling.

A first run requested the obsolete `BloodPH` channel and retained the tracker's default 1 Hz cadence. BioGears logged `FATAL Unknown Data Request : BloodPH`, emitted `-1.$` for that channel, and still returned zero. That failed attempt is retained in `biogears_native_run_failed_request/` and is excluded from passing results. The corrected run uses the source-defined `ArterialBloodPH` identifier and 50 Hz sampling. No invalid values were repaired, deleted or imputed.

## Completed run evidence

The corrected run exited zero after 41.27 wall seconds (including initialization), produced 3,000 samples from 0.02 to 60.00 seconds, and retained all 24 columns: time, 22 requested physiological channels and the engine diagnostic. Every value is finite and time strictly increases by 0.02 seconds. No fatal or unknown-request message appears. The log reports primary convergence of 402 simulated seconds, secondary convergence of 40 seconds, then approximately 30 seconds of final homeostasis before resetting simulation time to zero.

Arterial pressure spans 74.3238–114.7222 mmHg; total lung volume spans 2595.8481–3132.6704 mL; arterial pH spans 7.416431–7.417852. A sampling check detected 71 arterial-pressure peaks (median period 0.84 s, 42 samples/cycle) and 16 lung-volume peaks (median period 3.62 s, 181 samples/cycle). These checks show the beat and breathing waveforms are resolved; they are not clinical validation. Exact extrema for every channel and units in the native column names are preserved in `summary.json`; peak method details are in `waveform_checks.json`.

`native_build_manifest.json` records package versions, compiler/architecture, original release URL, source revision, and hashes of libraries, executable, adapter, original trajectory and engine state files (`states/native_stabilized.xml`, `states/native_final.xml`). The build script was re-executed successfully against the existing build, and all 8,850 acquired source-file hashes still pass verification.

## Configurable predictor API

`ihm.native.NativeConfig` validates an upstream patient identifier, optional existing XML state, duration (0.02–3600 seconds on the native 0.02-second grid), sample frequency (1, 2, 5, 10, 25 or 50 Hz), and up to 100 ordered interventions. `NativeConfig.from_dict` rejects unknown fields. `run_native(config, fresh_output_directory)` returns a checked summary, and `load_trajectory(directory_or_csv)` returns original columns, time and numeric channel arrays. State selection is intended for trusted local callers; the HTTP API must exclude arbitrary state paths. Existing execution directories are never overwritten.

```python
from ihm.native import NativeConfig, Intervention, run_native
config = NativeConfig(seconds=180, interventions=(
    Intervention(30, 'exercise', 0.3),
    Intervention(90, 'exercise', 0.0),
))
summary = run_native(config, 'data/derived/physiology/my_exercise_run')
```

Exercise uses upstream `SEExercise::SEGeneric.Intensity` (0–0.5); zero stops exercise. Hemorrhage uses upstream `SEHemorrhage`, a fixed `RightLeg` compartment and initial rate (0–100 mL/min); zero stops the hemorrhage. This is an initial flow used to derive the upstream bleeding resistance, so subsequent actual flow is pressure dependent. Saline uses upstream `SESubstanceCompoundInfusion`, 0–100 mL/min and a fresh 500 mL bag at each action. Zero stops the infusion. These are simulator workload limits, not clinical safety recommendations. Source examples are `projects/howto/Exercise/src/HowTo-Exercise.cpp` and `projects/howto/common-source/HowTo-ThreadedBioGears.cpp`. No physiological equations were copied or modified.

`--config path.json` on the runner accepts the same configuration fields with JSON intervention objects (`time_s`, `kind`, `value`). `--state path.xml`, `--patient`, and `--sample-hz` support direct command-line use. The engine stores complete initial and final states in `states/`. The summary records source SHA, executable/adapter/CSV/state SHA-256 hashes, explicit quantity/unit/system bindings, simulation evidence labels and unknown parameter uncertainty. Diagnostics have no asserted physiological binding. An ensemble confidence interval is not inferred from a deterministic run.

The adapter can be rebuilt quickly against existing upstream libraries with `.venv/bin/python scripts/compile_native_adapter.py`. `.venv/bin/python scripts/verify_native.py` verifies contract validation and CSV rejection. `.venv/bin/python scripts/native_experiments.py` materializes exercise/recovery and hemorrhage/saline scenarios from the same native stabilized state as `native_baseline_v2`.

An explicit case-insensitive source search for `supine|posture` in `projects/biogears/libBiogears/src` found no matches. The adapter therefore asserts no posture setting; any hydrostatic display adjunct must be described separately and cannot be treated as an upstream body-posture action.

## Paired experiment results and failure evidence

The new `native_baseline_v2`, `native_state_baseline`, `native_exercise_low` and `native_fluid_recovery` runs each produced 9,000 finite samples over 180 seconds at 50 Hz. A separate 2-second loaded-state run verified 10 Hz sampling (20 rows). `native_experiment_verification.json` preserves numerical assertions and all state-reload channel discrepancies. The baseline blood-volume range is 0.0764% of initial volume; this measures stability, not whole-body mass conservation.

State reload is **not bitwise identical**: maximum blood-volume discrepancy is 0.072754 mL and arterial CO2 discrepancy 0.065009 mmHg, while instantaneous lung-volume discrepancy reaches 118.243 mL (waveform timing). Accordingly perturbations are compared to the reloaded control. A 0.05 intensity exercise from 30 to 60 seconds increased late-exercise mean heart rate by 1.827 bpm and metabolic power by 0.856 W; late recovery power excess is 0.088 W. CO2 remains elevated, so this is partial recovery. The 60 mL/min initial hemorrhage from 30–90 seconds caused a 44.831 mL mean volume deficit near its end; saline at 60 mL/min from 90–150 seconds increased the corresponding volume difference by 55.811 mL. Pressure-dependent bleeding and native fluid exchanges explain why the nominal hemorrhage-rate integral is not an exact lost-volume measurement.

The first 0.3-intensity exercise trial, stopped at 90 seconds, entered irreversible hypercapnia at 130.32 seconds and returned exit 3. `native_exercise_recovery` retains that failed experiment; it is excluded from passing summaries. This is a material limitation of the source model/scenario, not repaired data. Bounds on accepted actions do not certify physiologically safe responses. Complete clinical validation, exact state restart equivalence and whole-body conservation remain unestablished.

### Fractional-time scheduling regression

A browser scenario requested 2 seconds at 10 Hz, with exercise actions at 0.4 and 1.4 seconds. The previous adapter called the upstream duration overload separately for each interval; floating-point division/truncation yielded only 1.98 seconds and the strict final-time check rejected it. That failed run is preserved under `data/derived/scenarios/20260904-231755-3a378e1e/output`.

The adapter now converts validated durations/action times to integer 0.02-second steps and calls the official single-step API exactly that many times. Automatic tracker scheduling is disabled, and original `TrackData` is called at integer sampling strides, avoiding floating tracker-clock offsets. Engine equations and tracker values remain unchanged. Final engine time is checked independently of the sampled horizon: a 2.04-second run at 10 Hz legitimately has 20 samples ending at 2 seconds, while its engine state must reach 2.04 seconds.

`.venv/bin/python scripts/verify_native.py --engine-clock` exercises fractional action intervals at all six allowed sample rates using the preserved stabilized state. All rates pass exact count/action/final-time checks. The original browser configuration also passes from fresh patient initialization in `native_fractional_clock_original` (20 rows, final engine time 2 seconds).

### Explicit thermal context and hourly drift

`NativeConfig(ambient_temperature_c=22, clothing_clo=1.5)` applies the original
`SEEnvironmentChange` action after initialization/state loading and before the first
recorded step. Omitted fields preserve the source environment. Ambient override sets
ambient, mean radiant, and respiration ambient temperatures together; humidity,
air velocity, pressure, gas fractions and emissivity remain source-defined. Allowed
ranges (10–35°C, 0–3 clo) are adapter experiment bounds, not clinical safety limits.
Original source and requested conditions are recorded in the runner log and summary.

The original StandardEnvironment is 22°C, 0.5 clo, 0.1 m/s air velocity. Its dynamic
stabilization tests core temperature over only a 15 s convergence window. This does
not establish an hour-long steady state. A matched stabilized-state 3600 s control
cooled from 36.561 to 33.817°C; increasing insulation to 1.5 clo at the same room
temperature ended at 35.638°C. A declared 26°C/1 clo comparison ended at 35.706°C. All three runs contained
3600 finite samples and reached exactly 3600 s; none establishes thermal equilibrium.
Blood volume rose approximately 229 mL in both 22°C runs.
The source stomach initially contains 500 mL water plus electrolytes; its absorption
and redistribution explain why resting blood volume is not fixed. The source-model
thermal drift remains independent of this input and is not a validated bed-rest prediction.
The preserved original independently initialized hour run is `native_hour_rest`;
matched experiments are `native_hour_thermal_*`. Run the declared comparisons with
`scripts/native_thermal_experiments.py --case control|insulated|warm` into fresh outputs.
No parameter was fitted to a desired physiological trajectory.

Source `Environment.cpp` uses a standing effective radiation area factor (0.73).
There is no native posture action; an external supine hydrostatic adjunct does not
modify that thermal assumption. Five original heat-loss channels are now exported.
A confirmed upstream telemetry defect assigns the convective accumulator to
`EvaporativeHeatLoss`, instead of the computed evaporation accumulator. Raw values
are retained with `usable_for_energy_accounting=false` in that channel's binding.
Do not sum these channels as an independent energy balance. The original thermal
equations and telemetry source are unmodified.

### Reproducible bounds and heat telemetry source variants

The original StandardFemale initialization abort is preserved in
`data/derived/scenarios/20260904-235212-40e8d665` and the no-action reproduction
`native_female_init_repro`. A debug/AddressSanitizer build identifies
`Saturation.cpp:536`: the negative-result correction reads/writes `x(3)` although
`Eigen::VectorXd x(3)` has only indices 0–2. The preceding blocks already handle all
three valid entries. `scripts/build_biogears_saturation_variant.py` removes only
that redundant invalid block in a separate source copy and shared library.
Original source and library are untouched. Explicit
`NativeConfig(engine_variant='saturation_bounds')` selects this variant; the default
remains `upstream`. Corrected Female rest and the previously failing 2 s fluid
protocol pass; matched Male original/corrected output values are bitwise identical.
These checks establish the regression correction, not clinical patient validation.

`--heatflux` builds `saturation_bounds_heatflux`, additionally assigning the computed
`eHeatLoss_W` to the evaporation telemetry scalar. The source scalar has no
physiological equation consumer. Matched 10 s original-telemetry/corrected-telemetry
runs have identical time and every other output column. Only this explicit variant
marks evaporation suitable for heat accounting; do not double-count aggregate skin
heat loss. Each variant directory under `data/runtime/physiology/variants` contains
exact unified patches, original/patched source hashes, and the selected library hash.
The runner checks library integrity and stores the variant manifest in each summary.
Run `scripts/verify_biogears_saturation_variant.py` to verify the preserved regressions.
