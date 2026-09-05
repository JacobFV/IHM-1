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
