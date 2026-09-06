# Opt-in lumbar muscle extension

`scripts/materialize_lumbar_muscle_variant.py` creates a separately hashed 98-muscle model and catalog by inserting the six retained Gait2392 lumbar force elements into the existing 92-muscle model. Removing the exact retained `insert.xml` bytes reproduces the original `.osim` byte-for-byte. The original 22 bodies, inertias, joints, original muscle order, parameters and paths are therefore preserved. No default model pointer is changed.

The six additions are bilateral erector spinae, internal oblique and external oblique. Their complete Thelen source properties and unit-scale homologous back-joint registration come from `LUMBAR_SHOULDER_MUSCLE_COVERAGE.md` and its hashed audit. These remain transferred population priors; successful native loading does not validate anatomical insertions or subject strength. No scapula/clavicle body or shoulder muscle is introduced.

Run the source-only tests:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_lumbar_muscle_variant.py
```

With the shared native build slot explicitly allocated, `scripts/accept_lumbar_muscle_variant.py` compiles one isolated probe against the retained OpenSim libraries. It keeps the existing 80 FunctionBasedPaths, loads all 98 muscles, and tests neutral plus each lumbar coordinate independently at −0.1 and +0.1 rad. The probe initializes equilibrium at unchanged engine-default activation, records that activation, and advances no time. This is a constitutive/geometry test, not a force-balanced body pose.

The receipt contains 126 native moment-arm observations, central finite differences of path length, path lengths, tendon/active/passive fiber forces, pennation and fiber velocity. Neutral native values are compared with the independent source-geometry audit. Every input source, compiled source and library is pinned, and checked again after the probe. Native tendon tension must equal projected active-plus-passive fiber tension within explicit numerical tolerance.

The resulting `variant/registration.json` is accepted by the existing `augmented_registration` model/catalog interface. Its six appended catalog entries supply exact muscle IDs, source law and parameters, side, pelvis/torso attachment stations, and explicit cortical assignment priors. The `native_control_ready` field remains false: a source label plus successful force evaluation does not establish measured recruitment, spinal motor-neuron ownership, or brain-driven dynamics. No excitation is selected to settle a support solver, and reserve actuators remain separately identified. The native receipt is separate from the immutable source materialization's `native_acceptance_complete: false`; callers must inspect the acceptance receipt rather than treating materialization itself as certification.

The next controller integration must retain the first 92 control indices, append these six IDs, preserve signed muscle energy/metabolic ownership, and verify actual 98-muscle unified runtime checkpoint/replay and motor-port control. This fixture establishes native load and geometry/constitutive acceptance only; it does not claim runtime brain-controller acceptance or a feasible supine equilibrium.

## Actual native acceptance

The bounded probe passed in **5.889 s**, peak child RSS **1,060,752 KiB**, against the retained OpenSim libraries. Its 98-muscle/22-body load retained all 80 fitted lower-limb paths. All seven poses and 126 observations passed; maximum moment-arm versus finite-difference error was **1.066×10⁻¹⁰ m**, and maximum tendon versus projected fiber force mismatch was **1.422×10⁻¹⁴ N**. The native default activation was **0.05** for every added muscle. No patient time advanced.

Committed evidence is in `data/research/lumbar_native_acceptance/`; `retained_files.json` pins the materialized model/catalog/registration and build logs. The actual selectable registration is `data/derived/lumbar-muscle-native-lb45uirs/variant/registration.json`. The earlier failed compile is retained at `data/derived/lumbar-muscle-native-v4f2tup7/compile.log`; the probe passed after replacing a raw C-string argument with explicit `Set<FunctionBasedPath>` construction. No shared adapter source changed.
