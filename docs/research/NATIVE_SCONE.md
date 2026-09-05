# Native reflex locomotion acceptance

SCONE source revision `39c4b7a52a90dff838ebe7a1cb8afe58af1aaae1` is acquired from the [primary repository](https://github.com/tgeijten/scone-core/tree/39c4b7a52a90dff838ebe7a1cb8afe58af1aaae1). The shipped H0914 planar model has nine degrees of freedom and fourteen muscles. Its `ControllerGH2010v9.scone` implements Geyer/Herr-derived reflex control, with shipped optimized values in `ResultH0914Gait10.par`. These values have not been replaced by guessed gains or claimed as newly calibrated human measurements.

## Reproduction

Run in repository root:

```sh
.venv/bin/python scripts/acquire_scone_controllers.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_native_scone.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_scone.py
```

The build requires the separately held native OpenSim4/Simbody installation. Acquisition retains all 1,036 tracked donor file hashes and does not rewrite an existing valid acquisition timestamp. Pinned xo/spot gitlinks and their file hashes are retained separately. Modified source/submodule bytes fail verification. `data/sources/scone-predictive-locomotion.json` indexes receipts; raw sources, builds, failed build logs and native trajectories remain local under ignored data directories.

Build compatibility changes occur only in `data/runtime/scone/source-opensim4`: explicit OpenSim headers, updated base-frame contact accessors, rejection of unsupported offset contact frames, removal of obsolete Manager halt from the fixed-control-step TimeStepper path, tclap constructor compatibility, the missing xo version-string return, and an isolated settings-directory override. `compatibility_patches.json` retains exact diffs and before/after hashes. The build refreshes donor bytes before applying these changes and rejects unexpected extra variant files. Its manifest records every variant file, linked library, executable, commands and upstream OpenSim receipt. Runner checks current variant, donor, submodule, patch and linked-library hashes; it hashes executed scenario/parameter/settings files and trajectory output.

## Causal intervention

Source `controllers/Controller.cpp` reads `disabled` and `stop_time`. `Controller.h::IsActive` gates control updates on these values. Disabling from time zero changes source initialization's activation equilibrium, so the accepted experiment sets `stop_time = 1` in a copied gait controller. Physical state matches exactly at initialization and all samples before one second; subsequent actuator input is zero. Source minimum excitation/activation and passive muscle mechanics remain. The original controller, subject and parameter files are unchanged.

Native GaitMeasure terminates when center-of-mass height drops below 0.85 of initial height (the shipped measure's setting). Controller removal reaches that condition at 1.29 s. This is a source termination event, not simulation of the body striking the ground. A normal runner call rejects early termination; the intervention deliberately uses `--allow-early-termination`, and its receipt still says the requested duration was not completed.

## Numerical checks and limits

The intact run completes 30 s with 3,001 rows and 581 finite channels, traveling 23.7148 m. All fourteen activation and excitation outputs remain in [0.01, 1]; fiber lengths stay positive and muscle force nonnegative. Ground-reaction vertical components are unilateral and nonzero. The controller-off trial travels 1.36176 m before termination.

Accuracy 0.001 gives 23.5962 m over 30 s, a 0.50% displacement difference from source accuracy 0.002. Successive first-second runs at 0.002, 0.001 and 0.0005 differ by at most 0.35 mm translation and 0.00434 rad joint angle. Predeclared engineering acceptance bounds are 2 cm translation, 0.1 rad angular coordinate difference over one second, and 5% 30-second distance difference. These are regression tolerances, not human validation thresholds, and do not establish long-horizon pointwise convergence. Control sampling remains the source 5 ms; output sampling is 10 ms. Refinement changes integrator accuracy only.

The verifier also rejects nonfinite duration, malformed/nonfinite/backward-time trajectory data, and unapproved early completion. It checks native/header time agreement and logs actual native step counts. This validates execution and a causal role of the shipped reflexes in this reduced plant. It does not establish canonical 3D walking, a personalized gait, cortical intention, whole-body energy closure, or human predictive accuracy. The native evaluation remains separate from the canonical body assembly.
