# Native OpenSim mechanics

The official OpenSim C++ engine and Simbody are built natively for ARM64 in `data/runtime/opensim`. They execute the unchanged acquired Rajagopal2016 model: real source wrapping, musculotendon lengths, independent rotational moment arms, muscle equilibrium and forces. The app now displays native wrapped paths and source-state mechanics. No Python binding, visualizer, paid compute or system package installation is needed.

## Reproduce

```sh
.venv/bin/python scripts/build_native_opensim.py
.venv/bin/python scripts/build_opensim_wrap_variant.py --iterations 8 --tolerance .0005 --cache-fix
.venv/bin/python scripts/verify_native_opensim.py --engine --corrected
.venv/bin/python scripts/build_opensim_display.py
.venv/bin/python scripts/verify_opensim.py
```

Build stages `dependencies`, `engine`, and `adapter` support incremental reproduction. OpenSim source revision is `86b30588374650fbaf012a345a836a64f6855522`; Simbody is the exact revision required by that OpenSim dependency recipe, `59c6e7b89b3bdf266a2f3d54c611599d964205f7`. Both source repositories stay clean under `data/raw/mechanics`. Ubuntu ARM64 LAPACK/BLAS/Fortran dependency packages are downloaded and extracted locally, with hashes in `build_manifest.json`. Libraries, commands, build logs, original source identity and adapter identity are retained. API-only installation needs an empty `bin` directory because the upstream installed CMake configuration requires it. The adapter RPATH includes the extracted Fortran runtime directory. Rebuild after moving the workspace.

The original installed engine is preserved. An **explicit, separate cache correction** is used for app mechanics. It is not silently described as untouched upstream execution. Its library resides in `variants/wrap_8_0.0005_cache`, with a reproducible unified patch, original/patched file hashes and resulting library hash. The original installation remains separately executable.

```python
from ihm.native.opensim_backend import OpenSimConfig, run_opensim
result = run_opensim(OpenSimConfig(
    coordinate='knee_angle_r', delta_rad=0.1, activation=0.05,
    engine_variant='wrap_8_0.0005_cache'),
    'data/derived/opensim/my_knee_prediction')
```

`delta_rad` offsets the selected source default, bounded to ±0.15 rad; activation is uniformly set explicitly in [0.01, 0.5]. The native wrapper rejects source range violations, dependent coordinates and an assembled coordinate that fails the requested tolerance. Locked coordinates are not unlocked. Source knee flexion has a lower bound of zero; negative requested knee offsets are invalid. The default variant is `upstream` for auditing; selecting the corrected variant is explicit. Output directories cannot overwrite prior executions.

## Upstream cache defect and isolated correction

The first native derivative check found knee moment-arm discrepancies of 16.014 mm for `gasmed_r` and 10.860 mm for `gaslat_r`. Central finite differences at 0.1 rad with steps 1e-4 and 1e-5 reproduced them. Three identical-state cache rebuilds reproduced the same lengths and moment arms, ruling out mere call ordering or one stale outer-state evaluation. A separate experiment raised the multiwrap iteration cap from 8 to 200 and tightened length convergence from 0.5 mm to 1e-10 m; it did not resolve the discrepancy. Those outputs remain in `native`, `cache_diagnostics`, `refined` and `native_initial_checks`.

The cause was then identified directly: `PathWrapPoint::setLocation` updates the wrap-local coordinate cache without invalidating its inherited `Point` caches. During iterative wrapping at the same simulation stage, `Point::getLocationInGround` can therefore reuse a previous tangent location. Original native cached Ground positions differed from the current local coordinate transformed through its body by **6.03 mm**. The isolated patch invalidates exactly the inherited `location`, `velocity` and `acceleration` caches when wrap-local location changes. No geometry, force law, model parameter, wrapping algorithm or default convergence threshold changes in this corrected variant.

With that correction, the cache discrepancy drops below 1.4e-17 m; the gastrocnemius derivative discrepancies fall to 0.01054 mm and 0.00428 mm. All 80 right-knee moment arms agree with finite differences within 0.052 mm. The patch is `data/runtime/opensim/variants/wrap_8_0.0005_cache/cache_fix.patch`; `scripts/build_opensim_wrap_variant.py` compiles only the patched translation units and relinks a separate library from the original build objects. Original source and libraries are untouched. No geodesic approximation or replacement cylinder equations were introduced.

## Numerical evidence and scope

`data/derived/opensim/native_corrected` contains default pose, 0.1-rad knee flexion, centered finite-difference perturbations and doubled activation. All 22 body transforms match the restricted XML evaluator within 3.33e-16; there are 80 muscles and 30 active wrap segments in the default pose. Maximum baseline tendon-versus-fiber-force equilibrium residual is 6.09e-7 N. Increasing activation from 0.05 to 0.1 increases total tendon force; knee flexion changes source musculotendon lengths. Additional coordinate checks verify left knee, right hip flexion/adduction and right ankle finite differences, with respective maximum discrepancies 0.052, 0.346, 0.349 and 0.0092 mm. Native wrapping has finite numerical precision; these residuals are retained rather than declared zero. A locked subtalar coordinate trial was rejected and its log preserved.

The source model is equilibrated at fixed coordinates and zero coordinate speeds. These are muscle-equilibrium forces at explicit activation, **not** an externally balanced standing pose, measured human forces or a subject-specific calibration. No force covariance or clinical validity is inferred. Moment arms have units m (equivalently m/rad for dimensionless radians); translational coordinates and dependent/locked coordinates are excluded. Current tests establish local derivative consistency at stated poses, not accuracy over all possible motions.

The native output preserves actual body transforms, coordinates, lengths, equilibrium forces, moment arms, source path nodes, wrapping tessellation and diagnostic cache checks. Geometry follows the original `GeometryPath::generateDecorations` traversal. The app retains 81 bones and 80 muscles under `opensim-rajagopal`; muscle geometry contains paired vertices for native wrapped polylines, `native_mechanics`, explicit engine variant and activation. Original straight attachment chords remain in `*-attachment-fallback.json.gz`. Source physiological parameters are distinct from computed forces. These mechanics remain a separate source specimen from BodyParts3D and BioGears.
