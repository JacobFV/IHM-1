# Coherent OpenSim default-pose geometry

`ihm.spatial.opensim` evaluates the acquired Rajagopal2016 XML default pose in its native meter/radian coordinates. `scripts/build_opensim_display.py` adds the separate `opensim-rajagopal` model to the local app: 81 bone meshes, 80 muscle attachment paths and 288 original attachment/via points. All 22 body frames resolve through 22 joints to Ground. The source model SHA-256 and repository revision remain those in `opensim__Rajagopal__Rajagopal2016.json`; raw assets are unchanged.

Run from the project root:

```sh
.venv/bin/python scripts/verify_opensim.py
.venv/bin/python scripts/build_opensim_display.py
```

The builder looks first in the model-local Geometry directory, then the acquired repository's common Geometry directory, recording each selected file and its SHA-256. Mesh coordinates are multiplied by their explicit source `scale_factors`, transformed through the resolved body frame, rotated into the display frame and centered without resizing. The model spans approximately 1.669 m vertically. The source Ground frame, all body transforms, default coordinates and display rotation/translation are preserved in `data/derived/opensim/default_pose.json` and the app manifest. Source bone meshes are decimated only for display if they exceed 2,500 triangles. Every structure keeps its source identity.

The exporter replaces only this model's manifest entries, uses an atomic manifest replacement and produces deterministic gzip payloads. Re-running does not duplicate structures. The model is a separate source specimen, not a registration to BodyParts3D or the BioGears patient. Its source default pose is standing, not a supine posture assumption.

## Exact supported evaluation

An offset frame uses the source body-fixed XYZ rotation convention. With `X_GB` mapping body coordinates to Ground, each body transform is `X_GP X_PF X_FM inverse(X_BM)`. The custom joint's first three axis rotations multiply in source order; translations are the sum along its three specified fixed translation axes. Axes are normalized as in Simbody. These conventions were checked against official [OpenSim OffsetFrame implementation](https://github.com/opensim-org/opensim-core/blob/main/OpenSim/Simulation/Model/OffsetFrame.h) and [Simbody FunctionBased implementation](https://github.com/simbody/simbody/blob/master/Simbody/src/MobilizedBodyImpl.h); downloaded reference files and hashes are retained under `data/derived/opensim/reference/`.

Supported functions are Constant, LinearFunction and **exact-knot SimmSpline evaluation**. Every one of the 14 SimmSpline evaluations in this model lands exactly on an original spline knot at default coordinate zero. No spline type, endpoint derivative or extrapolation is guessed. Off-knot spline evaluation fails explicitly. PinJoint uses its native Z rotation. UniversalJoint is supported only at its zero-coordinate identity transform, which covers both wrists in this source. The two enforced linear coordinate couplers are checked against default values; a violated or unsupported constraint is rejected instead of silently assembling a guessed pose. Missing files, unresolved frames, unsupported points/functions and incomplete body-frame closure fail explicitly.

The acquired model contains only fixed `PathPoint` muscle points: no moving or conditional points require approximation. All points are transformed from their declared parent body, preserving both original local coordinates and computed Ground coordinates. This evaluator is restricted to the default configuration and is **not** an installed OpenSim engine, generalized motion solver, forward dynamics model or muscle-force solver.

## Muscle path limitations

There are 46 source path-wrap references. The display explicitly shows straight chords between attachment/via points and does not solve wrapping around the source wrap objects. Consequently lines can cross bone surfaces; chord length is not presented as musculotendon length. The payload labels wrapping unsolved and leaves activation, force and moment arm null. Maximum isometric force, optimal fiber length and tendon slack length are source parameters, not newly measured or fitted values. Neither source-defined wrapping nor independent clinical calibration is claimed.

Tests verify known rigid transforms, a two-body frame chain, affine function evaluation, exact spline knots and off-knot rejection, complete real-model frame closure, the native pelvis height of 0.94 m, rotation orthogonality, all 80 muscles and 288 points. Generated-asset tests check all geometry coordinates and indices, units, provenance hashes and the 161-structure manifest contract.

## Native mechanics extension

When verified native output exists, the display builder now uses native wrapped paths instead of attachment chords, while archiving the chord fallback. `native_mechanics` exposes computed length, forces, rotational moment arms, activation, derivative checks and explicit engine variant on both structure metadata and geometry payloads. Bone frames were independently checked against the native engine. The exact cache correction, original failure evidence, build procedure and scope are documented in [NATIVE_OPENSIM.md](NATIVE_OPENSIM.md). The earlier restricted XML evaluator remains available as the source-preserving fallback and independent frame check.
