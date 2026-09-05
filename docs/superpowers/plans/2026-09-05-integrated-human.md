# Integrated Human Implementation Plan

> **For agentic workers:** Use superpowers:subagent-driven-development for independent domains and task reviews. User requests continuous autonomous execution.

**Goal:** Implement a source-grounded integrated human workbench with native physiology, physical interfaces, calibration, temporal predictors and real 3D visualization.
**Architecture:** Independent model families share explicit observation and physical exchange contracts. A local application exposes geometry, trajectories, spectra and evidence with provenance.
**Tech Stack:** Python/NumPy/SciPy, native BioGears C++, browser WebGL/Three.js; local HTTP API.
**Spec:** docs/superpowers/specs/2026-09-05-integrated-human.md

## Global constraints
- Never label illustrative or source-simulated coefficients as independently calibrated.
- Preserve raw assets and their hashes; write generated artifacts under data/derived or artifacts.
- No overlapping agent ownership; root integrates CLI/dependencies/API.
- No publishing, credential changes, paid compute or source-data replacement.
- Test numerical behavior, not implementation details; preserve failure evidence.

## Task 1: Native physiology materialization
Files: ihm/native/, scripts/native_biogears_rest.cpp, scripts/run_native_biogears.py, scripts/build_native_biogears.py, scripts/verify_native.py, docs/NATIVE_BACKEND.md.
- [ ] Test parameter bounds, scenario validation, request rejection and trajectory roundtrip.
- [ ] Add configurable upstream-native scenarios and patient/state selection; expose a Python native predictor API without copied equations.
- [ ] Bind native variables to physical quantities with explicit units; run baseline and perturbation/recovery scenarios; persist hashes and complete state.
- [ ] Verify finite outputs, temporal sampling, changes under perturbation and baseline conservation; inspect posture support from source.
Interface: native module exposes validated run configuration, run execution, trajectory loader; coordinate with root before shared API adoption.

## Task 2: Temporal evidence and explicit predictors
Files: ihm/temporal/, scripts/build_temporal_atlas.py, scripts/verify_temporal.py, docs/TEMPORAL_MODELS.md.
- [ ] Inspect IBM spectral basis and priors; test analytic exponential Laplace transform and known oscillatory frequency before implementation.
- [ ] Implement finite-time Laplace, normalized PSD/cross spectra, coherence and transfer resolvent with validation.
- [ ] Fit stable temporal predictors with chronological holdout, reconstruct/forecast artifacts and rank/conditioning diagnostics.
- [ ] Generate actual native-variable spectral atlas and clearly distinguish short-run simulated evidence from observed human recordings.
Interface: data/derived/temporal/index.json and per-run JSON outputs; share output schema with root/app.

## Task 3: Local 3D application
Files: app/ only (including package manifest and browser tests), docs/APP.md.
- [ ] Build local WebGL app from real geometry/API contract supplied by root; no invented anatomy.
- [ ] Add model-family/system filters, search/picking, opacity/clipping, source detail, real flow/vector playback and chart/spectral panels.
- [ ] Add bounded scenario submission and run status with explicit backend availability.
- [ ] Verify production build, keyboard/control behavior and browser rendering; document startup.
Interface: /api/manifest, /api/geometry/{id}, /api/physiology, /api/temporal, /api/evidence, /api/scenarios. Root owns API and generated assets; frontend must handle error/empty/loading states.

## Task 4: Spatial supports and conservative integration
Files: ihm/spatial/, ihm/coupling/, scripts/build_spatial_atlas.py, scripts/verify_spatial.py, scripts/verify_coupling.py.
- [ ] Test known-frame transforms, invalid registrations, constant-field mapping and conserved exchange cancellation.
- [ ] Implement frames, weighted landmark registration, residual/uncertainty metadata, conservative overlap maps and dimensional flux contracts.
- [ ] Build real atlas display meshes, OpenSim muscle path geometry in valid source frames and vascular flow frames with metadata.
- [ ] Implement executable vascular/interstitial/lymph and bioelectric exchange adapters; validate conservation and avoid double counting native circulation.

## Task 5: Calibration and system coverage
Files: ihm/calibration/, scripts/calibrate_skin.py, scripts/build_system_coverage.py, scripts/verify_calibration.py, docs/CALIBRATION.md.
- [ ] Test recovery of identifiable synthetic parameters and detection of rank-deficient parameter combinations.
- [ ] Implement bounded weighted fitting, Jacobian rank, covariance/uncertainty, source binding and separated holdout reports.
- [ ] Fit appropriate observation-level skin quantities to real clinical data, preserve excluded/missing groups, test holdout calibration.
- [ ] Link all major/minor system inventory to actual source mechanisms and measurement coverage; expose unresolved unknowns quantitatively.

## Task 6: API, integration and verification
Files: ihm/app/, ihm/cli.py, pyproject.toml, scripts/verify_app.py, README.md.
- [ ] Test HTTP traversal rejection, strict run arguments, asset manifest consistency and useful missing-backend errors.
- [ ] Implement local API, safe asset serving and bounded scenario worker; integrate CLI materializations and app launch.
- [ ] Run existing and new checks, native perturbation experiment, asset roundtrips, app production build and real browser flows.
- [ ] Review independent domains, fix findings, preserve reproducible commands and update actual remaining scientific constraints.

## Execution ledger
Ruling: Work in the existing project directory on a new feature branch after preserving a baseline commit; original project was not a Git repository. This avoids moving or duplicating gigabytes of source assets and compiled absolute paths.
Ruling: Independent domains use parallel agents per dispatching-parallel-agents; shared CLI/dependencies/API remain root-owned.
Ruling: 'All concerns' requires implementing mechanisms for uncertainty and unavailable evidence, not claiming impossible universal calibration. Unknown experimental quantities remain explicit and testable.
