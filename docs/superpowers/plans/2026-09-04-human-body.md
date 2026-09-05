# Human body implementation plan

Goal: executable IBM-style body substrate with heterogeneous evidence and explicit predictors.
Architecture: four registered primitives, Gaussian evidence conditioning, dependency-traced affine dynamics.
Tech stack: Python 3.11+, NumPy, SciPy; JSON/CSV/NPZ artifacts.
Spec: ../specs/2026-09-04-human-body-design.md

- [x] Define a verification script for inference, ingestion, tracing and export before implementation.
- [x] Implement immutable field/support, anatomical, topology and process declarations plus registry validation.
- [x] Implement body ontology and six named predictor requests with bounded dependency closure.
- [x] Implement normalized source cards and CSV/JSON/NPZ ingestion, unit conversion and evidence validation.
- [x] Implement full-covariance filtering, exact affine dynamics, forecasts, clamps and portable export.
- [x] Implement source-attributed affine fitting and held-out subject evaluation.
- [x] Add synthetic multiformat evidence and verified external source inventory, CLI, examples and documentation.
- [x] Execute numerical verification and every CLI workflow; review scientific limitations and package build.

- [x] User refinement: heterogeneous per-cell skin electrical networks and independent extracellular graph.
- [x] User refinement: all-system inventory with distinct blood, vascular, interstitial and lymphatic supports.
- [x] User refinement: conservative explicit vascular/interstitial/lymph drainage network.
- [x] User refinement: nonlinear supine heartbeats and tidal breaths, with numerical checks and plots.
- [x] Review corrections: scale-aware covariance, preserved source manifests, frozen anatomical maps, assembled conservation checks, independent held-out fixture and ventricular elastance limits.
