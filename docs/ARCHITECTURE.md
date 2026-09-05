# Architecture

IHM-1 adapts the architecture and evidence/provenance discipline of the local
IBM-1 project. It does not copy brain-specific geometry, cortical operators or
assume that body tissue has one common interaction topology.

```
ihm/
  registry.py                    validated namespace
  body.py                        body ontology assembly
  fields/                        components, supports, system inventory
  anatomy/                       independently overlapping memberships
  topologies/                    permitted interactions
  processes/                     affine declarations; cardiac/respiratory physics
  materialize/                   symbolic tracing and explicit predictor views
  runtime/                       Gaussian transition and covariance arithmetic
  forge/                         source binding, ingestion and process fitting
```

## Representations

The implicit ontology exists independently of any one discretization. Organ-scale
components offer a coarse initial inventory; explicit skin and fluid geometries
add states on finer supports. Their parameters are specified per cell,
compartment or edge. Distinct topologies remain distinct even when their spatial
regions overlap.

`body(skin=patch, fluids=network, cardiopulmonary=parameters)` assembles compatible
declarations into one registry. Registration rejects duplicate identifiers,
unknown supports/regions, invalid priors, dangling process selectors, illegal
edges and invalid anatomical memberships. Sealed anatomical maps are immutable.

`materialize(registry, Request(targets, subject, max_states))` computes upstream
closure including feedback cycles, then allocates only selected components and
processes. An explicit state budget bounds dense allocations. Scalar Gaussian
models use full cross-component covariance. These models have no spatial/spectral
refinement knobs: unsupported resolutions are not silently accepted.

The skin and conservative fluid materializations use the same affine runtime.
Cardiopulmonary cycles require a nonlinear executor; generic affine materialization
rejects those processes rather than treating nonlinear declarations as zero drift.
The nonlinear view retains the same registered state/process provenance but has
no quantified uncertainty or measurement assimilation yet.

## Evidence is not a force

Process pressures sum in the derivative. Observations instead condition the
state distribution through a measurement likelihood. Gaussian conditioning uses
a full measurement covariance and Joseph-form covariance update. Missing data is
absence of a constraint. Finite time, units, subject identity, source identity and
measurement error are required. An exported model includes the source manifest,
input file hashes, evidence history and chosen process forms.

Source cards distinguish synthetic, measured and derived evidence. Derived data
needs calibrated error and dependence accounting; it is not an independent
observation merely because its row has a different identifier. Cross-batch
measurement error dependence is outside the current implementation.

## Dynamics and uncertainty

Affine processes assemble `dx = (A x + b) dt + L dW`, `Q = L L^T`. Matrix
exponentials propagate both drift and diffusion. Repeated squaring avoids a
large unstable block exponential on long stable intervals. Covariance validation
normalizes by component scales before checking positive semidefiniteness.

Hard clamps replace the selected derivative rows and initial variance over an
interval; they are recorded as interventions. Forecasts use a copy of the current
posterior. Persistence uses numeric NPZ arrays and JSON metadata without pickle.

Affine process fitting requires physical inputs and process-specific derivatives.
It preserves input/output/topology declarations and records source/training
subjects. It does not infer a posterior over model-form/parameter uncertainty.
Prior-only and fitted process labels remain visible; biological validation remains
false until external performance evidence actually exists.

## Current limits

The inventory spans major body systems, but a field declaration is not a complete
physiological model. Many organ mechanisms are only illustrative local affine
couplings. Conservation is specifically enforced in the specialized fluid and
paired electrical networks, not claimed for the entire organ scaffold.

Dense covariance and explicit matrices suit small materializations. Whole-body
cell-resolution state, spectral beliefs, nonlinear electrodiffusion, automatic
image/omics operators, automatic species transfer and whole-body closed-loop
homeostasis remain unimplemented. Finer skin, drainage and cardiac models are
alternative views; they are not automatically coupled to every coarse component.
