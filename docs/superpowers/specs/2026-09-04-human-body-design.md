# Implicit human body model

Adapt IBM-1's four primitives: fields are uncertain physical quantities on named
supports; anatomy supplies overlapping partitions; topologies permit interactions;
processes supply dynamics. Evidence and materialization are operations on these
primitives. A separate package avoids importing IBM's brain-specific registry.
Copying the entire brain implementation would also copy head geometry assumptions;
a single universal latent predictor would discard the physical-state contract.

The first executable representation is organ/compartment scale, scalar Gaussian
state with full cross-component covariance, continuous-time affine dynamics and
exact linear Gaussian transitions. Named predictor requests trace all upstream
process dependencies. A hard state budget bounds dense allocation. Geometry
refinement, spectral state and nonlinear physiology are explicitly future work,
not silently accepted request settings. The ontology covers cardiovascular,
respiratory, metabolic, endocrine, renal, immune, thermal, musculoskeletal,
peripheral neural and device state. All default parameters are illustrative weak
priors, not claimed population norms or calibrated physiology.

CSV tables, JSON records and NPZ time series normalize through explicit source
cards into measurements with subject, time in seconds, units, variance and stable
record identity. Missing values are absent evidence. Unknown units, components,
subjects, nonfinite values and duplicate evidence fail. Simultaneous correlated
measurements use a full error covariance. Different subjects never fuse. Source
cards distinguish measured, derived and synthetic evidence; catalog-only external
sources never claim to have been ingested. Raw images and omics need a declared
measurement operator before they can constrain physical state.

Materializations support filtering, forecasting, hard state clamps, provenance,
and portable NPZ/JSON export without pickle. Regression fitting provides an
alternative affine process form with source and training-subject provenance;
held-out evaluation must exclude training subjects. Synthetic fixtures demonstrate
the complete workflow; they do not establish biological predictive validity.

Verification covers Gaussian conditioning against closed form, cross-source unit
conversion, isolation and deduplication, dependency closure, covariance validity,
forecast/clamp behavior, serialization, regression recovery and held-out leakage.


## Scope refinements from the user

Cover all body systems explicitly rather than treating the body as a uniform
extension of brain tissue. The inventory separates blood from vessel walls,
capillary exchange, interstitial fluid, lymph transport/nodes and hematopoiesis,
as well as the remaining organ systems. Declared states are not evidence of
implemented or calibrated physiology.

Integumentary voltage interactions are a first-class materialization, motivated
by the user's interest in Michael Levin. Separate cellular membrane and
extracellular networks support heterogeneous cells, conductances and geometry.
A conservative fluid network handles vascular/interstitial/lymphatic interfaces.

The user asks whether the body literally breathes and its heart beats, and
suggests lying flat in bed. Add a nonlinear supine cardiopulmonary view with a
four-chamber pulsatile circulation, pressure-operated valves and respiratory
mechanics coupled by pleural pressure. Defer movement. Fixed phase pacemaking,
uncalibrated parameters, no gas chemistry and no nonlinear assimilation must be
stated explicitly. Verify cycle generation, volume conservation and tidal airflow.
