# Whole-body causal simulation expansion

The user's breathing, eating, reflex and walking examples establish a realism bar, not a feature boundary. Continue the existing generic-body project toward coupled neural, mechanical, vascular, metabolic, endocrine, renal, thermal and gastrointestinal state. Preserve the existing anatomical/evidence substrate, working UI and prior numerical artifacts. Work is authorized on the current feature branch; no additional approval ceremony is required for reversible implementation.

## Architecture

Use the retained native BioGears and OpenSim engines as causal state owners where their equations exist, with explicit adapters rather than silently duplicating their state in Python. Add a persistent native stepping protocol so a coordinator can read physiological state, deliver interventions, advance mechanics/controllers and exchange bounded inputs during a run. Keep native chemoreceptor and autonomic feedback under native ownership until a specified replacement is implemented and verified. A neural population label is not evidence that its output drives a muscle.

The native nervous system already computes chemical respiratory feedback. Expose requested driver pressure/frequency and actual applied circuit pressure, airway flow, phase, blood gases and autonomic scales. Apnea severity provides a real intervention on active respiratory pressure. The native default pressure is retained during drive suppression; source timing/phase latching and checkpoint limitations must be tested. Full left/right phrenic, intercostal and abdominal actuation requires an explicit partition beyond this one-actuator source model.

Meals enter the source stomach as a named composition and pass through the source digestion, chyme, vascular, endocrine, hepatic, tissue and renal processes. Retain native compartment quantities and source fluxes separately: differences in vascular storage are not absorption flux because transport/metabolism also act. Confirmed source defects receive isolated, checksummed library variants and original-failure evidence. The observed calcium g/mg defect is corrected; ambiguous macronutrient rate naming does not justify blindly dividing empirically tuned rates by60. Missing colon/fiber/stool, mouth/esophagus, enzyme physiology and intestinal segmentation remain explicit until added with owned state and conservation.

Native OpenSim time integration replaces pose-only evaluation for locomotor experiments. Retained 3D walking resources provide an articulated 22-body/80-muscle model, observed kinematics/GRF/EMG and contact definitions. Distinguish observed tracking, computed muscle controls, endogenous ground contact, reflex feedback and autonomous gait. Acquire/build required native optimization or controller dependencies when needed. Do not promote the source's explicitly non-research 2D example to a validated human gait model.

Canonical soft-body coupling must retain anatomical attachment transforms, finite strains, equal/opposite internal reactions, external load/work and inferred-material provenance. Whole-body all-part contact is not established by an affine pose or independent reference patch. Progress from audited nonuniform multi-domain contact to anatomical material partitions and native articulated boundary motion. Any model reduction records its support and error; higher source mesh resolution remains retained.

## Required contracts

- SI for new cross-engine quantities; native units are explicit at adapters and converted once.
- Native dt remains0.02s; commands land on integer steps. Other engines subcycle with declared exchange latency.
- One owner per blood, gas, solute, energy, mechanical and neural state. No duplicate native storage or recoil.
- State-dependent controllers read computed state; event schedules and recorded tracking targets are labeled as external inputs.
- Each experiment records initial state, source/runtime/library hashes, action timeline, clocks, solver settings, outputs, numerical audits and predictive limits.
- A named body-system interaction requires an executing edge and a perturbation check; file presence, geometry or labels alone do not qualify.
- Preserve original native source checkouts and prior outputs. No brain microvascular/microcircuit acquisition here; IBM imports remain versioned.
- Use `.venv/bin/python` and one OpenBLAS thread for reproducible host execution. Keep the user's live workspace and data corpus intact.

## Acceptance and continued scope

Required causal experiments include drive suppression/restoration with blood-gas feedback, meal versus matched fasting with vascular/hormonal/storage outputs and nonnegative nutrient contents, reflex blockade/delay with mechanically measured response, native articulated/contact perturbations, and conservative heterogeneous soft-tissue response. Add combined perturbations to detect unintended isolation between systems. Numerical accuracy, empirical parameter support and anatomical completeness are distinct gates. Extend coverage to autonomic, immune/inflammatory, reproductive, special-sense, marrow/adipose/connective and excretory interactions through the same process rather than treating the examples above as the complete target.

The UI must expose actual engine state and causal chains on one body, with regional detail as a materialization of that substrate. Display must show breathing/heartbeat/motion from computed state, and food/substrate/nerve information from recorded owned quantities. Time scales from milliseconds to hours remain distinct; finite-window Laplace descriptors are not mistaken for mechanistic transfer functions or calibrated uncertainty.
