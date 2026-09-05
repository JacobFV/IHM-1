# Whole-body causal expansion implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Make physiological control, nutrient transport, reflexes and articulated/contact mechanics execute as observable causal body processes, extending the wider system interaction program.

**Architecture:** A persistent native stepping boundary preserves BioGears physiology ownership and enables online co-simulation. Native anatomical mechanics and controller integration feed explicitly typed body ports; source integrity corrections, conservation audits and material-coordinate visualization remain separate concerns.

**Tech Stack:** Python/NumPy/SciPy, C++20, retained BioGears8/OpenSim/Simbody/FEBio, Three.js and Playwright.

**Spec:** `docs/superpowers/specs/2026-09-05-whole-body-causal-design.md`

## Global constraints

- SI for new cross-engine quantities; native units are explicit at adapters and converted once.
- Native dt remains0.02s; commands land on integer steps. Other engines subcycle with declared exchange latency.
- One owner per blood, gas, solute, energy, mechanical and neural state. No duplicate native storage or recoil.
- Preserve original native source checkouts and prior outputs. No brain microvascular/microcircuit acquisition here; IBM imports remain versioned.
- Use `.venv/bin/python` and one OpenBLAS thread for reproducible host execution. Keep the user's live workspace and data corpus intact.

### Task 1: Isolated gastrointestinal mass-integrity correction

**Files:** Create `scripts/build_biogears_gi_variant.py`, `scripts/verify_native_gi_integrity.py`, `docs/research/GI_MASS_INTEGRITY.md`; generated variant under `data/runtime/physiology/variants/whole_body_integrity`; original-failure and corrected probes under `data/derived/audits/gi-integrity`.

**Interfaces:** Consumes retained source revision3f16a5fa1dade9c511b88d923606fa51cc35e95d and existing saturation/heatflux/thermal-unit variant build objects. Produces a standalone libbiogears.so.8.0.0 and manifest with library_sha256, exact original/patched source digests and patch text; root native configuration later registers this variant.

- [ ] Write a failing source/native calcium-transfer regression. A transfer amount calculated in mg must decrement stomach and increment chyme by the same mg; include tiny remainder, zero contents and the retained negative native stomach value.
- [ ] Execute the regression against retained original evidence, retain its failure, and inspect exact upstream branch before editing.
- [ ] Build an isolated variant that caps transferred calcium at available nonnegative stomach mass, uses mg for both endpoints and empties a positive remainder. Preserve calcium absorption-fraction semantics and all macronutrient rate equations; document their unresolved calibration.
- [ ] Compile/link the new object against existing native objects without changing donor source or older variants. Run a real source-state calcium experiment and verify finite/nonnegative contents and the local paired-transfer equality; do not infer global calcium conservation from a local branch check.
- [ ] Record build/runtime hashes, test results and source limitations; commit only owned files and request task review.

### Task 2: Persistent physiological stepping and organ ports

**Files:** Create `scripts/native_body_ports.h`, `scripts/native_biogears_stream.cpp`, `ihm/native/session.py`, `scripts/verify_native_session.py`; modify `scripts/build_native_adapter.py` and `ihm/native/__init__.py`.

**Interfaces:** `NativeSession(config, output_dir)` starts a pinned local engine; `snapshot()` returns typed current native quantities, `step(seconds)` advances integer0.02s ticks, `apnea(severity)` applies the actual native action, `meal(**composition)` delivers carbohydrate/protein/fat/sodium/calcium in g and water in mL, `exercise(intensity)` applies native demand, `save_state()` writes a native checkpoint, `close()` closes/terminates only its owned process. Native stdout responses have a dedicated prefix and sequence; unrelated engine logs are retained separately. A context manager guarantees process cleanup.

- [ ] Test validation before native startup: nonfinite/negative quantities, invalid steps, empty meals and unsupported variant names fail without creating a run. Test expected current APIs with an intentionally missing stream binary before building.
- [ ] Implement the C++ command loop and actual source API calls. Emit pressure/frequency requests, cycle phase, applied pressure/flow, blood gases, autonomic scales, stomach contents, selected chyme/vascular solute masses, plasma glucose/insulin, glycogen, metabolic rate and renal output with explicit units and source identifiers.
- [ ] Implement Python process/log handling, timeouts, action receipts, immutable source/library identity, sequence/clock validation and cleanup. Preserve full native state for restart; record known serialization gaps instead of claiming exact replay without a test.
- [ ] Run native batch/session parity from the same saved state with no actions, then drive-suppression/restoration and matched meal delivery. Assert exact clock advancement and changed physical observables, not just changed command labels.
- [ ] Keep legacy NativeConfig and adapter behavior compatible; run native/body regression checks and commit the new executable boundary.

### Task 3: Shared systemic experiments and source interaction audits

**Files:** Create `ihm/assembly/systemic.py`, `scripts/build_systemic_experiments.py`, `scripts/verify_systemic_experiments.py`, `docs/research/SYSTEMIC_CAUSAL_EXPERIMENTS.md`; modify `ihm/human.py` and `ihm/app/experiments.py`.

**Interfaces:** `run_systemic(root, output_dir, protocol, seconds, sample_hz)` uses NativeSession as sole physiological owner. Recorded schema contains frames, fields with units/support/owner, actions, mechanism_edges, source receipts and numerical limitations. `human.materialize('body-systemic', **options).run(output_dir)` creates the predictor; API exposes acquired completed experiment records.

- [ ] Specify and test protocols: matched resting/meal states, apnea/restoration, exertion/recovery and combined meal/exertion. Reject overlapping conflicting control actions and clocks outside native ticks.
- [ ] Run actual long enough native trajectories to observe digestion, vascular/hormonal response and energy storage, retaining native initial/final states. Record calcium integrity and finite/positive storage. Retain branch-counter gaps where source telemetry cannot establish absorption flux.
- [ ] Bind respiratory/vascular/thermal outputs to canonical body materialization using matching source states and explicit brainstem control records. Preserve old replay artifacts and output ownership.
- [ ] Add sensitivity comparisons that require downstream changes (ventilation and gas state, absorbed substrates and hormonal/storage state), and source evidence for physiological parameter meaning.
- [ ] Audit broader native interactions and unresolved systems with executing-edge status; commit results and scientific limitations.

### Task 4: Native articulation, ground contact and motor control

**Files:** Create `scripts/native_opensim_dynamics.cpp`, `scripts/build_native_dynamics.py`, `ihm/native/locomotion.py`, `scripts/verify_native_locomotion.py`, `docs/research/LOCOMOTION_DYNAMICS.md`.

**Interfaces:** `LocomotionConfig` records source model, controller mode, duration, initial state, perturbations and solver tolerance. `run_locomotion(config, output_dir)` returns native body poses/velocities, muscle excitations/activations/forces, ground forces and mechanical work with exact model/library/source receipts. Tracking targets are input data; integrated coordinates remain solved states.

- [ ] Use the retained3D walking model/contact sets and measured reference data identified in LOCOMOTION_CONTROL_AUDIT.md; verify hashes, body masses, state dimensions and attachment paths before simulation.
- [ ] Build the required native time integration/controller or optimization dependencies locally. Retain failed builds and original source. Execute a forward mechanical interval with endogenous contact and unprescribed coordinate evolution.
- [ ] Verify force/gravity balance, contact direction, finite strain-independent articulation, state clocks and tolerance refinement. Contrast controller removal or changed ground loading with the baseline. Keep prescribed/measured-GRF tracking as a separate labeled mode if used for initialization.
- [ ] Record computed muscle work/activation and source reference agreement. Distinguish a successful finite forward interval from sustained stable/autonomous walking; extend the controller and horizon until the declared walking scenario passes.
- [ ] Commit the native dynamics bridge and review its anatomical/control claims.

### Task 5: Heterogeneous mechanical interaction and canonical attachment

**Files:** Create `ihm/assembly/material_domains.py`, `ihm/assembly/contact_dynamics.py`, `scripts/build_body_material_domains.py`, `scripts/verify_material_domains.py`, `scripts/verify_contact_dynamics.py`; modify `ihm/assembly/mechanics_backend.py` only for audited shared primitives.

**Interfaces:** Material domains contain source entity IDs, reference nodes/elements, owned mass, constitutive cards, attachment boundaries and contact surfaces. Dynamic state contains positions/velocities and energy/work/impulse ledgers. Anatomical registration maps source-native body transforms to material attachments and records unsupported assignments.

- [ ] Test rigid motion invariance, nonoverlapping ownership, finite positive element Jacobians, equal/opposite internal reactions and source-frame consistency on analytic fixtures.
- [ ] Materialize source-attached heterogeneous volumes with explicit inferred partitions; reject invalid/overlapping material ownership rather than silently duplicating organ mass.
- [ ] Extend contact from quasistatic reference indentation to time-dependent coupled material domains, including external gravity/ground work and source articulation boundary reaction. Compare nonuniform contact to the retained FEBio backend before asserting accuracy.
- [ ] Run mechanical perturbations and mesh/time refinement; record source resolution and material-prior sensitivity separately. Scale across anatomical domains and list unsupported interfaces.
- [ ] Commit the mechanics materialization and reference evidence.

### Task 6: Causal sensory–motor and neuromechanical coupling

**Files:** Create `ihm/assembly/reflexes.py`, `scripts/verify_body_reflexes.py`, `scripts/build_neuromechanical_experiments.py`; modify `ihm/assembly/body_runtime.py`, `ihm/assembly/body_protocol.py`, `ihm/assembly/body_states.py` and checkpoint interfaces.

**Interfaces:** Reflex controllers consume computed sensory/force/length/velocity state and optional descending commands; produce bounded actuator commands after specified delays. State/checkpoint includes delay buffers, neural/controller state and mode transitions. Native muscle work and physiological demand use a declared observation law with work/energy units.

- [ ] Acquire primary controller equations/parameters and distinguish spinal, brainstem and cortical paths; keep transferred gains/thresholds explicit. No arbitrary unconditional touch-to-muscle assignment.
- [ ] Test afferent delay, efferent delay, selective block, release/decay, descending modulation, checkpoint equivalence and absence of output before sensory arrival.
- [ ] Close the mechanics→sensor→controller→actuator→mechanics loop and feed declared metabolic demand into the native physiological owner. Record the distinction between muscle mechanical work and metabolic cost.
- [ ] Run coupled perturbations demonstrating changed motion, force and physiological state, with one owner per state and no hidden recorded-trajectory forcing.
- [ ] Commit and review source fidelity and coupling claims.

### Task 7: One-body inspection, temporal analysis and acceptance

**Files:** Modify `app/src/main.js`, `app/src/regional.js`, `ihm/app/__init__.py`, `scripts/build_body_experiment_spectra.py`, `scripts/audit_body_coverage.py`, `scripts/verify_all.py`; add scoped browser tests and update research status.

- [ ] Show systemic/neuromechanical state, control requests versus applied drive, nutrients and blood-borne quantities on the canonical body. Use one materialization inventory, with explicit source-space inspection for unregistered geometry.
- [ ] Preserve 3D-first panels and physical camera scales; provide computed time playback and selected causal-chain readouts. No decorative peristalsis, walking or breathing in place of model state.
- [ ] Derive finite-window spectra from actual uniformly sampled signals with appropriate time units and Nyquist limits. Source freshness checks must include transitive inputs and native libraries.
- [ ] Run the appropriate complete numerical/source/native/browser suite and independent broad review; fix concrete findings before committing acceptance evidence.
- [ ] Record what became integrated and what remains scientifically unresolved across all body systems. Passing the current experiment suite is not universal calibration or complete anatomical realism.
