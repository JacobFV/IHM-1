# Human Simulation Exploration and Implementation Plan

> Planning only. The user has explicitly paused implementation. On later authorization, use the executing-plans workflow task-by-task; do not infer execution or delegation approval from this document. Every unchecked item is future work.

**Goal:** Build one generic human whose physical, physiological and neural interactions support the observable experiments in the concern register, with an immersive 3D workspace and explicit evidence/uncertainty at every materialization.

**Architecture:** One versioned implicit body owns identities, source constraints, state definitions and exchange interfaces. Reference mechanical, vascular, electrophysiological and neural models provide heterogeneous numerical realizations; reduced predictors and display geometry are derived materializations. A conservative co-simulation runtime coordinates state owners and event timing.

**Tech stack:** Existing Python/NumPy/SciPy, C++ BioGears adapter, IBM materialization/runtime, OpenSim source assets, Three.js/Vite and the existing HTTP workbench. FEBio, SOFA and SimVascular are candidate backends to evaluate before selection, not approved new dependencies.

**Spec:** [Requirements and concern register](../specs/2026-09-05-human-simulation-concerns.md).

## Global constraints

- Implementation is paused. This document and its companion roadmap are planning deliverables only.
- One generic adult human with coherent anatomy and physically meaningful interactions, initially resting supine on a bed.
- Detailed structures may be synthesized, but their assumptions, uncertainty and consequences must remain inspectable.
- The viewport is the primary workspace; information panels can appear/disappear and scroll independently.
- Single authoritative owner per physical state; source-derived variants cannot silently create independent subjects or duplicated mass, energy, fluid, charge or control.
- Preserve original source bytes and canonical IDs; all repaired, registered, synthesized and reduced outputs retain transformation lineage.
- No passing numerical test is described as whole-body empirical validation.
- Exact solver choice and all physiological tolerances are evidence decisions. Numerical acceptance targets below are proposed engineering gates, not biological confidence intervals.

## 1. Work organization and decision gates

Use an exploration gate before committing to a new domain solver or full implementation. Each exploration has a concrete question, input set, output artifact and decision rule. The exploration outputs include negative findings: missing data, contradictory sources, unsupported platforms and failed convergence are results, not reasons to create replacement numbers silently.

Dependency order:

```text
E0 baseline and coverage ── E1 evidence/geometry ── E2 mechanics feasibility
                    ├──── E3 IBM and pathways ──┐
                    ├──── E4 circulation/skin ───┼─ E6 architecture decision
                    └──── E5 compute/UI ────────┘

I1 contracts/evidence ── I3 material geometry ── I4 mechanics/contact
          ├──────────── I2 fixed workspace          ├─ I5 breathing
          ├──────────── I6 IBM integration ─────────┴─ I7 sensory/motor
          ├──────────── I8 vascular detail ── I9 perfusion/lymph/skin
          └──────────── I10 hair ────────────────┘
                               I11 bioelectricity and tissue response
                      I12 remaining systems and whole-body expansion
                      I13 calibration/reductions/spectra/release gates
```

I1 precedes every physical integration. I2 can be verified independently once identities and observation semantics are fixed. Reference mechanics and IBM compatibility can be investigated independently; their integrated sensorimotor experiment depends on both. Regional examples are milestones, not the definition of full coverage.

Do not begin by sprinkling procedural hairs, vascular branches and delayed nerve lines over the current view and call the result complete. First establish where each representation lives, what state it owns, how its interactions are computed and which observations would falsify it.

## 2. Exploration work packages

### E0 — Freeze the audit baseline and enumerate actual execution

**Questions:** Which artifacts exist? Which processes execute in the canonical run? Which later drafts/commits have not been integrated? Which runtime is the app actually serving?

**Inspect:** git revisions and dirty files; `ihm/human.py`; `ihm/assembly/`; native adapter/configuration; canonical manifests and source receipts; workbench build pipeline; app routes and loaded asset versions; IBM's revision and materialization contracts.

- [ ] Record the starting revision, all pre-existing edits, generated-artifact hashes and live-server code/build identity without overwriting them.
- [ ] Build a read-only inventory of anatomical entities, physical DOFs, solved state variables, state owners, active processes and actual nonzero exchanges.
- [ ] Classify each claimed capability as declared, implemented, numerically exercised, integrated, empirically constrained or empirically validated. These are separate columns, not a single score.
- [ ] Reconcile the early panel, respiration, peripheral and mechanics work against this plan. Retain, revise or reject each explicitly; no automatic acceptance because code exists.

**Outputs after authorization:** `data/derived/audits/execution-coverage.json`, `docs/research/EXECUTION_BASELINE.md`.

**Exit gate:** Every capability advertised in the UI has a corresponding artifact/runtime path and validation status. Seven moving surfaces cannot be summarized as full-body mechanics.

### E1 — Anatomical, material and parameter evidence audit

**Questions:** Which surfaces represent material, boundaries, lumens or assemblies? Which interfaces/attachments are known? What is missing in each body region? Which coefficient transfers are supportable?

- [ ] Compare canonical topology against BP/Z anatomy ontology, source OpenSim joints/paths, available vascular centerlines and anatomical reference datasets.
- [ ] Audit geometry at head/neck, thorax, abdomen/pelvis, each upper limb and each lower limb; distinguish duplicate parent surfaces from physical material.
- [ ] Measure seam duplication, open boundaries, intersections, connected components, local mesh quality and ambiguous tissue interfaces. Do not weld nearby anatomical structures solely because coordinates are close.
- [ ] Create an interface inventory for joints, tendon insertions, fascia, pleura, pericardium, organ suspensions, skin/subcutis and bed contact.
- [ ] Extract parameter cards: value/unit, tissue/region, species, subject cohort, loading rate/direction, temperature, boundary conditions, reference state, source location, uncertainty, dependence group and transfer assumptions.
- [ ] Design the minimum missing geometry for each first acceptance experiment; specify source-constrained synthesis and alternative realizations where data are insufficient.

**Outputs:** `anatomy-gaps.json`, `interfaces-evidence.json`, `parameter-cards.json`, `geometry-repair-proposals.json`, and a source-access/license ledger under `data/derived/audits/`.

**Exit gate:** Every first-experiment material domain and interface is either evidence-backed or explicitly modeled as an assumption with a sensitivity range. No unexplained global density multiplier substitutes for this audit.

### E2 — Mechanics and respiratory solver feasibility

**Questions:** Which backend can solve the required tissue/contact/joint behavior robustly on this machine? Can chest-wall loading alter ventilation through an actual supported BioGears interface?

- [ ] Compare FEBio, SOFA and the current reduced mechanics using the same reference units, material law, boundary conditions and test geometry.
- [ ] Specify three benchmark problems: layered skin indentation on a stiff substrate; diaphragm/chest cavity pressure-volume loading; a tendon-driven hinged limb with surrounding soft tissue contact.
- [ ] After execution approval, benchmark each candidate at increasing volume-element/constraint counts. Record build support on ARM64, CPU/GPU mode, license, wall time, memory, residuals, nonlinear iterations, contact error and mesh/time refinement.
- [ ] Inspect BioGears source setters/internal ownership to determine whether external respiratory muscle pressure, pleural pressure or mechanical compliance can participate in the integration loop. A telemetry getter alone is insufficient.
- [ ] Define a retained-snapshot replay mode and a distinct closed-coupling mode. Identify precisely which native respiratory states must be replaced or externally driven in the second mode.

**Outputs:** `docs/research/MECHANICS_BACKEND_DECISION.md`, `mechanics-benchmarks.json`, `native-respiratory-port-audit.json`.

**Exit gate:** One candidate passes the regional numerical gates and has a supported integration path. If the native model cannot accept the chosen feedback, the plan must include a source-level adapter or a replacement respiratory submodel with explicit boundary ownership; no animation fallback qualifies.

### E3 — IBM reuse and somatic pathway compatibility

**Questions:** What does IBM actually implement for receptors, motor units, afferent/efferent stages, spectral state and uncertainty? Which pieces can be materialized against IHM's body supports without rewriting their meaning?

**Read in IBM:** `docs/{ARCHITECTURE,CONTRACT,EVIDENCE}.md`; `ibm/materialize/{request,build,model,sites,geometry,provenance}.py`; `ibm/runtime/{state,step,ensemble,fuse,intervene}.py`; `ibm/fields/{transduction,effector,neural}.py`; `ibm/processes/{transduction,effector,neural}.py`; `ibm/topologies/{afferent,efferent}.py`.

- [ ] Enumerate implemented laws separately from declarations, aliases and unavailable data dependencies.
- [ ] Audit the actual somatotopy, receptor supports, relays, corticospinal/corticobulbar chains and fiber-class delay assumptions already present in IBM.
- [ ] Define donor-versus-integrated parity fixtures for tactile pulse, thermal step, delayed afferent relay and motor-unit drive; record required temporal representation and original outputs.
- [ ] Determine integration packaging: a version-pinned local package or a vendored dependency closure with source receipts. Do not copy an entire dataset tree indiscriminately or compile a few functions while claiming full IBM reuse.
- [ ] Audit every canonical muscle/effector and skin territory for innervation, mixed supply and missing roots/branches; keep geometry certainty separate from connectivity certainty.

**Outputs:** `docs/research/IBM_BODY_COMPATIBILITY.md`, `ibm-feature-parity.json`, `innervation-coverage.json`, `body-brain-port-contract.json`.

**Exit gate:** The proposed brain materialization names exactly which IBM capabilities execute and which remain unsupported. A stimulus-gated cortical-rate multiplier is not accepted as a physiological reflex controller.

### E4 — Vascular, lymphatic, skin, hair and bioelectric evidence

- [ ] Find source-supported named vessel topology/radii and terminal territories; identify acquisition resolution and missing branches before synthesizing them.
- [ ] Gather regional arterial/capillary/venous density, length/radius distributions, hematocrit/rheology, extraction and exchange evidence. Test applicability at each modeled scale; capillary behavior is not automatically described by a large-vessel Newtonian law.
- [ ] Audit lymph flow direction, valves, pumping, interstitial compliance/protein transport and venous return against source evidence.
- [ ] Gather regional skin layer thickness/material, appendage density, receptor distributions, innervation and terminal/vellus hair statistics. Mark glabrous/anatomical exclusion regions explicitly.
- [ ] Define which skin electrical quantities are available: cell membrane potentials, transepithelial potential, extracellular electric field, ion concentrations/conductivities, gap-junction conductances and wound boundary conditions.
- [ ] Distinguish measured human skin constraints from model-organism bioelectric findings and conceptual morphogenesis hypotheses. Record what any proposed cellular response law is allowed to predict.

**Outputs:** `docs/research/REGIONAL_TRANSPORT_AND_SKIN.md`, `vascular-territories.json`, `skin-region-priors.json`, `bioelectric-evidence-map.json`.

**Exit gate:** A specified forearm/skin reference region has arterial supply, venous return, interstitial/lymph boundaries, mechanical material, receptor and electrical boundary definitions. Unknown quantities are identifiable parameters or explicit model-form uncertainty.

### E5 — Compute, visualization and data movement

**Observed machine, not a performance promise:** 20 ARM64 CPU threads, about 121 GiB total reported RAM and an NVIDIA GB10. Reported GPU memory was unavailable in the inventory query; do not assume a separate additive VRAM pool.

- [ ] After approval, measure a fixed synthetic rendering workload at 1M, 5M, 10M and 20M triangles; repeat with shafts, vascular instances, picking and one deforming region.
- [ ] Independently scale solver variables, volume elements, contact pairs, vascular edges and neural spectral coefficients. Triangle count is not the simulation budget.
- [ ] Measure artifact load/decompression, CPU/GPU copies, browser memory, frame-time distribution and simulation-to-render lag.
- [ ] Define fidelity presets by observed error and cost: offline reference, interactive reduced mechanics, fine regional inspection and whole-body overview. A preset changes the materialization contract and exposes its validity envelope.
- [ ] Verify panel state/camera retention and no document scroll at 1440×900, 1024×768 and a narrow viewport. Explore resizable versus fixed panel dimensions without rebuilding anatomy.

**Outputs:** `docs/research/COMPUTE_AND_VIEWPORT_BUDGET.md`, `compute-benchmarks.json`, proposed workspace acceptance screenshots/wireframes.

**Exit gate:** Supported memory/performance envelopes are measured. Proposed display target: p95 frame time ≤33 ms for the agreed interactive preset on the target machine; a 60 Hz target is additional, not promised. Offline physical accuracy takes precedence over keeping the full reference solve real time.

### E6 — Architecture decision and executable backlog

- [ ] Review E0–E5 outputs against all C01–C36 concerns.
- [ ] Select mechanics/coupling backends and IBM packaging; freeze source versions and units.
- [ ] Select empirical datasets for each first acceptance experiment, reserve holdout subjects/conditions, and freeze scientific metrics before calibration.
- [ ] Replace provisional targets below with measured workload limits and justified experiment-specific error thresholds.
- [ ] Split implementation into independently reviewable work packages, estimate each from the spikes and identify high-risk dependencies. Do not promise an overnight delivery date for unbenchmarked multiphysics work.

**Exit gate:** User reviews the concrete design and authorizes implementation. No build is triggered automatically by finishing this plan.

## 3. Shared contracts to implement after approval

Proposed serialized records, not current APIs:

```text
BodyEntity:
  id, anatomical_concepts, region, laterality, material_domain_id,
  reference_frame, source_constraints[], uncertainty_id

MaterialPoint:
  body_version, entity_id, element_id, local_coordinates
  # Stable under deformation; independent of display vertex ordering.

Interaction:
  id, endpoint_a, endpoint_b, kind, orientation,
  law_id, parameters_id, active_domain, evidence_ids[], unresolved_reason
  # kind: joint, attachment, contact, sliding, flow, diffusion,
  #       neural_afferent, neural_efferent, receptor, gap_junction, thermal.

StateSpec:
  id, owner, support, unit, conserved_quantity_or_none,
  spatial_discretization, temporal_representation, validity_domain

Exchange:
  interface_id, interval_start_s, interval_end_s,
  quantity, amount_SI, source_owner, target_owner, work_J_or_none
  # Rate samples and integrated amounts are explicitly different records.

ObservationRequest:
  body_version, region, observables[], physical_resolution,
  temporal_bandwidth_Hz, error_target, compute_budget, intervention

MaterializationReceipt:
  data_hashes, code_hashes, solver_versions, configuration, seed,
  state_owners, transfers, error_estimates, uncertainty, validity_domain

Frame:
  model_time_s, receipt_id, transforms, deformation_fields,
  physiological_fields, neural_fields, observations, audit, playback_kind
```

Proposed execution interfaces:

```text
materialize(request: ObservationRequest) -> ExplicitBodyModel
initialize(model, reference_conditions) -> Checkpoint
advance(checkpoint, end_time_s, interventions) -> (Checkpoint, Observations, Audit)
observe(checkpoint, material_points, quantities) -> values with units/uncertainty
project_for_display(checkpoint, view_request) -> Frame
```

Each engine adapter must support checkpoint/restore for rejected coupling steps. Continuity of mass, species, charge, momentum and energy must be assessed on integrated exchanges. Prescribing displacement while also independently solving the same motion is an ownership error. Contact/attachment constraints transfer work; interpolation of values alone is not a conservative coupling rule.

## 4. Implementation work packages

Every package follows the same concrete review cycle: write the specified counterexample/regression fixture; confirm the expected failure in the old implementation; implement that package; run its named verification; inspect its physical/visual output; document measured errors and unresolved scope; commit only the reviewed package. None of these steps is being executed during this planning pass.

### I1 — Evidence, identities, state ownership and coupling infrastructure

**Concerns:** C01–C05, C27, C29–C31, C35–C36. **Depends on:** E6.

**Files:** extend `ihm/human.py`, `ihm/assembly/body.py`, `certainty.py`; create `ihm/assembly/contracts.py`, `interfaces.py`, `cosimulation.py`, `evidence.py`; create `scripts/audit_body_coverage.py`, `scripts/verify_body_contracts.py`, `scripts/verify_body_cosimulation.py`.

**Consumes:** audited entity/interface/evidence records. **Produces:** the shared records and execution interfaces above, field ownership registry, checkpointed coupling and machine-readable coverage.

- [ ] Implement stable material coordinates and typed interactions; reject missing IDs, dimension mismatch and dual state ownership before integration.
- [ ] Implement a source/repair/synthesis/fit/reduction dependency graph and parameter cards; represent unknown uncertainty without probability values.
- [ ] Implement interface interval accounting, explicit coupling order, nonlinear convergence/rollback and error estimates. Separate co-simulation from recorded source replay.
- [ ] Implement tests: a two-reservoir exchange preserves combined storage; interface energy/work signs cancel; a rejected step restores every engine and event queue; halved coupling step converges; duplicated owner is rejected.
- [ ] Require each downstream work package to register all states/ports and unresolved dependencies before enabling a scenario.

**Verification:** `.venv/bin/python scripts/verify_body_contracts.py` and `scripts/verify_body_cosimulation.py`.

### I2 — Fixed 3D workspace and inspection contract

**Concerns:** C33–C35. **Depends on:** E5 and I1 identity/observation semantics. Earlier UI changes must be reviewed before reuse.

**Files:** `app/src/main.js`, `style.css`, `state.js`; create `app/src/workspace.js`, `scene.js`, `body-playback.js`; extend `app/test/browser/body.spec.js`; create `app/test/browser/workspace-layout.spec.js`.

**Consumes:** Frame and Observation records. **Produces:** a fixed viewport with independently scrollable, toggleable and resizable left/right/bottom panels.

- [ ] Give the workspace the available viewport height; use constrained grid/flex tracks and `min-height: 0`; only panels own scrolling.
- [ ] Keep show/hide controls reachable with every panel hidden; preserve camera, selection, clipping, layer state and simulation state when panels are hidden or restored.
- [ ] On narrow screens, panels overlay or replace each other without removing the viewport's primary role. Hidden panels leave the focus order; controls expose expanded state.
- [ ] Tie inspector readouts and plot cursors to the visible simulation timestamp; distinguish recorded replay, interpolation, current simulation and prediction. Do not reset or advance physics through a layout event.
- [ ] Support anatomical selection separately from physical stimulation. A skin press has explicit position/material coordinates, area, force or indentation mode, duration and units; orbit dragging is never a stimulus.
- [ ] Verify eight combinations of three visible/hidden panels, independent wheel scroll, focus mode restoration, resize, missing WebGL fallback and stale-data messages.

**Acceptance examples:** document scroll remains zero during panel scrolling; hiding both sidebars increases canvas width; reopening them does not recreate the scene or change the camera; renderer size follows container size; a paused field stays at the same time while the inspector opens.

**Verification:** `npm --prefix app test`, `npm --prefix app run build`, and the targeted Playwright layout/body suites.

### I3 — Physical geometry, tissues and interfaces

**Concerns:** C02–C06, C17, C31. **Depends on:** I1 and E1.

**Files:** extend `ihm/assembly/anatomy.py`, `scripts/build_canonical_anatomy.py`; create `ihm/assembly/domains.py`, `scripts/build_body_domains.py`, `scripts/verify_body_domains.py`.

- [ ] Classify geometry into material volumes, shells, lumen/boundary surfaces and composite assemblies; resolve overlaps before mass integration.
- [ ] Derive simulation domains with local refinement, fibers, thickness and contact surfaces; retain full source surfaces and repair lineage.
- [ ] Add fascia, organ supports, joint/cartilage interfaces and bed contact surfaces for the first thorax/forearm experiments, then expand by region.
- [ ] Build material-to-simulation and simulation-to-display mappings; show error bounds and ambiguous regions.
- [ ] Test material mass exclusions, closed-domain volumes, element Jacobians, orientation consistency, attachment ownership and mapping under known rigid/affine deformations.

**Deliverable:** a meshed reference region with legitimate mechanical boundaries, not only more triangles. **Verification:** `scripts/verify_body_domains.py`.

### I4 — Deformable/contact mechanics and articulated force transfer

**Concerns:** C04–C09. **Depends on:** I1/I3 and selected E2 backend.

**Files:** retain `ihm/assembly/mechanics.py` as an explicitly bounded coarse model; create `ihm/assembly/mechanics_backend.py`, `mechanical_mapping.py`, `scripts/verify_body_contact.py`, `scripts/verify_body_articulation.py`; extend `scripts/build_body_mechanics.py`.

- [ ] Implement the backend adapter with rigid articulated bodies, tissue solid/shell domains, tendons/ligaments, contact/sliding and state checkpoints.
- [ ] Initialize gravity-loaded supine equilibrium with bed reaction, resting tone/prestress and a reconciled reference configuration.
- [ ] Implement source-consistent muscle/tendon behavior and dynamic wrapping or constrained path updates. Check moment arms and torque, not only line-force cancellation.
- [ ] Transfer equal/opposite forces and consistent work across coarse/fine and rigid/soft interfaces.
- [ ] Verify indentation against reference curves, joint moment arms, tissue contact separation, rigid-body objectivity, force/work balance and mesh/time refinement; run sustained loading and recovery.

**Deliverable:** pressing the forearm or loading a muscle causes local deformation, correctly constrained motion and reaction on neighboring tissue/bone. **Verification:** contact/articulation scripts and a recorded reference-solver comparison.

### I5 — Respiratory mechanics and closed pressure/flow coupling

**Concerns:** C06, C10–C11, C27. **Depends on:** I4 and native-port E2 gate.

**Files:** review/replace draft `ihm/assembly/respiration.py`; create `respiratory_coupling.py`; extend `scripts/native_biogears_rest.cpp`, `ihm/native/__init__.py`, `ihm/assembly/body.py`; extend `scripts/build_body_respiration.py`, `scripts/verify_body_respiration.py`.

- [ ] Create diaphragm, intercostal, rib/cartilage, sternum, pleural, lung and abdominal mechanical interfaces.
- [ ] First verify mechanical response to a known pressure/volume boundary. Label this prescribed-input mode accurately.
- [ ] Implement the approved physiological port so respiratory drive and mechanical muscle effort produce pressure/flow and evolving gas volume; remove conflicting ownership in the coupled mode.
- [ ] Include lung/chest recoil, airway resistance, gas transport and contact as required by the reference experiment. Account for pressure-volume work and boundary energy.
- [ ] Compare diaphragm excursion, chest marker displacement, tidal volume, airway/pleural pressure and phase lag against reserved targets.
- [ ] Test changed chest-wall stiffness, added external thoracic load and airway resistance; feedback must change pressure/flow/volume consistently. Rendering the skin and ribs must use this state, with no secondary oscillation.

**Verification:** existing respiratory unit checks plus a closed-loop load perturbation and a recorded chest/diaphragm visualization with simultaneous physical traces.

### I6 — Actual IBM materialization/runtime integration

**Concerns:** C21, C25, C28–C30. **Depends on:** I1 and E3.

**Files:** create `ihm/brain/ibm_backend.py`, `body_supports.py`, `port_mapping.py`; extend `ihm/assembly/brain.py`, `scripts/build_body_brain.py`; create `scripts/verify_ibm_body_parity.py`.

- [ ] Package the selected IBM dependency closure at a recorded revision and preserve data/code license/source receipts.
- [ ] Map canonical body skin, receptor, motor-unit, spinal/brainstem and brain supports to IBM's support/region semantics, including physical units, frames and laterality.
- [ ] Execute an IBM request through its actual materializer and runtime, including chosen temporal and uncertainty representations.
- [ ] Preserve the existing 80-population model as a separately declared reduced materialization until parity/error tests qualify replacements.
- [ ] Verify donor and integrated tactile/thermal/motor fixtures under identical inputs, parameters and temporal basis; reject silently empty region mappings or unsupported processes.

**Deliverable:** an IBM-backed neural materialization that can consume and emit the agreed body ports, with an explicit reuse/coverage report. **Verification:** `scripts/verify_ibm_body_parity.py`.

### I7 — Skin sensing, peripheral wiring, spinal circuits and motor feedback

**Concerns:** C09, C19–C22, C25. **Depends on:** I4/I6; skin domain from I3.

**Files:** review `ihm/assembly/peripheral.py`, `scripts/build_body_peripheral.py`; create `receptors.py`, `innervation.py`, `spinal_circuits.py`, `scripts/verify_body_sensorimotor.py`; add stimulation/observation routes in `ihm/app/__init__.py` and app controls.

- [ ] Build named motor/sensory pathways with cranial/spinal roots, plexuses, branches, fiber class, path length, delay and mixed-innervation rules. Label geometric routes separately from functional connectivity.
- [ ] Anchor spatial receptor populations to skin and muscle material points. Implement modality-specific transduction, adaptation/sensitization, spindle/tendon feedback and neuromuscular activation dynamics using audited IBM laws where available.
- [ ] Implement explicit local reflex circuits and separately specified descending commands. Innocuous touch does not automatically imply withdrawal or arbitrary motor activity.
- [ ] Couple physical contact/temperature to receptor inputs and motor output to actual muscle mechanics; return proprioception from the computed muscle/tendon state.
- [ ] Test pressure pulse/release, vibration-frequency response, thermal response, imposed tendon stretch and selective nerve block. Verify causality, laterality, path delays, loss of target response and preserved unrelated responses.
- [ ] Publish coverage separately for skin area/territories, receptor classes, nerves, motor targets and functioning closed loops. No “all nerves connected” claim until those inventories resolve.

**Verification:** `scripts/verify_body_peripheral.py` and `scripts/verify_body_sensorimotor.py`, including an interaction-to-neural-to-mechanical trace visible in the app.

### I8 — Named vascular networks and source-constrained microvascular detail

**Concerns:** C13–C15, C31–C34. **Depends on:** I1/I3/E4.

**Files:** create `ihm/assembly/vascular.py`, `microvascular.py`, `scripts/build_body_vascular.py`, `scripts/verify_body_vascular.py`; reuse audited regional CFD/circuit code through explicit adapters; add vascular rendering to the scene module.

- [ ] Reconstruct or acquire validated lumen centerlines, branch identities/radii and terminal perfusion territories; pair arterial supply and venous return with tissue beds.
- [ ] Define a hierarchical stochastic network conditioned on territory geometry and measured regional statistics. Keep a persistent seed and identities so zooming cannot invent a different vascular realization.
- [ ] Resolve explicit microvascular graphs in selected regions and homogenized exchange elsewhere. Branching trees alone are insufficient where capillary plexus architecture is required by the evidence.
- [ ] Solve pressure/flow and relevant transport using laws appropriate to each scale; connect regional boundary flow and storage to the global physiological owner.
- [ ] Verify branch conservation, known network solutions, positive/passive behavior where expected, species transport, oxygen extraction, return flow and occlusion effects. Compare regional statistics and flow targets against held-out constraints.
- [ ] Render vessels from the same radius/topology/state records; use physical radius by default and label any visibility exaggeration.

**Deliverable:** fine vessels with traceable supply/return and solved field provenance, rather than free-floating decorative branches. **Verification:** `scripts/verify_body_vascular.py`.

### I9 — Regional skin/interstitial/lymph transport and systemic coupling

**Concerns:** C15–C17, C26–C27. **Depends on:** I4/I8.

**Files:** create `ihm/assembly/skin_transport.py`, `lymph_transport.py`, `scripts/build_body_transport.py`, `scripts/verify_body_transport.py`; bridge existing `ihm/coupling/` native circuits.

- [ ] Implement skin/subcutis fluid storage, microvascular exchange, protein transport, initial lymph uptake, valve/pump transport and venous return.
- [ ] Couple deformation/pressure to perfusion and interstitial transport within the chosen constitutive domain; preserve the physical meaning of transmural pressure and compliance.
- [ ] Assign thermal, sweat and water transport to a single owner, with explicit handoff if replacing native/JOS3 regional behavior.
- [ ] Verify controlled inflow reduction, raised venous pressure and lymph obstruction; compare tissue fluid/oxygen/temperature trajectories and all integrated boundary balances.

**Verification:** transport conservation, refinement and source-target regression in `scripts/verify_body_transport.py`.

### I10 — Body hair and appendage materialization

**Concerns:** C17–C18, C19, C31/C34. **Depends on:** I3/E4; mechanical/sensory coupling depends on I4/I7.

**Files:** create `ihm/assembly/hair.py`, `scripts/build_body_hair.py`, `scripts/verify_body_hair.py`; add instanced shaft rendering and follicle inspection in the scene module.

- [ ] Generate terminal/vellus follicle populations conditioned on region, exclusions and uncertain density/orientation/radius/length distributions.
- [ ] Store follicles in material coordinates so they move with skin and survive remeshing/display LOD through a recorded mapping.
- [ ] Materialize shaft geometry independently of visual sampling density; preserve aggregate physical descriptors for insulation/contact where modeled.
- [ ] Add shaft bending and follicle-force receptor input in the selected detailed region. Distinguish decorative appearance, mechanical shaft and sensory follicle capability in coverage.
- [ ] Test no roots in excluded regions, statistical counts within the sampling model, reproducibility, deformation attachment and unchanged underlying follicle population when display density changes.

**Verification:** `scripts/verify_body_hair.py` and a zoom/deformation browser test.

### I11 — Integumentary bioelectricity and supported tissue response

**Concerns:** C17, C23–C24, C27–C30. **Depends on:** I1/I3/I9 and E4 evidence gate.

**Files:** create `ihm/assembly/skin_bioelectric.py`, `electrical_coupling.py`, `scripts/build_body_skin_electric.py`, `scripts/verify_body_skin_electric.py`; integrate the existing BETSE/human wound materializations through explicit boundary contracts.

- [ ] Define separate intracellular/extracellular potentials, membrane potential, epithelial voltage, ion concentrations and currents with unit/reference conventions.
- [ ] Implement or adapt channel/pump/gap-junction and transport laws for the chosen local skin domain; couple perfusion, temperature, hydration and wound boundaries only through supported mechanisms.
- [ ] Evaluate local electrical interventions and barrier disruption against conserved ion/charge balances and available human wound observations.
- [ ] Add cell migration/healing response only where an explicit supported law and target dataset exist; annotate organism/tissue transfer and uncertainty.
- [ ] Show membrane voltage and extracellular field as different selectable fields with time-linked provenance. Neural sensory activity must not be relabeled as non-neural bioelectric state.

**Verification:** equilibrium/perturbation tests, finite-domain boundary tests, solver refinement, source-model parity and holdout comparison in `scripts/verify_body_skin_electric.py`.

### I12 — Whole-body expansion, heart and remaining systems

**Concerns:** C02, C08, C12, C20–C22, C25–C27. **Depends on:** successful regional work packages and I1 coverage contracts.

**Files:** extend domain builders and `ihm/assembly/body.py`; add dedicated regional/system adapters only when their state ownership is defined; extend `scripts/audit_body_coverage.py` and system-specific verification.

Apply the same evidence→domain→law→coupling→perturbation cycle to head/neck, thorax, abdomen/pelvis, upper limbs and lower limbs. Each region has its own acceptance report; anatomical continuity at region boundaries is mandatory.

| System | Concrete expansion beyond a visual organ | Required perturbation / observable |
|---|---|---|
| Cardiac/pericardial | Chamber partition, active myocardial material/fibers, valve behavior, conduction/activation and pressure-volume work | Preload/afterload or activation change → altered contraction and stroke output with circulation balance. |
| Respiratory | Bronchial regional ventilation, gas exchange/perfusion, chest/abdomen interaction | Resistance/compliance change → pressure, volume, regional ventilation and blood-gas response. |
| Renal/urinary | Renal perfusion/filtration/tubular and urine-storage ownership; spatial detail where needed | Perfusion or fluid load → filtration/urine/electrolyte balance. |
| GI/hepatic/pancreatic | Organ support/contact, motility, luminal storage/transport, absorption, metabolism/secretion | Defined nutrient delivery → transport, portal/systemic metabolites and endocrine response. |
| Endocrine/metabolic | Hormone transport, secretion, receptor response and feedback timing | Defined hormonal/metabolic challenge → target response and feedback with species balance. |
| Immune/hematologic/lymph | Cell/protein transport, regional inflammation, hematology/coagulation where supported | Local inflammatory/vascular challenge → tissue/lymph/systemic response with declared model domain. |
| Integumentary | Regional layers/appendages, sensation, perfusion, thermoregulation and bioelectric state | Pressure, heat, electrical or wound perturbation → appropriate coupled local response. |
| Reproductive/pelvic | Chosen profile anatomy, organ mechanics and reproductive/endocrine ports | Model-appropriate endocrine or mechanical challenge; variants tracked explicitly. |
| Cranial and special senses | Eye/ear/vestibular/chemosensory transduction, cranial nerves and effectors | Physical modality stimulus → correct pathway response and declared effector behavior. |
| Autonomic/visceral/enteric | Regional ganglia, visceral afferents/effectors, reflex ownership and control interfaces | Baroreceptor/visceral challenge → appropriate autonomic response without duplicate native control. |
| Connective/adipose/marrow | Region-specific support, compliance, storage and marrow physiology where in scope | Mechanical or metabolic perturbation with an explicit corresponding observable. |

“All body systems” requires entries in this execution ledger even when a system remains lumped. Increasing geometric detail must not silently imply increased physiological resolution. Full articulated movement and gait require a separate validated control/contact milestone; supine local movement does not establish them.

### I13 — Calibration, uncertainty propagation, temporal predictors and final acceptance

**Concerns:** C28–C32, C35–C36; applies continuously to I1–I12, with a final integrated gate.

**Files:** extend `ihm/assembly/certainty.py`, `temporal.py`, `ihm/temporal/`, calibration modules and `ihm/human.py`; create `scripts/verify_body_crossscale.py`, `scripts/verify_body_experiments.py`; update `scripts/build_workbench.py` and `scripts/verify_all.py` only after dependencies are integrated.

- [ ] Perform sensitivity and structural/practical identifiability analysis before fitting large parameter sets. Fit only supported parameter combinations; preserve unresolved model-form alternatives.
- [ ] Separate subjects/conditions between fitting and validation; propagate correlated source uncertainty and registration uncertainty into scenario observables.
- [ ] Fit reduced models against reference solvers over declared state/input domains; report worst-case and held-out error, energy/passivity behavior and fallback triggers.
- [ ] Integrate IBM spectral-window state where supported; retain time-domain event/contact handling where necessary. Verify representation conversion, delay phase, nonlinear round-trip error and overlap/window boundary causality.
- [ ] Keep finite Laplace descriptors, mechanistic resolvents and identified predictors as distinct output types. Do not infer causal coupling from coherence alone.
- [ ] Test all six end-to-end scenarios, source/code mutation rejection, checkpoint restart, deterministic seeds, LOD invariance and intervention traceability.
- [ ] Publish the final supported-query/coverage report and explicit unsupported requests. The app must expose the active fidelity/validity domain when it affects interpretation.

**Verification:** targeted package checks, empirical holdout reports, then `.venv/bin/python scripts/verify_all.py --app` and the new integrated experiment suite. Native closed-loop experiments are mandatory for capabilities claiming feedback; recorded replay is insufficient.

## 5. Proposed numerical and usability gates

These are starting engineering thresholds to evaluate during E2/E5/E6. They are not fitted human accuracy claims and may be tightened after reference evidence is available.

| Gate | Proposed criterion |
|---|---|
| Identity / ownership | No dangling required entity/port references or duplicated authoritative state owners. Every unresolved relation explicitly enumerated. |
| Material geometry | No inverted active volume elements; orientable active boundaries; intersection/contact exceptions individually classified. |
| Conservation | For each controlled test, normalized storage-versus-integrated-exchange residual ≤1e-5 using a predeclared quantity-specific reference scale; also report dimensional residuals so tiny fluxes are not hidden. |
| Contact | Maximum penetration below 1% of local element/feature scale and decreasing under tolerance/refinement study; boundary/friction work accounted for. |
| Time/space convergence | Changes in target observables below 1% under final refinement for reference numerical tests; failure triggers further refinement or a narrower validity domain. |
| Reduced-model fidelity | Target observable error ≤5% over a declared tested input/state domain, with absolute tolerances for near-zero quantities. No claim outside that domain. |
| Causal neural propagation | No response before the declared path/synaptic delay; delay converges under temporal refinement; selective block removes only the specified path's contribution. |
| Rendering | Agreed interactive preset p95 frame time ≤33 ms on the measured machine; source truth/solver accuracy not reduced silently to meet it. |
| Layout | No document scroll; three independent scroll regions when visible; all eight panel combinations; camera/selection/checkpoint preserved. |
| Empirical validity | Scenario-specific metrics and acceptance bands frozen from data/measurement variability before fitting. Numerical thresholds above cannot substitute for these. |

Concrete test specifications:

```text
Given an equilibrated layered-skin reference patch and a specified contact tool,
when an indentation pulse is applied and released,
assert positive valid elements, bounded penetration, converged reaction force,
receptor-class-specific response, correct afferent arrival and accounted work.

Given a coupled supine respiratory model with an identified chest-wall parameter,
when that parameter or external chest load changes,
assert altered pressure/flow/volume and phase-consistent chest/diaphragm motion;
reject a run whose physiology is unchanged merely because motion was prescribed.

Given a named muscle, motor pathway, tendon and joint,
when the motor pathway is blocked under otherwise identical drive,
assert loss of its evoked activation/force contribution and preservation of
unrelated pathways; compare joint torque against the moment-arm/force relation.

Given a fixed body version, vascular/hair seed and simulation checkpoint,
when the camera zoom or display triangle/instance budget changes,
assert identical underlying topology, state, aggregate transport and follicle IDs.
```

## 6. First deliverables and review sequence

1. Review this concern register and roadmap, including generic subject, supine scope and the proposed six acceptance scenarios.
2. Authorize exploration separately from implementation if desired. E0–E6 then produce the source coverage, compatibility reports, backend benchmarks and final numerical/empirical criteria.
3. Implement I1–I4 plus the fixed workspace; demonstrate a mechanically meaningful regional experiment.
4. Implement I5–I7; demonstrate chest mechanics with feedback and the skin→nerve→brain/spinal→muscle→mechanics loop.
5. Implement I8–I11 and expand regions/systems through I12, with calibration and uncertainty work active throughout.
6. Accept only claims supported by I13's integrated and empirical evidence. Save every failed experiment and source conflict that changes the model design.

No new implementation, simulation, installation, deployment or source-data build was performed to produce this plan. Existing experimental code is preserved for later review, and implementation remains paused pending the user's direction.
