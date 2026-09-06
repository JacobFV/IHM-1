# IHM-1 remaining objectives and acceptance concerns

Consolidated from the user's requests, including the latest priority on full
mechanical feedback, brain-driven sensory–motor connectivity, whole-body fine
blood/lymph/integument and digestion. This is the continuing objective register,
not a claim that its unchecked items are implemented. Existing isolated tests,
source acquisitions or interfaces do not close a whole-body objective.

## 1. One definite implicit human

- [ ] One coherent generic body with shared geometry, identities, coordinates,
  clocks, units and state ownership; the app operates that body rather than a
  collection of mutually exclusive source models.
- [ ] Engineer a unification of heterogeneous sources, synthesizing missing
  anatomy and interactions with explicit assumptions. Do not wait for one
  complete atlas or treat incompatible donor atlases as one measured person.
- [ ] Materialize mechanical, physiological, neural, transport, electrical and
  visual predictors from the same substrate. Preserve the possibility of later
  personalization without calling the current generic body the user's twin.

## 2. Complete anatomical and material construction

- [ ] Bones with meaningful mass/inertia and joints; muscles with origins,
  insertions, paths, wrapping, fiber directions and force laws; tendons,
  ligaments, cartilage, fascia, fat, marrow and connective interfaces.
- [ ] Organs with walls, cavities, lumens, internal structure and realistic
  material properties, including smaller glands, ducts, valves and structures.
- [ ] Resolve overlapping source representations and exclusive material/mass
  ownership. Repair or conditionally reconstruct invalid source surfaces while
  retaining their original evidence. Of424 muscle bulk candidates,180 currently
  self-intersect; the244 passing screening are not yet certified tissue volumes.
- [ ] Complete upper-body, neck, hands/fingers and feet/toes actuation and
  articulation. The original80 native muscle model is lower-limb-only.

## 3. Whole-body mechanical feedback — immediate priority

- [ ] Actual articulated rigid and deformable soft-body dynamics with equal
  and opposite forces, moments, contact, friction, sliding and load transfer.
- [ ] Solve gravity-supported supine equilibrium on a bed, then changing
  posture and ground support. Initial contact or an assumed balancing preload
  does not establish equilibrium.
- [ ] Connect skin, fat, muscles, organs, skeleton and garments mechanically;
  no independent objects passing through one another or being teleported by
  recorded motion. Resolve body self-contact and body/environment contact.
- [ ] Close mechanics↔physiology feedback: breathing resistance, tissue
  compression/perfusion/lymph drainage, muscle work/metabolism and effects of
  physiological state on mechanical behavior. Regional successful ports are
  starting points, not full mechanical coupling.
- [ ] Verify momentum, angular momentum, energy/work, contact and timestep/
  spatial convergence at meaningful body scale.

## 4. Brain-driven sensory and motor integration — immediate priority

- [ ] Reuse and periodically update pinned IBM-1 imports with compatibility,
  provenance and numerical parity. Do not modify or duplicate its concurrent
  brain microstructure acquisition/training.
- [ ] Complete sensory, motor and autonomic peripheral pathways, with named
  targets, modalities, conduction delays and appropriate effector laws.
- [ ] Skin senses touch, pressure, vibration, temperature and relevant noxious
  stimuli; muscles/tendons/joints supply proprioception and load feedback.
- [ ] Actual mechanics→receptors→afferents→brain/spinal processing→efferents→
  muscles→mechanics loops, with selective blocks, release and descending drive.
- [ ] Coordinate diaphragm/intercostals and other respiratory muscles,
  posture, balance and whole-body movement. Distinguish native brainstem
  regulation, preserved IBM dynamics and inferred recruitment/control mappings.
- [ ] Extend special senses and visceral feedback: vision, hearing,
  vestibular balance, smell/taste and internal organ sensing where required.
- [ ] Autonomous movement learning remains future work; physically valid
  controlled movement must not depend on solving learning first.

## 5. Breathing, heartbeat and circulation

- [ ] Neural respiratory drive produces muscle force and thoracic/diaphragm
  motion; resulting pressure, airflow, lung volumes and gas exchange feed back.
  Include chest/abdominal loading and appropriate airway/alveolar structure.
- [ ] Heart electrical activation, myocardial mechanics, chambers, valves and
  pumping interact with preload, afterload and coronary supply.
- [ ] Couple visible heartbeat/chest motion to solved live state; do not use
  decorative animation as proof of integration.
- [ ] Account for pressure–volume work and avoid duplicate passive recoil,
  storage or metabolic demand when connecting engines.

## 6. Fine blood vessels and blood state — immediate priority

- [ ] Connected arterial, arteriolar, capillary, venular and venous structure
  across the body, including organ-specific beds and macro/micro boundaries.
- [ ] Blood pressure, flow, compliance, rheology, oxygen/CO2, nutrients,
  electrolytes, hormones, proteins, cells and relevant coagulation state.
- [ ] Conditional organ-specific synthesis from acquired imaging, histology,
  autopsy and derived graphs. Constrain geometry by perfusion, exchange area,
  transit times, oxygen extraction and pressure losses—not appearance alone.
- [ ] Measured versus transferred versus synthesized topology remains explicit.
  Sparse images cannot recover individual unobserved capillaries. Brain,
  muscle, liver, kidney, lung and skin require different structural priors.

## 7. Lymph, interstitium and immune transport

- [ ] Initial lymphatics, collectors, valves, nodes and return to circulation
  distributed throughout the body.
- [ ] Fluid/protein/cell exchange, oncotic/osmotic pressure, drainage, edema,
  muscle pumping and compression response use one physiological storage owner.
- [ ] Connect regional geometry and local mechanical loads to those processes;
  a uniform lumped Skin compression test does not complete local drainage.
- [ ] Integrate lymphoid organs, immune/inflammatory processes and relevant
  vascular/skin/endocrine interactions at declared model resolution.

## 8. Integument and bioelectricity

- [ ] Layered epidermis/dermis/subcutis, vessels, nerves, follicles, sweat and
  sebaceous glands, regional material properties, barrier and water/heat flux.
- [ ] Membrane and transepithelial voltages, ion reservoirs, pumps/channels,
  intercellular coupling and mechanically/chemically dependent behavior.
- [ ] Address the Levin-inspired interest through explicit hypotheses and
  testable bioelectric interactions. Bulk Nernst estimates or illustrative RC
  networks are not established keratinocyte voltages or regeneration models.

## 9. Digestion and nutrition — immediate priority

- [ ] End-to-end intake, oral processing/swallowing, esophageal transport,
  gastric mechanics/chemistry, intestinal motility, secretion and absorption,
  colon/water handling and elimination.
- [ ] Liver/portal circulation, bile/gallbladder, pancreatic digestion and
  endocrine control, nutrient storage/use, renal losses and appropriate enteric/
  autonomic feedback are integrated with the same body state.
- [ ] Timed food, water and other intakes change actual gastrointestinal,
  vascular, hormonal, thermal and mechanical states.
- [ ] Resolve observed integrated glucose and acid–base failures, sodium/
  counterion handling and missing meal thermogenesis; retain failed runs.
  Local mass-integrity corrections do not establish healthy homeostasis.
- [ ] Account for intake/excretion changes in total mass and mechanical inertia,
  rather than permanently distributing newly ingested mass across all segments.

## 10. Other major and minor systems

- [ ] Renal filtration, tubular handling, electrolyte/acid–base control, urinary
  storage/voiding and pressure interactions.
- [ ] Endocrine axes and relevant hypothalamic/pituitary, thyroid/parathyroid,
  adrenal, pancreatic and gonadal control.
- [ ] Hepatic metabolism/detoxification, hematology/marrow, spleen/immune state,
  adipose and connective-tissue physiology.
- [ ] Reproductive, pelvic and urogenital systems receive appropriate anatomical
  and mechanical integration as part of whole-body coverage; a localized example
  must not absorb effort intended for the entire body.
- [ ] Thermoregulation, sweating, environment/clothing heat exchange, fluid
  balance and longer-term homeostasis across the integrated system.

## 11. Realistic clothing

- [ ] Complete usable garments worn by default and independently switchable.
- [ ] Physically computed drape, stretch, bending, seams, openings, friction and
  contact keep clothes on the body during motion without penetration or hidden
  anatomy. Local cloth-panel tests are insufficient.
- [ ] Clothes affect and respond to soft tissue, movement, heat and moisture.
  Source coefficients must match fabric, region and conditions or remain priors.

## 12. Real hair

- [ ] Scalp and region-appropriate body hair as independent rooted strands,
  with believable density, thickness, lengths, orientation and appearance.
- [ ] Particle/rod dynamics, bending, damping, skin motion and relevant hair/
  body/garment/self-contact. Interpolated rendering retains its guide mapping.
- [ ] Efficient solver scheduling/off-thread execution without quietly turning
  a requested physical feature into painted skin or decorative animation.
  Current dynamics are opt-in because performance is not yet reliable.

## 13. Movement, environment and direct interaction

- [ ] The assembled body moves through its environment under solved forces:
  supported motion, transitions, reaching, crawling and ultimately walking.
  A separate planar walking model does not complete canonical3D movement.
- [ ] Choice of bed/floor/other environments and physical objects, with actual
  body/object contact and reaction forces.
- [ ] Select, translation/rotation gimbal and free-force interaction, correct
  material grab points, force/torque limits and visible applied-force readouts.
- [ ] A user poke/drag propagates through mechanics, sensing and physiology
  rather than affecting an isolated view-only object.

## 14. One live, 3D-first application

- [ ] Retain the implemented fixed viewport, clean left View sidebar and
  independently scrolling right monitor-card stack; no bottom pane/header/page
  scrolling or return to sprawling source-choice-driven workflows.
- [ ] Finish configurable hospital-like monitors: multiple independent signals,
  physiological input/setpoint controls, clickable body stimulus diagram,
  intake schedules and searchable overlapping categories in +New Pane.
- [ ] Show one coherent live state. Replays, experiments and source inspection
  must not overwrite or falsely masquerade as that live body's state.
- [ ] Keep anatomy selection, uncertainty, inputs, graphs and stimuli consistent,
  with clear physical units and useful local zoom for microscopic structures.

## 15. Temporal materializations

- [ ] Extend IBM-style temporal representations where mathematically appropriate,
  including finite Laplace spectra of actual body variables and causal predictor
  materializations when supported.
- [ ] Preserve clocks, units, sampling, window length, damping and numerical
  limits. Finite-window signal spectra are not physiological poles or an
  identified causal transfer model.

## 16. Data, calibration and precision

- [ ] Continue acquiring substantial real heterogeneous data, including raw
  foundation products when preprocessed atlases are inadequate. Record acquired
  bytes, versions, licenses, specimen identities and processing as work proceeds.
- [ ] Fit and validate interaction coefficients against appropriate observations;
  transfer source parameters with conditions and uncertainty, not blanket claims
  of human calibration. Avoid double-counting correlated evidence.
- [ ] Track source resolution, geometry registration, biological variation,
  inference/model discrepancy, numerical error and display approximation
  separately, through every explicit materialization and selected structure.
- [ ] Preserve high-resolution truth, including10million-plus triangles when
  useful. Lower display resolution only through documented projections/error
  limits, not by discarding the underlying data.

## 17. Integration acceptance and operating constraints

- [ ] Demonstrate integrated perturbations: load→breathing/gases; touch→neural
  response→muscle/force; compression→perfusion/lymph; intake→absorption→blood/
  hormones/stores; motion→work→physiology; clothes↔whole-body mechanics.
- [ ] Validate corresponding releases, selective blocks and matched controls,
  conservation, positivity, equilibrium and temporal/spatial convergence.
- [ ] Run long enough to expose failures. Short finite trajectories, passing
  tests, adapter existence and acquired datasets are distinct from whole-body
  predictive validation.
- [ ] Maintain reproducible snapshots, source receipts, actual/uncertain command
  outcomes and crash recovery. Never claim native rollback that is not exact.
- [ ] Parallelize independent research/code work while bounding local native
  computation, builds, downloads and browser load. Preserve IBM-1 training,
  single-thread numerical settings and low-priority execution on this machine.

Immediate work order: finish the shared articulated/neural/native runtime and
its real feedback ports; connect full-body material/contact and garments;
expand regional blood/lymph/skin and digestion integration; then expose those
same states and controls through the unified viewer. All other concerns remain
in scope and must not disappear when a narrower implementation is completed.


## Current integration checkpoint (2026-09-05, continuing wave)

All broad objectives above remain open. These accepted increments narrow the
implementation gaps without closing whole-body validation:

| Area | Current executable evidence | Next unresolved integration |
| --- | --- | --- |
| Signed mechanical energy | Native chemical/heat/signed-work ledger; corrected shared-donor library; six signed whole-engine boundary steps and native water intake accepted (`e998b4b`); latched signed sessions reject ordinary-step bypass (`efa4fc8`) | Actual mechanical startup asks −36.49 W versus allowed −8.69 W and is correctly rejected; need supported reference and causal energy-supply feedback |
| Support | Explicit force/moment-constrained native solve with immutable archived dependencies and assembly-error equivalence evidence; objective reduced from millions to ~75 while force/gauge checks pass | Latest maximum acceleration 4.77 rad/s²; descent slows near hip rotation bound (`247e82b`). Investigating rejected steps and physical joint-limit ownership; no accepted equilibrium, forward proof or startup promotion |
| Skin–brain–motor | Direct native sensor identity/registration and causal block/latency tests; immutable current IBM candidate shared by regional brain, receptor and checkpoint identity; actual candidate+regional+sensor factory initialization (`c3a8814`) | No accepted whole-body reflex trajectory; default source promotion and evidence-driven new microcircuit/vascular imports remain separate. Cortical recruitment and body mapping remain explicit priors |
| Fine vessels | Muscle and curved kidney patches available through bounded API and monitor; actual 122-edge kidney HTTP response, whole-kidney containment (`82f41b5`, `0fa42a9`) | Cortical kidney location and macro-vessel attachment unresolved; no full regional perfusion feedback or all-organ coverage |
| Lymph | Actual three-region native Skin circuit and signed/embodied initialization; separate inventories, local pressure response and sweating ownership. Exclusive source-skin territory validation and eight unregistered lower-leg annotation slots (`036f03a`) | Native engineering fractions still have no anatomical territory binding; no actual local contact-to-drainage feedback, complete lymph topology or long-run acceptance |
| Intake/GI | Authoritative consumed-intake receipts; actual native consumption adds 30 g to a registered mechanical owner at common endpoints (`beac491`); shared-runtime opt-in integration (`b5ed55f`). Finite epithelial transport with marginal-affinity bounds; conditional human colon perfusion fits and conservative transit (`128a3fc`) | Whole-body signed-step acceptance, intake beyond 0.5 kg validated payload range, excretion and internal regional mass redistribution remain open. Colon/membrane models not native inventory owners; bicarbonate/current, acid–base and long-run glucose behavior unresolved |
| Bioelectric integument | Finite three-ion epithelial inventories and two membrane capacitors; primary human TEP data acquired; bounded conservation/energy tests (`2e2a087`) | Native inventory/energy ownership, calibrated ion pathways, anatomical layers and wound/regeneration predictions remain unresolved |
| Mass and cervical mechanics | Opt-in native instance mass with source-state inertia, momentum/energy accounting, replay and post-mutation rollback accepted in 92-muscle model (`874e5f0`). Cervical source recipe preserves mass/moments (`9b9cfd6`) | Cervical recipe is not loadable anatomy: C1 registration residual 36.2 mm, missing exact C7 mesh and unresolved shoulder/contact mapping. No native activation or calibrated ligament/joint-limit behavior |
| Respiratory mechanics | External chest pressure affects native airflow/gases. Exact source-flow pressure–volume work observer passes isolated native circuit checks (`6cf3b29`) | Chest expansion is massless display deformation; no articulated diaphragm/ribs or applied IBM breathing drive (`3779682`). Whole-engine work observation and physical thoracic/recoil/contact ownership are next |
| Hair and clothes | Worker hair display dynamics; Python beam/contact/follicle reaction prototype with finite-body momentum, energy, replay and synchronization tests (`66de7ea`); exact source-face and actual segment-inertia factory preparation | Hair prototype not attached to live native body/viewer; mass debit and owner registration unresolved. Canonical hair proxy mass is descriptive, not an independently integrated native store. Whole-body garment/contact/CCD/self-contact and realistic clothing remain open |
| Visualization | Served mass monitor plus kidney inspection, finite Laplace spectra and intake monitor; current bundle `index-C5Q8H2-C.js`; backend exposes strict regional/intake mass mode validation with no active body required | Stable live body and browser acceptance remain open. Research modes are opt-in; passing component/factory checks does not make default body fully coupled. Remaining cursor/environment/force/gimbal and whole-body soft dynamics requirements stay open |

Cardiovascular signed-demand reader audit needs particular care: in the held
source it runs only inside a drug-resistance-change branch. Rest-fixture count
zero is reported honestly; it does not establish routine demand-driven vascular
feedback. Endocrine and both observed nervous readers consume signed effective
demand in the native fixture. Tests of reader access are not empirical validation.

Current parallel lanes include support/joint-limit diagnosis, physical thoracic
work interfaces, anatomical hair factory and mass-ownership preparation, and
conditional GI data/calibration. Previously accepted source and short-fixture
results remain distinct from integrated predictive validation. The latest local
server refresh receipt is `artifacts/workbench-mass-refresh-verification.json`;
no browser or native body was launched by that HTTP check.
Native compilation and runtime fixtures use one explicit shared heavy-job slot.
