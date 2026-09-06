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


## Current integration checkpoint (2026-09-06, continuing wave)

All broad objectives above remain open. These accepted increments narrow the
implementation gaps without closing whole-body validation:

| Area | Current executable evidence | Next unresolved integration |
| --- | --- | --- |
| Signed mechanical energy | Native chemical/heat/signed-work ledger; corrected shared-donor library; six signed whole-engine boundary steps and native water intake accepted (`e998b4b`); latched signed sessions reject ordinary-step bypass (`efa4fc8`) | Actual mechanical startup asks −36.49 W versus allowed −8.69 W and is correctly rejected; need supported reference and causal energy-supply feedback |
| Support | The 92-muscle constrained solve reached objective 71.82 with support checks passing but max acceleration 5.10 rad/s². Source-backed 98-muscle model loads and passes seven-pose geometry tests (`a145530`); fresh seam-omitted static diagnostic (`3e330d9`) | Old pose under 98-muscle physics has max acceleration 54.17 rad/s² and is not reusable as equilibrium. One unresolved contact cell is omitted; no old responses reused. Native Thelen passive-energy getter uses dimensional fiber length where normalized length is required; an observational correction is retained separately from raw values. Four original Arm26 wrapped paths violate length-gradient/force consistency. An isolated BRA-private cylinder candidate passes 208 native work checks and 9,984 unchanged observations of the other 96 muscles (`a42c62f`); BIClong ellipsoid routing remains nonstationary and no production wrap is promoted. Actual-force solve `ec51c33` failed with 139.18 N support imbalance and max acceleration 88.53 despite lower numerical cost. Native inverse-mass/sign checks pass; a new hard-support solve (`support-physical-root-nayktxy6`) reduces max acceleration 54.17→5.688 while preserving support/gauge/constraints. Equilibrium remains failed; isolated wrap corrections continue. No supported reference, forward proof or startup promotion |
| Skin–brain–motor | Direct skin sensor identity and causal latency/block tests; immutable IBM candidate and verified API selection (`7a6c83e`). All six new lumbar muscles pass controller/replay tests; actual 98-muscle regional factory initializes (`ae26353`) | No accepted whole-body reflex trajectory or default source/model promotion. Source geometry, inferred cortical recruitment and peripheral mapping remain distinct; factory readiness is zero-time, not integrated motor validation |
| Fine vessels | Muscle and curved kidney patches available through bounded API and monitor; actual 122-edge kidney HTTP response, whole-kidney containment (`82f41b5`, `0fa42a9`) | Human-evidence-conditioned skin plexus/papillary-loop patches now bind exact exterior faces in all eight engineered lower-leg territories (`a294ca9`), with 165,587 primary source bytes retained. Cortical kidney location, physical dermal containment and macro-vessel attachment remain unresolved; no full regional perfusion feedback or all-organ coverage |
| Lymph | Existing native Skin graph retains separate inventories, pressure response and sweat ownership. Eight exclusive engineered lower-leg masks use an explicit exterior component. Configured named regions preserve aggregate circuit state and pass local drainage/release fixture (`127a5cb`). Separate native 28-species transport and Albumin phases pass conservation, local pressure response/release and donor exhaustion (`fd5020a`) | Area/thickness/EC-fraction weights are explicit priors, not measured drainage volumes. Full-patient configured-region/transport acceptance, live contact-to-drainage binding, complete lymph topology and long-run validation remain open |
| Intake/GI | Native consumption adds 30 g to actual mechanical owner at common endpoints (`beac491`); opt-in shared runtime. Conditional colon and finite epithelial models. Actual mapped native fixture transports all 61 loaded substances through distinct lumen/fecal stores, commits NextVolume once and preserves settled state/epoch (`0a24655`) | Full-engine ordering and activation remain unverified. Original GI codecs reproduce omitted oral-drug state, nested writer crash and cleanup defects; isolated full-TU repair passes native roundtrip, distinct metabolized/excreted totals and ownership cleanup (`7c31bfb`). Shared-engine composition and continuation remain pending. Native wall current/bicarbonate, acid–base, long-run glucose, excretion-to-mechanics and intake beyond 0.5 kg remain open |
| Bioelectric integument | Finite three-ion epithelial inventories and two membrane capacitors; primary human TEP data acquired; bounded conservation/energy tests (`2e2a087`) | Exact exterior-area shadow partition now preserves existing IC/EC species and unresolved complements, with separate apical film and explicit capacitor/heat owners (`6151cd3`). Native exchange exclusion, epoch-atomic commits, calibrated ion pathways, anatomical layers and wound/regeneration predictions remain unresolved |
| Mass and cervical mechanics | Opt-in native instance mass with source-state inertia, momentum/energy accounting, replay and post-mutation rollback accepted in 92-muscle model (`874e5f0`). Cervical source recipe preserves mass/moments (`9b9cfd6`) | Cervical recipe is not loadable anatomy: C1 registration residual 36.2 mm, missing exact C7 mesh and unresolved shoulder/contact mapping. Restricted human hip-capsule research law and primary evidence (`94ef078`) apply only to one measured pose/load slice; they do not justify a general three-DOF force or ±40° behavior. No native activation or calibrated whole-body ligament/joint-limit behavior |
| Respiratory mechanics | Native external-pressure/airflow and source PV-work checks pass. Executable source thorax includes complete material mass/cross terms, moving attachments and a geometric cavity/Jacobian (`1a3ec1b`) | Native display chest remains massless; the new thorax is not coupled to native recoil or IBM drive. Its cervical-reduced mass basis needs explicit composition with actual torso. Correct body-frame dynamics and SO(3) pose rates pass energy/momentum rate checks (`1de3fa8`); a 50 ms midpoint trajectory passes refinement and the cervical+thorax composition plan closes parent mass/moments (`a1aed47`). Native coupled inertia and physical locked-joint reactions remain open. No calibrated full thoracic constitutive/contact model. Corrected sleep variant now passes 32 native respiratory observer/control steps with exact full-state parity (`5c2de08`); this does not establish save/reload continuation |
| Hair and clothes | Exact source eligibility retains 94 physical guide candidates and records excluded inner roots. All 22 residual native tensors and combined momentum/energy pass initialization (`396eb6d`). Two actual 0.1 ms clamp-coupled intervals, replay and rollback pass (`b3291c5`) | No collision occurred; 20 ms readiness, native total-energy closure and live integration remain unverified. Original hair population and ambiguity retained. Whole-body garments, cloth contact/CCD/self-contact, realistic clothing and source-consistent viewer projection remain open |
| Native state persistence | Exact source audit (`f0fa8ed`) identifies missing sleep state, vascular region metadata, oral-drug inventories and Tissue burn histories. Circuit IO native roundtrip passes all five regions (`3d79ed0`). Sleep owner persistence and paired-DSO observer/control acceptance pass (`804eee4`, `5c2de08`); isolated GI repair passes (`7c31bfb`), and Tissue burn-history repair remains source-staged | Compose and validate corrected native owners/libraries, then exact state roundtrip and next-tick continuation. Missing historical state cannot be silently restored; explicit fresh initial conditions must be labelled. Existing short output-port checks are not full-state equivalence |
| Skin source geometry | Two dominant opposing shells explain raw area 3.5026 m² versus selected exterior proxy 1.78046 m². Future layer generation uses explicit exterior support (`99685de`); contact thickness metadata migration checked (`803ddd3`) | An isolated migration candidate recalculates all proxy masses/inertias and 3,341 damping links while preserving shell thickness (`4481fbf`); live canonical artifacts remain unchanged. Finish dependent consumer acceptance before promotion; open boundaries, self-intersections and physical surface exclusivity remain unresolved |
| Visualization | Served mass, finite Laplace, intake and bulk Skin Vm monitors plus muscle/kidney/skin vessel inspection; current bundle `index-f-cYJevs.js`. Native Vm observational fixture reads −82.822 mV without advancing or changing 1,248 tissue fields (`871ef3a`). Commands use actual pressure topology (`3a1baf4`, `5de64bb`); actual Skin HTTP request passes with 25 edges and zero native volume allocation | Stable live body and browser acceptance remain open. Research modes are opt-in; passing component/factory checks does not make default body fully coupled. Remaining cursor/environment/force/gimbal and whole-body soft dynamics requirements stay open |

The original cardiovascular reader ran only within a drug-resistance branch;
a separate routine observation patch passes exact output parity. A transferred
vascular-law experiment failed because LoadState omitted region labels and left
regional resistance caches empty. IO repair passes its native codec test and the repaired-baseline experiment now
passes actual signed resistance/flow response, exact zero/inactive parity and
local release (`a2ed766`). Configured regional composition remains pending. Reader
access and conditional engineering laws are not empirical calibration.

Current parallel work includes supported equilibrium, thoracic dynamics,
hair/body contact, complete native state persistence, configured regional
exchange and GI ownership. The latest local server refresh receipt remains
`artifacts/workbench-skin-refresh-verification.json`; this wave has not launched
a browser or promoted the experimental physiology/mechanics into the default
viewer. Native compilation and runtime fixtures use one explicit shared slot;
`artifacts/heavy-job-queue.json` records the manual grant and waiting tasks for
recovery. All broad objectives remain in scope.
