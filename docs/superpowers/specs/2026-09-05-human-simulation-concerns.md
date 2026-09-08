# Human simulation: requirements and concern register

Status: planning proposal, September 5, 2026. Implementation is paused at the user's explicit request. This document and its companion roadmap are planning deliverables only. Earlier embodied-human drafts and commits are evidence to review, not authorization to proceed.

## Intended result

One generic adult human with coherent anatomy and physically meaningful interactions, initially resting supine on a bed. Its explicit visual, mechanical, vascular, electrophysiological and predictive representations must refer to the same underlying body. Detailed structures may be synthesized, but their assumptions, uncertainty and consequences must remain inspectable. The viewport is the primary workspace; information panels can appear/disappear and scroll independently.

The acceptance question is not merely whether a structure is present. It is whether a perturbation reaches the correct structures through the correct mechanisms, with observable consequences, conserved exchanges, justified parameters and a declared resolution/validity domain.

## What the audit actually establishes

The inspected canonical artifacts contain 2,408 entities, including 425 muscle surfaces, 639 arterial and 395 venous surfaces, 146 nerve surfaces, four cardiac cavity surfaces, five registered lung lobes, three synthesized skin layers and a registered lymphatic graph. Counts are representation counts, not counts of complete independent anatomical organs.

The recorded 30-second materialization contains 301 frames. Only seven entities have nonidentity transforms: five lung lobes and two cardiac cavities. This does not establish chest-wall motion, myocardial contraction, pleural contact or movement of all surrounding tissues.

The mechanical artifact contains 443 muscle actuators and 24 tendon/ligament paths. Eighty muscle paths inherit native OpenSim parameters; 363 use synthesis priors. Its runtime uses affine tissue elements and constrained bone orientations. A connected support graph is not a complete articulation/contact model. A total-mass normalization is not evidence of correct individual organ masses or absence of overlapping material volumes.

The brain artifact contains 80 populations and 640 inferred connections. The adapter executes selected preserved IBM neural functions. It does not execute the complete IBM materialization, spectral-state, inference and uncertainty runtime. IBM itself contains more extensive receptor, effector, afferent/efferent and spectral mechanisms that require a proper compatibility audit.

Earlier peripheral work produced an artifact with 52 schematic nerves, 10 relays, 16 receptor patches, 215 supported bindings and 252 unsupported targets. Its target inventory includes muscle actuators and additional paths; those numbers must not be presented as an anatomical muscle coverage percentage. This prototype is absent from the inspected recorded body frames. Its simple stimulus-gated motor readout is not an established sensorimotor policy or reflex circuit.

The canonical anatomy includes bulk head-hair and pubic-hair surfaces, but no demonstrated body-wide follicle/shaft model. The 1,034 vascular surfaces are not evidence of a connected perfused microcirculation. The lymphatic graph's source topology does not supply calibrated vessel radii, valves or drainage.

Source files also contain early fixed-viewport changes and a reduced respiratory prototype. Some earlier work completed during this audit. This planning pass has not integrated, rebuilt, deployed, reverted or otherwise modified those implementations. Their review and reconciliation belong to the first future work package.

Audit anchors: `ihm/assembly/body.py`, `mechanics.py`, `brain.py`, `peripheral.py`, `respiration.py`; `scripts/build_body_mechanics.py`, `build_body_peripheral.py`; `app/src/main.js`, `style.css`, `state.js`; `data/derived/canonical/{anatomy,mechanics,brain,peripheral,trajectory}.json`; IBM's `ibm/materialize/`, `ibm/runtime/`, `ibm/processes/{transduction,effector}.py` and `ibm/topologies/{afferent,efferent}.py`.

## Concern register

“Absent” below means not demonstrated in the integrated canonical execution inspected here. It does not mean that no source code or regional example exists anywhere in the corpus.

| ID | Concern | Current gap / failure mode | Required evidence or acceptance result |
|---|---|---|---|
| C01 | Scope and generic subject | Composite sources can silently mix size, sex, age and disease states. Generic is not personal. | One versioned reference profile; explicit source-to-target transforms, parameter transfer rules and anatomical variants. |
| C02 | Anatomical completeness | Surface inventory omits many minor structures; named parents may be assemblies rather than tissue. | Region-by-region inventory of tissue, cavities, interfaces and unresolved structures with denominators and no double-counted parents. |
| C03 | Geometric integrity | Open surfaces, attribute seams, overlapping organs, inconsistent normals and atlas intersections. | Repair provenance; oriented material domains; collision/contact classification; positive valid volume elements in simulated domains. |
| C04 | Mechanical anatomy | Nearest-bone supports are not joints, ligament insertions, fascia, tendon sheaths or mesentery. | Typed attachment and interface graph; source-backed constraints; explicit unresolved attachment regions. |
| C05 | Inertial/material consistency | Bbox centers are not centers of mass; source scaling and mass normalization can hide impossible density distributions. | Disjoint material accounting, cavity/storage exclusions, COM/inertia estimates, fiber directions and density uncertainty. |
| C06 | Supine rest | Display rotation does not introduce gravity, bed contact, resting tone, prestress or hydrostatic redistribution. | Mechanically equilibrated body/bed state; compatible physiological initialization and recorded preload. |
| C07 | Soft tissue dynamics | One affine element per surface cannot resolve local indentation, sliding, folds or regional strain. | Volumetric/shell reference mechanics in active regions, measured material laws where available and error-controlled reductions. |
| C08 | All-parts interaction | Connectivity alone does not enforce nonpenetration, friction, attachment or exchange. | Explicit adjacency/interface registry; contact only where physically warranted; bidirectional force/work transfer. |
| C09 | Musculoskeletal behavior | Rotations constrained, generic attachments, frozen wrapping, simplified force law. | Joint axes/limits, compliant tendons, updated wrapping, activation/force-length/force-velocity laws, proprioception and force/torque checks. |
| C10 | Visible breathing | Native gas volume changes only lobes; chest/ribs/diaphragm/skin do not respond in the recorded run. | Coupled diaphragm, ribs, cartilage, chest wall, pleural interface, lungs and abdomen; source-clock deformation and measured kinematic targets. |
| C11 | Endogenous breathing | Prescribed volume is not a mechanically generated breath or closed respiratory feedback. | Respiratory drive produces muscle force and pressure, airway flow changes volume; one owner for each state; loading/compliance perturbations change ventilation. |
| C12 | Heart mechanics | A cavity-volume projection is not active myocardium, valve action or a 3D electromechanical heart. | Chamber partition, wall/fiber mechanics, activation timing, valves, pressure-volume work and circulation coupling; explicit reduced versus regional detailed scope. |
| C13 | Named vessels | Rendered branches lack demonstrated lumen centerlines, radii, connectivity and tissue ownership. | Validated arterial/venous networks, terminal territory bindings, lumen/wall distinction, boundary conditions and lesion propagation. |
| C14 | Microcirculation | No connected body-wide arteriolar-capillary-venular transport model. | Hierarchical tissue-bed model with arterial supply and venous return, exchange surfaces, vessel statistics and explicit regional capillary realizations. |
| C15 | Blood state | Bulk chemistry does not establish local hematocrit, oxygen extraction, pH, viscosity or transit distributions. | Advected species and local exchange, appropriate rheology, no duplicate blood storage, mass/species accounting across resolution boundaries. |
| C16 | Lymph/interstitium | Undirected chords and regional node associations do not solve drainage. | Interstitial fluid/protein balance, initial lymph uptake, directional collecting vessels/valves, pumps and venous return with uncertainty. |
| C17 | Integumentary structure | Uniform shell thickness omits regional layers, glands, follicles, fat, dermal mechanics and local perfusion. | Regional layered skin domains and attachment, appendages, barrier, sweat/heat/water transport, perfusion and innervation. |
| C18 | Body hair | Bulk hair meshes do not provide individual shafts, follicles, regional distributions or follicle receptors. | Follicle placement tied to material coordinates, regional terminal/vellus priors and exclusions, shaft geometry/dynamics at selected resolution. |
| C19 | Skin sensation | Prototype patches and generic gains do not establish whole-skin sensing or distinct receptor classes. | Spatial receptor populations, deformation/temperature-driven transduction, adaptation, conduction and observable target activity; parameter fits by modality. |
| C20 | Peripheral wiring | 146 nerve surfaces and 52 schematic prototype routes do not constitute all nerves/axons or correct innervation. | Cranial/spinal roots, plexuses, mixed nerves, branches, receptor territories and motor targets; correct laterality, fiber classes, delays and variations. |
| C21 | IBM integration depth | Selected equations omit IBM's request compiler, heterogeneous supports, temporal representations and evidence handling. | Versioned reusable IBM runtime integration, body-to-brain coordinate/port contract and donor-versus-integrated parity tests. |
| C22 | Spinal/brainstem control | Generic cortical stimulation is not a reflex pathway, motor program or respiratory rhythm generator. | Explicit reflex relays and inhibition, ascending/descending pathways, brainstem/autonomic ownership, motor recruitment and sensory feedback. |
| C23 | Non-neural bioelectricity | Neural activity, membrane voltage, extracellular fields and transepithelial potential are different states. | Distinct state definitions and units; ion channels/pumps, gap junctions, charge transport, boundary conditions and wound/perfusion coupling. |
| C24 | Healing/regeneration interpretation | A voltage field does not itself establish migration, morphogenesis or human regenerative prediction. | Explicit supported biological response law and target-specific evidence; species and tissue transfer recorded; unsupported outcomes excluded from claims. |
| C25 | Other senses and effectors | Eyes, ears and sensory anatomy are not functional visual, auditory, vestibular, olfactory or gustatory systems. | Physical stimulus→receptor→pathway→regional response contracts per modality, plus ocular/facial/visceral effector bindings. |
| C26 | Remaining body systems | Renal, endocrine, GI, hepatic, immune and reproductive code may remain lumped or separate. | Coverage ledger for structure, state, law, exchange, solver ownership and validation across every system; gradual replacement with conserved boundary exchange. |
| C27 | Multiphysics coupling | One-way replay cannot test feedback; multiple engines can double-own autonomic, blood, thermal or fluid state. | Single authoritative owner per state, conservative transfers, causal timing, rollback, stable two-way iteration and no duplicate controllers. |
| C28 | Temporal/spectral semantics | Finite Laplace descriptors are not IBM's spectral state or identified causal predictors. | Explicit representation type; phase/lag, units, bandwidth, window boundaries, anti-aliasing, delays and nonlinear-event handling verified. |
| C29 | Parameter calibration | Literature lookup is not target calibration; stiffness depends on loading, scale and condition. | Parameter cards with units, conditions, provenance, source dependence, identifiability, fit/holdout data and joint posterior or explicit unknown status. |
| C30 | Error and uncertainty | Local confidence metadata does not propagate geometry, parameter and model-form errors through predictions. | Correlated uncertainty across sources and couplings; scenario-specific predictive uncertainty and validity envelopes; no fabricated confidence. |
| C31 | Simulation versus drawing resolution | More triangles do not supply constitutive laws or microvascular topology. | Separate truth, simulation and display discretizations; conservative field transfer and persistent entity/material identities. |
| C32 | Compute feasibility | Uniform whole-body cellular/FEM/CFD detail is not justified by triangle budget; ARM64/GPU solver compatibility untested. | Measured memory, solve rate, stability and rendering budgets at escalating resolutions on this machine; documented degradation strategy. |
| C33 | Interaction UI | Scrolling document, hidden anatomy and unrelated signal traces prevent inspecting consequences. | Fixed viewport; independent toggleable/resizable left/right/bottom panels; camera state retained; frame-linked inspector and traces. |
| C34 | Honest visual scale | Tiny vessels and hairs may be physically subpixel; exaggerated diameter can masquerade as anatomy. | Physical scale by default, explicit magnification/exaggeration controls, local zoom/LOD and labels distinguishing structural versus solved fields. |
| C35 | Reproducibility/runtime integrity | Partial builds, mutable jobs, changed source code and ambiguous snapshots produce stale or nonreproducible results. | Pinned data/code/solver/seed/config, atomic artifacts, consistent checkpoints and readable stale-result rejection. |
| C36 | Verification versus validity | Numerical tests and screenshots can pass while biological mechanisms are absent. | Separate geometry, numerical, causal, empirical and usability gates with falsifiable perturbations and published coverage. |

## Architecture recommendation and alternatives

Recommended: a single implicit body with heterogeneous, task-specific materializations. Use a physiological backend for global compartment states; a dedicated mechanics backend for solids/shells/contact; a network solver for regional vascular transport; the actual IBM machinery for neural/receptor/effector materializations; and detailed skin/electrodiffusion domains where the scientific question requires them. The coupling contract, not the viewer, determines how these pieces form one body.

Alternative 1: keep extending the present affine/body-overlay runtime. This is useful as a test harness and coarse surrogate, but it cannot be the final answer for local contact, all-part interaction or calibrated sensory mechanics. Retain it only as a declared reduced materialization with a tested operating envelope.

Alternative 2: construct one uniformly resolved, monolithic whole-body cell/FEM/CFD model. This would require far more unresolved data and solver work before the first credible experiment. Do not select it without a specific question and a compute feasibility result. Regional high-detail reference models can later be coupled and expanded.

Candidate reference mechanics tools are FEBio for nonlinear/multiphysics biomechanics and SOFA for interactive mechanical/contact models. This is a benchmark shortlist, not a selection or installation decision. Their official descriptions establish relevant capabilities, not fit to this body or ARM64/GPU performance: [FEBio](https://febio.org/), [SOFA](https://www.sofa-framework.org/). SimVascular's [reduced-order tools](https://simvascular.github.io/documentation/rom_simulation.html) and [svZeroDSolver](https://simvascular.github.io/svZeroDSolver/) are candidates for testing regional circulation coupling. Preserve OpenSim's actual [path/wrapping semantics](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53090145) when evaluating muscle transfer.

The bioelectric scope must distinguish non-neural membrane-potential networks from sensory nerve signaling. Levin's [non-neural network modeling paper with Manicka](https://arxiv.org/abs/1912.04246) is a candidate conceptual/computational source. It does not supply whole-human calibrated coefficients or establish that the current body can predict regeneration.

## Proposed end-to-end scientific acceptance scenarios

1. **Resting supine breathing:** measured respiratory inputs/targets, equilibrated bed contact and preload; observable chest, diaphragm, rib, lung and abdominal movement; gas/mass and mechanical work accounting. Alter chest-wall compliance and airway resistance and observe changed pressure/flow in the coupled mode.
2. **Skin indentation and release:** a tool presses a defined skin area; skin/subcutis deform against deeper tissue, receptor classes respond with their own dynamics, afferent activity arrives after a path-dependent delay, and release gives the appropriate recovery/adaptation.
3. **Muscle stimulation and loading:** activate an identified motor pool; muscle/tendon force loads the correct joint and adjacent tissue; external load changes motion and spindle/tendon feedback. Nerve block selectively abolishes the appropriate response. No generic rule that any touch contracts an arbitrary muscle.
4. **Regional perfusion and drainage:** modify inflow or venous pressure in a specified tissue bed; capillary exchange, oxygenation, interstitial fluid and lymph return respond, with all storage and boundary flows accounted for.
5. **Integumentary electrical perturbation:** change a local epithelial boundary/channel parameter; membrane and extracellular fields, ion flux and any enabled cellular response evolve with charge/species checks and explicit tissue-specific evidence.
6. **Resolution and materialization equivalence:** the same experiment at two spatial/time resolutions, and in reference/reduced solvers, agrees within an error envelope established for the requested observables. Display LOD alone does not change the prediction.

These are acceptance experiments, not assertions that their complete mechanisms are currently implemented. Full-body coverage additionally requires every in-scope region/system to have an explicit representation and every declared interaction to resolve; passing a forearm or thorax test alone cannot establish that.

## Planning assumptions for review

Retain the current generic adult male as the first reference body, with explicit future variant support. Use supine rest and controlled local movements first; arbitrary locomotion remains a later scope gate. Permit offline high-fidelity reference solves and interactive reduced playback/control, clearly labeled. Preserve full source geometry and permit more than ten million displayed triangles when measurements justify it. Synthesis is allowed but must be conditioned on source constraints and must not be presented as measurement. No software installation, new simulation, model integration or app deployment is authorized by this plan alone.
