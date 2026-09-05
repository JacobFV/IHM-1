# IHM-1 implementation and evidence status — 2026-09-05

The workbench opens one clothed generic body at `http://127.0.0.1:8765/`. It supports computed whole-body replay and body-anchored garment/tissue contact, systemic, sensory and electrical materializations. The implicit substrate binds anatomical, physiological, electrical, population and microscale evidence with explicit source identities. It is not a fully calibrated digital human: global all-organ contact, complete peripheral connectivity and organ-wide microstructure remain incomplete. The research roadmap's broad gates are not marked complete by passing software tests.

## Implemented behavior and retained evidence

| Domain | Executable or acquired result | Scope of the result |
| --- | --- | --- |
| Whole body | 2,408 anatomical representations; regenerated 30-second native-driven record with 301 frames and 2,323 changing reduced transforms | Canonical frame, explicit source fits and assumptions; representation count is not a count of independently simulated organs |
| Clothing and soft tissue | Independently visible default tank/shorts; source-derived three-tissue penile material domain with 5,796 tetrahedra at 4 mm, and a 46,890-tetrahedron 2 mm partition | No single-atlas requirement; shared interfaces and exclusive material ownership are engineered from retained surfaces. Elastic constants and pelvic friction remain explicit priors |
| Coupled garment contact | Actual elastic panel/tissue dynamics over 240 ms, paired contact forces and frictional work; 4% panel prestrain yields approximately 5.03 mm posterior glans displacement | Positive volume Jacobians; half-step mean glans displacement differs by 5.83 µm. Cloth trajectory differs by up to 1.53 mm; no whole-shorts containment or global collision claim |
| Measured textile friction | Held primary skin/cotton study, extracted chest and dorsal-forearm static/dynamic coefficients and experimental conditions | Region/fabric/condition-specific evidence; unsupported mappings rejected. These measurements do not calibrate genital clothing contact |
| Human tissue constitutive data | Two institutional full texts acquired; 2024 Table1 SI parameters execute as Ogden, HGO and split neo-Hookean laws through the implicit body | Three elderly ex vivo donors; layer-specific inverse fits, no glans assignment or living-body calibration. Stress/energy, objectivity and initial-modulus checks pass |
| Continuing systemic physiology | Persistent native 0.02 s process, 180 observed ports, typed meals/apnea/exercise, actual downstream respiratory/metabolic changes | Exact action ownership and source receipts; no closed whole-body elemental/energy ledger or mechanical work return |
| Predictive planar locomotion | Actual native SCONE feedback controller executes 30 s, 65 steps and 23.7148 m; matched controller removal causes a fall; timestep refinement retained | Planar 9-DOF, 14-muscle source model. Separate 3D 80-muscle contact plant executes short forward dynamics; converged contact-driven 3D walking remains a separate gate |
| Stretch reflex | Source Geyer/Herr tibialis-anterior length feedback with 20 ms delay, 10 ms activation and actual attachment-path mechanics; nerve block changes displacement | Stiff-tendon fiber proxy and selected reduced mechanics; no calibrated global posture controller or IBM descending recruitment |
| Breathing and heartbeat | Native cardiovascular/respiratory state drives the three-mode thorax, chest skin, diaphragm and affine organ/cavity state; optional native chest-compliance intervention | Chest compliance changes the actual native pressure/flow solution; full 3D chest reaction is not fed back into the native circuit |
| Somatic pathways | Timed pressure, muscle-command and selective nerve-block protocols, finite delays, checkpoint rollback and coordinated clocks | Inferred named peripheral support and reduced neural state; no implicit policy converting touch into voluntary contraction |
| IBM reuse | Verified-byte snapshot of 94 source files, real selected materialization, periodic source-runtime parity and exact causal rational/ZOH transfer realizations | Selected tactile/effector dynamics, not a full imported brain connectome or calibrated afferent firing model |
| Fine skin anatomy | 201,528 regional follicle samples, 30,000 physical-radius display shafts transported with computed skin deformation; eight paired forearm vascular units, 256 capillary links and 1,248 total edges | Regional density evidence plus geometry/rheology priors; incomplete hair coverage; vascular boundaries prescribed and native blood storage disconnected |
| Regional contact | 98,304 tetrahedra, 10×10×3 mm synthesized reference volume anchored to a source skin face; 0.1 mm indentation gives 0.0179365 N; computed displacement drives causal IBM receptors | Quasistatic neo-Hookean reference patch, transferred apparent shear prior, inferred bulk modulus/thickness; no global force feedback |
| Independent mechanics | Pinned native FEBio compression plus five rigid-wall contact runs; contact penalty refinement reaches 0.011225% force error against an analytical confined limit | Actual native execution, positive Jacobians and gap/load/release checks; not an equivalent nonuniform indentation benchmark or whole-body contact validation |
| Regional fluid/lymph | Native-law fluid and albumin transport, one-way lymph uptake/return, isolated finite reservoirs; venous pressure, inflow and lymph obstruction experiments | Conserved regional ledger and positive storage; generic regional scaling remains uncalibrated and does not duplicate native-owned blood volume |
| Non-neural electricity | Distinct membrane and apical RC networks, fixed basal reference, baseline/shunt/electrode/membrane perturbations, exact state stepping and independent charge checks | Fixed-reservoir illustrative RC parameters; field observations do not calibrate membrane voltage or establish wound-healing dynamics |
| Kidney raw foundation | Pinned arterial example with 148 nodes, 147 polylines and 19,924 points; 32 paired image/dense-label slices at verified 50 µm spacing | Same donor/right kidney established; graph-to-image and generic-body registration unknown; graph thickness semantics unresolved |
| Temporal analysis | Recorded body/source trajectories and regional Fourier/finite-Laplace views, with units, damping and source receipts | Finite-window signal descriptors, distinct from causal IBM transfer functions; no inferred physiological poles or manufactured confidence intervals |
| Visualization | Fixed 3D viewport; independently scrollable, toggleable and resizable panes; regional camera scales, interventions, trajectories and spectra | Display projection preserves exact recorded times and renderer quantities; source geometry, synthesis uncertainty and numerical precision remain separate |

The kidney image slab contains 490,087 labeled voxels, or 61.260875 mm³ at its documented nominal spacing. Labels cross both z cuts, so this is a slab-limited segmented volume. The box includes background; its labeled fraction is not a kidney tissue vessel fraction. See [the acquisition and metadata record](KIDNEY_MICROSTRUCTURE.md) for raw byte receipts, published/archive slice-count discrepancy and unavailable final-graph DOI. The source card never treats this one-donor example as a calibrated population prior.

## Entry points

```python
from ihm.human import ImplicitHuman
human = ImplicitHuman.open()
description = human.describe()
body = human.materialize('body')
evidence = human.microstructure_evidence()
measured_graph = human.materialize('kidney-arterial-geometry')
transport = human.materialize('body-skin-transport')
transport.step(0.1, lymph_obstruction=0.5)
receptor = human.materialize('ibm-causal', sites_m=[[0., 0., 0.]], response_kind='rapid')
```

The last example deliberately supplies a support coordinate; meaningful body experiments should use the retained anatomical anchor. `body-touch` and `body-skin-electric` produce deferred predictors with `.run()`. The contact default is a substantial reference solve; its builder sets one OpenBLAS thread to avoid oversubscription on this host.

In the UI, **Active materialization** also includes actual garment contact and accepted systemic records in the same canonical anatomical frame. The garment view replaces only the participating tissue meshes and selected shorts faces with recorded positions. Electrical interventions retain separate voltage scales. **Spectrum → Finite Laplace magnitude** displays the active regional record's damping values. **Explore anatomy** includes fine hair and microvascular layers. Those structures have physical dimensions and need regional zoom; their radii are not inflated into whole-body-visible anatomy.

`/api/body/coverage` separates canonical recorded integration from regional experiments and source acquisitions. `/api/body/trajectory?view=display` returns an exact field projection of the full trajectory. `/api/body/experiments/{forearm-touch,skin-transport,skin-electric,spectra,details}` rejects stale source/model receipts, including the spectra's transitive experiment dependencies.

## Reproduction and checks

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_all.py --app
.venv/bin/python scripts/verify_native_respiratory_port.py --resume
.venv/bin/python scripts/verify_reference_contact.py --self-test
.venv/bin/python scripts/audit_body_coverage.py
```

The broad runner writes actual exit codes and output logs to `artifacts/verification/report.json`. The native continuation check resumes an actual saved respiratory state; the reference contact mutation checks reject invalid Jacobians, nonfinite/missing states, residual release, altered energy and false convergence reports. Numerical checks test source parity, conservation, positivity, objectivity, contact/refinement, rollback and evidence integrity. Browser checks exercise computed chest motion, regional states, panel scrolling, controls and spectral/source transitions. None substitutes for empirical whole-body validation.

The latest broad run passes all75 command groups in380.25 seconds, including29 JavaScript tests and28 browser checks (one opt-in native test skipped). Its immutable logs are under `artifacts/verification/run-1788645421340890626/`. It includes native systemic receipt binding, detached-environment selection and acquired tissue constitutive laws. The preceding71- and73-group failures remain retained alongside their fixes. Source-volume CC/CS pulse/refinement checks also pass separately at `artifacts/verification/penile-volume-audit-14mcegde/verification.json`.

The latest persistent native regression is `data/derived/audits/session-vlbgd4t0/actions/verification.json`:180 ports, breathing recovery, batch parity, bounded exercise demand, first-step demand clearing, and actual prestart selection of detached resource files and copied initial state all pass. Historical systemic runs retain frozen numerical inputs and native dependency/resource archives. Their original capture phases remain explicit; new runs archive and select resources before starting the process. `respiratory_v4` reruns the full180-second rest/apnea pair using the depletion-corrected library and this startup path; both causal checks pass and the source-bound projections are now the defaults.

The six-hour meal/hydration pair finishes with nonnegative required stores and measurable glucose, insulin and glycogen contrasts, but core temperature falls to approximately 24 °C. It is excluded from ordinary meal-response visualization. A paired thermal-source audit has identified duplicated clothing/film factors and incorrect parallel-resistance allocation; correcting those improves one-hour core temperature substantially, but continued drift and skin-perfusion heat exchange remain under investigation. These findings are retained model failures, not calibrated whole-body predictions. See [the systemic evidence record](SYSTEMIC_CAUSAL_EXPERIMENTS.md).

Native source corrections also cover calcium unit accounting, renal transport debit caps, dry-gut sodium transfer, bounded exercise power and mixed-unit water depletion. Each correction has isolated source/library identity and native branch evidence; previous libraries and failed trajectories remain available.

Raw inputs remain under `data/raw/`; regional simulations and byte receipts are under `data/derived/canonical/`, kidney products under `data/derived/microstructure/kidney/`, and native FEBio inputs/outputs under `data/derived/mechanics-reference/`. These large/generated trees are largely ignored by Git but retained locally. Builders and source cards are versioned. `scripts/build_workbench.py --plan` gives dependency order; it does not redownload the corpus, rebuild external engines or silently repin IBM.

## Open scientific and integration gates

| Concern | Concrete remaining gate |
| --- | --- |
| Definite volumetric body | Extend the existing non-overlapping pelvic tissue partition to the rest of the body, with material interfaces, contact pairs and joint constraints; audit each inferred partition and thin structure |
| Clothing mechanics | Extend the coupled local panel to the complete garments with edge, self and continuous collision; identify regional constitutive/friction evidence and verify gravity-supported wear without hidden mesh suppression |
| Homeostasis | Correct and verify thermal boundary/perfusion accounting; establish long resting trajectories under explicit bed, ambient and clothing conditions before accepting ordinary long-horizon digestion/exercise predictions |
| Walking | Complete and independently forward-validate a 3D contact-driven controller, map source/canonical anatomy and mass differences, then couple it to soft tissues and physiological work ownership |
| Bed and all-organ mechanics | Benchmark nonuniform regional contact against FEBio, then curved/multiple contacts, friction and gravity; demonstrate supported supine equilibrium and load transfer before whole-body coupling |
| Respiratory feedback | Replace one specified native constitutive contribution with a 3D pressure/work port; close its energy ledger without double-counting recoil |
| Complete blood/lymph geometry | Acquire organ-specific microstructure and boundary observations; fit conditional topology and transport observables; replace regional replicas through explicit native storage ownership contracts |
| Kidney synthesis | Register graph/slab/organ, resolve radius semantics and collapse corrections, quantify segmentation/censoring errors; obtain glomerular, venous and capillary evidence before synthesizing those systems |
| Other organ microstructure | Follow the [organ-specific acquisition matrix](MULTISCALE_BODY_DATA.md) for muscle, liver, lung, heart, gut, lymph and peripheral nerves; preserve species/donor/preparation and measured-versus-synthesized fields |
| Brain and complete innervation | Import future IBM revisions with new hashes and compatibility/parity checks; identify all peripheral sensory/motor/autonomic observation and effector laws, units and target assignments |
| Hair and skin sensing | Add region-specific terminal/scalp hair evidence, follicle mechanics and sensory observations; extend beyond present vellus samples and selected skin supports |
| Integumentary bioelectricity | Identify epithelial layers, ion reservoirs/pumps and intercellular coupling from human data; test perturbation observables before asserting regenerative behavior |
| Minor systems and interactions | Audit endocrine, hepatic/GI, immune/hematologic, reproductive, special senses, autonomic/visceral, connective/adipose/marrow state ownership and perturbations separately; code presence does not establish every interaction |
| Calibration and uncertainty | Separate numerical error, source resolution, registration, biological variation and model discrepancy; require held-out observations for each claimed predictive domain and account for correlated same-donor evidence |
| Scale and visualization | Measure solver, loading and rendering budgets independently; introduce display reduction only with a recorded mapping/error and retain the higher-resolution substrate |

The next IBM vascular and cortical acquisitions remain owned by IBM-1, as requested. IHM's existing import stays pinned; the concurrently changing donor checkout cannot silently alter a saved body. The body-wide microstructure program follows the same evidence discipline, with organ-specific topology and functional validation rather than a universal cortex-derived generator.
