# IHM-1 implementation and evidence status — 2026-09-05

The workbench opens one generic body at `http://127.0.0.1:8765/`. It now supports computed whole-body replay and separate body-anchored contact and electrical materializations. The implicit substrate binds anatomical, physiological, electrical, population and microscale evidence with explicit source identities. It is not a fully calibrated digital human: global all-organ contact, complete peripheral connectivity and organ-wide microstructure remain incomplete. The research roadmap's broad gates are not marked complete by passing software tests.

## Implemented behavior and retained evidence

| Domain | Executable or acquired result | Scope of the result |
| --- | --- | --- |
| Whole body | 2,408 anatomical representations; regenerated 30-second native-driven record with 301 frames and 2,323 changing reduced transforms | Canonical frame, explicit source fits and assumptions; representation count is not a count of independently simulated organs |
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

In the UI, **Active materialization** selects whole-body replay, forearm contact, or skin electricity. Electrical interventions retain separate voltage scales. **Spectrum → Finite Laplace magnitude** displays the active regional record's damping values. **Explore anatomy** includes fine hair and microvascular layers. Those structures have physical dimensions and need regional zoom; their radii are not inflated into whole-body-visible anatomy.

`/api/body/coverage` separates canonical recorded integration from regional experiments and source acquisitions. `/api/body/trajectory?view=display` returns an exact field projection of the full trajectory. `/api/body/experiments/{forearm-touch,skin-transport,skin-electric,spectra,details}` rejects stale source/model receipts, including the spectra's transitive experiment dependencies.

## Reproduction and checks

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_all.py --app
.venv/bin/python scripts/verify_native_respiratory_port.py --resume
.venv/bin/python scripts/verify_reference_contact.py --self-test
.venv/bin/python scripts/audit_body_coverage.py
```

The broad runner writes actual exit codes and output logs to `artifacts/verification/report.json`. The native continuation check resumes an actual saved respiratory state; the reference contact mutation checks reject invalid Jacobians, nonfinite/missing states, residual release, altered energy and false convergence reports. Numerical checks test source parity, conservation, positivity, objectivity, contact/refinement, rollback and evidence integrity. Browser checks exercise computed chest motion, regional states, panel scrolling, controls and spectral/source transitions. None substitutes for empirical whole-body validation.

The final broad run passed all 63 command groups in 240.25 seconds, including 20 JavaScript unit tests and 20 browser checks. The additional opt-in browser test launched a real two-second native hemorrhage/saline scenario and passed separately through the direct Playwright CLI in 68.33 seconds, with exit code zero. Its native run is `20260905-124726-c7bbeef7`; `artifacts/verification/native-browser/report.json` records the result. The earlier npm wrapper printed success but returned143; that receipt is retained separately. Saved native-state continuation and contact rejection self-tests also passed. The earlier supplemental manifest-loading failure remains in its log: a binary NPZ graph had reached the JSON-only manifest reader. The corrected loader passes the final broad run.

Raw inputs remain under `data/raw/`; regional simulations and byte receipts are under `data/derived/canonical/`, kidney products under `data/derived/microstructure/kidney/`, and native FEBio inputs/outputs under `data/derived/mechanics-reference/`. These large/generated trees are largely ignored by Git but retained locally. Builders and source cards are versioned. `scripts/build_workbench.py --plan` gives dependency order; it does not redownload the corpus, rebuild external engines or silently repin IBM.

## Open scientific and integration gates

| Concern | Concrete remaining gate |
| --- | --- |
| Definite volumetric body | Establish non-overlapping material volumes, interfaces, contact pairs and joint constraints from registered surfaces; audit each inferred partition and thin structure |
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
