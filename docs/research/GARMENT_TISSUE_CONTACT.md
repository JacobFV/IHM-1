# Generated shorts and deformable tissue contact

The generated shorts now supply a real elastic front panel that transfers contact impulses into the source-shaped pelvic tissue volume. Both finite-mass owners advance on one physical clock. The local experiment demonstrates force-driven displacement, elastic preload, friction sensitivity, and load/release response. It does **not** establish whole-garment containment, calibrated human tucking, or a converged cloth wrinkle field.

## Physical experiment and provenance

`ihm/assembly/garment_contact.py` extracts 135 nodes and 225 triangles from the recorded generated shorts, using canonical coordinates `abs(x)<75 mm`, `-145<y<41 mm`, `z>35 mm`. It preserves exact node positions and records original garment node/triangle indices. Forty-three perimeter nodes remain fixed, representing the unresolved far-field garment support. Contracting reference edge lengths by 4% supplies elastic preload; no tissue positions are prescribed. The panel has an engineering spring coefficient of 10 N/m and areal density 0.18 kg/m². Static/kinetic friction coefficients 0.4/0.3 are explicitly uncalibrated genital textile priors.

The participating volume is the source-shaped 4 mm tetrahedral domain documented in [MATERIAL_DOMAIN_DYNAMICS.md](MATERIAL_DOMAIN_DYNAMICS.md): 1,708 nodes, 5,796 positive tetrahedra, three exclusive source material owners, and 0.0762282 kg inherited canonical mass. It combines corpus cavernosum, corpus spongiosum, and glans; it does not include an invented tunica or skin shell. Tissue stiffness remains the inherited generic compressible neo-Hookean prior. Canonical coordinates are metres, x left, y superior, z anterior. Posterior/proximal tissue supports remain fixed. Gravity is zero in this isolated fixture.

The clock advances 240 ms. A 0.02 N superior glans load acts from 120 through 180 ms, then releases. Elastic panel preload remains throughout. Two node/triangle directions are resolved twice each step, transferring equal/opposite impulses between the solids and their fixed supports. A conservative triangle-enclosing-sphere broadphase accelerates candidate search without fixed-count nearest-neighbour truncation. Face-interior contact is implemented; edge/vertex and swept collision remain unresolved.

Per-owner impulse work is the exact kinetic energy change caused by each contact resolution. These works sum to **negative physical dissipation**, not zero when friction or inelastic impact is present. The independent total-energy audit records integration/projection defect separately. Every 5 ms observation stores positions, velocities, cumulative work, and nodal contact impulse divided by the preceding observation interval. The initial force frame is zero. These are interval-mean contact forces, not instantaneous samples.

## Retained results

The complete causal run is `data/derived/garment-tissue-g5me2rlz/`; half-step integration is `data/derived/garment-tissue-refined-oufb46p6/coupled/`. Configuration, source code snapshots, trajectory NPZ, JSONL observations, reports and SHA-256 receipts are retained. The stronger 8% preload pilot at `data/derived/garment-tissue-pilot-dixve2n4/` remains an exploratory artifact and is not the selected display case.

| Case | dt, µs | Final mean glans z shift, mm | Minimum J | Contact dissipation, µJ | Numerical energy defect, µJ |
|---|---:|---:|---:|---:|---:|
| Coupled | 50 | −5.02670 | 0.813094 | 78.8644 | 9.59757 |
| Friction disabled | 50 | −4.44142 | 0.803147 | 47.6069 | 12.3223 |
| Contact disabled | 50 | +1.33462 | 0.985683 | 0 | 0.165537 |
| Coupled, half step | 25 | −5.02881 | 0.813423 | 73.8203 | 4.81370 |

At the refined step, tissue interface work is +58.6475 µJ and panel work is −132.468 µJ. Their sum plus 73.8203 µJ dissipation closes to floating-point tolerance. Maximum pair impulse residual is 3.39×10⁻²¹ N s; maximum pair work balance error is 1.22×10⁻¹⁹ J. Initial elastic energy is 905.134 µJ; final energy is 868.712 µJ after 32.5849 µJ external work. The separately reported 4.81370 µJ numerical defect is approximately 0.53% of initial energy; it is not physical damping.

Halving dt changes the final mean glans displacement by 5.83 µm and approximately halves the numerical energy defect. However, the largest individual tissue-node difference is 155 µm and the largest panel-node difference is **1.53 mm**. Contact dissipation differs by about 6.4%. Two temporal resolutions and one active spatial resolution do not establish full contact/cloth convergence. Further refinement, bending stiffness, and complete contact features are required before interpreting local wrinkles or pressure maxima.

The coupled run has zero final sampled face-interior penetration; disabling contact yields 13.81 mm. This is a causal contact check, not a watertight containment certificate: the refined run reports 34 unresolved edge-contact candidates over its steps and 326 final nearest-edge sample locations whose face projection is outside the triangle. Those sample locations can be safely distant; the number alone is not an intersection count. Full triangle/edge intersection and continuous collision tests remain required.

## Computed display export

`scripts/export_garment_tissue_display.py` exports all original volume boundary faces, with source owner labels, without changing the tetrahedral volume. The selected export is `data/derived/garment-tissue-display-v2/display.json.gz`, with `manifest.json`. It contains 49 actual recorded times, 1,327 boundary nodes, 2,676 boundary faces, the 135-node panel, positions and interval-mean contact forces. Float32 display conversion introduces at most 4.77 nm coordinate error; original float64 trajectories and full volume remain independently retained. This display precision is not anatomical accuracy: synthesized 4 mm boundary geometry has much larger uncertainty.

The payload records material priors, source hashes, runtime source hashes, original garment subset positions and triangle indices, and source garment topology hashes. A viewer must validate correspondence before masking the corresponding original shorts faces, and replace `body-bp3d-FJ3132`, `body-bp3d-FJ3133`, and `body-bp3d-FJ3134` only while this local experiment is shown. It must not superimpose duplicated physical owners or call the remaining garment simulated.

## Reproduction and verification

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_contact_dynamics.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_garment_tissue_contact.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_garment_tissue_contact.py --existing data/derived/garment-tissue-g5me2rlz --refined data/derived/garment-tissue-refined-oufb46p6/coupled
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_garment_tissue_display.py data/derived/garment-tissue-refined-oufb46p6/coupled
```

To create another half-step run, call `run_garment_tissue_experiment(fresh_directory, dt_s=.000025)`. All experiment builders and display exporters reject existing output directories. Rechecking retained results writes a fresh verification receipt. The enriched causal/energy/refinement audit is `data/derived/garment-tissue-audit-jzqjnynf/verification.json`; final display verification is `data/derived/garment-display-check-um2wk6l_/verification.json`.

Tests check contact on/off and friction on/off response, positive volume, actual force transfer, per-owner work balance, total energy closure, saved per-node force balance, temporal refinement, source hashes, source ownership, display positional fidelity, and refusal to overwrite retained artifacts. No scalar endpoint is used to conceal local mesh sensitivity.

## Material evidence still needed

[Khorshidi et al., Acta Biomaterialia 2024, doi:10.1016/j.actbio.2024.06.035](https://pubmed.ncbi.nlm.nih.gov/38945188/) reports human penile tissue experiments and inverse finite-element identification. Its numerical material tables have not been acquired here. They are a target for a tunica/cavernosum/spongiosum law, not justification for treating the current 3 kPa template as measured human mechanics.

[Berardo et al., Applied Sciences 2024, doi:10.3390/app14041357](https://www.mdpi.com/2076-3417/14/4/1357) reports anisotropic tensile and relaxation tests across several lower urinary tract tissues from the same 77-year-old fresh-frozen male donor, including tunica, Buck fascia, bladder, urethra and prostate. Prior prostate resection, preservation, specimen orientation and loading protocol must accompany any inferred material prior. This is a single-donor study, not population calibration. Download attempts against publisher and institution returned HTTP 403 in this environment; no paper bytes or numerical-table acquisition are claimed. High-strain fitted slopes must not be substituted for small-strain constitutive coefficients without identifying the actual law and fit interval.

Current missing physical interfaces include tunica and fascial/skin shells, pressure-dependent porous erectile tissue, true pelvis/thigh contact, textile bending and sew-line mechanics, complete collision features, and globally reconciled mass/density. Until these are implemented and calibrated, this is a reproducible coupled mechanics experiment with visible uncertainty, not a clinically predictive genital or clothing fit model.
