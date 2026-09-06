# Organ-specific microvascular evidence, 2026-09-05

This increment acquires **1,243,271 bytes** of primary articles and public metadata, including a 21,572-byte DERMA-OCTA workbook. It acquires **no new vascular image volumes or vessel graphs**. Exact URLs, acquisition timestamps, SHA-256, sizes and failed requests are in `data/research/organ_microvascular/acquisition.json`; `derive_evidence.py` verifies receipts and extracts tables without network access. `derived_statistics.json` preserves source rows and an inventory of every listed DERMA file, distinguishing held files. UTC acquisition dates are September 6, corresponding to September 5 in the workspace timezone. Brain evidence remains IBM-1's responsibility.

## What the body currently represents

`body_microstructure.anatomical_supports` weights anatomical surfaces by surface area. These weights are engineering allocation proxies, not measured perfused-volume fractions. `body_exchange` conservatively subdivides native BioGears quantities; it does not integrate another vascular reservoir, resolve capillary pressure or reconstruct solute flux.

The 181,082-byte canonical `microvascular.json` contains eight **synthetic** forearm units. `build_body_details.py` calls `microvascular_unit(..., .003, 32, seed=101+j)`, attaches four units to each forearm and prescribes 4500/1500 Pa isolated pressure boundaries. Its own assumptions explicitly say count, length and radius are not measured human morphology. Thus geometry read from a held file is not necessarily measured geometry. The misleading `source_derived_graph_with_inferred_native_territory_allocation` label in `body_exchange` was corrected to `synthetic_graph_with_inferred_native_territory_allocation`; allocation arithmetic was unchanged. Geometric lumen/native-reference fraction remains an inference, including its transfer to extracellular and intracellular views.

Existing kidney bytes are distinct: two Amira files totaling 3,620,045 bytes, a derived NPZ graph, and 64 TIFFs totaling 114,093,568 bytes (32 image/mask pairs). Existing graph/slab manifests remain authoritative; no kidney reacquisition occurred. The audit records selected prior-file hashes, so later source changes can be distinguished from this snapshot.

## Skin: a usable dataset with a cohort discrepancy

[Rotunno et al. 2025](https://doi.org/10.1038/s41597-025-05763-6) and [Zenodo 15088516](https://zenodo.org/records/15088516) provide DERMA-OCTA. The article describes 330 human OCTA volumes, 74 subjects, 76 healthy volumes, 235 chronic venous disease volumes and 19 lesion volumes. TIF volumes are 90×512×483 in Z,X,Y over 0.9×10.0×9.5 mm; 2D projections are PNG. These are flow-sensitive OCTA observations, not complete histological lumens. Do not identify voxel spacing with capillary diameter resolution.

The held XLSX instead contains **330 scan rows, 75 unique literal patient IDs and 51 `HEALTHY` rows from 15 IDs**, located on arm (18) and leg (33). Other labels and their exact counts are preserved. These workbook counts are not reconciled to the paper by guessing, merging IDs or silently relabeling CVI-C0. The mismatch blocks a claim of a resolved healthy population sample, but leaves metadata useful for explicit, provisional scan selection and patient-level grouping.

The API lists **16,899,739,453 bytes** across 14 files. Examples: `Manual segmentations - 3D.zip` 227,118,191 bytes; `Manual segmentations - 2D.zip` 12,026,407; `Norm - 3D.zip` 2,128,148,458; `Contrast - 3D.zip` 4,977,565,511. Only `DERMA-OCTA.xlsx` is held, with publisher MD5 verified. Image archives are listed, not downloaded. Dataset license is **CC-BY-NC-4.0**; article license is **CC-BY-NC-ND-4.0**. These are different grants.

Next conditional input: retain subject ID, disease label, arm/leg, age bin, acquisition depth, preprocessing version and mask status. A bounded future mask/image pair can establish apparent depth-conditioned vessel fraction, skeleton length, components and tortuosity with a tissue denominator and resolution censoring. Preserve projections and 3D masks separately. No current measured skin density, radius prior, or exchange coefficient is extracted from this metadata.

## Muscle: cross-sectional human density and supply domains

[Coupling between skeletal muscle fiber size and capillarization is maintained during healthy aging](https://pmc.ncbi.nlm.nih.gov/articles/PMC5566646/) supplies human vastus lateralis biopsy measurements. The recruited cohort is 47 healthy men/women; Table 3 capillary measurements cover 46. Older participants were unusually healthy and recreationally active. Table 3 values below are **mean ± SD**, not standard errors or uncertainty of a universal constant. Article XML is held (160,919 bytes), under the article's CC-BY license; raw image/graph repository and its size were not established.

| Group | Capillary n | Density, capillaries/mm² | Capillary/fiber ratio | Domain area, µm² | Log-domain SD |
|---|---:|---:|---:|---:|---:|
| Young men | 13 | 331 ± 94 | 1.74 ± 0.57 | 3233 ± 1252 | 0.171 ± 0.022 |
| Young women | 5 | 340 ± 120 | 1.44 ± 0.69 | 3149 ± 855 | 0.166 ± 0.016 |
| Older men | 22 | 286 ± 78 | 1.46 ± 0.41 | 3667 ± 969 | 0.172 ± 0.018 |
| Older women | 6 | 312 ± 139 | 1.08 ± 0.38 | 3434 ± 963 | 0.158 ± 0.019 |

Table 1 group ages are 22.1±2.9, 21.0±2.4, 73.5±3.9 and 74.5±3.7 years, respectively; young-men demographics use n=14. A first synthesis can condition capillary locations on muscle-fiber cross sections and compare all four Table 3 observables jointly. Capillaries/mm² is a numerical section density, not automatically 3D length/mm³ or vessel volume fraction. Longitudinal alignment, tortuosity, radius, branching and permeability remain unsupported by this table. Independent Gaussian draws for correlated density/domain quantities would not recover measured individual configurations.

## Lung: alveolar constraints, not an acquired capillary sheet

[Vasilescu et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7311688/) directly measured 13 single donor lungs: 3 left, 10 right, 6 women, 7 men, age range 25–77 years. Air-inflated lungs were frozen at 10 cmH₂O and systematically sampled by micro-CT. Three specimens showed mild emphysema on radiological review despite no known respiratory disease. Table 2's **present-study** column reports per-single-lung mean±SD: alveolar surface **67±20 m²**, septal wall thickness **12±3 µm**, alveolar number **106±41 million**, lung volume **2883±721 mL**, alveolar/parenchymal volume fraction **71±3%**, duct/parenchymal **13±4%**, tissue/parenchymal **16±4%**. Source table rows are held in derived JSON.

The article HTML (261,150 bytes) is held; Europe PMC full-text XML returned 404. No raw CT archive, capillary mesh, explicit raw-data license or archive size was established. The article is not treated as an unrestricted dataset license. These alveolar statistics constrain an acinar/septal scaffold at the source inflation state. Septal thickness is not blood-air barrier thickness; alveolar surface is not measured capillary endothelial area. Avoid pooling its current measurements with the historical comparison columns or treating each single lung as a pair. Recruitment, capillary-sheet geometry, pressure dependence, diffusion barrier and local ventilation/perfusion need separate evidence before functional gas-exchange synthesis.

## Kidney: preserve the retained evidence's actual scale

[The kidney acquisition audit](KIDNEY_MICROSTRUCTURE.md), [2025 arterial study](https://doi.org/10.1038/s44303-025-00090-2) and [2026 segmentation study](https://doi.org/10.1038/s41467-026-74050-8) establish one LADAF-2021-17 right kidney from a 63-year-old man, ex vivo and processed. The reduced graph has 148 nodes, 147 edges, 19,924 polyline points, mean segment length 1.103236 mm and path/chord tortuosity 1.416294. `thickness` remains uninterpreted; radius versus diameter and reuse license are unresolved. This is not capillary, glomerular or venous completeness.

The paired slab is nominally 50 µm isotropic, 32×1303×912 voxels. Existing measured labels occupy 61.260875 mm³, with both z-faces intersected; the 0.01288795 fraction uses the entire rectangular image box including background, not a kidney-tissue denominator. The CC-BY-4.0 segmentation ZIP is 51,177,854,763 bytes, with only the retained slab acquired. The final 2025 graph remains unavailable at its cited record; an older record is restricted. Do not substitute the reduced example for that graph.

Next inputs are boundary-censored arterial centerline/length constraints in donor-local coordinates. Acquire a tissue denominator before density calibration, establish radius semantics before resistance or Murray analysis, and register the graph/slab before claiming correspondence. Glomerular afferent/efferent serial topology and medullary vasa recta need additional evidence, not extrapolation from a 50 µm overview.

## Liver: human zonation measurements located, numerical curves still missing

[Segovia-Miranda et al. 2019](https://doi.org/10.1038/s41591-019-0660-7) studied human biopsies grouped as normal control, healthy obese, steatosis and early NASH. The total study has 25 patients; extended-data sinusoidal analysis uses NC=5, HO=3, STEA=5, eNASH=3, a different subset from some displayed reconstructions. Approximately 100 µm-thick fixed sections were optically cleared and imaged at **0.3×0.3×0.3 µm voxels**, with stitched central-to-portal vein axes. Fibronectin marks sinusoids. Reconstructed surfaces and centerline graphs were used to measure regional volume fraction, radius, junctions, length per tissue volume, connected fraction, connectivity density and branches crossing regions. Regional summaries are median±MAD; overall figures use box plots.

Held manuscript XML is 207,380 bytes and publisher HTML 481,008 bytes. The manuscript carries academic-research access terms, not a blanket CC-BY grant. No public raw graph/image package with a verified size/license was established. Publisher supplementary XLSX tables are linked, but not acquired; no numerical sinusoid distribution has been digitized or invented. The study's bile-canalicular radius/flow findings must not be assigned to blood sinusoids. Often cited 8.81±2.20 µm sinusoid diameters and 19.2% volume from [large-volume EM](https://pmc.ncbi.nlm.nih.gov/articles/PMC5105151/) are **rat** measurements, not human calibration, and were excluded.

A liver generator should accept portal and arterial inlet identity, central drainage, normalized CV–PV coordinate, tissue mask and the joint spatial observables above. Numerical normal-control curves and a reuse-permitted source graph are acquisition prerequisites for calibration; until then those inputs stay null. This architecture is not a paired binary tree cloned from the forearm.

## Conditional implementation contract

Each evidence row should preserve species, donor/cohort, specimen preparation, body site, disease, measured observable, units, denominator, uncertainty type, spatial scale, imaging resolution, boundary censoring, license, source hash and whether bytes are held. Keep observed geometry, transferred population constraints and generated realizations separate. Record a deterministic seed and unmet constraints for every synthesis. Validate observables through the measurement model (for example cross-sectional sampling or OCTA visibility), not just raw synthetic histograms.

The next useful implementation is an evidence-conditioned scaffold/materialization interface that keeps missing radii, transport laws and unobserved topology explicitly unresolved. It can allocate conservative views of a native owner, but geometric calibration alone does not validate perfusion, oxygen extraction, permeability or lymphatic exchange. This increment neither adds a physiological solver nor claims a whole-organ microvascular reconstruction.

Verification: `python3 data/research/organ_microvascular/derive_evidence.py` verifies seven acquired payloads, checks 330 workbook records and regenerates extracted tables and inventory. `.venv/bin/python scripts/verify_body_exchange.py` runs the existing small ownership checks; the retained-snapshot path now also asserts synthetic provenance for every fine skin partition.

The retained check also passed against `data/derived/audits/native-skin-compression-s725ktdm/zero` without running a native simulation, exercising the new provenance assertion and conservative geometry projection. System `python3` lacks NumPy; the six ownership unit tests passed using the existing `.venv` interpreter.
