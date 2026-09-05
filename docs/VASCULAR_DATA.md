# Primary vascular assets

This acquisition contains three separate human cases from the official [Vascular Model Repository](https://www.vascularmodel.com/). It provides real image-derived geometry, finite-element meshes, supplied flow boundary inputs and published CFD output. It does **not** provide a registered whole-body subject, a running calibrated solver, or a measured whole-body flow state.

## Local cases

| VMR case | Anatomy / source subject | Mesh points | Tetrahedral cells | Local flow evidence |
|---|---|---:|---:|---|
| `0001_H_AO_SVD` | Aorta; male, 3 years; single ventricle defect | 357,362 | 2,055,236 | Original inlet waveforms, job/preprocessor inputs and outlet resistance BCs |
| `0050_H_CERE_H` | Cerebral / vertebrobasilar; female, 31 years; healthy | 153,082 | 738,037 | Inlet waveforms, 12 resistance outlets, published volumetric pressure and velocity at 201 time-step suffixes |
| `0077_H_PULM_H` | Pulmonary; female, 67 years; healthy | 291,617 | 1,549,224 | Original inlet waveforms, job/preprocessor inputs and outlet resistance BCs |

The cerebral result is a single 1,229,955,608-byte VTU containing temporal arrays, not 201 separate VTU files. The suffixes run from `01010` to `02010` in steps of 5. Arrays include `pressure_01010`, `velocity_01010`, and `average_pressure_01010`, with equivalent fields for the other steps. A numerical export of the first stored state contains the original 153,082 × 3 velocity vectors, pressure values, point coordinates, tetrahedral connectivity and element IDs. All 201 states remain available in the unmodified VTU and original ZIP.

## Reproduce / inspect

```bash
.venv/bin/python scripts/collect_vascular.py
.venv/bin/python scripts/collect_vascular.py --index-only
```

The collector uses Python standard-library networking/archive/XML facilities plus NumPy. The offline pass reads every archive member, checks ZIP CRC integrity, computes SHA-256 hashes, indexes VTK array names and dimensions, and decodes original numerical arrays. It does not generate fields or alter coordinates. Existing downloads are reused; catalog lengths must match. The selected archives plus extracted content remain below the 5 GB collection cap.

- `data/raw/vascular/vmr/*.zip`: four original downloadable archives.
- `data/raw/vascular/vmr/extracted/`: original project files including medical images, models, meshes, boundaries, flow time series, clinical metadata, per-case license, and the full CFD VTU.
- `data/derived/vascular/vmr_index.json`: case identities, metadata, mesh counts, field names, temporal suffixes, decoded first-state ranges, BC definitions and numeric flow time series.
- `data/derived/vascular/vmr_archive_members.json`: exact member names, sizes, CRCs and SHA-256 hashes.
- `data/derived/vascular/*_mesh.npz`: numerical mesh arrays; keys preserve VTK association and name, e.g. `Points/Points` and `Cells/connectivity`.
- `data/derived/vascular/0050_H_CERE_H_first_cfd_frame.npz`: original first-state numerical fields; vector key `PointData/velocity_01010`.
- `data/derived/vascular/vmr_verification.json`: cross-checks of geometry, topology, fields, local file counts and bytes.

## Provenance and evidence limits

The official catalog and license are pinned to [SimVascular/vascularmodel revision `7a1c17dfc9ce4b34eec587fd12f709e8e7baea4e`](https://github.com/SimVascular/vascularmodel/tree/7a1c17dfc9ce4b34eec587fd12f709e8e7baea4e). The catalog supplies patient metadata, anatomy, imaging modality, study DOI and Stanford repository DOI. Its official JavaScript builds the direct archive links used by the collector. File-specific download URLs, sizes and SHA-256 values are in the index. The original license permits research and development with preserved copyright notices and acknowledgements; it is not labeled as an unrestricted public-domain dataset. Every project retains its original license.

The cerebral case links to Bockman et al., [Fluid Mechanics of Mixing in the Vertebrobasilar System: Comparison of Simulation and MRI](https://doi.org/10.1007/s13239-012-0112-8), and the case PDF identifies paper patient 1. The study used phase-contrast MRI to assign inflow boundary conditions and compared simulations with arterial spin-labeling measurements, as described in the [author institution's publication record](https://profiles.wustl.edu/en/publications/fluid-mechanics-of-mixing-in-the-vertebrobasilar-system-compariso/). The downloaded `.flow` files are **processed solver inputs**. Their precise transformation from the original scanner measurements has not been reconstructed. They are not mislabeled as raw MRI measurements. No independent measured pressure time series has been identified in this acquisition.

All CFD pressure and velocity arrays are **simulated**. All three saved `.sjb` jobs report `Simulation failed`; those statuses are retained verbatim. The separately published result archive is not proof that this exact saved job was successfully reproduced. The saved preprocessor uses a parabolic prescribed inlet, rigid no-slip wall, density `1.06`, viscosity `0.04`, period `1.0`, and the solver specifies time step `0.001`. These are original model settings, not newly calibrated parameters. Physical time is not assigned to result suffixes until the job/result correspondence is established.

Units in `.flow` and VTK files are not explicitly declared. Native numbers are preserved. CGS (cm, s, g; flow cm³/s; velocity cm/s; pressure dyn/cm²; resistance dyn·s/cm⁵) is consistent with the supplied density and viscosity, but this is recorded as an inference, not confirmed per-file metadata. Flow sign differs between the supplied `inflow_1d.flow` and `inflow_3d.flow`; neither is silently flipped. Initial-condition arrays and prescribed inlet profiles are labeled separately from simulation results.

## Lymphatic access gap

No paired measured human lymphatic 3D geometry and flow field was acquired. Searches found:

- Tretyakova et al., [Developing Computational Geometry and Network Graph Models of Human Lymphatic System](https://www.mdpi.com/2079-3197/6/1/1), and Savinkov et al., [Graph Theory for Modeling and Analysis of the Human Lymphatic System](https://www.mdpi.com/2227-7390/8/12/2236). These describe anatomical/network models, not a paired measured whole-body flow field. Direct publisher supplementary endpoints returned HTTP 403 in this environment; no agreement or account was used.
- [Human lymphatic vessel electrophysiology and contractility study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4532530/): isolated human vessel measurements exist, but the PMC page presented a browser challenge and its Europe PMC full-text endpoint returned 404. No machine-readable traces or subject-matched 3D geometry were acquired from it.
- [Human skin lymphatic 3D imaging study, GSE282417](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282417): GEO provides transcriptomics metadata; it is not established here as a downloadable paired geometry/flow dataset.

Exact attempted endpoints and their outcomes are recorded in `data/derived/vascular/lymphatic_access_gaps.json`. No lymphatic velocity, pressure, contractility coefficient or geometry was fabricated to fill this gap.

## Verified acquisition totals

The completed raw collection contains **452 local files, 3,167,384,406 bytes**. Four source ZIPs total **1,318,358,528 bytes**; their 442 non-directory members total **1,848,513,708 bytes**. Raw totals include retained archives and extracted copies, so these quantities must not be added together as unique data. All four archive lengths match the official catalog and all indexed members passed ZIP CRC checks. SHA-256 hashes cover source archives, catalog/license files, every indexed member, and numerical NPZ exports.

All three meshes have finite coordinates, valid connectivity indices, and exclusively four-vertex tetrahedra. Cerebral result coordinates and element IDs match the supplied cerebral mesh exactly. Its tetrahedra match after accounting for a different local vertex order; original connectivity in both files is preserved. The first result velocity/pressure arrays are finite with expected shapes `[153082, 3]` and `[153082]`. The first stored velocity magnitude spans `0` to `3.498440222864099` and pressure spans `-2949.832649189771` to `-2662.478614882564`, in unconverted native units. These values verify decoding, not physiological calibration or numerical convergence.
