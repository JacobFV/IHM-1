# Locally acquired anatomy data

This collection contains real published mesh geometry and musculoskeletal model parameters. It does not constitute a subject-specific digital human. Source coordinates are preserved; no registration between the atlas, OpenSim models, or IHM runtime has been performed.

## Contents

| Collection | Local data | Counts | Exact raw bytes |
|---|---|---|---:|
| OpenSim official models repository | `data/raw/anatomy/opensim-models/` | 1,038 files including archive; 343 mesh files; 34 models under `Models/` | 554,465,054 |
| BodyParts3D 4.0 official DBCLS atlas | `data/raw/anatomy/bodyparts3d/` | 2,234 OBJ meshes, 4,209,773 vertices, 6,681,030 faces; six ontology tables with 56,036 rows; archive and license evidence | 624,843,860 |
| **Total** | Includes archives and extracted source files | **3,282 files** | **1,179,308,914** |

The OpenSim collection includes Rajagopal 2016 and Lai/Uhlrich 2023, Hamner, gait2392/gait2354, Arm26, wrist/hand, leg and example models. Across the 34 overlapping model variants the parser records 393 bodies, 393 joints, 1,254 muscle definitions and 4,480 muscle path points. These totals are **not counts of distinct human structures**. Some downloaded examples are mechanical demonstrations. Tutorial and pipeline variants remain raw and are not included in the 34-model derived index.

Rajagopal2016 alone contains 22 bodies, 22 joints, 80 muscles and 288 path points. The wrist model contains 28 bodies, 28 joints, 25 muscles and 210 path points. BodyParts3D supplies a broad atlas including bones, organs and vessels, but not a complete physiological model of those systems. Its official distribution is 99% polygon-reduced, not full-resolution anatomy.

## Sources and rights evidence

OpenSim source: [official opensim-org/opensim-models repository](https://github.com/opensim-org/opensim-models), pinned to revision `d9b05d470b1a481c222372c85b75772faf8f7792`. The [official Rajagopal notes](https://github.com/opensim-org/opensim-models/blob/d9b05d470b1a481c222372c85b75772faf8f7792/Models/Rajagopal/README.txt) describe upstream muscle path and force-curve modifications and cite the associated papers. Original credits, publications and README files are retained. This repository has no root LICENSE file; this collection does **not** infer that the separate OpenSim software license licenses every model/geometry asset. Asset redistribution terms remain unverified.

BodyParts3D source: [official DBCLS archive downloads](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html), specifically `isa_BP3D_4.0_obj_99.zip`, last modified 2013-05-22 according to the server. The [official license page](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html), updated 2025-02-27, specifies **CC Attribution 4.0 International**. Local `lic.html` and `README_e.html` retain this evidence; older mirrors may still quote the previous license. Required attribution: **BodyParts3D, © The Database Center for Life Science licensed under CC Attribution 4.0 International**.

Downloaded on 2026-09-05 UTC (2026-09-04 Pacific). Archive integrity:

| Archive | Bytes | SHA-256 |
|---|---:|---|
| OpenSim `source.zip` | 119,531,629 | `377c0f2e972bd17e2f81d29074ea08a10af81205ba43df191bf5b86287274578` |
| BodyParts3D `isa_BP3D_4.0_obj_99.zip` | 142,903,898 | `40665852c49f218326590e204db91064a1ecfc3c6f8cbd7bbbcaac62c7cd409e` |

## Derived records and coordinate semantics

Run `.venv/bin/python scripts/collect_anatomy.py` to extract retained archives and regenerate the indexes. Use `--skip-extract` to index existing extracted files, or `--download` to retrieve absent source archives/metadata directly from the official sites. Downloads use temporary files and promote them only after completion. The manifest hashes pin the data observed during this acquisition; DBCLS's `LATEST` URL can change on a future download.

All derived records are under `data/derived/anatomy/`:

- `inventory.json`: source URL, path, exact byte count, SHA-256 and license status for every raw file; repository revision and archive version.
- `opensim_index.json`: per-model paths and counts, plus totals across variants.
- `opensim__<folder>__<model>.json`: bodies, inertial properties, joints, offset frames, geometry references, muscle coefficients, origin/via/insertion path definitions, conditional/moving-point functions and wrapping objects. Recursive source definitions preserve XML child order and identifiers. Defaults are retained separately and are not silently expanded. Legacy C++ `::` XML tag names are escaped only while parsing and restored in the derived representation.
- `bodyparts3d_index.json`: all mesh paths, hashes, bytes, vertex/face counts and bounding boxes in original source coordinates, plus ontology mappings and all six original tables. Every OBJ has an associated concept mapping; compound-organ membership is preserved rather than mistaken for an independent mesh.

OpenSim parameters are source model parameters, not new measurements. Coefficients record their source strings and units where known; unknown unit semantics are explicitly marked. Length/force declarations, credits and publications are retained at model level. A path's first/last point denotes model origin/insertion order, not a measured attachment footprint. Moving points must be evaluated using their original coordinate functions; wrap objects influence the actual path. The parser does not approximate these as fixed world-space points.

Atlas coordinates have not been assigned assumed units or transformed to OpenSim. Mesh bounding boxes preserve source values. No conversion, subject fitting, physical attachment matching, segmentation validation or biomechanics solving has been performed by this import.

## Verification and remaining gaps

All 34 primary model files parse, including the legacy XML tag syntax. Every one of the 4,480 muscle path points retains a parent body/frame reference. All 2,234 atlas meshes contain vertices and faces and map to concept records. Archive CRC and inventory SHA-256 checks pass. Geometry references resolve locally except two `Unassigned` placeholders and the upstream `Leg39/leg39.osim` reference `fibula_R.vtp`, which is not present in the pinned repository; it is explicitly unresolved in the output.

This pass provides real geometry and attachment/model data, but it does not add a calibrated detailed spine biomechanics model, subject-specific muscle architecture, tissue material measurements, comprehensive neural innervation, or a runnable coupled whole-body solver. BodyParts3D spine geometry is atlas anatomy; it is not a spine dynamics calibration. Upstream model validation/calibration described in source publications is distinct from calibration of the IHM target subject.
