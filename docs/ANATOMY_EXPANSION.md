# Anatomy expansion — 2026-09-05

Accuracy is the primary contract. Original source files, annotations, coordinates,
topology, evaluation settings and hashes remain authoritative. Display geometry
is a derived projection, never a replacement for the source. No geometry is
upsampled to inflate detail or silently repaired.

## Acquired and integrated

| Family | Geometry and identity | Evidence boundary |
|---|---|---|
| BodyParts3D | 2,234 surfaces, 4,209,773 vertices, 6,681,030 triangles | Reference anatomy; all source topology preserved |
| Z-Anatomy | 2,581 evaluated surfaces, 8,111,194 triangles; 495 muscle surfaces, 158 regional lymph-node groups, 256 external regions | Authored atlas partly derived from BodyParts3D; not another independent human |
| Published lymphatic graph | 996 vertices, 1,117 structural edges, 272 lymph-node flags | Model-derived anatomical estimate; no measured radii or flows |
| OpenSim and vascular cases | Original surface topology restored in the display; native wrapped muscle paths and full archived CFD retained | Independent source frames and validity domains |

The application contains 4,982 selectable structures in eight source families.
Node groups, mesh objects, graph vertices and named organs are different counting
units. They are not summed into a claim about anatomical completeness.

## Fidelity and integrity checks

Every BodyParts3D source mesh was compared with its display projection: identical
triangle indices, no decimation and zero coordinate discrepancy under the recorded
rigid transform/unit conversion. The 147 degenerate triangles in the source were
reported and retained. Original anatomy is not thereby certified watertight or
suitable for a physical volume solver.

The ontology crosswalk corrected 230 earlier name-based groupings. Every recorded
relationship path was checked. Four pancreatic primary-group ambiguities remain
explicit; contextual membership is separate from anatomical identity.

Z-Anatomy retains both original datablocks and evaluated surfaces. Original
viewport modifiers and procedural nodes are evaluated without changing authored
resolution. Blender version, world matrices, modifier settings, source hashes and
logs are retained. The source oesophagus dependency-cycle warning remains visible.
Build and integration reject changed raw, cached or projected files before writing
shared display assets. Nine synthetic corruption/integrity regressions pass.

All published lymphatic coordinates, endpoint references, flags and lengths match
the source supplement. Display segments are straight endpoint chords. Their
maximum discrepancy from the rounded published length is 0.000506 mm. Source
From/To ordering does not establish physiological direction. Integration rejects
changed graph or display identities.

## Data records

- `data/raw/anatomy/extended/provenance.json`: pinned official archive, files,
  blend hash, licensing and attribution.
- `data/derived/anatomy/extended/source_index.json`: original/evaluated surfaces,
  transforms and source-authored modifiers; excluded objects remain inventoried.
- `data/derived/anatomy/anatomy_coverage.json`: recorded ontology paths and ambiguities.
- `data/derived/anatomy/fidelity.json`: per-mesh source/display hashes and topology audit.
- `data/raw/lymphatic/provenance.json` and `data/derived/lymphatic/graph.json`:
  publisher supplement, coordinates, topology, source lengths and unknown fields.
- `data/sources/z-anatomy.json` and `data/sources/published-lymphatic-network.json`:
  tracked acquisition descriptors.

The raw corpus is now approximately 6.00 GB. The searchable catalog has 36,656
records, including 1,117 published model-geometry length constraints. Those
lengths are not experimentally calibrated conductances.

## Verification and remaining evidence gaps

The full model/data/API/frontend verification run passed all 35 commands. Its
browser suite passed 10 tests, including full-resolution muscle/skin discovery,
independent layer transparency and both new source families; the unchanged native
intervention test was intentionally not enabled in this anatomy-only run. The nine
new cache-integrity regressions and the additional lymphatic corruption check also
passed separately. `scripts/verify_all.py --app` includes these checks on subsequent
runs. Outputs are retained in `artifacts/verification/`.

The new geometry does not establish dermal layer thickness, glands/follicles,
capillary-resolved whole-body vasculature, measured lymphatic lumen geometry, or
cross-source registration. Unknown radii, directions and flows remain null.
The native lumped lymph model is unchanged; the new graph does not silently add
or double-count its flow paths.

Sources and reproduction: [extended atlas](EXTENDED_ANATOMY.md),
[ontology coverage](ANATOMY_COVERAGE.md), [lymphatic network](LYMPH_NETWORK.md).
