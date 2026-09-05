# Canonical surface volume-readiness audit

This read-only audit evaluates every canonical reference geometry and its recorded source identities. It merges only numerically identical coordinate triples; it does not move coordinates, remove faces, fill holes, orient surfaces, or construct volumes. Raw and welded topology are recorded separately so duplicated rendering vertices do not masquerade as anatomical holes.

Each unique source mesh is processed once. Tests cover finite coordinates and valid triangle indices, repeated-index and zero-area faces, scale-relative near-degenerate faces, duplicate triangles, undirected edge incidence, pairwise directed-edge winding consistency, vertex-connected components and edge-connected face components. Every component has a signed surface volume integral, both around canonical zero and around its bounding-box center. An open surface's integral depends on origin and must not be interpreted as occupied material volume.

A closed, consistently oriented, nondegenerate surface is only a topological candidate for subsequent volumetric work. Self-intersections, vertex-link manifoldness, nesting/overlap among shells, inter-organ overlaps, physical tissue occupancy and material calibration are not established. In particular, summing closed surfaces can double-count tissue, and skin layers referencing one envelope do not acquire independent solid volumes. Graph/chord and shell-layer references retain explicit non-solid status.

The script writes a new directory containing per-source topology, per-entity provenance and reasons, per-system counts, input snapshots and digest receipts. Source hash failures remain failures in the output rather than silently adopting current bytes. It never rewrites canonical anatomy, geometry, mechanics or previous audit artifacts.

## First complete audit

`data/derived/surface-volume-readiness-v1/` covers all 2,408 canonical entities. Its 2,404 unique triangular source meshes contain 6,836,998 faces; three skin layers reuse one already-audited skin envelope, and one additional reference is a structural graph. All 4,810 canonical/raw source-file identities passed. Unique mesh vertices decrease from 4,289,153 raw rendering vertices to 3,427,602 exact coordinate triples. The pass took 28.7 seconds while processing one mesh at a time.

Exact welding changes 2,063 previously open entity references to closed edge incidence. Overall, 2,381 references have closed incidence, but only 2,089 pass the stricter topological screen. The 319 blocked references include 306 entities with duplicate faces, 16 with zero-area faces, 21 with nonmanifold edges and six with boundary edges; these categories overlap. Three shell-layer references and one graph are explicitly excluded as independent solid candidates. There are 147 exact zero-area faces, one additional near-degenerate face and 3,111 duplicate faces across unique meshes. No faces were removed.

| System | Entities | Topological candidates | Blocked |
|---|---:|---:|---:|
| Arterial | 639 | 627 | 12 |
| Muscular | 425 | 189 | 236 |
| Venous | 395 | 389 | 6 |
| Skeletal | 285 | 272 | 13 |
| Lymphatic | 164 | 163 | 1 |
| Nervous | 146 | 141 | 5 |
| Digestive | 110 | 100 | 10 |
| Respiratory | 106 | 99 | 7 |
| Connective | 48 | 34 | 14 |
| Sensory | 33 | 26 | 7 |
| Cardiac | 23 | 21 | 2 |
| Reproductive | 12 | 10 | 2 |
| Endocrine | 9 | 9 | 0 |
| Integumentary | 7 | 3 | 4 |
| Urinary | 6 | 6 | 0 |

Stomach, spleen, both kidneys and urinary bladder pass the candidate screen. All five lung lobe surfaces retain genuine boundary edges; the right inferior lobe also has nonmanifold edges. Thus welding alone does not make these lung surfaces suitable for closed-volume occupancy. Many muscular surfaces have duplicate triangles even when edge incidence is two; they remain blocked pending an explicit interpretation of those source faces. A total of 578 unique meshes have multiple edge-connected face components, whose separate signed integrals remain available in the source records.

This identifies inputs for subsequent geometry validation and material ownership work. It does not certify self-intersection-free solids, assign tissue volume from the signed integrals, replace the canonical mass ledger or resolve inter-organ overlaps.

Closed arterial and venous atlas envelopes do not distinguish vascular wall from lumen, provide wall thickness, or authorize assigning their enclosed space to tissue mass. Native blood remains the exclusive blood-state owner. Hollow-organ closure likewise does not establish lumen nesting, wall occupancy or the physical interpretation of an interface. All candidate counts here are geometric only.

## Reproduction and receipts

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_surface_volume_readiness.py --self-test
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_surface_volume_readiness.py --output data/derived/surface-volume-readiness-new
```

`sources.jsonl` stores the metrics and component integrals once per source digest. `entities.jsonl` maps every body entity to that record and preserves its source files, raw/welded counts, prior closure flag, current blockers and unresolved tests. `summary.json` provides the system counts. `manifest.json` binds the result files, original source hashes, copied anatomy and exact audit script. Geometry and raw sources are hashed at their existing paths; this is explicitly not a detached archive of every source file.

Synthetic checks cover duplicated-coordinate seams, open boundaries, flipped face winding, repeated indices, collinear triangles, duplicate triangles, separated and vertex-touching components, and reversed closed surfaces with negative signed integral. The full artifact check confirmed all entity/source associations, per-component face counts, manifest digest closure and unset occupied/union volumes. Original source files and prior reports remain unchanged.
