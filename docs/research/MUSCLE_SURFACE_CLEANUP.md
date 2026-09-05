# Derived muscle duplicate-face cleanup

The permitted transformation is deliberately narrow: preserve every original coordinate, retain the original vertex array, and remove repeated triangle occurrences only when their coordinate bytes and cyclic orientation agree. If an unoriented triangle group contains both windings, retain the whole group as unresolved. Exact zero-area triangles can be excluded separately with their source-face IDs recorded. No holes are filled, near-degenerate faces removed, vertices moved, normals repaired or canonical files rewritten.

Each output must include a source-face→retained-face mapping, retained source-face IDs, excluded zero-area IDs, removed duplicate IDs and ambiguous opposite-winding groups. A duplicate maps to the first retained matching occurrence; zero-area faces map to −1. Derived topology is re-audited with the same checks as the full canonical audit. Triangle-set equivalence applies to nondegenerate duplicate-only geometry; surface integrals/areas that count multiplicity can change when duplicates are removed and are not evidence of a new physical volume.

Resulting candidates remain geometric. Self-intersections, vertex-link manifoldness, source shell nesting, inter-muscle overlap, muscle-fiber architecture, attachment sites, tissue material/activation laws and exclusive mass ownership remain unresolved by this cleanup. No generic atlas envelope is promoted to calibrated muscle mechanics.

## Complete canonical muscle result

The retained run `data/derived/muscle-surface-cleanup-v1/` processes all 425 muscle sources. There are **no same-cyclic duplicate faces eligible for removal**. Instead, it finds 1,902 opposite-winding groups, each containing exactly two faces. All are retained. Three exact-zero-area faces are excluded from the right aryepiglotticus, which still has boundary edges, nonmanifold edges and duplicate faces.

The candidate count consequently stays **189 of 425**, with 236 still blocked and **zero new candidates**. This negative result identifies why a duplicate-removal recipe cannot resolve the muscular volume problem. The outputs contain 1,866 edge-connected components of exactly two faces, suggesting that interpreting mixed surface/solid strata merits separate work. Some opposite pairs belong to larger connected structures; no blanket removal or physical reinterpretation is justified by this count.

The run took 120.2 seconds. Every original coordinate array survives a float64-byte comparison after derived JSON serialization and decompression. Original source indices determine every retained triangle. Per-source mappings retain all ambiguous faces, and all geometry/mapping/result files pass their manifest hashes. The original canonical geometry and metadata remain unchanged, including the original source face counts; derived geometry records its own retained face count and source identity. Vertex normals remain those of the source.

## Verification and reproduction

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/clean_muscle_surfaces.py --self-test
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/clean_muscle_surfaces.py --output data/derived/muscle-surface-cleanup-new
```

The synthetic tests prove equivalence of the oriented coordinate-triangle set when removing same-cyclic duplicates, retention of the entire ambiguous group when opposite winding appears, explicit zero-area exclusions, preservation of representably different coordinates and retention of a tiny nonzero-area triangle whose area norm would underflow. The implementation tests cross-product components directly for zero to avoid that norm-underflow error.

`sources.jsonl` records complete before/after topology and integral diagnostics; `entities.jsonl` joins muscle identities to their derived geometry and mappings. `manifest.json` binds all result bytes and the recorded raw/canonical source identities. Source anatomy and audit/cleanup code are copied into `inputs`. The full-run log and independent 425-source hash/coordinate/mapping checks are retained at `artifacts/verification/muscle-surface-cleanup-v1.log` and `muscle-surface-cleanup-v1-verification.json`.
