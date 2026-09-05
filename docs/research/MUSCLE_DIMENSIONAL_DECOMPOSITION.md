# Retained mixed-dimensional muscle geometry

This new derived representation separates only edge-connected components of exactly two nondegenerate triangles whose coordinate bytes are identical and cyclic windings opposite. Both faces remain in a zero-thickness surface sidecar with their original source-face indices and original coordinate array. Their physical role, thickness, material and mass remain unknown. Opposite-winding pairs inside larger edge-connected components stay in the remaining geometry.

Every original face belongs to exactly one representation: the remaining bulk-candidate surface or the retained paired-surface stratum. No coordinates, normals or face winding are changed. No holes are filled and no source faces are discarded. A muscle containing only paired strata is reported as having no volumetric bulk. The remaining geometry is re-audited for closure, degeneracy, winding and components; even a successful bulk candidate does not establish self-intersection freedom, shell nesting, occupied tissue volume, activation law or exclusive mass ownership.

The geometric basis is cancellation, not inferred anatomy: reversing a triangle's cyclic orientation reverses its signed solid angle at points away from the triangle and reverses its signed tetrahedral volume contribution about any origin. The two retained triangles therefore have zero combined winding and signed volume. Tests must verify that identity numerically, off the surfaces, and retain its floating-point residuals. The pair's geometric surface remains present and may eventually require a physical interpretation; cancellation does not authorize deleting it or assigning it a thickness.

For a triangle `(a,b,c)`, orientation reversal gives `Ω(a,b,c;q) = −Ω(a,c,b;q)` off the surface. Its signed volume contribution about origin `o` is `(a−o) · ((b−a) × (c−a)) / 6`; reversal negates the cross product. Cyclic rotation preserves both quantities. The implementation checks exact coordinate-byte/cyclic identity and samples each retained pair on both sides of its plane, at normal offsets of ±0.1, ±1 and +3 maximum edge lengths. It checks volume cancellation around canonical zero, the pair centroid and another displaced origin.

## Full retained decomposition

The successful run is `data/derived/muscle-dimensional-decomposition-v2/`. All 425 canonical muscles are represented, with every source face retained exactly once across the two outputs.

| Result | Count |
|---|---:|
| Original geometric candidates | 189 |
| Remaining bulk geometric candidates | 418 |
| New bulk candidates with unresolved sheet strata | 229 |
| Isolated opposite pairs retained as strata | 1,866 |
| Faces retained in surface strata | 3,732 |
| Faces retained in bulk-candidate geometry | 2,400,872 |
| Still-blocked bulk sources | 7 |
| Pair-only sources with no bulk | 0 |
| Faces discarded | 0 |

The seven blocked sources are the left/right sets of levatores costarum breves, left/right descending trapezius, right rectus capitis anterior, right aryepiglotticus and diaphragm. Their larger connected components still contain duplicate/nonmanifold topology; the right aryepiglotticus also retains exact-zero-area faces. Opposite faces connected to these larger components have not been separated or deleted.

The maximum sampled pair winding residual is 1.77e−16. The maximum stable signed-volume cancellation residual is 1.24e−24 m³, or 7.17e−15 normalized by the larger contribution/edge scale. These numerical residuals support geometric cancellation and do not quantify biological accuracy. Remaining bulk components still require self-intersection, nesting, overlap and material-ownership checks before volumetric mechanics. Every original coordinate remains in both coordinate arrays; unused coordinates do not imply material occupancy and consumers must use the retained triangle indices.

## Numerical failure retained

The initial `v1` attempt remains intact with `failure.json`. Direct triple products of absolute coordinates lost relative precision for tiny foot triangles far from canonical zero: pair volume residual reached 3.50e−19 m³ and its normalized value exceeded the chosen check. Exact opposite orientation and solid-angle cancellation passed. A regression reproduces this conditioning issue. The successful run uses the algebraically identical local-edge determinant above, preserving the naive residual as a separate diagnostic. No pair-classification criterion or tolerance was loosened.

## Verification and reproduction

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_muscle_dimensions.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_muscle_dimensions.py --existing data/derived/muscle-dimensional-decomposition-v2
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/decompose_muscle_dimensions.py --output data/derived/muscle-dimensional-decomposition-new
```

Synthetic cases check isolated opposite pairs, retention of a pair attached to a larger edge-connected component, same-winding rejection, a pair-only source with no bulk and translated small-triangle numerical stability. The independent artifact verifier checks source/result digests, exact coordinate bytes, a disjoint exhaustive partition of original faces, retained triangle indices, component size two, opposite-cyclic coordinate identity and the cancellation tests again. Reverification writes a new receipt under `artifacts/verification/muscle-dimensions-audit-*`, preserving all original results.

`sources.jsonl` records original and remaining-bulk topology, both geometry representations and cancellation results. `entities.jsonl` maps every muscle to that record and explicitly leaves mass unassigned. Per-source `mappings` retain both source-face lists and each classified pair. `manifest.json` binds all derived bytes, copied anatomy/code and raw/canonical source identities. No canonical geometry, metadata, mass ledger or runtime mechanics has been changed.
