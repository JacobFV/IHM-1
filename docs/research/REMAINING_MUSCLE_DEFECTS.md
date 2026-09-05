# Remaining muscle bulk defects

The seven sources blocked after isolated-pair separation are inspected directly against their original BodyParts3D OBJ files. Every additional transformation must preserve all source faces and coordinates in derived bulk/stratum representations. No canonical replacement, mass assignment, smoothing or hole filling is authorized by this report.

The new unambiguous geometric case is an entire edge-connected component made exclusively of opposite-coordinate triangle pairs whose unoriented support is a manifold disk: edge incidence at most two, one simple boundary loop, Euler characteristic one, and connected path/cycle vertex links. Such a component has no closed 3D enclosure in its unoriented support and its oriented pair contributions cancel. It can be retained as a sheet stratum without selecting an orientation, thickness or anatomical role. A doubled closed tetrahedron, a component partly containing an unpaired bulk face, or a component with a nonmanifold support must remain unresolved.

Degenerate triangles attached to remaining bulk are not deleted merely to obtain a passing audit. If exclusion creates boundary edges, repairing those boundaries changes topology and needs evidence or an explicitly conditioned synthesis proposal. Atlas geometry and numerical degeneracy alone do not establish missing anatomy.

## Retained result and original source evidence

`data/derived/remaining-muscle-bulk-v1/` preserves nine additional complete paired-disk components (72 source faces) as surface strata. Each component consists of four opposite triangle pairs over five coordinate vertices. These are micrometre-scale fans within much larger authored muscle meshes. They were already present in the raw OBJ: every original vertex, after the recorded canonical transform, matches exactly; face order and indices also match exactly. Reports identify the relevant original OBJ line numbers. This supports an authoring/export-artifact hypothesis, but does not establish its cause or the physical role of the sheets.

Six additional muscles now pass the geometric bulk screen: both sets of levatores costarum breves, both descending trapezius parts, right rectus capitis anterior and diaphragm. Combined with the parent decomposition, **424 of 425** canonical muscles have bulk candidates. Every source face and coordinate is retained across the derived bulk/stratum outputs, and none has a mass or physical thickness assigned.

The right aryepiglotticus remains blocked by three exactly collinear triangles: canonical source faces 3118, 3139 and 4079, corresponding to raw OBJ lines 13029, 13050 and 13990. Simply excluding them creates nine boundary edges. That experiment is recorded as an **unapplied proposal**, and the original faces remain in its bulk representation.

A future conservative investigation can test whether each collinear face represents a nonconforming edge subdivision: an existing vertex lying exactly on an opposite edge could permit source-coordinate-preserving retriangulation with explicit barycentric face mappings. This has not been established or applied. If no exact geometric correspondence exists, any reconstructed patch must instead be labeled conditioned synthesis, with explicit local boundary constraints and independent anatomical support. The atlas contains the left aryepiglotticus as a possible anatomical comparison, but contralateral symmetry would remain an assumption, not a measurement of the missing right-side geometry.

The original OBJ license headers are preserved verbatim alongside the canonical license label as source metadata. Legacy header wording alone does not show that the current archive label is wrong; current archive licensing is a separate dated source record.

## Verification

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_remaining_muscle_bulk.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_remaining_muscle_bulk.py --existing data/derived/remaining-muscle-bulk-v1
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/inspect_remaining_muscle_bulk.py --output data/derived/remaining-muscle-bulk-new
```

The full independent check covers all 100,008 source faces across the seven reviewed meshes, exact coordinate bytes, exhaustive disjoint face partitions, original OBJ face-line correspondence, entire paired-disk support, signed cancellation, re-audited remaining topology and retention of unapplied proposals. Synthetic checks reject a doubled closed tetrahedron, a sheet attached to unpaired bulk and nonmanifold support. The verification receipt is retained under `artifacts/verification/remaining-muscle-audit-3bw7x7nh/verification.json`; the run logs are `artifacts/verification/remaining-muscle-bulk-v1*.log`.
