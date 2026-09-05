# Cerebral CFD numerical audit

The original archived cerebral velocity/pressure states have good **global flux balance**, but do **not establish periodic convergence or physiologic validity**. The archived job remains marked `Simulation failed`. This audit retains that status, the source pressure gauge, and uncertain case-wide units.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_vascular_flow.py
.venv/bin/python scripts/verify_vascular_flow.py
```

The default audits all 201 original states (solver steps 1010–2010, every 5 steps), corresponding to 1.010–2.010 seconds using the saved solver step of 0.001 s. `--stride N --output DIRECTORY` supports a clearly labeled subset. Generated artifacts are `data/derived/vascular/audit.json`, `cap-flow.json`, `audit.svg`, and `audit.png`; source files remain untouched and their SHA-256 hashes are recorded.

## Geometry and integration

The 153,082 original CFD nodes are identical to the saved simulation mesh. All 738,037 tetrahedra have the same per-cell node sets, although the archive reverses their winding. Boundary orientation is determined geometrically: each triangle's area vector must point away from its tetrahedron's opposite vertex. Unordered connectivity identifies interior faces and cancels adjacent contributions. Overlapping, nonmanifold and degenerate tetrahedra are rejected.

All 14 cap files and `walls_combined.vtp` match **exactly** in original node coordinates and unordered triangle connectivity. Their 99,626 triangles form a complete, disjoint boundary partition. No nearest-triangle replacement or cross-specimen registration is performed. General surface matching retains incomplete-match diagnostics and excludes incomplete surfaces from claimed cap integrals.

For a triangular face, signed flow is its outward area vector dotted with mean nodal velocity. This exactly integrates the piecewise-linear nodal field. Pressure is likewise integrated with triangle-area weighting. Normal-velocity RMS uses the exact triangular integral of the squared linear normal field. Independently, each tetrahedron's linear shape-function gradients produce its constant velocity divergence, multiplied by its geometric volume.

## Results from all 201 states

| Check | Result |
|---|---:|
| Maximum relative discrete divergence-theorem residual | 2.02 × 10⁻¹⁶ |
| Maximum relative global net-flow imbalance | 3.69 × 10⁻⁵ (0.00369%) |
| Mean relative global net-flow imbalance | 5.23 × 10⁻⁶ |
| Maximum wall nodal speed | 0 in every state |
| Maximum volume-weighted local divergence RMS | 1.876 source inverse-time units |
| Prescribed cap flow versus saved BCT: relative RMS difference | 1.17 × 10⁻¹⁶ |
| Prescribed cap flow versus raw waveform: relative RMS difference | 1.71 × 10⁻⁴ |
| Worst cap pressure versus saved resistance × outward flow: relative RMS difference | 4.20 × 10⁻⁴ |

Relative global imbalance is absolute net full-boundary flux divided by the larger of total outward and inward cap flow. Constant density cancels in this ratio. The divergence-theorem result checks geometry and integration of the **same** nodal field; it is not an independent incompressibility test. Local divergence is explicitly nonzero despite strong global cancellation.

The actual saved `bct.dat` vectors are integrated with the same outward normals and interpolated on the original saved time clock with its declared 1 s period. No fitted phase offset or sign correction is applied. Instantaneous `pressure_STEP` fields are used, not cumulative `average_pressure_STEP` fields. Area-mean pressure is compared with the saved `p_ref + R Q_out`; normal viscous traction and finite-element discretization can prevent exact equality.

The cap named `inflow` has **positive outward flow** throughout this archive. The named outlets and `inflow_2` predominantly carry inward flow. Names remain original labels; the audit does not reverse values to match the labels. Moreover, `inflow_2` is a real matched boundary cap but has no assigned condition in the saved `.svpre` file. Thus the saved setup is not a fully reconciled description of every archived boundary.

Snapshots at 1.010 and 2.010 s have identical prescribed inlet phase, yet their velocity nodal RMS difference is **24.86%** of the initial state's RMS and pressure difference is **73.71%**. Mean pressure changes by −2051.26 source pressure units. These unadjusted endpoint comparisons provide direct evidence against treating the saved state as demonstrably periodic. No pressure-gauge shift is fitted away.

## Units and limits

The saved density of 1.06 and viscosity of 0.04 match documented CGS examples, and the image filename contains `-cm`. Those facts support a **CGS inference**, but do not independently confirm the unit consistency of this entire case. SimVascular's own documentation says the solver leaves unit consistency to the analyst. Accordingly, the audit applies no SI conversion and does not interpret negative archived gauge pressures as absolute physiological pressures. [Official SimVascular solver documentation](https://github.com/SimVascular/simvascular.github.io-archive/blob/master/docsFlowSolver.html).

No mesh-refinement study, time-refinement study, full momentum residual, or independent human validation is supplied by this audit. The tests cover analytic affine fields, constant fields, reversed winding, exact scalar-area weighting, corrupted/degenerate meshes, real cap matching, and the actual archived conservation identity. They cannot establish clinical validity or erase the failed-job provenance.
