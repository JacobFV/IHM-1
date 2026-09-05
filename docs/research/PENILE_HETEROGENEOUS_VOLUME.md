# Separate human ex vivo CC/CS volume experiment

## Design and primary review

This experiment assigns the published Khorshidi et al. (2024) point estimates to a new, separate tetrahedral owner containing the existing source-derived corpus cavernosum (CC) and corpus spongiosum (CS) cells. It excludes glans cells and orphan vertices. The parent voxel partition and its explicit overlap priority remain intact; this is a retained subdomain, not a new segmentation or a measured volume. The old garment experiment and all its material receipts remain unchanged.

Table 1, printed page 233, matches the acquired parameter JSON after kPa→Pa and kPa⁻¹→Pa⁻¹ conversion. Equations 1–3 on page 229 match the implemented energy and first Piola stress. Table 3 independently gives initial shear/bulk moduli: CC 22/1111.11 Pa, CS 600/20000 Pa, fascia 50/2500 Pa, and TA isotropic matrix 200000/200000000 Pa. The reported D2=0 cannot be evaluated as the literal reciprocal in equation 2; omitting that higher-order volumetric term is an explicit interpretation whose native Abaqus runtime parity is unresolved.

The new dynamic class dispatches a constitutive law per tetrahedral material index, integrates reference-volume energy, and assembles its exact nodal gradient. Its initial isotropic tangent supplies only a conservative starting step restriction; finite-strain stiffening is not bounded by that initial tangent. Positive determinants and numerical refinement must be checked on each actual trajectory. Density remains the parent's effective mass-allocation prior so a matched generic-law comparison changes only elasticity.

The experiment fixes the retained parent support nodes, applies a smooth, prescribed superior load to a distal source-surface node subset, then releases it. It records actual positions, velocities, interval-mean applied forces, interval support impulses, determinant minimum, work and numerical energy defect. Three time steps test the published-law trajectory and one matched run retains the parent's generic neo-Hookean law. Nodal finite differences test energy/force conjugacy independently of the time integrator.

This is a three-donor elderly, fresh-frozen ex vivo constitutive transfer, not a living whole-penis calibration. No glans, tunica/fiber architecture, fascia/skin shell, resolved urethral lumen, perfusion, active tone, poroelastic pressure, viscoelasticity, garment contact or whole-body material handoff is supplied here. Shared source-volume nodes impose bonded interfaces. Support placement, load and inertia remain explicit engineering assumptions. A 4 mm voxel spacing and 6.93 mm cell diagonal describe discretization, not a measured anatomy-error bound.

## Retained numerical experiment

`data/derived/penile-ex-vivo-volume-v1/` contains four separate immutable runs on the same 1,524-node, 5,172-tetrahedron subdomain: 3,570 CC elements and 1,602 CS elements. Its unchanged effective mass is 68.68 g. A 20 mN peak half-sine pulse lasts 15 ms, followed by 15 ms without applied force. All source inputs, the parent domain, canonical CC/CS surfaces, primary paper, parameter JSON and numerical code are copied and hashed per run. `configuration.law` selects the actual constitutive model; the collected published `constitutive_evidence` is a comparison reference and does not mean the `generic_prior` run used those coefficients.

| Material law | Step | Mean loaded-region superior displacement at 30 ms | Minimum J |
|---|---:|---:|---:|
| Published CC/CS | 50 µs | 0.763788 mm | 0.966229 |
| Published CC/CS | 25 µs | 0.763793 mm | 0.966229 |
| Published CC/CS | 12.5 µs | 0.763794 mm | 0.966229 |
| Parent generic prior | 25 µs | 0.299674 mm | 0.996181 |

Across all saved nodes and times, the maximum coarse-to-half-step difference is 13.42 nm and half-to-quarter-step difference is 3.36 nm. The matched published/generic difference reaches 1.76 mm. These are numerical differences within this prescribed fixture, not an anatomy-accuracy estimate. The published-law numerical energy defect decreases from 16.95 to 8.48 to 4.24 nJ as the step halves; finest applied work is 4.794 µJ. Momentum residuals are below 5.3e−21 N·s. Algebraic energy-audit closure is separate from the reported nonzero integration defect.

Both the two-element heterogeneous test and virtual-work finite differences at the actual final source-volume deformation pass. The latter's directional energy-gradient errors are below 1.92e−10 N. Verification also checks exact geometry and input equality across comparisons, all artifact hashes, positive J, a material-law effect exceeding numerical differences, and rejection of unmapped/glans material assignments. Space refinement, constitutive validation against the original test curves and Abaqus D2 parity remain unestablished.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_penile_volume.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_penile_volume.py --existing data/derived/penile-ex-vivo-volume-v1
# A new experiment must name a fresh directory:
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_penile_volume.py --experiment data/derived/penile-ex-vivo-volume-new
```

Reverification writes a new receipt beneath `artifacts/verification/penile-volume-audit-*`; it does not rewrite retained experiment data. The full run and subsequent extended source-volume checks are logged under `artifacts/verification/penile-constitutive/volume-v1.log` and `volume-audit-v1.log`.

## Explicit materialization

```python
from ihm.human import ImplicitHuman
volume = ImplicitHuman.open().materialize('penile-volume')
evidence = volume.evidence
```

This constructs the actual detached 5,172-element heterogeneous dynamic body,
using the opened body's evidence root. Both the parent partition and acquired
constitutive evidence are hash checked. Its evidence records the source frame,
discretization, effective density basis and source owner identities.
`canonical_handoff_applied` remains false: constructing this predictor does not
add a second copy of its mass to the canonical body or replace the older garment
fixture. The quick verifier exercises this public materialization in addition
to the constitutive and heterogeneous virtual-work checks.
