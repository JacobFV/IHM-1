# Source-shaped material domains and dynamic contact

Implemented 2026-09-05. The new material layer supports an engineered unification
of source anatomy and explicit synthesized volume partitions. Each existing
canonical entity retains a source and material identity. The current resolved
pelvic volume is a numerical research fixture; it does not establish whole-body
soft contact or biologically calibrated genital mechanics.

## Source volume and heterogeneous ownership

`scripts/build_body_material_domains.py` reads canonical anatomy and mechanical
material priors. It records all **2,408 entities**, including bones, muscles,
tendons, organs, nerves, vessels and skin, in `registry.json`. An unresolved
volumetric domain remains visible as unresolved; a surface record is not silently
promoted to a volume or a collision interface.

The first resolved partition uses three actual canonical BodyParts3D surfaces:
corpus cavernosum FJ3132, corpus spongiosum FJ3133, and glans FJ3134. Their earlier
metadata reported open surfaces because of duplicated mesh vertices. Welding
**exactly equal coordinates**, without moving any point or filling any hole,
produces closed, consistently oriented surfaces: respectively 648/1,272,
290/576 and 693/1,374 vertices/faces. The original bytes remain unchanged.

Cell-center generalized winding numbers condition a synthesized voxel volume.
Every occupied cube is split into six positively oriented tetrahedra; adjacent
material labels share nodes. Intersecting source cells have one owner, with the
explicit priority cavernosum → spongiosum → glans (later wins). This is a
documented partition assumption, not evidence that the overlapping source
surfaces are perfectly registered tissue interfaces.

| Resolution | Vertices | Tetrahedra | Union volume | Overlapping source cells |
|---|---:|---:|---:|---:|
| 4 mm | 1,708 | 5,796 | 61.824 mL | 16 |
| 2 mm | 10,726 | 46,890 | 62.520 mL | 170 |

The volume changes by 1.126% between these grids. This is a discretization
comparison, not anatomical validation. Boundary positions remain limited by the
voxel geometry; the cell diagonal is 6.93 mm or 3.46 mm respectively. Increasing
triangle count does not create additional source observations.

The three source owners allocate **0.0762281941 kg** from the existing canonical
mass ledger, once. `MaterialOwnership` rejects duplicate claims, and the dynamic
step rejects repeated physical owner IDs even when different Python objects are
supplied. These are detached replacement domains. The manifest explicitly says
the old affine owners must be deactivated before activation in the global body;
that runtime handoff has not been performed here.

Effective inertial densities are existing allocated mass divided by voxel
material volume. At 4 mm they are about 1,239, 1,259 and 1,134 kg/m³. These values
**are not independent tissue density measurements** and exceed the inherited
900–1,100 kg/m³ density prior range. This exposes the existing whole-body mass
allocation/geometry mismatch rather than concealing it. Resolve that allocation
before interpreting domain inertia as measured human inertia.

## Equations and solver ownership

`DynamicTetrahedra` reuses `DeformableRegion`'s compressible neo-Hookean energy,

`W = μ/2 (I₁ - 3) - μ log(J) + λ/2 log(J)²`.

The gradient supplies actual nodal elastic forces. Per-tetrahedron μ, λ and
density arrays allow distinct constitutive parameter fields. The current three
pelvic labels deliberately retain the same weak generic soft-organ elastic
template (E=3 kPa, ν=0.45), with separate mass allocations. Different labels do
not imply empirically measured differences that the available data do not yet
support. Their `calibration_status` explicitly says the penile material tables
have not been acquired.

`step_coupled` advances all elastic owners on one clock using symplectic Euler.
It includes inertia, body forces, fixed support reactions and optional finite
mass contact. It works directly with the existing `Cloth` object's position,
velocity, mass, triangle and `elastic_forces` interface. All proposed states
are checked before any owner is committed. Nonfinite inputs, invalid clocks,
duplicate owners, excessive steps and inverted tetrahedra fail without partial
state advancement. The initial elastic-wave step bound is conservative for the
tested small strains; it is not a universal nonlinear stability proof.

The contact kernel projects a node to its nearest interior triangle feature,
using barycentric target masses. Normal impact has zero restitution. Tangential
impulses use static/kinetic Coulomb bounds and finite relative effective mass.
Both owners move and receive equal/opposite impulses. Position projection also
moves both masses; its energy effect is recorded as a numerical defect, not
invented physical dissipation. Normal/friction impact loss is computed from the
actual relative velocity and inverse mass. Every receipt contains paired linear
and angular impulse residuals, contact counts, unresolved edge contacts and
preprojection penetration.

The solver audits `ΔP - external impulse - support impulse` and
`Δ(K + elastic + gravitational energy) - external work + contact dissipation`.
It does not force the latter residual to zero. Mesh/time refinement must reduce
the numerical error. Nearest edge/vertex contacts, continuous swept collision,
self-collision, arbitrary multiple-contact convergence and whole-body collision
coverage remain unresolved. These limitations preclude a global containment
claim from this kernel alone.

## Executed evidence

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_body_material_domains.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_body_material_domains.py --spacing-m .002
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_contact_dynamics.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_material_domains.py --canonical data/derived/material-domains/pelvis-0.004m
```

Builders require fresh output directories; use `--output-dir` to retain new runs
alongside previous evidence. Generated meshes, original source hashes, material
and literature cards, ownership contracts and registry are retained under
`data/derived/material-domains/`.

Tests were observed failing before the missing material API, dynamic API and
duplicate-owner rejection were implemented. They verify positive and unique
tetrahedra, explicitly overlapping source partition, exclusive mass ownership,
rigid-motion objectivity, analytical stress, physical parameter sensitivity,
finite mass contact impulse/dissipation, common clocks and failed-step rollback.

The regional constitutive stress test independently reproduces the pinned native
FEBio compression result −105.555555556 Pa. A fresh native rerun is retained at
`data/derived/material-domain-febio-prbibafj/compression/benchmark.json`; its
maximum stress error is 4.44487e-10 Pa. This verifies the constitutive convention,
not dynamic contact, anatomical material calibration or whole-body equivalence.

The inertial bent-beam refinement at dt=100/50/25 µs yields accumulated energy
defects 9.23817e-9 / 4.59525e-9 / 2.29165e-9 J from initial strain energy
3.28125e-6 J. The final position difference between the last two time grids is
1.03605 µm (Euclidean all-node diagnostic). Doubling stiffness in only the distal
material region or doubling density changes the computed trajectory. Maximum
momentum residual is below 5e-21 N s. Receipt:
`data/derived/contact-dynamics-_ivr7ien/verification.json`.

The actual canonical pelvic volume was then advanced for 50 ms with the declared
posterior support fixture and a total 0.02 N superior load distributed over glans
nodes. Mean glans displacement is **1.39919 mm superior**, maximum nodal
displacement 1.61540 mm, minimum J=0.9867135 and maximum momentum residual
8.11e-21 N s. Halving dt from 50 to 25 µs changes mean glans displacement by
0.499751 µm and halves the energy defect magnitude from 1.34691e-8 to
6.73902e-9 J. The unloaded control remains stationary within 1.74e-18 m.
Full computed trajectories and receipts are retained at
`data/derived/canonical-material-dynamics-0e96yt93/`. This establishes actual
source-shaped elastic deformation under an assumed force; it is not a calibrated
tucking response or a garment-driven experiment.

An independent clothing-agent probe also coupled a finite-mass cloth patch to
the dynamic solid for 100 one-microsecond steps: 93 contact resolutions, no
unresolved edges, maximum momentum residual 4.08e-24 N s, positive J≥0.9999711
and numerical energy defect 6.96e-15 J. Its source-hashed receipt is
`artifacts/verification/clothing/two-way-probe.json`. Full shirt/shorts and body
contact remain separate work.

## Evidence needed for anatomically realistic mechanics

The next parameter source is the **human fresh-frozen** study
[Khorshidi et al., 2024](https://pubmed.ncbi.nlm.nih.gov/38945188/), DOI
10.1016/j.actbio.2024.06.035. Its inverse-FE approach combines whole-segment and
individual-layer testing. Its numerical tables are not yet held. The preceding
[Bose et al., 2024](https://pubmed.ncbi.nlm.nih.gov/38494081/) study uses **horse**
tissue; its fitted parameters cannot silently become human coefficients.

[Fereidoonnezhad et al., 2023](https://research.tudelft.nl/en/publications/development-of-in-silico-models-to-guide-the-experimental-charact/)
identifies tunica albuginea as a major contributor to whole-penis load response
and distinguishes indentation from plate-compression information. This volume
currently lacks a resolved tunica, fascial/skin shells, urethral lumen and
vascular/poroelastic pressure state. Those omissions matter directly to large
deformation and clothing containment. Apparent modulus values from
[in-vivo shear-wave measurements](https://pubmed.ncbi.nlm.nih.gov/33945173/) are
measurement-context evidence, not interchangeable large-strain material laws.
The literature cards retain these distinctions rather than pretending acquired
calibration tables exist.
