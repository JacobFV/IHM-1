# Native FEBio frictionless planar contact reference

The pinned, unmodified FEBio 4.13 build at commit
`32ae206ff4881dfb54f62296cd1558e58ed9fcc6` passes five native contact runs.
This checks a deformable volume against a moving rigid plane, including open
gap, indentation, unloading, force, penetration, Jacobian, displacement and
native strain energy. It is a bounded contact solver reference, not calibration
or whole-body collision validation.

## Model and independent limit

All values use SI units. The block is 10 × 10 × 3 mm, with compressible
neo-Hookean μ = 2800 Pa and λ = 25200 Pa (FEBio E = 8120 Pa, ν = 0.45).
All nodes are fixed in x and y; bottom nodes are fixed in z. The full top face
contacts an infinite rigid plane with downward normal. There is no prescribed
top-node displacement: the unilateral penalty contact determines it.
The plane begins 0.05 mm above the undeformed top, advances 0.15 mm to produce
nominal 0.1 mm indentation, then returns to its initial gap. Twenty static
increments cover the loading/unloading cycle, with the initial state also saved.
Time labels parameterize quasistatic loading and are not physical dynamics.

The block dimensions, material parameters and indentation depth match the local
regional mechanics example. Its contact footprint and constraints do **not**:
this benchmark uses full-face loading and lateral confinement, whereas regional
circular indentation creates a nonuniform deformation field. No force equality
with that regional model is claimed. Confinement intentionally supplies a known
analytical limit independently of the native solver.

With `F = diag(1,1,s)`, area `A = 1e-4 m²`, height `H = .003 m`, and stretch
`s = 1 + u_top/H`, direct differentiation of the compressible neo-Hookean energy gives:

```text
W(s) = μ/2 (s² − 1) − μ ln(s) + λ/2 ln(s)²
compressive force = −A [μ (s − 1/s) + λ ln(s)/s]
strain energy = A H W(s)
J = s
```

Hard contact at 0.1 mm indentation has `s = 29/30` and force
**0.1073663439085578 N**. Finite penalty relaxes compression by a measurable
penetration; the checker evaluates the constitutive force and energy at the
actual displacement, then separately compares to the hard-contact limit.

## Measured results

| Hex8 elements | Penalty (Pa/m) | Peak force (N) | Hard-contact force error | Peak penetration (m) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1e10 | 0.1072459533864 | 0.112131% | 1.072459534e-7 |
| 8 | 1e10 | 0.1072459533864 | 0.112131% | 1.072459534e-7 |
| 64 | 1e10 | 0.10724595338645 | 0.112131% | 1.072459534e-7 |
| 8 | 1e9 | 0.10617497760316 | 1.109627% | 1.061749776e-6 |
| 8 | 1e11 | 0.1073542921664 | 0.011225% | 1.073542920e-8 |

The three mesh levels have a force spread of about 5e-14 N. Hex8 elements
represent this affine solution exactly, so this is mesh invariance, **not an
observed spatial convergence order**. Independent penalty refinement reduces
force error by approximately tenfold per tenfold penalty increase. At the
finest penalty, minimum native J is 0.966670245143 and peak native strain energy
is 5.28722899389e-6 J. The open-gap and final released states pass the zero-force
and zero-displacement checks.

## Checks and reproduction

After building the reference executable as described in
[the compression benchmark](REFERENCE_MECHANICS_BENCHMARK.md), run:

```bash
python3 scripts/verify_reference_contact.py
python3 scripts/verify_reference_contact.py --self-test
```

Both commands passed on 2026-09-05. There are no additional Python dependencies.
The checker verifies the executable and shared-library hashes against the build
manifest before running. Existing native outputs for each case are removed
before execution. Each run must exit successfully and report normal termination.
At every recorded state, the checker requires the exact time grid, node/element
IDs, row counts and finite values, and checks:

- Full nodal affine displacement within 2e-9 m; positive element-average J and
  consistency with that affine field within 2e-6. The affine field also establishes
  positive local Jacobians for this structured hex mesh.
- Nodal contact force versus analytical force within 2e-6 N, and element stress
  within 0.03 Pa. `Rz` is FEBio's `node.get_load(z)`; sum the loaded top nodes.
  Fixed bottom nodes do not provide the force in this output variable.
- Force versus penalty × reference area × measured penetration within 3e-6 N;
  penetration below 2e-6 m for every tested penalty.
- Native element `sed`, integrated over reference volume, versus analytical
  energy within 1e-10 J. The `.xplt` also retains strain energy density.
- Open gap and release force below 1e-8 N and displacement below 1e-10 m;
  positive loaded force and negative top displacement.
- Mesh force spread below 2e-6 N, recomputed penalty errors strictly decreasing
  by more than fivefold, and final hard-contact relative error below 0.0002.

Mutation checks reject negative J, NaN, missing time state, release residual,
corrupted native strain energy and a falsely reported convergence error. No
spatial order is manufactured from an already exact solution.

Generated evidence is under `data/derived/mechanics-reference/contact/`:
`benchmark.json` records commands, source/build identity, measured time histories,
refinement checks and SHA-256 hashes for each native input/output file. Each
`n*-p*` directory retains `contact.feb`, `contact.log`, `stdout.log`, `nodes.txt`,
`elements.txt` and `contact.xplt`. These are reproducible local generated
artifacts in the repository's ignored data tree.

## Official implementation references

Syntax and native output semantics were inspected in the pinned official source:

- [Rigid-wall implementation](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/FEBioMech/FERigidWallInterface.cpp): plane projection, moving offset, unilateral normal penalty and contact residual.
- [XML 2.5 contact parser](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/FEBioXML/FEBioContactSection.cpp): `ParseRigidWall` and named contact surface.
- [Native output registrations](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/FEBioMech/FEBioMechModule.cpp): `J`, `sz`, `sed`, `Rz` and plot variables.
- [Native nodal force logging](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/FEBioMech/FEBioMechData.cpp): `FENodeForceZ` returns the node's z load.

Further validation still requires nonuniform regional indentation, friction,
curved or multiple interacting surfaces, anatomical geometry, material fitting
and coupled sensor response. This benchmark establishes none of those claims.
