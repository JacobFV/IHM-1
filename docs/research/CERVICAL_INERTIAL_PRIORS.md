# Geometry-derived cervical inertia and exclusive torso partition

The new opt-in materialization is `data/research/cervical_inertia/v2/manifest.json`. It supplies physically admissible cervical/head/jaw inertial priors and a conservative debit from the **actual mass-scaled native torso**. It does not modify the native model or authorize a cervical run. Joint-frame registration, muscle-frame binding and reassignment of skin/contact owners remain separate integration work.

The first trial was `v1`; `v2` preserves the same physical result and additionally archives every execution input. The retained `v2` artifact is the reviewable result. No original donor tensor or canonical geometry was edited, no inertia eigenvalues were clipped, and no OpenSim `Body.cpp` fallback was used.

## Inputs and geometric scope

The frozen reference is `data/derived/audits/cutaneous-factory-ylro2d66/body/mechanics`: its initial registration, canonical mechanics, identity, assembled native model and input model are all captured. The assembled model contains 22 bodies and 92 muscles. Its torso mass is **27.654676965260336 kg**, not the unscaled source XML's 30.3832390808 kg. The transformation used is the retained canonical-to-torso embedding, checked against the inverse initial torso/world transforms; no moving or settled pose supplies the COM locations.

The donor is the original MASI v2 acquired in `CERVICAL_MODEL_SOURCE_AUDIT.md`. Its seven cervical masses, 3.8 kg head composite and 0.2 kg jaw are retained without an additional mass scale. Their sum is 7.418919568221101 kg. Original tensor values remain in each proposed body's record for comparison.

Canonical C1–C7 geometry is FJ3176, FJ3177, FJ3161, FJ3164, FJ3167, FJ3170 and FJ3172. Head geometry combines fifteen explicitly listed craniofacial bone surfaces; the separate mandible is FJ3289. These 23 source geometry files contain 83,298 vertices in total, occupy 3,838,096 compressed bytes and are copied unchanged. The code rejects display-decimated meshes and bounds the size of each input. Original coordinates, faces, normals and coordinate-transform metadata remain intact in the archived files.

The geometry comes from BodyParts3D, © The Database Center for Life Science, under [CC Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). The held provider README is retained as `data/research/cervical_inertia/BodyParts3D_README_e.html`; the corresponding [provider description](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/README_e.html) and original acquisition code identify the source. This work adds convex-envelope moments and explicit mass-distribution priors; it does not alter or relicense the source geometry. MASI attribution and its Figshare MIT metadata remain in `data/research/cervical` and the input archive.

## Why the source surface integral is not used

All selected source groups have unresolved open indexed boundaries. The C1–C7 meshes have 632–938 boundary edges each; an edge-incidence test is not proof of no self-intersections even when it passes. Thus the new code does **not** use their signed surface integrals as material volumes.

Instead, it constructs a convex hull of the actual source vertices and integrates that explicitly approximated solid. Every hull facet is oriented outwards; its tetrahedron with an interior reference point has positive volume; the hull edge manifold closes; the independently accumulated volume agrees with QHull. Convexity defines the closed-volume approximation. This fills canals, concavities and some gaps, so it is not a recovered bone segmentation.

Tetrahedral integration computes the volume centroid and full second moment analytically. For tetrahedron vertices `v0…v3`, the uniform first moment is `sum(v)/4` and second moment is `(outer(sum(v),sum(v)) + sum(outer(v,v)))/20`, multiplied by its volume. Coordinates are translated before integration to reduce cancellation. The final covariance yields the full COM inertia, including off-diagonal terms.

The artifact reports source-vertex and source-face-centroid depths inside the hull, hull-facet-center distances to the nearest source vertex, hull/AABB volume ratio and AABB-vs-hull covariance sensitivity. These are explicit finite-sample geometric discrepancies, **not** a global Hausdorff bound. Cervical source-face centroid depths reach 8.76–13.09 mm; the head union reaches 59.94 mm because the convex composite fills the intracranial and other concave spaces.

## Mass-distribution priors

For each cervical level, a uniform bone-proxy core occupies its convex hull with generic density 1900 kg/m³. Its mass is `rho_b * V`. The remaining donor segment mass is assigned to a concentric homothetic soft shell at the generic prior 1000 kg/m³. A core exceeding the donor mass is rejected, not clipped. The shell excludes that core within the segment.

If the hull covariance is `C`, its volume is `V`, and the shell expansion is `s`, then:

```text
m_soft = m_donor - rho_b * V
s³ = 1 + m_soft / (rho_s * V)
S = V * [rho_b + rho_s * (s⁵ - 1)] * C
I_COM = trace(S) * identity - S
```

This defines a positive mass distribution and derives its inertia physically. It is not a numerical projection of the invalid donor tensor. The bone hull already fills canals and uses one density, so its bone mass is a proxy, not a measured cortical/trabecular split.

The derived shell scales are 2.32–2.61 times the bone hull. C1's full transverse proxy width is about 200 mm. Eighteen segment pairs have overlapping envelope AABBs; those intersection volumes are reported only as upper bounds on actual envelope overlap. These broad shells represent distributions of uniquely debited lumped masses. They do **not** claim mutually exclusive tissue compartments, an anatomical neck surface, or a validated continuum density field. A future anatomical neck partition can replace them without hiding this prior.

The head is a separate homogeneous 3.8 kg composite inside the craniofacial convex envelope; brain, bone and other contents are not counted again. Its effective density is approximately 1301 kg/m³. The jaw separately assigns 0.2 kg to its convex envelope, effective density approximately 651 kg/m³. The jaw's filled mouth concavity and the head's unresolved density heterogeneity/outer soft surface remain explicit. No extra canonical tissue mass is added to these donor segment masses.

## Result in the actual native torso frame

All masses/COMs/tensors are transformed through the one proper rigid canonical-to-torso map. The full tensors are retained, not reduced to scalar spheres.

| Quantity | Result |
|---|---:|
| Original native torso | 27.6546769653 kg |
| Nine transferred neck/head/jaw bodies | 7.4189195682 kg |
| Residual torso | 20.2357573970 kg |
| Residual COM, torso frame | `[-0.02393148745, 0.25211890208, -0.00100757564]` m |
| Minimum residual second-mass-moment eigenvalue | 0.3501680277 kg·m² |
| Reconstructed mass residual | 0 kg |
| Largest first-moment residual component | 9.60×10⁻¹⁷ kg·m |
| Reconstructed inertia Frobenius residual | 1.11×10⁻¹⁵ kg·m² |

Residual torso COM inertia in kg·m²:

```text
[[ 1.124308893927, -0.030301303403, -0.000492531928],
 [-0.030301303403,  0.756340771974,  0.005261166017],
 [-0.000492531928,  0.005261166017,  1.073292933344]]
```

The debit preserves mass, first moment and inertia about a common origin using parallel-axis terms, then recenters the residual tensor at its new COM. Every transferred and residual tensor must have positive inertia eigenvalues and nonnegative eigenvalues of `trace(I)/2 * identity - I`. This latter check rejects the original cervical diagonal `[0.02,0.08,0.02]`, despite its positive definiteness.

The four density corners at bone 1500/2200 kg/m³ and soft 900/1100 kg/m³ all preserve a physically admissible residual torso; its minimum second-moment eigenvalue remains 0.34992–0.35035 kg·m². These are sensitivity scenarios based on generic prior ranges, not empirical confidence intervals. Registration and head-envelope uncertainty are **not** covered by that density range.

## Opt-in integration recipe and unresolved work

The manifest contains the nine proposed inertial records in the **current torso frame**, the replacement torso record, source identities and a proposed 31-body/170-muscle MASI augmentation count. These are inertial inputs, not an executable `.osim` augmentation. `native_activation_allowed` remains false because no native/body-frame integration was performed in this task.

Before generating a loadable model, register the donor root and each body's frame, then express each COM and tensor in its new body frame. Preserve the donor joint tree, 24 coordinates/18 couplers, 78 unique muscle routes and original 92 actuators. Resolve clavicle/scapula/thoracic muscle attachment frames without adding duplicate donor shoulder/spine masses. Rebind head/neck skin and contact to the new bodies, preserving the global force/gravity/bed frame, and verify the instantiated native tensors match this recipe exactly.

The inherited global registration still has 62.356 mm RMS correspondence residual. A conservative inertia ledger does not turn that into an anatomically accurate joint map. Head contact geometry, overlapping soft-shell priors and the separate bed/support feasibility issue need their own acceptance. No native solver, model load, equilibrium, activation or model-builder change was used to obtain this result.

## Bounded verification

```sh
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m unittest scripts.verify_cervical_inertia
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.build_cervical_inertia --verify data/research/cervical_inertia/v2/manifest.json
```

Six analytic tests cover exact box volume/moments, rigid-transform covariance, core/shell mass and second moments, exclusive torso reconstruction, rejection of positive-definite but impossible tensors, and reflection rejection. The second command checks all archived input/geometry identities, recomputes every convex geometry prior and transferred tensor, reconstructs the torso debit, and verifies the inertial-only activation guard. All eight execution inputs are gzip-archived byte-for-byte; source and generated code are retained with their raw-byte hashes. The materializer requires a fresh workspace-owned output directory.
