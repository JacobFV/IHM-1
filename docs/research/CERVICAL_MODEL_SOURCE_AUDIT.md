# Cervical/head donor evidence and compatibility

Two actual research models are retained in `data/research/cervical`, replacing the source-search gap left by the five-muscle `Neck3dof_point_constraint.osim` engine fixture. The **MASI model is an original author deposit**. The Mortensen model is a clearly identified **reference copy**, with the original project's metadata also retained. Their joint/path data provide a concrete cervical extension, but **neither retained model has physically admissible cervical inertia tensors as written**. This is a blocking defect for unchanged dynamic import, not a claim that cervical articulation cannot be built.

No native process, simulation, shared model edit, factory edit, source geometry modification, or full geometry load was performed. The XML-only audit uses less than 1 GiB address space. Acquisition retained about 1.27 MB of model and provenance inputs before local audit/source copies; the 172 MB original HYOID package was not downloaded.

## Acquired identities and provenance

| Source | Actual retained bytes | Identity / license | Use and limits |
|---|---:|---|---|
| [MASI v2, authors' Figshare deposit](https://doi.org/10.6084/m9.figshare.4249808.v2) | 631,598 | SHA-256 `191d1bfdd328273d644f5822d85457e5e8fefd5037459c17d7f24b27d8b2681a`; Figshare declares MIT | Complete `.osim`; file size and provider MD5 verified. Credits Cazzola/Trewartha. XML publication field is a placeholder, but the linked primary article explicitly identifies this DOI. |
| [Mortensen2018 reference copy](https://github.com/mjhmilla/kinematicPassengerModel/blob/b0eb96127ca07dea0266764e837faeaa397092b5/models/reference/Mortensen2018.osim) | 494,571 | SHA-256 `f85b25349a0132938769c18a9f1e6c3049b6e1395c12ffb693e03f38e51acf0b`; pinned Git blob `41fb3623066aff0f7ce9fc4f5c5f13ae07589fe1`; repository MIT notice retained | Model itself credits Jonathan Mortensen and the 2018 paper. Reference copy only, not the repository's assembled passenger model. Its equivalence to the original latest release is unverified. |
| [Original HYOID project](https://simtk.org/projects/neckdynamics) | Metadata/license only | HYOID MIT; project separately identifies modified MASI files in its publication subfolder as CC BY-NC 3.0 | Release 1.2 dated July 26, 2019 says skull inertia was updated. Original ZIP is 172 MB; its confirmation endpoint redirected to account login. No authentication bypass or package acquisition. |

Original XML files are unchanged. `sources.json` records each acquired file's URL, size and SHA-256; retained HTML includes the actual HYOID license and login response. Do not infer that the secondary copy includes the 2019 update, or apply its repository license to every upstream derivative. Referenced visual geometry was not acquired: 90 distinct geometry names in MASI, 15 in the reference HYOID file. Their muscle routes use fixed `PathPoint` objects and no wrapping objects; those XML path coordinates are present, but mesh-based anatomical registration still needs the geometry.

The [Cazzola et al. primary article](https://doi.org/10.1371/journal.pone.0169329) distinguishes generic MASI from a rugby-specific model. It describes cervical muscle strength scaling and tests of motion, moments and activation. This acquisition selects generic MASI, not the rugby player's personalized inertia or strength. The [Mortensen et al. primary article](https://doi.org/10.1371/journal.pone.0199912) identifies inadequate upper-cervical flexion in earlier models and adds hyoid muscles/passive elements. Those results do not establish supine equilibrium after transplantation into IHM.

## Exact XML findings

`audit.json` retains every selected body's mass, COM, tensor, joint frame, coordinate range, rotation axis, coupling function, muscle identity/parameters/path-body dependencies, and passive bushing constants. Counts below are measured from acquired XML, not inferred from publication summaries.

| Property | MASI | Mortensen reference |
|---|---:|---:|
| XML document version | 30000 | 30000 |
| Cervical bodies | `cerv7` through `cerv1` | same |
| Head bodies | `skull`, welded `jaw` | same |
| Neck coordinates / couplers / independent coordinates | 24 / 18 / 6 | 24 / 18 / 6 |
| Muscle law / path elements | 78 Thelen2003 / 198 | 72 Millard2012Equilibrium / 186 |
| Other force elements | 23 coordinate actuators | 8 bushings, 1 point actuator, 1 torque actuator |
| Neck + skull + jaw mass | 7.4189195682 kg | same |
| Skull / jaw mass | 3.8 / 0.2 kg | same |

The six master coordinates are lower-neck `pitch2`, `roll2`, `yaw2` and upper-neck `pitch1`, `roll1`, `yaw1`. The upper master coordinates are implemented at `aux2jnt` (C2–C1), with a coupled C1–skull joint; do not rename them as three independent anatomical rotations at every vertebra. Eighteen other rotations are explicitly dependent. Jaw is welded and no distinct hyoid body or swallowing DOF exists in either XML. The HYOID copy represents its named hyoid routes using points on the existing bodies; it does not provide a deformable hyoid or complete swallowing apparatus.

Units are meters, newtons and rotational radians; gravity is `[0,-9.80665,0]`. Source neutral orientation is Y superior, X anterior, Z right (the anatomical direction inference agrees with right/left attachment signs). The neck `SpatialTransform` lists rotations in order about Z, X and Y for pitch, roll and yaw. Preserve the source transformation sequence, functions, parent/child frames and radians; coordinate labels alone are insufficient.

Important frame discrepancy: MASI `auxt1jnt` anchors C7 at `[-0.0475,0.5343,0]` in its `spine`, while HYOID uses `[-0.0475,0.3613,0]` in its `spine`. The difference is 0.173 m and matches MASI's welded ribcage/torso offset. This is a frame-origin difference to resolve with anatomical registration, not a head translation to copy directly into the current torso. Both use child anchor `[0.00684,-0.005655,0]` in C7. All other cervical joint offsets are retained exactly in `audit.json`.

## Blocking inertia finding

All seven cervical bodies in **both** acquired files have diagonal tensor `[0.02,0.08,0.02] kg m²`, with zero off-diagonals. Although positive definite, this is impossible for any nonnegative mass distribution: `Iyy > Ixx + Izz`. Equivalently `S = trace(I)/2 identity - I` has eigenvalues `[-0.02,0.04,0.04] kg m²`. Rotating frames or scaling every inertia by one positive factor cannot fix this inequality. The source data are preserved, and no unreported repair is made.

The held OpenSim implementation, copied as `OpenSim_Body.cpp` with its license header, catches invalid tensors and substitutes a sphere with diagonal `norm([Ixx,Iyy,Izz])/sqrt(3)` (lines 146–163). For these numbers the replacement is about **0.048989795 kg m² on each axis**. This is source inspection, not an observed native load in this task. A successful native load can therefore conceal a changed inertial model; compare the instantiated tensor with the XML before accepting any donor.

The skull and jaw pass the algebraic tensor check. That check establishes realizability only, not anatomical accuracy. For example, the 0.2 kg jaw's `[0.04,0.02,0.04] kg m²` deserves independent anatomical review despite satisfying the inequalities. The later original HYOID release advertises a skull update, which is an additional reason not to promote this reference copy as the latest dynamic calibration.

## Concrete compatibility path for the current 92-muscle assembly

1. **Choose and freeze one donor version.** MASI provides first-party full neck geometry/path parameters using the same Thelen class as the existing twelve Arm26 actuators. HYOID provides a more complete flexion/passive-law candidate, but use an authenticated original release or keep its reference-copy limitation explicit. Do not combine both overlapping neck muscle sets. A MASI addition would total 170 muscles; the retained HYOID set would total 164, assuming each original 92 remains exactly once. These are intended catalog counts, not validated assemblies.
2. **Register T1, cervical joint centers and skull landmarks into the current native torso.** Solve and retain a proper local donor-to-torso transform against actual anatomical landmarks, including scale uncertainty. Then carry the result through the existing single world transform. Use the exact scaled `subject_with_arms.osim` identity, not unscaled Rajagopal geometry as a substitute. Keep world gravity, force points, skin and bed in that same frame. Donor display mesh scale alone does not scale joint offsets, COMs or muscle lengths.
3. **Partition exclusive inertial ownership before adding bodies.** Current torso is 30.3832390808 kg, COM `[-0.03,0.32,0] m`, COM inertia diagonal `[1.669987476,0.8556633016,1.6211733287] kg m²`. Its original model hash is in `audit.json`. Register and reidentify physically admissible head/neck tensors; subtract their mass and full origin inertia from the existing torso, then recover the residual torso COM/tensor. For donor masses used without rescaling, the arithmetic residual mass would be 22.9643195126 kg, **not an accepted partition**. Reject negative mass or invalid residual second moments. Do not import the donor ground, 70 kg HYOID spine, auxiliary 1 kg shoulder/ribcage bodies, or MASI's other full-body inertias on top of existing mass.
4. **Preserve the eight cervical joints, eighteen couplers and welded jaw.** Convert version-30000 body-local joints into version-40500 joints/offset frames using the held OpenSim migration semantics; bind every parent/child socket and geometry path explicitly. Do not replace the six master coordinates with ad hoc head contact rotations. Retain coordinate-dependent transform splines; zero/default coordinate values are not interchangeable with zero joint rotations.
5. **Resolve shoulder and thoracic muscle anchors explicitly.** Both donors require `spine`, bilateral clavicles/scapulae and cervical/skull bodies; MASI also references `torso`, HYOID `ribcage`. Current native skeleton has no separate scapula/clavicle owners. For a first fixed-shoulder reference, transform those source anchors into named massless offset frames on the current torso and state the welded-shoulder approximation. A full moving shoulder rhythm needs its own DOFs/couplers and mass partition. The donor scapular paths cannot follow current humeral elevation merely by renaming them `torso`.
6. **Keep native actuator and passive-law ownership.** Add donor muscle routes under unique names without overwriting original lower-limb fitted paths or twelve Arm26 paths. The current native adapter must fit only its original eighty muscles; new cervical paths must remain the acquired geometry paths. For HYOID, preserve the eight bushing laws and account for their stored/dissipated energy; reject unplanned reserve point/torque actuators. Neither MASI's 23 whole-body actuators nor HYOID's experiment actuators should silently become external support. Keep joint and muscle law changes separately receipted.
7. **Rebind skull/neck skin and contact to the new material owners.** Registration currently maps those canonical points to torso. Change ownership using the accepted anatomy map, preserving full skin coverage and exclusive mass. Rebuild current-pose surface quadrature and check passive generalized equilibrium, force/moment balance, compression domains and constraint residuals before a sustained native support run. Adding cervical DOFs does not remedy the separate rigid-bed travel and support-polygon limitations in `SUPINE_GEOMETRIC_FEASIBILITY.md`.

For step 3, in one registered torso frame let `P(r)=||r||² identity-r rᵀ`. Compute `m_r=m_0-Σm_i`, `c_r=(m_0 c_0-Σm_i c_i)/m_r`, and `I_r=I_0+m_0 P(c_0)-Σ(R_i I_i R_iᵀ+m_i P(c_i))-m_r P(c_r)`. This preserves mass, first moment and origin inertia at the reference configuration. Apply both positive-inertia and second-mass-moment checks; neither an arbitrary mass subtraction nor a nearest positive-definite matrix is anatomical calibration.

## Verification and remaining acquisition

Run:

```sh
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python data/research/cervical/audit_sources.py --check
```

This checks retained hashes, provider MD5, pinned Git blob identity, XML counts, six-master/18-coupler structure, the reproduced seven-tensor failures in both sources, current target identity, and exact extraction reproducibility. It performs no forward simulation and asserts no equilibrium. The next material inputs are corrected, source-supported cervical inertia/COM distributions and donor-to-current anatomical landmarks; the original HYOID release and geometry remain useful acquisition targets. No model value has been replaced or promoted by this report.
