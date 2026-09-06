# Source-only cervical registration recipe

The opt-in artifact is `data/research/cervical_registration/v2/recipe.json`.
It assembles nine new body descriptions, 24 coordinates, 18 cervical couplers,
and 78 muscles with all 198 source path points. It is a JSON component recipe,
**not a loadable OpenSim model**. `native_activation_allowed` remains false.
No default model, controller, skin owner, or native process was changed.

The recipe binds the exact initial 22-body/92-muscle factory model through its
SHA256 receipt, and the unchanged prior manifest from
[CERVICAL_INERTIAL_PRIORS.md](CERVICAL_INERTIAL_PRIORS.md). This is the assembled
77.6122029 kg model whose torso is 27.6546769653 kg. The 30.383239 kg torso in the
source input is before scaling; it is not a different production distribution.
No new normalization is applied. Validation rejects any other model bytes.

## Frames and geometric evidence

MASI v2 donor coordinates are meters/radians, X anterior, Y superior, Z right;
cervical rotations follow body-fixed Z/X/Y transform axes. The reference sets
all 24 cervical coordinates to zero, satisfying all 18 cervical couplers within
1e-14 rad numerical tolerance. Original joint XML, coordinate limits/defaults,
and constraint functions remain available verbatim as parsed XML text. Recipe
reference values explicitly override the tiny nonzero donor experiment defaults.

Pure Python reference kinematics accept CustomJoint and WeldJoint, affine
LinearFunction, Constant, and two-knot SimmSpline only. Unknown joints, reverse
joints, unrecognized axis order, nonunit axes, nonzero mobility translations,
and unsupported function forms fail. No module/model initialization is used.
The transformation convention still needs a native perturbation comparison.

Eight held OpenSim library meshes have exact MASI filenames: cerv1–cerv6,
skull, and jaw. Their raw bytes and hashes are retained. Filename matching does
not authenticate their mesh revision as the one used by MASI v2. Their separate
geometry license is not established by the donor model MIT license. These held
research inputs should not be represented as an authenticated source atlas.
MASI references `rotatedcerv7.vtp`, which is absent. The available `cerv7.vtp`
is deliberately not substituted. Nonidentity display scaling/transforms and
multiple body meshes are unsupported and rejected by this builder.

A single proper-rigid fit maps donor-spine coordinates to the frozen target
torso. Correspondences are convex-volume centroids of the eight held meshes
and the canonical bone hulls already retained in the inertia recipe. They are
engineering proxies, not measured homologous landmarks. Relative donor joint
frames and segment dimensions remain unchanged; there is no scale fitting.

| Body | Centroid residual (mm) |
| --- | ---: |
| C6 | 2.071 |
| C5 | 2.752 |
| C4 | 9.242 |
| C3 | 10.071 |
| C2 | 8.198 |
| C1 | 36.229 |
| Skull | 7.949 |
| Jaw | 11.557 |

RMS is 14.896 mm. The source correspondences are nearly sagittal: singular
values are 126.188, 69.387, and 0.138 mm. The fit enforces determinant +1 even
though unconstrained least squares prefers a reflection on the weak third
axis. The resulting donor-right axis aligns with target-right (dot product
0.99998). This is numerical observability evidence, not validation of joint
centers. The preexisting canonical global registration also carries its
62.356 mm COM/envelope proxy RMS. The new residual does not erase that error.

## Explicit component and mass mapping

Each `bodies` entry supplies the source joint XML, unique `cervical_` identity,
parent body, reference body-to-target-torso transform, local COM/inertia, and
canonical-geometry-to-new-body transform. At the root only, the recorded donor
spine-to-target-torso transform composes with the original C7 parent joint
frame. Every remaining joint keeps its original relative frames. A future XML
converter must rewrite source coordinate names inside axes/couplers, apply the
recipe reference coordinate values, and convert old XML socket/function
structure explicitly; simply appending old XML is not supported.

The prior mass distributions remain in exactly the same target-torso positions.
Only their coordinate expression changes using each new body's inverse rigid
transform. The replacement torso loses the entire 7.4189195682 kg donor
neck/head mass, leaving 20.2357573970 kg. No donor spine, torso, clavicle, or
scapula inertia is imported. Recombination gives mass error 0 kg, COM error
2.94e-19 m, and COM-inertia Frobenius error 1.26e-15 kg m². A future articulated
pose naturally changes the total COM/inertia; conservation here is the explicit
reference partition, not a claim of pose-invariant inertia.

## Muscle evidence and approximations

All 78 donor Thelen2003 muscle parameter records and fixed paths are retained.
Cervical/skull path points keep their source body-local coordinates. External
points are transformed to fixed target-torso locations: 48 spine, eight torso,
six clavicle, and four scapula points. Donor shoulder frames use source default
coordinates and fixed frame orientations; the recipe does **not** claim the
shoulder couplers have been assembled or that these anchors match target
shoulder anatomy. Their transformations and coordinate values are explicit.
No shoulder mass or new shoulder DOF is added. Moving/conditional points,
wrapping, unresolved owners, or different muscle types fail this recipe.

Unchanged donor strengths/rest lengths are source evidence, not a validated
operating point after registration. Existing 92 muscles may overlap in function
with some imported neck/shoulder forces. Moment arms, force-length states,
strength overlap, and shoulder motion need explicit acceptance. Canonical bone
ownership mappings do not automatically authorize skin/soft-tissue/contact
ownership. Current torso head/neck display and contacts must be repartitioned
before physical use, with no duplicate head mass, geometry, or contact forces.

## Verification and remaining work

`OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.verify_cervical_registration`

Five checks cover analytic rigid recovery, collinear/nonfinite rejection,
reflection handling, supported function/joint scope, raw receipt integrity,
independent XML torso conservation, changed-target/changed-mass rejection, and
muscle owner/location mapping. The source builder uses bounded geometry and
retains the exact input bytes. Rebuild into a fresh directory with
`python -m scripts.build_cervical_registration --output <fresh workspace path>`.

Before native acceptance: resolve C1 geometry/reference mismatch and missing C7
mesh, verify source geometry provenance, quantify surface/joint alignment,
convert XML with explicit sockets and names, validate fixed anchors or shoulder
mechanics, map new skin/contact owners, then run native frame/coupler/moment-arm
and supine support checks in an opt-in model. The artifact supplies concrete
inputs for those steps without concealing their unresolved accuracy.
