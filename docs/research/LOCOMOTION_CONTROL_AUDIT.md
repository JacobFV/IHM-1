# Native locomotion, control, and soft-tissue coupling audit

Audit date: 2026-09-05. Read-only inspection of retained sources, build configuration, binary symbols, manifests, and existing results. No new walking simulation was executed; no implementation or retained source was modified. Paths below are repository-relative.

## Decision

The repository has enough native OpenSim infrastructure and retained model/data assets to start a defensible forward-dynamics locomotion program. It does **not** currently have integrated, self-supporting neural walking or a whole-body deformable contact simulation. The strongest retained whole-body starting point is `data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/`, rather than the existing one-pose adapter or the smaller 2D demonstration.

Keep OpenSim responsible for articulated skeleton, muscle dynamics, and initially foot-ground contact. Couple it through explicit rigid transforms and equal/opposite wrenches to independently resolved deformable domains. Establish free forward integration and perturbation sensitivity before interpreting visual motion as controlled walking. Broader organ/tissue coverage should use the same exchange contract, while individually recording unresolved anatomy, material, boundary, and numerical coverage.

## What currently executes

`ihm/native/opensim_backend.py` launches `data/runtime/opensim/native_opensim_export`, implemented in `scripts/native_opensim_export.cpp`. It always loads `Models/Rajagopal/Rajagopal2016.osim` from retained opensim-models. It changes one independent rotational coordinate relative to its default, sets every muscle to a common activation, equilibrates muscle states, and calls `realizeDynamics`. It does not construct an OpenSim `Manager`, advance simulation time, solve standing balance, add ground contact, or run a neural controller.

The permitted adapter perturbation is +/-0.15 rad and activation 0.01–0.5. Those limits belong to this exporter, not the physical capabilities of OpenSim. Exported muscle length/fiber length/path points/moment arms are meters; force is N; fiber velocity is m/s; rotational coordinates are rad and translational coordinates meters. `body_transforms_ground` maps body-local coordinates into OpenSim Ground. The force scope explicitly says external skeletal force balance is not solved. A displayed walking trajectory assembled from successive prescribed poses would not change that status.

The source Rajagopal2016 XML has 22 bodies, 39 coordinates, 80 Millard2012 equilibrium muscles, and 17 coordinate actuators. Eight coordinates are locked (bilateral subtalar, MTP, wrist flexion/deviation); dependent patellar coordinates further distinguish coordinate count from independent DOFs. There is no ContactGeometrySet content and no controller. Gravity is `(0, -9.80665, 0)` m/s²; the pose export realizes dynamics but does not integrate gravity-driven movement.

Native installation evidence:

- `data/runtime/opensim/build_manifest.json` pins OpenSim `86b30588374650fbaf012a345a836a64f6855522` and Simbody `59c6e7b89b3bdf266a2f3d54c611599d964205f7`; architecture is aarch64.
- Read-only SHA-256 recomputation matched the manifest for the exporter and every listed OpenSim/Simbody `.so` library.
- `data/runtime/opensim/opensim-build/CMakeCache.txt` reports `BUILD_API_ONLY=ON`, `BUILD_PYTHON_WRAPPING=OFF`, `BUILD_JAVA_WRAPPING=OFF`, `OPENSIM_WITH_CASADI=OFF`.
- `install/opensim/bin` is empty. Do not propose `opensim-cmd run-tool ...` as an already available command. Native C++ APIs are available through `libosimTools.so`, `libosimSimulation.so`, `libosimActuators.so`, `libosimAnalyses.so`, and `libosimExampleComponents.so`.
- `libosimMoco.so` exists, but that does not establish an enabled CasADi/IPOPT optimization backend. The retained walking examples call `MocoCasADiSolver`; rebuilding/enabling its dependencies is needed for those optimization solves.
- `nm -DC .../libosimExampleComponents.so` exposes `OpenSim::ToyReflexController::computeControls(...)`.

Existing corrected static mechanics uses experimental variant `wrap_8_0.0005_cache`, not upstream by default. `data/derived/opensim/native_corrected/baseline/execution.json` identifies it. Its static verification reports body-transform maximum error `3.3306690738754696e-16`. The existing corrected wrapping tests do not establish dynamic wrapping stability across a full gait cycle. Any forward adapter must record which library it actually loads and repeat path/moment-arm checks over its full operating range.

## Retained locomotion resources and their actual meaning

### Primary whole-body starting point: native 3D walking example

Directory: `data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/`.

| Retained file | Usable content and limitation |
|---|---|
| `subject_walk_scaled.osim` | Rajagopal-derived model: 22 bodies, 33 coordinates, 80 Millard muscles. Separate scaled subject/model from both Rajagopal2016 baseline and canonical BodyParts3D anatomy. |
| `subject_walk_scaled_ContactGeometrySet.xml` | Twelve contact spheres (six per foot) and one floor half-space. |
| `subject_walk_scaled_ContactForceSet.xml` | Twelve SmoothSphereHalfSpaceForce components. Enables simulated foot-ground forces rather than supplying measured forces as external loads. |
| `subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml` | Added passive coordinate stiffness/damping; inspect and retain parameter provenance. |
| `subject_walk_scaled_FunctionBasedPathSet.xml` | Fitted polynomial muscle paths used by optimization. These are an alternative to geometric wrapping; they do not supply anatomical surface muscle geometry. |
| `coordinates.sto`, `grf_walk.mot`, `grf_walk.xml`, `markers_walk.trc`, `id_walk.sto`, `electromyography.sto` | Retained kinematic, force, marker, inverse-dynamics and EMG reference data. EMG normalization/processing is documented in README; amplitude is not direct physiological excitation calibration. |
| `example3DWalking.cpp` | Assembles contact model, then solves torque-driven and muscle-driven kinematics/GRF tracking for 0.48–1.61 s. Not already built into the IHM adapter. |
| `exampleMocoTrack.cpp`, `exampleMocoInverse.cpp` | Different workflows; for example, MocoTrack explicitly adds measured external loads in lieu of ground contact. Do not treat every similarly named example as a contact-driven forward model. |

The main example adds contact geometries and raises the spheres 0.02 m relative to the supplied locations. Copying XML alone therefore does not reproduce the prepared model. It initializes minimum muscle activation/excitation to zero, adds coordinate forces, and finalizes connections. Its optimization muscle branch ignores tendon compliance, replaces Millard muscles with DeGrooteFregly2016 muscles, ignores passive fiber forces, widens the active force-length curve by 1.5, and replaces geometry paths with fitted functions. These simplifications must be explicit in the resulting model identity; do not claim unchanged Rajagopal/Millard dynamics.

The source README is a useful provenance guide, but some descriptions (e.g. welded MTP coordinates) must be checked against actual model topology before assigning DOFs. Source comments are not a substitute for exported coordinate/constraint metadata from the constructed model.

This is the preferred basis for 3D contact gait and full-body articulation. It still lacks an established perturbation-rejecting reflex gait policy and volumetric tissue interaction model. Tracking optimization produces dynamically constrained trajectories; independent forward integration is a separate verification stage. The distinction between Moco optimization and a simulation controller is also described in the [official Moco developer guide](https://opensim-org.github.io/opensim-moco-site/docs/0.3.0/mocodevguide.html).

Additional native-read findings:

- No saved state/control solution for this 3D endogenous-contact example was found in the retained opensim-core or opensim-models trees. Its local `.sto` files are `coordinates.sto`, `id_walk.sto`, and `electromyography.sto`. The coordinates table has 1,350 rows, 40 columns, and `inDegrees=no`; it is not a full muscle state/control trajectory. The example writes solution files when run; those names in source code are not existing results.
- Total base-model segment mass is **85.26984854173146 kg**, compared with **75.337 kg** for retained Rajagopal2016 and the canonical mass target **77.1107029 kg**. Establish one explicit physical subject/mass partition before feedback coupling.
- Actual base model has 13 torso/arm CoordinateActuators with optimal force 50 (N m for their rotational coordinates); it does not contain a populated ControllerSet. Upper-body motion is therefore partly torque-actuated rather than driven by anatomically resolved upper-body muscles. Actual MTP coordinates exist and are unlocked despite the README's weld description. Two knee beta coordinates are coupled; export independence/constraint status natively.
- Contact spheres have radius 0.035 m and attach to `calcn_r/l` or `toes_r/l`. All twelve contact forces specify stiffness `1e6 N/m²` (a contact material modulus, **not N/m**), dissipation `2 s/m`, static/dynamic friction `0.8`, viscous parameter `0.5`, transition velocity `0.2 m/s`, Hertz smoothing `300`, and Hunt-Crossley smoothing `50`. Preserve each parameter's native definition rather than inferring all dimensions from the word coefficient. Report native forces in N, torques in N m, and penetration in m. The source's +0.02 m sphere-height preparation still applies.

No optimization is needed to execute a short physically free integration with this assembled model and chosen controls, or to develop a feedback controller that supplies controls online. Optimization or another valid control-identification method is needed to obtain defensible sustained tracking gait when no compatible saved policy exists. Do not transplant gait2392 CMC histories by muscle name: the models, muscle laws, paths, masses, contact support and initial states differ.

### Existing gait2392 forward reference: shortest integration regression

Directory: `data/raw/anatomy/opensim-models/source/Pipelines/Gait2392_Simbody/`.

`subject01_Setup_Forward.xml` defines an actual ForwardTool integration from 0.8 to 1.18 s. It uses:

- `subject01_simbody_adjusted.osim` (also retained under `Models/Gait2392_Simbody/` and `OutputReference/`): 12 bodies, 23 coordinates, 92 Thelen muscles, no built-in contact.
- `gait2392_CMC_Actuators.xml` appended to its force set.
- `ResultsCMC/subject01_walk1_controls.xml` through a ControlSetController.
- `ResultsCMC/subject01_walk1_states.sto` for initial conditions and correction targets.
- Enabled CorrectionController with `kp=16`, `kv=8`.
- `subject01_walk1_grf.xml`, supplying measured GRF/COP/torque from `subject01_walk1_grf.mot` to `calcn_r` and `calcn_l`, expressed in Ground.
- Retained reference states/controls under `OutputReference/ResultsCMC/` and output files under `OutputReference/ResultsForward/`.

The retained forward states file has 1,734 rows spanning 0.80000000–1.18000000 s. Outputs include states, controls, kinematics q/u/dudt, body positions/velocities/accelerations, and actuator force/speed/power. Kinematics analysis requests degrees; do not silently interpret every `.sto` rotational column as rad. Use each file's header and originating analysis. Native state variables should be exported with explicit units and full component paths.

This is more than prescribed coordinate replay: positions/speeds are integrated under supplied controls, measured loads, and corrective feedback. It is less than autonomous contact gait: ground reactions and desired correction targets come from prior data, and corrective/residual actuation can support the motion. First rerun with original settings, then disable correction and separately perturb/remove measured loads. Preserve all versions and label them accurately. Resolve relative paths in a fresh run directory; the setup expects `ResultsCMC`, whereas retained reference data sits in `OutputReference/ResultsCMC`. A wrapper must copy/rewrite references deliberately rather than writing into the source tree.

### Smaller contact examples: useful numerical fixtures

`data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example2DWalking/2D_gait.osim` has 12 bodies, 10 coordinates, 18 DeGrooteFregly muscles, one coordinate actuator, four SmoothSphereHalfSpaceForce components, four spheres and one half-space. Forces/muscles are serialized under generic `components`, so counting only ForceSet gives a false zero. `referenceCoordinates.sto` has 50 rows from 0 to 0.47008941 s, rad/m with `inDegrees=no`; accompanying `referenceGRF.sto` and XML support force comparison. It does not contain a solved excitation policy or controller.

**Known source defect:** `example2DWalking.cpp` states that the gastrocnemius path does not cross the knee and that the model should not be used for research. Use only as an explicitly limited integration/contact fixture unless repaired against source anatomy and verified for nonzero knee moment arm and anatomical routing. It is not the recommended final whole-body seed.

`Models/Gait10dof18musc/subject01.osim` retains 18 Millard muscles and 10 coordinates but no ground contact/controller. Tutorial metabolic CMC outputs are useful reference controls, not evidence of a transferable autonomous gait policy.

`Models/WalkerModel/{WalkerModel,WalkerModelTerrain}.osim` retains HuntCrossley contact and passive walking resources; the controller sets are empty. Python's strict XML parser rejects legacy tags such as `HuntCrossleyForce::ContactParametersSet` at lines 640/632 respectively. This is a strict-parser compatibility finding, not proof that OpenSim's legacy model loader rejects the files. Native loading needs a controlled probe. The source builder `Tutorials/Building_a_Passive_Dynamic_Walker/DynamicWalkerBuildModel.cpp` is another retained passive/contact construction reference. Passive downhill walking is not muscle-driven neural walking.

## Reflex and locomotor control

Actual retained native implementation: `data/raw/mechanics/opensim-core/OpenSim/ExampleComponents/ToyReflexController.{h,cpp}`. It computes additive control proportional to positive musculotendon lengthening speed:

`u_add = gain * max(lengthening_speed, 0) / (optimal_fiber_length * max_contraction_velocity)`.

The denominator is m/s (maximum contraction velocity is fiber lengths/s); gain and excitation are dimensionless. It accepts muscle actuators only. It adds to existing controls. The implementation itself has no delay buffer, force feedback, sensory noise, stance/swing state machine, foot placement, balance policy, learned/optimized gait parameters, or explicit saturation in this formula. The header calls it demonstration-only. Connecting it to every muscle is not a defensible claim of neural walking.

Shortest useful controller progression (proposed, not implemented):

1. Native forward plant with reproducible initial states and supplied excitation histories; record actual free/prescribed/locked coordinates and all externally imposed loads.
2. Contact gait tracking policy with feedback on state, finite torque/excitation bounds, and explicit correction/residual work. Use optimized feedforward controls as an initialization if a compatible Moco build is available.
3. Add delayed force/length/velocity feedback, contact sensing and gait phase logic; identify gains against retained data and evaluate perturbed trials excluded from fitting. Describe this as a reduced neuromechanical controller until validated, not a reconstructed human nervous system.
4. Progress from sagittal gait to lateral balance, foot placement, trunk/arm participation, variable speed, transitions and uneven contact. Short fixtures verify machinery but do not define the eventual anatomical scope.

Controller history must update on accepted integration steps or a well-defined discrete sampling clock. OpenSim `computeControls` may be called repeatedly on trial states by an adaptive integrator; advancing mutable delay history per evaluation breaks reproducibility and causality. Archive controller states/history in checkpoints, enforce signal units and bounds, and separate physiological delay from communication latency.

Installed controller headers expose `Controller::addActuator(const Actuator&)`, `setActuators(...)`, `setEnabled(bool)`, and `computeControls(const SimTK::State&, SimTK::Vector&) const`. Use a concrete Controller for online state feedback and add actuator controls using the actuator's own indexing. For supplied control functions, the current `PrescribedController::prescribeControlForActuator(const std::string&, const Function&)` copies the function; the old `(int, Function*)` overload now throws. Build against these installed headers rather than older examples. `Model::setControls` alone is not a persistent online policy. Tissue reactions should enter a correctly implemented OpenSim force component/native Simbody force interface; Moco's `DiscreteForces` wrapper explicitly says it is internal to Moco and is not a public general-purpose coupling recommendation.

## Registration to canonical anatomy and current mechanics

`scripts/build_opensim_display.py` applies the fixed rotation

`R = [[0,0,-1], [0,1,0], [1,0,0]]`

and a translation equal to minus the rotated rest-geometry bounding-box center. Thus `p_display = R p_ground + t`; source +Y remains vertical. Its manifest explicitly says the source subject differs from BodyParts3D. It stores body-frame IDs and original attachment/via points. `ihm/spatial/opensim.py` only evaluates exact default poses/source spline knots and is not suitable for walking-time joint evaluation; use native OpenSim transforms.

`scripts/build_body_mechanics.py` maps source bodies to canonical bone groups, fits per-body axis-aligned bounding-box scales/offsets, and transforms each rest path point as `p_canonical = scale_body * (R p_ground + t) + offset_body`. Endpoints are projected to the nearest canonical bone surface. Wrap supports are inferred from nearest exported native path nodes, and wrapped paths are frozen at rest. Fmax and pennation are retained; fiber/slack lengths are scaled by registered path-length ratio; passive force is a different assumed reduced law.

Audited `data/derived/canonical/mechanics.json` has:

- 2,408 entities; 257 rigid bones; 2,145 affine soft solids; six numerical boundary carriers.
- 443 muscle actuators, of which 80 source paths are mapped and 363 are synthesized; 24 tendon/ligament paths; 3,341 support links.
- 22 registered source-body groups. Torso, pelvis, hands and toes each collect multiple canonical bones.
- Registered path-length ratios from 0.8652118927 to 1.1550413846; maximum recorded endpoint surface-projection distance 0.0383586546 m. These are geometric transfer diagnostics, not anatomical accuracy estimates.

Current `ihm/assembly/mechanics.py` integrates translations with constrained orientations, 2 ms substeps, no gravity, no collision, and no articulated joint rotation. Soft shape is one affine compressible neo-Hookean deformation per surface, solved quasistatically. The minimum-distance skeletal support tree is a geometric prior, not a measured articulated skeleton. A tissue moving because its support translation was prescribed is not an independently resolved organ deformation field.

For forward coupling, define a persistent body-local attachment map and explicit rigid transform chain at rest and each time. Do not re-fit AABBs each frame or use anisotropic scaling as if it were a rigid rotation. Many canonical bones will initially share a source segment's transform; mark missing articulation (individual ribs, vertebrae, fingers, etc.) in coverage. A world-space anisotropic map and a rotating rigid body do not generally commute; choose a calibrated local geometry map and test rigid-rotation objectivity. Recompute source wrapping with native dynamics, or disclose a validated reduced path model. Register the 3D walking subject separately; current Rajagopal2016 mapping cannot silently serve as its anatomical calibration.

## Soft tissue and native reference capability

`ihm/assembly/mechanics_backend.py::DeformableRegion` is a quasistatic tetrahedral compressible neo-Hookean solver with per-element material arrays and lumped nodal mass. It supports a restricted frictionless planar indenter acting on selected top-surface nodes. The nodal mass supports gravitational loading; the solver has no transient inertia/time integration, arbitrary inter-organ collision, articulated rigid boundary motion API, cavity flow, sliding fascia, or resolved anisotropic muscle fibers. It can represent assigned spatially varying parameters, but current capability should not be described as dynamic heterogeneous all-organ tissue mechanics.

Native FEBio 4.13 is retained at `data/derived/mechanics-reference/build/bin/febio4`, source revision `32ae206ff4881dfb54f62296cd1558e58ed9fcc6`. Existing `compression/benchmark.json` reports analytical stress error `4.4448711378208827e-10 Pa` and displacement error `1.78792204428e-16 m` for a homogeneous two-hex compression test. `contact/benchmark.json` reports a passed full-face frictionless rigid-plane/confined-block test, with decreasing penalty errors. These are useful numerical reference tests; neither validates anatomical soft tissues or locomotion-scale transient coupling.

Native FE extension needs explicit volumetric domains, density and mass partition, suitable per-domain constitutive laws/fiber fields, dynamic time stepping, contact pairs and attachment surfaces, prestress/equilibrium initialization, cavity pressure/volume boundary treatment, and meshing/quality checks. Different tissues need distinct validated descriptions; a single homogeneous spring or affine law is not adequate evidence for whole-body heterogeneity. FEBio can be a reference/execution choice, but its exact chosen dynamic/contact/material configuration needs new native evidence.

Minimum exchange contract:

| OpenSim to tissue | Tissue to OpenSim |
|---|---|
| Time in s; body `X_GB`; body angular/linear velocity rad/s and m/s; optionally acceleration | Resultant force N and moment N m per attachment body, expressed in a declared frame and about a declared origin |
| Registered attachment positions/velocities and muscle excitation/activation with explicit semantics | Nodal traction integration, interface work J, contact impulse N s, solver residuals and domain validity |
| Environmental gravity, active stress inputs, pressure boundaries Pa, reference material parameters | Tissue state/checkpoint, deformed geometry m, element strain/Jacobian, stress Pa and strain energy J |

Use `F_B = R_GB^T F_G` and rotate/shift moments consistently; sum nodal moments using `(x_node - x_body_origin) × f_node`. Check equal/opposite interface wrench and power/work conservation. If the geometric map includes scaling, use its Jacobian-transpose force transfer for virtual-work consistency rather than blindly copying force vectors. Avoid double-counting tissue mass already present in source segment inertias and avoid duplicate foot contact in OpenSim and FE. When FE takes over sole contact, replace the corresponding lumped ground force contribution.

Begin with one-way imposed skeleton motion to establish registration and boundary correctness, explicitly labeled one-way. Genuine interaction requires soft-domain reaction forces to change the next free skeletal state. Stiff partitioned coupling may require subiterations, rollback/checkpoint support and smaller communication steps; simply exchanging at renderer frame rate has no established stability guarantee. Whole-body coverage is a staged set of resolved domains under this contract, not an extrapolation from one successful patch.

## Exact next execution work and honest outputs

### Enabling retained CasADi/Moco optimization (optional for forward plant)

Pinned source prerequisites are recorded in `data/raw/mechanics/opensim-core/dependencies/CMakeLists.txt`: CasADi 3.6.5 with `WITH_IPOPT=ON`, IPOPT 3.14.16, ThirdParty-Mumps `releases/3.0.5`, METIS 5.1.0, and a Fortran-capable native build with BLAS/LAPACK. OpenSim CMake requires `find_package(casadi 3.6.5 REQUIRED)` and on UNIX `pkg_check_modules(IPOPT REQUIRED ipopt IMPORTED_TARGET)`. CasADi's IPOPT plugin and every transitive library must be loadable; merely finding `libcasadi` is insufficient. CMake minimum is 3.15 and the installed OpenSim targets require C++20.

Host checks found cmake, c++, make, pkg-config, autoconf and automake, but `command -v gfortran libtool` produced no paths. `pkg-config --modversion ipopt` fails. There are no retained CasADi/IPOPT binaries/config files in the searched raw/runtime trees. The local `libgfortran.so.5` is a runtime library, not a Fortran compiler. Obtain a matching aarch64 compiler/toolchain and required BLAS/LAPACK build interfaces before starting the dependency superbuild; do not treat the current build flag as the only missing item. Network acquisition and compilation were not performed by this audit.

The following is a concrete **proposed**, unexecuted build sequence after those host prerequisites are present. Separate prefixes preserve the audited current runtime. Dependency source archives/commits and resulting libraries need recorded hashes; upstream's URL downloads alone do not create an immutable local provenance manifest.

```sh
task_root="$PWD"
cmake -S "$task_root/data/raw/mechanics/opensim-core/dependencies" -B "$task_root/data/runtime/opensim-moco/dependencies-build" -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$task_root/data/runtime/opensim-moco/dependencies" -DOPENSIM_WITH_CASADI=ON -DSUPERBUILD_simbody=OFF -DSUPERBUILD_catch2=OFF -DSUPERBUILD_ezc3d=OFF
cmake --build "$task_root/data/runtime/opensim-moco/dependencies-build" --parallel 4
PKG_CONFIG_PATH="$task_root/data/runtime/opensim-moco/dependencies/ipopt/lib/pkgconfig" pkg-config --modversion ipopt
PKG_CONFIG_PATH="$task_root/data/runtime/opensim-moco/dependencies/ipopt/lib/pkgconfig" cmake -S "$task_root/data/raw/mechanics/opensim-core" -B "$task_root/data/runtime/opensim-moco/engine-build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$task_root/data/runtime/opensim-moco/install" -DSIMBODY_HOME="$task_root/data/runtime/opensim/install/simbody" -DOPENSIM_DEPENDENCIES_DIR="$task_root/data/runtime/opensim-moco/dependencies" -DOPENSIM_WITH_CASADI=ON -DBUILD_API_ONLY=ON -DBUILD_TESTING=OFF -DBUILD_PYTHON_WRAPPING=OFF -DBUILD_JAVA_WRAPPING=OFF -DOPENSIM_COPY_DEPENDENCIES=OFF -DOPENSIM_C3D_PARSER=None
cmake --build "$task_root/data/runtime/opensim-moco/engine-build" --target install --parallel 4
```

The dependency superbuild's Make-based METIS/MUMPS/IPOPT steps may need environment flags for the already-local BLAS/LAPACK sysroot and newly supplied Fortran toolchain; confirm configuration output instead of assuming discovery. Locate generated `casadi-config.cmake` and `ipopt.pc` if a platform uses `lib64`, and point `casadi_DIR`/`PKG_CONFIG_PATH` accordingly. With dependency copying disabled, extend the runner's runtime search path to new OpenSim, CasADi, IPOPT, MUMPS/METIS and existing Simbody/BLAS/LAPACK libraries. Build a tiny native IPOPT/Moco smoke solve before the 3D optimization, then compile the example/runner explicitly against the new prefix; API-only install still does not promise a command-line application. This sequence has not been run or validated on this host.

### Raw identities and existing registration provenance

Read-only SHA-256 values for the primary 3D example directory:

| File | SHA-256 |
|---|---|
| `subject_walk_scaled.osim` | `806ad1b248c6b2bc48ab43b58865915d1664e144b4ff6846457d80adacc23ce5` |
| `subject_walk_scaled_ContactForceSet.xml` | `659d571d9f731c9e119d93371edcd46c7fc265796ef0681670d46bd875b1b54d` |
| `subject_walk_scaled_ContactGeometrySet.xml` | `b64f1f1c7ebc9da5008c14896dc877f16b958505d8c4adbb29f86bdd3f1be47f` |
| `subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml` | `ff09a63fa9c2036a3692214889c87bee28bea6f51a503233692b288c5b28aba0` |
| `subject_walk_scaled_FunctionBasedPathSet.xml` | `932ac403e4d6d0d10d9a9da551e43d0ded94dd2b99c370f9b25053ee9232f8a5` |
| `coordinates.sto` | `981e94b5d7bb6fcae403a381fee8c3011533d8c683b17460096b400d20044cdb` |
| `grf_walk.mot` | `372d59855f80d7de81303cb7d3a5c5878ef8e0cb82a53a27cfce44e93357e42c` |
| `grf_walk.xml` | `3c5313d5680c52fd130a0ed86620472fe3c576a5a9599ed08a53652b84dd82e6` |
| `electromyography.sto` | `0746a08b1107f816ef35183d7639ca46dc863b60e178404d822ddc42fbb58061` |
| `example3DWalking.cpp` | `a6512ebe7ca4e54fb567f24c814a72803c6bb16c6ec72e2135b0ded47def283c` |

The retained opensim-core checkout is clean at the build-manifest commit. Current canonical registration instead points to original Rajagopal2016 raw hash `3f5c5f23e486073f2ad2aa4a4967ffe2fcdd582b1e355512bc54f70c36376bf4`, matching the raw file and its derived source metadata. That metadata records opensim-models source revision `d9b05d470b1a481c222372c85b75772faf8f7792`. The opensim-models source folder has no `.git`; `git -C ... rev-parse` would resolve the enclosing IHM repository and must not be used as its source revision.

All three canonical mechanics `source_files` hashes matched current file bytes: canonical anatomy `ee4494e43893d81d9798071d9228e201be62f009ee3e1374e522d1670250f64c`, corrected baseline mechanics `39b30a273ceee28ca1d1ec1a0766daf516a4f3b87a2f7e17118a2ef4e126c069`, and derived Rajagopal metadata `629744e8cb54e26f6130cbcb9bfaf6a193a014ad2968432c1cbe407aa3595208`. This verifies current registration identity, not registration accuracy or compatibility with the different 3D walking subject.

### Forward runner work

No retained executable currently exposes a walking launch command. The following is a proposed implementation contract, not a command claimed to run today.

1. Add a separate native C++ forward runner linked through the existing `find_package(OpenSim REQUIRED)` configuration. Follow `scripts/build_native_opensim.py`'s adapter build conventions. Keep immutable model inputs and fresh output directories.
2. First implement `ForwardTool` execution of a copied gait2392 setup with corrected relative file paths. Preserve original measured-GRF/correction settings for reference comparison. This is the shortest defensible forward regression using already retained states/controls.
3. Add programmatic assembly of the primary 3D contact model following the preparation part of `example3DWalking.cpp`; export the exact assembled model and coordinate/constraint/control list. Initialize q, u and auxiliary muscle states consistently from matched references. Do not derive a full muscle state from a q-only file without explicitly solving/identifying it.
4. Use native `Manager::initialize(state)` and `Manager::integrate(t_final)` with explicit accuracy and step limits. A 50–100 ms free response verifies integration/contact; it is not sustained gait. For a gait interval, require controls fitted to this exact model, then independently integrate them. The retained `OpenSim/Moco/MocoUtilities.cpp::simulateTrajectoryWithTimeStepping` demonstrates prescribing **controls**, initializing the full state, and then using Manager; prescribing controls does not prescribe coordinates.
5. Persist SI state trajectories, controls/excitations/activations, body transforms/velocities, muscle path/fiber/tendon quantities, GRF/COP/contact force and moment, constraint residuals, actuator/residual work, energies, contact events and all termination reasons. Copied source references must have separate filenames from computed outputs.
6. Add the tissue coupling contract above and verify reaction feedback before broadening domains.

Current read-only reproducibility commands (run from repository root):

```sh
file data/runtime/opensim/native_opensim_export
rg -n 'BUILD_API_ONLY|BUILD_PYTHON_WRAPPING|BUILD_JAVA_WRAPPING|OPENSIM_WITH_CASADI' data/runtime/opensim/opensim-build/CMakeCache.txt
nm -DC data/runtime/opensim/install/opensim/lib/libosimExampleComponents.so | rg 'ToyReflexController::computeControls'
rg -n 'Manager|integrate|realizeDynamics|equilibrateMuscles' scripts/native_opensim_export.cpp
rg -n 'CorrectionController|controls_file|states_file|external_loads_file|initial_time|final_time' data/raw/anatomy/opensim-models/source/Pipelines/Gait2392_Simbody/subject01_Setup_Forward.xml
rg -n 'ModOp|ContactGeometrySet|ContactForceSet|location\[1\]|set_initial_time|set_final_time' data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/example3DWalking.cpp
```

Observed results: aarch64 native executable; API-only/wrapping-off/CasADi-off flags; linked toy-reflex symbol; exporter only equilibrate/realizeDynamics matches (no Manager/integrate); gait2392 settings and 3D contact assembly as described above. The shell offers `python3`, not `python`.

Existing adapter environment uses library paths under `data/runtime/opensim/install/{opensim,simbody}/lib` plus the local `sysroot/usr/lib/aarch64-linux-gnu/{lapack,blas}` and parent directory. A new runner must preserve this environment and explicitly place a selected wrap variant first when used. Record executable/source/model/controller/library hashes, source commits, solver options, output units, input-table interpolation policy, loaded model state count, and elapsed simulation versus wall time.

## Required falsification and validation runs

These are proposed acceptance experiments, not completed results. Set numerical tolerances from scale/convergence studies before claiming validation.

| Experiment | Evidence it should establish |
|---|---|
| Gait2392 original forward rerun | Compare all state labels and time-aligned q/u, forces and body transforms with retained forward reference; report current-version differences rather than silently overwriting the reference. |
| Disable CorrectionController; then disable/perturb measured GRF | Exposes dependence on correction/residual supports and prescribed external loads. A fall or drift is useful causal evidence; do not conceal it by prescribing q. |
| Same controls, initial pelvis speed +/-0.05 m/s or knee angle +/-0.01 rad | State trajectories must respond; identical geometry would expose replay/prescription. Check constraints and muscle equilibration policy. |
| Contact stiffness 0.5x/2x, friction perturbation, small floor-height change | Simulated GRF/penetration/slip must change and forces must follow contact geometry. Record contact work and distinguish smooth-force leakage from true stance. |
| Horizontal pelvis force pulse, e.g. 20 N for 50 ms | Record 1 N s imposed impulse, COM momentum change, contact impulse and controller response; reduce pulse if outside validated domain. It is an engineering perturbation, not a clinical protocol. |
| Reflex gain zero versus baseline, increased delay, sensory channel ablation | Demonstrates feedback's causal contribution and delay sensitivity, with bounded controls. Feedforward-only comparison is essential. |
| Reduce maximum integrator/communication step by 2 and 4; tighten accuracy | Quantify convergence in q/u, GRF peaks/impulses, interface work, constraint residuals and event times. Do not infer accuracy from a finite trajectory alone. |
| Rigid rotation/translation of registration fixture | No spurious strain or changed anatomical attachment lengths under a common rigid motion. Verify frame inversion and moment origin. |
| Tissue density/material perturbation; feedback off/on | Changes to genuine tissue dynamics should alter inertial/deformation response; turning on reaction feedback should change free skeletal motion. |
| Energy/momentum and mass allocation audit | Close internal force/moment exchange, account for contact and active muscle work, avoid duplicate segment/tissue mass and duplicate contact. |
| Mesh/contact and constitutive refinement per resolved domain | Positive element Jacobians; bounded penetration; convergence of traction/displacement/stress; physically appropriate material tests. Numerical convergence is distinct from subject/experimental validation. |

For whole-body claims also require sustained multi-step gait, lateral stability, speed transitions, multiple initial conditions and held-out perturbations, plus an explicit anatomy/material/interaction coverage table. Report which organs are merely carried by the skeleton, which have affine deformation, which resolve local dynamic strain, and which exchange forces through contact/sliding attachments. Those categories must remain visible even as the system expands.
