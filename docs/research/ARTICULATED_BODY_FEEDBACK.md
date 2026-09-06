# Continuing articulated body and garment feedback

The implementation uses one retained OpenSim model as the body inertia owner. It retains native joints, source muscle/tendon/activation states and complete in-process `SimTK::State` checkpoints. Canonical anatomy supplies attached geometry and explicit registration priors; this does not establish all-organ finite-element deformation.

## Interface and current acceptance

`ArticulatedBodyPlant(root, output, environment='supine', target_mass_kg=..., augmented_registration=None, enable_garments=False)` provides `snapshot`, `advance(dt_s, forces, actuation)`, `checkpoint`, `restore`, `release`, and `close`. Force records are canonical world `{id, point_m, force_n}`; muscle commands are excitations, not an additional activation state. A step is at most 20 ms.

The retained source model has 22 inertial bodies and 80 muscles. Passing the separately retained Arm26 registration selects its augmented 92-muscle model; exact model, source and catalog hashes are checked and copied. Native coverage must match the selected catalog. The corrected `whole_body_arm26_v2/registration.json` has passed the short native acceptance below. This does not validate its moment arms across the full three-dimensional motion range; XML coverage alone was demonstrably insufficient.

The initial source80 native smoke preceded the metabolism/Arm26 extension: `data/derived/native-stream-smoke-n6e0pvqi/supine/smoke.json`. A 2 ms advance gave a total linear momentum residual below 3e-13 N and exact transform replay after checkpoint restoration. Initial support was only about 20 N against about 756 N weight. **This is an initial falling/contact transient, not settled supine equilibrium.** The extension compiled in `data/runtime/mechanical-stream/build-09xtbmes` in 5.62 seconds with 1,168,724 KiB compiler peak RSS under the 4 GiB address-space ceiling. Its first 92-muscle initialization failed before dynamics because the donor pre40500 GeometryPath property name was not migrated into the current target XML schema. The failure is preserved in `data/derived/articulated-acceptance-h2tqwejy`; no new native acceptance is claimed by that failed initialization. Old executable identity cannot silently satisfy new source checks.

The corrected v2 native acceptance is retained at `data/derived/articulated-acceptance-noo34jds`. All 92 muscles and 2408 canonical entities were present. Four 2 ms branches completed in 1.32 seconds, with native peak RSS 101,556 KiB under the 4 GiB ceiling. Increasing right brachialis excitation to 0.8 raised its activation from 0.05 to 0.232323 and native total muscle metabolic power from 135.413681 W to 193.230602 W. Metabolic interval energies were 0.266094 J and 0.361585 J. A separate canonical hand-force branch changed the pose; all branch total-momentum residuals were below 3.44e-13 N. Restored baseline body poses and cumulative metabolic energy reproduced exactly. Garments were disabled, and approximately 20.5 N initial contact remains far below full body weight: this is not settled support.

## Shared physical coordinates and mass

One proper rigid least-squares fit maps source ground into canonical ground. The 22 correspondence pairs are native segment COMs and canonical named-bone bounding-box centers: approximate engineering anchors, not homologous measured landmarks. Fit RMS residual is 0.0623561 m and maximum 0.132032 m. No per-body world-frame translation conceals these errors. The same source bed plane, native joint points and forces have one canonical world location.

Canonical reference mesh points have fixed segment-local embeddings through the shared registration. Named bones use retained grouping; other tissues use an explicitly inferred nearest bone-envelope assignment. Initial geometry is preserved, but anatomical joint surfaces and continuous soft interfaces may disagree with source kinematics. Retained-state tests verify identity initialization, shared rigid motion, and wrench virtual-power invariance for all 2408 entities.

Native segment masses and inertias are uniformly scaled together to an explicitly selected generic-body mass; source proportions are retained as a prior. Original native total is 85.26985 kg. Canonical mechanics total is 77.11070 kg. Root orchestration may select the native physiological saved-state weight, which includes baseline GI contents. Canonical attached tissues add no second body inertia. Cloth has its own explicit finite mass. Metabolic analysis mass is not added to mechanics.

## Environment and energy

Supine gravity acts along native negative X. The single unilateral plane contacts 22 posterior spheres centered at retained source segment COMs, with radius inferred from an equivalent uniform ellipsoid's moments of inertia. These are contact proxies, not anatomical skin surfaces. Source foot `SmoothSphereHalfSpaceForce` parameters are transferred as an engineering prior; there is no calibrated mattress claim. Force and moment outputs are complete world-frame resultants about body origins. Their anatomical distribution cannot be recovered from a proxy contact area.

`body_environment.canonical_wrenches` retains force, moment and current body-origin point; this must not be treated as a point force on deforming skin or lungs. `body_environment.plane` provides the common canonical plane. Native contact action/reaction and whole-body momentum residuals remain observable.

External point-force work and positive active-fiber work use endpoint trapezoidal power quadrature. The frozen source-native Umberger/Uchida helper separately reports total muscle metabolic power, muscle heat, signed active-fiber work and analysis mass. Its cumulative metabolic energy enters checkpoints. Generic fiber composition and analysis-mass priors remain explicit. Positive mechanical work is not metabolic demand. Total muscle demand cannot simply be added atop overlapping BioGears basal demand; orchestration owns that budget.

## Whole garments

`GarmentFeedback.from_root` materializes the retained full shirt (2038 nodes) and shorts (2527 nodes), preserving source coordinates and shared-index topology. Explicit current engineering parameters are 0.18 kg/m² areal density, 12 N/m edge stiffness, 12 N/m elastic support stiffness, and Coulomb coefficients 0.4/0.3. These are not measured whole-fabric/skin calibration. Sixteen shirt and sixteen waistband source-triangle supports use barycentric elastic forces with explicit initial rest distances.

The retained source skin mesh is bound once per vertex to the nearest reference bone envelope. Each cloth interval subcycles below both edge and added-support stiffness bounds. A moving source-triangle boundary interpolates predicted native endpoint positions. Nearest face-interior unilateral impulses and Coulomb friction return equal/opposite body impulses. Exact resultant force and moment are represented by four native point forces per segment. Native and cloth trials restart from the same checkpoint until the maximum boundary-position mismatch is within the declared tolerance. Failure restores both owners and their clocks.

The report separates friction/normal dissipation, cloth and support-spring integration energy defects, native garment work, and interface quadrature mismatch. Tiny tests verify force/moment reduction, finite-body reaction, total momentum, common clocks and failed-convergence rollback. Whole-body garment drape has not yet been accepted. Garments remain opt-in pending the bounded source-native run. Discrete face contact does not certify continuous collision, self-contact, persistent friction, edge/vertex collision, global exterior/cavity interpretation, intersegment skin continuity or full containment. There is no genital tucking claim.

## Resource and provenance discipline

A fresh thin build retains its C++ source, frozen metabolic header, command and library hashes. New engine processes and compilation launch through `prlimit --as=4294967296 -- nice -n 10`, with single BLAS/OMP threads. Missing limit support fails clearly. This limits new IHM-owned subprocesses, not unrelated IBM work. The Python caller still needs its own process resource discipline. No whole-body native run is authorized merely by importing a verifier.

Pure small checks:

- `scripts/verify_articulated_registration.py`: retained state, no native launch.
- `scripts/verify_moving_surface_contact.py`: moving-boundary contact parity, impulse/work and bounded candidate allocation.
- `scripts/verify_garment_feedback.py`: finite-mass patch/plant fixture, conservation and transactional failure.

`verify_articulated_native.py --run-native` is a separately coordinated native acceptance, with no default native launch. Its optional `--garments` step is 50 microseconds and remains unrun; enabling that flag is not evidence of garment acceptance.

The implementation follows retained OpenSim/Simbody native APIs and source contact laws. Relevant primary references include [Simbody force API](https://simbody.github.io/3.8.0/classSimTK_1_1Force_1_1DiscreteForces.html), [OpenSim simulation concepts](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089017/SimTK%2BSimulation%2BConcepts), and the separately retained muscle-metabolism source receipts in `NATIVE_MUSCLE_METABOLISM.md`.
