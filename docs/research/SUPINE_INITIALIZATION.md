# Supine initialization: diagnosis and explicit static residual construction

The current 92-muscle v2 source pose is not a supported equilibrium. In retained run `supine-support-5ma720yd`, the initial plane X coordinate is -0.403499 m and touches only the torso's 0.25766 m COM/inertia-derived sphere. That sphere initially generates 21.276 N, about 2.79% of body weight. The other posterior proxies float 0.174–0.373 m above the plane; the pelvis gap is 0.23815 m. A common rigid translation cannot change these relative gaps. Raising the plane until it touches the pelvis would penetrate the torso proxy by 0.23815 m, far beyond the unchanged 0.05 m diagnostic stop bound. Arbitrary plane translation or preload is therefore not an acceptable initialization repair.

The native source retains zero joint angles apart from pelvis height (0.93 m), and turns gravity to source -X. `equilibrateMuscles` solves internal musculotendon state equilibrium at that pose; it does not solve generalized static body equilibrium under gravity and contact. Source80 Millard muscles retain 1% default activation. The added Thelen arm muscles use their native 5% default activation (that property is implicit in their XML). Initial biceps-long fiber length is 1.281 times optimal length, with tendon force about 14.5% of maximum isometric force. Those values diagnose a loaded retained pose; they do not establish passive strain, pathology, or a reason to disable muscle force. Passive and active tendon-force components are not separately exposed in the retained observation.

`build_supine_initial_state.py` first produces a source-only, identity-bound diagnosis and static-pose seed. The seed excludes dependent patellar coordinates, retains independent source bounds, default excitation, gravity and exact plane placement, and explicitly declares `accepted_initial_state=false`. Its fixture checks confirm torso-only support, pelvis clearance, rigid-translation invariance, dependent-coordinate exclusion and rejection of an unverified seed as equilibrium. The first such receipt is `data/derived/supine-initialization-4fcn25e9`.

## Native static-pose evaluation contract

The proposed repair is a bounded static residual solve of the existing physical model, followed by unchanged forward acceptance. A static candidate is evaluated on a **copy** of the continuing SimTK State. The `evaluate_static_pose` command in `native_static_pose.h` permits only finite, source-bounded independent free coordinates and rejects duplicates, locked/prescribed/dependent coordinates, external loads and non-supine environments. It assembles dependent coordinates, assigns zero generalized speeds only to the isolated static candidate, and equilibrates muscle states at unchanged controller excitation. It neither advances physical time nor changes work/metabolic accumulators, continuing state, or its immutable metabolic reference.

The response includes source-order q, u and udot; per-coordinate values, speeds and accelerations; position/velocity/acceleration constraint errors and multipliers; tree and constrained zero-acceleration residual mobility forces; COM acceleration; contact forces and penetration; and candidate state-variable values. The tree residual ignores constraint force; it can be nonzero even for a valid constrained equilibrium. The constrained residual uses the native dynamic multipliers and is a force/torque diagnostic in Simbody mobility order. These mixed units are reported explicitly. A single near-zero COM acceleration is insufficient: the full generalized acceleration must settle.

The Python solver uses bounded nonlinear least squares with at most **200 actual candidate calls**, including finite-difference calls, and a **60-second wall deadline** that kills and reaps its owned native process. Source bounds are enforced both before dispatch and after native assembly. Residual conditioning uses g for translations and g/0.3 m for angular accelerations; this is an optimizer scale only and does not alter native forces. No gravity cancellation, contact preload, hidden controller change, muscle disabling or coefficient tuning is permitted. Candidate state observations and every evaluation norm/cost are retained. The continuing state is observed and compared before restoring its checkpoint.

A statically converged candidate requires all native udot magnitudes ≤1e-4 in their source mixed units, constraint position error ≤1e-5 and posterior penetration ≤0.01 m. Optimizer termination alone is insufficient. Even then it remains an unaccepted initial state until an explicit new-run initialization establishes an appropriate reference and the unchanged sustained forward criteria in `SUPINE_SUPPORT_ACCEPTANCE.md` pass. This first implementation does not import or publish candidate state observations as a live supported body.

## Limits and remaining evidence

The source joint ranges include broad engineering values, and the COM/inertia spheres are not posterior anatomy. A static solution could still have an inappropriate folded pose; its supine morphology needs separate assessment. The solver therefore cannot by itself validate these proxies as a mattress or validate canonical tissue contact. Failure to find a solution within the finite budget is an unresolved initialization result, not permission to tune the convergence thresholds or extend the run automatically.

No native static-pose run is claimed in this initial implementation. Build and native validation require the shared heavy-resource queue; source-only fixture checks and Python syntax verification run without a native process.

## First authorized native static solve — 2026-09-05

The header and minimal dispatch compiled successfully as `build-z5d66e33` in **22.35 s**, with **1,172,432 KiB** peak compiler child RSS. One authorized static solve at `data/derived/supine-static-solve-zppmkd4p` stopped at its hard **200 actual evaluations** after **7.495 s**. No physical time advanced. An observation before checkpoint restore confirmed that the continuing state remained unchanged; restore matched too. Dependent-coordinate and source-range rejection checks passed. The process was reaped and the shared resource slot released.

The initial scaled residual norm was **168.0408**, and the best was **167.7376**: no static solution was established. Maximum generalized acceleration fell only from **2964.35** to **2958.04 rad/s²**, dominated by the elbows. The best candidate had zero kinetic energy by construction but still accelerated its COM at **9.6764 m/s²** toward the bed and had only **10.369 N** of support. Zero candidate velocity therefore provides no support evidence. The solver's best candidate is explicitly **not accepted or installed**. The retained report's `stopped` status denotes the evaluation budget, not physical divergence.

A source-only follow-up identifies an additional specific source of initial pose imbalance. The retained passive elbow law is `6.09*exp(-6.94*(q-0.30))-11.03*exp(11.33*(q-2.40))-1.0*qdot`. At q=0 and zero speed it supplies **48.845 N·m** of flexion. Its isolated zero-torque angle is **1.56979 rad**, far from the retained straight-elbow default. This follows directly from the retained source expression and does not require a coefficient change. It does not imply that this angle is the full muscle/contact equilibrium. A passive-law-informed seed may merit a separately bounded investigation; the first 200-call solve barely moved elbow flexion from zero. No second static solve was run.

Evidence SHA-256:

- `report.json`: `b39a163d736854b83c061f4f98a6a2fa20280510f7faf89204cf4d996e41f305`
- `evaluations.json`: `c5af52e49e4cb825d1577998bfec116142612947d18df66dba43d5d8190329ee`
- `best_candidate.json`: `6d601de3e292430d779549d8c23504bea16b446b888bc35504434e358cdf1266`
- `seed.json`: `4bf25da44252a0ef27eaec8d883bd705e1b53fe9fbc769a24d52b13a86a455e0`
- `native/execution.json`: `ddc5587611b322e950e75cc1306ab29d9cbe184caa83050768f2c44bb57c67cb`
- `native/inputs/subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml`: `ff09a63fa9c2036a3692214889c87bee28bea6f51a503233692b288c5b28aba0`

## Source passive-neutral seed construction

`--passive-neutral` derives isolated coordinate-force zeros from the retained `ExpressionBasedCoordinateForceSet`, at qdot=0 and within the model's declared coordinate bounds. A restricted arithmetic/`exp` AST evaluator avoids general expression execution. Each active coordinate law is sampled at 129 bounded points to find brackets; each sign-change bracket gets 80 bisection steps, and the nearest root to the original coordinate is selected. Identically zero laws retain the original coordinate. Multiple laws on one coordinate are summed. A missing bounded root or residual above 1e-6 fails closed. This procedure changes only candidate q; model coefficients, controller excitation, gravity and plane remain identical.

The held source produces bilateral elbow flexion 1.569788 rad, knee flexion 0.490078 rad, arm flexion 0.160698 rad, arm adduction -0.625613 rad, and hip flexion -0.060267 rad. Pure ankle damping has no preferred static angle and preserves zero. These are **isolated passive-coordinate-force zeros**, not muscle-force or whole-body equilibrium. The native residual solver must evaluate all unchanged muscle/gravity/contact forces at these poses. Initial source receipt: `data/derived/supine-initialization-n69xgzc9`; fixtures verify all resulting roots, unchanged controls/plane and rejection of unsafe expression syntax.

The proposed minimal future clock-zero initialization interface is an optional identity-bound startup q seed, applied and assembled after constructing the original contact plane but before initial muscle equilibration and the first observation. It must retain the original fixed bed geometry, zero starting speed, source coefficients and controls. Such a seed would be labeled `unsettled_seed`, with its startup metabolic observation explicitly **not verified as a resting reference**. This is a proposal, not an implemented state-import or accepted-baseline feature. Sustained forward verification remains necessary before claiming a supported reference state.

## Authorized passive-neutral static experiment — 2026-09-05

One run, `data/derived/supine-static-solve-qblxfx5z`, used the same native model/build/environment/controls as the default-seed experiment and changed only the documented starting q. It stopped at **200 evaluations / 2.472 wall seconds**, with zero physical-time advance. Continuing-state observation and checkpoint restore checks passed. No startup import or reference promotion occurred.

The raw isolated passive-neutral seed was **worse as a support pose**: maximum proxy penetration **90.881 mm**, maximum generalized acceleration **13,274.68 rad/s²**, and scaled residual norm **805.576**. Isolated coordinate-force neutrality clearly does not imply feasible muscle/contact geometry.

The optimizer reduced its aggregate residual to **26.9125**, compared with **167.7376** for the previous default-seed search. However, this numerical improvement is **not an improvement in supported equilibrium**. Its best candidate lifted the pelvis about **0.205 m** along +X, away from the fixed bed, leaving support of only **2.14e-6 N** and COM acceleration **9.81 m/s²** toward the bed. Maximum generalized acceleration remained **609.50 rad/s²**, dominated by ankles. Zero penetration reflected lost contact, not proper load support. The low residual therefore selects a largely unsupported pose while reducing larger limb accelerations.

This result does not satisfy the condition for implementing a clock-zero startup seed as a useful supported initialization. The next static formulation needs explicit support-force and feasible-contact constraints, with full generalized force balance, so limb residual reduction cannot compensate for losing bed support. No coefficient tuning, tolerance relaxation, reference reset, or additional native run was performed. All source-law seeds and native candidates remain marked unaccepted.

Evidence SHA-256:

- `report.json`: `d72055d8db96c3cbf4df2e84f0d65febd2cb035dc9a107db6cf2d05763fa40a9`
- `evaluations.json`: `114142d36f99eb6c4e96a56250ec19c864837c7dbf709f954588cc1bc8be45e8`
- `best_candidate.json`: `c46bbc0bcff41ec8a38e9d4c85d30080cf349d9fc5faa76097b297f574705ed9`
- `seed.json`: `4f368cee657ec5e016e7bbe8838a23b281dd16ee6709ff254a8125569a7b2875`
- `native/execution.json`: `cdaeeb2abfa4e142e202ccaa08585b9c3f7f305fc4a477cc785c025225de68ae`


## Measured-bed constrained statics — 2026-09-05

The measured MM bed/skin series law and posterior skin quadrature replace the inertia spheres in these opt-in experiments. Native MM and HM two-millisecond implementation fixtures passed (`native-surface-foundation-0hzzust1` and `native-surface-foundation-ch3sdupy`); their fresh default poses still had zero support. They verify the contact law, separate energies and indentation identities, not a settled body.

The rigid source-only seed `supported-rigid-seed-ljyewh9y` balances weight and both normal-support moments for MM within the declared material domain. The first full native constrained solve (`constrained-supine-k970xisy`) stopped at 200 evaluations in 41.008 s. It held plane translation/heading gauges without locking physical joints and constrained normal force, pitch and roll residuals. Maximum generalized acceleration remained 1525.55 rad/s². Uniform weight-times-height torque normalization understated distal acceleration, so its small objective was not equilibrium evidence.

The replacement objective is `udotᵀ M udot / (mass |gravity|²)`, cross-checked with `−rᵀ udot` and `r + M udot = 0` under the emitted dynamic constraint multiplier convention. Negative or nonfinite metrics fail; no cost floor is used. Native fixture `native-static-metric-tfbdb6ua` passed with identity norm 3.33e−13 and unchanged continuing state, using build `build-64ooel0c`.

One authorized resumed MM solve (`constrained-supine-672_i109`) exhausted 200 native calls in 53.763 s; no physical time advanced. It reset only toes to their isolated source passive-law neutral, retaining ankle q because pure ankle damping has no preferred static angle. Its best support-feasible dynamic metric decreased from 3.45733 to 3.37995, and maximum acceleration from 1526.61 to 1486.95 rad/s². Toe accelerations remained 1466.54/1486.95 rad/s², ankle accelerations −992.99/−989.60 rad/s², and elbows 447.55/441.06 rad/s². Skin compression stayed at 1.35524 mm and bed deflection at 64.0424 mm.

The support-feasible candidate's raw generalized residuals were pelvis translation X 0.0186719 N, pitch 0.0473951 N·m, and roll 0.000393530 N·m. Held translation gauges Y/Z had residuals 2.84e−14/7.11e−15 N and heading −3.71e−5 N·m. These nearly balanced external loads coexist with large internal accelerations. The lowest unrestricted objective 2.93574 violated the force constraint (normalized residual −0.000487793) and was not selected. Neither candidate meets the unchanged all-acceleration threshold of 1e−4 or has forward support acceptance. At the budget exception, the solver closed/reaped the native process; the separate identity fixture supplies continuing-state immutability evidence, not an unperformed post-budget observation.

No startup pose or metabolic reference was promoted. A future accepted startup artifact must be loaded before first metabolic sampling in a fresh process, with zero clocks/work and a separate immutable neutral `registration_reference` to preserve canonical mesh embedding. Supported-pose transforms must not become the source registration reference. Cervical model or contact-owner changes invalidate these 22-body candidates.

Receipt hashes:

- `data/derived/constrained-supine-672_i109/best_supported_candidate.json`: `43ea1233176f0099f3b3416f0f8766df1e2e07a44aa313d5b806ba921057f52d`
- `data/derived/constrained-supine-672_i109/initial_supported_seed_native.json`: `587c234662f5762226af0346c1cf028fa82fbee368959696c37b064cfdddc1a1`
- `data/derived/constrained-supine-672_i109/report.json`: `1e768c791ab25dbd62cf19f94f6d3906b57a3ff7d7e2eb37f6643aa9e02735dc`
- `data/derived/native-static-metric-tfbdb6ua/report.json`: `31897538b3d6166e7c598ab8bcf1c03b1e45bfc8e042b731f0f902370982cfa9`


## Source-bound evaluation recovery

The constrained harness now writes every complete native response, requested coordinate vector, objective and support/gauge residual to `candidates.jsonl`, with per-record SHA-256. `cache_identity.json` binds native input hashes, build/dependency hashes, protocol/helper source, material, mass, coordinate order/bounds, held gauges and geometry manifest. `--resume-cache DIRECTORY` fails on identity mismatch, truncated records or conflicting responses; keys use exact floating-point representations, with no interpolation. Cached responses do not consume the 200-new-call allowance. All native work, cache loading and optimization remain inside the existing hard 60-second alarm.

Optimizer callback iterates are retained in `optimizer_iterates.jsonl` and `last_optimizer_iterate.json`; all finite-difference samples remain recoverable from the evaluation journal. Restart uses a chosen recorded q and reconstructs matching finite differences from exact cache hits. It does **not** claim to restore SciPy's internal quasi-Newton/trust-region state. Pre-journal runs cannot recover evaluations that were never retained. Source-only fixtures verify full response recovery, exact-key separation, identity rejection and truncated-record rejection. No native resume has yet tested this increment.
