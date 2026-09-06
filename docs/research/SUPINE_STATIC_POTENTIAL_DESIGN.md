# Fixed-activation static potential: derivation and bounded solver proposal

No native run, integration, controller change or new static solve was performed for this design. The target remains the separately bound 98-muscle, seam-omitted reference. All physical equilibrium thresholds remain unchanged.

## What is conservative at zero speed

Gravity has potential −Σ m g·x_COM. Each zero-speed expression force Q_i(q_i,0) has a scalar primitive −∫Q_i dq_i. The 20 retained expressions are separable; their velocity terms vanish in this diagnostic. The native expression component does not provide a corresponding potential override, so simply querying native total potential energy would omit these primitives.

The existing skin/mattress energy is already source-defined. The series law solves equal pressure between bounded skin and mattress compression, and stores the sum of their energies. On an interior valid branch, differentiating with respect to total approach gives pressure; the two internal deformation derivatives cancel at equal pressure. The unilateral zero-contact branch is continuous, and the same source deformation limits must reject invalid trials. No energy is assigned to the omitted seam cell. Friction/damping at zero velocity contributes zero force here; this does not make a moving contact experiment conservative.

## Muscle derivation and its limits

For fixed activation a, zero fiber velocity and a smooth path length L(q), define fiber force

F_f(l_f,a) = F_iso [a f_AL(l_f/l_opt) + f_PE(l_f/l_opt)].

Let the constant pennation height be h and the fiber projection s(l_f)=sqrt(l_f²−h²), so tendon length l_t=L−s(l_f), and cos(alpha)=s/l_f. A candidate *effective* potential is

U_eff(L,l_f;a) = U_t(L−s(l_f)) + ∫ F_f(l_f,a) dl_f.

Its fiber derivative is F_f−T/cos(alpha). Thus an interior stationary fiber satisfies the native static force balance T=F_f cos(alpha). Along a differentiable stationary branch l_f*(L), the envelope derivative is dU_eff(L,l_f*(L);a)/dL=T. With source-consistent moment arm r_i=−∂L/∂q_i, the force is Q_i=−∂U_eff/∂q_i=T r_i.

This derivation requires a stable, locally single-valued internal branch, finite pennation geometry and no unaccounted minimum-fiber constraint. Positive internal second derivative is necessary for eliminating the fiber as a local energy minimum. Descending active force-length behavior can introduce nonconvexity/multiple branches; blindly taking a global fiber-energy minimum could switch away from the branch selected by native initialization. Native clamped/minimum-length warning branches need separate inequality/reaction treatment and must not silently pass the interior derivation.

The active primitive is a mathematical fixed-activation work potential, **not physiological stored energy**, chemical energy or metabolic credit. If activation changes, its parameter derivative supplies an additional term; in moving muscles, force-velocity and damping invalidate this static potential reduction. It must never be substituted into the runtime signed muscle energy ledger or used to reset M0.

Retained Thelen source initializes a bounded force-equilibrium solve and explicitly distinguishes minimum-fiber warnings. Its active curve is Gaussian, passive curve exponential above normalized length one, and tendon curve exponential/linear with zero tension below slack. The diagnostic analytically integrates those exact functions. Millard source likewise uses active + passive + damping fiber force and explicitly sets velocities to zero for static equilibrium. The same local derivation applies structurally, but the 80 actual Millard curves, lower-bound behavior and pennation have not yet passed a native full-body effective-energy gradient check.

All 80 fitted paths in the actual assembled model omit independent moment-arm functions. `FunctionBasedPath.cpp` consequently computes moment arms as the negative derivative of its length polynomial. The six lumbar geometry paths have retained native virtual-work validation. The 12 arm paths still require branch-aware wrapping/virtual-work coverage for this proposed whole-body merit function.

## Actual source-law evidence

`scripts/diagnose_fixed_activation_potential.py` reproduces the six zero-pennation lumbar Thelen branches at all 126 retained observations in `lumbar-muscle-native-lb45uirs/observations.csv`. It preserves activation 0.05 and exact model parameters. Maximum tendon-force discrepancy is 8.56965e−6 N, within each muscle's native initialization tolerance 1e−8 F_iso. The maximum effective-energy length-derivative error is 9.25511e−9 N with a central 1e−7 m difference. The minimum sampled internal stiffness is 72,861.66 N/m. These results validate these sampled branches, not every possible length or the other 92 muscles.

`verify_fixed_activation_potential.py` checks those per-muscle tolerances, work derivatives, positive sampled internal stiffness and all 20 zero-speed expression-force primitives. No runtime muscle evaluation was added. The retained native lumbar virtual-work test independently bounds moment-arm error around 1.07e−10 m.

## Concrete next solver and its prerequisite gate

The proposed solver is reduced-coordinate, source-bounded trust-region minimization of the complete *effective* potential, with a physical force gradient and a safeguarded BFGS/Newton model. Use a sensitivity-scaled coordinate metric, report its scales, and choose radius from actual versus predicted energy decrease. The effective energy replaces the poorly conditioned acceleration-squared merit; it does not replace the final acceleration test. The prototype should first remain a separate diagnostic, not a default adapter behavior.

Before any optimization, a bounded native fixture must explicitly emit/check gravity energy, skin/bed energy, the 20 joint primitives and all 98 fixed-activation muscle primitives, along with tendon force, path length, activation, fiber velocity, pennation and branch/clamp status. Check finite-difference energy gradients against actual generalized virtual work at the new reference and a few source-bounded perturbations. Small closed-loop work checks and repeated-length/branch checks must fail on inconsistent moment arms or history-dependent branch selection. Passive energy must be counted once; the active primitive and expression contributions cannot be assumed present in `calcPotentialEnergy()`. A scalar integral assembled from noisy finite differences of the whole-body force is not an adequate substitute.

True kinematic couplers are eliminated or enforced exactly. The coordinate-to-mobility work map must be explicit: pull back generalized force through the independent-coordinate velocity map, including both knee and dependent beta components, rather than dropping dependent mobility residuals. Verify the resulting gradient by virtual work. Held global translation/heading gauges remain symmetry choices whose actual residuals must still pass.

Force/pitch/roll balance should emerge from stationarity of the free root coordinates under gravity/contact energy. Adding those force balances as extra optimization equalities would produce mathematical multipliers that are not physical contact reactions; a constrained-energy stationary point with nonzero such multipliers is not equilibrium. Similarly, XML coordinate bounds are retained as trial-domain limits, but a bound KKT multiplier must not be counted as an anatomical stop torque. At any candidate, all actual accelerations, force/moment balances, gauges and native constraint errors must meet the original thresholds. A minimum resting against an unsupported coordinate bound fails that acceptance.

If the complete effective-gradient prerequisite fails or branch handling remains unresolved, retain a sensitivity-scaled support-tangent trust-region force-root method instead. Do not guess inertia weights from the emitted M·udot vector, which does not determine the full mass matrix. No physics coefficient, activation, ROM bound or equilibrium tolerance is adjusted by either numerical option.

Audit receipt SHA-256: `69e41d790e9e7e81fb80d8071cfa43bd84b5564aeb8b2e11039a4a97e9b4fdf3` (`data/derived/lumbar-supine-static-1g1q08u2/fixed_activation_potential_audit.json`).

## Prepared isolated native prerequisite probe

The source-only build manifest `data/derived/effective-potential-build-pqbqcago/manifest.json` copies the frozen engine source/headers and actual 98-muscle, seam-omitted input directory. The only engine additions are a new diagnostic header and command. It links the independently attested archived libraries; SDK headers and every copied/source input are hashed. It never changes `latest`, the old archive or physiological energy ledgers. Compilation is queued, with a 4GiB address-space limit, one thread/nice10 and a 60s process-group kill/reap cap. Python source/pullback fixtures pass; C++ compilation has not yet been attempted.

The new header evaluates separate gravity, joint, skin, mattress, native passive-muscle and active-effective primitives after the same copied-state muscle equilibration. Gravity/contact generalized forces are independently obtained from body spatial wrenches via the native system Jacobian transpose; passive expression forces are mapped into mobilities directly. It emits every muscle's native tendon force, length, moment arms, activation, fiber velocity, minimum-length margin, pennation and internal stiffness. The active Millard curve is integrated with bounded adaptive Simpson quadrature; the Thelen Gaussian uses its analytic primitive. The merit is explicitly labeled numerical only.

The native fixture is capped at 65 evaluations/45s. It plans one base, two source-bounded perturbations for each of 31 independent coordinates and one repeat (64 actual requests). It reconstructs a second-order directional gradient using the *actual assembled* coordinate-difference matrix, rejecting ill-conditioning. It independently checks muscle length/moment-arm derivatives, muscle energy/work, gravity, joint and combined skin/mattress work. It verifies native N=I instead of assuming coordinate rates equal mobilities, and explicitly pulls back both knee and dependent beta residuals. Ideal coupler reaction work must cancel; actuator forces outside the documented merit must be zero.

The preregistered diagnostic gradient gate is 1e−4 in each native generalized-force unit, with length/moment-arm error <=1e−6m, positive interior branch stiffness, fiber minimum margin >1e−10m and fiber velocity <=1e−5m/s across sampled states. These are numerical prerequisite tolerances, not relaxed physical equilibrium criteria. Actual errors and branch failures are retained regardless of gate outcome. A passed diagnostic alone does not activate an optimizer, certify global branch uniqueness, or promote an initial physiological reference. Further trust-region implementation remains conditional on actual results.


## Isolated compilation result

The first compile (`effective-potential-build-pqbqcago`) failed because OpenSim `PathActuator::computeMomentArm` requires a mutable Coordinate reference; its failed log is retained. The API-only correction obtains that reference through `model.updCoordinateSet().get(name)`, without changing q or the physical force model. One authorized retry (`effective-potential-build-3a9juno_`) compiled successfully under the existing4GiB/60s, nice10/single-thread limits. The compiler was reaped and the slot released; the native gradient fixture has not yet run.

Build manifest SHA-256: `8df1d4a894bf80ee52c1cc27a851d477b3ecfa41da6668e5166fc633c2b456da`; executable SHA-256: `f4d387bde698a4074fddf383d4df792bc798fd82f628e5a82bf37962415bb9a0`.


## Actual all98 gradient gate: failed, no optimizer

`native-effective-potential-_livswt7` completed64 copied-state evaluations in13.287s, with continuing state unchanged and owned process reaped. No physical time advanced. The actual assembled difference matrix condition was2.9387; native coordinate-rate/mobility mapping passed, ideal knee-coupler reaction work canceled to4.44e−16, and all sampled fiber branches met the declared interior/stiffness/velocity checks. Nevertheless the gradient gate **failed**: maximum effective-energy versus physical generalized-force error was10.17787. Component errors were gravity1.39862e−5, joint expressions1.30733e−8 and contact1.92377e−4 in their native force/torque units. A repeat changed actual q by up to1.43832e−6, so its energy difference−1.00884e−7J is not a same-q branch-equivalence certificate.

The largest defect is source-specific: retained `Thelen2003Muscle.cpp` lines420–428 passes `mli.fiberLength` (metres) into `calcfpefisoPE`, while that function at1468 onward expects normalized fiber length and activates only above1. All18Thelen fibers here are shorter than1m. Native passive fiber potential is therefore omitted despite nonzero passive fiber force. Using `getMusclePotentialEnergy()` uncritically in the proposed merit was wrong. The exact source-law passive primitive from the earlier analytical derivation can correct the *numerical merit*; no physiological library or energy ledger was modified.

The offline attribution retains every original response and adds the normalized-length primitive minus the actual metres-argument primitive. It reduces the full gradient error from10.17787 to0.01572835. The four remaining dominant muscle failures are bilateral Arm26 BIClong and BRA, all with `hybrid` path wraps. Corrected muscle energy versus tendon-force times measured path-length derivative agrees about1e−9 for those paths, but their measured length gradients disagree with native moment arms: BIClong errors0.00010719/0.00023356m, BRA0.00010915/0.00013442m. The corresponding worst torque discrepancy is0.01572835N·m. The other94muscles have corrected energy/work errors at most1.26690e−5. No wrapped muscle was excluded to make a balance pass.

The normalization fixture reproduces the erroneous zero energy for a0.14m fiber with0.1m optimal length, verifies the correct primitive derivative against passive force, and asserts that the remaining native gate still fails. Next work requires a targeted source/path-wrapping virtual-work and numerical precision audit, plus contact finite-difference refinement. The complete effective-potential solver is not enabled. Neither unchanged acceleration optimization nor a tolerance increase is justified by this failed prerequisite.

Evidence SHA-256:

- `data/derived/native-effective-potential-_livswt7/report.json`: `eed0fa36f75131ed8d71673a76453a23ebedac8d04bdfd4cc9c51fcd2bd2aa3e`
- `data/derived/native-effective-potential-_livswt7/observations.jsonl`: `b179b815e1fa9ca8265d56a7749cdff5cdb1caa14008a301812faa045b13e612`
- `data/derived/native-effective-potential-_livswt7/failure_attribution.json`: `fe99de385312521b37e98f131a5381b6e46a192805dcc481f56a3cb361b01a56`

### Targeted wrap gate preparation

The next diagnostic uses the unchanged attested `effective-potential-build-3a9juno_`
executable. It evaluates the original requested seed and both signs of 1e-3,
1e-4 and 1e-5 rad changes in all eight arm flexion/adduction/rotation and elbow
coordinates, then repeats the exact requested seed: 50 copied-state calls,
45 s process-group deadline, 4 GiB limit, single thread and nice 10. No trajectory
or optimization is authorized by this fixture. Every raw response and pending
request is retained; actual assembled coordinate changes are used in virtual
work. Endpoint-average moment arms distinguish ordinary first-order secant
error from a discrepancy that persists as the step shrinks. Sided continuity
is evidence about a branch, not an explicit wrap-branch identifier.

`GeometryPath.cpp`'s 0.0005 m outer wrapping iteration stop cannot explain these
single-wrap paths: source `maxIterations` is one when there is one wrap object.
BIClong uses `WrapEllipsoid`, whose hybrid algorithm blends axial and sampled
fan constructions. BRA uses `WrapCylinder`; its stored wrap length is the
spiral length computed from radial angle and axial separation. Native moment
arms use `MomentArmSolver`, rather than finite differences of stored length.
Neither a discontinuity nor a conservative gradient is assumed from these
source observations alone.

`scripts/effective_passive_energy_observation.py` exposes raw native passive
energy, the exact source-derived Thelen normalization correction, and their
sum, with model and retained source hashes. It leaves native getters, force
laws and physiological ledgers untouched. Millard observations have a zero
correction. This observation interface can be reused in independent energy
closure diagnostics without claiming that the corrected numerical merit has
passed the complete native virtual-work gate. Its fixture checks the analytic
passive-force derivative and preservation of all 98 original observations.

### Actual targeted wrap result: persistent force/path discrepancy

`data/derived/native-wrap-work-ovst3fef/report.json` retains the 50-call result
(8.934 s), raw and corrected observations, and actual assembled coordinates.
The native process closed/reaped, continuing state remained unchanged, and no
physical time advanced. Repeating the exact requested seed changed actual q
by at most 6.624e-10 and raw effective energy by 8.680e-11 J.

Both sides at 1e-4 and 1e-5 rad converge to a nonzero discrepancy between stored
path-length derivative and native force moment arm. Fine sided means change
by less than 1e-8 m/rad for all eight tested directions, while every mean remains
above the original 1e-6 m/rad gate. Representative fine means:

| Path / coordinate | Length–moment-arm discrepancy (m/rad) |
|---|---:|
| BIClong left / arm adduction | 0.0002335568 |
| BIClong right / arm adduction | 0.0001071940 |
| BRA left / elbow flexion | 0.0001344210 |
| BRA right / elbow flexion | 0.0001091474 |

The corrected muscle energy derivative still agrees with tendon force times
measured length derivative to roughly 1e-9 Nm at the fine steps. Therefore an
ordinary finite-difference truncation or muscle primitive error does not explain
the remaining native virtual-work failure. Right BIClong has larger sided
deviations at some 1e-3 rad steps; those may indicate wrap-construction transitions,
but the probe does not expose a native branch identifier and does not assert one.

`scripts/summarize_wrap_virtual_work.py` reproduces `convergence.json`. The complete
energy-gradient gate remains failed. A justified next numerical route is a
trust region on actual native generalized-force residuals, with source-backed
sensitivity scaling and unchanged physical acceleration/support acceptance.
A scalar path-length potential cannot supply an exact gradient for this native
force mapping until a separate wrap tangent-construction correction is validated.
No muscle, wrap, force law, source bound or acceptance threshold was changed.

### Actual-force trust region implementation (source-only preparation)

`scripts/solve_native_force_supine.py` uses the attested 98-muscle/seam-omitted
probe to solve the 31 independent generalized-force equations. Native tree
residuals are pulled back through the actual linear knee couplings; ideal
coupler reactions are not fitted as applied forces. The 28 optimization
coordinates exclude the same three plane translation/heading gauges as prior
work. All 31 force equations remain in the merit, and all native accelerations
and gauge reactions remain in physical acceptance.

Coordinate scales start from XML coordinate spans and are equilibrated using
actual Jacobian column sensitivity. Residual scales are the resulting row
sensitivity norms, with explicit N/Nm and m/rad units. Scales freeze within a
chunk. The source fixture verifies invariance to changing residual units from
N to mN. A exactly zero-sensitivity equation is retained in its native unit;
a zero-sensitivity free coordinate rejects the formulation. The actual retained
98-muscle Jacobian is rank 28, with scaled singular values 3.50355 to 0.00681764.
Its initial bounded model predicts numerical cost 1.78885e-5 to 1.04091e-5.
This is a local model prediction, not native improvement or equilibrium.

Each iteration obtains a fresh actual-coordinate Jacobian and solves bounded
linear least squares in the scaled trust box. The original XML bounds and
0.03 native-unit maximum coordinate step hold. Trial acceptance requires a
positive actual reduction and actual/predicted reduction at least 0.1; poor
agreement shrinks the box and recomputes its step. SVD/rank, linearized residual
floor, active bounds, gradient and prediction errors are retained. No scalar
mechanical energy or KKT-bound force is invented. Intermediate numerical
iterates may lose support; each retains its full unaccepted physical gate.

A source-bounded 200-call/60 s chunk performs at most six fresh-Jacobian
iterations with eight attempts each, counts domain rejections, records every
requested and actual q plus native response, and kills/reaps its owned process
group at the wall cap. Only the exact known skin/bed-domain rejection permits
backtracking; other errors stop. All native responses are from the new 98-muscle
identity; historical 92-muscle responses are never reused. Existing 98-muscle
samples are used for offline scaling design, not substituted for live responses.

Original final limits remain unchanged: every native acceleration <=1e-4,
normalized support and gauge residuals <=1e-4, and native position/velocity/
acceleration constraint errors <=1e-5. A bound-stationary nonzero residual
fails. Even a passed static gate requires separate forward acceptance before
an immutable startup reference can be adopted. No native run of this solver
has occurred at this preparation checkpoint.

### Actual first force-root chunk: numerical improvement, physical regression

`data/derived/native-force-root-919h2o1p` retains 179 calls in 30.941 s,
including all actual q/force/acceleration records, six Jacobians, scales and
trial ratios. The process closed/reaped and continuing state remained unchanged.
All six numerical iterates were accepted by their declared model-ratio test;
four additional attempts in iteration 2 were rejected. No physical iterate
passed equilibrium acceptance and no trajectory/reference was promoted.

The scaled force merit fell from 2.86598979e-5 to 8.15047330e-6, while maximum
native acceleration increased from 54.1717563 to 88.5339739. Normal support
residual increased from 1.98840e-5 to 0.182804410 of weight (139.182837 N).
Final pitch and roll normalized residuals were 0.00481840 and -0.00134266;
maximum gauge residual was 0.000255531. Native constraint error remained small,
2.13902e-14. Thus this is an explicitly failed physical result, not settlement.

The source-derived sensitivity scaling was unit-invariant but not an adequate
physical-error metric: the normal-support row scale was 409158.6 N, so losing
139.18 N of support cost little relative to reducing flexible-joint errors.
Lumbar extension torque decreased from 12.8272 to 4.64611 Nm and bending from
5.81625 to 1.65872 Nm, but left ankle acceleration reached -88.534 and its toe
75.6604. This exposes a meaningful tradeoff failure of residual equilibration,
not merely a missed optimizer tolerance. Accepted reduction ratios ranged
0.272 to 0.799; favorable prediction agreement does not imply useful physical
progress.

Do not continue this formulation unchanged. A subsequent formulation must
retain support feasibility (for example exact root-force equalities with native
nonlinear verification) and connect residual weighting to inertial acceleration
error rather than inverse stiffness alone. Any mass/inverse-mass observation
must be source-verified before use; no guessed inertia, artificial balancing
reaction or relaxed physical limit is justified by this failed chunk.

### Native physical-metric prerequisite (prepared, not yet compiled)

The retained Simbody `SimbodyMatterSubsystem.h` documents `calcMInv` as the
inverse free-mobility mass operator. It does not by itself enforce ideal
constraint reactions. The isolated `native_physical_metric_probe.h` therefore
emits native `calcM` and `calcMInv`, both tree and constrained force residuals,
and actual udot. The two-pose prerequisite must verify the sign relation
`-MInv * constrained_residual = udot` rather than assume it. It also reports
`-MInv * tree_residual - udot` to expose omission of constraint reactions,
checks inverse symmetry, positive mass eigenvalues and `M * MInv = I`.
Failure requires a properly constrained projection; it never licenses using
the unconstrained tree residual as acceleration.

The prepared isolated build is `effective-potential-build-_jnt2x88/manifest.json`.
Only a copied diagnostic header changes; forces, model, library, activations,
contact and prior frozen probes remain unchanged. The planned bounded identity
probe uses two copied states (the initial supported q seed and the explicitly
failed force-root q), 2 calls/15 s, with no optimization or trajectory. Numerical
identity tolerance is 1e-8 per native acceleration component, separately in
m/s² and rad/s². Physical acceptance remains 1e-4 in those same respective
units. No mass operator will be used by the solver before this gate passes.

### Native inverse-mass gate passed; support-feasible physical solver prepared

Actual receipt `native-physical-metric-04jdqqbn` contains two evaluations in
0.542 s after successful isolated compile. The constrained sign relation holds
within 7.75e-14 per acceleration component. Inverse symmetry and M*MInv errors
are below 5.69e-14; both mass matrices have strictly positive minimum eigenvalues
(1.4217e-4 and 1.3290e-4 in the native mixed-coordinate representation). Applying
MInv to the *tree* residual instead gives errors 102.079 and 70.175, confirming
that retaining native constraint reactions is essential. All processes reaped;
continuing state unchanged. This passes only the metric prerequisite.

`scripts/support_feasible_physical_step.py` and
`scripts/solve_support_physical_metric.py` now use
`a = -MInv * constrained_residual`, verified at every sampled pose. The local
model `B = -MInv * d(residual)/dq` freezes MInv for each trust-region iteration;
it does not silently claim to differentiate MInv. Actual retained data measures
its difference from the full acceleration Jacobian at 0.9702% in matrix norm.
The actual/predicted reduction ratio checks this approximation during trials.
Only columns are preconditioned; the residual merit is all native accelerations
in their respective m/s² and rad/s² units.

Every local QP enforces normal support, pitch and roll equalities. At most three
additional bounded root-coordinate corrections use the native support block
before accepting a trial. Actual support AND gauge residuals must remain within
1e-4, and actual acceleration cost must decrease with ratio at least 0.1.
Otherwise the trust box shrinks and its QP is recomputed. No constraint multiplier
or bound reaction can supply physical equilibrium. XML ranges, original .03
native-unit maximum coordinate box, unchanged muscles and original wrap physics
are retained. Full acceleration/constraint acceptance and forward promotion
remain separate gates.

The offline design receipt `physical_qp_design.json` in the native metric folder
uses an exactly matching q/force anchor. The .03 box predicts acceleration cost
6293.4602 ->1126.6866; using the separately measured full acceleration Jacobian
on that same step predicts957.4281. Boxes .015/.0075/.00375 also predict decreases;
linear support error is <=6.4e-17. These are unaccepted local predictions. Source
fixtures additionally reject a lower-cost candidate with .18 normalized support
error. An initial QP column scaling that left tiny normalized Hessian diagonals
hit its offline iteration cap; scaling columns by the common objective norm
restored unit column norms without changing the mathematical objective or bounds.

The prepared gated solver requests at most 200 calls/60 s, six refreshed Jacobians,
eight recomputed attempts and four total support evaluations per attempt. It
records raw native M/MInv, q, force and acceleration responses, rejection records,
scales/rank/model diagnostics and exact last accepted numerical candidate. It
starts from the original supported98 seed, never the unsupported failed candidate,
and requires the successful metric receipt to match the executable manifest hash.
No full native solve has occurred at this source checkpoint.

### Actual support-feasible chunk and bounded continuation design

`support-physical-root-nayktxy6` completed 190 calls in36.596 s. Native state
remained unchanged and the process reaped. Acceleration cost fell
6293.4602 ->1250.3931 ->190.3182 ->97.2130 ->64.1873 ->54.5911 ->53.0012.
Maximum native acceleration fell54.1718 ->5.68810. Final normalized support
was at most7.018e-7, gauge1.242e-7, and native constraint error3.450e-14.
Physical equilibrium still fails the1e-4 acceleration criterion. No forward
or reference promotion occurred.

The final ratio was0.1137, versus earlier accepted ratios0.440–0.972; the
solver's existing rule reduces the next radius from.015 to.0075. This weak
late prediction must not be hidden by the large earlier improvement. No actual
coordinate is within.005 native units of its XML limit. The final QP rank is28,
with singular values2.461 to.003105; active hips/arms/lumbar limits are local
trust boxes, not anatomical stops. Remaining dominant acceleration is right
arm adduction5.6881 (torque−2.7724 Nm), followed by right arm rotation3.6109
and flexion−3.3669; lumbar extension acceleration is−2.5564.

A single smaller-radius continuation is prepared through `--resume` using the
exact saved requested coordinates, verified executable-manifest identity and
SHA-bound previous candidate/report/response journal. Prior response files are
immutable and none is substituted for a fresh response. The last retained
Jacobian precedes the final accepted step, so a fresh28-column Jacobian at the
new point is required. Fresh resumed support/gauge feasibility is checked;
actual replay q/acceleration differences are reported. The radius is derived
from the prior accepted ratio by the same existing rule. If prediction quality
or descent remains weak, further unchanged continuation is not justified;
actual arm path/source-registration errors require renewed diagnosis.

### Smaller-radius continuation: support retained, physical outlier stalled

`support-physical-root-xgpiykt1` completed182 calls in36.358 s; native process
reaped and continuing state unchanged. Exact requested-q replay differed by
7.137e-11 in actual q and5.779e-7 in maximum acceleration. Previous receipts
remain immutable. Final support is<=7.961e-8 normalized, gauge<=2.609e-9 and
constraint error2.441e-14. The physical acceleration gate still fails:
maximum5.903845, compared with5.688101 before continuation.

Cost fell53.0012 ->42.9547 ->41.7605 ->41.1224 ->41.0010 ->40.4029 ->40.3499.
After the useful first smaller-radius step, several predicted decreases failed
actual merit tests. The last iteration rejected radii.00375/.001875/.0009375
(negative actual reductions) and.00046875 (ratio.0869), accepting only.000234375
with a0.05297 cost reduction and ratio.4511. This is a32-fold contraction from
the chunk's initial radius; the last gain is only0.131% of current cost.

Right arm adduction remains the dominant acceleration5.90385, followed by
right arm flexion−3.23334, left adduction3.08073 and lumbar extension−2.82789.
No equilibrium, forward run or reference promotion is claimed. Stop unchanged
native continuations: investigate local model/derivative agreement and actual
arm path/source-registration effects before spending another solve chunk.
The separately demonstrated native wrap force/path mismatch is relevant
background, but these solver results alone do not prove it causes the current
stagnation; candidate wrap corrections remain separate unpromoted physics.

### Directional derivative diagnosis and corrected model

`support-physical-root-xgpiykt1/derivative_diagnosis.json`, reproduced by
`scripts/diagnose_physical_metric_derivatives.py`, identifies a concrete cause
of the late prediction failure. For `a(q)=-MInv(q)r(q)`, the derivative is
`da/dq=-MInv dr/dq-(dMInv/dq)r`. The second term was omitted from the previous
local model. It is only0.783% of the full Jacobian in global norm, yet dominates
some chosen directions. Adding the observed mass derivative reduces relative
Jacobian error to1.5951e-6 (remaining finite-product/finite-difference terms).

For the retained.00375 step, frozen-mass acceleration prediction error is1.02037,
versus0.10503 using the full observed Jacobian. At.000234375 it is0.07938 versus
0.0008744. The full Jacobian correctly predicts larger rejected steps increase
cost (for.00375:40.4029 ->41.0016, actual41.1593). Thus the observed ratio collapse
cannot be assigned to wrap errors or muscle coverage without first fixing this
known derivative omission.

Recomputed hard-support QPs using the full observed Jacobian predict decreases
at all five tested boxes. At.00375 the predicted cost is39.78845 instead of an
increase; smaller boxes predict39.9964/40.1197/40.2001/40.2587. Active-manifold
stationarity norms are below2.9e-7 in the conditioned QP, with linear support
errors below1e-8 and no active source-ROM bounds. This supports a derivative
correction before blaming the bounded linear solver. These are offline
predictions at the retained pre-final point, not new native acceptance.

Shoulder force attribution at the final point closes to below6e-15 Nm:
all shoulder muscle moment comes from TRIlong, BIClong and BICshort per side.
For right adduction, muscle torque is+2.67090 Nm, gravity+1.12977, contact−0.983994
and passive joint−0.000292, yielding native tree residual−2.81639 Nm. The
limited dedicated shoulder coverage is real, but this accounting does not
prove absence of any support-feasible equilibrium. No forces are added or
controls adjusted to compensate.

The solver now uses the full observed native acceleration Jacobian, while
verifying the native force/inverse-mass sign identity at every sample. This
includes mass configuration and constraint-reaction dependence without extra
native evaluations. The former approximation and failed trial records remain
intact. A bounded comparison may explicitly restart the numerical box at.00375
after this derivative formulation change; the override is recorded and all
physical/source bounds remain unchanged. The proposed resume is the latest
supported saved q, with unchanged original98 wrap physics. No candidate wrap
or metadata epoch is combined with this numerical correction.
