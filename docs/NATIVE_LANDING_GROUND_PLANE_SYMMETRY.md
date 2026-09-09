# Why the native landing failed, and the symmetry the stance design was ignoring

Scope of this document: it explains a measured failure of the 98-muscle native
landing program, names two defects in the feedback design, and records what was
changed. It does **not** claim walking. Nothing here shows sustained forward
locomotion, and the receipts it points at say so in their own `scope` strings.

## 1. What the failed landing actually did

`data/research/locomotion_control/patient_right_landing_program/` ran
`lift_right` (9 s + 2 s hold), `swing_right_forward` (6 s + 2 s) and
`land_right_forward` (6 s + 4 s) and stopped at 24.13 s of a 30 s horizon with
`completed_horizon: false`. Reading the retained trajectory rather than the
report changes the picture the report alone suggests.

| t (s) | segment | blend | right clearance (m) | clipped | ‖x − x_ref‖ | hip_adduction_r |
|------:|---------|------:|--------------------:|--------:|------------:|----------------:|
| 20.0 | land (start) | 0.00 | 0.0309 | 0 | 0.033 | −0.042 |
| 21.4 | land | 0.09 | 0.0218 | 1 | 0.114 | −0.040 |
| 22.0 | land | 0.21 | 0.0073 | 8 | 0.395 | +0.034 |
| 22.6 | land | 0.38 | **0.0052** | 20 | 1.72 | +0.129 |
| 23.0 | land | 0.50 | 0.0473 | 31 | 3.99 | +0.252 |
| 24.1 | land | 0.81 | 0.2071 | 60 | 14.9 | +0.032 |

Three things follow directly from the trace.

* The reported `final_per_foot_clearance_m.r = 0.242` is **not** a lifted foot.
  It is the last frame of a fall: `pelvis_ty` is 0.874 m and `pelvis_tilt` is
  −0.963 rad there. The run ended on the evaluator's fall guard.
* The right foot got to within **5.2 mm** of the ground at t = 22.6 s and then
  went back up. The landing segment never produced a load transfer.
* The state that ran away is `hip_adduction_r`, and it ran away with the
  **opposite sign to its own reference** (the landing equilibrium asks for
  −0.121 rad; the plant went to +0.33 rad). `pelvis_list` and
  `hip_adduction_l` followed it. This is a frontal-plane collapse off a
  single-support base whose measured COM margin at t = 20 s was **11 mm**.

Clipping is a symptom, not the trigger: at t = 21.8 s the reference error was
already 0.28 with a single clipped command.

## 2. Defect one: the handoff has no stable blend

The evaluator interpolates both the reference and the gain on a clock:
`K(b) = K_old + b·(K_new − K_old)`, `x_ref(b) = x0_old + b·Δx`. The pre-existing
`landing_gain_interpolation_audit.json` reported the raw spectral radius of the
blended closed loop, which is dominated by three structural unit eigenvalues
(section 3) and is therefore uninformative. Recomputing it while excluding those
three (`landing_gain_interpolation_audit_v2.json`) gives:

| blend | swing plant, ρ excl. symmetry | landing plant, ρ excl. symmetry |
|------:|------------------------------:|--------------------------------:|
| 0.00 | 0.978638 | 0.999999999 |
| 0.20 | 0.987619 | 0.999999999 |
| 0.50 | 0.997830 | 0.996575 |
| 0.75 | 0.999998 | 0.985878 |
| 1.00 | 1.000000 | 0.982416 |

The swing plant (right foot airborne) is what was actually under the controller
for the whole landing segment, because the foot never touched down. On that
plant the slowest controllable mode's time constant goes from 0.47 s at blend 0
to effectively infinite at blend 1. Meanwhile the interpolated pair
(x_ref(b), u_ref(b)) is not an equilibrium of either plant, so it injects a
constant forcing term `b·(A_swing Δx + B_swing Δu)`. A constant forcing into a
mode with no restoring action integrates without bound. That is the ramp in the
table above.

There is no blend value at which both plants are well damped. The handoff is
only safe if the plant has already become the landing plant — that is, if the
foot is already in contact when the blend starts.

The repository already contains the fix for the geometry: the staged approach
(`patient_approach10mm` → `patient_approach1mm`) ends with the right foot at
0.89 mm clearance carrying 1.7 % of body weight, settled
(COM speed 1.2 × 10⁻⁵ m/s) with zero clipped commands. The failed landing
program simply did not use it: it went from a 31 mm swing straight to a loaded
landing equilibrium.

## 3. Defect two: the stage cost penalises an exact symmetry no muscle can reach

Every one of these linearisations has **three discrete eigenvalues at exactly
1.0 whose left eigenvectors are orthogonal to `B` to 10⁻¹⁰ or better**. Their
right eigenvectors are supported on `pelvis_tx/value`, `pelvis_tz/value` and
`pelvis_rotation/value`.

This is not numerical noise. On flat ground, translating or yawing the whole
body changes neither gravity nor the contact geometry, so the ground-plane
isometry group is an exact three-parameter symmetry of the plant. `A` has it in
its kernel and `B` has no component along it.

The stage cost in `linearize_native_stance.py` puts `Q = 1000` on exactly those
coordinates. Asking for a stabilising Riccati solution that regulates an
uncontrollable, marginally stable direction asks for something that does not
exist, and `scipy.linalg.solve_discrete_are` duly fails to produce one:

| policy | waypoint | `dare_residual_relative` |
|--------|----------|-------------------------:|
| `linearization_8_mijxgp/discrete_margin` | right lift | 0.0124 |
| `linearization_j0b9dci6/discrete_margin` | right swing | 0.0106 |
| `linearization_g_x24sqn/discrete_margin` | landing | **0.633** |
| `linearization_ufj51r6y/discrete_margin` | crouched landing, 10 % load | **0.591** |

Both loaded-landing designs — and only those — are effectively unsolved. The
continuous designs failed outright with "the associated Hamiltonian pencil has
eigenvalues too close to the imaginary axis", which is the same statement.

There is a second, more consequential cost. Because `Q` demands regulation of
absolute world position, the resulting gain has real feedback on it. For the
10 % landing design:

```
||K N||_2                        = 11.39      (N = the symmetry subspace)
1 cm of absolute pelvis_tx  ->  max |Δu| = 0.0475, 12 muscles shifted by > 0.01
1 cm of absolute pelvis_tz  ->  max |Δu| = 0.0629, 13 muscles shifted by > 0.01
```

Most baseline excitations in these equilibria sit between 0.01 and 0.43, so a
few centimetres of travel is a large fraction of the whole command. Nothing the
muscles do can reduce this error, so it is pure, accumulating disturbance. **Any
controller with nonzero feedback on absolute ground-plane position cannot be
used to walk**, because walking is precisely the act of increasing `pelvis_tx`.
This is a structural blocker on repeated stepping, independent of the landing
question.

## 4. What was changed

`scripts/design_native_deflated_lqr.py` (new) measures the symmetry subspace
from the retained `A`/`B` (smallest singular directions of `Ad − I`, required to
be unreachable from `Bd` and separated from the rest of the spectrum by a
factor of 100), deflates it out of both the design coordinates and the stage
cost, solves the reduced DARE, and lifts the gain back. No physical model
parameter and no equilibrium is touched: `x0`, `u0`, `A`, `B` are copied
through unchanged and only `K` differs. `K N = 0` holds to 10⁻⁹ by
construction, so the design has exactly zero feedback on absolute ground-plane
position and yaw, and the closed loop keeps those three unit eigenvalues.

| waypoint | deflated dim | reduced DARE residual | reduced ρ | max |K| |
|----------|-------------:|----------------------:|----------:|--------:|
| right lift (`8_mijxgp`) | 3 | 4.4 × 10⁻⁸ | 0.9787 | 354.8 |
| right swing (`j0b9dci6`) | 3 | 2.0 × 10⁻⁷ | 0.9786 | 394.2 |
| approach 10 mm (`rnwkf3bb`) | 3 | 1.2 × 10⁻⁸ | 0.9785 | 539.8 |
| approach 1 mm (`5j_icido`) | 3 | 1.9 × 10⁻⁹ | 0.9794 | 272.4 |
| landing, 10 % load (`ufj51r6y`) | 3 | 6.6 × 10⁻¹¹ | 0.9890 | 215.1 |
| landing, 25 % load (`landing25`) | 3 | 4.4 × 10⁻¹¹ | 0.9947 | 213.6 |
| landing, 45 % load (`landing45`) | 3 | 2.6 × 10⁻¹⁰ | 0.9969 | 242.3 |

For the 10 % landing waypoint this replaces a design whose closed loop had a
5.6 × 10⁶ s time constant on its slowest mode with one that has 0.9 s, at a
slightly *lower* peak gain.

Two new load-transfer equilibria were solved with the existing crouched-landing
solver, over the same footprint the swing produced:

| waypoint | right load | pelvis_tx (m) | pelvis_tz (m) | residual accel (rad/s²) |
|----------|-----------:|--------------:|--------------:|------------------------:|
| approach 1 mm | 0.017 | −0.0089 | −0.1378 | — |
| landing 10 % | 0.100 | 0.0057 | −0.1173 | — |
| landing 25 % | 0.250 | 0.0257 | −0.0733 | 1.9 × 10⁻⁵ |
| landing 45 % | 0.450 | 0.0555 | 0.0021 | 5.0 × 10⁻⁶ |

## 5. What the native runs showed

Both arms below ran the same three programs through one continuing native
process with `scripts/run_native_movement_chain.py`: the retained swing program
(lift 9 + 2 s, swing 6 + 4 s), the retained staged approach (10 mm then 1 mm,
5 + 2 s each), and a 10 % load landing (7 + 4 s). They differ only in the
landing gain.

| | `step_deflated_v1` | `step_dare_control_v1` |
|---|---|---|
| landing gain | `ufj51r6y/deflated_margin` | `ufj51r6y/discrete_margin` |
| all stages completed horizon | **yes** | **yes** |
| final simulated time | 47.0 s | 47.0 s |
| landing max clipped commands | 0 | 0 |
| final right-foot clearance | −0.00496 m | −0.00444 m |
| final active foot contacts | 12 of 12 (6 right) | 12 of 12 (6 right) |
| COM support margin, single → double support | 0.0164 → 0.0376 m | 0.0164 → 0.0394 m |
| final left (support) knee angle | 0.1415 rad | 0.1422 rad |
| final COM speed | 1.1 × 10⁻⁵ m/s | 6.4 × 10⁻⁸ m/s |
| **COM forward travel over the whole chain** | **0.0233 m** | **0.0241 m** |

The landing works. The body reaches a settled double-support stance with the
right foot 16 cm ahead of the left, all twelve contact spheres loaded, the
supporting knee flexed to 0.14 rad so the pelvis could lower, and no muscle
command clipped at any point in the landing segment.

**The gain-design ablation came out negative, and that is the honest result.**
Replacing the invalid Riccati gain with the deflated one did *not* decide the
landing: the pre-existing `discrete_margin` gain completes the same landing just
as well once the staged approach precedes it. The proximate cause of the
original failure was geometric, not a bad gain — the failed program blended
toward a loaded-contact equilibrium starting from a 31 mm airborne foot, so the
plant under the controller was never the plant the gain was designed for.
`step_no_approach_control_v1` is the matched arm for that variable: identical
swing, identical landing equilibrium, identical deflated gain, approach
waypoints removed. It fails, hard.

| landing segment measurement | with staged approach | approach removed |
|---|---:|---:|
| segment survived | 11.0 s of 11.0 s | **4.08 s of 11.0 s** |
| terminating error | none | `TimeoutError: Native mechanical stream response timed out` |
| max clipped commands | **0** | **69 of 98** |
| max ‖x − x_ref‖ | 0.135 | **80.6** |
| final COM speed | 1.1 × 10⁻⁵ m/s | **0.412 m/s** |
| final COM support margin | +0.0376 m | **−0.0548 m** (outside the base) |
| final right-foot contacts | 6 of 6 | **0 of 6** |
| final support-knee angle | +0.1415 rad (flexed, as intended) | **−0.109 rad (hyperextended)** |
| final foot forward separation | 0.161 m | 0.354 m (the leg flailed) |

This is the causal statement the landing needed: **the staged contact approach
is what makes the load transfer survivable**, because it makes the plant under
the controller the plant the landing gain was designed for before the load
blend begins. Everything else in the two arms is byte-identical, and both share
a bit-identical swing stage (final right clearance 0.030939745187 in both).

What the deflation *does* buy is stated in sections 3 and 4 and is not a claim
about this landing: a Riccati solution that actually solves its equation, real
closed-loop damping instead of a 5.6 × 10⁶ s time constant, exactly zero
feedback on absolute ground-plane position, and a design procedure that
succeeds on the 25 % and 45 % load waypoints where `solve_discrete_are` again
returns an unusable answer.

## 6. Carrying the step through to a load transfer

Two further static equilibria were solved over the same footprint with the
existing crouched-landing solver — 25 % and 45 % right-foot load — linearised,
and given deflated gains (`solve_discrete_are` again failed on both; the
continuous design reported "Hamiltonian pencil eigenvalues too close to the
imaginary axis"). `step_transfer_v1` chains all four programs through one
native process: swing, staged approach, 10 % landing, then 25 % and 45 %.

All four stages completed their horizons. 71.0 s simulated, 1568.8 s wall,
one native stream (`c76970c01a9240b0800545dbcf596d10`).

| segment ends at | right clr (m) | left clr (m) | contacts | COM margin (m) | clipped | **COM forward travel (m)** |
|---|---:|---:|---:|---:|---:|---:|
| initial stance, 1 s | −0.00919 | −0.00920 | 12 | 0.1160 | 0 | 0.0000 |
| lift, 12 s | +0.00526 | −0.01899 | 8 | 0.0125 | 5 | 0.0033 |
| swing, 22 s | +0.03094 | −0.01911 | 6 | 0.0112 | 0 | 0.0023 |
| approach 10 mm, 29 s | +0.00955 | −0.01908 | 6 | 0.0115 | 0 | 0.0024 |
| approach 1 mm, 36 s | +0.00089 | −0.01875 | 9 | 0.0164 | 0 | 0.0069 |
| landing 10 %, 47 s | −0.00496 | −0.01719 | 12 | 0.0376 | 0 | 0.0233 |
| transfer 25 %, 58 s | −0.00747 | −0.01379 | 12 | 0.0748 | 0 | 0.0467 |
| transfer 45 %, 71 s | −0.00920 | −0.01140 | 12 | **0.0856** | 0 | **0.0803** |

The right foot compresses monotonically while the left decompresses; the COM
support margin grows 7.6× from its single-support minimum; the run ends settled
at 5.1 × 10⁻⁵ m/s with a maximum reference error of 0.081 and not one clipped
muscle command after the lift.

**This is one half-step: the right foot placed 16.1 cm forward and 8.0 cm of the
body's centre of mass moved onto it, entirely through muscle excitation.** It is
the first forward centre-of-mass transport this line of work has produced. It
took 71 s of simulated time, which is 1.1 mm/s. It is not walking, and the
receipt says so.

## 7. What is still not shown

* No result in this document or its receipts is walking. The best forward
  transport shown is 8.0 cm of centre-of-mass travel in 71 s of simulated
  time — 1.1 mm/s, quasi-statically, ending stationary. Walking would require
  sustained travel at roughly three orders of magnitude that speed, with flight
  or single-support phases that exploit momentum rather than avoid it.
* The load-transfer ladder stops at 45 %. Lifting the trailing foot needs more
  than 50 % on the leading foot, and no such equilibrium has been solved; the
  crouched-landing solver refuses `--right-load-fraction >= 0.5` by construction.
* No mirrored (left-swing) waypoint exists, so a second step is not possible
  with the current asset set.
* The quasi-static waypoint architecture itself has not been shown capable of
  gait. Every waypoint is a static equilibrium and every transition is a
  clock-driven blend; nothing here has a stance phase that exploits momentum.
* Perturbation recovery has been shown for the standing controller only
  (`data/models/engineering_stance_v1/acceptance/`), not for anything in the
  stepping chain.

## 8. Receipts

`data/research/locomotion_control/landing_experiment_receipt.json` collects all
of it — the sha256 of every script, program, policy, target and manifest, each
arm's native stream identity and wall time, the solved load waypoints' native
acceptance, and each arm's own `scope` string. Regenerate it with
`scripts/build_landing_experiment_receipt.py`.

The individual runs, each with its retained `program.json`, `evaluator.py`,
policy and target `.npz` files, `trajectory.json`, `chain_report.json` and
`summary.json`:

| directory | result |
|---|---|
| `step_transfer_v1` | 4/4 stages, 71.0 s, 8.0 cm COM travel |
| `step_deflated_v1` | 3/3 stages, 47.0 s, 2.3 cm COM travel |
| `step_dare_control_v1` | 3/3 stages, 47.0 s, 2.4 cm COM travel (gain ablation, negative) |
| `step_no_approach_control_v1` | **fails at 26.08 s**, 69/98 clipped, COM outside support (approach ablation, positive) |
| `patient_right_landing_program` | the original failure this work started from |

`landing_gain_interpolation_audit_v2.json` re-reports the blended-gain spectral
radii of the original failing handoff with the three symmetry modes excluded.
