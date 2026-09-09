# Explicit supine initialization

The default 92-muscle supine plant starts from the walking source's zero-angle
pose with posterior collision spheres. At 77.6122029 kg, its initial contact
force is only 21.28 N against 761.38 N weight. Its first 20 ms include 1.64 mm
pelvis translation, 27.76 mm forearm-origin translation and substantial arm
rotation. Registered sampled skin motion reaches 164.50 mm. This is not a large
initial floor reaction; unsupported body weight and source joint/muscle forces
drive the initial motion.

`NativeMechanicalStream(initial_pose={coordinate: value})` now supports an
explicit alternative initial coordinate state. It copies and hashes the request
into `inputs/initial_pose.txt`. Native code applies it after contact-reference
construction, before muscle equilibrium and the metabolic reference sample.
Dependent, locked, prescribed, out-of-range and assembly-altered coordinates are
rejected. Time and cumulative work/energy ledgers begin at zero. The initializer
does not add an ongoing pose constraint, tether or support force.

Material registration remains attached to the retained source-default pose.
Native frames expose `registration_reference_bodies` captured before applying
the initial coordinates. Consumers must use that reference for their material
embedding and use actual `bodies` for the current pose. Rebuilding registration
from the altered pose would hide that pose in canonical identity geometry and
change force levers between comparison arms.

`scripts/verify_native_initial_pose.py` verifies copied-state static evaluation
and explicit initialization give the same coordinate/contact state, preserve
the original support plane and registration reference, and sample the metabolic
baseline after initialization. The reference-preservation check passed on
build `build-f2du33q4`.

The solver `scripts/solve_supine_sphere_initial_pose.py` optimizes all independent
source coordinates with source bounds and limits pelvis/lumbar rotations to
±0.15 rad. Its constrained artifact is:

`data/derived/supine-sphere-pose-7i99jsg8/initial_pose.json`

This candidate has 760.10 N initial support and maximum generalized acceleration
17.59 rad/s², down from 2964.35 rad/s². It **is not an equilibrium**: residual
lumbar and pelvis acceleration remains. The 92-muscle model has no lumbar muscle
actuators or source lumbar passive expression force. A one-microsecond forward
check reproduces the static lumbar acceleration, confirming this is a physical
model residual rather than a static/forward protocol discrepancy.

The corrected, unchanged-reference forward comparison
`data/derived/supine-initial-motion-3xx4d_td/report.json` measures:

| Measurement | Default | Constrained candidate |
|---|---:|---:|
| Maximum skin sample displacement, first 20 ms | 164.50 mm | 2.924 mm |
| Kinetic energy after 20 ms | 4.606 J | 0.649 J |
| Maximum skin sample displacement, 80–100 ms | 89.70 mm | 26.51 mm |
| Kinetic energy after 100 ms | 17.90 J | 13.32 J |

The onset improves substantially, but the body still drifts. This is a short
native body/skin observation test without world feedback; it does not establish
stable blanket contact, a settled human, or sustained operation. The initializer
remains optional and this candidate has not been promoted to the default.

Earlier unconstrained optimizer artifacts rotated the pelvis about 66°. They
have explicit `POSTURE_REJECTION.json` sidecars and must not be presented as flat
supine equilibria. An early unconstrained body/blanket test also failed its
environment integration envelope; it provides no stability evidence.

Reproduce the constrained comparison with:

```sh
.venv/bin/python scripts/verify_native_initial_pose.py
.venv/bin/python scripts/verify_supine_initial_motion.py \
  data/derived/supine-sphere-pose-7i99jsg8/initial_pose.json
```

## Measured-bed solver performance, 2026-09-08

The native skin/mattress equal-pressure solve now uses safeguarded Newton steps
with the analytic derivative of the unchanged skin pressure and piecewise linear
bed stress. It retains both original material-domain checks and the original
50-bisection algorithm as a fallback. The convergence threshold is 1e-10 Pa;
this does not change the native time integrator's tolerances.

`data/derived/bed-newton-production-parity/report.json` compares the retained
previous header with the production header over 20,021 approaches for each of
SM, MM, and HM. Domain acceptance matched exactly. Maximum pressure difference
was 1.06e-10 Pa, total energy difference 5.69e-13 J/m², and skin indentation
difference 5.13e-17 m. The compiled comparison source and both headers are
retained alongside the report.

`scripts/verify_native_bed_solver_parity.py` replays identical retained inputs
against the separately hashed historical and current native binaries. In
`data/derived/native-bed-solver-parity-eeho035r/report.json`, a 98-muscle body
advanced four 5 ms steps on MM. Maximum coordinate difference was 6.25e-16,
kinetic energy difference 7.55e-15 J, potential energy difference 2.70e-13 J,
and surface numeric difference 5.01e-12. Wall time including startup was
30.13 seconds for bisection and 13.82 seconds for Newton. This bounded comparison
supports equivalent mechanics and reduced runtime for this workload; it does
not establish a universal speedup or a stable supine equilibrium.

## External wrench is not pelvis angular acceleration

A later static search exposed an incorrect objective interpretation: reducing
pelvis angular acceleration does not establish balanced external support
moments. Internal muscle torques can hold the pelvis while other segments
accelerate. At fixed pose, a changed activation vector altered pelvis angular
acceleration while the root generalized static moment remained 262.51 N·m.
The six root static generalized forces were activation-invariant to 1.14e-13.
Evidence and the explicit correction are retained in
`data/derived/supine-hip-tones-w3_krvzv/root_wrench_activation_invariance.json`
and `OBJECTIVE_CORRECTION.md` in that directory.

Current equilibrium searches use actual root zero-acceleration generalized
residual forces and moments, with force scale mg and moment scale mg times
body length. Earlier reductions in skin onset, acceleration, or total residual
norm do not establish equilibrium. A measured support-wrench solution also
requires joint equilibrium and a successful actual initialization before it
can be used as a stable supine starting state. No candidate discussed above
has been promoted on the strength of reduced motion alone.

## Supported posture and tonic activation refinement

A constrained solve now balances the actual external root wrench before
addressing internal joint residuals. It uses the current 98-muscle model,
`supine-surface-contact-5jqy1juo/manifest.json`, the measured MM curve, and
26 explicitly selected tonic muscle activations. The native material reference
remains the original source geometry. Pose changes are initialization only;
no prescribed coordinates, root tether, or ongoing support controller is added.

The final static candidate is
`data/derived/supine-selected-defaults-lp6uv34z/`: maximum initial acceleration
2.996 rad/s², total acceleration norm 5.420, support 761.376 N, and normalized
root wrench residual 6.73e-8. Actual native initialization agrees with the
static pose, forces and all 26 activations in
`supine-lumbar-default-receipt-t435gyzt/report.json`. Remaining joint residuals
reach 1.65 N·m, so this is **not exact equilibrium**. The variant remains
explicit and is not a default or a learned cortical controller.

Two earlier, separately identified candidates passed finite-duration reciprocal
bedroom runs with cloth and held native tonic defaults:

| Candidate | Duration | Peak body KE | Final body KE | Report |
| --- | ---: | ---: | ---: | --- |
| 53 rad/s² initial maximum | 0.2 s | 0.577 J | 0.169 J | `supine-tonic-world-txit_xh3/report.json` |
| 26 rad/s² initial maximum | 0.5 s | 0.103 J | 0.0183 J | `supine-tonic-world-vcvr3qwx/report.json` |

The second run still changed a hip coordinate by 0.161 rad and moved the cloth
by up to 0.148 m. These are bounded forward runs, not proof of immobility or
long-term balance. Their different postures also initialize different cloth
contact configurations. All reports are under `data/derived/`. A brain adapter
that sends zero motor commands will replace these tonic defaults; integrating
the baseline requires an explicit, disclosed spinal policy and a severed-brain
control arm. Held-tone tests cannot establish cortical causation.

The final `lp6uv34z` candidate subsequently completed 1.0 s / 200 reciprocal
5 ms body–bedroom exchanges in `supine-tonic-world-r90qq71u/report.json`.
Peak body kinetic energy was 0.0405 J and final energy 0.00243 J; final sampled
skin motion was 0.316 mm per exchange. The cloth remained finite and moving
(final maximum node speed 0.587 m/s). Mattress indentation reached 80.9 mm.
However, both shoulder rotation coordinates changed about 0.399 rad during
relaxation. This is successful finite-duration mechanics with an appreciable
postural relaxation, not an initially motionless equilibrium. The actual
258.45 W initial muscle metabolic reference and all tonic activations are
retained in the report. The diagnostic omits brain and physiology integration;
it cannot establish their behavior under this tonic reference.

The bedroom runs above use the existing explicit projection cloth integrator.
A separate world-physics stationary-body audit identified possible energy
injection in that integrator. Therefore finite state and matched momentum in
these runs do **not** establish passive cloth contact or energy-conserving world
coupling. Validation with its passive replacement is a separate requirement.

A later joint pose-and-tone solve reached no-cloth static equilibrium in
`data/derived/supine-selected-defaults-ih90ww49/`. At the actual initialized
native state, maximum acceleration was 1.92e-7 rad/s² and maximum generalized
force residual 3.72e-8; an independent 1 µs forward step gave maximum speed
derivative 1.54e-7 rad/s². All 26 tonic activations match exactly and the source
material reference is unchanged. Evidence is retained in
`supine-equilibrium-initial-xfhqpo4e/report.json`. This proves the static
no-cloth initialization to the stated numerical tolerance, not sustained
loaded-body stability or cortical control.

The `ih90ww49` equilibrium then passed **1.0 s of native-only forward dynamics**
with no cloth or external ports in `supine-native-equilibrium-3gw30l4k/report.json`.
Across fifty 20 ms exchanges, maximum coordinate change was 1.34e-7, maximum
kinetic energy 7.12e-15 J, and maximum momentum residual 6.71e-12 N. The
predeclared limits were 1e-5 coordinate change, 1e-6 J and 1e-5 N respectively.
This establishes a numerically stable supported body for that duration with
its explicit tonic activations. `scripts/verify_supine_native_equilibrium.py`
retains initial/final native state and the exact verification source.

The same candidate completed 1.0 s with the existing projection cloth in
`supine-tonic-world-6u1wzfk9/report.json`. Its accompanying `cloth_scope.json`
and retained environment source explicitly exclude an energy-passivity claim;
validation with the passive cloth replacement remains separate.
