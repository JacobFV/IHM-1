# Brain, body, and world integration

The runtime connects a persistent IBM controller, segmental cord, native OpenSim
muscles and articulated body, BioGears physiology, and interactive world objects.
This is an engineering sensorimotor simulation. It is not a validated human brain,
and the evidence below does not establish walking.

## Controller ownership

| Selection | Actual motor owner | Scope |
|---|---|---|
| `regional` | Existing regional controller | Baseline integration |
| `implicit` | Original IBM implicit cortex, optionally followed by cord | Untrained motor behavior; motion alone is not success |
| `implicit_curriculum16` | Retained IBM 16-objective kernel, optionally followed by cord | Same untrained motor path, newer kernel; its sixteen objectives are all sensory |
| `implicit_ankle_primitive` | Trained eight-site motor primitive | Privileged ankle feedback; larger cortex is observational |
| `implicit_cortical_ankle` | Trained persistent 128-site IBM E/I cortex, optionally followed by cord | Privileged ankle error/velocity feedback; narrow learned task |
| `implicit_cortical_stance` | Trained persistent 1024-site IBM E/I cortex at 10 ms | Fixed identified mass/model; privileged full native state; bounded small-push balance, runtime integration under verification |
| `implicit_curriculum16_stance` | Same trained 1024-site policy, built from the 16-objective kernel | Same scope as above; retained so the two kernels can be compared on one task |
| `engineering_stance` | Discrete native LQR feedback at 10 ms | Fixed identified mass/model; regional neural output is observational; native physiology and workbench acceptance passed |

The original IBM checkpoint is unchanged. Copied trained artifacts and their
source are retained in `data/models/`; identities and ablations are exposed in
the workbench. None of these modes is a walking controller.

## What has been measured

- Reduced mechanics closes activation → deformation → spindle afference →
  stretch response. The isolated cord fixture has a 30 ms delay. The native
  runtime's 20 ms neural exchange quantizes sensing and actuation; the fixture
  does not prove a 30 ms physical response in that runtime.
- A separate actual native mechanics + IBM cord experiment at 1 ms exchange
  applies a 300 N foot pulse: evoked soleus Ia begins at 12 ms, named stretch
  and alpha response at 42 ms (30 ms later), next mechanical command interval
  at 43 ms, and native activation difference at 51 ms. The matched
  stretch-blocked arm has zero named stretch response. Receipt:
  `data/runtime/native-reflex-latency-20260908/receipt.json`. Cortex is held at
  zero and BioGears is not advanced at 1 ms in this experiment.
- Three native ankle targets, each lasting 0.3 seconds, show 88.7% lower mean
  squared error with the trained 128-site cortex than with its severed control.
  These are short free-fall ankle tasks, not balance or locomotion.
- With cord, physiology, and bedroom included, a paired 0.3-second run measured
  ankle target MSE of 0.063154 rad² full versus 0.262714 rad² severed, about 76%
  lower. Receipt: `data/derived/unified-world-2ww8inol/report.json`. That historical
  receipt retains an earlier conservative `cortical_control_demonstrated=false`
  classifier; its raw task metrics support only this narrow, short comparison.
- An actual browser test accepted nine native frames through 0.18 seconds with
  the trained cortex and cord bypassed, and submitted three mouse-generated ball
  force commands, up to 3.873 N. No HTTP or browser errors occurred. Receipt:
  `data/derived/browser-integration/drag-summary.json`.
- The trained cortex, cord, native body, physiology, and bedroom completed five
  seconds and 250 exchanges with synchronized clocks. Receipt:
  `data/derived/unified-world-uqzq1wgv/report.json`. This unpaired run establishes
  integration survival only; its ankle MSE is not a causal comparison. It predates
  the subsequent respiratory quadrature and intake refresh audit fixes.
- Ball ground impacts use restitution 0.75. A one-meter test rebounds 0.561993 m
  versus the analytic 0.5625 m prediction; the 0.5% shortfall is 20 ms apex
  sampling, not solver error. That single case has since been broadened. See
  "Rigid ball contact" below and `data/derived/ball-dynamics-v1/report.json`.
- Baseline regional physiology completed ten seconds and 500 exchanges with no
  unmet muscle energy. The metabolic reference is a moving zero-command
  counterfactual, not measured human resting metabolism. See
  `METABOLIC_REFERENCE.md`.
- IBM-1's current 16-objective consolidated kernel is retained at
  `data/models/ibm_curriculum16_kernel_v1/` and drives two new controller kinds.
  Rebuilding the whole stance pipeline on it, changing nothing but the kernel,
  reproduces the existing bundle's result and adds nothing to it: matched
  ten-second pushes recover with peak COM 1.7469 mm against the fused kernel's
  1.7844 mm, and 2.1287 mm against 2.2640 mm, with both severed arms falling.
  **A matched control kernel — the same kernel with its site rows permuted, so
  the learned association is destroyed and every marginal statistic preserved —
  goes through the same pipeline and recovers the same push at 1.7452 mm.** So
  the stance policy does not depend on what the IBM kernel learned; it depends on
  the E/I network being a well-conditioned filter that an offline decoder fit can
  invert. The severed control shows the cortical dynamics own the motor output;
  it was never evidence about the kernel's content, and now demonstrably is not.
  Receipt: `data/runtime/motor-learning/ibm-kernel-comparison-20260908.json`.
  The matched arms through the full runtime, with cord, physiology and the world
  in the loop, are retained at `data/derived/unified-world-ofb2dp7z/`: full holds
  five seconds and 250 exchanges without falling at 1.7344 mm peak COM, severed
  falls at 2.90 s at 499.12 mm. Both agree with the isolated harness, so the
  behaviour survives the full runtime. That is an integration result; it does not
  distinguish a learned kernel from the permuted control.
  Full account and limits in `IBM_CURRICULUM16_KERNEL.md`.

## Mechanics and remaining limits

The body and world exchange forces every 5 ms inside each 20 ms neural and
physiological step. Loads retain impulse quadrature and muscle work is summed
over the mechanical substeps. A same-checkpoint comparison against 20 ms world
exchange reduced severe contact transients; this is refinement evidence, not a
convergence study.

The production blanket previously showed nonpassive energy growth and constraint
violations: over a held-body half second its mechanical energy rose 0.94 J with
stationary boundaries, and the paired five-second audit reached 1.724 times rest
edge length in the severed arm against a 1.12 target. Its sharp folds were in
the simulated mesh, not the render.

Both defects were traced and repaired; see the fold and energy repair section of
[ENVIRONMENT_DYNAMICS.md](ENVIRONMENT_DYNAMICS.md). The folds were a zero-energy
mesh mode: with distance springs only, rotating half the sheet about a grid line
changes no edge length, so a crease cost exactly nothing. A discrete hinge
bending term now supplies the missing shell physics. The energy growth came
mostly from repairing a vertical penetration depth vertically instead of along
the contact normal, which inflates the correction by one over the normal's
vertical component and diverges at every skin sample sphere's rim. On the same
held-body fixture the blanket's mechanical energy change is now -0.00198 J
against +0.9355 J, its maximum hinge fold angle after preparation is 35.9 degrees
against 174.5, and no hinge exceeds 60 degrees where 85 did. Vertex/triangle and
edge/edge collision with a continuous swept test replaces the particle-only
guard, conserving linear and angular momentum exactly and never adding kinetic
energy. Overlapping contact constraints no longer make the implicit KKT solve
singular.

This is still not accepted as calibrated cloth. The production explicit path is
not a passive solver: its energy now oscillates without drifting rather than
growing, and the receipt records that gap instead of asserting a tuned threshold.
The implicit path is the one with a per-step energy and momentum acceptance test,
and it completes one second on the actual sampled skin with worst energy excess
-0.00091 J and worst momentum residual 1.2e-14 N*s. Neither implementation is
volumetric soft-body FEM or a validated textile model. Receipts:
`data/derived/cloth-energy-audit-20260908/receipt.json` and
`data/derived/passive-cloth-skin-yzwjp6iy/report.json`.

Whole-body floor contacts prevent the former feet-only falling-through-floor
failure. Initial supine support and unsupported standing remain separate
validation problems. A changed pose must use the fixed canonical material
registration; refitting the registration to hide its displacement is invalid.

Balance controllers, native equilibrium solves, and gait transfer experiments
are being evaluated separately. A static pose, local stable linearization, or
moving limb is not accepted as unsupported walking. Walking requires sustained
forward motion with actual muscle forces and contact, perturbation recovery,
and causal controller ablations.

The engineered controller passed a ten-second native trial on the 98-muscle,
77.6122029 kg model with a 100 ms [5, 0, 5] N push applied at the current pelvis
center of mass in ground coordinates. Peak COM displacement was 1.706 mm;
matched held excitation fell after 2.89 seconds. Both feet remained loaded.
This establishes bounded stance disturbance recovery, not stepping. Earlier
push trials used a local-coordinate tuple as a ground point and therefore had
an additional lever arm; their claims do not describe COM-centered pushes.
Corrected evidence is retained in `data/models/engineering_stance_v1/acceptance/`.

The cortical stance pilot now passes matched native ten-second small-push trials
in both directions. Its persistent 1024-site E/I network owns all 98 corrective
muscle commands; offline engineered control supplied training/readout targets,
but no live LQR runs on this motor path. Full native state feedback and a
same-kernel zero-input reference are explicit engineering interfaces. The
negative push gives peak COM displacement 2.264 mm and final 0.172 mm; severing
falls at 1.94 seconds. The positive push gives peak 1.784 mm and final 0.196 mm;
severing falls at 2.90 seconds. These are bounded native mechanics trials, not
yet whole-physiology acceptance or evidence of walking. Earlier failed neural
controllers and their receipts remain retained. See `IBM_CORTICAL_STANCE_PILOT.md`.

An engineered movement program has shifted about 94% of support load onto one
foot, lifted the opposite foot, and swung it forward 145.6 mm using actual
muscle commands and native contact. Final swing clearance is 30.94 mm. Landing
has not passed; this is not an accepted walking step or cycle.

Native gravity/support planes and the world's authored frame formerly differed
by roughly 4.5 degrees. A fixed proper rotation and translation now express the
world and body in one frame; native common-frame contact checks and actual
browser stance/ball-drag checks pass. Canonical anatomical registration is
unchanged. `scripts/verify_common_frame.py` now measures that correction on the
actual native body in both environments rather than asserting it: the
uncorrected mismatch is 4.4918 degrees supine and 4.4881 degrees upright, and
the residual after the shared frame is applied is exactly 0.0 rad in both. One
measured gravity vector crosses the boundary, with a magnitude spread of at most
8.9e-15 m/s^2 between the source, canonical and world expressions of it. Supine,
the body's published ideal plane and the world's bed template coincide to
1.7e-16 m. Upright, the body publishes no plane, so the twelve loaded native
contact spheres are measured instead: all push along the source floor normal to
within 1e-6 rad and are coplanar to 2.0 mm. At t=0 they sit 111.9 mm below that
plane carrying 38952 N, about 51 body weights -- an initial-pose contact
transient of the native model, not a frame disagreement. The frame map is a
proper rigid transform: determinant 1 to 1.8e-15, point round trip 2.2e-16 m,
pair distances preserved to 4.4e-16 m, and force-port power preserved to
7.1e-15 W over 1193 skin samples.

## Rigid ball contact

`scripts/verify_ball_dynamics.py` broadens the single one-meter drop into 18
checks; the receipt is `data/derived/ball-dynamics-v1/report.json`. Every
"analytic" comparison there is a closed form of the *declared* impulse model --
Newton restitution plus a clamped Coulomb tangential impulse on a uniform solid
sphere -- so agreement shows the solver implements what it declares. It is not a
measurement against a real ball; restitution 0.75 and friction 0.35 remain
uncalibrated engineering values.

- Normal impacts at drops of 0.05 to 5 m return a rebound speed ratio of exactly
  0.75 (to 1e-16), with the contact impulse and dissipation matching
  (1+e)m|v| and the closed-form energy loss to 12 decimal places.
- Oblique impacts match the closed form in both friction regimes, sticking and
  sliding, to 5.6e-17 m/s, including the post-impact spin. A sticking impact
  leaves exactly 5/7 of the tangential speed and zero contact-point slip.
- A sphere launched sliding at 4 m/s reaches rolling at 2.857142857 m/s versus
  the classical 5/7 v0 = 2.857142857 m/s, at 0.333 s versus the analytic
  0.33285 s, with contact slip decaying at the analytic 7/2 mu g = 12.017 m/s^2.
- Over 400 randomized impacts the contact never adds kinetic energy, the
  reported impulse equals the momentum change to 1e-12, and the tangential
  impulse stays inside the friction cone.
- Supine (axis 2) and upright (axis 1) produce bit-identical bounce histories,
  so the ball obeys the same shared frame the body does.

Two limits are now measured rather than assumed, and neither is fixed:

- Furniture is a 1500 N/m penalty collider sampled at six probe points on the
  sphere, far softer than the analytic ground plane. A 0.4 kg ball arriving at
  3.43 m/s sinks 47.3 mm into a table before being pushed back out (15.2 mm at
  0.99 m/s). It never tunnels through, and the depth follows the linear-spring
  scaling, but it is visibly compliant, not rigid.
- Total angular momentum is conserved only to first order in the environment
  substep. An off-axis prop-against-prop strike leaves a spurious net angular
  momentum of 2.54% of the angular momentum exchanged at the production 1 ms
  substep, halving to 1.08% at 0.5 ms and 0.26% at 0.125 ms. Linear momentum is
  exact to 2.3e-16 kg m/s.

Ball against the real 98-muscle body, against cloth, and against the pillow
lattice have no closed form and are not validated here. Ball-vs-skin is checked
only for equal and opposite impulse exchange, substep by substep, to 1e-12. The current paired five-second cortex/cord/body/physiology/bedroom
comparison completes 250 exchanges per arm with ankle MSE 0.0295883 rad² full
versus 0.2615622 rad² severed (88.69% lower), establishing a narrow cortical ankle
contribution while also exposing the cloth defects above. Receipt:
`data/derived/unified-world-b_xo1l9j/report.json`.

Engineered stance also completed five seconds/250 exchanges with native
physiology. Receipt: `data/derived/unified-world-kk8i4684/report.json`; its observer
source attribution correction is recorded separately, with original and
corrected reports retained. The browser accepted actual upright stance frames.
Skin still uses approximate segment attachment; distributed rendering and
shared continuous surface attachment are being corrected and do not establish
soft-tissue FEM.

## Reproduce the unified comparison

Install the project's `brain` optional dependencies, then run from IHM-1:

```bash
PYTHONPATH=. .venv/bin/python scripts/verify_unified_world.py \
  --seconds 0.3 --controller implicit_cortical_ankle --arms full sever
```

Native initialization and execution are substantially slower than real time.
The script retains source/model identities, individual traces, and failures.
Longer trials are required before extending any short-run claim.
