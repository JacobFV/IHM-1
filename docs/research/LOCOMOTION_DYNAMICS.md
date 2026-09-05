# Native articulated muscle/contact dynamics

Implemented 2026-09-05. This is a native forward plant for a separately identified
source subject. It is useful for building the body's motor and mechanical
interfaces, but it is not a claim of controlled walking or a replacement of the
canonical anatomy by a single atlas.

## What executes

`ihm.native.locomotion.run_locomotion(LocomotionConfig(...), fresh_directory)`
assembles the retained 3D walking source into a fresh native model and advances
it with OpenSim `Manager` / Simbody Runge–Kutta–Merson integration. Positions,
speeds, muscle activation and compliant tendon/fiber states evolve under the
source mechanical equations. No coordinate is prescribed after initialization.

- 22 rigid segments, total mass **85.26984854173146 kg**.
- 80 source Millard2012 muscles, retaining tendon compliance, active and passive
  fiber force laws. Muscle paths use the source polynomial FunctionBasedPath
  coefficients. This is an explicit path approximation, not regenerated anatomy
  or the known unstable geometric wrapping variant.
- Twelve SmoothSphereHalfSpaceForce foot contacts and the source ground plane.
  Contact sphere centers receive the source example's +0.02 m Y preparation.
  Native stiffness is 1e6 Pa, dissipation 2 s/m, static/dynamic friction 0.8,
  viscous friction 0.5 and transition velocity 0.2 m/s, with native smoothing.
  Stiffness and Coulomb friction scales are exposed separately for experiments.
- The source passive coordinate force set is applied. Thirteen retained
  torso/arm coordinate actuators receive zero controls. There are no imposed
  measured ground reaction forces, corrective pelvis actuators, or coordinate
  target tracking in this runner.
- Initial independent q are taken once from the measured/RRA-adjusted reference
  at 0.48 s. Initial u use the centered adjacent reference rows (0.479/0.481 s).
  Dependent coordinate reference values are intentionally not imported; native
  constraint assembly and velocity projection supply them. Unreferenced
  coordinates retain source defaults. Muscle activation starts at the requested
  excitation and native muscle equilibration supplies fiber states. This is
  **an assumed initial muscle state**, not a measured or optimized gait state.

The optional bounded instantaneous positive velocity controller is

`u = clamp(u0 + gain * max(v_MT, 0)/(l_opt * v_max), 0.01, 1)`.

It is derived from the retained OpenSim demonstration controller's dimensional
normalization, with an explicit baseline and saturation. Gains and the uniform
baseline are engineering assumptions. It has no biological delay, gait phase
logic, balance, muscle-specific identification, IBM neural input or calibrated
human reflex interpretation. Its callback is stateless, so adaptive trial
evaluations cannot advance fictitious controller history.

## Outputs and provenance

Each fresh run retains `execution.json`, `engine.log`, `assembled_model.osim`,
`model.json`, `frames.jsonl` and `termination.json`. The execution receipt hashes
the runner, its source, native OpenSim/Simbody libraries and every source file
in the retained example directory. Stale executable/source/library combinations
are rejected before creating a run directory. The build manifest records pinned
engine source identities.

Frames contain all native state variables, generalized q/u, all 22 rigid
transforms into source Ground, COM position/velocity, native kinetic/potential
energy, every muscle's activation/excitation/length/fiber state/tendon force and
mechanical power, and all 12 simulated contact forces/moments/penetrations.
Lengths are m, time s, forces N, moments N m, energy J and power W. Generalized
coordinates use their native rad/m units; generalized speeds use rad/s or m/s.
State labels are exported in native order. Rotations are not deformation fields.
Contact moments are expressed in Ground about their attached body's origin.
Reported penetration is **virtual sphere penetration against the source Y=0
plane**, not a resolved plantar tissue strain or measured skin indentation.

The force audit computes `mass * (COM_acceleration - gravity) - sum(GRF)` directly
from native acceleration and independently sums paired sphere/floor forces.
It does not infer conservation from animation. Active muscle mechanical power
is `-tendon_force * musculotendon_lengthening_speed`; it is not chemical energy
consumption and must not be copied directly into BioGears metabolic demand.
Native energy readouts alone do not establish a closed chemical/thermal ledger.

## Verification actually executed

Build and reproduce from repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_native_opensim_dynamics.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_locomotion.py
```

Acceptance receipt: `data/derived/locomotion-6lm0neja/verification.json`.
Nine independent native runs passed: baseline, changed initial pelvis speed,
contact stiffness, friction, excitation, feedback, two integration refinements,
and a longer fall trial. The primary comparison is 0.1 s at 200 Hz. These are
numerical/causal checks, **not experimental validation**.

| Check | Observed result |
|---|---:|
| Initial pelvis speed +0.05 m/s, final q difference norm | 0.0200059 |
| Contact stiffness ×0.5, final q difference norm | 0.138147 |
| Coulomb friction ×0.25, final q difference norm | 0.158927 |
| Uniform excitation 0.03 → 0.1, final q difference norm | 0.817041 |
| Velocity gain 0 → 0.5, final q difference norm | 0.572422 |
| Step cap 1 → 0.5 ms, accuracy 1e-6 → 1e-7, q difference | 1.19131e-5 |
| Step cap 0.5 → 0.25 ms, accuracy 1e-7 → 1e-8, q difference | 1.27842e-7 |
| Baseline maximum instantaneous momentum residual | 9.10195e-13 N |
| Maximum paired contact force residual | 0 N |
| Maximum native position / velocity constraint error | 1.11022e-16 / 1.58281e-15 |
| Baseline vertical GRF range | 297.482–1191.913 N |
| Baseline maximum virtual sphere penetration | 0.0203319 m |
| Feedback trial excitation range | 0.03–0.400395 |

The q difference norm mixes native rotational and translational entries. It is
an implementation sensitivity/convergence diagnostic, not a physical length,
clinical error, or statistical anatomical confidence interval. The large
virtual contact penetration and assumed muscle initialization need contact
identification/initialization work before any gait interpretation.

Tests were written and observed to fail before implementation: missing native
API, then missing force audit, then missing explicit termination receipt. The
completed acceptance also checks finite input rejection, frame counts/clocks,
actual q motion, all muscles' activation bounds, body/muscle/contact counts,
absence of prescribed coordinates, constraint residuals and causal differences.

## Fall and coverage boundary

The unbalanced one-second request leaves the upright plant's admissible region
at **0.415 s** and returns `completed: false` with
`reason: upright_plant_envelope_exceeded`. A COM-height guard of 0.55 m is an
explicit engineering restriction because this model has only foot collision
coverage. It is not a physiological fall threshold. A diagnostic longer trial
before this guard continued falling through the floor; those original bytes are
retained under `data/derived/locomotion-long-_ov6g7kd/` and its owned process was
terminated. A finite trajectory with missing body-ground contact is not success.

Remaining integration work is concrete: identify a model-compatible walking
policy and initial muscle state, add whole-body contact, characterize the
fitted-path validity envelope, run sustained and perturbed gait tests, and pass
physical wrenches to dynamic tissue domains with equal/opposite reactions.
This runner currently has **no deformable organ/tissue dynamics**, no canonical
attachment transfer, no online physiological coupling, and no neural controller
from IBM. It does not display its measured subject as the canonical generic body.

The canonical body is 77.1107029 kg and uses heterogeneous anatomical evidence.
Its existing Rajagopal registration was fitted to a different 75.337 kg model;
neither registration nor a rigid transform may silently be reused here. A next
unification should allocate each segment's mass across skeleton and deformable
domains once, preserve per-source attachment transforms and uncertainty, and
verify virtual work under wrench transfer. Additional muscles/tendons from other
sources or synthesis can participate in that shared model without requiring one
complete donor atlas. This native subject is a verified mechanical fixture for
that integration, not its anatomical scope limit.

## Primary retained sources

See [the pinned OpenSim 3D walking example](https://github.com/opensim-org/opensim-core/tree/86b30588374650fbaf012a345a836a64f6855522/OpenSim/Examples/Moco/example3DWalking)
for the exact model, contact, coordinate and path-function sources. Its README
links the original [full-body model](https://simtk.org/projects/full_body),
[high-flexion model](https://simtk.org/projects/model-high-flex) and
[passive calibration model](https://simtk.org/projects/fbmodpassivecal).
The broader source and registration audit is [LOCOMOTION_CONTROL_AUDIT.md](LOCOMOTION_CONTROL_AUDIT.md).
