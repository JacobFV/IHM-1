# Physical environments

Environment selections now cross the live body API boundary. `environment_selection`
contains catalogue scene/object/component IDs, which the server validates against the
base environment and dependency rules. Browser geometry is not a physics owner.

The server instantiates fixed compound furniture/wall colliders, spring cloth,
a three-dimensional pillow spring lattice, and six-degree-of-freedom movable balls
and blocks. It samples the canonical skin envelope at 45 mm voxel spacing and
attaches the samples to the existing named native bone owners. There is no second
human mass. Contact forces act at these material stations through
`ArticulatedBodyPlant.advance`, including the force moment about the native segment.
The same forces enter the embodied runtime's existing respiratory load projection;
mechanical observations continue into its neural and metabolic exchange.

Environment integration uses 1 ms substeps and a 5 ms explicit exchange with the
body. Cloth/body, prop/body, cloth/prop and prop/prop interactions exchange reciprocal
impulses. Fixed furniture and the anchored pillow foundation supply external support.
Environment checkpoints are restored alongside body checkpoints on a pre-native
failure. An uncertain native commit retains the existing fail-stop behavior.

The accepted `environment_state` reports the clock, force ports, impulse, meshes,
rigid transforms, selection, coefficients and asset hashes. The browser consumes
these positions; it performs no independent cloth or rigid-body integration.
Factory outputs retain source and selected asset copies. Mattress material and skin
foundation choices reach the existing native factory; mattress support is not
counted again by this environment solver.

## Appearance

Paint relief, wood grain, textile weave, stone, rubber and metal use deterministic
albedo/bump textures and physical materials. Furniture has beveled edges; trees use
instanced leaves within the declared canopy envelopes. Shadows follow the gravity
frame. Camera framing also follows that frame, including supine rooms. Cutaway walls
remain colliders. Bare native floors are drawn at their declared level, without a
second floor overlay on scene ground. The blanket is now a loose supine top cover
without pinned tucks. It drops for 1.5 s against held skin before time zero. A
sampled upper-skin envelope keeps the cover outside the body; mass-weighted
unilateral constraints limit stretch/shear edges to 112% of rest length.
Alternating contact and stretch projections report their actual residual and
whether they converged. Passive contact impulses do not turn positional
correction into cloth velocity. Preparation resets residual velocity as an explicit
initial condition and reports the discarded kinetic energy. Position projections
can change spring and gravitational energy; this nonconservative work is reported.
This is not equilibrium certification or a model of arbitrary wrapped fabric. Initial pillow geometry has a rounded surface and a volumetric lattice.

## Verified and unverified

- Twelve causal solver/rollback tests pass, including native force-port delivery through
  an embodied runtime fixture and reciprocal prop/cloth impulses.
- JavaScript tests verify server mesh updates, stable duplicate instance identities,
  rigid rotation about the correct center, and textured material families.
- All six catalogue scenes complete 3 s against a stationary canonical body fixture.
  This tests environment integration, not whole-body feedback stability.
- A real native OpenSim comparison with and without bedroom contacts changes body
  position by approximately 32.9 mm after the first 20 ms. The 0.2 s native run
  finishes with 82 contact force ports. Retained result:
  `data/derived/environment-native-1snh2tfi/report.json`. This demonstrates a causal
  load path, not realistic transient magnitude or calibrated biomechanical behavior.
- All six scene previews were inspected in Chrome. Screenshots and browser errors
  are retained under `test-results/environment-upgrade/`.
- Earlier full-runtime attempts failed signed-muscle budget checks. Subsequent
  runtime work is tracked in [BRAIN_BODY_WORLD.md](BRAIN_BODY_WORLD.md). Historical
  measurements here predate the current common-frame and loose-cover changes.


Contact stiffnesses, damping, friction and pillow/textile coefficients are engineering
values, not measured or clinically validated properties. The body surface attachment
is approximate. Fixed cylinders and canopies use bounding boxes; rigid manifolds use
sampled vertices, so edge-edge and high-speed collisions are not comprehensive.
Cloth now has approximate vertex self-contact at 12 mm separation, excluding nearby
material neighbors. It preserves free-pair momentum and dissipates kinetic energy
at impacts; position projection can alter spring potential energy. This is not
triangle/edge collision detection, so face intersections and fast tunneling remain
possible. There is no resolved cutaneous pressure, fabric thermal insulation, or
separate leaf/blade mechanics. The native support law remains the
owner of bed/floor support. Live temperature changes are unavailable in this adapter
and are disabled in the UI rather than silently accepted.

The cloth implementation is a damped spring model, **not XPBD**. The
[XPBD paper](https://mmacklin.com/xpbd.pdf) is a useful reference for a future compliant
constraint solver with more robust stiffness control; this implementation does not
claim its iteration/timestep properties.

## Reproduce

```sh
PYTHONPATH=. .venv/bin/python scripts/verify_environment_dynamics.py
PYTHONPATH=. .venv/bin/python scripts/verify_environment_scenes.py
PYTHONPATH=. .venv/bin/python scripts/verify_environment_native.py
PYTHONPATH=. .venv/bin/python scripts/verify_environment_embodied.py
npm --prefix app test
npm --prefix app run build
node app/test/environment-preview.mjs
```

The native and embodied scripts create fresh retained run directories and close
only their own owners. The preview script starts and closes its own local server on
port 8766. The commands exercise current source; historical reports retain their original scope.

## Current common-frame mechanics proof

Native body/world coupling now exchanges every 5 ms within the 20 ms neural and
physiology tick. World physics uses 1 ms substeps. The native body owns its support;
a separate bounded mattress collider supports cloth and props without duplicating
body support. Outside the bed, objects fall to the actual room floor.

Physics runs in a gravity-aligned world frame. Native material skin stations enter
through an explicit rigid transform; public force inputs and returned body reactions
are canonical. `environment_state.world_frame` retains both 4×4 transforms,
gravity vectors, and native support plane. Object geometry remains in
`gravity_aligned_world`. Rendering must transform the whole surround, including
fixed furniture, and convert picked targets back to canonical coordinates. The
transform preserves the canonical tangential origin while aligning support planes;
it changes no anatomy, mass, or registration.

A current 98-muscle native body + bedroom + ball run completed 0.1 s at 5 ms
exchange: canonical gravity agreed within 1e-12 and support plane discrepancy was
1.67e-16 m. A ball placed 1 cm above the actual room floor rebounded between 40 and
60 ms. Canonical applied-force and reaction-impulse checks passed. Blanket maximum
edge extension was 1.1201165 at 20 ms (reported not converged), then 1.1200794 at
100 ms (within tolerance); maximum cloth speed at 100 ms was 2.04 m/s.
This short causal mechanics proof does not establish equilibrium, calibrated fabric,
long-run stability, brain control, or walking. Report:
`data/derived/common-world-native-907cp3gk/report.json`.

```sh
PYTHONPATH=. .venv/bin/python scripts/verify_world_frame.py
PYTHONPATH=. .venv/bin/python scripts/verify_cloth_cover.py
PYTHONPATH=. .venv/bin/python scripts/verify_cloth_stretch.py
PYTHONPATH=. .venv/bin/python scripts/verify_common_world_native.py
```

The historical reports below predate this frame correction and loose-cover model.
They remain causal diagnostics, not validation of the current complete scene.

## Bouncy ball and cloth update (2026-09-08)

Both the interactive sphere and server environment ball now use the same analytic
plane time-of-impact solver. The uncalibrated rubber material has restitution
0.75 and friction 0.35. A 0.05 m/s low-speed threshold transitions to resting
contact. Normal and tangential impulses cannot add kinetic energy; any initial
penetration repair is reported separately and does not become velocity. Sphere
collisions with furniture, skin and other props retain the compliant contact law.

The environment drop test runs 10 simulated seconds, measures rebound clearance
against the theoretical restitution-squared ratio (0.5625), verifies final rest,
and bounds the mechanical energy plus dissipation residual below 1e-8 J. The
interactive scene check separately verifies its rebound and energy ledger. Five
additional tests exercise cloth self-contact, anchors, replay and momentum. These
are component and stationary-body fixture results, not evidence of walking or
brain control.

```sh
PYTHONPATH=. .venv/bin/python scripts/verify_cloth_contact.py
PYTHONPATH=. .venv/bin/python scripts/verify_interactive_scene.py
```

Native scene objects accept the same public `{id, force_n, point_m}` force ports
as the body. The runtime partitions ownership before mutation. Environment ports
are limited to 100 N each and 32 entries; stale/outside material points and fixed
cloth nodes reject. Rigid points bind in local coordinates and receive torque;
soft points bind to their nearest free vertex. These targets follow the material
through substeps. Applied object ports are reported in the environment frame and
roll back with its checkpoint. `fixed_nodes` lets the browser avoid anchored cloth.

An actual native OpenSim bedroom-plus-ball run completed 0.12 s with synchronized
clocks, 109–170 contact ports, and a beside-bed ball reversing from -0.3924 m/s
at 40 ms to +0.18655 m/s at 60 ms. Its initial floor clearance was explicitly
set to 1 cm. The blanket moved up to 0.3035 m during this short coupled transient;
this is evidence of coupling, not calibrated draping realism or long-run stability.
Retained result: `data/derived/native-ball-blanket-6gcn_wrx/report.json`.

## Moving skin barrier correction

The previous skin barrier converted penetration depth to velocity using the
1 ms environment substep. A 10 mm body move over a 20 ms exchange therefore
created 10 m/s cloth speed from a 0.5 m/s moving boundary. The replacement
projects position separately and removes only inward relative normal velocity,
returning the equal/opposite velocity impulse. Six focused tests check this
amplification regression, static passivity, momentum, anchors and diagnostics.
Spring-energy change caused by positional repair is explicitly reported; this
is still a discrete approximate contact model, not an energy-conserving surface
solver. Prescribed-boundary work is an estimate, not measured native work.

Diagnosis found canonical and initial native skin agree within 1.31e-16 m. The
large early transient includes native posture settling, amplified by the old
barrier. A runtime-patched passive comparison completed 2 s of actual native
body/bedroom coupling; final cloth speed was 0.907 m/s and body skin movement
was 0.010 m over the last 20 ms. The original-barrier comparison was interrupted
at 0.92 s after 276 s wall to free compute; this is not a solver-failure claim.
Reports: `data/derived/native-blanket-passive-f7q_e_7f/report.json` and
`data/derived/native-blanket-diagnose-9dvbf1u0/report.json`.

That historical implementation (including its corrected cloth preparation) also
completed a fresh 2 s actual native/bedroom run, in 57.95 s wall. Final cloth
speed was 1.636 m/s, maximum displacement 0.363 m, and native skin movement
over the final 20 ms was 0.00891 m. The body still underwent a large posture
transient; completion is not equilibrium or biomechanical calibration. Final
source result: `data/derived/native-blanket-final-emmvub3y/report.json`.

A same-initial-checkpoint 0.4 s comparison exposed strong mechanical exchange
timestep sensitivity. Reducing body/world exchange from 20 ms to 5 ms reduced
peak sampled-skin speed from 32.60 to 7.65 m/s and final skin displacement from
1.268 to 0.396 m. Native external work changed from +20.59 J to -2.50 J; blanket
displacement fell from 0.339 to 0.228 m. Peak blanket speed remained high
(12.34 versus 10.64 m/s), and its final speed did not decrease monotonically.
These measurements support reducing partitioned contact lag; they are not a
convergence or anatomical accuracy certificate for 5 ms. Reproduce with
`scripts/verify_environment_exchange.py`. Retained report:
`data/derived/native-environment-dt-kq0sbfhi/report.json`.

## Sustained current-frame result and remaining cloth defect

The current trained 128-site cortex + cord + native body + physiology + bedroom
and ball completed 5 s (250 synchronized neural ticks) in
`data/derived/unified-world-b_xo1l9j/full/report.json`. The ball rebounded against
the bounded mattress and settled. The blanket remained bounded, but maximum
sampled extension reached 1.18637 despite its configured 1.12 target; 41 of 250
sampled final substeps correctly reported nonconvergence. Its final maximum speed
was 2.15 m/s, so this is not settled fabric.

A held-body diagnostic confirms a remaining numerical defect: over 0.5 s the
current positional constraints supplied +11.06 J of unresolved work and total
cloth mechanical energy increased 0.94 J despite stationary boundaries.
Disabling self-contact did not remove the energy increase. Evidence and the
reproducer are retained as `cloth_jitter_diagnosis.json` and `.py` in that same
run directory. The contact/material projection is therefore not a passive cloth
solver. A replacement is being tested separately; successful runtime completion
must not be presented as calibrated or energy-consistent blanket physics.

These measurements predate the fold and energy repair recorded at the end of
this document, which reruns the same diagnostic fixture. The 1.724 extension and
the +0.94 J energy growth are the BEFORE column of that comparison.

The matched severed-kernel arm also completed 5 s. Ankle tracking MSE was
0.0295883 rad² for the trained cortex and 0.2615622 rad² when severed (88.69%
reduction in this one target/posture/duration). This is narrow cortical task
contribution, not walking. The severed arm's blanket reached 1.7240 maximum
sampled extension and 8.37 m/s maximum speed, with 78 of 250 nonconvergence
flags. This reinforces the cloth limitation; survival alone is insufficient.
Full matched report: `data/derived/unified-world-b_xo1l9j/report.json`.

## Blanket fold and energy repair (2026-09-08)

The sharp blanket folds were reproduced in the simulated mesh, diagnosed, and
traced to two separate defects. Neither was a rendering artifact and neither was
fixed by changing appearance.

### A crease was a zero-energy deformation

The sheet carried only distance springs, so it was a membrane, not a shell.
Rotating half the 437-node sheet about a grid line is a rigid motion of that
half: no material edge changes length and the material spring energy is exactly
zero. A fold therefore cost nothing at all. The only pre-existing resistance to a
crease about a grid line came from the two-apart springs, whose chord shortens as
the square of the fold angle, so their energy falls as its fourth power. Halving a 0.4 rad crease divides their
energy by 16 and a bending energy by 4; both ratios are asserted in
`verify_cloth_bending.py`. A smooth one-centimetre sag across the 0.84 m sheet
cost 9 microjoules of elastic energy against 0.147 J of gravitational energy, so
nothing in the model set the drape shape at all.

`ihm/assembly/cloth_bending.py` adds the missing hinge term: the quadratic
flat-rest bending model of Bergou, Wardetzky, Harmon, Zorin and Grinspun (SGP
2006), whose per-hinge energy agrees with the discrete-shells law
`3 B |edge|^2 theta^2 / (A0 + A1)` to leading order. The blanket has 1148
hinges. Its Hessian is constant and positive semidefinite, so the implicit
solver takes it unchanged and it cannot create energy. The rigidity is
2.5e-3 N*m, chosen so the fabric bending length `(B/w)^(1/3)` is 49.7 mm, at
least the 46.7 mm vertex spacing: a smaller value leaves creases the mesh cannot
resolve. It is an engineering value in the heavy-fabric range, **not** a measured
Kawabata bending modulus, and the quadratic model understates the restoring
torque at large fold angles.

### Contact repair displacement was inflated by 1/normal_z

`SupineClothCover.query` returns a VERTICAL penetration depth against the
sampled upper-sphere height field, and `project` used to repair that depth
vertically. The minimum-norm repair of `z >= height(x, y)` is along the analytic
surface normal by `depth * normal_z`; a vertical repair is `1/normal_z` times
larger, and `normal_z` goes to zero at every sample sphere's rim. Vertices near
the body silhouette were therefore displaced far more than their geometric
overlap, and the springs turned that displacement into motion. This is the
positional analogue of the moving-skin-barrier defect corrected earlier: the
correction was never converted into velocity, but its magnitude was wrong.

The envelope itself is genuinely rough at the cloth's own wavelength. Measured on
the held canonical body, its height field has second differences up to 0.0442 m
over a 1 mm sampling, that is slope discontinuities of order 40, at a cusp
spacing (45 mm skin samples, 48 mm spheres) matching the 46.7 mm cloth vertex
spacing. `cloth_cover.py` now applies the minimum-norm normal-direction repair,
so a cusp no longer throws a vertex. The envelope remains a sampled sphere union
with seams, not a continuous anatomical skin.

### Overlapping contact constraints no longer break the implicit solve

Reproduced: two vertices whose positions are already fully determined (by anchors
or by three independent contact tangents) with a material edge between them at
its extension ceiling. The active constraint Jacobian is rank deficient, the KKT
matrix is exactly singular, `spsolve` fails, and the step was rejected with the
advice to retry a smaller timestep. A smaller timestep cannot repair a rank
deficiency, so the step could never succeed and `ClothPassiveGeometry.advance`
bisected to its minimum and abandoned the whole interval.

`ClothPassiveStepper._kkt_solve` now solves such a system by dual regularization
refined against the ORIGINAL rows. Redundant rows are not dropped by inspection,
because choosing which of several physically real contacts is the duplicate is
arbitrary. Regularization only picks one member of a non-unique multiplier set;
the total reaction is unchanged, which `verify_cloth_passive_step.py` asserts by
comparing the summed reaction for one, two and five duplicated planes. Every
acceptance test still runs against the unregularized constraints, so a genuinely
inconsistent overlapping set is still rejected atomically with the state
unchanged.

### Triangle-level and continuous self-collision

`ihm/assembly/cloth_self_collision.py` adds the two primitive pairs a sheet
needs: vertex against triangle and edge against edge, as a discrete proximity
test at a declared 6 mm sheet thickness and as a straight-line swept test over
the step. The previous particle guard only kept vertices 12 mm apart, so two
coarse faces could pass through each other between their vertices and a fast
vertex could cross a face without ever approaching one of its corners.

Every response acts along the closest-point separation with barycentric weights
that sum to zero, so linear momentum is conserved exactly and the impulse acts at
one coincident material point, giving exactly zero net torque. The normal impulse
is perfectly inelastic and can only remove kinetic energy. Position repair is
separate and is never converted into velocity. A swept hit is repaired back to
the side the primitives started on, so a tunnelled vertex does not merely stop.

The pass is skipped only under an exact clearance argument: the distance between
any two primitives can shrink by at most twice the largest vertex displacement
since the last evaluation, so while twice the travel stays under the measured
clearance margin no contact is possible. Folding the sheet shrinks that margin
and the pass resumes every substep. Remaining limits: no friction, straight-line
vertex paths within a step, and a capped Gauss-Seidel sweep count whose residual
separation is reported rather than assumed.

### Measured result

Held canonical body, fixed supports, 500 one-millisecond substeps of the
production path. The before column is the retained pre-existing diagnosis
`data/derived/unified-world-b_xo1l9j/cloth_jitter_diagnosis.json` on the same
fixture.

| Held-body 0.5 s | Before | After |
|---|---|---|
| Mechanical energy growth rate (fitted trend) | +0.9761 +- 0.1073 W | +0.00512 +- 0.00106 W |
| Largest transient energy rise above the start | +1.4372 J | +0.00769 J |
| Net mechanical energy change | +0.9355 J | -0.00477 J |
| Elastic energy injected by position projection | +17.942 J | +2.376 J |
| Unresolved projection work | +11.064 J | +0.873 J |
| Peak speed, mean over the last 100 substeps | 3.423 m/s | 1.115 m/s |
| Maximum sampled edge extension | 1.12015 (1.724 in the 5 s runtime) | 1.12010 |
| Substeps reporting unconverged strain | 41 and 78 of 250 in the 5 s runtime | 0 of 500 |
| Maximum hinge fold angle after preparation | 174.5 deg | 35.9 deg |
| Hinge folds over 60 deg after preparation | 85 | 0 |
| Mean hinge fold angle after preparation | 28.0 deg | 7.9 deg |
| Preparation unresolved projection work | +26.96 J | +1.84 J |
| Preparation discarded kinetic energy | 0.2470 J | 0.00218 J |

Exactly satisfied in the same run: cloth/body contact impulse reciprocity to
1.04e-17 N*s, self-collision net force 0 N*s, net torque 0 N*m*s, self-collision
kinetic energy change 0 J, no cloth vertex inside the skin envelope (0 m), and a
minimum non-adjacent primitive separation of 0.0226 m against the 0.006 m
declared thickness. Receipt: `data/derived/cloth-energy-audit-20260908/receipt.json`.

The same experiment on the separately verified implicit path (437-node sheet,
actual sampled canonical skin, 1.0 s, `verify_cloth_passive_skin.py`) gives, with
and without the hinge term: maximum fold angle 162.5 deg against 42.3 deg, hinge
folds over 60 deg 31 against 0, mean fold angle 21.2 deg against 6.7 deg. Every
5 ms interval of that run is accepted against its own audit: worst energy excess
-0.00091 J (always dissipating), worst momentum residual 1.2e-14 N*s, worst
constraint violation 1.1e-15 m, maximum extension exactly 1.12. Retained
reports: `data/derived/passive-cloth-skin-s45eaf2g/report.json` without the
hinge term and `data/derived/passive-cloth-skin-yzwjp6iy/report.json` with it.

Reproduce:

```sh
PYTHONPATH=.:scripts .venv/bin/python -m unittest \
  scripts.verify_cloth_bending scripts.verify_cloth_self_collision \
  scripts.verify_cloth_contact scripts.verify_cloth_stretch scripts.verify_cloth_cover \
  scripts.verify_cloth_passive_geometry scripts.verify_cloth_passive_step
PYTHONPATH=. .venv/bin/python scripts/verify_cloth_energy_audit.py --seconds 0.5
PYTHONPATH=. .venv/bin/python scripts/verify_cloth_passive_skin.py --seconds 1.0
```

`verify_cloth_energy_audit.py` writes its receipt into a fresh temporary
directory under `data/derived`; never hand it a source root.

### What is still not true

- The blanket is still a triangulated spring-and-hinge shell. It is **not**
  volumetric soft-body FEM and not a validated textile model. Stiffness,
  bending rigidity, thickness and friction remain uncalibrated engineering values.
- The production explicit path is **not** a passive solver. Its fitted energy
  growth rate is still significantly positive, +0.00512 W against a standard
  error of 0.00106 W, and single substeps still rise up to 7.7 mJ. The receipt
  asserts only an explicit regression bound, one percent of the retained
  baseline growth rate, and records the remaining non-passivity under
  `not_yet_met` rather than moving a threshold to make passivity look true. The
  implicit path is the one with a per-step energy and momentum acceptance test.
- Twelve alternating passes per substep are needed for the strain ceiling to
  converge on every substep, against six before; at six passes 59 of 200 held-body
  substeps still reported nonconvergence at 1.12039. The blanket costs roughly
  twice the wall time it did; the correctness is not free. All six catalogue
  scenes still complete their 3 s fixed-body fixture, the bedroom in 182 s wall.
- Self-collision has no friction, assumes straight-line vertex motion within a
  step, and caps its Gauss-Seidel sweeps; a dense pile can end a step with a
  reported residual overlap.
- The skin contact geometry is still a sampled union of 48 mm spheres on 45 mm
  samples. Its seams are real geometry in this model, not a continuous body
  surface, and the cloth mesh only just resolves them.
- Cloth preparation remains an initial condition, not a certified equilibrium.
