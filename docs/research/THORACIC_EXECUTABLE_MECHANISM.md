# Executable thoracic material mechanism

`ihm/assembly/thoracic_mechanism.py` loads
`data/research/thoracic_mechanism/v2/manifest.json` and evaluates source material
positions, Jacobians, mass properties, generalized kinetic energy, a geometric
cavity volume, and work-conjugate loads. It is an executable Python static/
kinematic mechanism, not a native dynamic integrator or a physiological muscle
model. No recoil stiffness, damping, activation, or chemical work is added.

The source anatomy remains the [registered thoracic recipe](THORACIC_ANATOMY_RECIPE.md).
All 48 debited moving material shares are represented. The residual core remains
17.937704388 kg; together they recover the cervical-reduced 20.235757397 kg
subsystem. This prospective composition does not modify current native body
ownership or account for other independently proposed mass partitions.

## Moving attachment maps

Rib source vertices follow their proposed hinge rigid transforms. The left and
right twelfth-rib axes are ambiguous in the source recipe; their coordinates
11 and 23 are constrained to zero, and their mass moves only with the parent.
The sternum's three pieces share one anterior translation coordinate.

Cartilage and intercostal material surfaces use their retained source-node
attachment candidates. Diaphragm material combines moving attachments with its
previously retained inferior-descent field. For each material node, the three
nearest reference anchors supply normalized inverse-square displacement weights.
At an anchor the weight is exactly one, so its material node follows its chosen
target owner. The source-side anchor point, including its offset from the target
surface, is transformed by that owner; the original gap rotates without being
silently snapped closed. These are cardinal geometric interpolation maps, not
calibrated collagen fibers, central-tendon recruitment or tissue strain laws.

Nine unsupported intercostal pairings are excluded from attachment choices.
Their material mass is still included and moved by the remaining anchors.
Repeated source nodes retain the closest-gap target and record the discarded
alternatives: 51 repeated candidates are retained in the audit. This avoids
pretending one material node can satisfy conflicting independent attachments.
Some material points lie as far as 150.5 mm from their nearest selected anchor;
that extrapolation remains an explicit fidelity limitation.

The material interpolation is fixed in reference coordinates. It is not
recomputed after movement, which would introduce unrecorded changes in the
Jacobian and work. Positions inside each source triangle use linear nodal
interpolation. A complete native model must further validate attachment
constraints, tissue strains and contact before using these maps dynamically.

## Complete material kinetic metric

The generalized velocity has 32 entries: three parent linear velocities,
three parent angular velocities, 24 rib angular velocities, sternum translation
rate and diaphragm descent rate. With two constrained rib coordinates there
are **30 independent generalized velocities**.

For a material point in instantaneous torso axes,

```
v_point = v_parent + omega_parent × x + J_internal(q) qdot
A = [I, -skew(x), J_internal]
M(q) = integral rho A.T A dA + residual-core rigid contribution
T = 0.5 * velocity.T M velocity
```

Positive three-point degree-two quadrature integrates each triangular lamina.
Each triangle keeps its reference mass; deformation does not renormalize mass
by current area. This quadrature integrates first/second moments and kinetic
energy exactly for the chosen linear triangle map. The residual core adds its
COM-offset translational contribution and its full physical COM inertia.
The diaphragm's prior effective mass is **not** added as an independent body;
its kinetic coefficient and parent/rib cross terms arise from the same material
integration. Cartilage and intercostal kinetic energy are also retained.

The active mass matrix is positive definite in reference and tested perturbed
poses, with minimum eigenvalues approximately 6.76486e-5 in the declared mixed
units. This is a positivity/rank check, not a unit-invariant conditioning metric.
No stiffness is used to remove degeneracy. The inactive constrained rows are
not presented as independent degrees of freedom.

Reference recombination differs from the bound torso ledger by 1.78e-14 kg,
6.60e-17 m COM, and 1.00e-15 kg m² inertia. At the retained perturbed pose, direct
quadrature kinetic energy agrees with `0.5*v.T*M*v` to 1.39e-17 J. The parent
translation rows of M also contain the total first-moment Jacobian; an applied
body-frame gravity vector g gives generalized gravity force `M[:3,:].T @ g`.
No default gravity direction or supine equilibrium is inferred here.

## Interpretable geometric cavity candidate

Five retained, registered Z-Anatomy lung-lobe surfaces define a closed oriented
**convex envelope**, with 3,422 vertices and reference volume 0.006526046 m³.
This 6.526 L envelope fills interlobar/mediastinal gaps and concavities. It is
not native gas-exchange volume, physiological lung capacity or a new gas store.
The registered source snapshot has Z-Anatomy identity fields rather than an
explicit display-decimation field; its bytes and original geometry hashes are
preserved without upgrading that schema into an unsupported fidelity claim.
Z-Anatomy attribution and CC BY-SA 4.0 terms accompany the lobe-derived data;
BodyParts3D attribution remains retained upstream.

Each envelope vertex carries a fixed three-nearest-material-node displacement
map using ribs, sternum and diaphragm from the same mechanism. Maximum nearest
node gap is 23.24 mm. This is an engineered attachment between geometric
representations, not identified pleural contact. The same source-node motion
and Jacobian drive both material evaluation and cavity shape.

For the closed triangle surface, V is the signed tetrahedral sum and its nodal
gradient is evaluated analytically. `J_V = sum grad_x(V) * J_material` follows
by the chain rule. The API can apply a **caller-supplied** native reference
volume as a constant offset:

`V_candidate(q) = V_native_reference + V_geom(q) - V_geom(0)`.

The retained probe supplies no native gas reference. A separate test uses an
explicit illustrative value to verify the offset algebra. Neither value
calibrates native gas volume or licenses a second volume state.

The geometry exploration bounds are ±0.05 rad on unlocked ribs, ±5 mm sternum
translation and ±20 mm diaphragm descent. They are declared engineering bounds,
not physiological ROM. Cavity topology must stay closed and oriented, volume
positive and facets locally unflipped. The retained nonzero probe also checks
all material triangles and finds no reversals. These sampled/local checks do
not prove global absence of self-intersections throughout the entire box.

## Pressure work and parent reaction

Positive input p is **outward mechanical transmural pressure**:

`Q_internal = p * J_V`, so `Q_internal·qdot = p*Vdot`.

The probe compares this power with independent nodal pressure-force work and
finite-difference cavity-volume rates. Across all 24 active internal
coordinates, maximum analytic/finite-difference Jacobian difference is
1.09e-11; pressure-work difference is 5.07e-11 W for a 75 Pa test load.
Closed-envelope uniform pressure produces zero resultant force and moment to
floating-point accuracy. It must not create an unbalanced parent push.

`point_load` returns both the full parent/internal generalized load and the
opposite resultant needed **if the parent is held fixed**. Tests verify force,
moment and point/generalized power identities. That reaction is not an extra
cancelling force to apply to a freely moving parent. Caller forces and twists
must be expressed in the declared torso frame or transformed consistently;
there is no independent global anatomical registration in this module.

## Native compliance and work ownership

The probe retains exact held-source hashes and matching lines from native
`BioGears.cpp` and `Respiratory.cpp`. The two
`PleuralCavityToRespiratoryMuscle` branches are existing chest compliances;
`EnvironmentToRespiratoryMuscle` is the existing ideal pressure driver.
Separate lung recoil also remains native-owned. None is replaced by this
kinematic mechanism, and it adds no overlapping elastic/dissipative law.

The signed worker's [native pressure-source observer](NATIVE_RESPIRATORY_WORK_PORT.md)
uses source pressure difference P and solved source flow Q with work P Q dt.
Our outward p and geometric Vdot are a distinct mechanical convention. A later
connection must establish effort/flow signs and cavity ownership explicitly;
it cannot assume total lung-volume change always equals native source stroke.
No observed ideal-source work is credited to muscle ATP or thermal metabolism.

Before dynamic coupling: choose and validate recoil/actuator ownership,
explicitly transfer overlapping native chest compliance and ideal drive,
retain airway/gas/lung ownership as intended, derive the inertial equations
from the full M(q), and perform synchronized pressure-volume/contact acceptance.
Repeated irreversible native physiology steps cannot serve as a mechanical
fixed-point rollback loop. No native adapter, source library or default model
was changed in this lane.

## Retained checks

`data/research/thoracic_mechanism/probe_v2/probe.json` records the full generalized
mass matrix, active eigenvalues, exact source receipts, tested pose/velocity,
material orientation checks, cavity Jacobians and power residuals. Six tests
pass in roughly seven seconds under 1 GiB. No compilation or native job runs:

```
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.verify_thoracic_mechanism
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.probe_thoracic_mechanism --output <fresh workspace directory>
```
