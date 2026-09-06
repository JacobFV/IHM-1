# Local intake/excretion mass: native instance prototype

2026-09-05. The current public Simbody API cannot update a mobilized body's
mass in a continuing State: `MobilizedBody::getBodyMassProperties(state)` is
const, while `setDefaultMassProperties` explicitly invalidates topology and
requires realizing it again. Doing the latter inside a meal action would risk
losing muscle state, integration state and opaque checkpoint validity.

The retained source has an internal path:
`SimbodyMatterSubsystemRep::updBodyMassProperties(State&, MobilizedBodyIndex)`
updates `SBInstanceVars.bodyMassProperties` through `updDiscreteVariable`.
Instance realization recomputes an internal total-mass cache from those variables.
The public getter also ignores them: `MobilizedBodyImpl.h:257–260` contains
an explicit TODO and returns topology mass; public system mass sums those getters.
However, deeper source tracing shows dynamic spatial inertia is built from
`RigidBodyNode::massProps_B`, an immutable topology member (`RigidBodyNode.h:1071`).
`RigidBodyNode.cpp:75–83` obtains inertia, COM and mass from this member, not the
updated instance properties. Thus a private instance write can change reported
internal mass fields without changing public mass or the actual equations of motion. It is not a valid mass port.

## Executable falsification fixture

`native_local_mass_instance_fixture.cpp` modifies one independent free body's
instance mass, first moment and inertia by a localized 0.5 kg point payload.
A second free body checks locality. Before and after it measures the native mass
matrix and kinetic energy at identical configuration/velocity. A co-moving
payload must contribute `0.5*dm*|v_at_payload|²` kinetic energy, and linear/angular
momentum `dm*v` / `x cross dm*v`. The fixture checks whether the held instance
write produces those required changes and reports the mismatch. It preserves
time/configuration/speed and checks State-copy restore/replay of the instance
fields. It intentionally does not use an unchanged mass matrix to compute a
pretend mass-transfer velocity correction.

The expected source limitation is `variable_mass_supported=false`: the native
mass matrix and kinetic energy stay unchanged despite the internal instance mass rising.
A test result of `passed=true` means that this limitation was reproduced, not
that variable mass works. There is no runtime mass updater or body integration.
The fixture includes no muscles, constraints, contact or finite-duration flow.
The script retains the translation unit, hashes private headers and libraries,
and uses the coordinated 4 GiB native resource slot.

A supported implementation must route all rigid-body mass/COM/inertia consumers
through the State instance, including mass matrix, velocity-dependent forces,
spatial momentum, kinetic/gravitational energy, COM operators and constraint
impulses. Updating a report field alone is insufficient. Do not mutate const
node members or globally swap topology mass to evade this failure: that breaks
concurrent States and rollback. A separately built native source variant with
matched impulse/energy tests is required before ingestion changes mechanics.

## Native physiology ownership before integration

The held GI source adds `Nutrition.GetWeight` to patient Weight exactly when
ConsumeNutrients increments stomach contents (GI PreProcess 290–292). Nutrition
weight sums carbohydrate, fat, protein, calcium, sodium and water; water uses
the explicit native 1 g/mL convention (`SENutrition.cpp:281–297`). A second
mechanical addition from both the action and patient weight change would count
the meal twice. Initial native stomach contents are also already included in
the native patient weight during initialization.

Digestion transfers stomach pools to chyme; it is an internal redistribution,
not another external mass input. Absorption similarly transfers GI fluid and
solutes into vascular/tissue owners. To affect the mechanics locally, the
coupling needs a retained baseline inventory and nonoverlapping physiological
compartment-to-segment mass/COM/inertia maps. The incremental stomach payload
belongs to an explicitly registered abdominal site, not a uniform scaling of
all bodies. Chyme, absorbed fluid, bladder contents and tissue redistribution
must move that already-owned mass between registered sites without changing
system mass a second time. Assuming all liquid volume is water mass and then
adding separately counted dissolved solids needs a declared density convention.

Patient Weight is not a complete boundary ledger. Energy subtracts sweat mass
(`Energy.cpp:655–657`) and also changes Weight for a carried exercise pack
(325). The inspected Renal urination path resets bladder next-volume to 1 mL
and rebalances solutes (1563–1573); it has no matching Weight decrement in this
function. Its reported urine flow alone may not reproduce that instantaneous
volume reset. Therefore neither a Weight-only delta nor a flow-only urine
integral can certify complete excreted mass from the held source. A native exit
receipt must expose actual before/after fluid and solute mass, source/sink
identity, location and boundary velocity, bound to the same accepted sequence.
Respiratory gas, sweat and other exits need analogous ownership.

## Required integration boundary

Add a native per-owner material ledger with baseline inclusion flags, donor and
recipient compartments, signed external flux, local station/distribution,
transport velocity and energy. Reject duplicate sequence consumption; bind
native material snapshots and mechanical instance properties in joint rollback.
Require total material transfer to equal the sum of mechanical mass changes and
external exits. Validate central inertia and positive mass after each update.
A source-owned supported Simbody instance setter and an OpenSim mass-reporting
policy are preferable to shipping an unversioned private-header dependency.
Do not edit static OpenSim model properties and pretend those changed the
continuing System state, or silently call initSystem to make them take effect.

## Measured falsification

`data/derived/local-mass-instance-x97iyjps/report.json` reproduces the limitation:
internal mass 10→10.5 kg, public body mass 10 kg, public total mass 20 kg, mass
matrix change norm 0, kinetic energy change 0 instead of required +0.2938 J,
and momentum flux error norm 0.5516402813 (mixed angular/linear diagnostic norm,
not a physical single-unit quantity). Payload gravitational transport would be
+1.4709975 J. Maximum compiler/native child RSS was 484,980 KiB. State-copy
restore/replay passed after explicit Instance-cache invalidation. The retained
first attempts document a C++17 compile failure (held headers require C++20),
the public getter's topology behavior, and a cache-prerequisite failure before
that invalidation. No successful variable-mass claim follows from this receipt.

## Pinned native correction plan

Keep the current installation untouched; build a separately identified Simbody
variant from retained revision `59c6e7b89b3bdf266a2f3d54c611599d964205f7` and bind
all OpenSim loaded-library identities. Concrete source edits needed:

1. Add a validated public State-instance mass setter and change
   `MobilizedBodyImpl::getBodyMassProperties` to use `SBInstanceVars` rather
   than the topology node. Initialize from topology exactly once. Validate
   finite positive mass and admissible central inertia; invalidate Instance
   and all dependent caches. A State-copy must own its own mass properties.
2. Pass State/instance mass through
   `RigidBodyNode::calcJointIndependentKinematicsPos`; update spatial inertia,
   world COM and all position caches from it. Update the velocity gyroscopic
   force in `RigidBodyNode.cpp:121` to use the same instance mass. Its current
   signature accepts only tree caches; thread instance properties through the
   call chain instead of caching mutable mass on the shared node.
3. Audit specialized paths: `RigidBodyNode_Weld.cpp:412–420` repeats topology
   spatial inertia; `RigidBodyNode_LoneParticle.cpp:175,187–188,292,327,372`
   directly uses topology mass in inverse mass, momentum and acceleration.
   These must receive State instance properties too, even if the first fixture
   uses only Free bodies. Keep topology-only shape/mobilizer choices separate.
4. Verify gravity consumers (`Force_Gravity.cpp:546`, `Force.cpp:1049,1075`),
   system COM/momentum/kinetic energy, M/M-inverse, inverse dynamics and all
   articulated-body inertia cache prerequisites against the same properties.
   A public getter fix alone cannot repair dynamics; a dynamics-only fix
   leaves gravity/momentum/reporting inconsistent.
5. Upgrade this exact falsification into positive matched tests: co-moving
   payload gives required ΔK/ΔP and local inertia; different-velocity capture
   closes kinetic transport plus dissipation; passive outflow preserves carrier
   velocity; joint/constrained impulses and gravitational energy close. Verify
   two concurrent States with different mass, full checkpoint replay, and an
   actual OpenSim muscle state plus signed-work/metabolic ledger survives.
   Only then add a native material port to the continuing mechanical stream.

This plan is an engine implementation task, not permission to redistribute a
meal uniformly or add compensating forces to hide unchanged inertia.
