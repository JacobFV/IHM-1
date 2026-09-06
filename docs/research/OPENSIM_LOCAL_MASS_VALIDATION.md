# Actual OpenSim local mass validation fixture

The opt-in full-body fixture uses the retained 92-muscle source assembly and
exact continuing native stream implementation copied at build time. It adds two
fixture-only commands in the copied translation unit: a full M/G/q/u/z and
effective-instance-mass observation, and a bounded torso material impulse.
Shared stream files, latest pointers, model properties and installed libraries
are unchanged. The same executable runs serially against original and corrected
Simbody libraries; the corrected setter is resolved only for nonzero transfer.

The baseline comparison covers initialized, warm and continued states: full
mass/constraint matrices, kinematics, continuous muscle states and all existing
native muscle outputs/energy fields. An explicit zero transfer returns before
any setter, scalar rewrite or cache invalidation. Both runs must agree within
a scale-aware 2e-10 tolerance; the maximum observed absolute differences are
retained by field group.

A synthetic 0.5 kg payload is placed at the retained BP3D stomach entity
`body-bp3d-FJ2564`. Its canonical centroid is mapped through the existing common
rigid registration into the native torso frame. The fixture requires that the
existing canonical ownership prior assigns that entity to torso. The retained
site receipt identifies that mapping and its uncalibrated anatomical limits.
This is a local material experiment, not a nutrition action or a claim that a
physiological GI inventory has already been connected. No meal is spread across
segments and no native GI mass is duplicated.

At the impulse time, q, muscle z state and time stay fixed. Incoming spatial
momentum is mapped through the native Jacobian transpose, then new mass-matrix
velocities are solved. Native constraint impulses solve
`G M_new^-1 G^T lambda = -G (u_unconstrained-u_old)` and are applied through
`M_new^-1 G^T`. The two retained CoordinateCouplerConstraints are internal;
no artificial ground constraint or contact impulse is inserted. The fixture
checks generalized impulse incidence, linear/angular momentum, velocity
constraints, gravitational energy transport and constraint-impulse work.
Capture dissipation is recorded separately from active-muscle energy; neither
chemical energy nor heat is fabricated to balance a failed transfer.

A nonzero-energy warm checkpoint is restored before both matched transfer
trials. The instantaneous transfer must preserve all muscle continuous states,
frozen M0/W0/H0 reference, existing chemical/signed-work/heat/positive-work and
external-work ledgers. It then advances actual native mechanics by 0.002 s;
restoring and replaying must recover identical effective mass, full state,
muscle outputs and all ledgers. An Instance-cache invalidation on fixture-only
restore preserves the continuing State without reinitializing the model.

Run `scripts/verify_opensim_instance_mass_native.py --run-native` only in the
coordinated native compile/run slot. Sources, model inputs, exact payload site,
loader-selected libraries, baseline/variant frames and result are retained in
a fresh owned output directory. This fixture does not validate a settled
supine posture, anatomical registration precision or finite-duration fluid
transport. Native physiology inventory/exit receipts and the runtime material
transaction remain a separate integration task.

## Actual retained acceptance, 2026-09-05

`data/derived/opensim-instance-mass-d6hv6h4w/report.json` passed. Across original
and corrected Simbody libraries, initialized/warm/continued full M/G/q/u/z and
native outputs agreed within the stated scale-aware tolerance; maximum absolute
field difference was 1.5535e-10. The fixture exercised all 92 catalog muscles.

The localized +0.5 kg transfer required two nonzero native constraint impulses
(−4.4232e-7 and −4.1332e-7 in their native multiplier coordinates). Generalized
impulse residual was 8.591e-16, linear momentum residual 1.102e-15 kg·m/s, angular
momentum residual 3.732e-16 kg·m²/s and velocity-constraint error 1.001e-16.
Kinetic inflow was 0.0006215993 J, gravitational energy transport −0.3416786780 J
and capture dissipation 0.0006123406 J; constraint-impulse work was 6.617e-23 J.
Only torso effective mass changed. Time, q, muscle z, frozen reference and
already accumulated energy ledgers were unchanged by the impulse. Actual
0.002 s continuation and exact full-state/muscle/ledger restore/replay passed,
as did rejected out-of-scope payload rollback.

Compile plus serialized native acceptance took 8.198 s with 1,219,772 KiB maximum
compiler/native child RSS under the 4 GiB native cap. One earlier retained trial
cleanly rejected a fixture cache-stage mismatch: OpenSim energy/constraint
queries required Dynamics rather than Velocity. The fixture now realizes
Dynamics before these queries; it does not integrate or reset state to fix them.
No default native source, library, latest pointer or physiology intake inventory
was modified by this acceptance.
