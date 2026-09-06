# Explicit endpoint local material port

The production adapter has an opt-in `instance_mass_variant` argument. The current
validated library path is `data/runtime/opensim/variants/instance_mass_v1`.
The baseline Simbody installation is unchanged. The normal adapter is built with
`python scripts/build_native_mechanical_stream.py`; a separately selected adapter
is built with `--instance-mass-variant <path>` and gets its own pointer keyed by
library SHA. Selection requires a complete, unchanged owned manifest/library,
the matching adapter manifest, and an actual loaded-library `/proc` receipt.

Both `NativeMechanicalStream` and `ArticulatedBodyPlant` expose:

```python
point = plant.body_point(body='torso', station_m=station)
frame = plant.transfer_mass(
    sequence=1, owner='consumed-meal:1', body='torso',
    delta_mass_kg=0.25, station_m=station,
    velocity_source_m_s=point['velocity_source_m_s'])
```

`station_m` is torso-local. `point_source_m` and `velocity_source_m_s` from
`body_point` are in the OpenSim native ground frame, not the canonical display
frame. Query output includes time, body and exact local station. Co-moving inflow
is an explicit physical assumption supplied by the caller. No intake, absorption,
GI-fluid or excretion quantity is inferred by this API. Physiology mass receipt
wiring must identify the actual consumed boundary once and declare the mapping.
The isolated whole-model fixture's stomach mapping is a rigid registration prior,
not a deforming stomach model.

Each transaction has finite signed mass magnitude at most 0.5 kg. Total explicit
owned payload across all owners must remain in [0, 0.5] kg. Larger meals are
rejected; silently splitting them to bypass this domain is unsupported. Positive
transfer creates or increases an owner inventory. Its torso/station binding is
immutable. Removal is limited to that owner's added inventory and requires
co-moving station velocity within 1e-10 m/s; jets and original tissue/GI removal
are unsupported. At most 64 owner records persist, including depleted owners.
A zero delta commits a sequence/zero receipt but creates no owner or State mass
mutation. Each sequence must be the next positive uint64. The native reference
is bound to the live stream identity, never accepted from application input.

The native constrained impulse uses effective State mass properties, updates
local mass, first moment and inertia, and solves native constraint impulses.
The native receipt retains source position/velocity, local station, endpoint time,
body mass before/after, kinetic and gravitational-potential inflow, capture
loss, constraint impulse work and momentum/constraint residuals. A copyable
ledger exposes cumulative fluxes separately from the muscle energy/work/heat
ledger. Muscle State, q, time, and muscle ledger/reference are preserved at the
endpoint. Native snapshot `mass_kg`, body mass/COM/inertia are effective State
quantities; each body also exposes `model_baseline_mass_properties` for the
original assembled/scaled model. Canonical entity mass remains its reference
registration; the articulated frame exposes `effective_native_body_mass_kg`.

Opaque native checkpoints retain State instance mass, complete continuous State,
loads/excitation, muscle ledgers, and payload owner/sequence/flux ledger. Restore
invalidates from Instance. A native exception restores all of these before an
explicit rejection reply. Unknown transport outcome, short write, malformed reply,
or mismatched success receipt closes the native process. Such a session cannot
be retried or restored; application-level coupling must treat it as failed.
No remote/process-restart checkpoint format is implied.

`verify_native_local_mass_port.py` runs five light selection/transport tests.
Its coordinated `--run-native` mode compares ordinary/opt-in initialization and
body-point output, exercises cumulative inflow and co-moving outflow, owner/site
and total-domain rejection, nonzero State/sequence replay and muscle continuity,
and compiles an isolated deliberate post-setter throw to verify rollback and
subsequent continuation. The fault macro is compile-time only; there is no
production command to enable it. Prior positive whole92 physics evidence is
retained in `OPENSIM_LOCAL_MASS_VALIDATION.md`; production acceptance gets its
own receipt and does not replace that evidence with a source-only assertion.

## Measured production acceptance (2026-09-05)

Receipt: `data/derived/native-local-mass-port-heewq9qr/report.json` and
`states.json`. Ordinary adapter `build-3m_az4tv` compiled in 6.138 s and opt-in
`build-nov097u4` in 6.527 s (maximum compiler RSS about 1.173 GiB). The serialized
actual92/fault fixture passed in 8.393 s, maximum child RSS 1,227,268 KiB under
the 4 GiB ceiling. Zero-mode observable parity maximum absolute difference was
4.12115e-13, including effective/baseline properties and muscles.

Two +0.25 kg transfers raised total effective mass from 77.6122029 to
78.1122029 kg; torso mass became 28.154676965260336 kg while its baseline stayed
27.654676965260336 kg. Last-addition generalized impulse residual was
4.02e-16, linear/angular momentum residuals 3.34e-16 / 3.22e-16, constraint
velocity error 1.39e-17. Cumulative incoming kinetic energy was 0.00035 J and
capture loss 0.0003509746299801275 J (the latter also draws on preexisting body
kinetic energy when capture slows it). Two co-moving -0.25 kg transfers returned
mass and local inventory to baseline/zero; final outflow capture loss was
-2.78e-16 J roundoff and cumulative potential flux returned to zero. Removed
mass leaves with current material velocity, so original whole-body momentum is
not expected to be restored after loading then removal.

The fixture also passed owner/site/sequence/domain rejections, nonzero muscle
ledger preservation, continuing-State replay, and failure after the setter had
mutated Instance mass followed by successful rollback and continuation replay.
Five local tests separately passed selection and poisoned-transport behavior.
