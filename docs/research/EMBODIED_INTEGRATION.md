# Continuing mechanical, neural and physiological integration

This implementation replaces replay-only exchange with a demand-stepped live
runtime. It remains under integration; the presence of adapters and passing
regional tests does not establish a complete digital human.

`EmbodiedRuntime` uses one 20 ms exchange clock. Native OpenSim owns source
segment inertia, articulation, muscle activation, force and contact. The pinned
IBM regional model plus explicit reflex/decoder laws owns neural state. Native
BioGears owns blood, gases, nutrients, interstitial fluid, lymph and heat.
Spatial tissue allocations partition those observations without adding stores.

Neural output from an interval applies during the next mechanical interval,
in addition to the controller's explicit afferent/efferent delays. Native Umberger muscle metabolic energy increments are converted to interval
power and referenced to the fixed initial muscle metabolic power. The signed
difference drives the native generic exercise setpoint; negative differences
are rejected before native commands because that port cannot represent them.
BioGears retains its own basal demand and metabolic ramp, so this does not prove
instantaneous whole-body energy equality. Positive fiber work is reported separately. Initial native patient mass scales
source mechanical inertia; the resulting spatial distribution is a declared
engineering assumption, including how initial GI content is distributed.

The optional `CoupledNativeSession` uses a separately compiled thin adapter.
Original native executables and libraries are retained. A derived engine keeps
the pinned native preprocess/process/postprocess order, event updates and time
bookkeeping, adding external mechanical boundaries after native preprocessing.
Initialization uses the original engine stepping method.

## Executed respiratory feedback

External chest forces are projected through a work-conjugate, three-mode
anatomical load port. Only external load enters the native driver-pressure
boundary; native lung and chest-wall recoil are retained. This avoids adding
the spatial model's passive recoil a second time. It is not yet a fully
resolved pleural, rib-joint or thoracic tissue constitutive replacement.

Actual receipt:
`data/derived/audits/native-mechanical-feedback-rxgksg1c/verification.json`.
Three sequential one-second conditions verify exact original-port parity at
zero load, pressure application residual below1.5e-14 Pa, load release and
changed native lung volume (maximum difference101.2476 mL). The initial failed
attempt exposed a1e-15 unit round-trip difference; zero boundaries now preserve
the original scalar instead of rewriting it. All failed records are retained.

## Executed skin-fluid feedback

An optional dedicated external-reference node and pressure source shifts the
native Skin extracellular compliance boundary. Native oncotic pressure sources
remain unchanged. The original law object remains manager-owned for native
updates; the active replacement copies its laws using the native `Override`
protocol, preserving scalar units and read-only flags.

Actual receipt:
`data/derived/audits/native-skin-compression-s725ktdm/verification.json`.
Three0.4-second conditions verify zero-pressure topology parity (maximum
absolute observed residual1.87e-10), a0.9999987 mmHg skin interstitial pressure
change under a1 mmHg load, and a0.00787095 mL/s lymph-flow change. This is a
uniform whole-skin extracellular boundary; intracellular reference is unchanged.
It does not establish local contact-pressure transfer or long-horizon behavior.
Saving/reloading the added topology is deliberately unavailable until verified.

## Runtime and failure handling

`BodyActor` provides one native-controller thread per body. HTTP callers queue
commands; the body never advances while idle. One active body is permitted
while sharing this machine with IBM-1. Immutable compressed events retain
commands, frames and predecessor hashes. Strict sequence numbers and current
snapshot reads prevent replaying an uncertain force impulse.

Mechanical and neural checkpoints support local rollback. If any native
command becomes uncertain, the body session terminates; it does not pretend
that the native serializer provides exact whole-body rollback. Initial and
failed native outputs remain independently recorded.

Small checks are `verify_respiratory_feedback.py`, `verify_embodied_runtime.py`
and `verify_embodied_actor.py`. Actual native checks are explicit opt-ins:
`verify_coupled_native.py --native` and `verify_native_compression.py`. Run one
local native verification at a time, single-threaded and at nice10. Broad
browser tests and optimization jobs remain off while IBM-1 is training.

Native source equations and the held source revision govern these adapters.
Background methodology is available in the official
[respiratory](https://www.biogearsengine.com/documentation/_respiratory_methodology.html)
and [tissue](https://www.biogearsengine.com/documentation/_tissue_methodology.html)
reports; older tissue documentation does not describe every lymph mechanism
present in the retained executable source.

## Current integration checkpoint

The factory now selects the verified 92-muscle Arm26 v2 artifact and composes
live native lung-volume geometry through EmbodiedRespiration. Its point-load
Jacobian follows current articulation. Segment contact wrenches remain explicitly
unresolved for local thoracic/skin strain rather than guessing their distribution.
Native horizon validation precedes any owner mutation. Partial-startup cleanup
retains a close-capable owner when termination fails. Six small runtime fixtures
pass; the actual integrated native acceptance is still pending at this checkpoint.

Coupled native launches carry a 4 GiB address-space cap and nice 10. This caps
each native child, not aggregate system memory. Initial state and source model
choices remain research assumptions and retain known long-horizon limitations.

The default now verifies the retained initial state and exact library hashes from
`systemic/exertion_v3/exercise/native/manifest.json` (final thermal-corrected
research variant), avoiding the known old thermal state/library mismatch. This
does not accept that variant's failed glucose/acid-base trajectories.

Actual integrated attempt `data/derived/audits/embodied-native-kdfvcq_a` initialized
the paired physiology and 92-muscle plant but failed on the first 20 ms mechanical
advance with `Out-of-domain native muscle metabolism`. No successful integrated
trajectory or viewer promotion is claimed. The earlier source-checker failure
is retained separately at `embodied-native-kzihzxr1`; generated dataclass method
receipts are now handled by a tested source-verification correction.

## Instrumented follow-up

`data/derived/supine-support-tt0k8jes` identified the rejected `ehl_l` power as
-9.4635567759726085e-17 W at positive analysis mass 0.10655017692467443 kg.
The native nonnegative-total calculation can leave a cancellation residual.
The adapter now preserves raw signed power and publishes a per-muscle numerical
tolerance: 64 double epsilons times max(1 W, absolute active fiber power, analysis
mass times the source 1 W/kg minimum-heat scale). It still rejects nonfinite
values, nonpositive mass and materially negative power; it does not clamp energy.

With rebuilt adapter `build-xz0ktmvi`, integrated attempt
`data/derived/audits/embodied-native-bigr07yz` passed the mechanical advance but
rejected a negative metabolic increment relative to the fixed initial reference
before sending physiological boundary commands. A source-aware signed-demand
port is required; redefining the reference or silently clipping would hide this
missing coupling. Production viewer promotion remains blocked.

Independent support run `data/derived/supine-support-5ma720yd` advanced 0.2 s
then hit its predefined kinetic-growth divergence stop (10.002 wall seconds).
The numerical metabolism correction does not resolve physical support stability.

## Signed native variant and intake integration

The runtime now selects `whole_body_integrity_signed_muscle_v2`, verifies the
manifest/library ancestry through the shared-donor and substrate corrections
to the retained thermal base, and binds its mechanical reference to the native
execution receipt. Signed chemical, work and heat interval ledgers replace the
generic Exercise proxy. Local native port acceptance is recorded independently
in SIGNED_METABOLIC_IMPLEMENTATION.md.

Actual full-body attempt `data/derived/audits/embodied-native-4em_78hv` reached
the signed native port, which rejected DeltaM=-36.4893172371 W against an
8.6861257440 W eligible decrement. DeltaH=+0.5663103651 W and DeltaW=
-37.0556276022 W. The fixed mechanical initialization/reference is inconsistent
with that native budget; no reference adjustment or clipping was applied. Native
state may have changed before rejection, so the whole session terminates and
retains its failed log. This is not a successful integrated trajectory.

Intake scheduling now runs on the same actor/controller and sequence as body
steps. Scheduling changes sequence without advancing physical time. Due events
issue native Meal commands once before their next physiological step, with
explicit accepted/uncertain receipts. Nineteen small runtime/actor/HTTP checks
passed for this integration; native whole-body intake acceptance is still gated
on a viable mechanical/physiological initial state. Intake mass changes do not
yet update the mechanical segment inertia distribution.

Whole-engine adapter smoke `data/derived/audits/signed-coupled-k9bxch69` passed
six signed 20 ms steps (0.12 s total) with exact nonboundary-port equality for
disabled versus zero, positive/negative/return boundaries, eccentric heat, and
one acknowledged native water intake consumed by the following step. Wall time
was 1.963 s; caller peak RSS 587276 KiB. These boundary inputs were explicit
small fixture values, not the incompatible mechanical resting reference.


The continuing factory now selects `whole_body_integrity_gi_absorption`, whose
immutable lineage includes signed-muscle-v2 and the shared-donor correction.
The signed adapter was rebuilt at `signed-adapter-9q4tc17x` (4.51 s,
623456 KiB peak compiler RSS). Actual whole-engine acceptance is retained at
`data/derived/audits/signed-coupled-dbzx_jtg`: six signed steps, matched zero,
positive/negative/eccentric heat, return and consumed water intake all pass
(0.12 s simulated, 2.01 s wall). This is a boundary fixture, not a repaired
whole-body mechanical reference or healthy long-run digestion demonstration.
The Python session also retains explicit support for the original signed-v2
variant; executable/library receipts must still match the requested variant.

The articulated and embodied constructors now forward optional `bed_material`
to the native measured-bed law. Bed selection does not supply a supported pose;
its initial mechanical state remains subject to force/moment and joint residual
acceptance. Native receipts retain the chosen curve and separate skin/bed strains.
