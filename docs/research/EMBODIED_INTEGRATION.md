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
