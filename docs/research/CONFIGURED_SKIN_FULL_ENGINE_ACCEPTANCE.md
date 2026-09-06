# Named configured Skin: isolated full-engine acceptance

This lane uses the configured three-region header and the accepted graph-v2 native
library, with a matched unsplit control. It does not combine newer sleep/GI/CV
experiments, change shared adapters, promote defaults, or claim exact save/reload
restart equivalence. The held stabilized state is loaded explicitly; fresh patient
stabilization is not part of this test.

`native_configured_skin_engine_probe.cpp` derives the accepted full-engine probe and
uses the **exact held** `BioGearsEngine::AdvanceModelTime(bool)` method with only two
regional hooks, after native PreProcess and PostProcess. The native engine's normal
Process, including cardiovascular transport and solute/diffusion systems, runs once
per tick. No extra Albumin/transport call is inserted. Each mode advances twelve
0.02 s ticks, with local 133.322387415 Pa pressure at ticks 6–7 and release at tick 8.

The immutable configuration is `ad93aa3d...` and header `7401a91d...`; source preparation
retains complete hashes. The library is graph-v2 `bc91cbab...`, manifest `52de403e...`.
Frozen body/tissue ports, sweat helper and the native initial state come from the
previously accepted `native-regional-engine-cu9dn9or` run. Native initial-state hash
is `cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3`.

Every output records the configuration digest and semantic region names. The probe
checks the Skin parent's non-owning/read-only role, all three owning-leaf identities,
manager leaf-cache membership and absence of the parent from transport graphs.
The verifier expects the owning-leaf inventory count to increase by exactly two.
It observes all 28 Skin species in every region, shared aggregate quantities, native
sweat donor/waste residuals and per-region/aggregate fluid incidence residuals.
Every native finite liquid leaf contributes **every HasMass species**, with no size
threshold; nonfinite/negative internal mass is rejected. Unknown mass entries remain
counted as unknown. Aggregate parents never contribute again to the leaf sum.

Every tick also exports original-law-cache and configured-branch resistance,
compliance and flow-source values. The verifier checks `R_i*f_i = R`, summed
compliance/flow sources and exact configured fractions throughout native processing,
not only in the installation fixture. This checks configuration propagation into
actual native caches. Only the 16 original active Skin circuit ports known to be
removed by the accepted split may become unavailable; other lost observations fail.

## Initial observer failure and boundary classification

`configured-skin-engine-qsffqz8c` compiled successfully (2.20 s, peak 536,808 KiB),
loaded the held native state, then stopped before tick 0 because the new broad liquid
inventory guard rejected a nonfinite volume. No native ticks or configured-mode run
occurred. Logs and `boundary_diagnosis.json` preserve this first result.

Source/state inspection identified one declared infinite liquid/aerosol boundary:
`Ambient`, mapped to native node `Ambient`. `BioGears.cpp` lines 1164–1171 explicitly
sets infinite node volume and maps both gas and liquid environment compartments to
it. The held XML has that same mapping and no serialized Ambient species quantities;
it contains no negative current volumes. This diagnosis was initially a source/state
inference because the first guard did not include the failing leaf name.

Prepared_v3 accepts **only** `Ambient` with that exact single-node mapping and a
positive infinite volume as an external boundary. It records a positive-infinity
volume flag, every finite/positive-infinite/unknown boundary species observation,
and active aerosol carrier-link flows. These carrier flows are not liquid-water
fluxes; no solute flux is reconstructed from concentrations. The external boundary
is not added to finite internal-store totals, and its classification/observations
must match the baseline. Any other invalid liquid volume now reports its exact
owner name and value and still fails. This is not a blanket nonfinite exclusion.
Finite aerosol-carrier compartments remain part of the explicitly named native
liquid inventory diagnostic, so that sum is not claimed as physiological body water.

## Verification and pending gate

`.venv/bin/python scripts/verify_configured_skin_engine_checks.py` passes three
lightweight validator groups using explicitly synthetic test records: correct
configuration, ownership, law scaling and observed mass requirements; rejection of
duplicate owners, missing species/stale laws/unexpected missing ports; and rejection
of aggregate-mass/configuration drift. Those tests are not native acceptance.

`verify_configured_skin_engine.py` defaults to frozen source preparation verification.
`--run` requires the root's shared native slot: one 1 GiB/60 s CPU/75 s wall compiler,
then three serial native runs, each 2 GiB address space and ≤10 s wall, nice 10 and
one thread. Preparation versions and all failed receipts are preserved. Full-engine
acceptance remains pending until the revised observer's bounded native retry passes.

No whole-body species closure, independent regional chemical laws, clinical lymph
watersheds, pressure-transmission validation or production readiness is claimed.
