# Configured Skin species transport: bounded fixture preparation

The new `native_configured_skin_species_fixture.cpp` exercises actual native
species operators on the accepted configured three-region circuit. It is an
isolated operator/compartment fixture, not full-patient composition or evidence of
independent regional chemical laws. Existing headers, variants and adapters are
unchanged.

The fixture imports the accepted configured circuit setup from its frozen source.
It seeds **all 28 species** exported by the retained native regional installation
snapshot, using each observed aggregate mass divided by aggregate volume. These
concentrations seed the fixture's 100 mL extracellular inventory; no new measured
concentrations are asserted. Originally zero species remain zero. Species names,
source snapshot SHA and generated seed bytes are retained in
`data/research/configured_skin_species/prepared_v1/`.

Two separate operator phases prevent duplicate transfers:

1. The actual `SELiquidTransporter` acts on one selected extracellular-to-Lymph
   link for each owning leaf (three configured versus one unsplit). Links map
   the native `SkinE3ToSkinL1` paths. There is one shared 10 mL Lymph recipient.
   The aggregate Skin parent is explicitly excluded from graph vertices. No serial
   duplicate links, vascular filtration, sweating, external waste or second lymph
   store is added. The fixture locally instantiates the exact pinned native
   `SESubstanceTransport.inl` specialization because the held CDM library does not
   export it. This is the native matrix transport algorithm, not a replacement.
2. A freshly initialized independent fixture pair invokes the actual native
   `DiffusionCalculator::CalculatePassiveLymphDiffusion` for **Albumin only**, as
   the production nonlinear diffusion caller does. Access uses the header's
   existing `BioGearsEngineTest` friend; no private-access macro or copied diffusion
   body is used. The method runs in the pinned graph-v2 library and consumes the
   original aggregate flow cache, retaining shared native chemical-law behavior.
   The other 27 species remain in the complete conservation inventory but are not
   falsely treated as substrates of that Albumin-specific production pathway.

Both phases prepare ten matched zero-pressure circuit steps, followed by local
pressure and release. Every step checks finite, nonnegative owning masses, unchanged
sum of all owning extracellular leaves plus Lymph, read-only aggregate identity and
zero-load aggregate transported-mass parity. The generic phase additionally requires
pressure-dependent **local Albumin mass transfer**, not merely different pressure.
The passive phase then supplies a deliberately huge operator timestep to request
more Albumin than available and check the actual native donor cap and exhausted
repeated call. No patient or circuit is advanced over that stress interval.

This closed species ledger does not imply a complete physiological fluid ledger:
only selected drainage transport is represented; Lymph volume is a finite fixture
recipient volume, not a patient-resolved lymph pump. Generic convection of all 28
quantity types is an operator test, not a claim that native skin lymph transports
all those types in vivo. Aggregate native Albumin diffusion remains a shared law;
this test does not validate regional oncotic/permeability feedback or separate
regional reaction kinetics.

Sources remain pinned to configured header `7401a91d...`, configuration
`ad93aa3d...`, parent graph-v2 manifest `52de403e...`, library `bc91cbab...`, and
actual seed snapshot `cb37d4e4...`. The staged receipt additionally hashes both
native operator sources and all fixture sources. Full ancestry is retained in the
parent configured preparation. The graph-v2 detached-path lookup correction is
necessary for the native passive Albumin operator to find its aggregate law cache.

Commands:

- `.venv/bin/python scripts/verify_configured_skin_species.py` prepares/verifies
  frozen source and seed receipts without compiling or running native code.
- `.venv/bin/python scripts/verify_configured_skin_species_source.py` checks complete
  actual species inventory/byte identities and rejection of a changed frozen header.
- `--compile-and-run` on the first command is the queued native gate: 1 GiB address
  space, 60 s CPU/75 s wall compiler, then ≤10 s native fixture, serial, threads 1,
  nice 10 under the root's explicit shared-slot grant. Every failure and resource
  log is retained in a new audit directory.

At source preparation, native execution is **pending**. Full-patient acceptance and
any production composition with sleep/GI/cardiovascular changes remain separate.
