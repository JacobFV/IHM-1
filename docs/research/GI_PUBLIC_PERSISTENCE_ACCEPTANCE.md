# GI public persistence acceptance preparation

Status: source-prepared helper for the composition lead's public whole-engine harness. No new native build or patient advance was performed for this preparation. The isolated three-translation-unit repair already passed its native codec fixtures; that is not whole-engine continuation evidence.

## Harness contract

Include `scripts/native_gi_persistence_assertions.h`. On an actual engine loaded from the declared state:

1. `request_tiny_oral_dose(engine)` submits a native gastrointestinal Fentanyl action with 0.001 microgram dose. The fixture rejects preexisting oral transit records. Fentanyl is selected because the held native CAT implementation has actual physicochemical data and a drug-specific metabolism path.
2. Advance one actual 0.02 second tick. Native Drugs PreProcess consumes the action and initializes its transit owner before Gastrointestinal Process runs CAT. No `SetUp` override, manual active-substance registration, or reconstructed dosing history is used.
3. `seed_distinct_transit(engine)` explicitly replaces this test inventory with nine solid, nine dissolved and eight enterocyte masses through existing setters. Values are `(1..9)`, `(21..29)`, `(41..48)` times 0.000001 microgram; cumulative metabolism/excretion are separately 0.000123/0.000456 microgram. This is a synthetic persistence challenge, not a physically reconstructed dose history or conservation test of that seeding operation.
4. Capture, public `SaveStateToFile`, fresh engine public `LoadState`, capture and compare exactly. The observation contains all 28 transit values, record count and SmallIntestineVasculature Fentanyl mass. Missing, extra and unequal entries are reported; no numerical tolerance is hidden. Short decimal fixture selection does not guarantee lossless serialization: any decimal roundtrip difference must remain visible.
5. Advance both uninterrupted control and loaded engine one 0.02 second tick; compare all observations exactly, and call `require_evolved` to reject a nonexecuting CAT path and decreasing cumulative counters. Keep the composition harness's broader state comparison: agreement of these GI observations alone cannot establish whole-engine equivalence.

The lead owns execution, bounds, source/library identities and all raw comparison receipts. Compare pre-save versus immediate reload separately from next-tick comparison so codec loss is distinguished from reconstructed runtime pointer/cache/ordering failures. The local SI vascular mass is an observation at the CAT portal source, not a sum of all redistributed drug or proof of whole-body drug conservation.

## Exact composition input

Use `data/derived/audits/gi-codec-full-tu-3juo1cj5/objects_manifest.json` and `data/research/gi_serialization/source_hashes.json` with `repair.patch`. Ownership, nested and outer patched source SHA256 are respectively:

- `f31d38899a58a515b4539af34d0fe74922e23be326d7d481e7ddc3bf59305ee0`
- `d383eef348ababd828974da7d973d0eb5db9a31d2227c1cf858f0dfbe22156c3`
- `fd581437bd71843138910d68b6243e2832cf08e13108b820de37586c22fe3d73`

The outer `src/io/biogears/BioGearsPhysiology.cpp` also contains sleep and Tissue serialization. Compose source patches into this translation unit and compile it once. Replacing the entire object independently would discard another repair. Original failures remain in `native-gi-codec-_e7npnsp`; corrected isolated codec outcomes remain in `native-gi-codec-o32kvmex`. Unknown historically omitted drug mass is not recoverable and must not be synthesized.

## Mapped lumen timing correction and remaining integration gap

This corrects the prior mapped fixture preparation's rationale for delaying mass balance until circuit PostProcess. In held `include/biogears/cdm/compartment/SECompartmentNodes.inl` lines 63–82, mapped `GetQuantity` returns **NextQuantity**, including the single-node case. Thus mapped liquid `GetVolume` already observes NextVolume, not current volume. Baseline-before-mapping remains necessary for the native quantity-node selection.

Actual ordering is GI PreProcess (including AbsorbNutrients), Cardiovascular Process (circuit Process then substance Transport), GI Process (CAT), and subsequently Cardiovascular PostProcess (next-to-current circuit commit). `SESubstanceTransport.inl` constructs its RHS from extensive mass only when the substance has a valid intensive quantity, then `SELiquidCompartmentGraph::BalanceByIntensive` balances by concentration. A newly populated lumen compartment with an invalid concentration can therefore be omitted from the transport RHS. The prior mapped probe did not call the native solute transporter and cannot exclude that loss pathway. This is a source-derived integration hazard, not an observed patient failure.

Next acceptance must use real mapped NextVolume, full substance quantities, graph transport and PostProcess in actual order. Before graph transport, any manual paired mass/NextVolume update needs valid intensive state against that prospective mapped volume. Explicitly test dry/zero-volume residues, isolated graph vertices, existing SI membership, conservation before/after transport, and repeated epoch rejection. Do not simply balance a dry residue by concentration and assume it survives: the transport matrix substitutes a tiny volume for zero-volume nodes and then performs a concentration-based balance.

A concrete ownership option is to keep proposed manually transported colon/rectum/feces outside the cardiovascular **solute** graph while mapping their fluid nodes into the circuit for volume ownership. That option still requires proving SI boundary mass is balanced before its existing graph transport and that circuit Process preserves manually assigned next volumes. The existing custom mapped-state serialization is not whole-engine schema integration; complete owner state and consumed epoch need a composed native persistence field before activation. No lumen correction is promoted by this header.
