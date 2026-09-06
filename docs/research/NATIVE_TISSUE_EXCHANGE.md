# Native-owned tissue exchange and mechanical pressure boundary

`NativeTissueExchange` exposes 13 native tissue regions (Fat, Bone, Brain, Gut, left/right kidney, Liver, left/right lung, Muscle, Myocardium, Skin and Spleen), each with vascular, extracellular and intracellular ownership. Global Lymph and VenaCava are separate native owners. No Python blood, lymph, glucose, albumin or water reservoir is integrated. This replaces the ownership pattern of the older isolated `RegionalSkinTransport` experiment for integrated consumers; that older experiment remains unchanged and must not be added to native mass totals.

## Existing native feedback, preserved

The held BioGears source implements albumin concentration → colloid osmotic pressure → cardiovascular filtration → albumin exchange and lymph return → updated concentration. `Tissue.cpp` 2075–2116 computes the native Landis–Pappenheimer polynomial after a 1.6 albumin-to-protein concentration conversion. Brain is excluded from this update; kidney COP uses renal arterial blood, and Gut has three vascular COP input paths. The helper exposes those source distinctions rather than inventing one Gut pressure source.

`Diffusion.cpp` 300–309 actually calls macromolecule exchange, extracellular-to-lymph albumin transfer, and lymph-to-VenaCava transfer, balancing the compartments afterward. The older `Tissue.cpp` comment near 1996 that lymph is not used is stale relative to this executable code. `CalculatePassiveLymphDiffusion` at 627–660 uses current path flow, source concentration and native dt, caps transfer at available source mass, and ignores negative albumin transfer. These operations remain untouched. Snapshot concentration times flow is **not** reported as an exact integrated solute flux because update order and the availability cap matter.

The native engine also owns osmotic extracellular/intracellular water exchange, nutrition absorption, metabolic reactions, and heat regulation. The previously retained GI depletion, substrate and counterion audits are not resolved by this observer. Free Oxygen/CarbonDioxide scalar masses are not total hemoglobin-bound gas content or oxygen delivery. Temperatures are global native core/skin boundaries, not independently measured organ temperatures.

## Interfaces and ownership

- `scripts/native_tissue_ports.h`: `native_tissue_ports(BioGearsEngine&)` returns native volumes, pressures, albumin/glucose/free gas/electrolyte mass and concentration, Na/K/Cl molarity, gas partial pressure, circuit flow/COP/compliance observations and temperatures. Missing optional scalars remain missing; they are not replaced by zero.
- `NativeTissueExchange.from_workspace(root, reference_snapshot, native_identity)` freezes the source anatomy, reference native allocation, code hashes and caller-supplied native receipt. The caller must supply snapshots from that owned session; plain snapshot values do not independently authenticate their generating process.
- `observe(snapshot)` returns 41 native owners and compact spatial partitions. `project_networks(view)` and `project_anatomical_supports(view)` expand spatial detail only on demand.
- `native_skin_ionic_boundary(view)` in `skin_bioelectric.py` exposes observed bulk native Skin ionic concentrations and corresponding Nernst potentials. It does not assign these to keratinocytes, infer membrane voltage, invent pump current, or write ions into native state.

Every partition has `independent_store=False`. Allocated quantities plus an explicit residual reconstruct each native owner. An incidence ledger applies equal and opposite signs for selected vascular→extracellular, extracellular→intracellular, extracellular→lymph and lymph→VenaCava flows. Its internal sum cancels; it is **not** whole-body dV/dt or total mass closure because vascular perfusion, sweat, renal loss, GI input, metabolic reactions and unselected boundaries also exist. `whole_body_mass_closure_claimed=False` and empty `solute_fluxes` prevent that interpretation.

## Spatial representation

Eight held forearm microvascular graphs retain their actual positions, radii, edges and skin material attachments. Their reference geometric lumen volumes divided by native reference Skin vascular volume provide explicit allocation fractions. The same fractions are transferred to interstitial and intracellular spatial views as an engineering assumption. Edge quantities subdivide their parent allocated volume and mass; the old isolated graph's prescribed pressure/flow state is never imported. Native terminal vessel correspondence and fine pressure gradients remain unresolved.

All 13 native organs receive support records. Existing canonical surfaces provide 257 bone, 425 muscle, 24 selected brain, 60 gut, 10 liver, three myocardial wall, two left lung, three right lung, one per kidney, skin and spleen support entries. Their normalized surface areas are **allocation proxies, not measured perfused volume fractions**. They may be incomplete or overlapping atlas representations and are not certified tissue volumes. In particular the separate muscle surface audit documents self-intersections; these allocations do not repair them or establish mechanical mass ownership. Fat remains a native unlocalized owner because canonical anatomy contains no registered adipose surface. No fat geometry or organ capillary network is invented.

Default observations contain 102 partition rows; the 2,367 anatomical surface partition rows are expanded only on request. This avoids large per-tick display materialization during resource contention. Existing primary source geometry and receipts are reused; no bulk imagery or new download was needed.

## Optional compression topology and actual result

`NativeTissueCompression(engine,"Skin")` installs only after state load. It removes `SkinE3ToGround` from the active circuit while retaining its manager-owned object for native law pointers. A replacement compliance connects SkinE3 to a dedicated external reference; a Ground→reference pressure source imposes the requested nonnegative pressure in Pa. Original native baseline/current/next compliance laws are copied after PreProcess. Existing vascular/interstitial COP sources and pressures are untouched. The intracellular reference stays unchanged: this is an extracellular-only boundary assumption over native **whole Skin**, not a regional patch.

The root-owned adapter calls `after_preprocess()` before native cardiovascular Process and `after_postprocess()` afterward. Native `SECircuitCalculator.inl` 815–855 uses per-scalar `Override<Unit>` to commit read-only current values and resets next laws to baseline. The helper now uses that same protocol when synchronizing the detached original law object, preserving its read-only flags and stored units. The earlier failing `native-skin-compression-u0ofb8kp` read-only exception is retained; it was not bypassed by disabling circuit protection.

Passed actual receipt: `data/derived/audits/native-skin-compression-s725ktdm/verification.json`:

- Duration 0.4 s; maximum zero-boundary absolute residual across the tested source ports 1.8633e-10 in their respective native units, within the root verifier's recorded tolerance.
- Positive boundary changed skin interstitial pressure by 0.9999987337 mmHg.
- Skin lymph flow changed by +0.0078709484 mL/s.

This establishes a short native pressure→drainage causal response. It does not establish local force/area projection, long-term edema behavior, total mechanical/thermal energy closure, or a skin constitutive calibration. Save/reload of altered topology is not accepted; the integrated adapter must reject checkpointing while this helper is installed until a reload path is validated. A local cloth/bed contact pressure cannot be applied to all Skin without an explicit area/distribution model; current contact proxies do not supply that validated mapping.

Frozen headers used for the accepted probe:

- observational helper SHA256 `16b67da259c6228814ec110472634afe1382f5654d6eb1b3aa5560c1b322a42b`
- compression helper SHA256 `7152c37ef968d644cebc8d039f6c5ddabb8985def527d8d36b51a5829959ad4d`

## Verification

Six tiny tests in `scripts/verify_body_exchange.py` cover partition conservation, internal incidence cancellation, missing required fluid values, negative available mass rejection, absent species propagation, overlapping allocation rejection, and refusal to claim unobserved solute flux/global closure.

The retained-data mode was run against `native-mechanical-feedback-s_il1bnp/coupled_zero`, without another native job. Report: `data/derived/audits/body-exchange-35x_9nvk/verification.json`. All observed owner volume and species mass partition residuals were exactly zero; internal flow incidence residual was 1.76e-17 mL/s. Native mass-versus-concentration residual was at most 1.42e-14 g. At 1 s, bulk Skin Nernst potentials were Na +59.26, K −85.28 and Cl −81.03 mV; these are not measured membrane voltages.

All work follows `INCIDENT_2026-09-05_RESOURCE_CONTENTION.md`: no additional native jobs, builds, browsers or downloads from this worker after the incident; tiny Python checks use one OpenBLAS/OMP thread and nice 10. The broader unclosed objectives remain in `WHOLE_BODY_REMAINING_CONCERNS.md`.
