# Meal transport ownership and native absorption boundary failures

2026-09-05. A new six-case, one-call native GI fixture reproduced strict remainder-gate failures in the actual inherited `AbsorbNutrients` method. This is a transport availability discontinuity with local mass conservation; it is not another demonstrated mass-creation defect and does not establish the cause of exertion_v3 glucose collapse.

## Native evidence

`python3 scripts/verify_native_gi_absorption_boundary.py` passed **19/19 diagnostic checks**, including required reproductions of five failing transfer cases. Receipt: `data/derived/audits/gi-absorption-boundary-vlie6g6g/report.json`. The executable calls the actual shared-library method; no copied replacement algorithm is executed. The selected `whole_body_integrity_shared_donor` library preserves the inherited GI object from `whole_body_integrity_depletion`, verified against the substrate parent's object manifest. Source/object/library/executable/dependency hashes and resource receipts are retained.

The fixture creates two native 10 mL compartments with actual substance definitions and a 20 ms GI timestep. It calls `AbsorbNutrients` directly without secretion or digestion. BioGears is inactive, so the method's active-only fluid-circuit update is skipped; its nutrient/ion branches and native concentration balancing execute. No full patient state is loaded or advanced. This is not evidence of active fluid-circuit closure.

| Case | Available donor | Native vascular credit | Finding |
| --- | --- | --- | --- |
| Glucose tail; 1 g TAG fixes nutrient-rate driver | 1 µg glucose, 1 g sodium | 0 glucose | Finite partial increment is blocked. |
| Glucose equality; same driver | 33.3333333333 µg glucose, 1 g sodium | 0 glucose | Exactly one requested increment is blocked. |
| Sodium-limited glucose | 1 g glucose, 0.1 µg sodium | 0 glucose and sodium | Finite paired partial extent is blocked. |
| Fat tail | 1 µg TAG | 0 TAG | Finite remainder is blocked. |
| Fat equality | 38 µg TAG | 0 TAG | Exactly one requested increment is blocked. |
| Ample glucose control | 1 g glucose, 1 g sodium | 33.3333333333 µg glucose and 16.6666666667 µg sodium | Full transfer proceeds. |

All six species' combined chyme+vascular masses close within 1e−12 g in every case, and ending masses are finite/nonnegative. The source gates require donor mass to be **strictly greater** than the full requested amount; there is no bounded partial extent. The same strict structure exists in the amino-acid branch, which needs added native cases before its correction is accepted.

A separate source concern is that independent sodium absorption is gated by the **proposed** cotransporter amount rather than an actual transferred amount. The native fat-tail case transfers some independent sodium while the fat-equality case transfers none because the proposed nutrient-driven amount crosses that gate's threshold. Do not silently bundle a new sodium-routing policy into the remainder correction. Missing counterion chemistry also remains unresolved.

The first fixture compile failed on the native template header's undeclared `assert`, retained at `gi-absorption-boundary-j7cbuszf`. Adding `<cassert>` before the native headers fixed the fixture compilation. The successful compile took 7.36 s and 516,848 KiB peak RSS; individual calls took 0.06–0.35 s and at most 40,444 KiB. All child processes had 1 GiB address-space, 120 CPU-second and 150 wall-second limits. The heavy slot was released after verification.

## Concrete next correction

Preserve native coefficient and ordering semantics for this numerical availability correction. For glucose, define the sodium-equivalent extent as `min(requestedNa, currentChymeNa, currentChymeGlucose/2)`; debit that sodium and twice that glucose and credit exactly those same masses to vascular owners. For amino acids, use `min(requestedNa, currentChymeNa, currentChymeAA)` **after glucose has updated the live sodium donor**. For TAG, use `min(requestedTAG, currentChymeTAG)`. Transfer positive extents including exact equality; retain empty/nonpositive-input behavior explicitly. This is a sequential allocation model because that is the held source's order, unlike the simultaneously computed diffusion transfers.

The coefficient `2 g glucose/g sodium` is a source mass ratio, not a declaration of molecular transporter stoichiometry. This task does not alter it or infer any counterion identity. A subsequent chemistry correction requires its own species/molar model and validation.

A corrected native variant must preserve the signed-muscle/shared-donor/thermal/substrate lineage, replace only GI, and retain the original failure receipt. Require original/corrected ample and zero parity, scarce/equality cases for glucose/AA/TAG, multiple sequential withdrawals from sodium, and paired species conservation before linked-library acceptance. No glucose normalization claim follows from that gate.

## Retained meal incidence and portal ownership

The script also reads the retained two 62-row, 600 s meal/rest CSVs, without rerunning them. Under the held source's gram-for-gram stomach-carbohydrate-to-chyme-glucose convention, net gut glucose export is `initial stomach carbohydrate + initial chyme glucose + added carbohydrate − current stomach carbohydrate − current chyme glucose`. This is a source mass-equivalent net lumen ledger, not elemental carbon closure or a directly measured portal flux. Sodium net export uses the corresponding stomach+chyme mass ledger and includes secretion minus reabsorption.

For the nutrients-only case, the first 20 ms has glucose net export **6.856121429296052e−9 g** and sodium net export **−1.217691507524689e−8 g**: sodium initially enters the lumen by secretion. At 600 s, glucose export is **1.295810297240572 g**, with **55.94902318115451 g** carbohydrate still in the stomach and **2.755166521604919 g** glucose in chyme. Sodium net export is **0.5003875364304876 g**, versus **0.0194670143637159 g** at rest. Chyme chloride stays zero. These values confirm slow source meal progression and the known unpaired sodium pathway, but do not isolate transporter throughput from later hepatic/tissue consumption.

The first sampled nutrient-minus-rest sodium net-export excess occurs at 140 s, SID excess at 150 s, and pH excess at 120 s (threshold 1e−9 in the respective reported units). Sampling is every 10 s after the first 20 ms; these are sampled crossings, not exact event times. Different ordering of pH and SID reinforces that the pH response cannot be assigned entirely to sodium without CO2/buffer inputs.

The actual source ownership chain is:

1. `Gastrointestinal::PreProcess`: consume action increments stomach contents, `DigestNutrient` credits `SmallIntestineChyme`, then `ChymeSecretion`, then `AbsorbNutrients` credits `SmallIntestineVasculature`.
2. `Cardiovascular::Process` transports the native circulatory graph. `SmallIntestineToLiver` is the substance graph link mapped to circuit path `SmallIntestineToPortalVein` (`BioGears.cpp:2137–2138`).
3. PortalVein1 is a circuit node mapped into **LiverVasculature** (`BioGears.cpp:1910`), not a separately owned portal substance compartment. Do not invent or double-count a separate portal glucose pool.
4. `LiverToVenaCava` maps the liver return path (`BioGears.cpp:2077–2078`); native circulation carries subsequent systemic distribution. Native hepatic processing and Tissue reaction/diffusion remain the owners of chemistry and tissue exchange.

## Missing observation identifiers

Existing ports expose stomach nutrient masses; `SmallIntestineChyme`, `SmallIntestineVasculature`, `LiverVasculature`, Aorta and MuscleVasculature masses/concentrations; and liver/muscle glycogen. They do not expose actual capped transporter extents or per-link substance fluxes. Add bounded source-owned counters with step/time/phase for:

- `GI.digest.glucose.credit_g`, `GI.secretion.sodium.actual_g`, `GI.absorb.glucose.actual_g`, `GI.absorb.glucose.sodium_actual_g`, `GI.absorb.amino_acid.actual_g`, `GI.absorb.amino_acid.sodium_actual_g`, `GI.absorb.TAG.actual_g`, plus requested extent and each live donor bound.
- Actual circulatory transported glucose/sodium/chloride mass on `SmallIntestineToLiver` and `LiverToVenaCava`, rather than uncapped flow×concentration inference.
- Hepatic reaction glucose sources/sinks, muscle V/E/I delivery and use, and renal glucose export, using non-overlapping native owners.
- GI sodium/chloride boundary incidence and exact VC-derived SID, pCO2/CO2 and buffer solver inputs for acid–base attribution.

The next chemistry correction remains explicitly declared dietary/secretion/absorption ions or compounds with paired molar accounting. Do not reinterpret elemental sodium as NaCl or patch SID/pH to compensate. The shared-donor correction and this absorption-boundary defect each require event evidence in the failing interval before either is assigned responsibility for systemic collapse.

## Implemented native correction and linked-library regression

The isolated correction is now built as `data/runtime/physiology/variants/whole_body_integrity_gi_absorption`. It derives from exact **`whole_body_integrity_signed_muscle_v2`**, library SHA-256 `5a02e2942e6ac5ff3cdce352246e5073de7871b500a572b4ab2d8f50eb87b695`, and replaces only the inherited GI object. All **359 other objects** remain hash-identical, preserving signed muscle demand, shared-donor transport, corrected substrate availability, renal, GI digestion/depletion and thermal lineage. The signed-port ABI `header_sha256` is inherited as `28e3fb59b2680d406bcb0f394b40ec4885bab7d43dbda7e44e1368132f96f545` for adapter build compatibility.

| Identity | SHA-256 |
| --- | --- |
| Corrected library | `9792d857c47a5907f571a03495fe9f4f1144f114afd72c0451869e1a7049588b` |
| Corrected variant manifest | `de1aa254b868b4b51e3fbd4f90323e09370a409957d83ac7552a09bb0f3e7125` |
| Corrected Gastrointestinal.cpp | `804cba46cea17e9b05be3c45c3d20b2e4683dfdb12a5ea8ea7016de415cb6235` |
| Original inherited Gastrointestinal.cpp | `13af1a9350d9da3c4344b6cae187c7ecde44d11cb57358dc39ba7fe5cc1d333e` |

Builder: `scripts/build_biogears_gi_absorption_variant.py`. It verifies the full signed-v2 object inventory before and after build, freezes source/parent identities, refuses to overwrite an existing variant, and retains compile/link/resource receipts. Source outside the no-argument `AbsorbNutrients` method is byte-identical. The independent sodium/water fallback policy, original rate coefficients and glucose-before-AA order are unchanged.

`python3 scripts/verify_native_gi_absorption_correction.py --variant original` first ran the expanded acceptance gate against the old actual linked library. The retained 14-case run at `data/derived/audits/gi-absorption-correction-aph3_yk3/report.json` exits **1**, with **nine expected bounded-transfer failures** out of 42 checks. It adds AA tail/equality, AA sodium scarcity, sequential glucose-plus-AA sodium use and ample controls. The earlier 11-case red run is retained separately at `gi-absorption-correction-z5723dn_`.

The exact same compiled thin fixture was then reused against the corrected library:

```sh
python3 scripts/verify_native_gi_absorption_correction.py --variant corrected \
  --probe data/derived/audits/gi-absorption-correction-aph3_yk3/probe \
  --original-report data/derived/audits/gi-absorption-correction-aph3_yk3/report.json
```

Receipt `data/derived/audits/gi-absorption-correction-4rrrl_ru/report.json` reports **47/47 checks passed**. All cases require finite/nonnegative native ending masses, paired species conservation within 1e−12 g and independently specified bounded transfers. Five controls—ample glucose, ample AA, ample fat, ample mixed glucose+AA and all-empty—have exact original/corrected native-output parity.

Representative corrected actual credits:

- A 1 µg glucose tail credits 1 µg glucose with 0.5 µg sodium; exact 33.3333333333 µg glucose credits the complete requested pair.
- A 0.1 µg sodium donor with ample glucose credits 0.2 µg glucose and exactly 0.1 µg sodium.
- AA tails/equality and sodium-limited AA credit matching 1:1 source mass amounts.
- With 25 µg sodium and ample glucose+AA, glucose consumes 16.6666666667 µg sodium first; AA then consumes the remaining 8.3333333333 µg. Total sodium credit is exactly 25 µg and the donor ends at zero.
- TAG tails of 1 µg and exact 38 µg both transfer completely. The existing independent sodium behavior in the fat-tail case is preserved and not reinterpreted as corrected counterion chemistry.

The production GI compile took 3.18 s with peak RSS 575,756 KiB; linking took 0.51 s and 116,500 KiB. The expanded thin fixture compiled in 2.15 s with peak 517,164 KiB; original/corrected native calls peaked at 40,648 KiB. Production children used 4 GiB address-space/180 CPU-second/210 wall-second bounds; fixtures retained 1 GiB/120/150 bounds. The slot was released after linked verification. No existing library, upstream donor source or signed-metabolic header was changed.

These tests establish the local native remainder correction, including its composable library lineage. Active fluid-circuit integration, whole-body meal trajectories, counterion chemistry and exertion_v3 cause attribution remain separate acceptance gates.
