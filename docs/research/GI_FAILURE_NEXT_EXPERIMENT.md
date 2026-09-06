# GI glucose and acid–base failures: next discriminating experiments

2026-09-05. This is a read-only plan based on retained native receipts and pinned source. No native run, build, engine edit, large trajectory load or new physiological calibration was performed for this document. Resource limits in [the incident record](INCIDENT_2026-09-05_RESOURCE_CONTENTION.md) apply to every proposed experiment; none is authorized to occupy a second heavy execution slot.

The next implementation should be **bounded, per-step reaction and transfer counters**, followed by isolated native branch probes. Another uninstrumented hour would establish persistence without identifying the mechanism. The already corrected stale glucose availability must remain corrected, but it has not been shown to resolve the systemic failure. A chemically specified nutrition/counterion interface is a separate necessary correction; do not infer salt identity from sodium mass.

## Evidence and immutable comparison identities

The five retained exertion_v3 records used the same initial state and RH library. Exercise added 150.12 W of modeled metabolic demand from 1800 to 2400 s; this was not measured mechanical work. Exercise glucose reached 32.896 mg/dL at 2460 s; combined meal/exercise reached 44.127 mg/dL at 2620 s. Exercise pH reached 7.2987; combined pH reached 7.5456. Muscle vascular glucose approached zero, and arterial lactate was still rising at 3600 s, 20 minutes after exercise stopped. These observations establish unacceptable trajectories, not a causal diagnosis. See [the readiness audit](EXERTION_V3_PHYSIOLOGICAL_READINESS.md) and `data/derived/audits/exertion-v3-physiology-audit-j33kwere/audit.json`.

| Identity | SHA-256 |
| --- | --- |
| RH predecessor library | `ab485812968b17e16915b04a6e5144bb1e1ec7f5bf78213936c3334d087c0a5a` |
| Common initial state | `cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3` |
| Substrate availability library | `429c63730b3be08b8cf158b16c34817185424cadf7a240abdcd49b4b98d687f7` |
| Corrected Tissue.cpp | `18e98ca5ca6c966f0d065b57f8a5ffca632b6d682dc760db87d1cba9b6dd4fae` |

The substrate variant manifest is `data/runtime/physiology/variants/whole_body_integrity_substrate_availability/manifest.json`; donor revision is `3f16a5fa1dade9c511b88d923606fa51cc35e95d`. Source line numbers below refer to that donor unless explicitly identified as the isolated Tissue replacement. These files were inspected directly, not inferred from their comments:

| Source under `data/raw/physiology/biogears/projects/biogears/libBiogears/src/engine/Systems/` | SHA-256 |
| --- | --- |
| Diffusion.cpp | `9f90e429a87ae1ba3b7071fc14525905e768d23a27a0f21d15a2581a49b7d48b` |
| Hepatic.cpp | `48f3addfeaf0e84630d6d6bf8a390fd53a71f4b73df671d0d1f47d3ce5364e47` |
| BloodChemistry.cpp | `738ecd4a24fd45390943f3233bfacf3ceee3cd71466d94355627696db6add460` |
| Saturation.cpp | `a4bc82a1d63151b53aeed1b55195af04aa0942885f8de702b12842e3e9e574be` |

Before any future executable comparison, bind these sources to its actual object inventory, build manifest, executable, dependency hashes, loaded library and frozen input state. A present donor file alone does not prove that every historical binary executed it.

## Ranked hypotheses and smallest discriminating controls

### 1. Stale glucose availability contributes during exhaustion, but the relevant interval has not been observed

**Proven local defect:** Tissue.cpp:1296–1344 debits aerobic glucose mass; line 1406 formerly reused pre-balance molarity for anaerobic availability. The native scarce-glucose probe turned 6e-8 mol carbon into 1.2e-7 mol carbon. The isolated mass/MW correction closes this glucose-carbon ledger, with scarce oxygen and other substrates covered by 33 checks. See [the substrate audit](GLUCOSE_SUBSTRATE_INTEGRITY.md) and `substrate-integrity-7g17hof9/report.json` under `data/derived/audits/`.

**Not proven:** that this branch produced the hour-long nadir or lactate accumulation. The paired 60 s continuation from the identical final 3600 s exercise state had exactly equal seven sampled trajectories across predecessor and corrected libraries (`substrate-replay-t8cmviiu/report.json`). That continuation cannot repair an earlier excess production event; it also provides no evidence of a benefit during those 60 s.

**Cheapest next experiment:** extend the existing single-tissue branch fixture with read-only counters for pre-aerobic glucose mol, post-aerobic mass/MW mol, pre-balance cached molarity-derived mol, aerobic and anaerobic reaction extents, glycogen extents, and the actual branch predicate. Compare predecessor and corrected execution for one native step with ample glucose, scarce glucose and half-required O2. This is a diagnostic instrumentation parity gate, not a new fix. Then use an actually retained near-depletion checkpoint for a two-library, 1–5 s matched continuation. Do not reconstruct a checkpoint from CSV concentrations. If no suitable checkpoint exists, explicitly budget one serialized instrumented precursor run from the frozen initial state; stop at the first stale-availability event and retain its preceding full state. Event-free short windows are inconclusive about the hour.

**Discriminator:** a nonzero stale availability excess and excess glucose-derived carbon in the predecessor, absent in the correction, proves contribution at that event. Only a later full matched trajectory can quantify its contribution to the clinical extrema. Do not promise glucose normalization.

### 2. Glucose delivery may be limited or misallocated by facilitated diffusion

**Source-supported concern, no native reproduction here:** Diffusion.cpp:313–374 computes each transfer as `Cmax * dC/(Km+dC) * dt`, with concentrations in g/mL and an untyped Michaelis coefficient. Negative gradients are not sign-symmetric; at `dC=-Km` the denominator is zero, and below it flux reverses sign. The held `share/data/substances/Glucose.xml` specifies Km=0.8 and maximum diffusion flux=1 g/cm²/s. A −0.8 g/mL gradient is far outside the reported blood glucose concentrations: the algebraic singularity is **not evidence that this failure occurred in exertion_v3**. The concentration convention and source parameter interpretation need explicit accounting before any unit correction is proposed.

A more directly actionable availability issue is that vascular→extracellular (VE) and extracellular→intracellular (EI) transfers are capped independently using the same starting extracellular pool. With VE<0 and EI>0, both withdraw from that pool; their combined withdrawal can exceed availability despite each individual cap. Whether native gradients and coefficients reached this case remains unobserved. Diffusion.cpp:280 additionally assigns membrane area as tissue mass times an explicit 1 cm²/g prior, not measured regional transport capacity.

**Cheapest next experiment:** one native diffusion call on a three-compartment fixture with nonzero volumes and a held glucose definition. Use matched positive/negative small gradients, zero gradient, extracellular concentration above both neighbors, and a deliberately scarce shared donor. Record raw and capped VE/EI masses, all three starting/ending masses, Km and its concentration convention, dt and combined coefficient. A separate adversarial near-−Km fixture tests finite rejection; label it outside the retained physiological range. Compare the source call with an independently computed mass-incidence ledger, not an invented replacement transport law.

**Discriminator:** negative extracellular remainder or transfer against the intended gradient establishes a local defect. A realistic delivery hypothesis requires actual native V/E/I gradients and fluxes around muscle glucose collapse. If no defect fires, compare integrated glucose delivery with reaction demand; do not tune permeability until its parameter units and evidence are resolved. A replacement signed flux or shared-donor allocation rule requires a separate design because simultaneous proportional allocation and sequential allocation are different models.

### 3. Nutrient-stimulated sodium uptake raises SID without a represented dietary counterion

**Proven source omission and matched physiological effect:** SENutrition has sodium, calcium, water and macronutrients but no chloride or compound identity. Native secretion/absorption pairs sodium mass transfers without a corresponding chloride transfer. In the retained five-case 600 s experiment, nutrients with no added sodium mobilized the initial gastric sodium: net GI→blood Na was 0.500388 g versus 0.019467 g in rest. Vena-cava SID was 42.457903 versus 40.531028 mmol/L, and pH 7.446243 versus 7.411568. Chyme chloride stayed zero. See [the counterion audit](MEAL_COUNTERION_AUDIT.md), `meal-electrolyte-native-02rzs870`, and `meal-counterion-audit-y1aj1a8q/audit.json`.

**Not proven:** the fraction of the combined hour's alkalemia attributable to this omission versus ventilation, lactate, water redistribution or buffers. Initial rest already contains 500 mL gastric water, 1 g sodium and 500 mg calcium; it is not an empty-gut control. Sodium mass does not identify NaCl, and calcium requires its own declared counterion chemistry.

**Cheapest next experiment:** first add observational counters at native GI secretion and absorption, then replay a single digestion step of the retained nutrient/rest fixtures and close paired sodium, glucose and water transfers. No need to repeat the existing uninstrumented five-case 600 s run. Independently test the unchanged Saturation calculator with recorded body-state inputs: change only measured SID, then only pCO2/CO2 state through the calculator's supported interface, then albumin/phosphate. This local sensitivity decomposition is not a dynamically closed whole-body intervention; incompatible solver inputs must not be forced simultaneously.

**Next implementation gate:** extend nutrition to explicitly declared compounds/ions with molar equivalents and conservative secretion/absorption bookkeeping. The first new matched chemistry experiment should use a declared NaCl input, water-matched to its control, with exactly equal Na and Cl mmol added and all other inputs held fixed. Preserve the legacy sodium-only case as an explicitly incomplete chemistry negative control. Paired counterion transfer is an implementation under test, not something the present native nutrition API can already deliver. Do not patch pH or offset SID to conceal the missing chemistry.

### 4. Persistent lactate reflects production, delivery/clearance or demand partition rather than a missing glucose clamp

**Established source structure:** Tissue.cpp:1140–1165 sends 80% of exercise demand to muscle and 20% to fat, adds a basal mandatory muscle anaerobic fraction, and consumes substrates in a fixed order. The fraction starts at 0.028 (line 942); the lactate multiplier starts at 1 and changes under hemorrhage (lines 945–964). Anaerobic glycogen conversion explicitly omits low-pH inhibition (line 1431). Energy::PreProcess leaves ManageEnergyDeficit disabled, so its zero deficit port cannot establish adequate supply. Hepatic.cpp:347–355 converts all available liver extracellular lactate to glucose each call with factor 1; the accompanying O2-cost TODO is not an implemented energy debit.

**Hypotheses:** sustained post-stop lactate can reflect continuing basal/deficit-driven production, slow delivery to liver, redistribution of a large tissue pool, or incorrect reported versus actual demand. The liver conversion itself being fast does not establish whole-body clearance. Missing gluconeogenesis energy cost is a demonstrated structural omission, but does not directly explain hypoglycemia; making glucose without its cost can bias in the opposite direction.

**Cheapest next experiment:** on the already frozen final exercise state, compare continued rest with an otherwise identical zero-exercise command over 1 s, logging actual exercise action/intensity, total and basal demand, exerciseEnergyRequested_kcal, otherEnergyDemandAboveBasal_kcal, tissueNeededEnergy_kcal before/after each branch, mandatory anaerobic demand, lactateScale and hepatic lactate conversion. This chiefly tests stop semantics and post-stop accounting; if both action states already equal zero, equivalence is the expected negative control. In the same interval close total lactate production minus hepatic consumption, renal loss and compartment transfers. A falling muscle pool with rising arterial concentration is redistribution; increasing whole-body lactate with continuing reaction production is not.

Do not immediately add pH inhibition, recalibrate substrate order, or debit arbitrary ATP/O2 for the Cori cycle. Those require a separately sourced reaction/kinetic model after flux ownership is known.

### 5. Other reaction chemistry is incomplete, but must not be conflated with the observed failure

The retained scarce-AA fixture found a −1.0000000000463443e-9 mol carbon residual when the lumped AA is interpreted literally as alanine; glucose carbon and source energy still closed. Hepatic.cpp:370 also declares one glucose per TAG glycerol backbone, while its comment describes a three-carbon glycerol becoming glucose. This deserves a separate reaction audit including the source TAG and ketone formulas and every carbon product; it is not established here as an executed failing branch. Tissue.cpp's previously identified TAG mol/gram comparison has an unreachable pool-deletion branch under held positive inputs and its clamped hormone multiplier; do not edit it merely because its syntax looks wrong.

**Cheapest next experiment:** one isolated hepatic call with only a finite TAG pool, then scarce O2, recording all TAG, glucose, ketone, O2 and other product changes. Independently count carbon using the actual source formulas (not a presumed generic fat). Separately preserve the already failing AA fixture. Reject a claim of complete elemental conservation while lumped chemistry or omitted products remain unresolved. Prioritize these after the transport and glucose branch probes unless event counters show material flux in the failed interval.

## Required ledgers and native observation contract

All new observations are counters on native-owned transfers/reactions, not additional blood or tissue stores. Sample immediately before/after each relevant operation with native time, step index, phase and compartment owner. A once-per-second final snapshot alone cannot reveal reuse of stale pre-balance scalars.

| Ledger | Exact quantities and closure |
| --- | --- |
| Glucose delivery and use | Native V/E/I glucose masses in g, volumes in mL, molarity and MW; actual `massToMoveVE_ug`, `massToMoveEI_ug`; Tissue aerobic/anaerobic glucose reaction extents in mol; GI glucose addition; hepatic glycogenesis/glycogenolysis and glucose from lactate/AA/TAG; renal glucose excretion. For each owner: ending−starting mass = signed transfers + production−consumption. Do not infer actual capped flux from flow×concentration. |
| Glucose-derived carbon | For isolated branches: 6 mol C per glucose/glucose-equivalent glycogen, 3 per lactate, 1 per CO2. `ΔC + C_out − C_in = 0`. Include glycogen debits once, using the native glucose-MW convention rather than silently treating stored glycogen grams as anhydrous polymer grams. Distinguish source glucose-equivalent bookkeeping from a chemical polymer model. |
| Oxygen and source energy | Capture actual O2 debit and extents per substrate. Existing source branch ledger uses 6 O2/glucose, 1.875 O2/AA and 72.5 O2/TAG; energy is 686 kcal/aerobic-glucose mol, 2×energyPerMolATP for anaerobic glucose, 387.189 kcal/AA mol and 7554 kcal/TAG mol. Include separate glycogen ATP yield, brain ketones and unmet demand when those branches execute. Capture heat versus ATP accounting separately. Source energy closure is not a validated human calorimetry result. |
| Lactate | All non-overlapping liquid leaf owners' lactate masses, Tissue production extents, hepatic removal, renal/external loss; actual lymph/vascular/V–E–I transfers. Include CO2/urea/ketones and all known carbon-bearing pools for broader carbon claims. Omitted water/product chemistry prevents full elemental H/O closure. |
| Ions and acid–base | Na/K/Cl/lactate masses and volumes in gut, vascular, E/I, lymph and urine; dietary compound inventory; paired native transfer counters. Moles = grams/MW; charge equivalents = valence×mol. Record VC Na/K/Cl/lactate mmol/L, albumin g/L, phosphate mmol/L, hematocrit, temperature, pH, bicarbonate, dissolved/bound CO2 and pCO2 at solver input/output. Track respiratory CO2 loss and renal acid/base outputs if claiming a total acid–base budget. |

BloodChemistry.cpp:388–390 executes `SID = VC.Na + VC.K − VC.Cl − VC.lactate − 1.02` in mmol/L, despite a stale −5.3 comment. Saturation.cpp:115 uses `SID − bicarbonate − albumin*(0.123*pH−0.631) − phosphate*(0.309*pH−0.469)` in its acid–base residual; `SetBodyState` at lines 191–197 supplies these fields. Record these exact inputs rather than substitute arterial concentrations for the global VC-derived SID. A residual of this reduced equation is not a complete body charge budget, and missing counterions do not justify claiming literal macroscopic Coulomb charge accumulation.

Do not sum parent vascular compartments with their children. Enumerate native non-overlapping owners and external boundaries first. Missing optional species remain absent, not zero. Record renal and lung exits as exports, not unexplained disappearance; preserve unmodeled/lumped species as explicit incomplete closure. Set numeric tolerances from units, scale and native rounding, with the existing 1e-12 mol local branch checks as precedent, not a blanket tolerance for an entire hour.

## Execution order and acceptance

1. Implement observational counters in a separately named diagnostic build, retaining original source and binaries. Bound logs to selected muscle/liver/gut owners and a handful of steps or first-event capture. Require instrumentation-on/off parity before interpreting physiology.
2. Run serialized one-call branch fixtures: stale glucose, shared-donor diffusion, signed gradients, GI sodium transfer, then hepatic TAG if indicated. Preserve failing inputs and outputs before changing numerical source. No new clinical parameter fit or global clamp.
3. Run the smallest matched saved-state windows that actually activate the identified branches. Freeze action history, all resources, library and execution inputs; distinguish checkpoint continuation from an uninterrupted trajectory. Use the failed hour only to select intervals, not to manufacture states.
4. Only after a causal defect is reproduced and an isolated correction passes branch ledgers, schedule a matched initial-state rest/exercise pair through onset, 2400 s stop and recovery. Add meal/exercise after compound-aware chemistry is separately demonstrated. Keep inherited GI, renal, energy, depletion and substrate acceptance gates; never replace the failed hour's artifacts.

The immediate recommendation is **instrumentation plus the shared-extracellular-donor diffusion fixture**, alongside incidence counters for the already corrected stale-glucose branch. That is cheaper and more discriminating than another hour, and can distinguish transport availability from metabolic reaction bookkeeping. Nutrition counterion support should proceed as an independent, explicitly specified chemistry change. None of these proposals establishes that the existing whole-body glucose or acid–base trajectories are physiologically acceptable.
