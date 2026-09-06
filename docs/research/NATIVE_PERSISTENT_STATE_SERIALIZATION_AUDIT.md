# Native physiology persistent-state serialization audit

Source-only audit, 2026-09-06. Scope: Gastrointestinal, Renal, Endocrine, Energy and Tissue in the retained BioGears source. No native execution, full build, download or source patch was performed. This complements the independently owned Nervous sleep-history and FluidCircuitPath cardiovascular-region defects; it does not duplicate their fixes.

The new [audit script](../../scripts/audit_native_persistent_state.py) verifies specific reviewed function bodies and emits source hashes, exact line anchors and small branch witnesses. The retained [JSON receipt](../../data/research/native_state_serialization/source-audit.json) records defects present in the audited source, not a passing native replay test. After fixes, its assertions deliberately require review rather than silently blessing changed source.

Paths below share prefix `data/raw/physiology/biogears/projects/biogears/libBiogears/`.

## Confirmed high-risk gaps

### 1. Oral drug transit inventories are read on load but never written

`src/engine/Systems/Gastrointestinal.cpp:789` calls `ProcessDrugCAT()` each Process. The method obtains the drug-transit map at line 819 and evolves nine solid lumen masses, nine dissolved lumen masses, eight enterocyte masses, cumulative metabolism and cumulative excretion; lines 967–972 commit those states and transfer drug to the vascular small-intestine compartment.

`src/io/biogears/BioGearsPhysiology.cpp:475`, `Marshall(const Gastrointestinal&, ...)`, calls only the base writer. `src/io/cdm/Physiology.cpp:410`, `Marshall(const SEGastrointestinalSystem&, ...)`, writes only chyme absorption rate and stomach contents. Neither enumerates `DrugTransitStates`. Conversely, the engine-level reader explicitly iterates `in.DrugTransitStates()` at line 462. The schema supports the optional repeated field (`share/xsd/cdm/Physiology.xsd:413`). This is omitted live inventory, not a recomputable pH/volume parameter table.

Full state load calls `out.BioGears::SetUp()` (`src/io/biogears/BioGears.cpp:76`), recreating the GI object (`src/engine/Controller/BioGears.cpp:121`). A saved state with a partially absorbed oral drug therefore resumes with an empty transit map. The next `ProcessDrugCAT()` returns immediately when that map is empty. Already absorbed vascular drug may remain in serialized compartments, but future delivery from the omitted dose reservoir and accumulated accounting are lost.

**Required native reproduction:** administer an oral drug with CAT parameters, advance until nonzero lumen and enterocyte inventories exist, save, reload into a fresh engine, compare every vector and cumulative counter before advancing, then compare the next tick's portal delivery. Use a settled full-tick checkpoint; a plain nutrition meal is not enough to activate this path. The source-only witness establishes map write omission; it does not claim that this native experiment was run.

### 2. Burn and escharotomy history resets during Tissue load setup

`Tissue::PreProcess()` actively calls `CalculateCompartmentalBurn()` (`src/engine/Systems/Tissue.cpp:506`). `Tissue::SetUp()` resets the following twelve fields at lines 382–396:

- Five accumulated regional resistance increments (`m_*DeltaResistance_mmHg_s_Per_mL`).
- `m_baselineECFluidVolume_mL` and `m_compartmentSyndromeCount`.
- Five `m_*Escharotomy` completion flags.

The burn method increments regional resistance history, captures ECF baseline on its first burn call, computes fluid-creep factors relative to that baseline, counts syndrome events, and uses the completion flags to suppress subsequent regional accumulation. After a successful escharotomy it zeroes the local history, sets the completion flag and removes the one-shot action. These are evolved state, not pointers or geometry-derived caches.

Tissue IO (`src/io/biogears/BioGearsPhysiology.cpp:747–773`) saves metabolic averages and resting masses but none of these twelve fields. Its reader finishes with `BioGearsSystem::LoadState()`, which calls `SetUp()`. Even retained circuit resistance cannot reconstruct which fraction came from this burn accumulator, the pre-resuscitation ECF baseline, or whether a removed escharotomy action previously completed.

**Reproducible branch consequences:** with burn-onset baseline ECF 10 L and current ECF 18 L, the source's fluid-creep factor is 5.4. Load resets the baseline to zero; the next burn call captures 18 L and computes factor 1.0. For an accumulated resistance increment 2.0 and next addition 0.1, the continuous branch contributes 2.1 while the reset branch contributes 0.1 under otherwise matched instantaneous inputs. A completed escharotomy flag returns to false, allowing history to accumulate again. These numbers exercise the source branches; they are not measured patient states or native trajectories.

**Required native reproduction:** capture a settled active-burn state after nonzero regional accumulation and changed ECF, then compare all twelve fields, next resistance contribution and action/event behavior across save/load. Add a separate post-escharotomy checkpoint. Fixing cardiovascular-region structural metadata alone cannot restore this lost history.

## Additional concrete codec/ownership defects

### 3. Nested drug-transit writer copies metabolized mass into excreted mass

`src/io/cdm/Physiology.cpp:459–461` checks that `m_TotalMassExcreted` is valid, but passes `*in.m_TotalMassMetabolized` into `out.MassExcreted()`. With metabolized 3 µg and excreted 7 µg, the emitted excretion is 3 µg. This latent corruption becomes relevant when the missing outer map writer is repaired, or when the nested codec is invoked directly. A direct object with valid excretion but absent metabolism also violates the dereference precondition. A fixture must use deliberately different nonzero counters.

The same nested vector writer pushes default-constructed `unique_ptr<ScalarMassData>` before accessing `.back()`. Whether the generated sequence API accepts this safely was not established here; the repair owner should include nonempty 9/9/8 vectors in a real codec test. This is a follow-up suspicion, not a confirmed additional runtime failure.

### 4. Direct GI invalidation/decode does not clear its owned transit map

`SEGastrointestinalSystem::Invalidate()` (`src/cdm/system/physiology/SEGastrointestinalSystem.cpp:48`) deletes stomach/absorption scalars but never deletes or clears `m_DrugTransitStates`. Its destructor calls that method. `NewDrugTransitState()` at line 155 allocates a new state and overwrites the map entry without deleting a previous object. Direct same-system decode of a state with no transit entries can retain old entries; redecoding the same substance overwrites an owned pointer.

**Scope matters:** full engine load recreates the GI system before decode, so this is not evidence that full engine load retains stale GI physiology. Full engine load instead suffers the inventory loss in finding 1. This omission still affects ownership cleanup and direct codec reuse. Required direct tests: nonempty→empty, same-substance replacement, and repeated decode with distinct objects. Do not conflate a map of pointers comparing equal with value-equivalent serialized state.

## Reviewed apparent omissions that are not demonstrated replay defects

| System/state | Why it was not reported as lost physiological history |
| --- | --- |
| GI transit pH/volume/surface/transit-rate/bile tables | Reassigned constants in `SetUp`; actual evolving drug masses are separate map state above |
| GI stomach contents and digestion phase | Stomach contents pass through base IO; engine reader explicitly sets `m_DecrementNutrients=true` after setup for post-stabilization load |
| GI initial substance masses | Used to restore configured contents at final stabilization, not ongoing active replay history |
| Renal sodium averaging and afferent feedback | Running averages, afferent resistances, sodium-flow setpoints and urination state have explicit two-way engine IO |
| Renal permeability multipliers | Recomputed in `PreProcess → CalculateReabsorptionFeedback → CalculateOsmoreceptorFeedback` before `Process → CalculateActiveTransport` reads them; sodium average is serialized |
| Renal transport scratchpad | Four glucose/lactate temporary masses reset at each `CalculateActiveTransport` before the per-substance loop and gluconeogenesis |
| Renal configuration/pointer members | Bound from serialized configuration/circuit/compartment state in `SetUp`; no additional private history demonstrated |
| Energy pack state and blood averages | `m_packOn`, prior pack weight, blood pH and bicarbonate averages explicitly serialized; `m_Test` is only a debug probe |
| Endocrine pointers/molar masses | Rebound from native substance manager/configuration; hormone quantities live in native compartments; `m_AverageBiologicalDebt` has an assignment but no active read in this source |
| Tissue nutrient stores and metabolic averages | Liver/muscle glycogen, protein and fat have base IO; O2/CO2/RQ/fatigue averages have engine IO |
| Hepatic→Tissue static O2/CO2 scratch | Hepatic Process precedes Tissue Process; Tissue consumes and resets scratch in the same full tick. No gap demonstrated for settled full-tick checkpoints; arbitrary mid-step saving is a different contract |
| Tissue unused/debug accumulators | `m_lastFatigueTime` unused; cumulative fat only appears in a disabled probe; anaerobic-tissue text is diagnostic, not a physiological driver |

These negative controls apply to the reviewed paths, not a completeness proof for every base-class field or every engine subsystem. No new high-risk private-history loss was established in the reviewed renal, endocrine or energy paths.

## Reproduction and repair boundary

```sh
.venv/bin/python scripts/audit_native_persistent_state.py
```

Optional `--output PATH` writes the JSON report. The script neither runs the native engine nor edits sources. It verifies concrete update/read/reset/IO branches, with source hashes and line anchors. Native save/load validation remains required for each isolated repair.

The GI owner is preparing isolated map/vector/counter codec fixtures and cleanup; sleep and Circuit repairs remain with their existing owners. Tissue burn-history repair needs explicit serialization and backward compatibility handling. Old states contain no information from which these histories can be uniquely recovered: fresh initialization is a documented fallback, not faithful continuation. Reconstructible caches should remain setup responsibilities; live inventories and history must cross the state boundary.
