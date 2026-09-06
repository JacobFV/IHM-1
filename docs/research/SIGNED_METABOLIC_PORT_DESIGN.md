# Signed incremental muscle metabolic port: native source design

2026-09-05. Design only. No native source, engine header or library is modified here. The immediate integration failure is a real decrease in modeled muscle demand relative to the frozen initial mechanical reference, after the separate floating-point domain issue was diagnosed. Do not erase that decrease by clamping it, moving the reference, or encoding a negative Exercise intensity.

## Contract and ownership

The coupling represents a signed **increment about a frozen reference**, not replacement of all BioGears basal metabolism with the absolute OpenSim muscle total. Keep reference model/state/control hashes and fixed reference muscle metabolic power `M0`, signed active-fiber power `W0` and heat `H0=M0−W0`. For each accepted exchange interval of length dt, provide:

```text
Mbar = interval native muscle metabolic energy / dt
Wbar = interval signed active-fiber mechanical work / dt
Hbar = Mbar − Wbar
DeltaM = Mbar − M0
DeltaW = Wbar − W0
DeltaH = Hbar − H0 = DeltaM − DeltaW
```

All quantities are signed increments in W. A negative DeltaM means reduced chemical demand; it does not reverse reactions or export newly synthesized fuel. A negative Wbar is legitimate eccentric absorption and can increase heat relative to chemical demand. Preserve raw roundoff audit metadata independently of these materially signed increments.

Use the same interval quadrature for Mbar and Wbar. The held stream currently integrates metabolic energy and **positive-only** active work, but reports signed work only instantaneously (`scripts/native_mechanical_stream.cpp:90–96,146–148`). Add a cumulative **signed** work/heat ledger with the same checkpoint/restore semantics before claiming interval energy equality. Positive-only work is not a substitute. Native tendon/passive elastic storage and external boundary work remain mechanical owners; signed active-fiber work must not be charged a second time as external work.

BioGears continues to own basal/thermal regulation, tissue reaction choice, oxygen/substrate availability, blood transport, renal/lung exits and fatigue. The new port requests a muscle-only demand increment; it does not mutate glucose, oxygen, ATP or lactate stores. Its heat boundary modifies the native thermal source once. Native chemical supply may fail to meet the request; expose unmet energy rather than silently asserting that requested OpenSim work was supplied.

## Exact held ordering and why the existing hook is insufficient

`BioGears.cpp:971–987` preprocesses Environment, Cardiovascular, Inhaler, Respiratory, Anesthesia, GI, Hepatic, Renal, **Energy**, Endocrine, Drug, **Tissue**, BloodChemistry, Nervous and ECG. `Tissue::PreProcess` at Tissue.cpp:506 handles fluid fluxes, oncotic pressure and burns; it does not execute metabolism.

`Energy::PreProcess` at Energy.cpp:200 calls `CalculateMetabolicHeatGeneration`, sweat, heat resistance and then `Exercise`; `ManageEnergyDeficit` is disabled. `CalculateMetabolicHeatGeneration` at 541–591 updates native TMR for temperature regulation and writes `m_temperatureGroundToCorePath->GetNextHeatSource()` from that rate. **This thermal write precedes Exercise's TMR increment.**

`BioGears::Process` at 990–1005 runs Energy thermal Process before Tissue Process. Tissue.cpp:518–526 performs storage/release, **CalculateMetabolicConsumptionAndProduction**, pulmonary exchange, diffusion and saturation. Therefore an after-all-PreProcess hook could change the pending thermal source before the solve, but cannot replace the source-owned demand calculation or provide correctly ordered observations to earlier controllers. The existing `CoupledBioGearsEngine::AdvanceModelTime` hook is suitable for its respiratory/compression boundary; a single extra scalar setter there is not a complete metabolic port.

## Why existing scalars cannot safely stand in for the port

Generic Exercise at Energy.cpp:243–258 maps an intensity or desired work rate to a bounded intensity; lines 339–347 ramp both ExerciseEnergyDemand and TMR toward BMR plus work. It is not an instantaneous signed chemical-demand interface. At Tissue.cpp:896–897, other demand is `max((TMR−BMR)*dt*(1−biologicalDebt)−exerciseEnergy, 0)`. A lower TMR alone cannot request reduced basal muscle consumption through that floor. Exercise energy is apportioned **80% muscle, 20% fat** (1146–1166), inappropriate for a measured muscle-only delta.

Mutating TMR and restoring a cached base next step is also unsafe. BloodChemistry.cpp:558 increments TMR for radiation stress and line 1217 sets it for transfusion reaction; Energy overrides can set it too. Restoring a saved scalar can discard those native-owned changes, while failing to remove the previous coupling contribution compounds it through the controller. Preserve native TMR as the native-owned base and introduce an explicitly separate effective-demand view.

## Proposed native implementation boundaries

The smallest defensible implementation is a separately named diagnostic native variant with a transactional per-interval port record, plus narrow Energy/Tissue accessors. It is not a Python blood-state correction and not a production shared-header edit without branch tests.

1. **Latch once before native PreProcess.** Record sequence, exact native start/end, dt, fixed reference identity, DeltaM/DeltaH/DeltaW, raw interval integrals and residual bounds. Reject duplicate consumption, missing/overlapping intervals, invalid numbers and conflicting generic Exercise/strength/running/cycling demand for the same muscles. Initial version should reject simultaneous Exercise rather than count both. A zero port must follow untouched native code without scalar round trips.
2. **Native Energy thermal hook.** Immediately after `CalculateMetabolicHeatGeneration` has produced its native source, and before Energy::Process solves it, add DeltaH once to that *next* heat-source rate. Preserve the native source separately as an observation. Never add both DeltaM and DeltaH to heat; `DeltaM = DeltaH + DeltaW`. Sweat/heat-resistance remain native. Placing the narrow hook immediately after heat generation in Energy::PreProcess makes its ownership explicit and precedes later native processing. The effective whole-body source must stay finite and within an explicitly accepted domain; a negative increment is allowed, an unexplained negative absolute source is a failed boundary contract.
3. **Native Tissue demand hook.** Keep original BMR, biological debt, hypoperfusion, brain and other-demand calculations intact. Add `DeltaM * dt` (converted to kcal with native units) **only to the muscle demand allocation**, after its original allocation and mandatory-anaerobic split have been determined and before substrate selection consumes that increment. Do not route this field through ExerciseEnergyDemand or the above-basal `max`. Add the identical amount once to `totalEnergyRequested_kcal` and its per-owner check. Brain/fat/nonmuscle requested budgets must remain unchanged in the isolated branch fixture.
4. **Preflight the decrease against native muscle obligations.** A signed increment may reduce only the eligible muscle budget; it cannot make that budget negative or undo native obligatory chemistry. The source's obligatory AA consumption precedes discretionary substrate branches and can reduce the remaining budget; mandatory anaerobic glycogen demand is accounted separately. Compute/observe the pending native obligation and resulting discretionary budget before committing reactions. Reject an increment below the available decrement with its requested/allowed amounts; do not cap it. Test the threshold at the actual frozen state. A whole-body positive TMR is not sufficient to establish this local domain.
5. **Effective-demand readers, without changing base writers.** Introduce a named effective-demand accessor combining native-owned TMR and the latched DeltaM for consumers that require coupled demand. Audit Cardiovascular.cpp:1655, Endocrine.cpp:292, Nervous.cpp:775/899 and public reporting. Keep native Energy/BloodChemistry regulation writers on their base-owned field. Keep Tissue's original nonmuscle/base partition on the base field; it uses the separate muscle delta once. Do not substitute the effective accessor globally, which would feed the increment back into native other-demand allocation. Preserve existing phase lag: Cardiovascular preprocesses before Energy's current-step update, while Endocrine follows it. Document the latched effective value each observes rather than reordering native systems.
6. **Commit once or restore everything.** Port state, reference identity, consumed sequence, interval energy audit, both simulators and native clocks participate in rollback/checkpointing. Rejecting the domain must not leave a partially consumed muscle pool, thermal step, mechanical interval, or action behind. Native XML alone cannot resume this adapter state unless the port is serialized with a bound sidecar.

The unchanged native heat source currently reflects demanded TMR; Tissue's locally accumulated `heatGenerated_kcal` is only a commented probe at Tissue.cpp:1568 and is not fed to that thermal source. Consequently even correct incremental `DeltaM=DeltaH+DeltaW` accounting does not establish complete absolute native chemical-to-thermal closure. Preserve this limitation. The first coupled acceptance tests must have sufficient substrate/O2 to meet demand; a depleted run requires explicit unmet-demand feedback to mechanics and thermal accounting before claiming energy conservation.

The source's mandatory AA behavior and anaerobic-glycogen branch need particular care for negative remaining demand. At Tissue.cpp:1431–1440, the condition allows entry if mandatory anaerobic demand is positive and then sums it with `tissueNeededEnergy_kcal`. A large decrement must not create a negative combined glycogen extent. A preflight domain check and branch counters are required; simply permitting negative ExerciseEnergyDemand would expose this hazard.

## Frozen-reference consistency gate

The existing arbitrary initial mechanical reference is not automatically a measured resting metabolic baseline. Keep it unchanged for replay, identify that modeling assumption, and compare the requested decrement with the actual native eligible muscle budget. If the observed decrease exceeds that budget, the port should report the inconsistency. It is not permission to borrow brain/fat energy, lower the reference until the request fits, or silently drop the negative part. A later reference change requires a separately specified, retained reference-state experiment and cannot repair the current failed receipt retroactively.

## Smallest matched tests, in serialized execution order

| Fixture | Required observation/acceptance |
| --- | --- |
| Disabled versus explicit zero port, identical frozen state | Exact untouched-source parity where already required; no scalar rewriting, no extra Exercise action, no changed native clock or reference. |
| Muscle-only +epsilon demand with sufficient fuel/O2 | Muscle requested energy changes by `+epsilon*dt`; brain/fat requested budgets unchanged; native source energy, substrate/O2/CO2 and glucose-carbon ledgers close; no artificial 80/20 allocation. |
| Muscle-only −epsilon within measured eligible budget | Muscle request changes by `−epsilon*dt`; native reaction extents remain nonnegative and measured fuel use falls as appropriate to the active branch; no reverse reaction or fuel credit. Fixed reference is unchanged. |
| Demand sign pair with zero signed-work increment | Native thermal-source change is exactly DeltaM because DeltaH=DeltaM; zero return removes the port contribution without drift. |
| Same DeltaM, different signed DeltaW | Thermal-source differences follow DeltaH=DeltaM−DeltaW, including eccentric absorption; chemical demand is unchanged; positive-only work cannot pass this ledger. |
| Lower bound, just-inside and just-outside | Accept admissible decrement; reject excessive decrement before committed effects. Preserve exact failing payload and allowed native budget. |
| Zero/positive/negative/zero sequence | No accumulation of old deltas, no stale Exercise demand, per-interval chemical/heat/work incidence closes. |
| Native stress/thermoregulation with zero port, then small signed port | Native base-owned TMR and stress increments survive; readers observe the documented effective value and phase lag; no cached-base overwrite. |
| Scarce O2/glucose and mandatory AA/glycogen controls | Record unmet demand separately; preserve source substrate correction and reject negative extents. Do not claim mechanical energy supply from requested demand alone. |
| Restore/replay and invalid/duplicate sequence | Exact reference/port ledger rollback; same accepted interval produces the same native result; rejected input consumes neither clock nor energy. |

Start with one-call Tissue budget/reaction and Energy next-source fixtures; do not use full-body clinical extrema to tune this port. Bind source bodies, compiled object inventory, loaded libraries, frozen state and port record. Then repeat the failing integrated 0.02 s exchange and a short matched zero/positive/negative sequence, only after the native local ledgers pass. This design changes a demand boundary; it does not establish physiological validity of the retained whole-body model.
