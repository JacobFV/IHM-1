# Muscle supply coupling: native feasibility and admission contract

2026-09-05. Source inspection and local contract tests only; no native run or
native source change. The held models do not expose an ATP/PCr store or an
energy-limited force constitutive law. A physiological force limiter cannot be
inferred from the existing signed port's unmet demand scalar.

## What the native source actually owns

In `data/runtime/physiology/variants/whole_body_integrity_signed_muscle_v2/Tissue.cpp`,
lines 920–941 use ATP yields and configured energy per mole ATP to convert
substrate reaction extents into energy. They do not integrate ATP concentration.
Aerobic glucose uses 29.85 ATP per glucose with a cellular efficiency whose
numerator contains that same yield; its accounted energy therefore includes
chemical heat, rather than representing an independently stored ATP balance.
Anaerobic glucose and glycogen use yields 2 and 3 with no corresponding aerobic
efficiency division. These branch conventions must be preserved and audited,
not replaced by a single assumed muscle efficiency.

The native owners are intracellular O2, glucose, amino acids, TAG, and the
muscle glycogen scalar, together with CO2/lactate/urea exits. Mandatory AA
consumption precedes discretionary TAG, glucose and glycogen branches. TAG
availability has a hormone-dependent tuning factor (1250–1263); aerobic
branches share finite oxygen; anaerobic branches use residual glucose/glycogen
(1429–1480). Summing caloric stock values would count mutually competing
oxygen/substrates twice and ignore native priority/rate rules. Whole-body blood
concentrations cannot certify intracellular supply over the next interval.

The port's read-only preflight at 1004–1025 measures a lower bound on signed
**demand decrease**, preserving mandatory AA and anaerobic obligations. It is
not a maximum energy supply quote. `unmet_muscle_kcal` at 1483–1484 is an
observation after reactions. Tissue::Process first performs protein/fat storage
and release, then metabolism, then diffusion (522–530). A beginning-of-step
snapshot does not reproduce this phase-local supply state. Native Energy has
already solved the demand-derived heat boundary before Tissue metabolism.

`scripts/native_muscle_metabolism.h` installs a value probe from muscle
activation, fiber states and force/velocity. It does not alter force or own
fuel. The held Millard implementation registers activation and fiber-length
states (`OpenSim/Actuators/Millard2012EquilibriumMuscle.cpp:1284,1287`); its
activation/deactivation time constants and minimum-activation restrictions
(261–263,336–344) make setting excitation to zero neither instantaneous force
removal nor a proof of zero demand. Passive/tendon work remains a separate
mechanical owner. Reducing all recorded work by a supplied/requested ratio
would erase performed work and corrupt eccentric heat accounting.

## Smallest defensible next native port

1. Bind a trial to complete pre-step physiology and mechanics checkpoints,
   source/build identity, native start/end/sequence, frozen M0/W0/H0 reference
   and exact candidate controls. A preview must run the same native phase
   ordering, including prior storage/release and controllers, on a restorable
   transaction. A Python duplicate of substrate kinetics is not that preview.
2. Produce native base/requested/provided/unmet **muscle chemical** energy,
   obligatory allocations, per-branch reaction extents and starting/ending
   substrate/O2 amounts. Preserve brain/fat/nonmuscle budgets. Emit native unit
   conversion identity; convert kcal using the native unit definition, not a
   rounded independent constant.
3. Integrate candidate mechanics tentatively, compute its same-quadrature
   signed M/H/W increments, and native-preview precisely that candidate. This
   is a pre-commit solve; the candidate is not yet accepted physical work.
   Restore both engines and all adapter ledgers after an insufficient trial.
   Never edit the candidate's work or append heat to repair an already accepted
   interval. All preview mutations, thermal effects, counters and reactions
   must be restorable before this loop is valid.
4. Search candidate controls only using actual mechanical reintegration and
   native supply evaluation. No monotonic supply-to-excitation law is assumed:
   activation memory, minimum activation, posture and eccentric motion matter.
   If even admissible minimum controls cannot satisfy obligatory demand, report
   the physical/model domain failure. A force law supporting ATP-limited
   crossbridges requires a separately sourced native model extension.
5. If an ATP/PCr buffer is introduced, give it exactly one native owner, measured
   initial inventory, capacity and kinetic rates, stoichiometric ATP/ADP/Pi/PCr
   transfers, and checkpoint state. Separate chemical release, ATP synthesis,
   ATP hydrolysis/work and heat. Never add the same substrate energy both to a
   buffer and to delivered muscle work. The existing Umberger chemical total
   and Tissue ATP-yield constants alone cannot initialize this extension.

## Implemented local boundary

`ihm/coupling/muscle_supply.py` supplies immutable trial/preview records and a
pure admission check. It validates exact state/control/reference/build/interval
binding, signed M=H+W, nonnegative absolute native budgets, native-base-plus-
increment incidence, and provided-plus-unmet incidence. Missing preview returns
`native_supply_preview_unavailable`; supply shortfall requests reintegration
from the joint checkpoint. It writes no force, work, heat, substrate or ATP.
It is deliberately not wired into embodied runtime: the native preview producer
and complete preview rollback do not yet exist. A supplied chemical verdict
alone is not an absolute thermal closure claim or physiological force model.
The record is an internal producer contract, not cryptographic attestation.

`verify_muscle_supply_contract.py` passes eight local tests: eccentric admission,
missing preview, scarce supply, signed decreases without fuel credit, each stale
binding dimension, malformed ledgers, heat/work mismatch and invalid intervals.
No force limitation or physiological supply validation is claimed by these tests.

## Source identities inspected

| Source | SHA-256 |
| --- | --- |
| signed_muscle_v2/Tissue.cpp | `8ca1eb1c9a09a477cd685adc7f9df2bbb45fdac67bde994490616640788add1d` |
| scripts/native_signed_muscle_port.h | `28e3fb59b2680d406bcb0f394b40ec4885bab7d43dbda7e44e1368132f96f545` |
| scripts/native_muscle_metabolism.h | `5f461b32aee13640a6982c0469745d5d8597dc23287641d22e0058a0afead33d` |
| OpenSim/Actuators/Millard2012EquilibriumMuscle.cpp | `164638c0655d2831fb98a95ea92e53363f4144a9fcd9b412b06bd2686cb7d163` |
