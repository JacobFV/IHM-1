# Isolated Tissue burn/escharotomy history repair preparation

The source-only audit found twelve evolved Tissue members reset by `SetUp()` and absent from native IO. This preparation preserves that owner state through a versioned `IHMBurnHistory` payload. It changes no held source, default library, physiology coefficients or circuit behavior. No schema regeneration, C++ compilation or native execution has been performed for this repair.

## Exact owner contract

The payload is canonical ASCII:

`IHM_BURN_V1:<two lowercase hex flag digits>:<seven 16-digit IEEE754 binary64 hex words>`

| Encoded state | Native owner members | Units/order |
| --- | --- | --- |
| Five resistance histories | `m_trunkDeltaResistance_mmHg_s_Per_mL`, `m_leftArmDeltaResistance_mmHg_s_Per_mL`, `m_rightArmDeltaResistance_mmHg_s_Per_mL`, `m_leftLegDeltaResistance_mmHg_s_Per_mL`, `m_rightLegDeltaResistance_mmHg_s_Per_mL` | mmHg·s/mL, in that order |
| Syndrome counter | `m_compartmentSyndromeCount` | Nonnegative integer-valued native double; preserves actual counter, not recomputed event cardinality |
| Historical ECF baseline | `m_baselineECFluidVolume_mL` | mL; captured at first burn call; not current ECF |
| Five completion flags | Corresponding `m_*Escharotomy` booleans | Bits 0–4, same region order |

Seven doubles and five booleans preserve all twelve fields. Binary64 bit encoding avoids native XML decimal precision loss. Values must be finite and nonnegative; count must be integral and at most 2^53; flags must fit five bits. Any positive accumulated resistance, nonzero counter or completion flag requires positive historical ECF baseline. An entirely fresh state has all zeros, including baseline. No sodium/chloride, pressure or resistance is invented to reconstruct missing history.

The C++ helper is [native_tissue_burn_state.h](../../scripts/native_tissue_burn_state.h). The Python preparation tool implements the same wire contract for explicit offline initialization and light checks. C++/Python runtime parity is staged for compilation, not claimed verified by Python tests.

## Initialization and load order

The local Tissue patch separates the lifecycle:

1. `Invalidate()` makes the seven numeric history values NaN and clears flags. Such an owner cannot be saved or processed as initialized.
2. Fresh `Initialize()` calls normal setup and then explicitly establishes zero no-burn history.
3. `SetUp()` continues rebuilding native pointers/caches but no longer zeroes burn history.
4. Tissue `UnMarshall()` requires, decodes and validates the payload **before** invalidating the target. It restores base/scalar state, runs `BioGearsSystem::LoadState()`/`SetUp()`, then restores all twelve private histories.
5. Tissue `Marshall()` validates its private owner and emits the exact payload. The burn preprocessing method validates history even when no current burn action exists.

This fixes loss of pre-resuscitation ECF reference, accumulated resistance and completed-escharotomy flags. It deliberately preserves existing burn formulas and the existing counter semantics, including any independent defects in those formulas. It does not derive count from patient event flags or infer prior procedures from absence of an action.

The payload is an optional `xs:string` only for XML/schema parsing compatibility. The runtime reader rejects its absence; optional syntax does not mean zero is a valid restore default. Malformed/version-unknown payload rejection is checked before mutation. Errors later in base IO/native setup are not claimed transactional for the whole engine.

## Legacy states and explicit fresh initialization

Old burn histories are unavailable from old XML. Existing circuit resistance, current ECF, event flags and removed one-shot actions do not uniquely identify those histories. Missing payload therefore fails closed, even if no burn action is currently present.

An independently justified **new no-burn initialization** can be prepared with a declaration JSON:

```json
{
  "schema": "ihm.explicit-fresh-tissue-burn.v1",
  "declaration": "fresh_no_prior_burn_or_escharotomy",
  "provenance": "Describe why this is a fresh no-burn preparation"
}
```

The migration refuses burn-wound, escharotomy or compartment-syndrome markers and refuses to overwrite an existing payload. Its declaration is explicit external evidence; marker absence alone is never treated as recovered history. It inserts only one namespace-correct Tissue child into the original bytes, preserving all other XML bytes and namespace declarations. A separate receipt labels this as new zero-history initialization, not recovery.

```sh
.venv/bin/python scripts/prepare_tissue_burn_state_patch.py \
  --seed-fresh-legacy INPUT.xml --declaration-json DECLARATION.json \
  --output NEW-STATE.xml
```

No migration of a real patient state was performed in this task. For a legacy active/prior burn, actual history must come from an external retained owner record; no generic zero-fill migration is supplied.

## Composable isolated preparation

[prepare_tissue_burn_state_patch.py](../../scripts/prepare_tissue_burn_state_patch.py) exports `patch_tissue`, `patch_io` and `patch_schema`. Transformations are scoped to Tissue functions or the `BioGearsTissueSystemData` complex type. Duplicate application and missing anchors fail. It never replaces another owner's Nervous/renal/GI blocks.

Default inputs are guarded by the retained [source hash manifest](../../data/research/tissue_burn_serialization/source_hashes.json). For a parent already containing sleep/GI fixes, supply an isolated directory with `Tissue.cpp`, `BioGearsPhysiology.cpp`, `BioGearsPhysiology.xsd` and an explicit expected-hashes JSON (`source_sha256` with `tissue`, `io`, `schema` keys). All three input hashes are checked before output creation.

```sh
.venv/bin/python scripts/prepare_tissue_burn_state_patch.py --output NEW-DIRECTORY
# Composed parent, including existing independent repairs:
.venv/bin/python scripts/prepare_tissue_burn_state_patch.py \
  --parent-dir PARENT-DIRECTORY --expected-hashes PARENT-HASHES.json \
  --output NEW-DIRECTORY
```

The new directory contains patched source copies, unified diffs, helper, staged fixture and a receipt. Existing source/variant files remain unchanged. Light tests confirm that sleep and Tissue schema patches commute exactly and that both orders of IO patching retain both payloads and correct Tissue restoration order.

**ABI constraint:** adding this schema field requires matching schema regeneration and dependent CDM/core rebuild. Do not combine objects generated against original, sleep-only and combined schemas. The sleep worker reports three dependent CDM and nine core objects for its schema change; the combined dependency closure and Tissue translation unit still need coordinated build verification. No such build was started here.

## Verification completed and remaining

Run the light checks:

```sh
.venv/bin/python scripts/verify_tissue_burn_state_patch.py
```

Six test groups pass: original zero-baseline/first-burn capture timing, exact Python codec/finite-domain validation, lifecycle source transformation, both sleep compositions, source-hash guarding/missing-history rejection, and byte-local explicit fresh-state insertion with marker/overwrite rejection.

[native_tissue_burn_state_fixture.cpp](../../scripts/native_tissue_burn_state_fixture.cpp) is staged but **uncompiled**. It exercises actual native Tissue Marshall/UnMarshall on new and reused owners with all twelve nonzero/mixed histories; an overridden setup deliberately zeroes the fields to verify restoration occurs afterward. Distinct sentinel doubles exercise exact bit preservation. Missing and malformed payloads must leave the live target untouched. The fixture initializes only base fields necessary for codec use, does not initialize a patient and does not advance physiology.

That fixture tests native owner codecs, not full-engine XML ordering. Completion requires, in order: matched combined schema/code build; C++ helper/native owner fixture; actual XML save/load with normal Tissue setup; settled no-burn, ongoing burn with changed ECF, and post-escharotomy checkpoint comparisons; next-tick resistance/event/action and fluid-creep comparisons. Include a negative legacy-burn case and rejected missing/malformed payload cases. Whole-engine replay must preserve both this payload and independently fixed sleep/Circuit/GI state; a successful isolated codec test cannot establish that integration.
