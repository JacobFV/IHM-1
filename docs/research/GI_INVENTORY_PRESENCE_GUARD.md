# GI inventory presence and continuation guard

Source preparation only; no new compilation or patient execution. This extends the accepted isolated three-TU codec repair without changing its original receipts or pretending it distinguished historical omissions.

`prepare_gi_inventory_presence.py` provides composable `patch_header`, `patch_gi`, and `patch_io` functions. The composition lead owns the matching optional `GITransitInventoryVersion` (`xs:unsignedInt`, at most one) schema field, schema generation, dependency closure, and builds. Missing is parseable at schema level so the native loader can explicitly reject it; version 1 declares the saved inventory authoritative even when its record list is empty.

| Owner transition | Required behavior |
| --- | --- |
| Constructor / Invalidate | Inventory initialized flag false |
| Actual native fresh Initialize | Require empty drug map; mark true after initialization succeeds |
| Continuation decode | Reject absent/non-1 marker before mutation; invalidate; decode every record; mark true only after successful completion |
| Encode unknown owner | Throw; do not turn an unrecorded empty map into an authoritative empty map |
| Encode initialized owner | Write all repaired drug records and marker 1 |
| Legacy fresh-empty preparation | Require explicit declaration/provenance; reject existing transit records, oral actions, markers, or ambiguous GI owner; insert marker and retain input/output hashes |

The runtime boolean is in the derived native Gastrointestinal owner, not the base CDM owner. Its addition changes object layout. Rebuild all dependent translation units, including allocation sites, against one composed header; do not mix an old constructor/object layout with the new IO code. Apply the GI.cpp transform to the accepted absorption-corrected source (`804cba46cea17e9b05be3c45c3d20b2e4683dfdb12a5ea8ea7016de415cb6235`), not pristine upstream, so this guard cannot discard the donor correction. The IO transform edits only GI reader/writer bodies and requires the accepted outer record writer already present; sleep and burn changes in the same translation unit remain intact.

Explicit fresh legacy declaration JSON has exactly these keys:

```json
{"schema":"ihm.explicit-fresh-gi-inventory.v1","declaration":"fresh_no_prior_oral_drug_inventory","provenance":"A concrete authorized initialization rationale"}
```

This declaration starts a new declared history. It does not establish that the old state's missing records were historically empty, and absence of an oral action is not evidence of absence of prior doses. The migration checks reject contradictory visible records but cannot recover omitted mass. It inserts only the marker bytes into one explicit GI XML element, preserving unrelated scalar text; output and `.initialization.json` paths must be new. There is no raw-load compatibility bypass and no automatic migration during save.

Light verification: `scripts/verify_gi_inventory_presence.py` first failed because the new transformer did not exist, then passed source transform/idempotence checks, explicit declaration validation, byte-local XML insertion, and rejection of missing declarations, blank provenance, existing drug records/actions, and duplicate markers. This verifies preparation logic only. Native acceptance remains required: unknown owner save rejects; missing/0/2 marker loads reject; declared empty reload remains empty; populated 9/9/8 vectors and separate cumulative totals survive public save/load; failed decode cannot be saved as initialized; accepted continuation advances matching next ticks. Keep legacy rejection and fresh-declaration receipts separate from the populated persistence receipt.
