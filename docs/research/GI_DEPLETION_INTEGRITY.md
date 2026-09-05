# Known stomach-water depletion across storage units

The `whole_body_integrity_energy` baseline fails at 1199.0599999992 s because stomach water becomes −5.2939559203393777e−20 mL after the final paired transfer. At the exception, the native digestion rate is still 0.417 mL/s and timestep 0.02 s: the water value, not either rate input, violates the GI preflight.

The original failed systemic group remains under `data/derived/systemic/exertion_v2`. Independent native reproduction and source receipts are retained in `data/derived/audits/depletion-integrity-7h6dabcw/report.json`.

## Native mechanism and narrow correction

The retained GI-water source debits `digestedAmount` in mL using `SEScalarVolume::IncrementValue`. When the scalar internally stores litres, conversion of the entire pool to mL and back can overshoot by one rounding unit. A direct native example starts with 1e−9 L, reads 1.0000000000000002e−6 mL, transfers that known available amount, and leaves −2.0679515313825694e−22 mL. The next digestion call correctly rejects that negative value. The corresponding mL-stored pool empties exactly.

`whole_body_integrity_depletion` changes only the already-capped water debit. When `digestedAmount > 0` and `digestedAmount == waterContent_mL`, the paired transfer establishes an exact empty pool, so the debit uses `SetValue(0, mL)`. Partial transfers retain the original increment operation. Chyme water credit and sodium co-transfer are unchanged. The preflight still rejects arbitrary negative or missing water before transferring anything; there is no tolerance clipping, guessed correction, or general negative-state repair.

The builder reads the calcium/dry-water-corrected GI source, verifies its receipt, compiles a new GI object and replaces exactly the inherited GI-water object in the energy variant's complete 360-object link list. It verifies all parent object hashes and retains the source patch, commands, object and library hashes. Earlier libraries and donor source remain unchanged.

## Verification

```sh
.venv/bin/python scripts/build_biogears_depletion_variant.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_depletion_integrity.py
```

The verifier runs actual native mixed-unit remainder, mL remainder, partial-transfer and arbitrary-negative-input cases. It checks paired water/sodium transfers and a subsequent digestion call, then advances the original canonical baseline through 1205 s with explicit native 0.02 s steps. Parent and corrected variants execute under individually checked library linkage, with source/state/probe/dependency hashes retained. `--parent-only` and `--fixed-only` permit independent reproduction while preserving existing artifacts.

Inherited calcium, dry-GI and renal probes can select this variant with their existing `--variant` arguments. For the energy checks, the retained test runner under `data/derived/audits/depletion-inherited/verify_energy_on_depletion.py` substitutes only the corrected variant name and repository path in the existing energy verifier; production source and existing verifier files are unchanged.

This correction concerns the exact representation of an exhausted known pool. It does not validate gastric emptying rates or make a long-duration whole-body physiological accuracy claim.

## Retained acceptance

The corrected baseline completes 1205 s with exact zero stomach water. Its direct and long replay checks pass in `data/derived/audits/depletion-integrity-3ymeeb97/report.json`. Inherited checks pass: energy 26/26, dry-GI 11/11, calcium 6/6, renal 20/20. `data/derived/audits/depletion-inherited/acceptance.json` binds all six reports by hash, including the parent failure reproduction. New library SHA-256: `5b8840947b3bb964ba06f257a2874cc37efd3bd67f866565d34acbbcbb0b381d`.
