# Native dry-stomach sodium integrity

## Investigation and correction plan

The retained BioGears `Gastrointestinal::DigestNutrient` computes aqueous sodium
delivery as stomach sodium divided by stomach water, multiplied by the water
digestion rate and native timestep. It evaluates that division even when the
water pool is zero. It also invalidates the sodium scalar after complete
depletion, conflating a computed empty pool with an unknown quantity.

The native friend-access probe exercises dry empty contents, dry positive
sodium, complete depletion, ordinary hydrated transfer, hydrated zero sodium,
unknown sodium and unknown water. Every case reloads the same native saved
state, isolating the executed branch. It records stomach validity, paired chyme
credits, native exceptions and exact library/dependency hashes. The original
regression remains retained separately from the corrected result.

The new variant layers over `whole_body_integrity_renal`, preserving
calcium, renal Tmax, saturation and thermal corrections and every parent's
object identity. The scope is an aqueous sodium availability cap and zero-water
guard. No secretion, solubility coefficient or full dry-food digestive process
is fabricated. A numerical zero is permitted only when the executing transfer
establishes an empty pool; pre-existing missing values remain unknown.

## Executed native regressions

Both runs execute the same C++ friend-access probe against the same saved native
state (`ca5b165c25118c1aae3f66296cafd063a73f87dcf3f46d8d635217b18b3df8a1`).
Each case starts from a fresh load, forces the ordinary decrementing mode, and
uses the native 0.02 s timestep. The 10 mL/s water rate is an explicit test input,
not a newly calibrated biological coefficient.

| Case | Parent renal variant | Corrected variant |
|---|---|---|
| 0 g sodium, 0 mL water | NaN reaches chyme sodium; native exception | Both sodium pools remain valid; zero transfer |
| 1 g sodium, 0 mL water | Transfers 1 g without water; invalidates stomach sodium | Retains 1 g; zero transfer |
| 0.001 g sodium, 0.1 mL water | Paired depletion, but stomach sodium becomes unknown | Transfers 0.001 g and 0.1 mL; stomach sodium remains numeric zero |
| 1 g sodium, 500 mL water | Transfers 0.0004 g and 0.2 mL | Identical measured credit and stomach remainder |
| 0 g sodium, 500 mL water | Invalidates known zero sodium | Preserves known zero while transferring water |
| Unknown sodium, 500 mL water | Sodium remains unknown | Same behavior; no inferred sodium |
| Known sodium, unknown water | Poisons stomach/chyme sodium before throwing | Throws before nutrient or ion mutations |
| Zero flow, dry or very small positive water | NaN contaminates sodium | Zero transfer; valid known pools preserved |
| Zero flow, ordinary positive water | Zero transfer | Same behavior |
| Subnormal positive water, half-volume request | Transfers all sodium after concentration overflow | Transfers half, retaining the rest |

The last case uses 1e-310 mL water and 2.5e-309 mL/s deliberately as a floating
point stress test, with no physiological interpretation. Its sodium credit is
0.5000000000000248 g, reflecting representable subnormal inputs. That minute
water credit is below the receiving compartment's floating point resolution;
the probe verifies its source remainder and bounds the receiver residual, and
does not claim the receiver resolves a subnormal volume increment.

The final native library SHA-256 is
`f2f67b8bef0f8bbad25d393af9545f79f288d241802501c3a12163f26be3c162`.
The dry-GI suite passes all 11 cases; the parent fails 8 of the same cases.
The inherited native calcium suite passes 6 cases and bilateral renal suite
passes 20 cases against this library. Ordinary hydrated sodium credit matches
the parent exactly at 0.00039999999999999996 g in this probe.

## Correction and remaining limits

For known nonnegative sodium mass `M`, positive known carrier water `V`, and
finite nonnegative requested water `dV`, sodium delivery is
`M * min(1, dV / V)`; zero `V` or zero `dV` gives zero sodium delivery. This is
the existing hydrated kinetic law with a bounded water fraction, avoiding
concentration overflow followed by multiplication by zero. A verified complete
transfer sets sodium to zero in grams atomically. Pre-existing unknown sodium
remains unknown. Unknown, negative or nonfinite water, invalid water rate or
timestep, and known negative/nonfinite sodium are rejected before this method
mutates nutrient pools. The unknown-water regression includes 10 g carbohydrate
to detect premature nutrient mutation.

The source's `m_DecrementNutrients=false` stabilization reservoir mode is
preserved. Its source pool does not deplete and its water is not credited into
the receiving compartment; that mode is not a closed-body conservation claim.
The transfer-budget probes explicitly use decrementing mode. Existing protein,
fat, carbohydrate, gastric secretion and downstream absorption behavior is not
replaced or calibrated here. In particular, preventing sodium transfer without
carrier water does not establish realistic handling of dry food, dissolution,
secretions, osmotic transport or complete digestion.

## Reproduction and retained data

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_biogears_dry_gi_variant.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_dry_gi.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_dry_gi.py --variant whole_body_integrity_renal
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_gi_integrity.py --variant whole_body_integrity_gi_water
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_renal_integrity.py --variant whole_body_integrity_gi_water
```

The parent dry-GI command is expected to exit 1. Native stdout, errors, executable,
probe source, exact state identity, and resolved library/dependency hashes are
retained under `data/derived/audits/dry-gi-integrity/<variant>/`. The initial
seven-case negative probe is separately retained in
`whole_body_integrity_renal_initial_seven`. Calcium and renal receipts are under
their respective `gi-integrity` and `renal-integrity` audit directories.

The generated variant manifest retains the complete parent object hash map,
new object map, donor revision and source hash, calcium-corrected source hash,
replaced object hash, compiler and link commands, full patch, and resulting
library hash. The builder checks all parent objects before and after compilation,
and checks donor and parent library bytes remain unchanged. Only its copied GI
source/object and new variant library are rebuilt. This is a numerical/native
data integrity correction, with no new biological calibration claim.
