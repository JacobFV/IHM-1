# Isolated native gastrointestinal calcium integrity

2026-09-05. The `whole_body_integrity` library fixes the retained native calcium stomach-to-chyme branch while preserving the saturation-bounds, evaporative telemetry and thermal-resistance corrections in `saturation_bounds_heatflux_thermal_units`. It does not establish whole-body nutrient or calcium conservation.

## Defect and correction

At donor revision `3f16a5fa1dade9c511b88d923606fa51cc35e95d`, `Gastrointestinal::DigestNutrient()` computes calcium transfer in mg, multiplies by the configured calcium absorption fraction, then debits stomach using g and credits chyme using mg. Its strict comparison also strands remainders smaller than or equal to the proposed transfer.

The isolated patch retains the native rate and absorption-fraction multiplier. It caps the proposal between zero and available nonnegative stomach calcium, and uses mg for both endpoints. A positive remainder is transferred completely. A pre-existing negative stomach value remains negative and credits nothing: correcting historical state would require a separately justified initialization procedure. The macro rate expressions, including their commented `/60` conversions, remain byte-for-byte unchanged. Their physiological calibration remains unresolved. Absorption fraction continues to scale transfer speed here; the patch does not reinterpret it as fecal loss or implement missing stool/colon physiology.

## Native regression

`scripts/verify_native_gi_integrity.py` compiles a small C++ probe against retained headers. It loads the existing `native_baseline_v2/states/native_stabilized.xml`, then uses upstream's existing `BioGearsEngineTest` friend access to set only test stomach calcium and call the real library's `DigestNutrient()`. This exercises real native scalars, quantities, units and balancing. The source timestep is 0.02 s; configured 2.7 mg/min and fraction 0.25 produce 0.000225 mg per call. Other digestion branches execute; absorption and whole-engine advancement are deliberately outside this local branch experiment. Chyme credits are measured directly before and after the call, not inferred from blood storage.

| Initial stomach mg | Original final mg | Corrected final mg | Corrected chyme credit mg |
|---:|---:|---:|---:|
| 500 | 499.775 | 499.999775 | 0.000225 |
| 0.001 | -0.224 | 0.000775 | 0.000225 |
| 0.0001 | 0.0001 | 0 | 0.0001 |
| 1e-14 | 1e-14 | 0 | approximately 1e-14 |
| 0 | 0 | 0 | 0 |
| -0.175000000020085 | unchanged | unchanged | 0 |

Original upstream and prior thermal variant regressions exit 1 with four failed cases; corrected regression exits 0 with six passing cases. The original failing regression was executed before the correction was implemented. Assertions require finite outputs, matching expected debit and credit within 1e-12 mg absolute tolerance, nonnegative output from nonnegative input, and exact zero when a positive remainder is exhausted. The negative fixture asserts no transfer and no silent repair. Tiny chyme differences show ordinary floating-point cancellation.

The retained `native_hour_rest/states/native_final.xml` independently contains stomach calcium -0.175000000020085 mg (SHA-256 `e9982d93dcab94f884288ee1e3c94f3fbc366ea74c289f2f04518417918a4337`). The fixture reproduces this value; that historical file is not edited. Initial regression state SHA-256: `42aa7af3ec40a4bfab582a6e2658d1869e28545cd9c831d709a89599e00117b1`.

## Build and provenance

Run from the repository:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_gi_integrity.py --variant upstream
# Expected regression failure above.
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_biogears_gi_variant.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_gi_integrity.py
```

The builder checks revision and GI source equality against the committed donor blob, verifies prior library/source manifest digests, compiles a new GI object, then relinks the 360-object parent response list with exactly that object replaced. It checks that prior engine Saturation and Environment replacement objects remain selected. The donor source, donor objects and older variant are unchanged. This requires the retained build tree and parent variant objects; it is not a fresh upstream bootstrap.

Generated `data/runtime/physiology/variants/whole_body_integrity/manifest.json` records exact source/patch digests and patch text, compiler/compile/link commands, original GI object hash, every linked object hash, parent manifest/library hashes and resulting library hash. Runtime audit reports record initial state, compiled probe, resolved dependency hashes and explicit `LD_LIBRARY_PATH`; `ldd` resolution is checked against the requested variant.

| Artifact | SHA-256 |
|---|---|
| New library | `68b6188236819d45f61f828051f88f4497b1f2ad45416cfaa889ed8724944a3e` |
| GI original source | `7a3f00270e1a18e3402cf0d4276f0143b39a554695edc5dbe78fd692e0278758` |
| GI patched source | `20ed9995adf3b9aff43a850670ebd9dd9b4c43c533605afdb932ab3f0f1cfbcb` |
| GI patch | `6765e69b7c5ff5ed463ca7ad421b057babd037a8bfa2567c789285ef9420d351` |
| Parent library | `ec9a52913887de4d823b34d7b328ea43fdf379324d105c9ef370f041bbd0665a` |

Native stdout/stderr, compiled probe source/binary, `ldd.txt`, JSON regression reports and build logs are retained locally under `data/derived/audits/gi-integrity/`. These generated artifacts are not committed. Systemic calcium absorption/removal, long meal trajectories, other nutrient nonnegativity and global elemental or energy conservation need separate checks. A library name is not a physiological validation claim.
