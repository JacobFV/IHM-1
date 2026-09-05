# Isolated native renal transport-maximum integrity

2026-09-05. `whole_body_integrity_renal` corrects `Renal::CalculateReabsorptionTransport()` from retained BioGears revision `3f16a5fa1dade9c511b88d923606fa51cc35e95d`. It is a new library layered over `whole_body_integrity`; donor source and earlier variant bytes are preserved.

## Executing defect and exact change

The source first computes a proposed tubule-to-peritubular mass, caps it at available tubular mass, and derives its rate by dividing by `m_dt` (seconds). For a finite transport maximum it caps only that rate. Both kidney compartment transfers and the saved glucose reabsorption mass still use the larger uncapped amount. Thus reported rate and actual transferred mass disagree.

The correction applies the existing maximum in mass units before downstream bookkeeping:

```cpp
massToMove_mg = std::min(massToMove_mg, transportMaximum_mg_Per_s * m_dt);
reabsorptionRate_mg_Per_s = massToMove_mg / m_dt;
```

The prior available-mass cap remains in place. The same resulting `massToMove_mg` is debited from tubules, credited to that kidney's peritubular compartment, and stored in `leftGlucoseReabsorptionMass_mg` or `rightGlucoseReabsorptionMass_mg`. Their sum divided by the timestep agrees with the reported total reabsorption rate. No other source statements change.

Native semantics remain: finite-ratio proposals use concentration (mg/mL), next flow (mL/s), dt, ratio and permeability modification; negative flow skips transport. Infinite reabsorption ratio proposes all available mass even at zero flow; this source behavior is not changed. Infinite transport maximum bypasses the cap. The source applies the configured transport maximum independently to each kidney; no inferred division by two is introduced. Valid nonnegative quantities/configuration and positive native timestep are the scope of this correction, not a new invalid-state sanitization policy.

## Native regression evidence

`scripts/verify_native_renal_integrity.py` compiles a native C++ friend-access probe and loads the existing `native_baseline_v2/states/native_stabilized.xml` (SHA-256 `42aa7af3ec40a4bfab582a6e2658d1869e28545cd9c831d709a89599e00117b1`). It tests the actual dynamically linked renal routine, scalars, balancing and compartment quantities. The test overrides only its declared local test inputs: glucose masses, transport maximum, reabsorption ratio, permeability factor and next reabsorption flow. Test-only read-only flags are cleared to set the native test controls. Glucose bookkeeping is reset as in the source active-transport entry point. Left/right initial masses differ to verify bilateral handling.

The original failure was established before writing the correction. At dt 0.02 s:

| Case | Left/right proposed mass mg | Permitted mg per kidney | Original transferred mg | Corrected transferred mg |
|---|---:|---:|---:|---:|
| Above cap: 1 mg/s maximum | 10 / 5 | 0.02 | 10 / 5 | 0.02 / 0.02 |
| Zero maximum | 10 / 5 | 0 | 10 / 5 | 0 / 0 |
| Finite ratio: 0.001 mg/s maximum | 0.016259770839691 / 0.008129885419846 | 0.00002 | proposed amounts | 0.00002 / 0.00002 |
| Below cap: 1 mg/s maximum | 0.01 / 0.005 | 0.02 | proposed amounts | proposed amounts |

Additional cases cover zero mass, infinite maximum, zero ratio, zero flow with finite ratio, backflow, and proposals bounded by available mass. There are ten cases and two kidneys: original native execution exits 0 but verifier exits 1, with six of twenty rows failing; corrected native execution and verifier both exit 0, twenty of twenty rows passing. Each row checks finite output, nonnegative remaining mass, availability bound, equal actual debit/credit, stored glucose mass and summed reported rate, within 1e-10 mg absolute tolerance (ordinary scalar subtraction can incur cancellation). The unchanged inherited GI regression passes six of six cases using this new library.

## Downstream bookkeeping and limits

`Renal::Gluconeogenesis()` reads the saved left/right glucose reabsorption mass and combines it with lactate excretion before another native glucose transport-maximum calculation. An uncapped stored glucose amount can make its computed lactate conversion negative. The patch ensures this input represents the actual capped glucose transfer. The regression checks this input directly; it does not execute or validate the complete gluconeogenesis pathway. Its other lactate/excretion accounting is unchanged and remains separately auditable.

This experiment invokes the reabsorption branch without complete engine advancement. It proves the local paired mass cap and matching bookkeeping, not renal physiology calibration, glycemic thresholds, full urine mass balance, whole-body elemental conservation or long-run stability. The native glucose substance's transport parameters are not independently validated by this numerical correction. Source glomerular reabsorption and other unrelated equations are untouched.

## Reproduction and provenance

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_renal_integrity.py --variant whole_body_integrity
# Expected regression failure above.
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/build_biogears_renal_variant.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_renal_integrity.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_gi_integrity.py --variant whole_body_integrity_renal
```

Builder verifies the donor committed blob, parent manifest/library/source and every parent object hash. It compiles a new renal object and relinks the parent's 360-object response list replacing only engine Renal. Saturation, Environment and GI replacement objects remain selected. Generated manifest records exact patch text/digest, source digests, commands/compiler, original renal object, every selected object and resulting library hashes. This depends on the retained native build tree and parent variant; it is not a bootstrap build.

| Artifact | SHA-256 |
|---|---|
| New library | `01205c01fb96bf711a975e937886db1ff7b8fa1dae913711e97c294aa41e19de` |
| Renal original source | `2eeb59b030def647a51497b724c2b078f01e26aac8c0ce1afd415dc3a3d2c0fb` |
| Renal patched source | `4306870b72afd2a98109826733365b6643df1775a153a53c96cf6eae1ab39a65` |
| Renal patch | `daea5348dabe5b01a00e1919d3d0e84a4ffbb7c6c856e136c2a9776cccb3f66e` |
| Parent GI library | `68b6188236819d45f61f828051f88f4497b1f2ad45416cfaa889ed8724944a3e` |

Generated variant files are in `data/runtime/physiology/variants/whole_body_integrity_renal/`. Native probe source/binary, stdout/stderr, `ldd.txt`, JSON reports and build logs are retained under `data/derived/audits/renal-integrity/`; inherited calcium evidence is under `data/derived/audits/gi-integrity/whole_body_integrity_renal/`. Runtime reports record library/state/probe/dependency hashes and explicit library search path; the verifier checks that `ldd` resolves the requested library. Generated artifacts remain local.
