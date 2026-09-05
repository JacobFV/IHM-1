# Native meal counterion audit

The native nutrition model can transfer sodium without any represented counterion. This affects meal acid–base trajectories even when the added meal has no sodium: nutrient cotransport mobilizes sodium already present in the initial stomach. The current meal response must remain a qualified research result, not normal calibrated digestion.

Actual five-case evidence: `data/derived/audits/meal-electrolyte-native-02rzs870`. Independent retained audit: `data/derived/audits/meal-counterion-audit-y1aj1a8q/audit.json`. All cases execute 600 seconds at the native 20 ms step, recording the first step and then every ten seconds. They share final-RH library `ab485812968b17e16915b04a6e5144bb1e1ec7f5bf78213936c3334d087c0a5a` and fresh-state SHA `cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3`. No source equations, initial pools or pH values were edited.

The initial stomach contains 500 mL water, 1 g sodium and 500 mg calcium. Initial intestinal chyme sodium and chloride are both zero. Consequently “rest” here is not an empty-gut control. Water-only adds 500 mL; sodium-only adds 1 g sodium using the baseline carrier water; nutrients-only adds 60 g carbohydrate, 20 g protein and 20 g fat; sodium+water adds 1 g and 500 mL. No added calcium or chloride in any case.

| Case | Final arterial pH | Final VC SID (mmol/L) | Net GI→vascular Na (g) | Final chyme Na (g) |
|---|---:|---:|---:|---:|
| Rest | 7.411568 | 40.531028 | 0.019467 | 0.480933 |
| Water only | 7.410610 | 40.489675 | 0.009734 | 0.240466 |
| Sodium only | 7.412979 | 40.613730 | 0.038934 | 0.961866 |
| Nutrients only | 7.446243 | 42.457903 | 0.500388 | 0.000012464 |
| Sodium + water | 7.411594 | 40.531028 | 0.019467 | 0.480933 |

Net GI sodium is calculated from initial stomach+chyme sodium, declared addition, and final stomach+chyme sodium. It includes secretion minus reabsorption; it is not gross transporter throughput or a whole-body charge ledger. All local recorded masses remain finite and nonnegative. Native SID and independently reconstructed `VC.Na + VC.K − VC.Cl − VC.lactate − 1.02` agree within 1e−10 mmol/L. The source’s stale nearby comment describes −5.3, but the executing constant is −1.02.

Nutrients-only versus rest increases vena-cava sodium by 1.74798 mmol/L and decreases chloride by 0.16899 mmol/L. Arterial CO₂ is 39.7303 versus 40.3497 mmHg; bicarbonate is 27.6817 versus 25.9559 mmol/L. Thus the sodium/SID pathway is demonstrated, while the pH change includes respiratory and buffer responses too. Chyme chloride stays exactly zero throughout every case. Sodium+water closely matches rest at ten minutes because the gastric sodium/water concentration is unchanged and the larger carrier reservoir has not yet emptied.

## Source ownership and missing chemistry

Pinned donor revision `3f16a5fa1dade9c511b88d923606fa51cc35e95d`; retained executing GI translation unit is `data/runtime/physiology/variants/whole_body_integrity_depletion/Gastrointestinal.cpp`, inherited unchanged by the final-RH library.

* `SENutrition.h` exposes carbohydrate, fat, protein, calcium, sodium and water. It has no chloride, potassium, salt identity or counterion field. The native action bridge preserves this limitation.
* `DigestNutrient` transfers gastric sodium to chyme with water; `ChymeSecretion` transfers vascular sodium to chyme; `AbsorbNutrients` transfers sodium alongside carbohydrate/amino acids, or alone with carrier water when nutrients are absent. These are paired sodium-mass transfers. None transports chloride. There is explicitly no chyme–vascular substance graph link to advect all ions with water.
* `ChymeSecretion` itself comments that other ions and substances are needed for correct pH balance. Its “salt” comment is implemented as sodium alone. Sodium secretion and absorption change distinct owner pools; this is not evidence of duplicated sodium mass.
* `BloodChemistry` publishes the above vena-cava SID. `Saturation` uses strong-ion balance in its acid–base solve. Calcium is not a dynamic term in that published SID expression; the fixed other-ion offset does not constitute a measured complete electrolyte budget.

The held sodium molar mass is 22.9898 g/mol: 1 g elemental sodium corresponds to 43.49756 mmol positive equivalents. It does **not** identify sodium chloride, bicarbonate, citrate or another salt. Inferring chloride would silently select a chemistry and acid–base load. Likewise, a nutrition calcium mass does not identify its counterion. This audit does not claim literal bulk Coulomb charge accumulation: buffers, bicarbonate and proton balances adjust inside the model, and a full charge ledger was not instrumented.

Required next schema: explicit intake chemical species, amount/unit, valence and counterion/compound identity; paired luminal secretion, absorption and excretion ownership; electroneutrality/buffer terms and unmodeled charge made explicit; no automatic Na→NaCl conversion. Existing nutrition inputs must remain qualified as incomplete chemistry. A replacement needs source-supported transporter stoichiometry and matched sodium-salt/counterion perturbations, not a pH reset or arbitrary chloride fit. These probes isolate a missing pathway; they do not establish long-term calibration.

## Reproduction and evidence limits

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/run_native_meal_electrolytes.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_native_meal_electrolytes.py data/derived/audits/meal-electrolyte-native-02rzs870
```

The custom diagnostic archives source, executable, library, actual `ldd` dependencies, initial state and runtime resources before execution and detaches its resource tree. The independent audit verifies archive bytes and detached-tree hashes. It predates the stronger `NativeSession` execution-selection receipt and therefore does not have `execution-inputs.json`; actual file-open/loader access was not traced. This record is diagnostic evidence, not a default-publication candidate. Failed audit attempts remain in separate directories; accepted audit and native outputs are not overwritten.
