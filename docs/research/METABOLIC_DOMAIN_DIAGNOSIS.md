# Native muscle metabolic domain failure: source diagnosis

2026-09-05. Read-only investigation of `data/derived/audits/embodied-native-kdfvcq_a/failure.json`, which records `Out-of-domain native muscle metabolism` after 1.948446433 s wall time. No native process, compilation, model change or guard weakening was performed for this diagnosis.

**Conclusion:** signed negative active-fiber mechanical power is valid under the selected native model. Negative **total metabolic power** is not its intended output with the enabled options. However, the held source's nonnegative-power correction can leave a tiny negative floating-point cancellation remainder that the adapter's exact `total < 0` guard rejects. That is a concrete candidate mechanism, not yet the measured cause of this event. The receipt lacks the offending muscle, power, mass and simulation time.

## Bound identities and configuration

The failed execution receipt `body/mechanics/native/execution.json` binds adapter build `data/runtime/mechanical-stream/build-09xtbmes`. Its frozen and current `native_muscle_metabolism.h` both hash to `0eabf9a9268f2cd5e8a30a79c3c4e548a8d40cca8ff631e2838e980c3627763d`; the executable hash is `351fd9e0578e67595401bc0ef4e1aae0a83159fcbe5eccd456b7ea829a3992a4`. The receipt records `libosimSimulation.so` hash `3a041f37914bbfece5ff0c7409220977cfa48d25ff5086cabe8c5456b98e3e66` and the remaining native dependencies.

Inspected held sources under `data/raw/mechanics/opensim-core/OpenSim/Simulation/Model/`:

| File | SHA-256 |
| --- | --- |
| Umberger2010MuscleMetabolicsProbe.cpp | `4a69981e2fdb6c9fc7625cabc1d170280988f45ab1dd2f98157d69829375c879` |
| Umberger2010MuscleMetabolicsProbe.h | `14563d98023821934db773cddc6faf5dab2a5e59b64f5d3e04ed35218f3aa07e` |

These source hashes identify the source inspected; they are not independently reconstructed object provenance for the previously linked OpenSim library.

The retained `assembled_model.osim` confirms 92 metabolic muscle parameters, `probe_operation=value`, gain=1, activation/maintenance=true, shortening=true, basal=false, work=true, minimum-heat=true, include-negative-work=true and forbid-negative-total=true. Thus the constructor's boolean order is correct; the third false disables whole-body basal rate and does not disable shortening heat. The adapter uses effort scaling 1 and Bhargava recruitment=true.

All 92 retained parameters have `use_provided_muscle_mass=false`. The native formula at cpp:821–831 is `Fmax / specific_tension * density * optimal_fiber_length`; defaults at cpp:852–853 are 250000 N/m² and 1059.7 kg/m³. Applying this formula to the retained XML gives masses from **0.07381873423371656 kg** (`glmin3_l`) to **2.589550403664654 kg** (`vaslat_r`), all positive. This is a static calculation, not a native `getMuscleMass` trace. The adapter also rejects nonpositive Fmax and optimal length during construction. Mass failure is therefore less supported by the static evidence than negative total power, but the combined exception does not itself identify which predicate fired.

## Signed work versus total metabolic power

At cpp:396–403 the native work term is `Wdot = -active_force * fiber_velocity / muscle_mass`, after flooring negative active force to zero. Positive fiber velocity denotes lengthening, so Wdot can be negative. The header's mechanical-work documentation at h:153–176 expressly allows this absorbed work and explains why it increases shortening heat when necessary to avoid negative metabolic demand.

With the adapter's settings, the native path computes activation/maintenance heat AMdot and shortening heat Sdot. At cpp:409–412 it evaluates:

```cpp
const double Edot_Wkg_beforeClamp = AMdot + Sdot + Wdot;
if (Edot_Wkg_beforeClamp < 0)
    Sdot -= Edot_Wkg_beforeClamp;
```

It then recomputes `totalHeatRate = AMdot + Sdot`, applies the optional 1 W/kg heat floor, adds Wdot, and multiplies by muscle mass (cpp:430–455). There is no final exact-zero operation after this recombination. In exact arithmetic, finite correctly formed terms and positive mass imply nonnegative total. In floating-point arithmetic, cancellation can leave a negative remainder.

A small Python binary64 arithmetic check, using the same operation order but **not executing the native probe**, demonstrates the issue:

| Quantity, W/kg | Binary64 value |
| --- | --- |
| AMdot | 0.1 |
| Sdot before correction | 0.1 |
| Wdot | −100.1 |
| `AMdot + Sdot + Wdot` | −99.89999999999999 |
| Sdot after subtraction | 99.99999999999999 |
| Recomputed total | **−1.4210854715202004e−14** |

The heat floor does not activate here because recomputed heat is about 100.1 W/kg. Multiplying by a positive mass preserves the tiny negative sign. This synthetic example demonstrates numerical susceptibility of the source expression; it does not establish that these AM/S/W terms correspond to any actual retained muscle or explain the failing native sample.

`native_muscle_metabolism.h:74–79` checks individual total, active force, fiber velocity and mass for finiteness before throwing on `total<0 || mass<=0`. Therefore this particular exception implies the values passing those preceding checks were finite. The adapter subsequently calculates signed work and heat as `total-work`; it does not reject negative work. Changing the signed-work convention would not address the guard's actual predicate and would change the energy model.

The generic stream calls the same sampler before integration, after integration, and when emitting snapshots. The retained failure does not identify which sample raised the exception or its simulated time. The command handler restores its pre-command state on exceptions, so no successful advanced state should be inferred from the failed receipt.

## Smallest discriminating observation

Preserve the strict guard until it reports the evidence. In a separately identified diagnostic build, expand only its exception/observation payload with:

- Muscle name, probe output label/index, simulation time and sampling phase.
- Individual total W and analysis mass kg at full round-trip precision, plus probe aggregate total and basal values.
- Native active force, fiber velocity, signed work W, activation, excitation/control, normalized fiber length and active force-length multiplier.
- Enabled probe flags and, if accessible through a diagnostic source copy, AMdot, Sdot before/after native correction, Wdot and final W/kg before multiplication.

Repeat the smallest already failing articulated command with identical held model, controls and inputs; no physiology advance is needed to diagnose this guard. A healthy preceding sample plus a constant positive measured mass would isolate the power predicate. Do not infer the branch from aggregate total alone.

If the measured negative is a roundoff-sized recombination residue, justify a tolerance against the scale of the contributing heat/work terms and the floating-point operations, retain raw and accepted power separately, and bound the resulting energy-ledger correction. Preserve rejection for materially negative or nonfinite total, and for nonpositive mass. Require a native regression at the original failing state and aggregate/per-muscle energy closure before adopting such a tolerance. An arbitrary fixed clamp, disabling `forbid_negative_total_power`, or interpreting negative total as ordinary signed mechanical work is not supported by this source diagnosis.

If the measured value is materially negative, inspect the actual loaded source/build, per-muscle label mapping, enabled terms and input state before altering acceptance rules. The current observations do not license a metabolic calibration or a claim that the articulated physiology is valid.

## Follow-up: instrumented native value and bounded residual rule

The coordinating native run subsequently reported diagnostic build `build-l2gezl04`, receipt `data/derived/supine-support-tt0k8jes/report.json`, and an actual failing `ehl_l` sample. Direct inspection of that receipt confirms last successful time 0.005 s, failed stage target 0.02 s, and these native diagnostic values:

| Field | Value |
| --- | --- |
| Total metabolic power, W | −9.4635567759726085e−17 |
| Analysis mass, kg | 0.10655017692467443 |
| Active fiber force, N | 3.6009271834261525 |
| Fiber velocity, m/s | 0.16601119770192224 |
| Activation and excitation | 0.01 each |

These values isolate the negative-total predicate with a positive mass. Independently recomputing the adapter's signed work gives **−0.597794234557985 W**; `total−work` gives heat **0.5977942345579849 W**. The total is about 0.713 binary64 epsilon times the heat/work scale, consistent with the source recombination mechanism. AMdot and Sdot at this actual native event were not separately observed here, so the synthetic algebra example above must not be relabeled as those actual reaction terms.

A proposed acceptance bound, preserving the raw native output, is:

```text
scale_W = max(abs(signed_work_W), abs(raw_total_W − signed_work_W),
              analysis_mass_kg * 1 W/kg)
tolerance_W = 32 * binary64_epsilon * scale_W
accept only finite raw_total_W >= −tolerance_W and finite analysis_mass_kg > 0
```

The mass factor supplies the native minimum-heat scale with correct units; it is not a universal 1 W absolute tolerance. The factor 32 is a conservative numerical operation margin around the source's correction/recombination, not a physiological parameter. It should be validated against retained native regression cases and tightened if evidence permits. For this sample the bound is **4.247583508444457e−15 W**, about 45 times the observed residual. A −1e−6 W result at the same mass/work scale remains rejected.

Do not overwrite or clamp the raw per-muscle total. Add explicit audit count and signed sum for accepted negative residuals, along with each residual's numerical bound; preserve the raw aggregate/per-muscle sum ledger. If energy is accumulated over time, quantify the corresponding residual contribution using the same integration ownership. Finite checks and the positive-mass requirement remain unconditional. The tolerance addresses a numerical-domain predicate; it does not validate the muscle model or whole-body physiological coupling.

A source-order binary64 regression was executed as a small pure Python check. It follows the original AM/S/W correction, minimum-heat step and final mass multiplication. With AM=.1, S=.1 and W=−100.1 W/kg, and the positive mass above, the original source-order algebra yields **−1.5141690841556174e−15 W**: the original exact-sign predicate fails, while the bounded residual predicate accepts. Six additional predicate checks passed: accept the measured native residual, reject −1e−6 W at its scale, accept zero, reject zero mass, reject NaN and reject negative infinity. This verifies the proposed acceptance algebra, not a changed native executable. The regression intentionally leaves the total untouched:

```python
import math, sys

def accept(total, mass, work):
    if not all(map(math.isfinite, [total, mass, work])) or mass <= 0:
        return False
    scale = max(abs(work), abs(total - work), mass * 1.0)
    return total >= -32 * sys.float_info.epsilon * scale

am, short, wkg = .1, .1, -100.1
before = am + short + wkg
if before < 0:
    short -= before
heat = am + short
if heat < 1.0:
    heat = 1.0
mass = .10655017692467443
raw_total = (heat + wkg) * mass
assert raw_total < 0 and accept(raw_total, mass, wkg * mass)
assert accept(-9.4635567759726085e-17, mass,
              -3.6009271834261525 * .16601119770192224)
assert not accept(-1e-6, mass, -.597794234557985)
assert accept(0., mass, -.597794234557985)
assert not accept(raw_total, 0., wkg * mass)
assert not accept(float('nan'), mass, wkg * mass)
assert not accept(float('-inf'), mass, wkg * mass)
```

The coordinating implementation subsequently selected the more conservative bound `64 * epsilon * max(1 W, abs(work_W), mass_kg * 1 W/kg)` and exposed per-muscle tolerances while preserving raw totals. That implemented policy is distinct from the narrower, mass-scaled 32-epsilon recommendation above. The coordinator reports that it removed this guard failure and exposed the next independent signed-demand interface failure. This document does not independently claim full native regression completion for that later implementation.
