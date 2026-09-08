# Long-horizon native thermal audit

The historical six-hour meal and hydration experiments are **not normal-rest
physiology**. Their native core temperatures reach 24.16 and 24.19 °C. Local
mass conservation and intervention contrasts did not establish physiological
validity. These outputs remain retained as failed thermal experiments.

This audit changes neither donor source nor historical states. It builds two
isolated, source-attached corrections and tests the actual native circuit.
No clothing level, environment temperature, controller gain, or subject
parameter was fitted to obtain a desired core temperature.

## Inputs and initialization

BioGears donor revision is `3f16a5fa1dade9c511b88d923606fa51cc35e95d`.
The generic male profile is 70.7713 kg (77.1107029 kg when this audit was written), 1.7194712 m, age 44, fat fraction
0.21; these are generic priors, not measurements of the user. Native skin
area is 1.9012784297 m². The actual environment is 22 °C air and mean radiant
temperature, 0.1 m/s air speed, 60% relative humidity and 0.5 clo.

The native radiation correlation assumes standing. Naming the canonical
reference pose supine does not introduce mattress conduction, bedding,
posture-dependent exposure, or bed contact into this circuit.

`native_baseline_v1` was initialized using `saturation_bounds_heatflux`, before
the later heat-transfer unit correction. Its core thermal node is
36.5750397 °C, already losing 258.64 W from core storage and 278.79 W from
total storage. Its reported output core temperature lags the thermal node by
about 21 µK. Loading it under a different constitutive variant is an explicit
initial-state/variant mismatch. A fresh initialization under the later energy
variant still loses 488.63 W: core 36.02955 °C, metabolic plus fixed skin
sources 116.90 W, boundary loss 605.53 W. Reinitialization alone did not fix
the boundary.

The circuit heat capacities total 258118.13 J/K, including core 243683.33 J/K.
The audit counts unique thermal paths, not duplicate compartment rollups, and
converts the skin sources from kcal/hour to watts. Instantaneous storage equals
metabolic input minus boundary loss to within micro-watts in these states.
The excessive cooling was encoded in the circuit; it was not unaccounted
numerical energy loss.

## Clothing and surface-film correction

The source documentation defines ensemble clothing insulation as a whole-body
average. See [BioGears environment methodology](https://www.biogearsengine.com/documentation/_environment_methodology.html)
and [CDM environment parameters](https://www.biogearsengine.com/documentation/_c_d_m_tables.html).
The following equivalent-circuit result follows directly from that definition.

For total area A and regional area fractions fᵢ summing to one, uniform local
clothing insulation I in m² K/W requires regional resistances
Rᵢ = I/(A fᵢ). At equal skin temperatures the six parallel branches therefore
have equivalent resistance I/A. One clo converts to 0.155 m² K/W.

The inherited code multiplied I/A by each regional fraction instead. For its
fractions and input 0.5 clo, the parallel equivalent represents 0.00332084 clo.
It additionally multiplied convective film resistance by `0.1 * clo` and
radiative film resistance by `5 * clo`, although clothing is already a separate
series layer. At 0.5 clo these multiply nominal convective conductance by 20
and radiative resistance by 2.5.

The isolated correction uses Rᵢ = I/(A fᵢ), retaining the source skin-area
fractions 0.36, 0.07, 0.092, 0.092, 0.193, 0.193. Film resistances are
1/(h A). The source heat-transfer correlations remain unchanged. Uniform
regional insulation is an explicit materialization prior, not a measured
regional clothing map.

`native_thermal_boundary.h` is used both by the compiled native patch and the
algebra tests. Tests cover unequal areas, subdivision invariance, the known
parallel/series equivalent, two boundary temperatures, heat balance, and zero
clothing. The original circuit short resistance 1e−100 K/W produced a native
NaN when clothing was zero. Version 2 retains a documented 1e−8 K/W numerical
short floor. This is not physical clothing: it changes resting film resistance
by less than 2e−7 relative. Compared with a smaller positive-insulation probe,
the two-second native core difference is 6.7 nK and boundary heat difference
32 µW. At 0.5 clo, versions 1 and 2 have identical sampled outputs.

## Total skin perfusion was counted six times

The producer/consumer chain in the pinned source is decisive:

1. `Cardiovascular::SetUp` binds `m_pAortaToSkin` to the single
   `Aorta1ToSkin1` circuit path.
2. `CalculateVitalSigns` samples that path's flow into
   `m_CardiacCycleSkinFlow_mL_Per_s`.
3. The cardiac-cycle boundary publishes the accumulator mean as
   `GetMeanSkinFlow()`, then resets it.
4. `Energy::UpdateHeatResistance` partitions this value over six thermal
   regions, but multiplies every regional share by six.

The mean is temporal, of total skin inflow. It is not a spatial mean of six
regional flows. In the actual saved native state, inverting the six thermal
resistances gives a summed regional flow exactly six times native skin inflow.
The new variant removes only those six multipliers. It retains source regional
weights (2.10, 1.44, 1.25, 1.25, 2.925, 2.925)/11.89, blood density and heat
capacity, the empirical transfer fraction alpha = 0.42, and resistance bounds.
The corrected saved state gives a sum ratio 1.0000000000000002.

This proves the partition identity; it does not independently validate the
empirical alpha or regional vascular distribution against human measurements.

## Native results and retained identities

| Variant / condition | Initial core °C | One-hour core °C | Final core slope °C/hour | Final storage W |
|---|---:|---:|---:|---:|
| Source boundary, fresh initialization | 36.02955 | 30.96511 | −3.89280 | −265.46546 |
| Corrected boundary only, fresh initialization | 36.86619 | 36.26845 | −0.47039 | −33.15253 |
| Corrected boundary and perfusion, fresh initialization | 37.04819 | 37.10022 | −0.01272 | −0.32184 |

The last row has 96.91918 W metabolic plus skin input and 97.24102 W boundary
loss. One-hour near-balance is encouraging but does not establish all-system
equilibrium, clinical calibration, or six-hour validity.

The subsequent uninterrupted six-hour paired experiment completed. Both runs
loaded the same retained state (core 36.86619 °C), and their object inventories
differ only in the perfusion-corrected Energy compilation unit. With correct
perfusion the final core is **37.10135 °C**, versus **35.38694 °C** with the
sixfold flow error. Corrected total storage is +0.85121 W and final core slope
+0.01875 °C/hour; the last half-hour linear slope is +0.01491 °C/hour. The
control still loses 3.19542 W and cools at −0.04378 °C/hour. Corrected core
stays within 36.86639–37.10135 °C over the recorded run. These results establish
bounded normothermic native rest for this variant and boundary, not clinical
calibration or exact thermoregulatory equilibrium.

Corrected stored-heat change is −7593.969 J (skin cools while core warms),
with sampled integral error 1.179 J. Control stored-heat change is −376991.077 J,
with error 0.572 J. Final instantaneous circuit residuals are −0.168 µW and
+1.826 µW. The verifier resolves every retained execution input and confirms
21600 uninterrupted one-second observations and the actual 1× versus 6×
regional perfusion identities. Receipt: `thermal-rest-verification-2f6e_3wd`.

For the paired boundary experiment, one-second integration of net heat differs
from actual stored-heat change by 5.716 J for the source boundary and 0.0893 J
for the corrected boundary. These are sampling errors, separately reported
from the much smaller instantaneous circuit residual.

Variant chain and exact library SHA-256:

| Variant | SHA-256 |
|---|---|
| `whole_body_integrity_depletion` | `5b8840947b3bb964ba06f257a2874cc37efd3bd67f866565d34acbbcbb0b381d` |
| `whole_body_integrity_thermal_boundary` | `141dbc368a70040705874e782108b39541f4218a21002b2f81b243d683ebda0e` |
| `whole_body_integrity_thermal_boundary_v2` | `14ac1f6e28d7ee020a18659ba1808d23e0b86c2fad1084df2781308e6e5cf81f` |
| `whole_body_integrity_skin_perfusion` | `c26e7ede45bfb0058a4baaf931daef7d64ffe22b68802d71ce93104382e66b80` |

Each builder verifies inherited object identities, retains patched source and
diff, and relinks only the affected compilation unit. Existing variants are
immutable. The perfusion parent is boundary version 2.

Retained audit directories under `data/derived/audits/`:

- `long-horizon-thermal-xnzladb2`: original-state and failed six-hour audit.
- `thermal-boundary-native-ambu9748`: fresh paired one-hour boundary experiment.
- `thermal-boundary-v2-probes-x5la17f3`: actual nude and clothed finite-circuit checks.
- `thermal-boundary-runtime-verification-1dd58arx`: passed algebra/runtime ledger.
- `thermal-perfusion-fresh-hour`: fresh perfusion-corrected initialization and hour.
- `skin-perfusion-verification-6_wn2uj_`: actual corrected flow partition.
- `thermal-six-hour-perfusion-corrected` and `thermal-six-hour-perfusion-control`:
  completed uninterrupted paired six-hour rest from the **same initial state**.
- `thermal-rest-verification-2f6e_3wd`: passed paired six-hour source, clock,
  constitutive and thermal-energy checks.

The private long-rest runner changes only the existing batch runner's permitted
duration from 3600 to 21600 seconds. It is explicitly an audit batch executable,
not a persistent-session adapter or an accepted systemic publication. It retains
the private C++ source, diff, binary, actual dynamic dependency identities,
one-second outputs, initial/final states, and post-start runtime-resource archive.

## Independent JOS-3 comparison

The held [authors' JOS-3 implementation](https://github.com/TanabeLab/JOS-3),
revision `3c74ee2af2f79aa360093cc517e38bf465ec8c5b`, was executed with the same
generic age, size, sex and fat fraction, air/radiant temperatures, humidity,
speed and 0.5-clo uniform garment input. Physical activity ratio is 1.0. At one
hour, central blood temperature is 36.62860 °C standing or 36.62596 °C lying;
mean skin temperatures are 32.86091 and 33.12841 °C. Results are retained in
`thermal-boundary-jos3-azo23xyu`.

The reproducible `compare_native_thermal_jos3.py` subsequently ran both postures
for six hours with a 15-second step, retaining its own source and all held
solver inputs in `thermal-matched-jos3-bt7l13p0`. Final central blood temperature
is 36.60386 °C standing or 36.61475 °C lying; mean skin is 31.30697 or
31.54208 °C. The independently reconstructed implicit thermal step balance
closes within 2.0e−10 W. This is a six-hour model execution, not extrapolation
from the one-hour endpoint.

JOS-3 uses independent physiological setpoints, compartments and initialization;
its central blood temperature is not identical to the BioGears core observable.
Agreement in temperature scale is a cross-model check, not measurement-based
calibration or evidence that the two models can be substituted without mapping.

Held bedding cases also exist in `data/derived/thermal/index.json`, based on
[Akimoto et al., 2025](https://doi.org/10.1016/j.buildenv.2025.113074). They use
a different female profile and whole-body *total* insulation, which includes
the boundary layer. Their mattress/pillow, blanket and duvet values must not
be copied directly into the native garment-clo input.

## Remaining scope

The perfusion-corrected six-hour rest passed the bounded thermal audit. Subsequent
sweat-area and humidity corrections are documented in
[NATIVE_EVAPORATION_AUDIT.md](NATIVE_EVAPORATION_AUDIT.md); their one-hour native
rest checks pass, but their final variant needs long integrated acceptance.
Meal, hydration and exertion protocols need fresh initial states
under the final accepted variant and integrated multisystem checks; old six-hour
outputs cannot be relabeled as corrected. Native supine/bedding/contact heat
transfer is absent. Regional perfusion and alpha remain source assumptions.
The isolated evaporation corrections still need exertional and garment-moisture
validation before heat-strain interpretation. The corrections here do not
claim that one whole-body model has now been clinically calibrated.

Reproduction (builders require a new destination and preserve existing variants):

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_long_horizon_thermal.py
.venv/bin/python scripts/verify_native_thermal_boundary.py --runtime data/derived/audits/thermal-boundary-native-ambu9748 --nude data/derived/audits/thermal-boundary-v2-probes-x5la17f3/nude
.venv/bin/python scripts/verify_native_skin_perfusion.py data/derived/audits/thermal-perfusion-fresh-hour/states/native_final.xml
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_thermal_rest.py data/derived/audits/thermal-six-hour-perfusion-corrected --control data/derived/audits/thermal-six-hour-perfusion-control
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/compare_native_thermal_jos3.py --seconds 21600 --dt 15
```
