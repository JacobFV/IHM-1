# Native sweat, humidity and regional evaporation

Two isolated variants correct the native evaporative boundary while preserving
the historical source, dry-boundary and skin-perfusion variants. The corrections
follow dimensional and producer/consumer evidence; no thermal coefficient was
fitted to a desired body temperature.

## Whole-body sweat is a mass rate, not a local heat flux

At pinned BioGears revision `3f16a5fa1dade9c511b88d923606fa51cc35e95d`,
`Energy::CalculateSweatRate` multiplies `sweatRate_kg_Per_s` by timestep, removes
that mass from the patient, updates lost electrolytes, and drives the single
`SkinSweating` fluid path. `GetSweatRate()` therefore represents the whole body.

The Environment consumer multiplied that rate by latent heat (J/kg), obtaining
watts, but divided the result by an evaporative potential in W/m² as though
the result were dimensionless wettedness. It added this watt quantity to a
diffusion flux in W/m², then multiplied by each regional skin area. The sweat
term consequently gained an extra numerical factor of whole-body area.
Unbounded wettedness also allowed heat loss above local evaporative capacity.

`whole_body_integrity_sweat_evaporation` distributes total sweat latent power S
with an explicit uniform sweat-flux prior over area A. For region i:

```text
local available sweat flux = S/A                                  [W/m²]
local capacity p_i = max(0, vapor-pressure difference / resistance) [W/m²]
evaporated sweat flux e_i = min(S/A, p_i)                          [W/m²]
sweat wetted fraction w_i = e_i/p_i (zero if p_i is zero)           [1]
regional heat = A f_i [e_i + (1 − w_i) d_i p_i]                    [W]
```

The source diffusion fraction dᵢ is retained (0.06 for ordinary skin). Thus
total sweat evaporation cannot exceed S, no region exceeds its capacity, and
subdividing an area cannot create latent energy. Unevaporated sweat has already
left the fluid circuit; it supplies no additional latent cooling. This does
not yet model liquid retained in clothing, runoff, condensation, sensible heat
carried by sweat, or later re-evaporation. The helper exposes the unused latent
budget to diagnostics without removing water a second time.

The [authors' JOS-3 evaporation implementation](https://github.com/TanabeLab/JOS-3/blob/3c74ee2af2f79aa360093cc517e38bf465ec8c5b/src/jos3/thermoregulation.py)
independently uses local area, vapor-pressure differences and a capped wetted
fraction. Its physiological sweat distribution is different; no JOS-3 control
coefficients were copied into the native BioGears sweat producer.

## Relative humidity and local skin temperature

`Environment::CalculateSupplementalValues` assigns saturation pressure at air
temperature to `m_dWaterVaporPressureInAmbientAir_Pa`. It multiplies that field
by relative humidity only in a separate local variable used for air density.
The field itself remains saturation pressure. The evaporative consumer used
that unscaled field. A direct native probe confirms that its evaporative
boundary is unchanged between 0, 50 and 100% RH at fixed temperatures.

`whole_body_integrity_evaporation_humidity` applies RH once at that consumer:

```text
ambient actual vapor pressure = RH × saturation pressure(T_air)
regional skin saturation = saturation pressure(T_skin_i)
```

Each `External*SkinToGround` path starts at its actual regional skin thermal
node, so the consumer now uses that node's current temperature. The old
consumer reused torso skin vapor pressure for all six regions. The source
Antoine saturation function, transfer coefficient, diffusion fraction and
latent heat polynomial are preserved. The original ambient saturation field
is unchanged, avoiding double application of RH in air density or respiration.

The corrected gradient uses the same dimensional structure as the held JOS-3
source. This is an equation/implementation check, not a human humidity-response
calibration. Native regional garment evaporative permeability still uses its
inherited coarse coverage assumptions; clothing wetness is not coupled yet.

## Actual native checks

`native_evaporation_probe.cpp` loads a retained state, explicitly unlocks only
the thermal circuit for a direct constitutive probe, supplies controlled
temperatures, area and sweat rate, then calls actual native `PreProcess`.
It does **not** integrate or save this artificial test state as physiology.
High sweat inputs exercise capacity limits and may exceed the native
physiological controller's permitted production rate.

Each variant executes 432 regional cases: four total areas (0.7, 1, 1.9012784,
3 m²), three humidities, uniform 34 °C or unequal 24–34 °C skin, three sweat
rates (0, 1e−6, 1e−3 kg/s), and six regions. The independent Python verifier
reconstructs local budgets from native coefficients and pressures. Small-sweat
increments are independent of total area, after accounting for displaced
diffusive cooling.

At area 1.9012784 m², uniform skin 34 °C, air 22 °C and zero clothing:

| Native variant | RH | Zero-sweat heat W | 1e−6 kg/s sweat heat W | 1e−3 kg/s sweat heat W |
|---|---:|---:|---:|---:|
| Before sweat correction | any tested | 19.71804 | 24.10022 | 4401.90171 |
| Sweat-area/capacity correction | any tested | 19.71804 | 22.02290 | 328.63397 |
| Humidity/local-temperature correction | 0% | 39.16644 | 41.47130 | 652.77395 |
| Humidity/local-temperature correction | 50% | 29.44224 | 31.74710 | 490.70396 |
| Humidity/local-temperature correction | 100% | 19.71804 | 22.02290 | 328.63397 |

Source-law tests were observed failing on the predecessor before each patch.
The area/capacity predecessor misses the expected regional budget by up to
2466.95 W in the artificial high-sweat test. The area-corrected variant passes
to 2.84e−14 W. It still fails the RH/local-temperature test by 186.07 W;
the humidity-corrected variant passes to 5.68e−14 W.

The actual helper also passes zero/negative evaporative capacity, zero sweat,
unequal subdivision, conservation of sweat latent power, and finite wettedness
tests. Numerical conservation is distinct from physiological calibration.

## Immutable lineage and receipts

| Variant | Native library SHA-256 |
|---|---|
| Skin perfusion parent | `c26e7ede45bfb0058a4baaf931daef7d64ffe22b68802d71ce93104382e66b80` |
| `whole_body_integrity_sweat_evaporation` | `c1b6616c924532fd542458e29b7ffe07e9cc4f5d3061daf2f78f3a7832b88ae5` |
| `whole_body_integrity_evaporation_humidity` | `ab485812968b17e16915b04a6e5144bb1e1ec7f5bf78213936c3334d087c0a5a` |

Each builder retains parent library and object identities, source diff, copied
headers, compiler/linker commands, and result identity. Only Environment.cpp
is recompiled; the perfusion-corrected Energy object remains unchanged.

Receipts under `data/derived/audits/`:

- `sweat-evaporation-algebra-rto0xcd8`: actual native helper checks.
- `native-evaporation-skin_perfusion-vl32lbsq`: source-law negative control.
- `native-evaporation-sweat_evaporation-lhxkf6c6`: area/capacity native probe.
- `native-evaporation-verification-24mlxpqm`: passed area/capacity verification.
- `native-evaporation-evaporation_humidity-ogvv_8zk`: humidity/local-temperature probe.
- `native-evaporation-verification-e941emd9`: passed humidity verification.
- `thermal-sweat-evaporation-fresh-hour` and
  `thermal-evaporation-humidity-fresh-hour`: fresh initialization and integrated
  one-hour rest checks, both completed and passed.
- `thermal-rest-verification-nl86eakv` and `thermal-rest-verification-u6sclow6`:
  inherited dry-boundary, actual 1× perfusion partition, finite 3600-second
  clock, retained execution inputs and full thermal energy checks.

Fresh sweat-area rest starts at 37.04819 °C and ends at 37.10034 °C, with
final core slope −0.01254 °C/hour and total storage −0.33533 W. The humidity
variant starts at 37.04812 °C and ends at 37.05208 °C, with final core slope
−0.06705 °C/hour and total storage −4.15987 W. One-second integrated heat versus
stored-heat errors are 0.566 J and 1.091 J. The humidity correction modestly
increases resting heat loss, as expected from the restored vapor gradient;
its one-hour endpoint is not a settled thermal equilibrium. Longer integrated
rest and intervention acceptance remain necessary for this final variant.

For matched downstream protocols, the fresh initial state and execution lineage
were frozen separately in
`data/derived/audits/thermal-final-initial-state-ctkr8bap/initial-state.json`.
Its `state_path` resolves to a retained copy with SHA-256
`cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3`.
This is the state before the one-hour rest, not the final endpoint. Its initial
total heat-storage rate is −26.38209 W, so native source stabilization must not
be called thermal equilibrium. Every intervention/control should load the same
frozen state under the same final library and retain its own runtime inputs.

Probe inputs include retained source, binary, initial state, exact linked
dependencies and post-start runtime-resource bytes. Failed initial probe
attempts were retained; they failed on intentional native read-only protection
before the explicit diagnostic thermal unlock was added. No production adapter
was weakened.

```sh
.venv/bin/python scripts/verify_native_sweat_evaporation.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_evaporation_probe.py data/derived/audits/native-evaporation-sweat_evaporation-lhxkf6c6
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_native_evaporation_probe.py data/derived/audits/native-evaporation-evaporation_humidity-ogvv_8zk --humidity
```

These variants are experimental until their fresh-state integrated rest and
intervention protocols pass the whole-body acceptance checks. They cannot
retroactively repair previously published source-run histories.

## Inherited source regressions

`verify_final_thermal_inheritance.py` executes the held calcium digestion,
bilateral renal transport, dry aqueous GI, energy and depletion regressions
against the final humidity library. Original verifier files are untouched.
The calcium, renal and dry-GI probes use their existing parameterized API. The
energy and depletion drivers are copied with retained diffs to select the final
library and classify the inherited fixes correctly. The energy driver omits
only historical *cross-variant* resting-output parity: corrected thermal laws
are expected to change resting outputs. Its cold/neutral partition, thermal
source, exercise bound, stop behavior, store positivity and within-variant
matched causal checks remain active.

GI, renal, dry-GI and energy use the final fresh state. Depletion deliberately
retains its original known water-store/clock fixture for the 1205-second empty
stomach regression. `archive_thermal_inheritance.py` checks that fixture and
the linked dependencies against execution receipts and retains their exact
bytes. Its runtime-resource capture is explicitly after execution, not an
assertion of historical startup identity.

All inherited regressions passed in
`data/derived/audits/thermal-final-inheritance-tudiahnb/verification.json`:
6 calcium cases, 20 bilateral renal cases, 11 dry-GI cases, 15 energy checks,
and 5 depletion checks. The actual uninterrupted depletion replay reaches
1205 seconds with exactly zero stomach water, including a safe subsequent
step after exhaustion; its mixed litre/millilitre remainders and invalid
negative-input rejection also pass. At 120 seconds, continuous exercise uses
150.12 W additional demand and raises oxygen consumption from resting
236.13 to 479.12 mL/min. Stopping at 30 seconds clears that demand and gives
236.13 mL/min at 120 seconds. These are bounded native causal regressions,
not sustained exercise or six-hour multisystem acceptance.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_final_thermal_inheritance.py
# Archive the fresh directory printed by the verifier:
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/archive_thermal_inheritance.py data/derived/audits/thermal-final-inheritance-tudiahnb
```
