# Independent source thermoregulation

The workbench executes the authors' [JOS-3 repository](https://github.com/TanabeLab/JOS-3) at commit `3c74ee2af2f79aa360093cc517e38bf465ec8c5b` (package version 0.5.0, MIT). It imports the checked source directly using the existing NumPy environment; no globally installed JOS-3 package or modified upstream equations are used. Every source module, license, paper PDF, and the local adapter is hashed in `data/raw/thermal/source.json`.

The actual model contains **85 thermal nodes and 17 body regions**, including central blood, regional arteries/veins, core, skin, and the source-defined subset of muscle/fat/superficial-vein nodes. `lying` selects the published source supine heat-transfer coefficients. This is physical source posture behavior, distinct from the anatomical viewer's rigid display rotation. It does not add venous hydrostatics, sleep stages, or lateral sleeping mechanics.

## Published bedding boundary

[Akimoto et al., Building and Environment 279 (2025), 113074](https://doi.org/10.1016/j.buildenv.2025.113074) supplies the measured total-insulation inputs. The author-hosted PDF is retained at `data/raw/thermal/bedding.pdf`; Table 9 was checked visually. Table 6 specifies a hypothetical female, height 1.68 m, weight 61.5 kg, age 20, fat 15%, cardiac index 2.59 L/min/m², and PAR 1.0. These describe a separate source model, not the BioGears patient.

| Case | Air/radiant temperature | Total insulation IT | Source condition |
|---|---:|---:|---|
| `lying_default` | Actual constructor neutral environment, approximately 28.7654 °C | Source nude clothing model | All constructor body/environment defaults, then lying |
| `supine_mattress` | 22.6 °C | 1.15 clo | Nude, mattress/pillow, 23.3% contact coverage |
| `supine_blanket` | 22.6 °C | 2.78 clo | Pajamas/blanket, 94.1% coverage |
| `supine_duvet` | 22.6 °C | 5.46 clo | Pajamas/duvet, 94.1% coverage |

For bedding, the explicit adapter follows paper Eq. 7/8:

`RT = 0.155 IT` in m² K/W; `ReT = RT / (0.38 × 16.5)` in m² kPa/W.

**IT is total insulation, not ordinary clothing Icl.** The paper's prose calls RT the reciprocal of IT, but its Eq. 2/7 and dimensions imply the conversion above. The adapter replaces only the source's dry/evaporative resistance functions during the run and restores them afterward. This is an identified paper-boundary implementation. It is not an undocumented change to JOS-3 or BioGears.

Whole-body measured IT is applied uniformly to all regions, one of the paper's comparison methods. Detailed regional manikin measurements were not supplied, so no local mattress/contact map is invented. The paper used more than 2000 minutes to approach steady state; these runs intentionally record the first hour after changing the source's neutral standing initialization to lying and the specified boundary. Their endpoints do not reproduce the paper's equilibrium skin temperatures.

## Conservation and timestep evidence

At every step the adapter observes the actual source matrix without changing it. With heat capacities C (J/K), directed transfer W (W/K), K = diag(row-sum W) − W, boundary conductance B, operative temperature To, and net source heat Q, it checks:

`C (Tnew − Told)/dt + K Tnew = B (To − Tnew) + Q`.

Column sums of K test internal energy conservation, including the directed circulation. The summed equation checks storage against external and source heat. The source's displayed sensible heat loss uses the previous temperature; the audit uses the implicit endpoint temperature required by its solved equation. Confusing these evaluation times would manufacture a heat-budget error.

Runs compare 60, 30, 15, 7.5, 3.75, 1.875, and 0.9375-second steps at common physical times. Early maximum differences are not uniformly monotonic because the first sampled instant moves through the rapid initial vascular transient. Endpoint and post-60-second discrepancies are reported separately; endpoint refinement is verified below 0.001 °C at the finest pair. This is numerical agreement, not human measurement accuracy.

The index records measured residuals, refinement errors, and final temperatures. Per-case coefficient files expose capacities, conductances, source node indices, body surface areas, and the final evaluated circulation/boundary matrices. `comparison.png` and `.svg` show central blood and mean skin temperature; central blood is not labeled as every organ's core temperature.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/collect_thermal_model.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_thermal_model.py --artifacts
```

Acquisition requires the pinned clean checkout. Outputs live under `data/derived/thermal/`. The original JOS-3 [model paper](https://doi.org/10.1016/j.enbuild.2020.110575) is retained as a separate hashed PDF. Neither the measured manikin insulation nor a conservative source solver establishes calibrated human sleep or bed-rest predictions.

At the finest recorded step (0.9375 s), the actual one-hour results are:

| Case | Central blood °C | Mean skin °C | Largest endpoint difference from 1.875 s, across all 85 nodes |
|---|---:|---:|---:|
| Source default lying | 36.75919 | 34.73876 | 0.000028 °C |
| Nude mattress/pillow | 36.76711 | 33.00866 | 0.000561 °C |
| Pajamas/blanket | 36.94693 | 35.39174 | 0.000058 °C |
| Pajamas/duvet | 37.00398 | 35.82317 | 0.000059 °C |

Maximum whole-body discrete heat residual across these runs is below 3.3×10⁻⁹ W. The final-node refinement errors decrease with halving the final steps. Tests also verify positive capacities, nonnegative symmetric conductive couplings, restoration of the original boundary functions, deterministic repeated source initialization, and source-file hashes. These checks concern the implemented numerical model and its evidence record.

## Independent runtime review

The wrapper validates clocks before loading or initializing JOS-3: duration is
positive and at most 604800 s; timesteps are 0.01–3600 s; the duration must contain
1–100000 whole timesteps. Boolean, nonfinite, fractional-step, and excessive-work
requests are rejected explicitly. These are computational bounds, not physiological
safety limits. The heat ledger rejects invalid coefficients and arithmetic overflow
instead of returning infinite diagnostic values.

Runtime source verification checks the actual clean Git checkout and complete source
hash inventory, rather than trusting the revision string alone. Historical acquisition
and generated-artifact `adapter_sha256` fields remain intact; new executions additionally
record `runtime_adapter_sha256` for the code actually executing them. Boundary substitution
holds the shared reentrant lock and restores both source functions even after an exception.
The observation hook preserves a source exception rather than replacing it with a missing
local-variable error during exception unwinding. None of these changes modifies the
source thermal equations or the recorded initialization procedure.
