# Source-derived cerebral circulation and CSF dynamics

`ihm.native.csf` executes the two-state Ursino–Lodi model from *A simple mathematical model of the interaction between intracranial pressure and cerebral hemodynamics*, J Appl Physiol 82:1256–1269 (1997), [DOI](https://doi.org/10.1152/jappl.1997.82.4.1256). The [original author-upload full text](https://www.researchgate.net/publication/14112614_A_simple_mathematical_model_of_the_interaction_between_intracranial_pressure_and_cerebral_hemodynamics) provides Table 1 and Appendices A–C. This is an independently implemented source simulation, not official executable author code or a patient calibration.

Publisher requests returned access challenges. The locally retained indexed author-upload retrieval is hash-pinned (`69a989fb59cbca3e2317a6553ded045dbed340ff4558177ea6ae523d23b386ea`) in `data/raw/csf/ursino_lodi_1997/indexed-author-source-review.json`. This hash authenticates the reviewed retrieval record, not original PDF bytes. The acquisition manifest retains failed HTTP responses separately. Copyright is American Physiological Society, 1997; repository availability does not establish an open reuse license. Source mathematical facts are implemented without redistributing the full article in the repository.

## Equations and coefficient interpretation

Native source units are mmHg, mL and seconds. States are intracranial pressure `Pic` and arterial compliance `Ca`. The prescribed arterial pressure port supplies `Pa` and its derivative. With basal sinus pressure `Pvs`, arterial volume is `Va=Ca(Pa-Pic)`, arterial resistance is `Ra=kR Can²/Va²`, and craniospinal compliance is `Cic=1/(kE Pic)`. Exact Appendix A3 determines capillary pressure:

```
(Pa-Pc)/Ra = (Pc-Pic)/Rf + (Pc-Pic)/Rpv
q = (Pa-Pc)/Ra
qf = max(0,(Pc-Pic)/Rf)
qo = max(0,(Pic-Pvs)/Ro)
x = (q-qn)/qn
sigma = Can - delta/2 tanh(2 G x/delta)
dCa/dt = (sigma-Ca)/tau
dPic/dt = [Ca dPa/dt + (Pa-Pic) dCa/dt + qf-qo+Ii]/(Cic+Ca)
```

Here `delta=delta_Ca1` for negative `x` and `delta_Ca2` otherwise. The tanh form is algebraically identical to A9 with `ks=delta/4`. `Ii` is an explicit injection rate available in the equation evaluator; shipped trajectories use zero. The implementation rejects states outside the source's collapsed terminal-vein regime `Pa > Pic > Pvs`, or nonpositive compliance. It uses exact capillary conservation rather than the optional B1 approximation that neglects the CSF branch.

| Parameter | Source basal value | Unit |
|---|---:|---|
| Ro | 526.3 | mmHg s/mL |
| Rpv | 1.24 | mmHg s/mL |
| Rf | 2380 | mmHg s/mL |
| Can | 0.15 | mL/mmHg |
| delta_Ca1 | 0.75 | mL/mmHg |
| delta_Ca2 | 0.075 | mL/mmHg |
| kE | 0.11 | 1/mL |
| kR | 49100 | mmHg³ s/mL |
| tau | 20 | s |
| qn | 12.5 | mL/s |
| G | 1.5 | mL/mmHg per fractional CBF change |
| Pvs | 6 | mmHg |

Initial `Pic=9.5`, `Ca=0.15`, and basal `Pa=100` follow rounded Table 1 values. Equation A9 and the tabulated deltas imply compliance saturation at 0.525 and 0.1125 mL/mmHg. The paper's prose instead describes sixfold/half basal values (0.9/0.075). This implementation preserves the literal equation/table, documents the discrepancy, and does not silently adjust coefficients. Rounded initialization is not an exact equilibrium: the computed basal equilibrium ICP is approximately 9.4264 mmHg, with local poles about −0.00215345 and −0.701211 per second. These are local mathematical modes, not measured human spectra.

## Execution and one-way coupling

```
.venv/bin/python scripts/collect_csf.py --offline
.venv/bin/python scripts/verify_csf.py
```

Offline building requires the hash-verified review artifact. Ordinary collection retries bounded HTTP acquisition and preserves responses; if the reviewed original is absent or altered, it writes an unavailable status instead of running undocumented coefficients. A fresh checkout must obtain the reviewed source artifact separately; ignored research assets are not embedded in source control.

`PressureInput.constant()` uses published 100 mmHg. `PressureInput.from_samples(time_s, values, unit='Pa'|'mmHg', provenance=...)` performs unit-checked piecewise-linear interpolation and analytic segment derivatives without extrapolation. `PressureInput.from_native_map(path)` requires explicit `#Time(s)` and `MeanArterialPressure(mmHg)` CSV headers, records the complete CSV SHA-256, and subtracts only the initial source timestamp. The solver resolves input intervals and evaluates `dPa/dt` in mmHg/s.

The native baseline's 3000-row, 50 Hz, 60-second MAP trace drives a separate CSF model for its available 59.98-second elapsed interval. MAP is appropriate to this low-frequency mean-pressure model; it is not an instantaneous arterial waveform. Native and literature subjects are unmatched. Published initial ICP/compliance remain in place when replacing the arterial boundary, so the resulting transient is not a jointly initialized patient. No ICP feedback enters BioGears.

`run_csf(root, pressure_input=..., duration_s=..., samples=...)` returns serializable channels and parameters. `build_csf_status(root)` exports `data/derived/csf/index.json` with `models[]`, plus `baseline.json` (600 seconds), `native_map_driven.json`, and `hypotension.json` (100→90 mmHg ramp at 10–20 seconds, basal parameters). The ramp boundary follows Figure 7; the run does not claim to reproduce the figure's altered parameter cases.

Trajectories include `time_s`, `state_values[2][n]`, `channels[{id,unit,si_unit,si_scale,values}]`, source/hash evidence, input provenance, solver tolerances, equation residuals, and limitations. Channels cover ICP, compliance, arterial/capillary pressures, arterial volume, CBF, CSF production and reabsorption. Units convert explicitly using 133.322387415 Pa/mmHg and 10⁻⁶ m³/mL. The reported linearization always identifies its constant arterial operating pressure.

Verification checks source-scale flow targets, exact capillary and intracranial storage conservation, sigmoid slope/saturation, SI pressure equivalence, interpolation bounds, independent DOP853/Radau integrations and local stability. Basal 600-second solvers agree within 8.1×10⁻¹⁰ mmHg; capillary balance residuals are below 9×10⁻¹⁵ mL/s. Numerical agreement does not establish clinical accuracy. The model excludes venous compliance, rapid ICP pulsatility, posture-specific outflow and glymphatic transport; a minute of native MAP cannot validate its slow dynamics.
