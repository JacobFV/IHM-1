# Native skin–lymph circuit predictor

`ihm.coupling.circuit.FluidCircuit` compiles a selected coherent native fluid circuit into a hydraulic descriptor system with native coefficients. This is a small-signal, frozen-operating-point predictor. It does not replace BioGears' changing elastance, oncotic pressures, osmotic transport, feedback or valve logic.

Run `.venv/bin/python scripts/verify_circuit_predictor.py`. To generate the actual native skin artifact:

```python
from pathlib import Path
from ihm.coupling.circuit import build_skin_circuit
result = build_skin_circuit(Path.cwd())
```

The builder reads `data/derived/native-circuits/graph.json`, verifies its underlying native XML hash, and writes `data/derived/coupling/native-skin-circuit.json`. It preserves the graph hash and raw native node/path properties, including source units, SI values and original baseline/current distinctions.

## Native selection and boundary conditions

The selection comes exclusively from the coherent `FullCardiovascular` view: Skin1, Skin2, SkinE1, SkinE2, SkinE3, SkinI, SkinL1, SkinL2, Lymph, Aorta1, VenaCava and Ground. Aorta1, VenaCava and Ground are fixed pressure boundaries. Their initial pressures come from the saved state. Cross-cut native paths into other selected free nodes, chiefly other tissues draining into Lymph, become fixed snapshot inflows. Paths crossing a prescribed-pressure boundary are recorded and absorbed into its necessary boundary inflow; no extra organ network is duplicated.

All resistances, compliances, pressure sources and flow sources are current native SI coefficients. No illustrative physiological coefficient is supplied. Initial pressures and available volumes are native state values. The 100 Pa perturbation used in the exported experiment is an explicitly hypothetical input, not a measured intervention or fitted physiology.

## Equations and numerical form

Use incidence +1 at a path source and −1 at its target. A resistance carries `q=(p_source-p_target+S)/R`. A standalone native pressure source constrains `p_target-p_source=S`; a plain conducting wire constrains zero pressure drop. A prescribed flow source retains its signed native flow. A compliance carries `q=C*d(p_source-p_target)/dt` and contributes `C*a*a.T` to the nodal storage matrix. Compound element forms not supported by this compiler and inertances requiring extra branch states fail explicitly.

Native **Closed** valves and switches conduct; **Open** valves and switches block. These states are frozen. An ideal conducting valve is a modified-nodal-analysis pressure constraint with its flow as a Lagrange multiplier. Reverse flow through a frozen conducting valve is reported in `gate_violations`; the code does not pretend it remains a valid diode regime or silently change topology.

The compiled equations are `M p_dot + K p + A j = b`, together with ideal pressure constraints. Prescribed boundary pressures are eliminated, including their compliance derivative terms. Numerical pressure unknowns are scaled by 10,000 Pa and ideal-flow unknowns by 1e-6 m³/s. Flow balance equations and pressure constraints are scaled separately, preventing raw SI magnitudes from creating avoidable conditioning problems. These numerical scales are not physiological parameters. Singular selections or redundant fixed-to-fixed ideal constraints fail compilation.

`step(dt, boundary_pressures_pa=None)` uses backward Euler and mutates the state. Boundary changes are interpreted over that numerical timestep; an ideal abrupt pressure step across compliance can create a finite-step approximation to an impulse. Smaller timesteps should be used when resolving fast transients. The result contains physical pressures, signed path flows, native volume changes, current volumes, free-node KCL residuals, required boundary inflows, valve validity flags and the scaled solve condition number.

Compliance flux is storage, not physical transport to Ground. Volume changes are its signed nodal accumulation integrated over the step, reported only where the native node actually has a volume state. A negative predicted volume is flagged as invalid extrapolation rather than clipped. Arbitrary large perturbations can invalidate frozen coefficients or gates; the caller must inspect these diagnostics.

## Transfer functions and poles

`response(s, boundary_name)` solves the complex small-signal descriptor pencil for a unit pressure perturbation at one declared boundary, with zero initial perturbation and other boundary/source perturbations held at zero. It returns per-node `pressures_pa_per_pa` and per-path `flows_m3_s_per_pa` arrays. `s` has inverse-second units; a Fourier frequency f corresponds to `2πif`. This is a causal LTI resolvent of the explicitly specified approximation, distinct from a finite-window Fourier or Laplace transform of a trace. Evaluation at a singular pole raises rather than inventing a finite response.

`finite_poles_per_s()` uses generalized descriptor eigenvalues. Infinite eigenvalues represent algebraic ideal constraints and are excluded from its finite-mode output. The native skin selection has finite poles approximately:

| Pole (s⁻¹) | Decay time |
|---:|---:|
| −0.07655154 | 13.06 s |
| −0.0001553642 | 1.79 h |
| −0.00005630455 | 4.93 h |
| 0 | Integrator |

The zero mode is the isolated intracellular compliance under a frozen native osmotic flow source. Its volume can drift secularly; this approximation does **not** establish long-horizon biological stability. The hour-scale modes arise from the native coefficients and boundary choices, not from fitting a 60-second trace or validating those timescales against human observations. The artifact evaluates response at σ=1e-8 s⁻¹ plus imaginary frequencies to avoid the zero pole; σ is recorded explicitly.

## API and artifacts

```python
model = FluidCircuit.from_native(graph, node_names, fixed_pressure_names,
                                 circuit_name='FullCardiovascular')
step = model.step(0.02, {'Aorta1': 15225.3})
response = model.response([1e-4 + 0j, 1j], 'Aorta1')
restored = FluidCircuit.from_dict(model.to_dict())
```

`to_dict` preserves state, volumes, ideal flows, selected raw native objects, cross-cut source fluxes, descriptor matrices/scales and source provenance. `from_dict` rebuilds matrices from serialized native coefficients and restores the state, rather than trusting arbitrary supplied matrix entries. Unknown units, wrong dimensions, nonpositive R/C, nonfinite source coefficients and unknown boundary names are rejected.

The generated JSON contains the executable initial `model`, one native-baseline frozen step, per-path comparison with the saved native snapshot, complex Laplace responses encoded as real/imag arrays, and a 300-second +100 Pa aortic boundary experiment sampled every five seconds. The latter is a prediction from the frozen circuit, not a native simulated patient scenario. Conditions and caveats are embedded in the file.

At the 0.02-second baseline step, maximum free-node KCL residual is approximately 6.0e-18 m³/s. The skin-to-lymph valve flow is 1.00885438e-8 m³/s versus the saved native 1.00885410e-8 m³/s; the skin vascular-to-interstitial flow differs by about 0.037%. This is an operating-point consistency check, not comparison with a withheld native future trajectory. The 300-second experiment has no frozen conducting-valve reversal or negative volume in its saved samples.

Tests cover an analytic RC time constant/pole and Laplace response, the exact backward-Euler step sequence, compliance storage/KCL, source pressure polarity, conducting/blocked gate semantics, reverse-valve diagnostics, physical dimensions, serialization continuation and real native one-step consistency. No new physiological coefficient is introduced to make those tests pass.
