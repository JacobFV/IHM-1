# Temporal evidence and predictors

Run `.venv/bin/python scripts/verify_temporal.py`, optionally acquire public measured evidence with `.venv/bin/python scripts/collect_temporal_evidence.py`, then build with `.venv/bin/python scripts/build_temporal_atlas.py`. Outputs are `data/derived/temporal/index.json` and one JSON per source run. Raw waveform bytes, source license, upstream SHA256SUMS, recording header and acquisition provenance are preserved in `data/raw/temporal/bidmc/`. No additional Python dependency is required beyond NumPy/SciPy.

## IBM interpretation

Read-only inspection covered the sibling IBM-1 files `ibm/fields/uncertainty/spectral.py`, `ibm/fields/priors.py`, `ibm/forge/spectra.py`, `ibm/processes/base.py` and `ibm/processes/mechanical.py`; the generated index records their SHA256 hashes. IBM's `TemporalBasis` is an orthonormal real Fourier transform on a finite cyclic time graph, with graph eigenvalues `4 sin²(πk/N)`. This is a temporal **Laplacian** basis, distinct from a causal **Laplace** transform. `SpectralGaussian` represents complex coefficient uncertainty; its variance is not automatically a signal-unit²/Hz density. Its low-pass, resonator and mechanical dynamics separately evaluate explicit LTI transfer functions at imaginary frequency. Finite periodic differentiation/delay assumes boundary continuity; it is not a general causal initial-value solver.

This implementation keeps those concepts separate:

* `spectral_estimate(x, fs, nperseg=...)` returns a one-sided physical PSD via Hann Welch normalization, in x-unit²/Hz. Segments do not overlap, segment means are removed, and incomplete trailing segments are reported and excluded. Density integrates by bin width to tapered variance. Graph eigenvalues and coefficient variances are not mislabeled physical power density.
* `cross_spectrum` uses `conj(X)*Y`, so cross phase is Y relative to X. Coherence is undefined (`null`) where either power is zero or only one segment exists. No independent degrees of freedom or confidence interval is claimed for correlated segments. Coherence is association, not causation.
* `finite_laplace(time_s,x,s)` integrates `x(t) exp(-s*t)` from the first to last sample by trapezoidal quadrature, using elapsed seconds. The function preserves complex amplitude and physical x-unit·s. There is no tail extrapolation. A finite-window transform is entire and its peaks are not system poles. The atlas separately stores real decay σ in s⁻¹ and frequency f in Hz, with `s=σ+2πif`; whole-record mean subtraction is explicit. Users can call the function on undetrended data when DC is needed.
* `resolvent(A,B,C,s,D)` returns `C(sI-A)⁻¹B+D` for a specified continuous linear state model. A must use inverse seconds. This is a causal system transfer only under the stated LTI model and zero initial conditions, not a transfer inferred from coherence. Singular evaluations at poles fail explicitly.

## Stable empirical dynamics

`fit_predictor` standardizes using the training prefix only, finds its SVD subspace, fits a regularized autonomous reduced VAR(1), and clips transition singular values below one. This contraction bound prevents both unstable eigenvalues and unbounded transient amplification in standardized reduced coordinates. It can bias real growth and poorly capture oscillations; its holdout error is retained, including failures against the persistence baseline.

The default split is 75% chronological training and 25% uninterrupted holdout. No holdout samples affect centering, scaling, rank or coefficients. The forecast starts from the last training sample and is never reinitialized with future truth. Reports include raw singular values, data/retained ranks, reduced design condition (null if deficient), constant training channels, transition eigenvalues, contraction, per-variable training one-step errors, physical-unit holdout RMSE, standardized RMSE and persistence RMSE. Hyperparameters are fixed a priori rather than selected on this holdout. JSON serialization preserves the model, cadence and diagnostics; `ReducedPredictor.from_dict(...).forecast(initial,steps)` executes it. Discrete eigenvalues are not reported as uniquely identified continuous physiological poles.

An autonomous trajectory fit does not learn intervention effects, microscopic parameters, causal coupling or a full human model. A poor holdout forecast is evidence against that empirical approximation. Native scenario simulation remains the mechanistic predictor.

## Sources and sampling

The native run contains 3,000 samples at 50 Hz and 22 physiological variables; the instrumentation-only CTSresistance column is excluded. Sample times span 0.02–60 s (59.98 s between endpoints). Units come from native CSV request headers. Constant channels have no reported dominant frequency. Ten-second segments yield 0.1 Hz bins and six segments. A 60-second trace cannot establish slow endocrine, renal, immune or thermal time constants.

The optional measured run is record 01 of the [BIDMC PPG and Respiration Dataset 1.0.0](https://physionet.org/content/bidmc/1.0.0/), DOI [10.13026/C2208R](https://doi.org/10.13026/C2208R), credited to Pimentel, Johnson, Charlton and Clifton. Cite Pimentel et al., *Towards a Robust Estimation of Respiratory Rate from Pulse Oximeters*, DOI [10.1109/TBME.2016.2613124](https://doi.org/10.1109/TBME.2016.2613124). License: Open Data Commons Attribution License v1.0. This open subset is from hospital recordings of critically ill patients.

Its verified WFDB header declares 60,001 samples at 125 Hz, spanning 480 s, with RESP (`pm`, impedance respiratory surrogate), PLETH (`NU`) and three ECG leads (`mV`). Units are read from that header rather than guessed from channel names. The CSV prints timestamps to only 0.01 s above 100 s, producing duplicate timestamps despite 0.008-second sampling. The importer reconstructs the sample clock from the header, checks that the discrepancy is at most 0.005 s and verifies sample count; measured discrepancy is 0.004 s. It does not resample or alter waveform values. This reconstruction is recorded in `time_basis`.

One hospital record is not independent validation of a different simulated patient. Impedance respiration is not measured lung volume, and ECG voltage is not arterial pressure. The atlas deliberately retains separate source runs and units; it performs no amplitude calibration or cross-subject coherence calculation. Raw acquisition is bounded and checksum-verified. Missing measured downloads leave the native atlas useful.

## JSON contract

`index.json` has `schema_version`, `runs`, `limitations` and `ibm_reference`. A run contains `id`, `source_kind` (`source_simulation` or `human_measurement`), source path/hash, time basis, cadence and horizon, plus:

* `variables`: id, label, unit, moments, constant flag, frequency/PSD arrays and estimator settings.
* `trace`: display-only decimated `time_s`, `variable_ids` and variable-major `values`.
* `cross_spectra`: x/y ids, common frequency axis, real/imag arrays, coherence, units and convention.
* `laplace`: σ array, frequency array and per-variable real/imag matrices shaped σ-by-frequency, units, detrending and quadrature metadata.
* `predictor`: executable serialized `model`, ordered variable ids/units, display-decimated holdout clock, variable-major actual/predicted traces, training boundary and diagnostics.

All fitting and spectral calculations use full-rate samples; stride decimation affects plotted traces only. `null` denotes undefined values, never zero evidence. JSON forbids NaN/Infinity. The analytic suite checks exponential Laplace integrals, absolute-clock invariance, sinusoidal frequency and integrated power, coherence/zero-power behavior, continuous transfer response, known stable dynamics, serialization, holdout leakage, deficient rank, input rejection and rounded-clock reconstruction.

## Materialized verification

The initial atlas is approximately 1.44 MB. Native arterial pressure has its strongest non-DC PSD bin at 1.2 Hz and lung volume at 0.3 Hz. The training matrix has rank 20, reduced to eight; spectral radius is 0.99910. Native holdout RMSE improves on constant persistence in only 7 of 22 channels. This is a substantive limitation, preserved in the report rather than concealed by one-step fit error. The measured recording has rank five and spectral radius 0.99917; its mean-reverting forecast beats last-sample persistence on five channels, which does not establish accurate beat-by-beat forecasting. Its RESP peak is 0.4 Hz. ECG's strongest bin (4.6 Hz in lead II) may be a waveform harmonic and must not be converted directly into heart rate. Spectral shape alone does not identify clinical rhythms.

## Hour-scale evidence and causal hydraulic response

The atlas also discovers completed `native_hour_*` experiments whose source configuration specifies 3600 seconds and 1 Hz output. Each retains the CSV hash and the complete condition configuration/hash. The native rest run represents the source environment at 22°C and 0.5 clo: core temperature falls from about 36.56°C to 33.76°C. It is not a stationary healthy baseline. Initial redistribution of approximately 500 mL gastrointestinal water contributes to fluid shifts; inspect the native fluid-budget report. Controlled thermal runs retain their individual ambient temperature and clothing values rather than sharing an assumed condition.

Hour-run Welch windows span 900 seconds (frequency-bin width 1/900 Hz). Finite-transform sigma values are 0, 1/3600, 1/600 and 1/60 per second, with frequency samples through 0.05 Hz. These resolve more of the slow source drift while retaining the distinction between finite-horizon transforms and an infinite-time system transfer function. Mean removal does not make drift stationary, and a spectral peak in cooling does not establish a physiological cycle.

At 1 Hz, Nyquist is 0.5 Hz. Instantaneous `ArterialPressure` and `TotalLungVolume` are excluded from the hour run's variables, spectra and empirical predictor because pulse/breath oscillations can alias. Their original columns remain in the hashed source CSV and their names appear in `omitted_undersampled_channels`. No source anti-alias filtering is assumed for remaining outputs. The original 60-second 50 Hz native run and 125 Hz BIDMC measured record remain unchanged.

The top-level `causal_responses[]` is a separate atlas collection, not a measured spectral run. It evaluates the actual frozen native skin/interstitial/lymph descriptor against an Aorta1 pressure perturbation at frequencies from 10⁻⁷ to 0.1 Hz and sigma 0 or 1/3600 per second. Complex pressure gains have units Pa/Pa; flow gains have units m³/s/Pa. Coefficients, source snapshot and circuit-artifact hashes remain attached. The intracellular storage integrator prevents evaluation at s=0; positive frequencies are explicit. Frozen gates, fixed external flows and local small-signal validity limit this causal response. Observed coherence cannot establish these causal gains.

CSF and reproductive source trajectories are linked as additional evidence, without inferring periodic physiology from their short nonstationary runs. CSF uses prescribed pressure responses and source-equation local poles; the reproductive example uses nonperiodic prescribed hormonal inputs and does not identify an autonomous menstrual cycle.
