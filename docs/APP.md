# Local integrated-human workbench

The app uses Three.js/WebGL to render the actual preprocessed source geometries served by the local Python API. It never substitutes generated anatomy. Each source family has its own coordinate frame; switching families removes the previous geometry. Upright/supine controls are rigid display rotations and do not establish native simulation posture or cross-specimen registration.

Start the Python API from the repository root on `127.0.0.1:8765`:

```sh
.venv/bin/python -c 'from ihm.app import serve; serve()'
```

Then, in a second terminal:

```sh
cd app
npm ci
npm run dev
```

Open the printed localhost URL. Vite proxies `/api` to port 8765. `npm run build` produces `app/dist` for the Python server to serve; no publication is required. All JavaScript and styles are local (no CDN runtime dependencies).

The anatomy library supports model selection, independent system layers, structure search, and picking. Drag to orbit, scroll to zoom, right-drag to pan. The inspector shows source, frame, specimen, units, calibration status and hashes where supplied. Opacity and section clipping apply to visible structures. Geometry loading uses six concurrent requests with cancellation when filters or families change. Source geometry remains selectable if optional physiological datasets are absent.

Flow playback appears only when decoded geometry has actual `frames` and `times`. Vector lengths use the supplied `vector_scale`; the UI does not invent or animate synthetic flow. Physiology and spectral panels show recorded API trajectories and the temporal atlas. Spectral sources are selectable (including BIDMC human recordings); choose Fourier power density or finite Laplace magnitude and its damping coefficient. A fixed archive-wide pressure color scale is retained across playback frames. Vascular walls start translucent so internal vectors remain visible. The human-measurement foldout reports held-out skin-field error and identifiability, with limits on what the fit establishes. Scenario submission uses bounded duration (1–600 s), StandardMale, 10 Hz sampling, and either baseline or low-intensity exercise with recovery; server validation is authoritative. Run status is polled every five seconds.

API contracts: `/api/manifest` exposes models and structures; `/api/geometry/{id}` exposes flat positions/indices (or line segments/vectors); `/api/physiology` exposes `time_s` and named `values`; `/api/temporal` exposes runs with per-variable frequency and PSD arrays; `/api/evidence` exposes status; `/api/scenarios` lists runs and accepts the native scenario configuration. Optional data has explicit empty/error states. Geometry failures do not fabricate replacements.

Verification:

```sh
cd app
npm test
npm run build
npm run test:browser
```

The browser test requires live API assets and the development server. It uses installed Google Chrome at `/usr/bin/google-chrome`; set `CHROME_PATH` or `APP_URL` to override. It exercises source loading, structure selection, empty search, opacity, clipping, posture and chart tabs, captures `app/test-results/workbench.png`, and rejects browser runtime errors. By default no native scenario is launched. Set `RUN_NATIVE_BROWSER=1` to additionally execute a real two-second exercise/recovery scenario through the form; this checks submission and completion against the live backend.

Additional evidence views use the local coverage, conservative coupling, CFD audit, and reproductive APIs. The source reproductive example retains its 0–10 day clock and labels prescribed endocrine inputs. Native OpenSim inspection exposes engine-evaluated wrapped paths and forces; external force balance is not solved.

The BETSE source model renders 212 actual planar solver cells with 34 recorded frames. Select membrane voltage (V), junction gating (dimensionless), or native ion concentrations (mol/m³); fixed archive ranges preserve comparisons across time. Click a cell to inspect its values. These are generic computational tissue cells, not registered human skin. Transmembrane voltage differs from the extracellular wound field (V/m).

The scenario form supports upstream patient selection and a bounded hemorrhage–saline protocol with explicit stop actions. StandardMale completed the 2-second browser integration check. The original upstream StandardFemale build fails during native initialization, including without interventions; failed runs remain visible. This limitation is labeled when the upstream implementation is selected. Explicit corrected implementations repair the source bounds defect; their provenance remains distinct from upstream. Short protocol checks validate transport and execution, not physiological recovery.

CSF trajectory choices expose the literature-model baseline, hypotension ramp, and native MAP-driven run. The MAP coupling is one-way, with unmatched source subjects and published rounded initial ICP/compliance; no ICP feedback enters BioGears. The chart retains the trajectory's declared pressure, volume, compliance, and flow units.

Optional ambient temperature (10–35 °C) and clothing (0–3 clo) inputs are omitted from requests when blank. Explicit zero clothing is preserved. These are environment overrides; they do not establish thermal comfort or posture-specific physiology.

The native implementation selector uses available API descriptors and defaults to the verified corrected implementation when installed. Requests include `engine_variant` explicitly. The bounds correction repairs an out-of-range source vector access; the heatflux variant additionally corrects diagnostic evaporation output. Neither constitutes clinical calibration.

Browser verification completed the corrected StandardFemale 2-second hemorrhage/saline protocol with `saturation_bounds_heatflux` (run `20260905-000949-80412525`), including completed status and recorded trajectory. Full production browser suite: nine flows passed; unit suite: ten tests passed.
