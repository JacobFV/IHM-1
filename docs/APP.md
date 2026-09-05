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
