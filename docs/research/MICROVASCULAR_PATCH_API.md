# Read-only local muscle vascular API

`ihm.app.microvascular_patch.MicrovascularPatchService(root).materialize(data)` returns a bounded on-demand vascular patch in retained canonical body coordinates. It does not create a native actor, issue native commands, write files or add blood stores. Shared-server route wiring is owned separately; the adapter is intended for `POST /api/body/microvascular-patch`.

Minimal request:

```json
{"entity_id":"body-bp3d-FJ1442"}
```

Supported IDs are `body-bp3d-FJ1442` (right vastus lateralis) and `body-bp3d-FJ1442M` (left). Entity-only requests use explicit retained-atlas fixture points `[-0.14,-0.23,0]` and `[0.141294,-0.23,0]` meters. Every call rechecks winding and full-patch surface clearance; these are not assumed centroids or arbitrary UI-space points. Other selected muscles return structured `unsupported_region` errors. The response states whether a point came from this fixture or the caller.

Optional fields are `position_m`, `resolution_m`, `section_mm`, `extent_mm`, `seed`, `cohort`, `scenario_id` and `zoom_bounds_m`. Unknown fields are rejected, including attempted native volume or ownership overrides. Defaults are 50 µm sampling resolution, 0.2×0.2 mm transverse section, 0.5 mm longitudinal extent, seed 19 and `young_men`. All four measured density cohorts are available. Query position and zoom bounds use canonical meters, not current posed renderer coordinates; pose-to-reference mapping is a separate unresolved UI/runtime responsibility.

## Declared scenario and response

The default named scenario is `human_quadriceps_1988_corrected_v1`. Its human evidence is the shrinkage-corrected pair of quadriceps diameter means; `human_quadriceps_1988_uncorrected_v1` explicitly selects the uncorrected histological state. Both use the documented area-equivalent radius transform. Their fixed scenario values are returned verbatim: radius-shape CV 0.2, supply/return radius multipliers 2/2.4, Newtonian viscosity 0.003 Pa·s and supply/return pressures 4000/1000 Pa. Those values are declared engineering choices, not measurements from the diameter source. A simple renderer request thus produces definite geometry/flow while retaining its conditional interpretation. No silent radius-distribution or in-vivo calibration claim is introduced.

Response schema `ihm.microvascular-patch-response.v1` contains:

- `patch`: complete source-conditioned graph, stable patch identity, registration and containment evidence, plus `zoom_edges`.
- `scenario`: ID, human-readable label, complete hydraulic input and interpretation.
- `selection`: canonical entity, name, point and point-selection basis.
- `limits`, `response_units`, `read_only=true`, and `native_commands_issued=0`.

Each zoom edge exposes `edge_id`, `edge_index`, `source_nodes`, `source_t`, `endpoints_m`, `radius_m`, `endpoint_pressure_pa`, `flow_m3_per_s`, and `centerline_samples_m`. A local renderer can project these coordinates and draw physical-radius tubes. Zoom clipping is an analytic centerline view retaining the original edge's flow and pressure interpolation; no new physiological boundary is created. Changing resolution or ROI does not alter full graph physics or patch identity. The default right patch contains 13 capillaries and 61 graph edges.

Source/preparation/uncertainty is available through `patch.graph.conditioning`, `radius_conditioning`, `hydraulic_model`, source receipts and `patch.registration`. `native_owner` remains Muscle, `independent_store=false`, `native_volume_allocated_ml=null`, and macro-boundary correspondence remains null. Source priors condition the selected canonical body; they are not an alternate anatomical body.

## Resource and error contract

A shared nonblocking lock permits one query per process. Busy requests return status **503** immediately instead of queueing CPU work. There is no cache: every call revalidates source hashes, avoiding stale evidence after source changes.

The adapter enforces 32 KiB request JSON, at most 64 capillary crossings, at most 256 graph/zoom edges, at most 20,000 sampled points, and at most 2 MiB serialized response. The tighter edge count may reject a request below the capillary bound. Preflight uses a conservative full-patch length bound to reject excessive sampling before geometry generation; postflight independently checks actual counts and serialized bytes. Resolution is 5 µm–1 mm and longitudinal extent is at most 2 mm. This is local detail, never whole-body generation.

`PatchRequestError` exposes `.status`, `.code`, `str(error)` and `.as_dict()` for HTTP integration. Codes distinguish invalid JSON/query (400), unsupported region or containment (422), request/geometry/response budgets (413), changed source identity/schema (409), unavailable evidence (503), and busy (503). Returned source errors do not expose raw filesystem paths or tracebacks. The HTTP server must still enforce its own request body-size limit before JSON parsing and its existing local-origin guard.

## Verification

`OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_microvascular_patch_api.py` passes six tests: actual right/left canonical responses, unsupported/invalid/oversized preflight without invoking geometry, structured source/schema/containment failures, nonblocking busy behavior and output geometry caps. The inherited five embedding tests and eleven prior/hydraulic tests cover containment, stable identity, units, moments and network conservation. No native engine or browser is required for these service tests.
