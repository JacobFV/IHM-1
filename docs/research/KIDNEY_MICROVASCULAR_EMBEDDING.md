# Kidney patches in the bounded local-detail API

The existing `MicrovascularPatchService` now accepts kidney selections through the same read-only request/response contract as muscle. Existing muscle requests, scenario IDs and numeric budgets remain supported.

```json
{"entity_id":"body-bp3d-FJ3147","scenario_id":"human_kidney_control_v1","resolution_m":0.00005}
```

Right kidney is `body-bp3d-FJ3147`; left kidney is `body-bp3d-FJ3145`. Supported kidney scenarios are `human_kidney_control_v1` (default) and `human_kidney_injury_v1`, selecting the respective human biopsy-region PTC mean-area transform. “Control” refers to preserved regions in the source biopsies, not a separate healthy donor body. Pressure, viscosity and hydraulic-exchange assumptions remain declared in the response. Kidney requests reject muscle-specific cohort, section and extent fields instead of applying muscle parameters to kidney.

The selected region remains the same retained canonical body. The anatomical source hash is checked, and exact-coordinate welding/edge incidence validates geometric surface closure. Explicit default query points are right `[-0.07,0.23,-0.012]` m and left `[0.075,0.25,-0.02]` m. They have approximately 11.266 and 11.501 mm nearest-surface clearance respectively. The left kidney's bounding-box centroid is outside its winding-defined interior and is not used. Every request independently recomputes containment; a caller may supply another canonical `position_m`.

`materialize_registered_kidney` in `ihm/assembly/kidney_microvascular_embedding.py` rotates the kidney-specific source-conditioned graph using the atlas principal axis and translates to the query point. The surface-clearance certificate bounds **all polyline points plus the maximum vessel radius**, so it covers the winding filtration capillaries, not just their endpoint chords. Every polyline segment lies in the same conservative enclosing ball. Outside or uncertifiable requests fail with a structured API error.

This establishes geometric containment in the authored whole-kidney surface only. **Cortical location is not verified**, and the API explicitly returns `cortical_location_verified=false`. Native terminal correspondence, glomerular afferent registration and the actual donor's nephron orientation are unresolved. `body_registered=true` describes the declared rigid placement in canonical coordinates, not recovered individual microanatomy or a validated cortex/medulla localization.

## One polyline-aware graph contract

The response wrapper remains `ihm.microvascular-patch-response.v1`, with `patch`, `scenario`, `selection`, `limits` and units. `patch.graph` includes the original kidney geometry, physical-node pressures, physical-edge flows and separate full hydraulic/exchange diagnostic. The common `graph.flow_solution` exposes physical blood-node/edge arrays for generic consumers; named Bowman/interstitial boundaries remain in the kidney-specific fields.

`zoom_edges` preserves original bends in `centerline_samples_m`. A kidney edge must be rendered as that polyline; connecting only its two endpoint coordinates would discard its length-conditioned geometry. Additional samples subdivide long straight polyline segments to the requested maximum length without deleting original bends. Resolution changes leave full graph positions, radius, lengths, topology, pressures and flows identical.

Optional axis-aligned zoom clipping applies analytic segment-box intersections to every polyline segment. Reentry creates separate display fragments; contiguous clipped segments merge. `edge_index` and `source_edge_id` retain the physiological edge identity; `fragment_index` distinguishes views. `sample_source_t` stores fractions of **original cumulative arc length**. Clipped endpoint pressures interpolate using arc length, consistent with a constant-radius Poiseuille edge. Every fragment carries the original flow and radius. Display cuts create no new physiological boundary or storage pool.

Patch IDs include entity, position, scenario, seed and geometry/anatomy/prior/generator hashes, excluding resolution and ROI. Thus source or physical-conditioning changes invalidate identity; display changes do not.

## Bounds and ownership

Kidney generation is fixed at 76 physical nodes, 122 physical edges and 77 filtration edges. The service declares a separate `max_kidney_filtration_edges=96` budget; the existing muscle 64-capillary bound remains unchanged. Both organs share limits of 256 output edges/fragments, 20,000 samples and 2 MiB serialized JSON, with one nonblocking process-wide query lock. Excessive clip fragmentation fails instead of truncating geometry. The actual default right-kidney response contains 122 zoom rows and 10,061 samples, about 1.68 MB under ordinary JSON serialization.

`native_owner` is LeftKidney or RightKidney, `independent_store=false`, `native_volume_allocated_ml=null`, and `macro_boundary_correspondence=null`. No native perfusion, filtration, solute/charge transport or duplicate blood volume is created. Source preparation, modeled filtration radius and engineering exchange assumptions remain inherited from [the kidney prototype](KIDNEY_MICROVASCULAR_MATERIALIZATION.md).

## Verification and integration

Three direct embedding tests verify polyline exit/reentry clipping, arc-length parameters, rigid length/physics preservation, stable identity and unsupported-position rejection. Nine service tests retain all existing muscle tests and exercise both kidney sides, preserved/injured source conditioning, output size and cross-organ-field rejection. Six kidney-model tests continue to check source constraints, serial beds and hydraulic conservation.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_kidney_microvascular_embedding.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_microvascular_patch_api.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_kidney_microvascular_patch.py
```

The HTTP route remains the existing `POST /api/body/microvascular-patch`. The frontend worker received a real default response in `/tmp/ihm-kidney-microvascular-ui-fixture.json` for polyline renderer verification. This increment does not edit the shared HTTP server or frontend files.
