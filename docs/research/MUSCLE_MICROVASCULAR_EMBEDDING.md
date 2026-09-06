# On-demand muscle vascular geometry in canonical body coordinates

`ihm.assembly.muscle_microvascular_embedding.materialize_patch` now takes a canonical muscle entity, physical query position and display sampling resolution, and returns an embedded vascular graph with a local zoom view. It reuses the human-prior density/radius materializer and its declared hydraulic scenario. No whole-body graph, texture or replacement anatomical body is generated.

```python
from ihm.assembly.muscle_microvascular_embedding import materialize_patch
scenario = {
    'radius_mode': 'human_quadriceps_1988_equal_area',
    'diameter_state': 'shrinkage_corrected',
    'capillary_radius_cv': .2,
    'supply_radius_multiplier': 2.,
    'return_radius_multiplier': 2.4,
    'viscosity_pa_s': .003,
    'pressure_boundaries_pa': {0: 4000., 1: 1000.},
}
patch = materialize_patch(
    '.', entity_id='body-bp3d-FJ1442',
    position_m=[-.14, -.23, 0.], resolution_m=50e-6,
    seed=19, cohort='young_men', section_mm=(.2, .2),
    extent_mm=.5, hydraulic_scenario=scenario,
)
```

`patch['graph']` contains the complete graph: canonical-meter endpoints, original edge indices/kinds, physical radii/lengths, geometric volume, resistances and predicted pressures/flows. `patch['zoom_edges']` contains stable edge IDs, endpoint arrays, physical radius, interpolated endpoint pressures, original-edge parameters and centerline samples. The app can render tubes directly. Physical radius is never enlarged to match screen pixels. An omitted hydraulic scenario retains explicitly unresolved radii/flow, as in the source materializer.

## Registration and containment

Only left/right vastus lateralis are currently accepted. Priors remain selectable conditioning on this same authored body structure, with no implicit cross-muscle or cross-species transfer. BodyParts3D geometry is read from the retained canonical meter-frame source and verified against its SHA-256. The source license/specimen, anatomy hash, source/generator hashes, requested position and rigid transformation accompany the result. The local longitudinal axis aligns to the atlas principal axis; that alignment is an engineering inference, **not measured muscle-fiber architecture**. Rotation and translation preserve radius, edge length, transverse density and topology. The physics is solved before rigid embedding and remains invariant under that transformation.

Canonical vastus-lateralis metadata says raw edge incidence is not watertight. Exact-coordinate welding resolves its duplicate-coordinate seams without moving a vertex or filling a hole. The new module independently verifies that every welded edge has exactly two incident faces with opposite directions. It requires nondegenerate triangles. It does not use a bounding box as occupied tissue: the right VL's recorded bounding-box centroid is actually outside its winding-defined interior. The example query `[-.14,-.23,0]` m has winding 1 and **0.01246640805435 m** nearest-triangle clearance.

Containment uses signed triangle solid angles for winding and analytic point-to-triangle distance for surface clearance. The patch's maximum local node distance plus its maximum realized radius bounds every straight edge cylinder inside a sphere. Requiring that sphere to lie strictly within the center's surface clearance, with a 1 nm numerical guard, certifies that the entire patch cannot cross the source triangles. Boundary-near or oversized requests are rejected; the physiological graph is never silently cut to fit. This sufficient condition is conservative: some admissible elongated patches can be rejected.

These are numerical geometric checks against authored atlas surfaces, not exact-arithmetic proofs or confirmation of exclusive biological tissue occupancy. Global self-intersections, neighboring-muscle overlap, internal material composition and actual terminal-vessel correspondence remain unresolved. A nonunit winding or surface point fails closed.

## Zoom is a view of one graph

An optional `zoom_bounds_m=[lower_corner, upper_corner]` clips **centerline segments** using analytic slab intersections in original-edge parameters. It does not voxelize the graph or generate new branches. Radius remains attached to the original edge; a renderer that requires clipping cylinder surfaces at the view boundary must apply its own geometric clipping plane. No cylinder-solid intersection claim is made by centerline clipping.

Predicted pressure at a clipped endpoint is interpolated linearly along the original constant-radius/resistance edge, consistent with the declared steady cylindrical model. Original edge flow is retained. A zoom cut is a display boundary, not a new physiological inlet or blood store. `source_t` preserves the exact segment parameter interval up to floating-point arithmetic.

Changing resolution or zoom bounds does not change patch identity, topology, radii, source conditioning or solved pressure/flow. Resolution controls only centerline sampling, bounded to 100,000 returned points. Patch IDs include geometry/anatomy/prior/generator receipts, position, seed and all conditioning; registry or physical-request changes produce a new identity. The underlying capillary count remains bounded to 128 by the source materializer.

## Ownership and limits

The result is a native-owner view with `native_owner='Muscle'`, `independent_store=false` and `native_volume_allocated_ml=null`. Geometric lumen volume describes the generated demand; no volume/mass/perfusion is added to or subtracted from native compartments. Actual conservative allocation requires a later explicit quantity-budget contract. `macro_boundary_correspondence=null` records that the scenario's supply/return nodes are not measured links to native or atlas macrovessels.

The diameter source remains human quadriceps histology, with selected shrinkage state and the uncertainty interpretation documented in [the prior materializer](ORGAN_MICROVASCULAR_MATERIALIZATION.md). CV, supply/return multipliers, pressure boundaries and Newtonian viscosity remain declared engineering choices. No local tissue exchange or native-flow feedback is added by embedding.

## Verification

`OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_muscle_microvascular_embedding.py` passes five bounded checks: analytic segment-box clipping, tetrahedron winding/distance and open-surface rejection, retained canonical registration with rigid physical invariance, zoom identity and pressure interpolation, and outside/oversized-patch rejection. The tests use one 6,688-face retained muscle and a 13-capillary graph, and finish in about one second. No native simulation, build, browser or bulk generation is required.
