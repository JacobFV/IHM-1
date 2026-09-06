# Whole-garment exterior selection and tunneling guard

The default shirt and shorts are visible, but `app/src/clothing.js` updates them using respiration deformation and a rigid entity transform. It does not consume the optional Python `garment_mechanics` frames. This change improves that optional mechanical solver; it does not establish an accepted full-body trajectory or change the viewport.

## Native and display ownership

`Cloth` owns finite lumped masses, objective edge-spring energy and symplectic Euler integration. It has no bending or calibrated fabric law. `GarmentFeedback` subcycles those garments against a moving source surface, computes barycentric tether/contact reactions and iterates work-conjugate native body loads with a common checkpoint and rollback on failure. `ArticulatedPlant` enables it only when explicitly requested. Current small feedback tests use a finite toy body; they do not validate a complete supported human trajectory.

Previously `GarmentFeedback.from_root` passed every retained raw skin triangle to contact and chose tethers from that raw surface. Raw inner components could therefore influence contact. Shared `MovingSurfaceContact` accepts only discrete nearest face-interior nodal contact, counts unsupported nearest-edge events, and limits candidates by an endpoint distance. A node can cross the skin during a substep and end deeper than its 4 mm search threshold, producing no endpoint correction.

## Exact exterior materialization

`garment_exterior_contact.exterior_source` calls the existing `physical_skin_support` contract on exact source geometry bytes/hash and reviewed component evidence. It filters triangle connectivity by the evidence's original `contact_eligible_triangle_ids`. Source vertex coordinates and IDs remain unchanged; `contact_source_face_indices` maps every local contact face back to the source. The full native body masses and ownership are unchanged.

Existing support/sample garment node identities are retained, but their nearest source locations are requeried against exterior faces. This is explicitly a **new engineering tether candidate**, with its complete previous raw-source pairing retained, an explicit mapping-status label, and a new mapping hash. It is not claimed to preserve the former tether's attachment, rest distance or energy. This occurs at materialization; it cannot silently retether an already running garment.

The inferred exterior component still lacks surface-exclusivity, self-intersection and closure validation. Nearest-bone assignment and skin continuity across articulated segments remain unresolved. Exact filtering does not certify those properties or prove wearable containment.

## Fail-closed linear swept guard

`garment_swept_vertex_face.swept_vertex_face` solves the cubic coplanarity polynomial for a linearly moving point and linearly moving triangle. At an interior root it returns normalized interval time, barycentric contact coordinates, oriented normal, world point and crossing direction. Edge/vertex contact, coplanar intervals, tangencies, near-multiple roots, intermediate triangle degeneracy and multiple crossings are unresolved. Tolerances are explicit numerical guards; this is not interval-arithmetic CCD certification.

The new `StrictGarmentSurfaceContact` reconstructs the actual `Cloth.step` kick/drift endpoints, interpolates body endpoints over that substep and uses a conservative swept-sphere/AABB candidate search with a bounded pair budget. Any detected crossing or unresolved event **raises before contact acceptance**. The feedback transaction restores cloth and native state. Error records include the source face ID for source-bound materializations. A nearest-edge event returned by the shared discrete contact solver also raises. The shared solver used by hair is unchanged.

This is a tunneling guard, not a completed impact integrator. It deliberately cannot proceed through first contact until an event schedule resolves the impact and remaining substep. It does not certify cloth edge-edge, triangle-triangle, self-contact, finite thickness, initial containment, persistent static friction or non-linear substep trajectories. Initial small penetration still follows the existing endpoint projection and reports its numerical energy defect; no invisible heat correction is introduced.

## Bounded validation

`scripts/verify_garment_exterior_contact.py` runs nine tiny tests covering:

* byte/hash-bound two-component source filtering, original face mapping, preserved prior tether record and a distinct candidate mapping;
* actual `GarmentFeedback.from_root` wiring with mocked tiny source owners;
* fixed and moving face crossings, barycentric reconstruction and rigid-transform invariance;
* outside-face clearance, edge contact, coplanarity, interior tangency and intermediate triangle degeneration;
* deep endpoint tunneling rejection, candidate-budget rejection and a clear interval;
* strict feedback rejection with exact cloth/native checkpoint restoration.

The existing `verify_garment_feedback.py` remains passing. No native execution, browser action, full garment materialization, large geometry scan or physical full-body contact acceptance was performed for this change.
