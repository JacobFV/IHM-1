# Next bounded clothing-mechanics step

The next implementation should close one demonstrable collision gap before activating the entire shorts mesh: **continuous finite-mass edge/edge impact with Coulomb friction**, integrated into a separately versioned garment fixture. The present face-interior node/triangle solver remains the retained baseline. Its validated results must not be silently regenerated or relabeled.

## Current evidence and missing case

`contact_dynamics.py` computes nearest edge/vertex distances but explicitly counts such cases as unresolved. It only resolves a nearest interior triangle projection, after the candidate step, and has no swept or self collision. `garment_contact.py` uses two-way face-interior contacts between a 135-node extracted front panel and the source tissue, with prescribed panel boundary support. The full shorts' shared crotch indices are present in the generated mesh, but that alone does not establish full-garment contact, waistband support or containment.

A real negative probe is retained at `artifacts/verification/clothing-next-step/edge-miss.json`: a 1 g node approaching just outside a triangle edge, 0.141 mm from that edge, returns zero contact impulses and one unresolved-edge contact. The proposed 0.3 mm proximity band in that probe is an engineering collision parameter, not a measured textile thickness. It demonstrates the missing feature; proximity by itself does not authorize an early impact or a physical-thickness claim.

## First implementation contract

1. Add a pure segment/segment swept-event kernel in a new module. Inputs are both edges' old positions, constant trial velocities, nodal masses/mobility and the time interval. Find the first **zero-thickness** intersection over the interval, require both segment parameters in `[0,1]`, and return the time and barycentric endpoint weights. Use conservative candidate bounds and verified roots; explicitly reject or report unresolved degenerate, parallel/coplanar and simultaneous-event cases. An unresolved case must not count as certified containment.
2. At the coincident contact point, apply finite-mass impulses to both edges. For barycentric weights `a_i,b_j`, the relative inverse mass is `w = Σa_i²/m_i + Σb_j²/m_j`, excluding fixed nodes. For zero restitution, `j_n = max(−v_n,0)/w`. A tangential impulse first tests the static Coulomb cone; a sliding impulse uses the kinetic coefficient and cannot reverse tangential motion. Return equal/opposite nodal impulses, fixed-support reactions, kinetic transfer and dissipation. Contact occurs at a common point, making angular-momentum checks meaningful without inventing rotational degrees of freedom for point masses.
3. Add a separately versioned transactional common-clock stepper or explicit optional port after reviewing the kernel. Advance to the event, apply the impulse, then advance the remaining interval; recompute or consistently subdivide elastic forces. Record numerical integration work/energy defect separately from collision dissipation. Do not fake the event by pulling noncontacting edges together inside a search band.
4. Run two small dynamic cloth patches whose crossing edges intersect while their vertices miss the other patch's interior projections. Compare contact-on, contact-off and friction-off, with two step refinements and actual stored trajectories. Then run a local shorts-panel edge/grazing fixture against a mobile source-tissue boundary using the same new contact path. Retain source geometry, law, configuration and code hashes for each fresh experiment.

Proposed files: `ihm/assembly/swept_edge_contact.py`, a separate versioned coupled-step module if needed, `scripts/verify_swept_edge_contact.py`, and a fresh patch/garment experiment script. The first deliverable is the mechanical kernel and missing-crossing regression; full-shorts activation follows only after its receipts are accepted.

## Exact mechanical acceptance example

Two perpendicular edges meet at their midpoints. Each of their four endpoint masses is 2 g; the first edge approaches normally at 0.1 m/s and the second rests. The effective inverse mass is 500 kg⁻¹, normal impulse is 0.0002 N·s, all endpoints leave with common normal speed 0.05 m/s, and normal kinetic loss is 10 µJ. This is an analytic synthetic fixture, not a measured garment coefficient.

Tests must also cover endpoint contacts, a high-speed crossing missed by endpoint-only sampling, reordered endpoints, both orientations, unequal masses, a fixed target with recorded reaction, static sticking and kinetic sliding. Required checks are nonnegative dissipation, equality of kinetic loss and computed contact loss, equal/opposite total impulse, angular impulse at the common event point, correct event time and barycentric weights, no event outside the interval, and atomic rejection of unsupported geometry. A no-contact case must preserve state and report no contact, rather than project to an arbitrary plane.

## Subsequent gates for actual full shorts

- Add swept vertex/triangle events and their edge/vertex limits, and a unified earliest-event schedule. Edge/edge alone does not close every collision case.
- Add self-contact for nonadjacent cloth features, with explicit topology exclusions for shared vertices, sewn seams and immediate material neighbors. Test a fold and two crossing nonadjacent garment panels; seam neighbors must not repel themselves.
- Audit the full 2,527-node, 4,864-face shorts for three intended openings, seam connectivity, positive masses and mechanical edge/rest-length correspondence. Crotch connectivity exists; elastic waistband rest lengths, thickness, bending and material anisotropy remain explicit measurements or assumptions.
- Replace fixed front-panel boundaries with a full-garment support model: elastic waistband and actual pelvis/hip/thigh contact, measured or explicitly assumed skin/textile friction, finite-mass tissue reactions and source tissue supports. Maintain exclusive body/material ownership.
- Only then test load/release containment over all relevant surfaces, with no unexplained intersections, positive tissue determinants, friction-off and support-off contrasts, force/work receipts and resolution checks. Local penile deformation does not by itself establish full shorts containment or physical tucking.

This sequence follows the collision/contact/friction separation emphasized by Bridson, Fedkiw and Anderson's [primary cloth-collision paper](https://graphics.stanford.edu/papers/cloth-sig02/). The proposed first kernel is narrower than that complete algorithm; this reference supplies a method, not calibrated fabric coefficients or a validation of the IHM body.
