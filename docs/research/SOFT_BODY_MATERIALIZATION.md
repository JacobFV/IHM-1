# Soft-body materialization direction

Owner decision, 2026-09-06: the declared articulation topology is the wrong
abstraction. The body should be true soft body interacting throughout, with
bones as rigid inclusions or very stiff soft material, and joints emergent from
tissue, cartilage contact and ligament constraint rather than declared degrees
of freedom.

## Why the declared topology fails, measured

The plant is one 27.7 kg rigid `torso` carrying skull, jaw, ribcage, scapulae
and the whole spine, joined to the pelvis by a single 3-DOF `back` joint that no
muscle in the default 92-muscle set crosses and no passive law restrains. The
supine equilibrium search fails for that reason and not a numerical one: whole
body balance is already 761.360 N support against 761.376 N weight, while a
5.94 N.m gravitational moment sits on a coordinate with nothing to react
against. Worst joint acceleration 5.688-5.904 rad/s2 against a 1e-4 criterion.

Canonical deformation is one affine transform per entity driven by seven volume
bindings; 8 of 2408 entities ever leave identity, and `mechanics.py` sets
rotation to identity and never assigns it again.

## Implicit model versus materialization

The implicit model holds every entity, its geometry and its evidence at full
fidelity. A materialization (`ImplicitHuman.materialize`) is a specific
executable model derived from it for a specific question, selecting its own
subset, resolution and physics. Being in the implicit model does not mean being
in every materialization.

Consequences: there is no global element budget, no single representation
crossover, and no materialization that instantiates the whole atlas at full
resolution. The retained regional domains are materializations, not failed
integrations.

The real gaps are narrower. There is no whole-body *soft* materialization, and
materializations do not compose: every regional domain is detached and records
it (`whole_body_force_feedback: false`). A local implant domain is only
meaningful if a coarse global one can drive its boundary.

Materialization classes in view:
- whole-body soft, coarse, major structures, no capillary bed, emergent joints
- regional high fidelity, a limb segment or organ with local vasculature
- integument bioelectric, epithelial sheet, Vmem over gap junctions, minutes to
  days, no musculoskeletal system
- implant / foreign body, an implant plus a few cm3 of tissue and its capillary
  bed, graded from micron scale at the implant surface out to millimetres

## Conforming meshes, not voxels

Per-entity conforming tetrahedral meshes are the representation. Measured
reasons:

- Thin structures survive. 856 surfaces have a 2V/A radius proxy under 1 mm,
  reaching 0.125 mm, and 855 tetrahedralize; face areas to 6.7e-34 m2 and edges
  to 1.4e-17 m. 1016 entities claim no voxel at 8 mm. Vanishing is a property of
  the grid, not of the anatomy.
- Boundaries are where the physics is. Median conforming volume drift 3.9e-15.
- Sliding is only expressible per entity. A monolithic mesh with shared nodes
  welds every tissue to its neighbour, forbidding muscle-on-muscle sliding and
  skin gliding over fascia. That reproduces the rigid-block failure.

Coupling therefore needs an explicit answer. A conforming interstitial matrix
binding tied tissue, with genuine contact only where anatomy licenses sliding.
The extended atlas marks those interfaces: 78 bursae, 40 tendon sheaths, 38
fascia.

## State as of 2026-09-06

Proven:
- 2396 of 2408 entities carry a volume-gated TetGen proof. Muscle 424/425 at
  4,643,624 tets; the rest 1971/1978 at 11,840,747 tets.
- Repair order: exact weld, drop repeated-index and collinear faces, collapse
  coincident opposite-winding fin pairs, orient per patch and flip inward ones,
  CGAL remesh_self_intersections, then CGAL self-union where non-manifold edges
  remain. Self-intersection and fin faces explained every pre-repair failure.
- Muscle is the pathological tissue: 44.7 percent self-intersecting against 3.4
  percent elsewhere, control 33.9 percent against 94.6 percent. The defect
  concentrates in sheet-like fascial and ligamentous tissue, so ligaments
  (10/20 to 20/20) and tendons (1/4 to 4/4) are where repair earns most. This
  matters because those classes are the ones intended to grow.
- Watertight genus-0 outer envelope, 69.720 L, 1.78125 m2. The canonical skin is
  a double-sided slab of 3.5026 m2, about twice a body surface.
- 46.08 percent of the interior belongs to no entity: 32.117 L, of which one
  connected interstitial component holds 30.907 L. No adipose geometry exists.
- Entity overlap is 1.080 per occupied voxel with 0.46 percent claimed by three
  or more, so exclusive ownership is a small obstacle.
- No genuine anatomical cavity exists to destroy: BodyParts3D pairs a wall
  entity with a lumen entity rather than nesting shells.
- Joint constraint network exists and is registered. 301 ligaments, 36 capsules,
  32 discs, 12 menisci, 4 cruciates, 4 labra; affine plus thin-plate spline at
  6.18 mm never-fitted centroid RMS against a 36.2 mm rejection precedent.

Open:
- Zero articular cartilage in either atlas, and no capsule encloses a volume.
  1012 offset shells over 0.334 m2 would be synthesized.
- The bone pose is not a consistent articular configuration: gaps run from
  interpenetration to 5.03 mm, so no single cartilage thickness works.
- The 30.907 L interstitium must be synthesized. It is the coupling matrix, not
  merely missing volume.
- Materials: two literature sources cover 2408 entities, nu 0.45 and rho 1000
  are blanket assumptions, bone has no modulus.
- libigl cannot enter the main venv on this aarch64 host; upstream ships no
  aarch64 abi3 wheel, so TetGen and CGAL are out of process.

## Intended uses driving the requirements

- Michael Levin style bioelectricity theses on the integument. A conforming
  mesh's node adjacency is a gap-junction network for free. Every per-area
  bioelectric quantity currently uses the 2x-too-large slab area.
- Microvasculature integration with an implant or other fine structure, needing
  graded local refinement and materialization composition.
- Growth in tendons, ligaments, vessels, integument and nerve fibres, all using
  the same conforming approach.
