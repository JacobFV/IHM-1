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
  merely missing volume. Composition is now sourced: adipose is 60.7 percent of
  it, the fill weighs 31.7461 kg, and the body totals 70.7713 kg against a
  declared 77.1107029 that no soft-tissue density can reach.
- Two infrapatellar fat pads in the joint promotion candidate are the first
  adipose geometry in the model, superseding the statement above that none
  exists.
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

## Cross-structure conflict and the mesher decision

Fusing several atlases means two sources can claim the same volume. Testing all
62,973 bounding-box-overlapping pairs of the 2,403 repaired surfaces with CGAL
finds 11,047 genuinely intersecting, 2,312,750 intersecting face pairs, and
2,345 of 2,403 entities in conflict with only 58 clean. Total pairwise overlap
is 3.1367 L, 8.31 percent of summed entity volume, and thin: median 12.9 mm3,
thickness proxy median 0.489 mm.

The conflict graph has two connected components and one holds 2,340 entities.
There is no decomposition into independent neighbourhoods.

Resolution is by declared role priority (rigid_bone, cartilage, tendon,
ligament, fluid_cavity, vascular, nerve, muscle, soft_organ, connective_tissue,
lymph_node_group), ties to larger volume, with conflicts graded rather than
uniformly differenced: 4,269 segmentation noise at 0.112 L, 4,250 substantive,
2,528 flagged as genuine inter-source modelling conflicts at 2.196 L, and 362
neurovascular peers where neither displacement is justified. Five pairs are
byte-identical surfaces authored twice under synonymous BodyParts3D ids and are
duplicate entities rather than conflicts.

Method: one exact CGAL arrangement imprints every crossing curve into both
surfaces, coincident facets collapse to a single shared interface, and tets are
labelled by winding number with disputed tets going to the priority owner. Six
thigh structures that crashed TetGen now mesh with zero warnings at 112,576
tets, box closure 0.0, per-structure volume error at most 1.83e-07 and residual
pairwise overlap of exactly 0.0 over all 15 pairs; 48 entities reach 1,060,964
tets with 0.0 residual overlap over 1,128 pairs.

### The ceiling is a TetGen defect, not geometry

An earlier diagnosis blaming sub-tolerance features was wrong, and wrong because
arrange() overwrites its report each snap-loop iteration so the stored PLC was
post-snap. Measured on the exact arrangement, 48, 50 and 60 entities are
indistinguishable: minimum edge 8.897890141411458e-09 m, exactly 4 edges under
tolerance, minimum facet area 1.062084714469123e-15 m2, zero residual
self-intersections in every case, all present already at 6 entities. The exact
route meshed 60 entities once and failed on a byte-identical PLC on re-run and
in six repeats, with glibc corrupted-size aborts and SIGSEGV. Deterministic
within a process, not across them. No PLC-cleaning route can lift it.

Snap rounding is rejected on measurement: it lowers the ceiling from 48 to 12.
CGAL ships no 3D snap rounding, and one quantisation at 1e-7 m drives minimum
facet altitude from 2.100e-09 to 3.271e-18 m while creating crossings in PLCs
that had none; the fixpoint loop diverges 57 to 10,130.

fTetWild is built at data/runtime/tolerant-mesher/ after three aarch64 blockers:
GMP 6.3.0 from source, a patch selecting geogram's Linux64-gcc-aarch64 instead
of the hardcoded Linux64-gcc with its -m64, and -fsigned-char for Mesh.hpp's
negative sentinels under an unsigned char.

Both routes preserve conforming shared nodes at interfaces, measured
independently of construction, so that property does not decide. Element quality
does. TetGen -pY leaves 38 percent of tets under 5 degrees, a minimum dihedral
of 0.000 degrees and one negative-volume tet, and any quality pass that fixes
that moves the nodes and surrenders the exactness that was its only advantage.
fTetWild reaches 100 entities and 526,645 tets with no tet under 5 degrees,
paying 2.08e-04 m of interface position.

Decision: fTetWild for the whole-body coarse materialization, the exact route
for regional domains of at most 48 entities. Open: fTetWild sizes elements
against the bounding-box diagonal, so at the 1.866 m whole-body diagonal an
epsr of 1e-3 is a 1.87 mm envelope, wider than the cartilage the atlas carries;
4.89e-06 m would need 2.6e-6, outside its design regime and unmeasured. Nothing
above 100 entities has been measured.
