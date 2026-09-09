# what the body touches the world with, measured

The upright plant's non-foot contact is `fall_proxy_<body>`: one sphere per
segment whose radius is inscribed in that segment's **inertia ellipsoid** and
whose centre is its mass centre. A femur represented as a ball. This file
records what replacing that costs, what the solver does with a concave surface,
and — the part that matters most — **why skin-mediated ground contact cannot be
switched on yet, with the number that says so.**

Everything here is from `scripts/build_segment_contact_meshes.py`,
`scripts/build_skin_contact_meshes.py` and
`scripts/measure_segment_contact_meshes.py`, over the
`engineering_stance_v1` pose at 77.6122029 kg.

## The mechanism

`ContactMesh` over `SimTK::ContactGeometry::TriangleMesh`, carried by
`ElasticFoundationForce`, one force per (mesh, floor) **pair** — never one set
over all the meshes, because an OpenSim contact set tests every pair inside it
and a shared set would run mesh-against-mesh between neighbouring segments and
invent a bone-on-bone force at every joint the source model lets overlap.

A bundle declares its **layer**. `bone` is the skeleton, and is a collider
against other bones and its own soft tissue. `skin` is what actually meets the
floor. The layer travels into the native emit as `segment_contact_mesh_layer`
so no report can say the body stood on its skin about a run that stood on its
femurs.

## Concavity: nothing here takes a convex hull

Read out of the sources rather than assumed.
`CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::processObjects` walks the
mesh's **own OBB tree** and returns the real face indices below the plane;
`ElasticFoundationForceImpl::processContact` then loops those faces and puts an
independent spring at each face centroid, with that face's own area. Concavity
survives end to end.

That matters here by measurement:

| surface | volume / its own convex hull | max depth inside the hull |
|---|---:|---:|
| rib cage + scapulae (`hat_ribs_scap`) | **0.043** | 60.4 mm |
| jaw | 0.168 | 22.0 mm |
| spine | 0.189 | 44.0 mm |
| femur | 0.213 | 28.6 mm |
| foot (`r_foot`) | 0.322 | 26.7 mm |
| pisiform (roundest bone in the body) | 0.997 | 0.2 mm |

All 81 bone meshes are concave. **The controls that make those numbers mean
something**: an icosphere and a box, measured by the same code path, print
1.000 and 0.0 mm.

## SimTK refuses 13 of the 81 bones OpenSim ships with this model

`ContactGeometry::TriangleMesh`'s constructor demands a **closed, consistently
oriented, non-degenerate edge-2-manifold**: no repeated vertex in a face, no
zero-area face, no two faces sharing the same *directed* edge, and every forward
half-edge matched by its backward twin. `simtk_precondition()` replicates that
check exactly from `ContactGeometry_TriangleMesh.cpp`. A surface that fails it
is not contact geometry at all — the engine throws while building the model.

61 of 81 pass untouched. Welding, rewinding and hole-filling recovers 7 more.
**13 cannot be admitted**: both pelvis halves, both tibiae, the spine, the skull,
the jaw, and six hand bones, all for `two faces share the same directed edge`.

Every segment keeps at least one bone, but **the torso keeps only its rib cage**
— the skull and the spine are absent from contact. That is said here rather than
discovered later as a body that falls through its own head.

(The first version of this gate compared only the *counts* of forward and
backward half-edges, which SimTK also does — and then SimTK checks each edge
individually. Meshes passed the gate and the engine threw on them anyway. The
gate now compares the half-edge **sets**.)

## What it costs

0.5 s from the stance pose, 50 advances of 10 ms, single-threaded, on a shared
and busy machine. `s / advance` is wall clock for one error-controlled Simbody
integration; its cost is not bounded by dt.

| arm | contact elements | faces carried | vertical contact force | pelvis_ty drift | s / advance median | worst |
|---|---:|---:|---:|---:|---:|---:|
| `spheres` — today | 28 | 0 | **761.38 N** | 0 | 0.058 | 0.070 |
| `mesh_proxy` — 16 inertia proxies → 62 bone surfaces | 74 | 27,734 | **761.38 N** | 0 | 0.067 | 0.114 |
| `skin_carried` — all 20 skin surfaces, feet still on spheres | 32 | 111,986 | **761.38 N** | 0 | 0.058 | 0.064 |
| `mesh_all` — every segment on bone, source feet gone | 68 | 35,754 | 206 N, falling | −191 mm | 0.223 | 1.060 |
| `skin` — every segment on skin, source feet gone | 20 | 111,986 | 302 N, falling | −190 mm | 1.341 | 5.036 |

**761.38 N is m·g to the last digit** (77.6122029 × 9.81), and the
momentum-balance residual `m(a_com − g) − contact − external` stays at 1e-14 in
every arm. That residual is the gate that says the new forces are *summed*
correctly rather than merely printed; a contact term bookkept wrong shows up
there and nowhere else.

**Mesh contact geometry is free until it touches.** `skin_carried` holds 112,000
triangles — 4× the bone bundle, 4000× the sphere count — and costs 0.058 s
against the baseline's 0.058 s. The OBB broad phase prunes every mesh that is
not near the plane. Cost scales with faces **in contact**, not faces carried.
The 0.223 s and 1.341 s rows are not the price of the geometry; they are the
price of a plant that is collapsing, and should not be quoted as a mesh-contact
cost.

**Replacing the 16 inertia proxies with 62 real bone surfaces changes standing
statics by exactly nothing** — those bones stand 76 mm clear of the floor — for
1.15× the wall clock. On this gate, mesh contact is affordable.

## Standing the body on its bones does not hold it up

`mesh_all` is the arm that answers the question the brief was originally written
around, and the answer is negative. The body starts **15.4 mm above the floor**,
because a skeleton is not a foot; sinks 191 mm in half a second; ends 6.2 mm
*through* the floor with 206 N of support; and what load there is arrives
through the metatarsals rather than the heel.

The ~15 mm that the source foot spheres' radii encode is the heel pad and the
plantar soft tissue. Taking it away is not a smaller simplification than a
sphere — it is a different one. Contact with the world is skin, over fat and
muscle, over bone.

## The skin exists, is cut per segment, and is 96 mm too high

`scripts/build_skin_contact_meshes.py` cuts the canonical exterior skin —
109,183 triangles, 1.7805 m² — into one closed surface per segment:

* partitioned by the repo's own `continuous_surface_binding`, a graph-diffused
  skinning weight per skin vertex per segment; a triangle goes to the argmax of
  its three vertices' mean weight;
* mapped into each segment's own frame through the reference run's t=0 body
  transforms (the supine environment rotates **gravity**, not the body, so those
  are the neutral standing transforms);
* **capped**, because cutting an open surface leaves open pieces and SimTK
  refuses them. `trimesh.fill_holes` closes between 0 and 16 triangles on these
  cuts and leaves the loop open, so the caps are built explicitly: chain the
  unmatched directed edges into loops, fan each loop to its own centroid. 20 of
  22 segments come out admissible; `talus_l` and `talus_r` carry 8 and 49
  exterior triangles respectively, which is correct — the talus is an interior
  bone with no skin of its own.

The caps are invented surface and are reported as such: 1.0946 m² of cap on
1.7793 m² of real skin, concentrated at the waist (pelvis 0.243 m², torso
0.229 m²) where two segments meet and nothing outside the body can reach.

And then the arm fails, for a reason that is nothing to do with contact:

> **The sole of the foot sits 96 mm above the floor, and 84–95 mm above the foot
> bone inside the same segment frame.**

Per segment, skin minimum minus bone minimum along the segment's own y:

| segment | skin−bone (mm) | skin/bone centroid gap (mm) |
|---|---:|---:|
| toes_l / toes_r | +94.6 / +93.8 | 107 / 105 |
| calcn_l / calcn_r | +84.0 / +83.5 | 111 / 112 |
| tibia_l / tibia_r | +72.5 / +69.8 | 90 / 96 |
| hand_l / hand_r | +58.6 / +56.7 | 100 / 101 |
| pelvis | −141.1 | 78 |
| torso | −165.3 | 92 |

The toe **skin** runs from +0.083 m to +0.164 m in the segment frame while the
toe **bone** runs from −0.011 m to +0.009 m. The skin is entirely above the
bone, which is not a thing a body can do.

The cause is in the registration's own report and needs no inference.
`CanonicalRegistration.global_fit` says:

> *"Unweighted proper-rigid least-squares fit of 22 approximate source
> COM / canonical bone-envelope center correspondences; no scale fit"*,
> `rms_landmark_residual_m` **0.1234**, `maximum_landmark_residual_m` **0.3528**.

**123 mm RMS.** The canonical skin is not registered to this skeleton well enough
to touch a floor whose height is fixed by the model.

**Why nothing caught this before.** The supine surface foundation is built
through the same map, and it sets its support plane to *the skin's own minimum*
— `plane = min source-x of the eligible faces`. The floor follows the skin, so
the error is invisible there by construction. Upright, the floor is at y = 0
because the model says so, and the same error becomes a 96 mm hover.
`docs/SUPINE_SURFACE_SUPPORT_CURRENT.md` already recorded its shadow — "pelvis
39 mm, femurs 88–89 mm, tibias 112 mm above the reference plane" — as a property
of the geometry rather than of the registration.

**This is the blocker for skin-mediated contact, and it is separable.** The map
to fix is one 4×4; `binding.json` carries a second, *scaled* similarity
(`scale = 0.96303`) fitted differently, and the two are not the same map. Which
one is right, or whether either is, is a measurement nobody has made.

## What this engine can and cannot express

Asked directly, because it decides whether the *fully present participant* mode
is weeks or years.

**There is no deformable continuum anywhere in this stack.** Simbody is a rigid
multibody engine. Every compliant element in the plant is a **one-dimensional
spring law attached to a rigid carrier**:

| element | what deforms | what carries it |
|---|---|---|
| `SmoothSphereHalfSpaceForce` | a Hertz/Hunt-Crossley indentation scalar | a rigid sphere on a rigid body |
| `ihm_surface::Foundation` | a confined compressible neo-Hookean **column** per quadrature point | 21,382 rigid stations on rigid bodies |
| `ElasticFoundationForce` | one spring per triangle, normal only | a rigid `TriangleMesh` |
| muscle / tendon | fibre and tendon length | a path between rigid stations |

`surface_contact_manifest` is therefore **a rigid quadrature over a fixed skin
shape**, not a soft body: the skin's stations never move relative to their
segment, the columns never couple to each other, no volume is conserved, and
nothing shears. It is a Winkler bed, and it is also axis-locked — the header
carries a single `plane` scalar and the force is hard-coded to `(normal,0,0)`,
so it only works against the supine x-plane.

So "soft tissue deforms and mediates contact" is **not expressible in this engine
at all**, at any cost. It needs a deformable solver — FEM, MPM or position-based
— coupled to Simbody at the body-force level, or a different engine. That is a
new solver, not a parameter. The honest schedule note is that the *participant*
mode's soft-tissue clause is gated on that decision, and no amount of work on the
existing contact path reaches it.

What the existing path *can* reach, and what this branch now has the machinery
for: **skin as a rigid surface over a normal-compliant layer**, with the layer's
stiffness derived from the body's own declared skin (`E` = 3000 Pa, `ν` = 0.45,
`h` = 6.6 mm from epidermis + dermis + hypodermis, giving
`k = (1−ν)E/((1+ν)(1−2ν)h)` = 1.72 MPa/m). That is a real improvement on a
sphere and it is one registration away from working. It is not a soft body, and
must not be reported as one.

Note the layer's own arithmetic while it is here: 761 N spread over ~0.02 m² of
plantar skin at 1.72 MPa/m needs ~22 mm of indentation, against a declared layer
6.6 mm thick. **The declared skin alone cannot hold the body up**, which is the
same thing `docs/SUPINE_SURFACE_SUPPORT_CURRENT.md` found. The missing stiffness
is the fat and muscle between skin and bone — the interlayering the direction
names, and which nothing in the model currently carries a thickness for.

## Muscles, tendons and ligaments at real points

**Muscles.** The model declares 98 muscles and **364 `PathPoint`s** on 18 of the
22 segments (pelvis 68, tibia 37 per side, femur 35 per side, calcn 23 per side;
`talus_l/r` and `hand_l/r` carry none). Those are real attachment stations on
the segment frames.

But the engine calls `replacePathsWithFunctionBasedPaths` at load, and the
`FunctionBasedPathSet` holds **80** entries. So 80 of the 98 muscles run as
**polynomials in the coordinates**, not as paths through those points; the points
survive only through whatever the fit captured. The remaining 18 keep their point
paths. "Anchored at real points" is true of the model's declaration and only
indirectly true of the running plant.

**Tendons** are not separate objects here: each muscle is a musculotendon unit
with a tendon slack length and a series-elastic curve. There is no tendon
geometry and no tendon attachment distinct from the muscle's.

**Ligaments are a gap, and the shape of the gap is specific.**

* The anatomy carries **328 ligament entities**, every one of them with
  `reference_geometry` — real meshes — bound across all 22 segments (torso 63,
  hand 36 per side, calcn 29 per side, pelvis 20, tibia 19/18, …).
* The plant carries **zero** ligament force elements. `ForceSet` is 80
  `Millard2012EquilibriumMuscle`, 18 `Thelen2003Muscle`, 13 `CoordinateActuator`.
  No `Ligament`, no `BushingForce`, no `ExpressionBasedBushingForce`.
* What stands in for them is the opt-in `CoordinateLimitForce` joint stops, whose
  *limits* are the model's own declared coordinate ranges but whose stiffness,
  damping and transition width are stated engineering constants — 30 N·m/rad in
  `crawl.py` — and not measured ligament properties.

The reason the anatomy cannot supply attachment sites as it stands is structural,
not a matter of effort: **the binding assigns each whole entity to exactly one
segment**, by nearest-bone-group vertex vote, and a ligament by definition spans
two bones. 123 of the 328 have a *runner-up* segment different from the one they
were assigned, which is the binding itself reporting that they straddle a joint.
Median assignment coherence over the ligaments is 0.68, against 1.0 for an entity
wholly inside one segment.

So the two attachment ends could be **constructed** — split each ligament mesh's
vertices by which segment they vote for, take a centroid per side, and there is
your origin and insertion — but that is a construction, and nothing in the repo
would validate it. What exists is 328 ligament surfaces and a one-segment
binding; what would have to be authored is the two-ended attachment, the
stiffness and the slack length, per ligament. Stating that plainly is the
deliverable here, not a number.
