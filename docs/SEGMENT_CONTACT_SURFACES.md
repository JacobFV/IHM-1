# what the body touches the world with, measured

The upright plant's non-foot contact is `fall_proxy_<body>`: one sphere per
segment whose radius is inscribed in that segment's **inertia ellipsoid** and
whose centre is its mass centre. A femur represented as a ball. This file
records what replacing that costs, what the solver does with a concave surface,
and — the part that matters most — **what it takes for the SKIN to be the thing
that meets the floor, which is where contact actually happens.**

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

The two skin rows are on the **worse** of the two registrations below, where the
skin never reaches the floor; they are here for the geometry-cost point and the
working numbers are in *Skin-mediated ground contact, on the better map*.

**761.38 N is m·g to the last digit** (77.6122029 × 9.81), and the
momentum-balance residual `m(a_com − g) − contact − external` stays at 1e-14 in
every arm. That residual is the gate that says the new forces are *summed*
correctly rather than merely printed; a contact term bookkept wrong shows up
there and nowhere else.

**Mesh contact geometry is nearly free until it touches.** `skin_carried` holds
112,000 triangles — 4× the bone bundle, 4000× the sphere count — and costs
0.058 s against the baseline's 0.058 s when nothing is near the plane, 0.079 s
on the registration where the toes do reach it. The OBB broad phase prunes every mesh that is
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

## The skin, cut per segment — and the registration that decides whether it touches

`scripts/build_skin_contact_meshes.py` cuts the canonical exterior skin —
109,183 triangles, 1.7805 m² — into one closed surface per segment:

* partitioned by the repo's own `continuous_surface_binding`, a graph-diffused
  skinning weight per skin vertex per segment; a triangle goes to the argmax of
  its three vertices' mean weight;
* mapped into each segment's own frame through a canonical→source map and that
  map's own reference pose;
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

### The gate: a body's bones are inside its skin

That sentence has an answer everybody knows, so it is the gate. `enclosure()`
ray-parity-tests every segment's own OpenSim bone-mesh vertices against its skin
piece (Möller–Trumbore written out, because trimesh's `contains` needs an rtree
this environment does not have). Controls: points at radius 0.05 inside a
0.1 m icosphere print **1.0**, points at radius 0.5 print **0.0**.

**This repo carries two different canonical↔skeleton registrations and they are
not the same map.**

| | `CanonicalRegistration.global_map` | `binding.json` similarity |
|---|---|---|
| fit | unweighted proper-rigid, 22 approximate COM / bone-envelope-centre pairs | 33 model coordinates **and** one similarity, jointly, on 22 bone-group centroids plus principal axes |
| scale | none | 0.96303 |
| self-reported residual | **123.4 mm RMS**, 352.8 mm max | 24.7 mm RMS on bone-group centroids |
| used by | the supine skin foundation | the anatomy→segment binding |
| **bone vertices inside their own skin** | **0.273** | **0.445** |
| segments enclosing their own bone (≥0.99) | **0 / 20** | **0 / 20** |

Per segment, skin minimum minus bone minimum along the segment's own y — a
positive number means the skin is *above* the bone, which is not a thing a body
can do:

| segment | canonical map | binding map |
|---|---:|---:|
| toes_l | **+94.6 mm** | −9.5 mm |
| calcn_l | **+84.0 mm** | +20.1 mm |
| tibia_l | **+101.3 mm** | +42.3 mm |
| hand_l | **+58.6 mm** | −18.1 mm |
| femur_l | +4.2 mm | −39.9 mm |
| torso | −37.2 mm | −73.1 mm |
| pelvis | −79.5 mm | −117.0 mm |

Under the map the supine foundation uses, the sole of the foot sits **96 mm above
the floor** and 84–95 mm above the foot bone inside the same segment frame: the
toe skin runs +0.083 to +0.164 m while the toe bone runs −0.011 to +0.009 m, so
the skin is entirely above the bone. Under the binding map the foot skin comes
down onto the floor and contact works.

**Neither map passes the gate.** The better one leaves 55% of the skeleton
outside its own skin. That is not a global-fit problem any more — the anatomical
body and the Rajagopal skeleton are different subjects, and one rigid similarity
cannot make a different person's bones fit inside this person's skin. The fix
is per-segment geometric transformation of the anatomy, which is exactly the
"bones and muscles taken from the source to become geometrically parametrized,
transformed entities" that `docs/DIRECTION.md` already asks for.

**Why nothing caught this before.** The supine surface foundation is built
through the worse map, and it sets its support plane to *the skin's own minimum*
— `plane = min source-x of the eligible faces`. The floor follows the skin, so
the error is invisible there by construction. Upright, the floor is at y = 0
because the model says so, and the same error becomes a 96 mm hover.
`docs/SUPINE_SURFACE_SUPPORT_CURRENT.md` already recorded its shadow — "pelvis
39 mm, femurs 88–89 mm, tibias 112 mm above the reference plane" — as a property
of the geometry rather than of the registration.

### The gate was measuring the partition, not only the registration

The paragraph above says the failure "is not a global-fit problem any more"
because the two bodies are different subjects. **That is at most part of it, and
the gate as built cannot tell.** Test the body's OWN anatomical bones against its
OWN skin -- one acquired body, one canonical frame, no registration anywhere --
through the same cut, the same caps and the same `enclosure()`
(`scripts/measure_skin_enclosure_premise.py`):

| segment | own bones inside own skin piece |
|---|---:|
| hand, pelvis | **1.000** |
| torso, toes | 0.97-0.98 |
| tibia, humerus | 0.86-0.92 |
| femur | 0.62-0.64 |
| calcn | 0.04 |
| ulna, radius | **0.01-0.04** |
| patella | **0.000** |
| **mean** | **0.547, 3 of 20 >= 0.99** |

So under this partition even a PERFECT registration tops out at 0.547, and the
binding map's 0.445 is already 81% of it.

Which piece does contain them (`measure_skin_enclosure_cross.py`)? The
controls hold -- hand bones read 1.00 in their own piece and 0 in every other,
and every piece far from a bone reads 0. Two different defects:

* **boundary misassignment.** calcn bones sit **0.78 inside the toes piece**;
  femur reads 0.62 in its own piece and 0.28 in the pelvis piece, where the
  femoral head is.
* **segments that own no closed region at all.** radius and ulna split ONE
  forearm into two strips, and the patella's piece is a patch on the front of the
  knee. A capped strip encloses nothing deeper than itself, and these bones are
  barely inside any piece (row sums 0.23, 0.23, 0.08).

Merging the regions those segments share is the decisive test
(`measure_skin_enclosure_merged.py`), and it passes:

| skin pieces merged | bones | inside |
|---|---|---:|
| radius | radius | 0.008 |
| radius + ulna | radius / ulna | 0.774 / 0.812 |
| radius + ulna + hand | radius | **1.000** |
| patella | patella | 0.000 |
| patella + tibia (+ femur) | patella | 0.917 (**1.000**) |
| calcn + toes | calcn | **1.000** |
| femur + pelvis | femur | 0.900 |

**This body's bones are inside its skin; the per-segment hard partition is what
fails.** The consequences:

1. **The per-segment enclosure gate is structurally unpassable** for radius, ulna
   and patella under ANY registration, so it cannot be the acceptance test for
   one. Registration quality has to be measured against the whole skin, where the
   partition cannot intervene (`measure_skin_enclosure_whole.py`).
2. **A skin contact piece is a SURFACE carried by the segment under it, not a
   volume that encloses that segment's bone.** A forearm is one tube carried by
   two rigid bodies that rotate relative to each other (pronation). Nothing about
   a hard partition can represent that; the real fix is the deformable skin that
   `DIRECTION.md` already requires, with rigid per-segment pieces as the
   scaffold that stands in for it.
3. **"Per-segment geometric transformation of the anatomy" is still needed, but
   it must be judged by the whole-skin gate**, or it will be tuned against a
   ceiling of 0.547 that no transformation can move.

### Skin-mediated ground contact, on the better map

25 advances of 10 ms from the stance pose, skin bundle built on the binding
registration, layer stiffness from the body's own declared skin.

| arm | loaded elements | vertical contact force | pelvis_ty drift | s / advance median | worst |
|---|---:|---:|---:|---:|---:|
| `skin_carried` — skin present, source foot spheres still doing the work | 14 | 761.81 N | +0.1 mm | 0.079 | 0.249 |
| `skin` — source foot spheres removed, the body stands on its skin | 4 | 1019 N, oscillating | −55 mm | 0.528 | 1.396 |

It works, and it is not yet right:

* the plantar skin is **not level** with the floor. The toe skin starts 8.5 mm
  *below* it while the heel is above, so the body rocks forward and 431 N per
  side arrives through the toes against 79 N through the heel. Under the source
  spheres the same pose loads midfoot 80 N, heel 74 N, rearfoot 69 N.
* it overshoots weight by 34% and is still oscillating at 0.25 s.
* it costs 9× the sphere baseline **while in contact**, and 1.4× while merely
  carried.

Both faults are downstream of the same 55%: a skeleton that does not fit inside
its skin cannot put that skin flat on the floor.

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

**Ligaments were a gap. 117 of them are now force elements, and the gap that
remains is a different one.** The full account is `docs/TISSUE_MECHANICS.md`;
what belongs here is the correction to what this file used to say.

* The anatomy carries **300 ligament entities** and **36 joint capsules**, every
  one with `reference_geometry`, bound across all 22 segments.
* The plant now carries **105 ligaments and 12 joint capsules** as
  `Blankevoort1991Ligament` force elements, derived by
  `scripts/build_tissue_force_elements.py`. Standing weight is unchanged at
  761.3757 N and the momentum balance stays at its relative floor, because these
  are internal forces.
* The remaining 195 ligaments are **one-segment**: the two bones they join are
  the same rigid body on this scaffold — 59 in `torso`, 29 per hand, 21 per
  `calcn`. That is the scaffold reporting its own resolution, not a derivation
  failure, and the sacrotuberous, sacrospinous and inguinal ligaments coming out
  one-segment is one of the gates.

This file previously said the two attachment ends "could be constructed ... but
that is a construction, and nothing in the repo would validate it." The
construction was made and it does have gates: 30/30 named bone pairs including
three negative controls, six published lengths at 0.70–1.04x, three published
Blankevoort stiffnesses at 0.44–2.03x, and two independent implementations of the
path length agreeing to 3.3e-16 m.

What that buys is measured on three different prone drops, because one drop is
one drop. **The derived set as a whole makes the plant worse on all three**: the
worst excursion past the model's declared ranges goes 17.4/30.7/30.3 deg bare to
32.8/35.3/32.9 deg with all 105 ligaments. The reason is measured per
coordinate — a real cruciate is near-isometric because it wraps and its femoral
footprint sits near the flexion axis, while a straight line between two
attachment centroids sits 22 mm off that axis, so the derived ACL reads **77%
strain at 90 deg of knee flexion** against a 17.1% ultimate.

**The 66 elements that never pass ultimate strain inside a spanned joint's own
declared range are never worse than the bare plant, and added to the joint stops
they improve every drop**: 5.30/6.62/5.96 deg becomes 4.05/6.39/5.35, mean 5.96
to 5.27. They do **not** replace the stops — on `prone` alone they hold 6.08 deg
with no stops at all, but on `prone_high` the same 66 give 30.47 deg against a
bare 30.69, which is no restraint. The 51 that fail the check are the cruciates,
the collaterals and the ankle ligaments, and what they need is a wrap surface per
joint.

`docs/TISSUE_MECHANICS.md` is the full account.
