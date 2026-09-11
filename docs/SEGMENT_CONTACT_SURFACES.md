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

### The registration, measured without the partition

`scripts/measure_skin_enclosure_whole.py` tests each segment's registered OpenSim
bones against the WHOLE capped exterior skin, where no partition can intervene.
The ceiling is the body's own anatomical bones against the same surface.

| | own bones (ceiling) | canonical map | binding map |
|---|---:|---:|---:|
| **mean** | **0.997** | 0.416 | **0.888** |
| segments >= 0.99 | 20 / 22 | 2 / 22 | 9 / 22 |

The ceiling reads 0.997 (toes 0.966-0.968, everything else 1.000), so the skin
does enclose its own skeleton and the metric is sound on this surface.

**The binding map is far better than the per-segment gate said: 0.888, not
0.445.** Femur, patella, radius, ulna and talus read 1.000 and tibia 0.99. What
genuinely fails is concentrated and specific:

| segment | binding map |
|---|---:|
| hand | 0.57-0.63 |
| calcn | 0.75 |
| toes | 0.78-0.80 |
| humerus | 0.67-0.82 |
| torso | 0.84 |
| pelvis | 0.95 |

That is the extremities and the shoulder girdle -- where a different specimen's
proportions differ most -- and it is where per-segment registration is still
needed. The canonical map reads 0.000 for every foot and hand segment, which is
the 96 mm hover above, and should not be used for anything that touches the
world.

### Per-segment registration does not carry the skin

`scripts/fit_segment_registration.py` fits one similarity per segment (atlas bone
group -> the scaffold's own bone mesh, trimmed symmetric ICP from the global map).
It fits every bone better (radius 4.4 -> 1.7 mm, patella 28.3 -> 2.5, humerus up
to 32.9 -> 4.9). Carrying the SKIN with those maps by linear blend skinning over the
existing skinning weights does not follow:

| | global map | per-segment, blended |
|---|---:|---:|
| **mean enclosure** | **0.888** | **0.873** |
| segments >= 0.99 | 9 / 22 | 10 / 22 |
| hand | 0.63 / 0.57 | 0.81 / 0.81 |
| humerus | 0.82 / 0.67 | 1.000 / 1.000 |
| torso | 0.84 | 0.94 |
| **toes** | 0.78 / 0.80 | **0.16 / 0.10** |
| tibia, radius | 0.99-1.00 | 0.93-0.98 |

(gate: blending with the GLOBAL map for every segment reproduces the global column
to within its print precision.)

**The toes are not a bone-group mismatch.** That was the first explanation offered,
and the names refute it: the atlas `toes` group is 28 phalanges and nothing else,
the five metatarsals sit in `calcn` with the tarsals and sesamoids, and OpenSim
splits the foot the same way (`l_foot.vtp` for calcn, `l_bofoot.vtp` for toes).

**What the rotations actually are.** Split each per-segment rotation relative to
the global map into twist about the bone's own long axis and off-axis swing:

| segment | total | twist | swing | elongation |
|---|---:|---:|---:|---:|
| femur | 19.6-19.8 | 19.5-19.8 | 1.3-1.5 | 8.3-8.7x |
| tibia | 15.0-16.6 | 14.6-16.3 | 3.4 | 7.3-7.4x |
| radius | 31.9-36.1 | 31.6-36.0 | 2.9-4.1 | 10.9-11.9x |
| **toes** | **16.9-17.0** | **1.0-1.6** | **16.8-16.9** | 1.5x |
| talus | 22.9-29.6 | 5.0-8.9 | 22.3-28.3 | 1.4-1.5x |
| patella | 21.1-25.0 | 5.9-6.6 | 20.1-24.4 | 1.1-1.2x |

The long bones rotate almost entirely about their own axis, which point-to-point
ICP cannot determine: those twists are unconstrained, not fitted. The toes are
the one segment whose rotation is almost pure SWING -- a 17 deg tilt of the
forefoot, which lifts the toe skin off the phalanges once it is carried by that
map. That is the likeliest reading of the collapse and it is not yet tested.

So per-segment rigid maps are the wrong instrument for skin. Skin spans joints; a
rigid map per segment, blended, puts seams exactly where the partition did.

The rotation controls settle whether the twists were the problem: they were not.
Blended skin, by how each segment's rotation is fitted:

| rotation | mean enclosure | segments >= 0.99 |
|---|---:|---:|
| global map, one similarity | **0.888** | 9 / 22 |
| per segment, free | 0.873 | 10 / 22 |
| per segment, twist about the long axis removed | 0.857 | 12 / 22 |
| per segment, rotation held at the global map's | 0.830 | 12 / 22 |

Every per-segment variant is worse than the one global similarity on the mean,
including the one with no per-segment rotation at all -- so it is not the
rotations. Per-segment scale and translation, blended across a joint, already
distort the surface. The skin has to be carried by something that is continuous
across joints.

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


## One layer thickness, and the map it should be (2026-09-10)

The skin's elastic foundation derives its stiffness from ONE layer thickness for the whole
body (`skin_layers`, from the three skin-layer entities). A real body has a few millimetres
of soft tissue over the shin and scalp and centimetres over the buttocks, and the programme
requires contact to be mediated by that tissue, never by bone.
`scripts/measure_soft_tissue_depth.py` reads the map out of this body's own geometry: for
40,000 points on its exterior skin, the distance to the nearest of 658 bone and muscle
surfaces.

| region (by nearest deep structure) | median depth | 10-90% |
|---|---:|---:|
| sternum | 6.6 mm | 2.4-11.1 |
| scalp | 7.4 mm | 5.1-9.8 |
| anterior shin | 10.3 mm | 6.6-13.8 |
| thigh | 18.4 mm | 8.5-29.7 |
| buttock | **23.2 mm** | 14.7-27.5 |
| whole exterior skin | 11.0 mm | 5.3-24.6 |

Gate, an anatomical known answer: the shin and the scalp must each be at most half the
buttock. PASS. The depth is to the NEAREST surface, not along the inward normal, so it
underestimates in folds (axilla, groin).

A threefold range: a single thickness makes the foundation too stiff over the buttocks or too
soft over the shin (stiffness goes as 1/h). ~~Next, not done: give each skin patch the local
depth instead of one h.~~ **Built and judged the same day (2428731), below:** each segment
carries its own depth and an in vivo modulus through the V2 engine path. The stance FAILED
never-bone on the toes, and the cause turned out to be the skin's registration on the foot,
not the layer, so the drops wait on the continuous skin carrier.
The per-point map is saved in `data/derived/soft-tissue-depth-v1/depth.npz`.

### Depth alone would make it worse: the modulus is the missing half (pre-registered 2026-09-10, before any run)

Per segment of the skin partition (a vertex belongs to its argmax binding weight), the median
depth is: calcn 18.5/18.6 mm, toes 7.7, tibia 12.1/12.5, patella 18.4/19.2, femur 15.6/16.3,
pelvis 23.4, torso 10.2, humerus 14.9/15.4, radius 9.4/10.2, ulna 11.9/11.6, hand 7.8/7.7.
Known answer: the in vivo unloaded heel pad is 16.0 mm (median, 9.6-17.7; Teng 2022) and
14.85 +/- 2.81 mm (Yang 2022); this body's heel reads 18.5 mm, slightly above, as a whole-
segment median of nearest-structure distance should be (it includes the heel's sides).

Substituting that h into the declared layer (E = 3 kPa) makes the feet SOFTER, not stiffer: k
goes as 1/h, and 761 N on ~0.02 m2 of heel at 18.5 mm would need ~62 mm of compression. That
is bone contact, by arithmetic, so the run is not worth making. The declared 3 kPa is a
source_informed_prior shared by epidermis, dermis and hypodermis, not a loaded-tissue value.

Loaded soft tissue is an order of magnitude stiffer, and it has been measured in vivo
(`data/sources/in-vivo-soft-tissue-compression.json`):

| site | apparent modulus | source |
|---|---:|---|
| heel pad, gait | 192.6 kPa (median; 130-266) | Teng et al. 2022, doi:10.1186/s12891-022-05197-w |
| heel pad, gait, non-diabetic | 265.5 kPa (median; 155-306) | Yang et al. 2022, doi:10.3389/fendo.2022.894383 |
| heel pad | up to 175 kPa | Gefen et al. 2001, as cited by Teng |
| buttock fat, sitting | 39 kPa secant (18 kPa at 46%) | Linder-Ganz et al. 2007, doi:10.1016/j.jbiomech.2006.06.020 |

The heel values are pressure over THICKNESS strain. Yang's own pair (144.8 kPa at 0.523
strain) gives 277 kPa against their reported 265.5, so they are layer moduli and enter as
`k = E/h` directly; the confined-layer factor `(1-v)/((1+v)(1-2v))` (3.8 at v = 0.45) would count
the confinement twice.

**The rule, fixed now:**

* `h` per segment = that segment's measured median depth (above).
* `E` for calcn = the median of the three in vivo heel values, **192.6 kPa**.
* `E` for pelvis = the Linder-Ganz secant, **39 kPa**, labelled transferred (it is a peak
  principal value from an FE model, weaker evidence than the heel).
* `E` for every other segment = the same 39 kPa, labelled **unsourced for that site**. No in vivo
  value was retrieved for the forefoot, knee, hand, forearm, elbow or trunk; the toes in
  particular carry plantar pad, which is probably stiffer than fat, and this rule makes them soft.
* `k = E/h`. Damping and friction unchanged. The foundation stays linear: both heel fits are
  linear-elastic plus viscous, so it is linear inside the range they were fitted over.

Standing arithmetic, as the expectation: 38 kPa on the heel at 192.6 kPa is ~20% strain, 3.2 mm
of an 18.5 mm layer.

**Gates, fixed before the engine carries per-segment stiffness:**

1. **Momentum balance** (the plant's arithmetic identity): the maximum residual over the stance
   run <= 1e-5 N, the bound `verify_native_fall_contact.py` uses. If this fails, the per-row
   stiffness was bookkept wrong and nothing else is read.
2. **Never bone** (the programme's requirement): on every contacting segment, the maximum
   compression over the run (depth of its lowest skin vertex below the floor) is less than
   that segment's `h`. One segment bottoming out is a FAIL.
3. **Heel strain inside the in vivo range**: the calcn maximum compression over `h` <= 0.73, the
   top of Teng's gait range. Standing should sit well under it.

Reported, not gated: vertical contact force against weight (the pose is not an equilibrium), the
declared-layer `skin` arm beside it (-55 mm, oscillating, in the stance table above), cost. The
three prone drops follow the stance with gates 1 and 2 unchanged.

### Result: FAIL, on the toes (2026-09-10)

The engine now carries a stiffness per contact mesh (`IHM_SEGMENT_CONTACT_MESHES_V2`: the
stiffness is the last field of every row; V1 bundles read as before), the layer map is
`scripts/apply_soft_tissue_layer_map.py` into `data/derived/segment-contact-meshes/skin-layer-map-v1`,
and the gates are computed by `scripts/score_skin_layer_map.py`. Stance, 50 x 10 ms, source foot
spheres removed, both arms in one run:

| | declared layer (E 3 kPa, h 6.6 mm) | per-segment layer map |
|---|---:|---:|
| 1 momentum, worst residual | 7.6e-13 N | **6.8e-13 N, PASS** |
| toes worst compression (L/R) | 22.1 / 22.0 mm = 3.35 h | 12.0 / 12.0 mm = **1.56 h, FAIL** |
| heel worst compression (L/R) | 14.8 / 15.0 mm = 2.24 h | 6.0 / 6.3 mm = 0.32 / 0.34 h |
| 2 never bone | FAIL | **FAIL** |
| 3 heel strain <= 0.73 | FAIL | **PASS** |
| final vertical contact vs weight 761.4 N | 769.9 N | 777.4 N |

The verdict is FAIL, by the gate as fixed. The heel is where the in vivo modulus was sourced and
it behaves like a heel: a third of its layer at standing, against more than twice under the
declared skin. The toes are the segment the pre-registration said the rule would make soft.

What the failure is made of, measured after the verdict and reported, not used to rescore:

* **The pose starts inside the floor.** The lowest toes skin vertex is 8.5 mm below the floor at
  t = 0, before any load: 110% of the toes' 7.7 mm layer. The engineering stance pose was solved
  for the source foot spheres, not for this skin, so the never-bone gate on the toes fails at the
  first sample whatever the stiffness.
* **The load sits on the toes.** At the end the toes carry 353 N each and the heels 29-43 N: 93%
  of body weight on the forefoot. Quiet standing puts roughly half on the heel. The same pose
  fact: the toes skin hangs lower than the heel skin in this registration.
* **The plantar patch is not the segment median.** Taking only the lowest 5-10 mm of each
  segment's skin, the heel reads 14.8-16.5 mm -- the in vivo pad is 16.0 (Teng) and
  14.85 +/- 2.81 mm (Yang), so the plantar band reproduces the known answer better than the 18.5 mm
  whole-segment median does. The toes' plantar band is 6.6-6.8 mm over the phalanges and flexor
  digitorum brevis, slightly under their 7.7 mm median. Heel strain at the plantar h would be
  0.42, still inside gate 3.

~~Next: a stance pose solved against the skin, and an in vivo forefoot modulus.~~ **Withdrawn
the same day, before anything was run on it: the pose is not the cause, the foot registration
is, and a pose solved against this skin would put bone through it.**

A probe at t = 0 (one stream, 17.9 s to start, pose applied exactly): the heel skin's lowest
vertex is **35.7 mm ABOVE** the floor while the toes skin is 8.0-8.5 mm below it. The foot is
tipped onto its toes before anything moves; the heel compression in the table came later, as
the body rotated down. And in each segment's own frame the heel skin's lowest point sits
**20 mm above** the scaffold calcaneus's lowest point (8% of that bone's vertices inside its
skin piece), where the canonical anatomy puts the plantar skin ~15 mm BELOW its own calcaneus
(the plantar-band depth above). Levelling the skin on the floor would therefore drive the
scaffold's calcaneus ~20 mm through the heel skin -- the one contact the programme forbids.

This is the registration defect already recorded above (whole-skin enclosure calcn 0.75, toes
0.78-0.80), now with its size. Against the scaffold's own foot bones, the global map leaves the
atlas calcaneus group at 15.0 mm RMS (toes 8.0); fitting them alone takes a similarity with
**scale 1.22** and a **20 deg** rotation from the global map (toes 1.21 and 17 deg). The
scaffold's foot is about a fifth larger than this specimen's, relative to the rest of the body,
which one global similarity cannot absorb; and carrying the skin with per-segment rigid maps
is already measured to make it worse (toes 0.10-0.16).

So skin-mediated stance is blocked on the item this file already names -- a skin carried by
something continuous across joints -- not on the contact layer. The layer map (per-segment
depth, in vivo modulus, V2 engine path) is built, gated and ready for it; its heel behaves like
a heel. The toe pulp has no in vivo modulus in the literature retrieved
(`data/sources/in-vivo-soft-tissue-compression.json`, `not_found`).

### A skin carrier continuous across joints: gates fixed before it is built (2026-09-10)

The instrument: one smooth space warp from atlas space to the scaffold's ground, fitted to
correspondences between the atlas bone groups and the scaffold's own bone meshes (the
per-segment fits above supply them), then applied to the WHOLE canonical skin before it is cut
per segment. A single warp has no seams, so it cannot do what blended per-segment maps did at
the joints. It can fold, which gate 4 is for.

1. **Bones (known answer):** on every segment, the warped atlas bone group lies no further from
   the scaffold's bone mesh (RMS nearest-surface) than the per-segment similarity leaves it, plus
   1 mm. A warp with more freedom than a similarity that fits the bones worse is fitted wrong.
2. **Enclosure (the goal):** whole-skin enclosure (`measure_skin_enclosure_whole.py`, the
   partition-free measure) mean **>= 0.95**, and calcn and toes each **>= 0.95** -- against
   0.888 / 0.75 / 0.78-0.80 on the global map and a ceiling of 0.997 (toes 0.966-0.968).
3. **The heel sits on its pad:** in the calcn frame, heel skin's lowest point minus the scaffold
   calcaneus's lowest point in **[-25, -5] mm** (in vivo pad 9.6-17.7 mm, Teng 2022; this
   specimen's plantar band 14.8-16.5 mm). Today it is +20 mm.
4. **No folding:** the warp's Jacobian determinant > 0 at every skin vertex, and no skin triangle
   inverts. One fold is a FAIL.

Reported, not gated: per-segment skin area change, the warp's bending energy, and the stance
arm (layer map, this skin) scored by `score_skin_layer_map.py` with its three gates unchanged.

#### The warp as fitted, before gates 2-4 were computed (2026-09-10)

`scripts/fit_skin_warp.py --stage fit`; the warp itself is `scripts/skin_warp.py`. `W(x) = G x +
d(G x)`: `G` is the binding similarity -- the exact map the 0.888 column is measured through --
and `d` a regularised 3D thin-plate spline (kernel `-r`). Correspondences: 200 area-weighted
samples per atlas bone group (torso 600, pelvis 300); source `G a`, target the exact nearest point
ON the scaffold's bone mesh to the per-segment similarity's image `M_seg a`; the farthest 10% per
segment dropped, the per-segment fit's own trim. 4,410 kept. They ask for 7-45 mm of displacement
from the global map (median by segment; up to 80 mm at humerus_r).

**Regularisation rule, fixed in the script before the fit ran and computed from the bone
correspondences alone:** 5-fold cross-validation over `lambda in {0} U logspace(-8, 0, 17)`; take the
LARGEST lambda whose held-out RMS is within 1% of the minimum. Minimum 4.570 mm (flat up to 1e-4);
chosen **lambda = 1e-3**, 4.588 mm, not at the grid's edge. Fitted residual at the correspondences
0.571 mm RMS; bending energy 8.69.

Known-answer controls, all before any gate: `G` in `registration.json` equals the bundle's binding
map bit for bit; zero displacement reproduces the binding map on all 102,467 skin vertices bit for
bit; gate 1's instrument, replaying the per-segment fit's own evaluation samples, returns its
stored RMS exactly (0.0 m difference, 22 segments); the nearest-point projector returns 2.5e-16 m
for points on the mesh (a first, k-nearest-candidate version returned 0.13 mm and was replaced by
an exact bounded search); the spline interpolates at lambda 0 (2e-15 m) and reproduces an affine
field with zero bending; its Jacobian matches central differences to 1.8e-9. **Zero warp through
the two gate instruments:** whole-skin enclosure reads **0.888** per segment identical to the
binding column (ceiling **0.997**), and the bundle built through a zero warp is identical to
`skin-binding` in every record -- calcn **+20.1 / +20.0 mm**, same mesh sha256s.

**Gate 1: PASS, 22 / 22.** The warped bone group is at or under the per-segment similarity on 21
segments (calcn 6.09 -> 5.37 / 6.00 -> 5.10 mm, torso 15.65 -> 11.28, humerus 5.38 -> 3.77); femur_l
is 0.04 mm over it (3.92 vs 3.88), inside the 1 mm margin.

#### Result: FAIL -- the heel reaches its pad, but the forefoot stays outside and the toe skin folds (2026-09-10)

The warp committed in 96e5f1b, unchanged, through the pre-registered instruments
(`fit_skin_warp.py --stage score`; gate 2 is `measure_skin_enclosure_whole.py --warp`, gate 3 the
calcn records of `build_skin_contact_meshes.py --warp` into
`data/derived/segment-contact-meshes/skin-warp-v1`):

| gate | value | threshold | |
|---|---|---|---|
| zero-warp control | enclosure 0.888 identical per segment, ceiling 0.997; calcn +20.1 / +20.0 mm; bundle identical to `skin-binding` | reproduce to print precision | **PASS** |
| 1 bones | 22 / 22; closest to its limit femur_l 3.92 mm against 4.88 | warped <= per-segment + 1 mm | **PASS** |
| 2 enclosure | mean **0.954**; calcn **0.926 / 0.920**, toes **0.868 / 0.880** | mean >= 0.95; calcn and toes each >= 0.95 | **FAIL** |
| 3 heel on its pad | **-6.7 / -6.6 mm** | [-25, -5] mm | **PASS** |
| 4 no folding | det J <= 0 at **45** of 54,949 skin vertices (min -0.68; global map 1.12); **72** of 109,183 skin triangles inverted | none | **FAIL** |

**The verdict is FAIL, on gates 2 and 4.** The stance was not run: the pre-registration makes it
conditional on all four gates, so no layer-map bundle was built from this skin and
`measure_segment_contact_meshes.py` is unchanged.

Whole-skin enclosure per segment, the segments that moved:

| segment | binding map | warped |
|---|---:|---:|
| calcn | 0.754 / 0.756 | **0.926 / 0.920** |
| toes (ceiling 0.966 / 0.968) | 0.784 / 0.796 | **0.868 / 0.880** |
| hand | 0.634 / 0.572 | 0.794 / 0.810 |
| humerus | 0.821 / 0.674 | 0.979 / 0.988 |
| torso | 0.844 | 0.922 |
| pelvis | 0.948 | **0.912** (worse) |
| talus, tibia | 0.987-0.994 | 0.994-1.000 |
| **mean** / >= 0.99 | 0.888 / 9 of 22 | **0.954** / 12 of 22 |

Reported, not gated: bending energy 8.69; skin area 1.9198 -> 1.9114 m2 overall, but calcn
x1.30, talus x1.22, toes x1.09-1.10, tibia x1.06-1.07, femur x0.89-0.90, hand x0.85-0.86 -- the warp
carries the scaffold's larger foot and smaller hand into the skin, as it was built to. Caps: 0 of
132 inverted.

What the failures are made of, measured after the verdict and not used to refit (`--stage
diagnose`, which reproduces `enclosure()`'s value on every segment before reading its points):

* **Every fold is in the toe skin**: 20 + 25 vertices and 33 + 39 triangles, all in the toes
  partition, 7.5-15.3 mm from the nearest spline centre. The toe skin is thin and sits between
  phalanx correspondences (toes map: scale 1.21, a 17 deg swing) and metatarsal ones a few mm
  behind them (calcn map: scale 1.22, 20 deg); where those two ask for different displacements
  over less than the skin's own thickness, the spline turns the skin over. That is the joint
  conflict the pre-registration named, at the MTP; it is the likeliest reading and it is not
  separately tested.
* **The foot bone still outside is the far forefoot, not the heel.** Under the binding map calcn's
  outside points spanned the whole bone (x -11 to +212 mm in its frame); under the warp they are
  the distal lateral end (median x +199 of 212 mm, z +-40 mm: the lateral metatarsal heads). The toes'
  are the phalanx tips (median x +73 on a bone running -41 to +87 mm). The warp brought the heel
  down onto its pad (gate 3) but the skin does not reach the ends of a forefoot a fifth longer
  than this specimen's: skin with no correspondence under it is extrapolated, not carried.
* **Pelvis got worse** (0.948 -> 0.912; 39 of 442 points outside, median 10.7 mm, against 23 at
  6.5 mm): the points outside are the most posterior bone (x -207 to -153 mm in the pelvis frame),
  the sacrum and posterior ilium, under the thinnest skin of the pelvis.

So one smooth warp on bone correspondences alone fixes what it has correspondences under -- heel
height, humerus, torso, most of the hand -- and fails at the two places where skin extends beyond
the bones it is fitted to (the toe tips, the forefoot's lateral edge) and where the per-segment
maps it interpolates disagree across a joint thinner than the skin over it (the MTP).

#### Next attempt: a warp that cannot fold. ONE change, and the gates do not move (2026-09-10)

A thin-plate spline is free to turn a surface over, and gate 4 says it did. The instrument is
therefore replaced by one that cannot: `W(x) = G x + phi(G x)`, where `phi` is the flow of a
STATIONARY VELOCITY FIELD integrated by scaling-and-squaring. A smooth velocity field's flow is a
diffeomorphism, so det J > 0 holds by construction rather than by luck, and the check becomes
numerical (enough squaring steps that the per-step displacement is small) instead of a gate the
fit can fail.

**Exactly one thing changes.** The correspondences are the same 4,410, built the same way; the
regularisation is chosen by the same rule (5-fold CV on the bone correspondences alone, largest
value within 1% of the minimum held-out RMS); the zero-warp control and gates 1-4 and their
thresholds are unchanged -- including calcn and toes at 0.95, which this attempt missed at 0.926
and 0.868-0.880. Changing the instrument after a FAIL is allowed; moving the line it failed
against is not, and it has not moved. Additionally required, because a flow is not a spline: the
velocity field's integration must be verified against its own known answer -- a constant velocity
field flows to a pure translation, and the composed forward and inverse flows must return every
skin vertex to itself within 1e-9 m.

**Predicted here, before it runs, because it separates the two causes the diagnosis names:**
folds go to ZERO (the instrument forbids them), and enclosure at calcn and toes is expected to
IMPROVE BUT STILL FAIL, because those failures are skin extrapolated beyond the bones it is
fitted to, which a different smoothness class does not supply. If the folds go and toes still
read ~0.87, the remaining problem is correspondence coverage, not the warp family, and the next
step is anchors for skin that has no bone under it. If the toes instead reach 0.95, the fold and
the coverage story were one thing, and this note was wrong about it.

#### Result: FAIL on all four. The fold-free family is not a trade, it is worse (2026-09-10)

`fit_skin_warp.py --stage fit --family flow`, then the same three instruments. The velocity field
is the same RBF family on the same 4,410 correspondences (identical per-segment kept counts and
medians to the spline's, so only the warp family changed), integrated to its time-1 flow.

**The flow's own known answers, both required by the pre-registration:**

* a **constant** velocity field flows to a pure translation: **1.7e-15 m** (threshold 1e-12, the
  double-precision accumulation bound for 2^7 additions over metre-scale coordinates, not a
  tolerance picked after seeing the number); its own inverse returns 2.8e-17 m;
* **forward then inverse returns every one of the 102,467 skin vertices to itself: 1.31e-14 m**
  (threshold 1e-9).
* Integration: **128 steps, 2^7 squarings**, max |v| 96.7 mm, **max per-step displacement 0.755 mm**
  against a 1.0 mm budget; smallest single-step determinant 0.9727.

The squaring is done by composing the analytic half-step map, not by interpolating a stored field
on a grid: the field is parametric, so squaring N times is exactly 2^N applications of it. That is
also why the inverse is exact to 1e-14 instead of to a grid's interpolation error.

| gate | spline (v1) | flow (v2) | threshold |
|---|---|---|---|
| zero-warp control | 0.888 / ceiling 0.997 / calcn +20.1, +20.0 mm | identical | reproduce exactly |
| 1 bones | 22 / 22 **PASS** | **20 / 22 FAIL** -- patella_l 4.02 mm vs limit 3.47, patella_r 4.10 vs 3.52 | <= per-segment + 1 mm |
| 2 enclosure | mean 0.954; calcn 0.926 / 0.920, toes 0.868 / 0.880 **FAIL** | mean **0.959**; calcn **0.922 / 0.914**, toes **0.872 / 0.878** **FAIL** | mean >= 0.95; calcn, toes each >= 0.95 |
| 3 heel | -6.7 / -6.6 mm **PASS** | **-5.7 / -4.1 mm FAIL** (calcn_r 0.9 mm short of the -5 mm edge) | both in [-25, -5] mm |
| 4 no folding | 45 vertices det J <= 0, 72 triangles **FAIL** | **0 of 54,949 vertices** (min det 0.108); **2 of 109,183 triangles FAIL** | none |
| | **2 of 4** | **0 of 4** | |

**The verdict is FAIL on all four, and the two instruments are not a trade-off.** The spline passes
gates 1 and 3; the flow passes none. Outside vertex-level det J the flow is worse or equal
everywhere the gates look.

**Gate 4 is the informative failure.** The flow delivered exactly what a diffeomorphism promises:
not one negative Jacobian determinant anywhere on the skin, against 45 for the spline, and every
single integration step orientation-preserving (worst 0.973). **Two finite triangles still
inverted**, because det J > 0 is local invertibility AT A POINT and a triangle has size. Both are
slivers in the source mesh: aspect ratio (longest edge / 2x inradius; 1 is equilateral) **12.4 on
toes_r and 16.6 on tibia_r, ranks 430 and 123 worst of 109,183, against a mesh median of 2.42**. A
sliver's normal is the cross product of two nearly parallel edges, so it is ill-conditioned by
construction, and these two span 15 and 52 mm -- far enough for the flow's rotation to vary across
them. The turns are not marginal: **156 deg and 147 deg**, and the tibia_r triangle also collapses
to 0.20 of its area. So this is two degenerate triangles in the canonical mesh, not a folding warp
(`--stage slivers`).

**The prediction in 3b1526c was half wrong, twice.** It said folds go to zero and calcn and toes
improve but still fail. Folds went to zero at the VERTICES and not at the triangles. And calcn and
toes did not improve: 0.922 / 0.914 and 0.872 / 0.878 against the spline's 0.926 / 0.920 and
0.868 / 0.880 -- calcn slightly worse, toes a wash. The mean moved (0.954 -> 0.959, 12 -> 14
segments at >= 0.99, on hands 0.79-0.81 -> 0.81-0.84, humerus to 1.000, torso 0.922 -> 0.940) and
the two segments the gate actually names did not.

**What the fold-free family cost, mechanically.** The same CV rule chose **lambda = 1e-1**, two
orders of magnitude smoother than the spline's 1e-3, because the flow's held-out error is
minimised there: 6.585 mm flat from 0 to 1e-4, then 6.432 (1e-3), 5.843 (1e-2), **5.522 (1e-1,
chosen)**, 5.720 (3.16e-1), 6.335 (1). The spline's minimum was 4.570 mm. The fit is
correspondingly loose -- residual at the correspondences **4.043 mm RMS, 34.11 mm max**, against the
spline's 0.571 and 8.45 -- and the greedy correction of the velocity field **did not converge**:
worst |flow - target| per pass 67.11, 26.94, 31.29, 30.25, 28.37, 28.51 mm. That looseness is
exactly what gates 1 and 3 then caught: the patellae need the largest local change of any segment
(per-segment scale 1.35) and come out 0.55-0.58 mm past their margin, and the heel is pushed 2.5 mm
less far down than the spline pushed it, which leaves calcn_r 0.9 mm outside the pad window.

**No variant in the sweep was scored on enclosure**, and none can be quoted as reaching 0.95: the
CV curve above is held-out error on the bone correspondences only, which is what the rule allows it
to see, and only the rule's own lambda was carried to gate 2. Running the other lambdas through the
enclosure gate would be choosing the instrument on the gate it is judged by.

**Both families now fail gate 2 at the same two segments and nearly the same numbers**, with
smoothness classes as different as a spline and a diffeomorphic flow. That is evidence about the
correspondences, not the warp: the spline's diagnosis said the toe tips and the lateral forefoot
are skin with no bone under them, extrapolated rather than carried, and a different smoothness
class has now confirmed it does not supply what is missing. The visible fix is anchors for skin
that has no correspondence beneath it -- and that is a new instrument, so it is a new
pre-registration and not a refit of this one.

#### Third attempt: anchors where the skin has no bone under it. Fixed before it is built (2026-09-10)

Two smoothness classes, one that fits the bones to 0.571 mm and one that cannot fold, fail gate 2
at the same two segments within 0.006. The warp family is not what is missing; correspondences
are. Every correspondence today is a bone sample, so skin over the toe tips, the lateral forefoot
and the hands is extrapolated from bone that is elsewhere.

**The instrument.** The SPLINE of 96e5f1b, unchanged -- it passed gates 1 and 3, and the flow is
set aside as strictly worse on this evidence -- with its correspondence set extended by ANCHORS:

* an anchor is added for every exterior skin vertex whose distance in atlas space to the nearest
  existing bone correspondence exceeds **20 mm**, fixed now and not tuned: this body's median
  skin-to-bone/muscle depth is 11.0 mm, so a skin vertex 20 mm from every bone sample has no bone
  under it in any sense the warp can use;
* its target is that vertex's image under **its own segment's per-segment similarity** -- the map
  that already fits that segment's bones best, and the only local statement available where no
  bone is beneath the skin;
* anchors carry the same weight as bone correspondences, stated rather than fitted;
* `lambda` is chosen by the SAME rule (5-fold CV, largest within 1% of the minimum), computed on
  the BONE correspondences only. Anchors are a modelling assumption, not measurements, and a rule
  that cross-validates on them would be scoring the assumption against itself.

**Known answers, before any gate is read:** with the anchor set empty the pipeline must reproduce
the spline warp of 96e5f1b on all 102,467 skin vertices to floating point, and the zero-warp
control must still return the binding map's 0.888 and calcn +20.1 / +20.0 mm.

**The gates do not move.** 1-4 exactly as in 93d3d57, including gate 4, which the spline failed with
72 inverted triangles -- of which the flow has now shown at least two are slivers in the canonical
mesh (aspect 12.4 and 16.6) that a fold-free map still inverts. Sliver triangles are reported
alongside the count; they are not excused from it.

**Predicted here, before it runs:** calcn and toes rise above 0.95 and gate 2 passes; gate 1 holds,
because anchors add constraints where there were none rather than moving the bone targets; and
folds fall well below 72 without necessarily reaching zero, since the two slivers are a property of
the mesh. If calcn and toes stay near 0.92 and 0.87 with anchors under them, extrapolation is NOT
what is wrong, and what remains is the hard partition itself -- the seam at the MTP where the toes'
and calcn's maps disagree across skin thinner than the joint.

**A flaw in that rule, recorded BEFORE any gate of it was computed.** The fit reports **34,411
anchors of 54,949 exterior skin vertices -- 62.6%**, with a median skin-to-nearest-correspondence
distance of 27.3 mm. That is not what the rule was for. The correspondences are sparse SAMPLES
(4,410 over the whole skeleton, 200 per segment), so a skin vertex can be 27 mm from the nearest
sample with bone squarely beneath it; I wrote "far from a bone sample" and meant "no bone under
it". At 62.6% the warp is largely driven by per-segment similarities, which is close to the blended
per-segment maps this file already measured as WORSE than the global map (0.873 against 0.888).

The run is not stopped and its verdict stands as the rule's own result: the rule was fixed, and
amending it because an intermediate count looked wrong would be exactly the move this structure
forbids -- the more so since no gate of it has been read. **A second instrument is pre-registered
here instead, to be built after that verdict is recorded and reported beside it:** identical in
every respect except that a vertex is anchored when it is more than 20 mm from the nearest point on
any atlas BONE SURFACE, rather than from the nearest sampled correspondence. This body's median
skin-to-bone/muscle depth is 11.0 mm, so 20 mm from bone surface is genuinely skin with no bone
under it. Gates 1-4 unchanged; the anchor count is reported for both. **Predicted:** the corrected
measure anchors a small fraction -- the toe tips, the lateral forefoot, the hands and the fleshy
trunk -- rather than two thirds of the body.

**Both anchor rules are wrong, and the second one refutes the idea rather than the threshold
(2026-09-10).** The corrected measure anchors **27,209 vertices, 49.5%** -- not a small fraction.
The reason is arithmetic: this body's median skin-to-bone-SURFACE distance is **19.7 mm**, so a
20 mm threshold selects about half the skin by construction. The threshold was doing the selecting
in both versions; only the direction of the error changed.

And it anchors the wrong places. Hands fall from 1,360 anchors to **47**, toes from 243 to 142,
while pelvis keeps 6,608 and torso 6,249. **Gate 2 fails at calcn, toes and hands** -- exactly the
segments the corrected rule strips anchors from, because hand and toe skin sits CLOSE to bone.
Anchoring "skin with no bone under it" was never a description of the failure: those segments fail
because the scaffold's foot is a fifth larger than this specimen's, not because their skin is
unsupported. No anchor rule of this kind addresses that, and none is pre-registered.

**What the failure probably is, now that the missing control has a first number.** The recovery
control of `d626af3` reports that for a known field displacing this body's bones **1.50 mm** on
average, the nearest-point targets miss the truth by **mean 1.318 mm, p90 2.372, max 4.052** --
roughly twice the chest-wall line's 0.671 mm on the equivalent test, and **of the same order as the
displacement being fitted**. A warp fitted to targets that wrong is fitting noise the size of its
signal, which is a sufficient explanation for calcn and toes sitting at 0.92 and 0.87 under two
unrelated smoothness classes.

**So the next instrument here is a correspondence, as it was for the chest wall.** Fixed before it
is built: replace nearest-point targets with `ihm/anatomy/normal_shooting.py` -- shooting along the
source surface's own normal, with the normal-agreement filter ON and declared, which the chest-wall
line adopted on the argument that a hit whose surface faces away is not that point's partner. Then
**re-run the recovery control first**, and only judge a warp family afterwards. **Predicted:** the
target error falls from 1.318 mm to a few tenths, as it did there (0.671 to 0.078); calcn and toes
improve but do not reach 0.95, because a fifth of a foot is a shape difference and not a
correspondence error. The thin-sheet failure applies squarely -- these correspondences include ribs,
scapulae and the sternum -- so the agreement filter is not optional here.

**Caught before acting on it: 1.318 mm is measured at the WRONG SCALE for this line, and I was
about to kill the v3 fit over it.** The control above displaces bones 1.50 mm because that is the
chest wall's scale, where 1.318 mm of target error is 88% of the signal and damning. This line's
skin moves about **20 mm**. At 20 mm the same absolute error would be 6.6%, which is not obviously
disqualifying. The programme already learned this twice -- *score a correspondence at the
separation it will be used at*, and *a tolerance calibrated at zero separation is a different test
at 8 mm* -- and I still read a 1.5 mm-scale number as a verdict on a 20 mm-scale instrument.

So **nothing is killed and nothing is concluded until the 20 mm-scale recovery number lands**
(running now, `data/derived/skin-warp-v1/recovery.log`). Fixed before it does, so this is a
prediction and not a reading:

- if target error at 20 mm stays near 1.3 mm absolute, the bias is a FIXED OFFSET, it is 6.6% of
  this line's signal, and the anchored-fit comparison in v3 is readable as a warp-family
  comparison -- while the anchor rule above stays misdescribed and gate 1's 1 mm margin stays
  inside the bias, because THAT margin is a millimetre quantity whatever the displacement is;
- if it scales with the displacement -- roughly 18 mm of error on a 20 mm field -- the targets
  carry no signal at this scale at all, the v3 fit is fitting noise, and it is stopped on the spot.

**Predicted: it stays near 1.3 mm absolute.** The d^2/R argument makes nearest-point bias a
function of surface CURVATURE and separation, not of how far the truth happens to move, so the
error should be governed by the skin-to-bone gap (median 19.7 mm here) rather than by the field
amplitude. If it instead scales with the field, that argument is wrong and the chest-wall line's
correspondence work rests on it too.

The v3 anchored fit is meanwhile at **7,401 s with no progress output** -- the iteration logging
the fit gained applies only to the next run, and the 2 h budget with it. It is left alone rather
than restarted for logging, since restarting costs more than the missing lines are worth; its
successor carries both.

**A selection effect caught in the next instrument before it ran.**
`normal_shooting.shoot_pairs` has `cap_m` defaulting to **20 mm** -- a shot that finds nothing
inside 20 mm is dropped as `no_hit`. That default was set for the chest wall, whose separations
are millimetres. On this line the per-segment RMS reaches 15 mm and correspondence displacements
reach **95 mm**, so the default would have silently discarded exactly the long correspondences
that carry the failing segments, and the survivors would have been the easy ones. The result would
have looked like an improvement and would have been a filtered population -- the same shape as
every row in this programme's ledger: a quantity computed correctly, compared against the wrong
population.

So three settings are declared in the artefact rather than defaulted, fixed here before the
instrument is built: `cap_m` raised past the largest real displacement, `min_normal_agreement`
set explicitly (it defaults to OFF and this line needs it ON), and the shot run **from the
per-segment similarity's image to the scaffold surface, both in ground space** -- not from atlas
space, where the two surfaces are not yet in correspondence at all. `shoot_pairs` returns its drop
reasons (`no_hit`, `no_return`, `return_too_far`, `normal_disagreed`), and the kept fraction with
its reason breakdown is reported beside every number this instrument produces. A correspondence
result without its drop census is not readable.

## The stop rule fired, and I am not honouring it. Why that is not gate-loosening.

| known field moves bones | target error mean | p90 | max | as % of signal |
|---|---:|---:|---:|---:|
| 1.50 mm (chest-wall scale) | 1.318 mm | 2.372 | 4.052 | 88% |
| **20.00 mm (called "this line's scale")** | **16.816 mm** | 31.016 | 43.783 | **84%** |

**Predicted: near 1.3 mm absolute. Got: 16.816 mm.** The fraction is near-constant, so the error
is PROPORTIONAL to the displacement. `0355f77` fixed the consequence in advance -- "the targets
carry no signal at this scale, the v3 fit is fitting noise, and it is stopped on the spot."

**The v3 fit is not being stopped, because the rule's premise is false as a matter of measured
fact.** I wrote "this line's skin moves about 20 mm" and treated that as the separation the
correspondence operates at. It is not. In `fit_skin_warp.py` the target is the nearest scaffold
point to **M_seg.a**, the per-segment similarity's image, which has already removed most of the
displacement. The separations the pipeline actually presents are the `|t - M_seg a|` medians:
**0.78 mm (radius) to 6.88 mm (torso)**. The 20 mm run measures a separation this line never
operates at, so the rule fired on a quantity outside the pipeline's range.

**This is the same error a third time, inside the pre-registration written to avoid it.** Score a
correspondence at the separation it will be used at. I conflated DISPLACEMENT with SEPARATION: the
skin does move ~20 mm, and the correspondence still runs at 1-7 mm, because a similarity transform
sits between them. Twice today I caught this in someone else's instrument and then wrote it into my
own gate.

**So the rule is retired, not quietly dropped, and its replacement is post-hoc and labelled.**
Fixed before the replacement number is read: re-run the recovery control with a starting map that
leaves **1-7 mm** of separation, the range measured above. If target error there is a large
fraction of the residual displacement, the line stops as the retired rule intended. The retired
rule's threshold is not reused, reweighted, or applied to the new instrument -- a threshold set for
20 mm has no meaning at 1-7 mm, which is the whole content of the error above.

**Two further findings in the agent's own control, one of them a flaw it reported against itself.**

1. **The dominant term is tangential blindness, and it scales linearly BY CONSTRUCTION.** The known
   field is a random spline, so it slides each bone surface along itself as well as normal to it,
   and the tangential part is invisible to *any* surface-based correspondence. The split shows at
   the chest-wall scale: pointwise recovery **1.510 mm RMS** against to-surface recovery
   **0.510 mm**. The identifiable component is recovered three times better than the pointwise
   number says.
2. **Therefore the d^2/R argument is NOT refuted here, and the alarm does not propagate.**
   `0355f77` said that if the error scaled, "that argument is wrong and the chest-wall line rests
   on it too." That inference does not hold: the control contains a linear term by construction, so
   it cannot separate a d^2/R curvature bias from tangential blindness. The control is silent on
   d^2/R rather than against it, and **nothing is withdrawn from the chest-wall line.** Retracting
   that inference matters as much as the retired rule -- it would have been a withdrawal made on a
   measurement that does not address the claim.

**And it caps what the next instrument can deliver.** Normal shooting removes nearest-point bias;
it does not see tangential motion either. So the predicted improvement of `95c1e94` (0.671 ->
0.078 there, "a few tenths" here) applies to the **to-surface** measure and NOT to the pointwise
one, where an identifiability floor sits underneath any surface method. Both are reported
separately, and the pointwise number is never quoted as a correspondence quality.

**v3's disposition:** left running -- it is nine hours deep, killing is irreversible, and the
anchored family is already closed as an idea (`95c1e94`), so it can only return a negative, which
is worth having and costs nothing further. **Everything downstream stays frozen**: no gates 2-4 on
v3, no normal-shooting instrument, and no verdict read from v3 until the 1-7 mm control lands. The
agent froze these on its own initiative when my two instructions conflicted, and that was right.

### The invariant 85-90% has a KNOWN ANSWER, and it is pi/4

The fraction held across every scale tested -- 88% at 1.5 mm, 84% at 20 mm, 89.2% at 1.17 mm of
residual displacement. An invariant fraction is the signature of a quantity set by construction
rather than by the thing being measured, so the construction was asked for its value directly.

For a random isotropic field, the component tangent to the surface is invisible to **any**
surface-based correspondence. That floor is computable in closed form, with no reference to any
mesh, correspondence or warp:

| statistic | analytic | Monte Carlo, N = 4e6 |
|---|---:|---:|
| mean \|v_t\| / mean \|v\| | **pi/4 = 0.7854** | 0.7854 |
| rms \|v_t\| / rms \|v\| | sqrt(2/3) = 0.8165 | 0.8166 |
| mean \|v_n\| / mean \|v\| | 1/2 | 0.4999 |

A random spline field has iid Gaussian components, so its direction is isotropic and `pi/4` is the
floor this control cannot go below however good the correspondence is.

| run | measured fraction | floor | **excess** | excess, absolute |
|---|---:|---:|---:|---:|
| 1.17 mm residual | 0.8920 | 0.7854 | **+0.1066** | 0.125 mm |
| 1.50 mm field | 0.8787 | 0.7854 | **+0.0933** | 0.140 mm |
| 20.0 mm field | 0.8408 | 0.7854 | **+0.0554** | 1.11 mm |

**So roughly nine tenths of the alarming number is the control's own construction, and the
correspondence's actual contribution is the excess.** Reporting 84-89% as "target error" attributes
to nearest-point matching a quantity that a perfect correspondence would also incur. This is the
programme's standing rule applied to my own control: check the metric against a case whose answer
you know, and the answer here was available analytically before any run.

**What the instrument must report from now on**, fixed before the replacement control's numbers are
read: the error vector decomposed per point into components NORMAL and TANGENTIAL to the scaffold
surface. The tangential part is the identifiability floor and is reported as such; the **normal
part is the only component a correspondence can be blamed for**, and it is the only one that may be
called target error. The floor subtraction in the table above is the crude version of this and is
labelled as crude -- a per-point decomposition does not assume isotropy, and a real bone
displacement field is not isotropic.

**Not concluded:** the absolute excess grows 0.140 -> 1.11 mm from 1.5 to 20 mm, which is
sublinear in a 13.3x change, and d^2/R would be superlinear. That is one more reason the d^2/R
question stays open rather than answered in either direction, and the denominators in those two
rows do not mean the same thing (residual displacement against field amplitude), so the trend is
not read as a result at all.

### The decomposition lands, and "gate 1's margin is inside the bias" was overstated

The decomposition code was verified against the isotropic known answer first: it returns **0.7849**
against pi/4 = 0.7854, rms **0.8162** against 0.8165, normal **0.5005** against 0.5. At the first
amplitude (achieved separation 0.39 mm, residual displacement 1.17 mm):

| component | mean | share |
|---|---:|---:|
| whole error vector | 2.146 mm | 89.2% of residual displacement -- **not the correspondence's error** |
| tangential (identifiability floor) | 1.986 mm | **0.9255** of the error vector |
| **normal -- the correspondence's own error** | **0.494 mm** | p90 **0.966**, max **32.5** |

**The error vector is more tangential than isotropic** -- 0.9255 against pi/4 = 0.7854 -- which is
the case I told the agent to watch for, and the mechanism it gives is convincing: a best-fit
similarity starting map preferentially absorbs scale and translation, which is normal-direction
content, leaving a residual that is more tangential than a random field would be. **So the true
floor here is above pi/4 and my crude subtraction would have mis-credited it.** The per-point
decomposition was the right instrument and the closed-form floor was only ever a sanity check on it.

**Correcting myself, for the third time on this line and in the same direction.** `95c1e94` said,
and `d626af3` predicted, that "gate 1's margin of 1 mm is inside the bias". That rested on
**1.318 mm**, which is the whole error vector -- nine tenths of it a floor that a perfect
correspondence incurs too. The correspondence's own error is **0.494 mm mean**. So:

* the claim as written is **wrong at the mean**: the margin is about twice the bias, not inside it;
* it is **very nearly right at p90**, where the bias is 0.966 mm against a 1 mm margin;
* the honest statement is that **a result turning on a millimetre is unsafe for the worst tenth of
  points and defensible at the median**, which is a much weaker and much more useful claim than the
  blanket one I committed.

Three times today I read a composite number as a correspondence error: displacement for separation,
error vector for correspondence error, and a scaled control for a fixed one. The common shape is
the ledger's own -- a quantity computed correctly and compared against the wrong thing.

**The max of 32.5 mm in the NORMAL component is the thin-sheet failure**, arriving exactly where
predicted: a nearest-point target landing on the far wall of a rib. It is the pathology normal
shooting shares and the agreement filter exists to remove, so it strengthens rather than weakens
the queued instrument -- and it is now visible as a correspondence fault rather than buried inside
a number that was 93% floor.

### The retired control's real finding: fitting the targets is not recovering the truth

The retired 20 mm control ran to completion. Its verdict on this line is **not read** -- it is
retired for the reason in `73dfc35`, it operates at separations the pipeline never presents, and
the separation control supersedes it. But it measured one thing that does not depend on separation
at all, and that thing is the most useful number the control produced:

| | 1.50 mm field | 20.00 mm field |
|---|---:|---:|
| residual at the correspondences (**agreement with the TARGETS**) | 0.149 mm RMS | 0.437 mm RMS |
| recovery vs the TRUTH, pointwise | 1.510 mm RMS | 19.478 mm RMS |
| recovery to the SURFACE (the identifiable part) | 0.510 mm RMS | 3.345 mm RMS |
| **understatement factor, residual -> to-surface truth** | **3.4x** | **7.7x** |
| understatement factor, residual -> pointwise truth | 10.1x | 44.6x |

**The warp fits its targets beautifully and recovers the truth badly, at both scales.** It agrees
with the targets to 0.149 mm while sitting 0.510 mm from the truth on the identifiable component
and 1.510 mm away pointwise. The fit residual is not a small version of the error; it is a
different quantity, and it is optimistic by 3.4x at best and 45x at worst.

**This applies to a number this line has been quoting.** The correspondence residual of **0.571 mm**
has been carried as though it described fit quality. On this evidence it understates the
to-surface error by something like 3-8x, so the honest reading of that 0.571 mm is a to-surface
accuracy of **roughly 2-4 mm**, not half a millimetre. Every place a correspondence residual is
quoted as an accuracy on this line is affected. Nothing is restated as a corrected number yet,
because the understatement factor is measured here at the wrong separations; the separation
control reports the factor at the operating point, and only then is a corrected figure written.

**What it does NOT touch.** The chest-wall line's 0.078 mm came from normal shooting with the
agreement filter, a different instrument from nearest-point TPS, and this control says nothing
about it. And the factor is not a constant to divide by -- it moves 3.4x to 7.7x across the two
scales measured, so it is a warning that the quantity is wrong, not a conversion.

### At the operating point the margin IS inside the bias, and my retraction was made out of range

The separation control's second amplitude lands at a median separation of **1.05 mm**
(0.35-5.07 across segments), squarely inside the pipeline's measured 0.78-6.88 mm range:

| achieved separation | residual displacement | tangential fraction | **normal error (the correspondence's own)** | p90 | max |
|---:|---:|---:|---:|---:|---:|
| 0.39 mm | 1.17 mm | 0.9255 | **0.494 mm** | 0.966 | 32.5 |
| **1.05 mm** | 3.53 mm | 0.9009 | **1.933 mm** | 5.470 | 54.8 |

**Gate 1's 1 mm margin is inside the bias after all -- by a factor of about two.** The full
accounting of this claim, which I have now had wrong twice in opposite directions:

1. `d626af3` predicted and `95c1e94` asserted that the margin sits inside the bias, resting on
   **1.318 mm** -- the whole error vector, nine tenths of it an identifiability floor. Right
   conclusion, wrong quantity.
2. `8fc58bd` retracted it on the correctly-attributed **0.494 mm** normal error. Right quantity --
   measured at **0.39 mm separation, below the pipeline's own minimum of 0.78 mm**. Wrong
   separation, so the retraction was made out of range.
3. The claim is reinstated on **1.933 mm at 1.05 mm separation**: right quantity, in-range
   separation. It stands for a reason neither earlier version had.

That is the fourth time today this line has produced a number read at the wrong separation, and the
second time I have been the one to read it. The lesson has now cost a prediction, a retraction and
a stop rule.

**A structural prediction, fixed before the third amplitude lands.** The normal error grows from
0.494 to 1.933 mm as separation goes 0.39 to 1.05 mm -- an exponent of **1.38**, superlinear. Two
points cannot establish an exponent, so the third amplitude (achieved separation ~2.7 mm) is a test
rather than a confirmation, bracketed now:

| hypothesis | exponent | predicted normal error at 2.7 mm |
|---|---:|---:|
| linear: a fixed angular error | 1.00 | 4.97 mm |
| the measured 1-2 exponent, continued | 1.38 | 7.10 mm |
| **d^2/R curvature bias** | **2.00** | **12.78 mm** |

**Predicted: nearer 7 than 13** -- the 1-2 exponent continues rather than steepening to the
curvature law. This is the first test this line has run that can speak to d^2/R at all, after
`ae34bd2` correctly found the earlier control silent on it. If the third amplitude lands near
12.8 mm, d^2/R is supported and the separations must be driven down rather than the correspondence
improved; near 5 mm and the bias is an angular error that a better correspondence can remove. No
gate attaches to this -- it is a structural question, and it decides which repair is worth building.

### VERDICT: the bias is an angular error, not a curvature law. Normal shooting is the repair.

The third amplitude landed at an achieved separation of **2.20 mm** (0.61-12.32 across segments),
residual displacement 8.21 mm. Normal error **4.724 mm**, p90 12.710.

| achieved separation | normal error | fit residual | to-surface recovery | understatement |
|---:|---:|---:|---:|---:|
| 0.39 mm *(below the operating range)* | 0.494 mm | 0.242 mm | 0.984 mm | 4.1x |
| 1.05 mm | 1.933 mm | 0.258 mm | 2.138 mm | **8.3x** |
| 2.20 mm | **4.724 mm** | pending | pending | pending |

**Exponents, computed and reported separately rather than fitted:** 1->2 is **1.378**, 2->3 is
**1.208**, and across the full 5.6x span 1.305. Against the brackets of `bd3151c`, recomputed at
the achieved ratio because the run reached 2.20 mm rather than the 2.7 mm I assumed:

| hypothesis | exponent | predicted | measured 4.724 is |
|---|---:|---:|---|
| fixed angular error | 1.00 | 4.05 mm | +17% |
| the 1->2 exponent continuing | 1.378 | 5.36 mm | **-12%** |
| **d^2/R curvature bias** | 2.00 | 8.49 mm | **-44%** |

**By the rule fixed before the number existed, this decides the repair: normal shooting, not a
better starting map.** Two exponents of 1.208 and 1.378 across a 5.6x span of separation are not a
quadratic. `bd3151c` said "near 12.8 and d^2/R is supported... near 5 and the bias is an angular
error that a better correspondence can remove" -- it is 4.724, so the correspondence is the thing
to fix.

**Scoring my own prediction honestly: directionally right, centrally wrong.** I predicted "nearer
7 than 13", i.e. that the 1->2 exponent would continue rather than steepen. The non-quadratic call
is confirmed. But the exponent did not continue -- it **fell**, 1.378 to 1.208 -- and the measured
value sits 12% **below** my central bracket and only 17% above the linear one. The honest summary
is that I called the hypothesis correctly and the magnitude slightly wrong, in the direction of the
bias being even more benign than I expected. This is also the first prediction of mine today that
was not simply refused.

### The fit residual carries almost no information about accuracy

This is the sharpest form of `d4ce65e`'s finding, and it is stronger than "optimistic by a factor":

| between amplitudes 1 and 2 | change |
|---|---:|
| achieved separation | **+169%** |
| to-surface error against the truth | **+117%** |
| **fit residual at the correspondences** | **+6.6%** |

**The error more than doubles while the residual moves by a fifteenth.** The residual is not a
compressed version of the error, nor a constant multiple of it -- the understatement factor itself
goes 4.1x to 8.3x. It is close to uninformative about accuracy across the range this pipeline
operates in.

**Consequence for this line's quoted numbers.** The **0.571 mm** correspondence residual cannot be
converted into an accuracy by any factor, because the relationship is not a factor. It should stop
being quoted as fit quality entirely, and be replaced by the to-surface recovery where a truth is
available and by nothing at all where one is not. A residual this flat is a statement about how
well the spline interpolates its own targets, which was never in doubt.

### Normal shooting: 0.444 mm against 1.933 -- and the comparison is not yet like-for-like

| at the operating point | nearest point | normal shooting |
|---|---:|---:|
| achieved separation | 1.05 mm | 0.96 mm |
| **normal error** | **1.933 mm** | **0.444 mm** |
| as % of residual displacement | 27.4% | 10.7% |
| thin-sheet hits | pending | 5 of 3,989 (0.13%), torso |
| **kept fraction** | **100%** (10% trim) | **27.1%** |

**The repair predicted in `5239b2a` delivers, and it puts the error under gate 1's margin for the
first time** -- 0.444 mm against a 1 mm margin, where the nearest rule sat at nearly twice it. The
shooting arm's 1->2 exponent is **1.059**: it did **not** collapse toward zero, which is what I
asked to be watched for, but it sits just below my predicted 1.2-1.4 band. Same shape as the last
two predictions of mine -- hypothesis right, magnitude slightly wrong, erring toward the benign.

**But 27.1% kept is not a detail, it is the whole question.** Nearest point is scored on 100% of
its samples; shooting is scored on the 27.1% that survived `no_hit` (1,031) and `return_too_far`
(3,821). **Those are different populations, so 1.933 against 0.444 is not a measurement of the
two rules.** This is exactly what `fec11d2` pre-registered the drop census to catch -- a
correspondence that discards three quarters of its samples and reports an improvement may simply be
reporting the easy quarter. The nearest rule's own pathologies (normal-component maxima of 32.5 and
54.8 mm, far-wall hits through thin ribs) are precisely the points shooting drops, so the kept set
is plausibly enriched for points nearest-point also handles well.

**The control that settles it, ordered now and fixed before it runs: score the NEAREST-POINT rule
on exactly the subset shooting kept.** Same points, same amplitude, same everything else. That is
the only like-for-like number, and it is cheap because the subset is already known.

* **If nearest-point on the kept subset is still near 1.9 mm**, the improvement is real and
  entirely attributable to the rule.
* **If it falls to near 0.44 mm**, the improvement is selection and normal shooting's advantage at
  the operating point is an artefact of which points it declines to answer for.
* **Predicted: it lands between 0.8 and 1.5 mm** -- a genuine improvement over 1.933 because the
  subset excludes the far-wall pathologies that inflate the nearest rule's mean, but still two to
  three times worse than shooting's 0.444. If it lands at or below 0.6, I am wrong and most of the
  gain is selection.

**And the drop census is now a cost in its own right, not only a diagnostic.** It worsens with
separation -- 27.1% kept at the operating point, `return_too_far` up from 1,281 to 3,821 and
`no_hit` from 138 to 1,031. A correspondence that answers for a quarter of the skin cannot drive a
warp over the segments that fail unless those segments survive the filter. **Report the kept
fraction PER SEGMENT at the largest amplitude**, because a 27% average that is 60% on the torso and
5% on the calcaneus would leave exactly the failing segments unconstrained, and the headline number
would not show it.

### FOUR TIMES BETTER TARGETS PRODUCE A WORSE WARP. Target error is not a proxy for warp quality either.

| amplitude | normal error, nearest -> shooting | **to-surface recovery, nearest -> shooting** |
|---|---:|---:|
| 1 (below the operating range) | 0.494 -> 0.199 mm | 0.984 -> 0.908 mm |
| **2 (the operating point)** | 1.933 -> **0.444 mm (4.4x better)** | 2.138 -> **2.697 mm (26% WORSE)** |
| 3 | 4.724 -> 0.901 mm | 2.985 -> pending |

At amplitude 2 the separations are comparable (1.05 against 0.96 mm), the targets are **4.4 times
more accurate**, and the warp recovers the truth **26% worse**. **This inverts `5239b2a` and
`ded65bf`.** I authorised normal shooting on the strength of target error, and target error turns
out to point the wrong way.

**So the line now has three quantities that do not measure what they look like:**

1. the **fit residual** is nearly flat while the error doubles (`5239b2a`) -- uninformative;
2. the **target error** is 4.4x better while recovery is 26% worse -- **anti-correlated** at the
   operating point;
3. only **recovery to the surface against a known truth** has tracked anything, and it exists only
   because a control manufactures the truth.

Every number this line quoted before the recovery control was one of the first two. That is the
strongest possible argument for the order that was imposed -- *run the recovery control first,
judge no warp family until it lands* -- and it is the only reason this was caught rather than
shipped as a 4.4x improvement.

**The understatement ratio makes it worse, and in the informative direction.** It is **larger** for
shooting (10.7x, 15.6x) than for nearest point (4.1x, 8.3x, 6.0x): better targets lower the
residual without lowering the error, so the residual's optimism *grows* exactly when the
correspondence improves. A better instrument makes the bad metric look better and the warp worse.

**The proposed mechanism is coverage, and it is a hypothesis with a control, not a conclusion.**
Shooting keeps 27.1% and the fit is subsampled to an equal count, so the two fits have the same
number of constraints but not the same *spatial* coverage: whole regions where the surfaces are
awkward drop out, and those are where a warp most needs constraining.

**Control C, fixed before it runs.** Hold the correspondence rule constant and vary only coverage:
fit the warp from **nearest-point targets subsampled to shooting's kept SET** -- same rule, same
count, same spatial pattern as the shooting fit. Compare its to-surface recovery against the
full-coverage nearest fit's 2.138 mm.

* **Degrades to near 2.7 mm** -> coverage is the mechanism, the correspondence rule is exonerated,
  and the repair is to raise coverage (cap and filter are declared settings, not thresholds, so
  they may be changed after a FAIL; no gate attaches to recovery).
* **Stays near 2.1 mm** -> coverage is not the mechanism and something about the shooting targets
  themselves harms the fit, which would be a far more interesting and more troubling result.
* **Predicted: it degrades to 2.5-2.9 mm** and coverage carries it. Recorded so a confirmation is
  not read as more than the gate allows.

**A separate flaw in the shooting arm that must be fixed before its exponent is quoted again.** The
kept separation barely moves across amplitudes -- 0.45, 0.96, **1.08 mm** -- while the field
amplitude more than doubles between aims 3 and 7 (62.7 to 146.3 mm mean displacement), because
`return_too_far` and the cap preferentially discard the large-separation samples. **The kept set is
separation-biased**, so the 2->3 exponent of 6.01 is an artefact of a denominator that did not move
and must not be quoted as a scaling. Only the 1->2 value of 1.059 is defensible, and the nearest
rule's exponents stand unaffected because it keeps everything.

**Scoring my predictions, which is the uncomfortable part.** All three of my shooting predictions
landed: normal error fell substantially at every amplitude (2.5x, 4.4x, 5.2x), the exponent did not
collapse to zero, and the thin-sheet maxima improved most (32.5 -> 15.8, 54.8 -> 34.3). **I was
right about everything except whether it would help.** I forecast the correspondence accurately and
never asked the question that mattered, because I had accepted target error as the thing to
improve. Being right about three sub-quantities of a metric that turned out to be
anti-correlated with the outcome is not a good record; it is a demonstration that forecasting skill
on the wrong quantity is worth nothing.

### The 4.4x was mostly selection: on the same population it is 1.19x. My criterion fired against me.

`ded65bf` predicted nearest-point on shooting's kept subset would land at **0.8-1.5 mm**, and fixed
the refutation in advance: *"If it lands at or below 0.6, I am wrong and most of the gain is
selection."* **It is 0.529 mm.** The criterion fired, and the conclusion it forced is the one it
was written to force.

| amplitude 2 | nearest point | normal shooting | ratio |
|---|---:|---:|---:|
| target error, **all** samples | 1.933 mm | 0.444 mm | 4.4x "better" |
| **target error, same population** | **0.529 mm** | **0.444 mm** | **1.19x** |
| to-surface recovery | 2.138 mm | 2.697 mm | **26% worse** |
| fit residual | 0.258 mm | 0.173 mm | "better" while the error worsened |

**Three quarters of the headline improvement was the filter declining to answer for the hard
points.** The nearest rule is not 4.4 times worse than shooting; it is 1.19 times worse, and its
apparent gap was the far-wall pathologies that shooting drops rather than solves.

**So normal shooting, as constituted, is strictly worse.** It buys 1.19x on targets and pays 26% on
the warp -- the only quantity with a known truth behind it. The repair authorised in `5239b2a` does
not repair this line. Whether coverage explains the 26% is Control C's question and is still open;
what is closed is that the 4.4x was never real.

**The fit residual inverted a third time here**, and it is worth stating as a completed pattern:
0.258 -> 0.173 mm, *better*, on the arm whose warp is worse. Across this line the residual has now
moved independently of the error (flat while it doubled), and opposite to it (down while it rose).
It is not a degraded measure of accuracy. It is not a measure of accuracy.

**My scoreboard on this line today, since it bears on how much weight my calls should carry:** the
anchor rule's "small fraction" -- refused. Target error "near 1.3 mm absolute" -- refused. The
d^2/R bracket -- directionally right, centrally wrong. The shooting exponent's 1.2-1.4 band --
direction right, below the band. Nearest-on-kept-subset at 0.8-1.5 -- refused, at 0.529. Three
shooting sub-predictions -- all landed, all worthless, because the metric was anti-correlated with
the outcome.

**What did hold was the sequencing, not the forecasting:** recovery control first, no warp family
judged until it lands, every criterion written before its number existed. That order caught all
three inversions, and each one was invisible in the quantity the line had been quoting. The
discipline earned its keep today; the predictions did not, and the honest reading is that my
intuitions about this instrument are not calibrated and should be used to generate controls rather
than to anticipate their results.

### Control C: coverage is the mechanism, and it outweighs the rule ten to one

**Predicted 2.5-2.9 mm. Got 2.758 mm.** The first prediction of mine on this line to land inside
its band, and it settles the question `0213dfc` opened:

| at the operating amplitude | full-coverage nearest | **Control C** (nearest, shooting's coverage) | shooting |
|---|---:|---:|---:|
| to-surface recovery | **2.138 mm** | **2.758 mm** | 2.697 mm |
| same-population normal error | 0.529 mm | 0.529 mm | 0.444 mm |

**The whole trade, measured: coverage costs 0.620 mm of recovery; the rule buys 0.061 mm back.**
Shooting is very slightly *better* than the nearest rule at matched coverage (2.697 against 2.758),
so the 26% degradation in `0213dfc` is explained entirely by coverage and not at all by the
targets. The correspondence rule was never the problem and is not the solution either -- it is a
tenth the size of the effect that is.

**Correcting a number I propagated into two commits.** `ded65bf` and `0213dfc` say the shooting rule
keeps **27.1%**. It does not: the rule keeps **62.8%** at the operating amplitude and **33.5%** at
the largest. The 27.1% was measured *after* count-matched subsampling -- an artefact of the fitting
procedure, not the correspondence's behaviour -- and I wrote it up as the rule's loss rate. The
conclusions in those entries survive because they turned on coverage being reduced, which it is;
the magnitude of the loss was overstated by more than a factor of two.

**And my specific worry was exactly backwards.** `ded65bf` warned that "a 27% average that is 60% on
the torso and 5% on the calcaneus would leave exactly the failing segments unconstrained." The
measured per-segment rates invert that: at the operating amplitude **torso 6.8%** and **pelvis
29.7%** are the starved ones, while **calcn 51.8%**, **toes 40.8%** and **patella 79.8%** stay well
covered. The segments that fail gate 2 are **not** the ones shooting starves. Torso and pelvis are,
and at the largest amplitude the starvation spreads to the long bones (tibia_l 6.8%, humerus_r 10%,
femur_l 10.7%).

**That inversion is coherent rather than puzzling.** The torso is where every thin-sheet hit
occurred. Ribs and scapulae both defeat the shot -- producing `no_hit` and `return_too_far` -- and
inflate the nearest rule's tail through far-wall matches. One geometry causes both failures, which
is why the filter that fixes the tail is the same thing that starves the coverage.

**So the repair is neither rule, and it is pre-registered here before it is built.** Since coverage
outweighs targets ten to one, a **hybrid** -- shooting's target where the shot is kept, the nearest
rule's target where it is dropped -- should beat both, *despite* reintroducing the nearest rule's
bad targets on exactly the worst geometry. That is a counterintuitive consequence of the measured
trade and is worth testing precisely because it is counterintuitive.

* **Predicted: the hybrid recovers 2.00-2.13 mm** -- at or slightly better than full-coverage
  nearest's 2.138, since it restores 100% coverage and improves the targets on 62.8% of it.
* **If it lands below 2.00**, the two effects combine better than additively and the line has a
  real repair.
* **If it lands above 2.20**, reintroducing far-wall targets on the torso costs more than the
  coverage it buys, and coverage and target quality are not separable the way this table implies.
* Either way, **the entire normal-shooting exercise is then worth under 0.15 mm**, which is the
  honest scale of what a correspondence change can do for this line, and is worth stating before
  anyone builds a third one.

### The hybrid: 1.686 mm. Predicted below 2.00 means a real repair, and it is.

| arm, amplitude 2 | mean target error (normal) | **to-surface recovery** |
|---|---:|---:|
| nearest, full coverage | 1.933 mm | 2.138 mm |
| Control C (nearest, reduced coverage) | 0.529 mm* | 2.758 mm |
| shooting | **0.444 mm** (best targets) | 2.697 mm |
| **hybrid** | 1.962 mm (**worst targets**) | **1.686 mm** (best recovery) |

\*same-population figure. Predicted 2.00-2.13; **below 2.00 means the effects combine better than
additively and the line has a real repair.** It is 1.686 -- 21% better than the best previous arm.

**The arm with the WORST mean target error has the BEST recovery, and the arm with the best targets
is third of four.** Ranked by mean target error the order is very nearly the reverse of the ranking
by recovery. This is the day's theme at its sharpest: **mean target error does not predict recovery
on this line, and never did.**

**The proposed mechanism sharpens it further, and it is mine to test rather than accept.** The
hybrid gives shot-quality targets to 62.8% of points and matches the nearest rule only on the hard
torso geometry. Since the torso tail dominates the mean, the mean hides an improvement spread
across two thirds of the body. If that is right, **mean target error was the wrong SUMMARY
STATISTIC, not merely the wrong quantity** -- a distinction worth having, because the first is
repairable and the second is not.

**Pre-registered before it is measured: per-segment recovery for the hybrid against full-coverage
nearest.** The two differ only in which targets the shot supplies, so the comparison is clean.

* **Predicted: the hybrid's advantage is concentrated where its shot share is high, and is near
  zero on the torso** (shot share 6.8%), where it is using nearest targets anyway. If the torso
  improves as much as the well-covered segments, the mechanism is wrong and something other than
  target provenance is driving the gain.
* This is a mechanism test with no gate attached, and it is the cheap version of the question that
  would otherwise be answered by building a fourth correspondence.

**Correcting a ceiling I set one entry ago.** `380df28` concluded "the entire normal-shooting
exercise is then worth under 0.15 mm". That was wrong, and wrong for a reason worth naming: it
assumed shooting would be used as a **replacement** for the nearest rule. As a **component**, with
its drops backfilled rather than discarded, it is worth **0.452 mm** of recovery (2.138 -> 1.686) --
three times the ceiling I set. The measurement that produced the ceiling was sound; the framing
around it was not, because I compared two rules when the answer was to use both.

### The mechanism is refuted, the hybrid's win is confounded, and the ceiling correction is suspended

**Predicted: the hybrid's advantage concentrates where shot share is high and is near zero on the
torso. Refuted.** Per-segment, against full-coverage nearest:

| | measured |
|---|---|
| Pearson r (gain vs shot share) | **-0.353**, p = 0.108, **r^2 = 0.124** |
| Spearman | -0.401, p = 0.064 (n = 22) |
| segments improving | **22 of 22**, range -6.2% to -38.2% |
| mean gain | **-22.7%**, sd 8.2% |
| torso | **-18.2%** at 6.5% shot share |

The correlation runs in the predicted direction, is weak, is not significant, and explains about an
eighth of the variance. The dominant feature is a **near-uniform ~23% gain across every segment**,
including those taking almost no shot targets. Target provenance is at best a minor term and is
**not** what drives the hybrid's win.

**And the comparison is confounded, which the agent found in its own construction.** The
full-coverage nearest arm discards its largest-separation **10% per segment**; the hybrid discards
none. The two arms therefore differ in targets *and* in trim, and a trim that removes the hardest
points from one arm and not the other would produce exactly this signature -- a uniform improvement
independent of provenance. Given that coverage has already been shown to dominate targets ten to
one (`380df28`), a 10% difference concentrated on the hardest points is a live explanation and
possibly the whole one.

**Withdrawn: the ceiling correction in `80afd28`.** I revised "the whole normal-shooting exercise is
worth under 0.15 mm" up to **0.452 mm as a component**, on the strength of 2.138 -> 1.686. That
comparison cannot carry the attribution, so the revision is suspended. The defensible position
returns to `380df28`'s **under 0.15 mm**, which was measured at matched coverage and is not touched
by the trim. I corrected a sound figure on the basis of a confounded one, within an hour of setting
it.

**What survives untouched**, because none of it turns on the hybrid's cause:

* Control C: coverage is the mechanism for shooting's deficit (2.138 -> 2.758);
* the subset control: 1.19x on a matched population, most of the 4.4x was selection;
* **the hybrid is genuinely better, across all three amplitudes** (-24%, -21%, -13%). The *result*
  stands; only the *explanation* is missing.

**The isolating design, pre-registered before either fit runs.** It completes a 2x2 whose diagonal
is already filled -- {nearest, hybrid} x {trimmed, untrimmed}, at lambda = 1e-3 with no CV needed:

| | trimmed 10% | untrimmed |
|---|---:|---:|
| nearest targets | **2.138 mm** (have) | **cell A** |
| hybrid targets | **cell B** | **1.686 mm** (have) |

* **Trim explains it:** A near 1.69, B near 2.14.
* **Targets explain it:** A near 2.14, B near 1.69.
* **Predicted: trim carries the majority but not all** -- A lands **1.75-1.95**, B lands
  **1.95-2.15**. Coverage has beaten target quality on this line three times now (Control C, the
  selection control, the coverage restoration), and the trim is a coverage manipulation aimed at
  precisely the hardest points.
* Stated with low confidence and recorded anyway: my predictions on this line stand at two hits in
  eight, and the design is worth more than the guess. **The 2x2 is decisive whichever way it falls**,
  which is the property being bought here.

### THE TRIM IS THE WHOLE EFFECT. The correspondence rule is worth 0.02 mm.

| to-surface recovery, amplitude 2 | trimmed 10% | untrimmed | **trim effect** |
|---|---:|---:|---:|
| nearest targets | 2.138 mm | **1.666 mm** | **-0.472 mm** |
| hybrid targets | 2.113 mm | 1.686 mm | **-0.427 mm** |
| **target effect** | -0.025 mm | +0.020 mm | |

Both rows agree on the trim to within 0.045 mm. Both columns agree that swapping nearest targets
for hybrid targets moves recovery by **±0.02 mm, which is nothing.** The factor this line has spent
the day investigating does not matter; a parameter nobody was looking at accounts for all of it.

**Scoring the prediction:** cell B landed inside its band (2.113 against 1.95-2.15); cell A fell
below it (1.666 against 1.75-1.95). I said trim would carry "the majority but not all". It carries
effectively all. Directionally right, magnitude under-called -- which puts my record on this line at
roughly two and a half hits in nine.

**This is not a property of the control. It is in the production pipeline.** `CORR_TRIM = 0.10` is
applied to the correspondences in the real instrument: **v1 spline, v2 flow and v3 anchored all
discard the largest-separation 10% per segment.** On this control that choice costs **22% of
recovery accuracy** while barely touching target error (1.933 -> 2.036 mm untrimmed). The points it
removes are the ones carrying the most information about where the surfaces disagree.

**How it got there is the lesson.** It was adopted in the v1 pre-registration as "the same 10% trim
the per-segment fit itself uses", reasoning by analogy with ICP's outlier rejection. It was never
measured. **A pre-registration protects against choosing a threshold after seeing a result; it does
nothing about a parameter adopted by analogy before any result exists.** Those are a different
failure and this programme had no guard against them. Recorded in `CLAUDE.md`.

**Consequences, stated plainly:**

* **Normal shooting's value on this line is ~0.02 mm** -- below even the "under 0.15 mm" figure
  restored in `269f1bb`. As a replacement or as a component, the correspondence rule is not where
  the accuracy is, and the whole normal-shooting branch was a well-run investigation of a term that
  does not matter.
* **The trim is where the accuracy is**, and removing it is an instrument change requiring its own
  pre-registration, below.
* **Every gate-1 number was produced with the trim in place.** If the trim goes, those fits must be
  **re-run, not re-read**. That includes v3, currently five hours into a solve that would become
  obsolete -- a cost worth paying, since the anchored family is closed as an idea anyway.

**Pre-registered before the sweep runs.** The 10% was arbitrary as well as unmeasured, so the
question is not only whether to remove it but whether its cost is monotone. Sweep
`CORR_TRIM` over **{0, 0.02, 0.05, 0.10}** on the control at amplitude 2, nearest targets, all else
fixed.

* **Monotone** -- recovery improves all the way to 0 -- **adopt 0** and delete the parameter rather
  than retune it.
* **Interior optimum** -- some small trim beats both 0 and 0.10 -- adopt it, label it **post-hoc**,
  and state that it was chosen on a synthetic control.
* **Predicted: monotone to 0.** The 22% cost at 10% and the ±0.02 mm target-error insensitivity
  both say the trim is removing information rather than noise.

**The caveat that must travel with any decision here, and it is not small.** This control's truth is
a known field applied to *this same body*, so its largest separations are hard but genuine. In a
real subject-to-scaffold registration the largest separations may be **wrong** correspondences
rather than merely hard ones, and a trim that costs accuracy here could be protective there. There
is no ground truth on real data -- that is the reason this control exists at all -- so the only
evidence available says remove the trim, and no evidence says keep it. That is a weaker warrant
than it looks, and what would overturn it is a real registration with independent truth, which this
programme does not have.

### The sweep is monotone: CORR_TRIM goes to 0 and the parameter is deleted

| CORR_TRIM | points fitted | **to-surface recovery** | normal target error |
|---:|---:|---:|---:|
| **0%** | 4,410 | **1.666 mm** | 2.036 mm |
| 2% | 4,410 | 1.729 mm | 2.027 mm |
| 5% | 4,410 | 1.847 mm | 1.968 mm |
| 10% | 4,410 | **2.138 mm** | **1.925 mm** |

**Predicted monotone to 0. It is monotone to 0**, with no interior optimum, and the fitted count is
identical at every step -- so this is purely *which* points are kept, not how many. The shipped
setting is the worst of the four, costing **0.472 mm, 28%**. The 10% row reproduces the original
nearest arm's 2.138 mm to the last digit, confirming the sweep is the same instrument.

**The two columns move in opposite directions**, and that is the day's theme in its final form:
trimming makes target error look **5.8% better** while making the recovered map **28% worse**. The
trim removes the points where the correspondence is least certain, which are exactly the points
carrying the information about where the surfaces disagree. **A metric that improves as the map
degrades** -- the same failure found earlier in the fit residual, now in a second quantity.

**Adopted: `CORR_TRIM = 0`, and the parameter is deleted rather than retuned.** A tuned constant
invites the next person to tune it; there is nothing here to tune.

### The gate-1 re-runs, and why their outcome is NOT evidence about the trim

Fixed before a single fit is re-run, because this is exactly where a favourable reading could be
taken later. **Gate 1 measures the warped bone group's RMS against the per-segment similarity, on
the per-segment fit's own evaluation samples.** Structurally that is *agreement with what the fit
aimed at* -- the same class of quantity as the fit residual and the target error, both of which
this line has now shown can move **opposite** to accuracy.

* The re-runs are a **compliance check on the shipped instrument**: whatever ships must satisfy its
  gates. They are **not** a test of the trim decision, which rests on the control, which is the only
  place a truth exists.
* **A gate-1 improvement is not confirmation** that removing the trim was right, and **a gate-1
  degradation is not refutation.** Both are stated now so neither can be claimed later.
* Gate 1's threshold and its 1 mm margin are **unchanged**. Only the instrument changes, which is
  the permitted move.
* Old and new numbers are reported **side by side, all 22 segments**, not as a summary.

**The conflict case, named in advance and not resolved by preference.** If gate 1 **fails** at
`CORR_TRIM = 0` where it passed at 0.10, the honest position is that this programme has an
instrument that is **more accurate on the only data with a truth** and **fails a gate on the data
without one**. That conflict is reported as a conflict. It is not resolved by shipping whichever
setting passes, and it is not resolved by declaring gate 1 unfit because we dislike its answer --
gate 1 caught real things and its 22/22 pass was earned. **Predicted: gate 1 holds at 0**, with low
confidence stated; my record on this line is two and a half hits in ten.

**The caveat carried on every artefact from this**, restated because the monotonicity makes it look
stronger than it is: the control's truth is a known field on the *same body*, so its largest
separations are hard but genuine. On a real subject-to-scaffold registration they may be **wrong**
correspondences, where a trim could be protective. No ground truth exists on real data -- that is
why the control exists. The only evidence available says remove the trim and none says keep it.

### Gate 1 at CORR_TRIM = 0: PASS 22/22, better on every segment — and that is not evidence

| segment | per-seg similarity | trimmed | **untrimmed** | change | | segment | per-seg | trimmed | **untrimmed** | change |
|---|---:|---:|---:|---:|---|---|---:|---:|---:|---:|
| calcn_l | 6.09 | 5.37 | **4.59** | −14.5% | | pelvis | 10.41 | 9.05 | **7.84** | −13.4% |
| calcn_r | 6.00 | 5.10 | **4.69** | −8.0% | | radius_l | 1.77 | 1.40 | **1.17** | −16.4% |
| femur_l | 3.88 | 3.92 | **3.19** | −18.6% | | radius_r | 1.65 | 1.58 | **1.37** | −13.3% |
| femur_r | 3.94 | 3.90 | **3.51** | −10.0% | | talus_l | 3.23 | 2.62 | **2.09** | −20.2% |
| hand_l | 3.96 | 3.39 | **2.99** | −11.8% | | talus_r | 3.12 | 1.78 | **1.67** | −6.2% |
| hand_r | 3.98 | 3.20 | **2.73** | −14.7% | | tibia_l | 5.28 | 4.59 | **4.30** | −6.3% |
| humerus_l | 5.38 | 3.77 | **3.21** | −14.9% | | tibia_r | 5.29 | 4.46 | **4.16** | −6.7% |
| humerus_r | 4.91 | 3.41 | **2.71** | −20.5% | | toes_l | 4.15 | 3.89 | **3.80** | −2.3% |
| patella_l | 2.47 | 1.26 | **1.04** | −17.5% | | toes_r | 4.15 | 3.85 | **3.75** | −2.6% |
| patella_r | 2.52 | 1.32 | **0.95** | −28.0% | | torso | 15.65 | 11.28 | **9.84** | −12.8% |
| ulna_l | 1.97 | 1.81 | **1.56** | −13.8% | | ulna_r | 2.01 | 1.81 | **1.60** | −11.6% |

**PASS 22/22, mean −12.9%, and not one segment worse.** `femur_l`, which was 0.04 mm *over* the
per-segment similarity at trim 0.10 and passed only on the 1 mm margin, is now 0.69 mm under it.

**And by `8caf15b` this is not evidence that removing the trim was right.** That entry fixed, before
any fit ran, that gate 1 measures agreement with what the fit aimed at and that neither an
improvement nor a degradation counts. A clean sweep across 22 segments is precisely the result that
would tempt a reclassification of gate 1 as more meaningful than it was declared to be, and the
reclassification is refused. The trim decision rests on the control, which has a truth. **Gate 1's
pass means the shipped instrument satisfies its gate — nothing more.**

**One observation recorded without being interpreted.** Two agreement-type quantities moved in
*opposite* directions on the same change: cross-validated correspondence RMS **worsened** 4.588 →
5.284 mm (+15%), while gate 1 **improved** on every segment. If both were simply "agreement with
targets" they should have moved together. They did not, so they are not measuring the same thing --
but working out which is which needs its own control, and no conclusion is drawn here. Noting it
because it is the kind of fact that gets quietly used later as an argument that gate 1 is
trustworthy after all.

**A limitation the agent stated and I am carrying forward.** In the 2×2, the *rows* isolate the trim
exactly, but the *columns* also differ in sampler, which was stated and not controlled. So
"the correspondence rule is worth ±0.02 mm" carries an uncontrolled difference and should be read as
"no effect was detectable against a sampler difference of unknown size", not as a measured zero. The
trim result, which is the row comparison, is unaffected.

### Gate 2 on the untrimmed fit: unchanged. The trim gain does not show on real data.

| | ceiling | binding map | **trimmed** | **untrimmed** |
|---|---:|---:|---:|---:|
| mean enclosure | 0.997 | 0.888 | **0.954** | **0.951** |
| segments ≥ 0.99 | 20/22 | 9/22 | **12/22** | **12/22** |

**Gate 2 still FAILS, and removing the trim moved it by −0.003 — marginally worse.** The worst
segments are unchanged in identity and nearly in value: toes **0.854 / 0.876**, hands 0.792 / 0.818,
pelvis 0.905, calcn 0.910 / 0.912, torso 0.918.

**This is the first real-data test of the trim change, and it is negative.** Set against what the
same change bought elsewhere:

| quantity | effect of removing the trim |
|---|---|
| to-surface recovery on the control *(has a truth)* | **−28%**, 2.138 → 1.666 mm |
| gate 1 *(agreement with what the fit aimed at)* | better on **22/22**, mean −12.9% |
| **gate 2 *(skin enclosing bone, real data, no truth)*** | **−0.003 mean, 12/22 either way** |

Three quantities, three behaviours. Gate 1's improvement I had already ruled inadmissible in
`8caf15b`, before it was measured. Gate 2 is not of that kind — it is a geometric containment of
bone by skin, not an agreement with any target the fit chased, and the warp moves skin by up to
56.33 mm, so it is sensitive to exactly what changed. **It shows no benefit.**

**So the practical value of removing the trim is, on today's evidence, unestablished on real data.**
The 28% is real and measured against a truth; it simply has not appeared in the one real-data
quantity that could show it. That does not reinstate the trim — nothing here argues for putting it
back, and the control remains the only place a truth exists — but it does bound the claim. Anything
quoting "removing the trim is worth 28%" should say *on a synthetic field on the same body*, and
should say that the one independent real-data gate moved by −0.003.

**And it sharpens what the toes are.** They have now survived: the anchor rule being wrong in both
its forms, the correspondence rule turning out not to matter, the trim removal that improved
everything else, and a 28% better warp. Gate 3 adds the same finding from another direction — toe
skin sits **29.3 / 27.6 mm below** the toe bone's lowest point, reported there because gate 3 is
defined on calcn only. Nothing this line has tried has moved the toes, and the standing explanation
survives untouched: the scaffold's forefoot is a fifth longer than this specimen's, which is a shape
difference that no correspondence, regularisation or sampling change can reach.

### The toes cannot pass gate 2. Their own ceiling is 0.966 and the bar is 0.99.

The enclosure table's first column is the ceiling -- *"the body's own anatomical bones against the
same whole skin, same frame"*, the same `enclosure()` and the same samples, with only the skin's
vertex positions differing between columns. It reads **1.000 for twenty segments and 0.966 / 0.968
for toes_l / toes_r.**

**Gate 2's bar is ≥ 0.99. The toes' ceiling is 0.966. They cannot pass, and never could.** The
script has been printing this all along, in its own summary row: **`>= 0.99   20/22`** *on the
ceiling column*. Two segments fail at ceiling, and the two are the toes. Nobody read that row --
including me, through four rounds in which the toes were the headline failure.

**So every round of work aimed at the toes was aimed at a gate they are barred from.** The anchor
rule, the correspondence rule, the trim removal, a 28% better warp: all of it was measured against
a bar that the canonical anatomy itself misses by 2.4 points. The toes' 0.854 / 0.876 is a real
shortfall against their own 0.966 ceiling -- 0.11 and 0.09 of genuine gap -- but "toes FAIL gate 2"
was never the informative statement, and the count that matters is **12 of 20 reachable**, not 12
of 22.

**This is the third gate on this programme found to be unreachable by construction**, after gate 3's
"inside ≤ 1%" against a 42-53% ceiling (`61e5d1c`) and gate A's 0.5 mm bar on a 0.51 mm noise floor
(`de88ef4`). The shape is identical each time: **a bar set in absolute terms, against a quantity
whose achievable maximum was measured but never compared to the bar.** The ceiling column existed
from the start; it simply sat next to the bar without anyone subtracting.

**Fixed here, and NOT applied retroactively.** Gate 2's bar stays at 0.99 and every recorded verdict
stands as recorded -- the toes' FAILs are not rescored into passes, because a gate is not loosened
after a result. What changes is the *reporting*: the ceiling is quoted beside the bar wherever gate 2
appears, the pass count is given as **n of 20 reachable** with the two barred segments named, and a
segment whose ceiling is below the bar is reported as **BARRED** rather than FAILED. That is a
statement about what the number means, not a change to the number.

**And the scaffold defect underneath it is now the finding.** The canonical body's own toe bones are
**3.4% outside its own toe skin**. That is not a registration error and no warp can be blamed for
it: this body's toe skin does not contain its toe bones. It belongs with the forefoot measurement
already recorded -- per-segment similarity scale **1.22** at calcn and **1.21** at toes against the
global map -- as a property of the scaffold's foot rather than of anything fitted to it.

**Why this was worth having from a retired control.** I retired it for firing a threshold at the
wrong separation, and it then answered a question I had not asked: whether this line's headline
quality metric measures accuracy at all. It does not. A control kept running after its gate was
withdrawn is how that surfaced.

#### A known answer this line has never had, and what the chest wall found without it

Every gate above compares a fit to another fit. Gate 1 asks whether the warped bone group sits no
further from the scaffold than the per-segment similarity leaves it; the zero-warp control asks
whether the pipeline reproduces the binding map. Neither asks the question a known answer asks:
**displace this body's own bones by a warp we chose, and does the fit return it?** The chest-wall
line (`docs/BODY_PARAMETERS.md`) built exactly that test and **failed it at 1.194 mm RMS against a
1 mm bar**, with its own fit residual at the correspondences reading 0.01 mm.

The cause is the correspondence, and this line uses the same one. Targets here are "the nearest
point ON the scaffold's bone mesh" (`fit_skin_warp.py`); measured against a known true preimage,
a nearest-point target is off by **mean 0.671 mm, p90 2.281 mm, max 7.152 mm** for displacements
averaging 1.53 mm -- the d²/R bias nearest-point matching carries on a curved surface, and on a
rib a few millimetres thick that is most of the signal. A field applied purely along the surface
normals, identifiable by construction, recovered no better, so it is not the well-known tangential
blindness alone.

**Additional control, fixed now, to run after the anchored fit reports** -- an addition, not a
change to any gate: displace this body's own bone groups by a known smooth field of the fitted
family, fit that field with this line's own pipeline and correspondences, and report the recovery
RMS both pointwise and to the surface. **Predicted:** it lands near the chest wall's 1.2 mm rather
than near this line's 0.571 mm correspondence residual, because the residual measures agreement
with targets and the control measures agreement with the truth. If it does, gate 1's margin of
1 mm is inside the bias, and no result here that turns on a millimetre can be read as a
transform's doing.