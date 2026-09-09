# the body's tissue as mechanics, and where it stops being possible

645 of the bound body's 4,000 entities are tissue structures — ligaments,
capsules, menisci, discs, cartilage, bursae, fascia, retinacula, adipose — every
one of them with real mesh geometry, and the plant carried force elements for
**none** of them. They existed as geometry and not as mechanics.

**117 of them now carry force**, and **66 of those hold the plant inside its own
declared joint ranges better than the engineering joint stops they replace** —
6.08° of worst excursion against the stops' 5.30° and a bare plant's 17.43°, with
no stops present at all. The other 51 make the plant worse, for a reason that is
measured rather than guessed.

This file records that, and — the part that matters most — **the specific reason
each remaining class cannot be closed on this plant**, because they are blocked
for three different reasons and only one of the three is a missing solver.

Everything here is from `scripts/build_tissue_force_elements.py`,
`scripts/measure_tissue_mechanics.py`, `scripts/measure_ligament_moments.py`,
`scripts/measure_articular_substrate.py` and
`scripts/measure_muscle_path_constraint.py`, over the `engineering_stance_v1`
pose at 77.6122029 kg.

---

## Where the 645 classified structures went

Classified by name from the atlas' own labels, with the same first-match order
`scripts/audit_joint_substrate.py` uses so both lanes bucket a name the same way.
*Two-segment* means the structure's own surface votes for two different scaffold
rigid bodies; *one-segment* means the two bones it joins are the **same** rigid
body here, so no relative motion exists for it to resist.

| class | entities | two-segment | one-segment | force elements written |
|---|---:|---:|---:|---:|
| ligament | 300 | 105 | 195 | **105** |
| bursa / synovium | 80 | 12 | 68 | 0 |
| meniscus / disc | 69 | 10 | 59 | 0 |
| cartilage | 51 | 2 | 49 | 0 |
| tendon (mostly sheaths) | 44 | 36 | 8 | 0 |
| fascia / aponeurosis | 44 | 25 | 19 | 0 |
| joint capsule | 36 | 12 | 24 | **12** |
| retinaculum | 18 | 10 | 8 | 0 |
| adipose | 2 | 2 | 0 | 0 |
| skin | 1 | 1 | 0 | 0 |
| **total** | **645** | **215** | **430** | **117** |

**117 force elements where there were zero.** The plant now carries
`Blankevoort1991Ligament` — tension only, zero below slack, a quadratic toe to
6% strain and linear above, damped only while stretched and lengthening.

The 430 one-segment structures are not a failure of the derivation. **They are
the scaffold reporting its own resolution.** 59 of the 195 one-segment ligaments
are in `torso`, because the scaffold has one torso body and the whole spine and
rib cage are inside it; 29 per side are in `hand`, which is one rigid body
carrying 27 bones. A sacrotuberous ligament runs from sacrum to ischium and both
are `pelvis` here — so it correctly comes out one-segment, and that is a gate,
not a gap.

---

## Ligaments: how an attachment was derived, and what proves it

The recorded blocker was structural. `data/derived/anatomy-segment-binding/`
assigns each whole entity to exactly **one** segment by nearest-bone-group vertex
vote, and a ligament spans two, so the binding artifact cannot supply a
two-ended attachment. It cannot, **but the vote it is built from can**: the same
per-vertex nearest-bone-group distances that pick the winner also partition the
surface between two segments.

1. Every vertex takes the segment whose anatomical bone group is nearest. The
   two segments with the most votes are the two ends.
2. Each end's attachment is the centroid of the quarter of that end's vertices
   **furthest from the other end's centroid** — the tip, not the middle. Taking
   the centroid of each half gives an ACL 20.9 mm long against a published 32;
   the tip quantile gives 33.4.
3. Points are carried into the two segments' own frames through the binding's
   own recorded similarity and reference pose. **No new registration is fitted**,
   so this inherits the binding's 24.7 mm RMS residual and adds nothing.
4. Slack length is the separation at that reference pose, so every element is
   exactly slack in the pose its geometry was registered in.
5. Linear stiffness is `E·A`, with `E` the body's own declared ligament
   along-fibre modulus — 332.2 MPa, Quapp and Weiss 1998, human MCL — and `A`
   the structure's own tissue volume over its own derived length.

### The volume rule matters more here than where it came from

An **open** surface's signed integral is not a volume; it is the volume of what
the sheet wraps. `ihm/assembly/tissue_materials.py` owns that argument already
(84 open connective sheets were carrying 29.6 L of the muscle they wrap). Here
it decides a stiffness, because a cross-section is a volume divided by a length:
the hip capsule is an open sleeve whose signed integral is the femoral head and
neck inside it, and taking that as capsule tissue gave a **2,781 mm² cross-
section and a 924 kN ligament**. Routed through `material_volume` it gives
101 mm² and 33.6 kN.

### Gates — against cases whose answer is known

**Named pairs, 30/30.** A ligament everybody can name has two bones everybody
knows, and all of these came out right without the derivation ever seeing a name:

| structure | derived | | structure | derived |
|---|---|---|---|---|
| anterior cruciate | femur + tibia | | iliofemoral (descending) | pelvis + femur |
| posterior cruciate | femur + tibia | | pubofemoral | pelvis + femur |
| tibial collateral (both parts) | femur + tibia | | ischiofemoral | pelvis + femur |
| fibular collateral | femur + tibia | | anterior/posterior talofibular | tibia + talus |
| annular ligament of radius | radius + ulna | | calcaneofibular | tibia + calcn |

And the **negative controls**, which are the half that could have gone wrong
silently: sacrotuberous, sacrospinous and inguinal must come out **one** segment,
because both their bones are inside `pelvis`. They do. A two-segment answer there
would have been the derivation inventing a joint.

**A structure inside one bone.** Six molar and premolar teeth, run through the
same vote, return a second-segment share of **exactly 0.0**.

**Published lengths**, reported as ratios rather than asserted as agreement:

| ligament | derived | published | ratio |
|---|---:|---:|---:|
| anterior cruciate | 33.4 mm | 32 | 1.04 |
| posterior cruciate | 37.5 mm | 38 | 0.99 |
| fibular collateral | 52.9 mm | 60 | 0.88 |
| calcaneofibular | 19.1 mm | 25 | 0.76 |
| superficial tibial collateral | 70.8 mm | 95 | 0.75 |
| anterior talofibular | 14.1 mm | 20 | 0.70 |

Systematically short, by about 20% at the median: the tip quantile pulls in from
the true insertion. **Published Blankevoort stiffnesses**: ACL 10.1 kN against
5 kN (2.03×), PCL 11.6 kN against 9 kN (1.29×), LCL 0.87 kN against 2 kN (0.44×).
Within a factor of two either way, and *reported*, not tuned to.

**Cross-implementation length.** The engine's `Blankevoort1991Ligament` path
length against the same length computed from the same coordinates through the
model XML's own joint definitions: **max 3.3e-16 m over 117 elements**. Two
independent implementations of the same geometry.

**Reference-pose slack.** Maximum error **exactly 0.0 m** — every element is
slack in the pose its geometry was registered in.

---

## What they do to the plant, and what they cost

25 advances of 10 ms from the stance pose, same mass, same excitation, only the
force set differing. `admissible` is the 66-element subset the filter two
sections below selects.

| arm | elements | weight at t=0 | momentum residual | s/advance median |
|---|---:|---:|---:|---:|
| `bare` | 0 | **761.3757 N** | 4.9e-14 | 0.073 |
| `stops` — what the plant does today | 0 | **761.3757 N** | 4.9e-14 | 0.064 |
| `ligaments` | 105 | **761.3757 N** | 6.4e-13 | 0.215 |
| `ligaments_capsules` | 117 | **761.3757 N** | 3.2e-12 | 0.303 |
| `stops_ligaments` | 105 | **761.3757 N** | 6.4e-13 | 0.162 |
| `admissible` | 66 | **761.3757 N** | 1.1e-12 | 0.134 |
| `stops_admissible` | 66 | **761.3757 N** | 1.1e-12 | 0.106 |

**Standing weight is 761.3757 N in every arm, to the last digit** — 77.6122029 ×
9.81. That is the gate that says the new elements did not invent or lose external
force, and they are internal forces so they were never allowed to.

The momentum residual rises from 4.9e-14 to between 6.4e-13 and 3.2e-12. That is
not drift: the ligament arm carries **12 kN** of internal tension at t=0 and the
capsule arm **27 kN**, so the residual sits at about 1e-16 of the forces being
summed — the same relative floor the bare arm sits at. An absolute threshold here
would be comparing against the wrong thing.

**Cost is 1.8× the bare plant for the admissible set and 2.9–4.1× for the full
one, on the median.** The machine was shared and loaded while these ran, so the
worst-advance column is not reported as a tissue cost: the bare arm's own worst
advance was 1.04 s against a 0.073 s median in the same run.

---

## The gate the brief set — it fails on the full set and passes on the filtered one

> *with ligaments carrying load, the joint stops should become less necessary —
> measure the range excursion with ligaments and without, on the same trajectory*

2.0 s prone drop, nothing driving, the same seed configuration the withdrawn
973 mm crawl came from. Worst excursion past the model's **own declared**
coordinate ranges:

| arm | elements | worst excursion | coordinate |
|---|---:|---:|---|
| `bare` | 0 | 17.43° | hip_rotation_l |
| `stops` at 30 N·m/rad — today's plant | 0 | 5.30° | hip_rotation_l |
| `ligaments` | 105 | **32.78°** | hip_rotation_r |
| `ligaments_capsules` | 117 | 34.22° | hip_rotation_r |
| `stops_ligaments` | 105 | 10.63° | hip_rotation_l |
| `stops_ligaments_capsules` | 117 | 11.69° | hip_rotation_l |
| **`admissible`, no stops at all** | 66 | **6.08°** | hip_rotation_l |
| **`stops_admissible`** | 66 | **4.05°** | hip_rotation_l |

Two results, and they are opposite.

**The derived set as a whole makes the plant worse.** 105 ligaments roughly
double the excursion, 17.4° → 32.8°, and adding them to the stops takes 5.3° to
10.6°. Said plainly because it is the result.

**The kinematically admissible 66, with no joint stops at all, hold the plant to
6.08°** — within a degree of what the 30 N·m/rad engineering stops achieve, and
achieved by 66 tension elements whose stiffness came from the body's own declared
ligament modulus and whose attachments came from the structures' own surfaces
rather than from a stated constant. Both together give **4.05°**, better than
either. And at the end of that run **no element is past ultimate strain**: max
strain 0.136, 20 of 66 loaded, 781 N of total tension.

So the brief's gate is answered: *ligaments can carry what the joint stops are
standing in for, but only the ones whose attachment geometry survives a
kinematic check, and that is 66 of 117.*

### Why the other 51 make it worse

A line element under tension `F` makes a generalized force `−F·dl/dq`. Summing
that over the derived set gives the passive moment the ligaments make about each
coordinate, as a function of that coordinate — arithmetic, not a rollout, so it
is exact.

On the **full** set, at `hip_rotation_l` the ligament moment is **−26 N·m at
neutral**, reaching −153 N·m at +56°, and in the direction the plant actually
drifts it never exceeds **+4.4 N·m**. Those elements *drive* the coordinate out
of its range and never pull it back. At `knee_angle_r` they oppose flexion with
**−648 N·m at 90°**, which is not a knee.

The dominant elements at 90° of knee flexion say what is wrong:

| element | length | slack | strain | force |
|---|---:|---:|---:|---:|
| posterior cruciate | 68.1 mm | 37.5 | **0.82** | 9.1 kN |
| anterior cruciate | 59.1 mm | 33.4 | **0.77** | 7.5 kN |
| transverse ligament of knee | 75.2 mm | 26.9 | **1.79** | 3.7 kN |

A ligament whose ultimate strain is 17.1% is being asked for 77%. **A real
cruciate is near-isometric because its femoral footprint sits close to the
flexion axis and the ligament wraps; a straight line between two tip centroids
sits 22 mm off that axis and lengthens with the full lever arm.** The missing
machinery has a name: a wrap surface per joint. `PathWrap` and `WrapCylinder`
exist in OpenSim and this model already carries 52 wrap objects, but a
`Blankevoort1991Ligament` built from a two-point path has none, and fitting one
per joint from bone geometry is a piece of work that has not been done.

### The filter, and what it is not

An element that passes ligament ultimate strain **inside the declared range of a
joint it spans** is an attachment in the wrong place, not a ligament. Sweeping
each element over every coordinate on the chain between its two bodies, one at a
time with the rest at the reference pose, **66 of 117 are admissible** — 57
ligaments and 9 capsules. The threshold is not a knob tuned until a count looked
right: the same 66 pass at 10% strain and at 17.1%, then 70 at 25%, 80 at 40%,
95 at 60%. There is a real gap there.

It is a **screen, not a proof**, and it shows: at the stance pose five of the
admissible 66 are past ultimate strain anyway, all of them right-shoulder
elements, because the stance pose differs from the reference in `arm_rot_r` by
34° *and* `arm_flex_r` by 11° at once and the sweep moves one coordinate at a
time. Two coordinates together can be worse than either alone. That is stated in
the script and is what happened.

**The single-coordinate moment sweep understates the filtered set.** Measured
about one coordinate with everything else at the reference pose, the admissible
66 make ≤ 2.2 N·m outside the declared range at every coordinate except the
metatarsophalangeal joints, and six coordinates have no admissible element
spanning them at all — which reads as a set that can do nothing. In the prone
rollout the same 66 hold the plant better than the joint stops, because the
rollout is far from the reference pose in several coordinates at once and the
elements load there. The two measurements do not disagree; the sweep is a
one-dimensional slice of a 33-dimensional restraint, and quoting it alone would
have been the "compared against the wrong thing" error again.

### One thing that only became visible because something carried load

**The anatomy's registered reference pose is outside the model's own declared
range on four coordinates** — both knees by 5° (the model declares
`knee_angle ∈ [0, 2.44]` and the registration puts it at −0.09) and both
`pro_sup` by 4°. The ligaments hold those joints at the anatomy's pose and the
stops push them back toward the model's, so the two pull in opposite directions
there. Neither artifact was wrong on its own terms; they disagree, and nothing
had ever asked them to agree.

### How to switch them on

`NativeMechanicalStream(..., tissue_ligaments='data/derived/tissue-force-elements-v1',
tissue_ligament_admissible_only=True)`. The flag is not defaulted on: the bundle
carries all 117 elements and which subset a run wants is the caller's statement,
not a silent default. **Every measurement above says to pass it.**

## Cartilage, menisci and discs: a DATA gap, not a mechanics gap

The direction asks for cartilage as the thing that makes bone-on-bone contact
physical. Two measurements decide whether that is reachable.

**There is no articular cartilage in this body.** 51 cartilage entities: 34
costal, 8 laryngotracheal, 7 nasal, 2 triradiate growth plate. **Not one is a
joint surface.** `data/derived/joint-substrate-assessment-v1/` reached the same
zero independently against the *unregistered* extended atlas, and counted what
closing it would take: **1,012 cartilage layers over 506 bone contact pairs,
0.334 m² of patch area** — all of which would have to be synthesised, which is
inventing anatomy.

**And bone-on-bone contact is not a substitute.** A synovial joint's two bone
surfaces *overlap by construction* — the femoral head is inside the acetabulum —
and the measurement says so. Sweeping each scaffold joint's own declared range
and ray-parity-testing one body's bone vertices against the other's surfaces:

| joint | closest approach | configurations interpenetrating |
|---|---:|---|
| hip (l / r) | 1.2 / 0.9 mm | **21/21 and 19/21** |
| back (pelvis–torso) | 0.8 mm | 14/21 |
| ankle (l / r) | 1.0 mm | 6/7 each |
| elbow (l / r) | 0.5 mm | 5/7 each |
| metatarsophalangeal (l / r) | 0.4 mm | **7/7 each** |
| subtalar (l / r) | 2.0 mm | 1/1 each |
| knee (l / r) | 2.5 / 4.1 mm | 0/7 |
| patellofemoral, radioulnar, acromial, radiocarpal | 2.2 – 7.0 mm | 0 |

**11 of 21 joints interpenetrate somewhere in their own declared range.** An
`ElasticFoundationForce` between the two bones of such a joint would fire
permanently and invent a load the body never carries — which is exactly what the
engine's own comment warns a shared contact set would do.

Geometry is read from the model's **own** attached Body meshes, all 81 bones with
each Body's scale applied, *not* from `data/derived/segment-contact-meshes`,
which is the SimTK-admissible subset and drops both hip bones, both tibiae, the
skull, the jaw and the spine. The first run of this measurement used that subset
and reported an 82 mm hip "joint gap" — it was measuring the sacrum against the
femur.

Of the 69 menisci and discs, **59 are inside one rigid body**: 23 intervertebral
discs and 23 nuclei pulposi in `torso`, the interpubic disc in `pelvis`, the
acromioclavicular and sternoclavicular discs in `torso`. The four knee menisci
and the distal radio-ulnar discs are the ones that sit across a scaffold joint —
and the knee is the only joint in the body with intra-articular geometry to build
an articular layer from **and** bones that stay clear through its whole declared
range. That is the one reachable articular-surface job in this body.

---

## Fascia, aponeurosis and retinacula: blocked on the muscle path, not on a solver

A retinaculum is not a joint restraint. It is a strap that stops a tendon
bowstringing, and what it does mechanically is **change a muscle's path**.

Read out of the model the engine actually assembled — not the model on disk,
because `replacePathsWithFunctionBasedPaths` makes them different objects:

* **80 of 98 muscles run as `FunctionBasedPath`** — polynomials in the
  coordinates, with no via point, no wrap surface, no geometry a strap could
  touch. The 18 that keep a `GeometryPath` are the twelve `arm26` muscles and
  six trunk muscles (erector spinae and the obliques).
* The model declares **52 wrap objects** and 14 `PathWrap`s survive on those 18
  paths, so the wrap *mechanism* is present and serves only what is left.

**Not one of the 18 retinacula can constrain a muscle path in this plant.** Four
hold tendons of muscles the plant does not have at all — there is no hand or
wrist musculature, and `hand_l`/`hand_r` carry zero `PathPoint`s. The other
fourteen hold muscles the plant *does* have — `tibant`, `edl`, `ehl`, `perlong`,
`perbrev`, `tibpost`, `fdl`, `fhl`, `vasmed`, `vaslat`, `recfem` — and every one
of those is a polynomial. The blocker is not the retinaculum and not a missing
force class; it is the path substitution, and that substitution is what made the
plant affordable in the first place.

**Fascia is a different block again.** A deep fascial plane's mechanics is a
glide surface between muscle groups and a constraint on muscle bulging, and both
need a deformable continuum. `docs/SEGMENT_CONTACT_SURFACES.md` establishes there
is none anywhere in this stack: every compliant element is a 1-D spring law on a
rigid carrier. The one part of fascial mechanics a line element *could* carry is
longitudinal load transfer between two segments — the fascia lata and the three
layers of thoracolumbar fascia are the cases, and 25 of the 44 fascia entities do
span two segments. It is not modelled here because a fascia is a sheet whose
tension is distributed over its width, and a single line between two centroids
concentrates it: the same error that gave the hip capsule a 924 kN cross-section.

---

## Adipose and periosteum: acquisition, not modelling

**Periosteum: 0 entities.** Not thin, not badly bound — absent. So is
perichondrium, and so is bone marrow as geometry. Nothing can be given a force
element because there is nothing there.

**Adipose: 2 entities, both infrapatellar fat pads, 0.630 mL each — 1.26 mL of
adipose geometry in the whole body.** The only other structure that can carry
body fat is the `hypodermis` skin layer at 8.902 L, which at 950 kg/m³ is 8.46 kg.
There is **no visceral, no intermuscular and no intramuscular adipose geometry at
all**. `body_fat_fraction` is in the body-parameter schema and
`docs/BODY_PARAMETERS.md` already records that it *reaches nothing*; this is why.

That matters mechanically and not only cosmetically.
`docs/SEGMENT_CONTACT_SURFACES.md` measured that the declared skin alone cannot
hold the body up — 761 N over ~0.02 m² of plantar skin at 1.72 MPa/m needs
~22 mm of indentation against a 6.6 mm layer — and named the missing stiffness as
the fat and muscle between skin and bone. **This is that gap, counted.**

What would have to be acquired, stated so nobody synthesises it and calls it
measured:

* **articular cartilage** — a per-joint-surface segmentation. 1,012 layers over
  506 bone contact pairs by the existing assessment's own count. No catalogued
  source in this repository ships one; BodyParts3D and Z-Anatomy both annotate
  articular *facets* as named features of the bone, not as separate surfaces.
* **periosteum** — a bone-enveloping sheet, one per bone. Derivable as an offset
  of each bone surface, which is synthesis, and must be labelled as such.
* **adipose** — subcutaneous depth per skin region, plus visceral and
  intermuscular compartments. The hypodermis carries a single thickness for the
  whole body; regional depth is what a real fat model needs and no source here
  has it.

---

## What still carries no force, and why — the honest list

| class | entities | blocked by |
|---|---:|---|
| bursa / synovium | 80 | a bursa is a **glide interface**, and the plant has no contact between two soft surfaces. 68 of 80 are inside one rigid body anyway. |
| meniscus / disc | 69 | 59 are inside one rigid body (23 IVDs and 23 nuclei in `torso`). The rest need an intra-articular body between two bones, which is a new mobilized body, not a force element. |
| cartilage | 51 | **no articular cartilage exists.** Data gap. |
| tendon | 44 | 36 of them are tendon **sheaths** — glide interfaces, same block as bursae. The two real tendons (calcaneal, EDL) are already inside the model's musculotendon units and adding them would double-count. |
| fascia / aponeurosis | 44 | glide surface and bulge constraint need a deformable continuum; the load-transfer part needs a sheet, not a line. |
| retinaculum | 18 | 80 of 98 muscle paths are polynomials with no geometry to constrain; the other 4 straps hold muscles this plant does not have. |
| adipose | 2 | 1.26 mL of geometry. Data gap. |
| periosteum | 0 | absent entirely. Data gap. |
| **ligament (inadmissible)** | 48 of 105 | a straight-line path is not a ligament's path; needs a wrap surface per joint. |
| **joint capsule (inadmissible)** | 3 of 12 | same. |

And the part that *is* now mechanics: **105 ligaments and 12 joint capsules
carry tension between two rigid bodies, and the 66 of them that stay within their
own failure strain across a spanned joint's whole declared range hold the plant
to 6.08° past its declared ranges with no joint stops at all** — against 5.30°
for the 30 N·m/rad engineering stops and 17.43° for a bare plant. Running both
gives 4.05°.

---

## What would move this next, in order

1. **A wrap surface per scaffold joint**, fitted from the bone geometry, so a
   ligament path bends where a real one does. This is the single change that
   would turn the 48 inadmissible ligaments from wrong into useful, and those 48
   are the cruciates, the collaterals and the ankle ligaments — the ones that
   carry a real joint's restraint. The admissible 66 already do the stops' job;
   these would be the ones that do a knee's.
2. **Articular cartilage geometry.** Acquisition, not modelling. Until it exists
   there is no articular surface to make physical, and the knee is the only joint
   with intra-articular geometry to start from.
3. **Reversing the path substitution for the ankle and knee muscles**, at a
   measured wall-clock cost, so the six ankle retinacula have something to hold.
4. **Regional subcutaneous depth**, which is simultaneously the adipose gap and
   the reason the declared skin cannot hold the body up.
