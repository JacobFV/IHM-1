# the body's tissue as mechanics, and where it stops being possible

645 of the bound body's 4,000 entities are tissue structures — ligaments,
capsules, menisci, discs, cartilage, bursae, fascia, retinacula, adipose — every
one of them with real mesh geometry, and the plant carried force elements for
**none** of them. They existed as geometry and not as mechanics.

**117 of them now carry force.** What that buys, measured on three different
drops rather than one: the derived set as a whole makes the plant **worse**, and
the 66 elements that survive a kinematic check make it **no worse and sometimes
much better, but not reliably so**. Added to the joint stops they improve the
worst excursion past the model's declared ranges on every drop tested, by 4–24%.
They do **not** replace those stops, and the first version of this file said they
did on the strength of one trajectory — the withdrawal is below and is the point.

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

## The gate the brief set — and the claim it took three drops to get right

> *with ligaments carrying load, the joint stops should become less necessary —
> measure the range excursion with ligaments and without, on the same trajectory*

2.0 s prone drop, nothing driving but a PD hold on the declared torque ports.
Worst excursion past the model's **own declared** coordinate ranges. Three
different initial conditions, because a rollout is deterministic and repeating it
proves nothing, but running the same arms from different drops is what says an
ordering is a property of the force set rather than of one drop.

| arm | elements | `prone` | `prone_high` | `prone_rolled` | mean |
|---|---:|---:|---:|---:|---:|
| `bare` | 0 | 17.43° | 30.69° | 30.27° | 26.13° |
| `stops` at 30 N·m/rad — today's plant | 0 | 5.30° | 6.62° | 5.96° | **5.96°** |
| `ligaments` | 105 | 32.78° | 35.33° | 32.88° | 33.66° |
| `ligaments_capsules` | 117 | 34.22° | — | — | — |
| `stops_ligaments` | 105 | 10.63° | — | — | — |
| `admissible` | 66 | 6.08° | 30.47° | 18.50° | 18.35° |
| `stops_admissible` | 66 | **4.05°** | **6.39°** | **5.35°** | **5.27°** |

Three results, and the middle one is a withdrawal.

**1. The derived set as a whole makes the plant worse, on every drop.** 105
ligaments take 17.4° → 32.8°, 30.7° → 35.3°, 30.3° → 32.9°. Adding them to the
stops takes 5.3° → 10.6°. Consistent in sign three times out of three.

**2. The kinematically admissible 66 do NOT replace the joint stops.** On the
`prone` drop they hold 6.08° with no stops at all, against the stops' 5.30° and a
bare plant's 17.43°, and *this file's first version reported that as ligaments
replacing the stated engineering constant.* It does not survive a second drop:
on `prone_high` the same 66 give **30.47° against a bare plant's 30.69°** — no
restraint at all — and on `prone_rolled` 18.50° against 30.27°, a real but partial
one. **The 6.08° was a lucky draw**, which is the failure mode
`docs/LOG.md` in IBM-1 records over and over; one trajectory is one trajectory
even when nothing in it is random.

What survives is weaker and true: **the admissible set is never worse than the
bare plant on any drop**, where the unfiltered set is always worse. The filter
buys safety, not restraint.

**3. Added to the stops, the admissible set helps on every drop.** 5.30 → 4.05,
6.62 → 6.39, 5.96 → 5.35. Between 4% and 24%, mean 5.96° → 5.27°, and the sign is
the same three times out of three. That is the defensible positive: the body's
own tissue, with attachments from its own surfaces and stiffness from its own
declared modulus, takes a little of the load off a stated engineering constant.
It does not take it over.

And at the end of the `prone` admissible run **no element is past ultimate
strain**: max strain 0.136, 20 of 66 loaded, 781 N of total tension.

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

**The single-coordinate moment sweep and the rollouts do not disagree, and both
are one-sided.** Measured about one coordinate with everything else at the
reference pose, the admissible 66 make ≤ 2.2 N·m outside the declared range at
every coordinate except the metatarsophalangeal joints, and six coordinates have
no admissible element spanning them at all — a set that on that slice can do
almost nothing. Three drops later that is the better predictor of the mean
behaviour (18.35° against a bare 26.13°) than the one drop where they held.

### One thing that only became visible because something carried load

**The anatomy's registered reference pose is outside the model's own declared
range on four coordinates** — both knees by 5° (the model declares
`knee_angle ∈ [0, 2.44]` and the registration puts it at −0.09) and both
`pro_sup` by 4°. The ligaments hold those joints at the anatomy's pose and the
stops push them back toward the model's, so the two pull in opposite directions
there. Neither artifact was wrong on its own terms; they disagree, and nothing
had ever asked them to agree.

### Rendered

`ibm1-progress-archive/video/tissue_prone_drop_three_arms.mp4` — the same 2 s
prone drop three times on 429 bone surfaces of the real anatomy: bare, with all
105 ligaments, and with the joint stops plus the admissible 66. The centre panel
is a negative result rendered on purpose. The anatomy is posed from the plant, as
always; nothing there is a body moving itself.

### How to switch them on

`NativeMechanicalStream(..., tissue_ligaments='data/derived/tissue-force-elements-v1',
tissue_ligament_admissible_only=True)`, **with the joint stops still installed**.
The flag is not defaulted on: the bundle carries all 117 elements and which
subset a run wants is the caller's statement, not a silent default. Every
measurement above says to pass it, and every measurement above says not to drop
the stops.

---

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
| adipose | 2 | 1.26 mL of geometry. Data gap -- see "Where adipose can come from" below. |
| periosteum | 0 | absent entirely. Data gap. |
| **ligament (inadmissible)** | 48 of 105 | a straight-line path is not a ligament's path; needs a wrap surface per joint. |
| **joint capsule (inadmissible)** | 3 of 12 | same. |

And the part that *is* now mechanics: **105 ligaments and 12 joint capsules
carry tension between two rigid bodies.** The 66 that stay within their own
failure strain across a spanned joint's declared range are never worse than the
bare plant on any of three drops, and added to the joint stops they improve the
worst excursion on all three — 5.96° mean to 5.27°. They do not replace the
stops; a first reading that said they did came from one drop and is withdrawn
above.

---

## The 51 are a REGISTRATION failure first, and a path failure second

The section above blames the straight path and prescribes a wrap surface per
joint. Two measurements say that was the wrong order.

**Moving attachments does not fix them**
(`scripts/measure_ligament_attachment_vs_path.py`, gate: recomputed tip peaks
match the stored ones to 1.6e-14). Even the most isometric fibre inside each
ligament's own footprint rescues only **13 of 51**, and only 4 of those are taut
over more than a quarter of the joint's range, so almost all of the rescued ones
restrain nothing.

**Correcting the joint's registration does**
(`scripts/measure_ligament_registration_probe.py`). The global atlas->scaffold
similarity displaces the atlas's articular centres from the scaffold's by
**20-26 mm** at knee, ankle and hip, measured by the same closest-approach rule on
both skeletons (`scripts/measure_joint_centre_registration.py`; the same rule on
the scaffold's own bones is the control). A cruciate is 33-37 mm long, so it was
sweeping about a centre displaced by most of its own length. Translate each
failing element rigidly by its joint's displacement and re-derive its attachments
exactly as the builder does (gate: a zero translation reproduces the stored peaks
to 7.8e-15):

| | before | after one translation per joint |
|---|---:|---:|
| admissible, of the 37 crossing a measured joint | 0 | **19** |
| peak strain fell | | 33 / 37 |
| median peak strain | 0.579 | **0.170** |
| anterior cruciate | 1.216 / 1.226 | **-0.022 / -0.033** |
| posterior cruciate | 1.184 / 1.191 | 0.170 / 0.326 |
| knee collaterals and capsules | 0.44-0.78 | **-0.02 to 0.00** |

A translation is the crudest correction there is, and it is not uniformly right:
the two intra-articular hip structures get WORSE (ligament of head of femur 0.73
-> 1.15, transverse acetabular 0.26 -> 0.58), the anterior talofibulars stay at
1.1, and 14 of the 51 cross joints not measured here (wrist, radio-ulnar,
iliolumbar, MTP).

### Per-segment registration, and the control that says it is not an artefact

`scripts/fit_segment_registration.py` fits one similarity per segment and carries
each ligament END through its own bone's map. On the 51 elements v1 rejects:

| registration | admissible of 51 | median peak strain |
|---|---:|---:|
| global map (v1) | 0 | 0.521 |
| one translation per joint | 19 of 37 covered | 0.170 |
| per segment, free rotation | 30 | 0.122 |
| **per segment, twist removed** | **29** | **0.136** |
| per segment, rotation held at the global map's | 28 | 0.143 |

(gates a-d pass in every mode; gate c reproduces the stored v1 strains to 7.8e-15.)

The free fit's rotations are suspect on the long bones -- femur 19.8 deg relative
to the global map of which 19.8 is twist about its own axis, radius 32-36 deg,
all spin -- because point-to-point ICP cannot determine an elongated bone's spin,
and a spin moves a ligament's attachment round the bone. **The controls say the
rescue does not ride on it**: removing the twist costs one element and holding
the rotation fixed entirely costs two. What rescues these ligaments is each bone's
per-segment position and scale, i.e. putting the joints where the scaffold's are.

`tissue-force-elements-v2` is built with
`--registration per-segment --segment-registration .../anatomy-segment-registration-rot-no-twist/registration.json`
(gate: `--registration global` reproduces v1's `ligaments.json` exactly):

| | v1 | v2 |
|---|---:|---:|
| tensile elements admissible, of 117 | 66 | **95** |
| gained / lost | | 29 / **0** |
| cruciates | 0 / 4 | **4 / 4**, peak 0.143 |
| collaterals | | **16 / 16**, peak 0.123 |

Still failing in v2: the anterior talofibulars, the ligament of the head of the
femur, the transverse ligament of the knee, the iliolumbars, and wrist and midfoot
elements that cross joints this scaffold does not have. Those are the ones a wrap
surface or a new joint has to answer, and they are now a short list.

**v2 is not yet the plant's.** v1 was judged on three drops against the joint
stops; v2 has to be judged the same way before it replaces it.

#### Direction, not magnitude

`scripts/measure_ligament_moments.py --tissue ...` gives each coordinate's passive
ligament moment outside its declared range. Its summary line is a MAGNITUDE, and
v1's pathology was a DIRECTION -- ligaments driving joints out of range instead of
pulling them back. So every out-of-range sample is classed restoring or driving.

On the **full** element set, v2 fixes the knee (v1 drove it at up to 101 N.m; v2
restores on all 36 samples) and appears to introduce a hip-flexion drive, 0 / 18
samples on the left and up to 37 N.m on the right.
`scripts/measure_hip_flexion_drive.py` (gate: per-element sums reproduce the moment
check to 4.3e-4 N.m) puts that drive on **one element per side, the ligament of the
head of the femur**, carrying 2.4-7.1 kN in v2 and 12.9-14.7 kN in v1. It is
inadmissible in both builds: an intra-articular ligament ~30 mm long, slack most of
the time in a real hip, crossing a hip centre still 11-15 mm misregistered. No
straight line from registered ends represents it.

The plant only ever receives the **admissible** set, so that is the comparison
that decides (66 elements in v1, 95 in v2):

| coordinate | v1 restoring / driving (worst drive), samples past the stop | v2 |
|---|---|---|
| knee | 5 / 11 (0.3 N.m), 1 | **18 / 0, 13-14** |
| ankle | 0 / 0 -- no restraint | 12 / 12 (0.1 N.m), **12** |
| forearm pronation | 11-13 / 0, 2-4 | 20 / 2, **11** |
| hip rotation | 26 / 0, 13 | 26 / 0, 13 |
| hip flexion | ~0 | ~0 (the femoral-head drive is excluded) |
| hip adduction | 5-6 / 16-18 (1.6 N.m) | 4-5 / 18-21 (**2.6 N.m**) |
| **toe MTP** | 30 / 0, **15** | 15-16 / 13-15, **0-1** |
| elbow flexion R | 11 / 0, 2 | 0 / 0 |
| **total** | 169 / 70 / 77 | **192 / 107 / 109** |

Net: v2's admissible set restrains the knee harder than the joint stop, adds an
ankle restraint where v1 had none, and improves pronation. It LOSES the toe
restraint -- consistent with the 17 deg forefoot swing the toe registration took
(docs/SEGMENT_CONTACT_SURFACES.md) -- and the right elbow's, and its largest drive
(hip adduction) grows from 1.6 to 2.6 N.m, still small against a 30 N.m/rad stop.
That is arithmetic at a held pose; the drops are the dynamic test.

#### The drops: v2 halves the worst excursion when added to the stops

`scripts/measure_tissue_mechanics.py --tissue data/derived/tissue-force-elements-v2`,
the same three drops that judged v1. Gate: the arms that carry no tissue at all
reproduce v1's recorded excursions on every drop to 0.01 deg (bare 17.43 / 30.69 /
30.27, stops 5.30 / 6.62 / 5.96).

Worst excursion past the declared ranges, deg:

| arm | prone | prone_high | prone_rolled | mean |
|---|---:|---:|---:|---:|
| bare | 17.43 | 30.69 | 30.27 | 26.13 |
| stops | 5.30 | 6.62 | 5.96 | 5.96 |
| v1 admissible | 6.08 | 30.47 | 18.50 | 18.35 |
| **v2 admissible** | **5.15** | **17.86** | **7.65** | **10.22** |
| v1 stops + admissible | 4.07 | 6.39 | 5.35 | 5.27 |
| **v2 stops + admissible** | **2.12** | **2.91** | **3.36** | **2.80** |
| v2 full set | 31.41 | 33.30 | 31.90 | 32.20 |

(v1's prone stops + admissible is not recorded separately; 4.07 is the documented
5.27 mean less the two recorded drops.)

* **Added to the stops, v2's admissible ligaments halve the worst excursion on
  every drop**, 5.96 -> 2.80 deg mean, where v1's took it 5.96 -> 5.27.
* **Alone, v2 does not replace the stops** -- 17.86 against 6.62 on the high drop
  -- though it roughly halves what v1 did alone (18.35 -> 10.22 mean).
* **The full set is worse than no tissue on every drop** (32.20 against 26.13),
  because it carries the inadmissible elements, including the femoral-head ligament
  that drives hip flexion. Only the admissible set is ever a candidate for the plant.

So v2's admissible set is the recommended tissue set for the participant plant,
added to the joint stops, not in place of them. The measurement scripts still
default to v1 so every result already recorded reproduces; switching the plant's
default is a separate, deliberate step.

#### ...and in the crawl it is neutral

`scripts/crawl.py --tissue DIR` loads a set's admissible elements into the plant alongside
the stops. Search B's best crawl parameters replayed for 16 s, three ways (gate: the
stops-only replay reproduces the recorded crawl, 0.9200 m and 0.2003 rad at the right
ankle -- PASS):

| arm | reaches 16 s | travel | worst excursion | peak fibre velocity (of 10) | wall clock |
|---|---|---:|---:|---:|---:|
| stops | yes | 0.9200 m | 0.2003 rad, ankle_angle_r | 9.69 | 707 s |
| stops + v1 admissible | **no -- 0.01 s** | 0 | -- | 5.20 | 8 s |
| stops + v2 admissible | yes | 0.9146 m | 0.1987 rad, ankle_angle_r | **9.99** | 863 s |

* **v1 cannot start the crawl.** A joint-speed guard stopped it on the first advance, with
  arm_rot_r past 25 rad/s.
* **v2 does no harm and no good here.** Travel and worst excursion are within 1% of the
  stops alone. The drop tests' halving was unactuated falls; in actuated crawling the stops
  already hold the worst excursion to 0.2 rad, at an ankle v2 barely restrains.
* **Its costs:** 22% more wall clock, and peak fibre velocity at 9.99 of the muscle model's
  10 -- at the edge of the force-velocity domain, where the stops-only crawl had 3% to spare.

One trajectory per arm: v2 is safe to carry in locomotion; its benefit is not shown there.

## What would move this next, in order

1. **Per-segment registration of the atlas onto the scaffold**, before any wrap
   surface. One translation per joint already rescues 19 of 37 and brings both
   cruciates inside their ultimate strain; a proper per-segment map is the same
   fix the skin needs (`docs/SEGMENT_CONTACT_SURFACES.md`). Only what survives
   THAT should be judged against the item below.
2. **A wrap surface per scaffold joint**, fitted from the bone geometry, so a
   ligament path bends where a real one does. This is the single change that
   would turn the 51 inadmissible elements from wrong into useful, and those 51
   are the cruciates, the collaterals and the ankle ligaments — precisely the
   ones that carry a real joint's restraint. The 66 that survive the filter
   survive it largely because they barely change length with the joint, which is
   the same thing as saying they have little to restrain.
3. **Articular cartilage geometry.** Acquisition, not modelling. Until it exists
   there is no articular surface to make physical, and the knee is the only joint
   with intra-articular geometry to start from.
4. **Reversing the path substitution for the ankle and knee muscles**, at a
   measured wall-clock cost, so the six ankle retinacula have something to hold.
4. **Regional subcutaneous depth**, which is simultaneously the adipose gap and
   the reason the declared skin cannot hold the body up.


## Where adipose can come from (2026-09-10)

This body carries two adipose entities, 1.26 mL, and the programme requires contact to be
mediated by fat, muscle and skin rather than bone. Two routes were checked:

* **TotalSegmentator `tissue_types` / `tissue_4_types`** (subcutaneous_fat, torso_fat,
  skeletal_muscle, intermuscular_fat): **licensed**, listed under "Available with a license"
  in the tool's README. A licence is free for non-commercial use but is obtained through a web
  registration form (`backend.totalsegmentator.com/license-academic/`) that only the programme
  owner can submit; the key is then set with `totalseg_set_license -l <key>`. Not obtained.
* **CT attenuation, no model and no licence.** Fat is defined physically on CT: -190 to -30 HU
  is the standard adipose window, and it is the threshold TotalSegmentator's own
  `tissue_4_types` post-processing uses. Applied inside the Apache-2.0 `body` mask, it gives an
  adipose segmentation for the CT subjects already on disk (s0790, s1159, s0970), split into
  subcutaneous and deep by distance from the skin. Tested on the three CT subjects that pass
  every extraction gate:

    subject   body L  fat L  fat % breast mL breast fat % glandular %     T6    T12  T12-120  T12-250   (subcutaneous depth, 95th pct, mm)
    s0790       37.9  12.74   33.6      1161         82.6        13.8   36.5   18.1     51.0     37.3
    s1159       34.3   7.35   21.4       956         67.0        29.7   16.8   10.6     34.0     28.1
    s0970       35.3  10.58   30.0      1437         87.8        11.5   43.5   14.8     31.7      nan

  Two checks, neither a known answer in the strict sense, both passed. The breast is fat plus
  glandular tissue, and glandular tissue (~30-50 HU) falls OUTSIDE the window, so the share of
  breast voxels in the window is the breast's fat fraction; it lands in or just above the range a
  rough prior expects for these ages, with the glandular remainder in the -30..100 HU band. And
  the skin-connected fat is deepest over the lower abdomen and hips, thinner over the upper
  abdomen -- the female, lower-body distribution. Limitation: "subcutaneous" is fat connected to
  the outer 3 mm of the body, so where the abdominal wall thins some deep fat can join it; the
  lower-abdominal depths may include deep fat. This gives a fat layer for three female subjects,
  not for this body; carrying one onto this body's skin is the skin-envelope problem.


## Articular cartilage: the knee, and its gates fixed before any data is seen (2026-09-10)

This body has no articular cartilage. OAIZIB-CM (cc-by-nc-4.0; 507 OAI knee MRIs; labels femur,
femoral cartilage, tibia, medial and lateral tibial cartilage) supplies it for the knee. Only the
label archives are taken: the femur and tibia carry the registration, the cartilage is the
geometry. The registration follows the pelvic one -- a knee MRI crops both bones, so one-way ICP
pulls the MRI's partial bone surfaces onto this body's COMPLETE femur and tibia, never the reverse.

Gates, committed before a single subject is registered:
* **known answer:** this body's own femur and tibia, truncated the way a knee MRI field of view
  crops them, moved by a known similarity, are recovered within 2 mm, 2 deg and 1%.
* **which knee:** OAI images one knee per scan. Each subject is fitted to BOTH of this body's
  knees; the lower one-way residual must win by at least 20%, or the side is undetermined and
  the subject is excluded.
* **placement:** at least 95% of the mapped femoral cartilage lies within 3 mm of this body's
  femoral surface and at most 1% inside the femur; the same for tibial cartilage and the tibia.
* **joint space:** mapped femoral and tibial cartilage overlap by at most 1% of the smaller
  volume at this body's reference pose.
Reported, not judged: cartilage thickness (OAI knees are older and often osteoarthritic, so
thin or eroded cartilage is expected and is the cohort's, not the registration's).