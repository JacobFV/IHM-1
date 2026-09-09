# the body as a parameter, and where the parametrization stops

Ask for a tall person and get one. Ask for a female person and get told, with
numbers, exactly what is missing.

    ihm/body_parameters.py                        the schema
    ihm/body_scaling.py                           the exponent table: no exponent is typed
    ihm/native/model_scaling.py                   measure and scale the mechanical SCAFFOLD
    scripts/describe_body_parameters.py           print the schema; resolve a request
    scripts/materialize_body_variant.py           build a whole body: scaffold AND body
    scripts/materialize_stature_variant.py        the scaffold half of that
    scripts/scale_anatomical_body.py              the 4,000 anatomical entities
    scripts/scale_nerve_conduction.py             146 routes, 1,743 conduction delays
    scripts/scale_skin_patches.py                 1,326 patches; count fixed, area s**2
    scripts/scale_muscle_and_tissue.py            PCSA, force, mass, 117 ligaments
    scripts/measure_stature_allometry.py          is geometric similarity even true
    scripts/verify_stature_scaling.py             gate it, including on purpose broken
    scripts/audit_sex_specific_anatomy.py         what sex anatomy exists, measured
    scripts/index_anthropometry.py                sex-stratified proportions from NHANES
    scripts/measure_anisotropic_pelvis_error.py   why proportions are not a knob

`data/derived/` is gitignored throughout. The scripts are the artifacts.

**Section 4 is the one to read if you only read one.** Sections 1-3 describe the
parametrization as it was when it reached only the scaffold. Section 4 is where
it reaches the body, and it is also where the single most important fact about
the whole thing lives: the native engine renormalises the scaled mass straight
back out.

## 1. The schema

Seven parameters, and the status field is the load-bearing part:

| parameter | unit | range | status |
|---|---|---|---|
| `stature_m` | m | 1.40 – 2.05 | **implemented** by this work |
| `mass_kg` | kg | 35 – 160 | **surfaced**; was already a knob |
| `muscle_force_scale` | – | 0.25 – 4.0 | **implemented**; defaults to `stature_scale²` |
| `environment` | – | free / supine / upright | surfaced; scenario, not body |
| `sex` | – | **{male}** | **declared**; reaches nothing |
| `age_years` | y | 18 – 90 | declared; reaches nothing |
| `body_fat_fraction` | – | 0.05 – 0.50 | declared; reaches nothing |

*surfaced* means the knob existed and the schema only names, bounds and
attributes it. *implemented* means this work made it one. *declared* means it is
written down and reaches the running body through nothing at all, and the
schema says which of the three every parameter is rather than presenting a flat
list in which they look alike.

`sex` has a domain of one member and `resolve({'sex': 'female'})` raises. That
is deliberate. Accepting the value would change a string in a JSON file and
nothing else.

### What surfacing found

`target_mass_kg` was a knob before this work. It was also a bare literal —
**77.6122029** — written out by hand in **67 places across 58 files**, and no
derivation for it exists anywhere in the repository. Four candidates were
checked and rejected:

| candidate | value | why not |
|---|---|---|
| the source model's own mass | 85.26984854 kg | the engine divides *into* this |
| the anatomical body's composed mass | 70.7713 kg | different body, 9.7% apart |
| BioGears StandardMale | 77.1107029 kg | superseded, and not equal |
| mean vertical GRF over the trial ÷ g | 60.83 kg | a duty-cycle artefact, not a weight |

It is most likely the AddBiomechanics subject's measured mass carried in by
hand. It is recorded as unattributed rather than guessed at.

### Two bodies, and they are not the same size

|  | mechanical | anatomical |
|---|---|---|
| what it is | OpenSim Rajagopal subject | BodyParts3D atlas + BioGears |
| mass | 77.6122029 kg (asked for) | 70.7713 kg (composed) |
| stature | 1.7972725 m (measured, below) | 1.7194712 m (skin extent) |

They disagree by **9.7% in mass**. `docs/DISCONNECTS.md` item 1 is the same seam
seen from the geometry side. The schema holds both rather than picking one, and
tags every parameter with which body it addresses.

### Measuring the stature the model does not record

An `.osim` file records no height. `head_marker_height_m` measures it as the
vertical distance, at the model's default pose, from the plane of the
AddBiomechanics floor-level virtual markers (`R/L.HeelGround`, `MT5Ground`,
`ToeGround`) to the `Head` marker. Both are properties of *this subject*; a mesh
bounding box would be a property of the generic donor mesh refracted through
nine anisotropic per-body scale triples, and the head rides `torso` with no neck
joint.

**1.7972725 m**, and it has a check with an answer fixed in advance. Carried
onto the atlas by the independently fitted registration scale 0.963 it predicts
**1.7308 m** against the atlas skin extent of **1.7195 m** — **0.66%**. Two
measurements sharing no input agree to within a centimetre, which also bounds
from above how far the `Head` marker sits above the vertex.

The forward-kinematic walk evaluates each joint's transform functions at its
coordinates' declared defaults. The first draft assumed the default pose put
every joint at identity; `patellofemoral_r` sits at 1.14e-3 rad with the knee at
zero, and the assumption was wrong by 7.6e-7 m.

### The ranges are measured now

They were engineering guesses in the first draft of this schema. `BMX_J.parquet`
— standing height, weight, upper leg length, upper arm length, waist and hip for
8,704 NHANES participants — has been on disk since the population index was
built, and **nothing has ever read it**. `scripts/index_nhanes.py` loads the
table into `tables` and then uses only its 23 laboratory and vital bindings.

Survey-weighted over 4,883 adults aged 20–79, 0.5th–99.5th percentile:

| | stature | weight |
|---|---|---|
| men | 1.560 – 1.955 m | 52.3 – 169.8 kg |
| women | 1.437 – 1.801 m | 43.2 – 162.3 kg |
| both | 1.437 – 1.929 m | 45.1 – 164.5 kg |

The weighting is checked against an answer known before the code ran — the
published CDC/NCHS means for adults 20+ — and lands on **175.4 cm against 175.3**
for men and **161.3 against 161.3** for women.

## 2. Stature, and what scaling actually changes

`materialize_stature_variant.py --stature-m 1.90` writes a variant directory with
a scaled model, a scaled fitted-path set, a scaled contact geometry set, a scaled
muscle catalog and a `registration.json` the native stream accepts.

**Scaled**, by one isotropic linear factor `s`:

- every `PhysicalOffsetFrame` translation — the joint frames, so segment lengths;
- every `PathPoint` location, every `Marker` location;
- every wrap object's radius, length, translation and ellipsoid dimensions;
- every `Mesh` and `FrameGeometry` scale factor, so the display geometry follows;
- body `mass_center`; body `mass` by `s³` and `inertia` by `s⁵`, which keeps the
  mass/inertia ratio right whatever the native engine's own uniform
  `target_mass_kg` renormalisation then does to both;
- the **translational** transform axes of the `CustomJoint`s — the knee's three
  translation functions are lengths and move with the skeleton; its rotations are
  angles and must not;
- the translational coordinates' `default_value` and clamp `range`;
- every muscle's `optimal_fiber_length` and `tendon_slack_length`;
- **every coefficient of all 80 fitted path polynomials**;
- the 12 foot contact spheres' locations and radii;
- `max_isometric_force` by `muscle_force_scale`, and `optimal_force` on the
  reserve torque actuators by `muscle_force_scale · s`.

**Not scaled, by design**, and stated in every variant's registration:

- rotational coordinate ranges — angles are dimensionless;
- `ExpressionBasedCoordinateForceSet`, the passive shoulder/elbow/hip stops.
  Under geometric similarity these moments should go as `s³`, but they are
  authored range-of-motion stops rather than measured tissue properties and
  nothing in this repository fixes their scaling. **A scaled body has the source
  subject's passive joint stops.** This is a known gap, not a decision;
- `ContactForceSet` Hunt–Crossley stiffness and dissipation — a contact-model
  question, not a geometry one;
- the muscle curve shapes (`ActiveForceLengthCurve`, `ForceVelocityCurve`,
  `FiberForceLengthCurve`, `TendonForceLengthCurve`), `pennation_angle_at_optimal`
  and `max_contraction_velocity` — all already in normalised or angular units and
  correctly invariant.

### The trap, and the mechanism that closes it

**The engine does not run the model's muscle paths.**
`native_mechanical_stream.cpp` calls
`ModelFactory::replacePathsWithFunctionBasedPaths` and substitutes all 80 with
fitted multivariate polynomials from a *separate file*, and that file was
hard-wired to the source subject's. A scaled model with it is a body whose
skeleton was resized and whose muscles were not — silently.

`augmented_registration` therefore now carries `geometric_scale` and
`source_overrides` together, and `NativeMechanicalStream` **refuses** the first
without the second.

Scaling the polynomials is exact, not approximate, and the scaler checks the fact
that makes it so: all 80 take **rotational coordinates only**. Path length is then
homogeneous of degree one in the geometry at fixed pose, so `L'(q) = s·L(q)`
identically in `q`, and multiplying every coefficient by `s` realises exactly
that. The moment arms follow, because OpenSim takes them as `−dL/dq` when no
explicit `moment_arm_functions` are given — none are here — and a moment arm is a
length. If any path took a translational coordinate the scaler raises instead of
producing a plausible wrong number.

### The gates, and the one that would not have caught it

`verify_stature_scaling.py` runs four arms. The one that matters is **sabotage**:
scale the model to 1.13×, leave the fitted path set at the source subject's.

| gate | on the sabotaged body |
|---|---|
| stature scales | **PASSES** — 2.030918 m, exactly right |
| path length / optimal fibre length in 0.2–20 | **PASSES** — drifted 11.5%, stayed in [0.922, 8.412] |
| every fitted path length scales | catches it |
| every moment arm scales | catches it |
| polynomial length / path-point polyline length does not move | catches it |

**The 0.2–20 band gate is not sufficient.** It is the gate this work was told to
use, and on its own it would have passed a 13% mis-scaled body. It is a corpus
quality check — it catches a motion file holding accelerations under coordinate
names, where excursions come out at 1e22 — and it was never a scaling check.
Stature alone is no better: exactly right on a body whose every lower-limb muscle
was wrong.

Of the three that catch it, two are consistency checks: they read back what the
scaler wrote through independent code, which catches a file scaled by the wrong
factor or not at all but cannot catch both files being scaled consistently and
both being wrong. The **polyline** gate can. It compares the fitted polynomial
length against a straight-line walk of the model's own `PathPoint` locations — a
different file, different elements, different code path — and requires their
ratio not to move. Before any scaling the two representations of the same 80
musculotendon paths agree to within **5.1%** at the reference pose, which is the
wrap geometry the polyline does not model.

The other arms: identity at `s = 1.0` reproduces the base to relative error
**0.0** in both stature and every path length; five scales from 0.80 to 1.13 are
monotone in stature and in soleus path length; and the refusal arm confirms the
stream raises on a scaled registration with no scaled path set.

### What stature scaling is not

Isotropic. One factor everywhere. It changes how **big** the subject is and not
what **shape** it is: segment length ratios, segment mass fractions and radii of
gyration are exactly the source subject's at every stature. And no native run has
been accepted on a scaled body — `native_acceptance_complete` is `false` in every
variant. The gates above are measured on XML, not on a simulation.

## 3. Sex

### What the body actually has

`audit_sex_specific_anatomy.py`, measured over three entity sets:

| | simulated body | display atlas |
|---|---|---|
| entities | 4,000 | 8,979 |
| **male-specific** | **32** | **78** |
| **female-specific** | **0** | **0** |
| mammary gland / nipple / areola | **0** | **0** |
| surface region only ("mammary region") | 4 | 8 |
| excluded as not sex-specific | 24 | 48 |

**The body is not sexless.** It carries a complete male genital tract: the three
penile bodies, both testes, epididymides, deferent ducts, seminal vesicles and
ejaculatory ducts, the prostate, and their arteries and veins — all bound to the
`pelvis` segment at coherence ≥ 0.578, so they move with the body like every other
entity.

An earlier search reported zero because it was run against
`bodyparts3d_index.json`, which **has no name field**. Its searchable text is FMA
concept labels, so a name search over it returns 0 for every term including
`femur`.

Every match is classified rather than summed, because the raw regex overcounts by
48 in the display atlas alone: cavernous sinuses are in the skull and the labial
arteries are of the lips. The one that mattered was **`salpingopharyngeus`** — the
*salpinx* there is the auditory tube, not the uterine tube. Taking the regex at
its word would have reported a body with two fallopian tubes in its throat and no
uterus.

**There is no breast of either kind.** Not a female one and not a male one. What
exists is Z-Anatomy's `mammary region` and `inframammary region`: named patches of
chest skin over a chest wall carrying no gland.

### Does any catalogued source ship a female variant

**No.** Surveyed across every source card and anatomy doc in the repository:

| source | what it is | sex |
|---|---|---|
| BodyParts3D 4.0 | the whole-body atlas | male only, stated in `docs/CANONICAL_ANATOMY.md` |
| Z-Anatomy | the extended surfaces | male; its card sets `independent_subject: false`, and it partly derives from BodyParts3D |
| OpenSim models (Rajagopal, gait2392, Arm26) | the mechanics | male 50th-percentile |
| MOBL_ARMS shoulder | shoulder complement | "50th-percentile male model" |
| HiPCT kidney | microvascular prior | male donor, 63 y |
| Vascular Model Repository | CFD cases | mixed; two female subjects, but vessel trees, not a body |
| NHANES | population | both sexes, and it is a survey — no geometry |
| Visible Human Project | — | **not catalogued at all**; the only mentions in the tree are inside retained third-party HTML, and both refer to the male dataset |

So a female body is **new source geometry**, not new code. This is the fact that
determines everything below it.

### The proportional half, which is measurable and is measured

`index_anthropometry.py`, survey-weighted, 4,883 adults 20–79, male minus female,
`d` the weighted standardised difference:

| | male | female | d |
|---|---|---|---|
| standing height (cm) | 175.37 | 161.33 | **1.98** |
| weight (kg) | 91.87 | 78.23 | 0.62 |
| BMI (kg/m²) | 29.78 | 30.03 | **−0.03** |
| upper leg / stature | 0.2376 | 0.2316 | 0.44 |
| upper arm / stature | 0.2244 | 0.2230 | 0.14 |
| hip circ / stature | 0.6041 | 0.6805 | **−0.91** |
| waist / hip | 0.9730 | 0.8946 | **1.07** |
| waist / stature | 0.5893 | 0.6107 | −0.21 |

The shape of this is not what a list of "sexually dimorphic proportions" suggests.
Sexual dimorphism in this sample is **overwhelmingly size** — `d = 1.98` in stature
and `d = −0.03` in BMI — then **pelvic breadth relative to stature** and waist:hip,
and limb-to-trunk proportion differs by under 3%.

The size term is exactly what `stature_m` already carries. The pelvic term is what
it cannot.

### Why the pelvic term is not a knob

`measure_anisotropic_pelvis_error.py`. Widen the pelvis medio-laterally by the
measured female/male hip-over-stature ratio **1.1265**, hold everything else, and
measure how far each muscle's path length moves:

- **56 of 98** paths change; **42 do not**;
- among those that change the span is **0.01% to 6.23%** — `addmagProx` 6.2%,
  piriformis 4.8%, `addbrev` 4.0%, `glmax3` 2.0%, and nothing below the knee.

A polynomial coefficient scale admits **one** factor. Any single factor is wrong
for every muscle but at most one. Anisotropic proportions therefore need the
fitted paths **refitted**, not rescaled — and refitting means OpenSim's
`PolynomialPathFitter`, sampling the real `GeometryPath` with its wrap objects
over the coordinate ranges. The OpenSim Python bindings are not installed here,
and the native adapter that is built exposes the stream protocol, not the fitter.

That is the honest boundary of stage 2, and it was found by measurement rather
than assumed.

### The hormone cluster

`ihm/native/reproductive.py` executes the hash-pinned Schlosser & Selgrade
Physiome model. Four ODE states — releasable and circulating LH, releasable and
circulating FSH. Estradiol, progesterone and inhibin are **prescribed analytic
functions of time**, not states, and its own limitation list says so: *"not an
autonomous ovarian cycle generator"*, *"extending its nonperiodic inputs does not
generate repeated menstrual cycles"*.

The sharper point is one its limitation list does not make. **It is a menstrual
cycle model, and it is running beside a body with testes and a prostate.** It is
also coupled to nothing: `ihm/assembly/body.py` groups it with the retained domain
materializations, and it is reachable only as a standalone trajectory through
`/api/reproductive`.

Coupling it to anatomy is not blocked by the prescribed inputs. It is blocked by
the same thing as everything else in this section: closing the loop means E2 and
P4 come *from* an ovary, and there is no ovary. On the body that exists, the
coupling that could be built is the male one — a testis is present, and a
testosterone axis is not the model that is implemented.

Separately, `ihm/fields/systems.py` declares a symbolic scaffold with
`reproductive.estradiol / progesterone / testosterone / gamete_production`,
`uterine.endometrial_thickness`, `uterine.contractile_pressure` and a `placental`
group. These are prior centre/spread placeholders in the ontology with **no
anatomical entity to attach to** — `uterine` on a body with no uterus.

### Hair

`ihm/assembly/hair_fields.py` partitions 244 Terminologia Anatomica surface
regions into **21 hair fields** over the real skin mesh, each with follicle
density, shaft radius and length, growth rate, anagen fraction and shaft-bearing
fraction; the built candidate holds 652,101.75 follicles over 1.7813 m², 20 of 21
fields resolved.

Four fields are labelled `terminal_androgen_dependent` (beard, axillary, pubic,
perineal) and three `mixed_androgen_dependent` (chest, abdomen, back). **Those
labels are strings in a metadata tuple.** There is no androgen variable, no input
from `reproductive`, and no coupling of any kind — and `reproductive` does not
model testosterone anyway.

So the structure for a sex-dependent hair distribution exists and is well shaped:
the seven fields that would move are already identified and separated from the
fourteen that would not. What is missing is measured female values for them. The
evidence file carries sex-stratified cohort data where the sources have it — and
it also records that the beard density in this male body is **transferred from
female cheek vellus hair**, because that is what the retrieved source measured.

## 4. What was delivered, and what a female body would take

**Delivered.** One declared schema with provenance and consumers for seven
parameters. Stature as a working, gated, geometrically consistent knob, together
with the mechanism that makes a mis-scaled body a hard failure instead of a
silent one. A measured account of what sex anatomy this body has. Measured
sex-stratified proportions from data that was already on disk. A measured reason
why those proportions cannot be applied to the current geometry.

**Not delivered, and not attempted.** A female body. Nothing here is labelled
female, because nothing here would be entitled to the label.

What it would take, in order of what blocks what:

1. **A female whole-body mesh source.** Nothing catalogued has one. This gates
   internal genitalia, external genitalia and breast, all of which are *additional
   entities* and none of which can be produced by transforming male ones. The
   Visible Human Project female dataset is the obvious candidate and is not
   catalogued; whether its licence permits redistribution of derived meshes has
   not been checked here.
2. **A path refitting tool.** Anisotropic segment scaling — the pelvis, the
   shoulder:hip ratio, Q-angle — is measurable but unusable until the fitted path
   polynomials can be regenerated. This needs OpenSim in the environment. It is
   the smallest piece of work on this list and it unblocks every proportional
   parameter, for both sexes, including the ones that have nothing to do with sex.
3. **An anthropometric correspondence.** NHANES `BMXLEG` is inguinal crease to
   proximal tibia and `BMXARML` is acromion to olecranon. Neither is an OpenSim
   segment length, and substituting one for the other without a stated
   correspondence is exactly the class of error `docs/LOG.md` catalogues.
   Bi-iliac breadth, pelvic inlet shape and Q-angle are measured by **no**
   catalogued source, in either sex.
4. **A gonadal axis attached to a gonad.** The implemented hormone model is a
   female cycle model on a male body, coupled to nothing.
5. **Measured female values for the seven androgen-dependent hair fields.**

An honest summary of the current state, for anything that quotes it: *male
anatomy including a complete genital tract, no breast of either sex, one
implemented isotropic size parameter with five gates, seven measured proportional
sex differences that cannot yet be applied, and no female-specific geometry
anywhere.*

---

# 5. Where the parametrization reaches now: the whole body

Sections 1–4 describe a knob that scaled **the scaffold**. Ask for a tall person
and you got a tall 22-segment skeleton posing a fixed-size anatomy: 4,000
entities, 1.7805 m² of skin, 146 nerve routes and every organ, vessel and
ligament at the source subject's size. This section is that closed, and what
closing it measured.

    .venv/bin/python scripts/materialize_body_variant.py --stature-m 2.03

builds the scaffold and the body together and runs **48 gates** across five
stages, refusing to leave a variant on disk if one fails. At the identity
`--scale 1.0` it runs 52 and every quantity comes back **exactly**, not to a
tolerance.

## 5.1 No exponent is typed

`ihm/body_scaling.py`. The failure mode this file exists to prevent is an
exponent chosen because it looks about right. Areas go as `s²` and everyone
knows it; what gets got wrong are the mixed quantities. Is a ligament stiffness
`s`? `s²`? Invariant? A conduction delay? A receptor density? Typing any of them
is guessing, and a guess that is right cannot be told from a guess that is wrong.

So every exponent is **computed from a dimensional formula**, and the only
hand-written numbers are the primitives — each either 1 (a length) or 0 (a
material property). `exponent('ligament_stiffness')` returns 1 because `E·A/L`
expands to `0 + 2 − 1`.

**The primitives, which are the whole argument:**

| primitive | exp | why |
|---|---:|---|
| `length` | 1 | the isotropic factor itself — a modelling assumption, see 5.6 |
| `angle`, `fraction` | 0 | a uniform scale is conformal; a ratio is a ratio |
| `count` | 0 | a taller person has the same 31 spinal roots and the same 206 bones. **Every density derived from a count therefore carries a negative exponent** |
| `density` | 0 | composition. The anatomical 70.7713 kg is composed at per-constituent densities, which are properties of fat and bone, not of the specimen |
| `elastic_modulus` | 0 | 332.2 MPa (Quapp & Weiss 1998) is a property of collagen |
| `specific_tension` | 0 | cross-bridge density per unit myofilament area — molecular |
| `conduction_velocity` | 0 | **the load-bearing one.** Axon diameter and internodal myelin length are cellular dimensions |
| `viscosity`, `regulated_pressure`, `gravity` | 0 | blood composition; a baroreflex set-point; the world |
| `authored_damping_time` | 0 | an authored ratio in seconds, not a measured property. Held, and flagged — see 5.6 |

**Every exponent used anywhere, and its formula:**

| quantity | formula | exp |
|---|---|---:|
| length, centroid, moment arm, fibre length, tendon slack, route length | `length` | **1** |
| area, PCSA, patch area | `length²` | **2** |
| volume | `length³` | **3** |
| mass, weight, muscle mass | `density · volume` | **3** |
| inertia | `mass · length²` | **5** |
| max isometric force | `specific_tension · PCSA` | **2** |
| muscle volume | `PCSA · fibre_length` | **3** |
| joint moment | `force · moment_arm` | **3** |
| **conduction delay** | `route_length / conduction_velocity` | **1** |
| ligament stiffness (N/m) | `E · A / L` | **1** |
| ligament stiffness per strain (N) | `E · A` | **2** |
| ligament damping | `stiffness_per_strain · damping_time` | **2** |
| strain, stress, normalised fibre length, pennation, joint range | dimensionless | **0** |
| receptor density | `count / area` | **−2** |
| self-weight stress | `weight / area` | **+1** |
| strength-to-weight | `force / weight` | **−1** |
| vascular resistance | `viscosity · L / r⁴` | **−3** |
| volumetric flow | `pressure / resistance` | **3** |
| BMI | `mass / length²` | **1** |
| surface-to-volume | `area / volume` | **−1** |

The table checks itself three ways and will not load otherwise: against answers
fixed **outside** it (`mass = 3` and `inertia = 5` are `model_scaling.py`'s own
literals; `max_isometric_force = 2` is `body_parameters.py`'s default), against
**itself** where two formulas describe one quantity (`PCSA · fibre_length` must
come out at the same exponent as `length³`, and it only does if PCSA is really
an area), and structurally (every primitive used, every primitive 0 or 1).

## 5.2 The gate that matters: a taller body has a slower periphery

Conduction velocity is invariant. Route length is anatomy. So the delay is a
length divided by a constant, and **every conduction delay in the body lengthens
in exact proportion to stature**. At 2.03 m against 1.7973 m that is +12.95% on
all 1,743 recorded delays:

| route | length | v | delay | at 2.03 m |
|---|---:|---:|---:|---:|
| vagal C fibre | 507.8 mm | 1 m/s | 507.84 ms | **573.60 ms** |
| optic, magnocellular | 65.7 mm | 20 m/s | 3.28 ms | 3.71 ms |
| optic, parvocellular | 65.7 mm | 12 m/s | 5.47 ms | 6.18 ms |
| optic, koniocellular | 65.7 mm | 6 m/s | 10.95 ms | 12.36 ms |
| greater splanchnic C | 168.4 mm | 1 m/s | 168.36 ms | 190.16 ms |
| median A-beta | 610.0 mm | 55 m/s | 11.09 ms | 12.53 ms |
| sciatic Ia | 200.8 mm | 100 m/s | 2.01 ms | 2.27 ms |

The three retinal populations ride **one** optic nerve, so their spread is a
delay too and separates in proportion: 7.66 ms → 8.65 ms. All seven are named in
advance in `HEADLINE` so the report cannot be read as whatever came out.

This is the consequence that says the parametrization reached the body.
`docs/MILESTONES.md` records that lumping peripheral delays costs as much as
deleting an entire fibre group, so delays are load-bearing for the brain and a
body of a different size genuinely has a different periphery.

**Central and synaptic delays deliberately do not scale.** They are peripheral
route lengths in neither name nor fact, and nothing here fixes their scaling.
They are in the table as explicit nulls, so that forgetting them and excluding
them do not look alike. A whole reflex latency is therefore *not* 12.95% longer;
only its peripheral part is.

## 5.3 What scales, and what does not

**Scales, all of it gated:**

- **4,000 anatomical entities** — centroid, bounds, surface area (`s²`), volume
  (`s³`), the `connections[].distance_m` on 3,996 of them, the three skin
  layers' `shell.thickness_m` and depth intervals, and the segment binding's
  centroids and registration translation. Meshes are **not** rewritten; the
  variant carries a uniform-scale transform, exactly as the OpenSim model
  carries `scale_factors` rather than resized geometry.
- **146 nerve routes, 249 muscle bindings, 16 receptor patches, 20 relays** —
  positions, route lengths, anchors, and 1,743 delays.
- **1,326 skin patches** — position, area (`s²`), afferent route length, and the
  area budget the bisection rule uses.
- **98 actuators** — optimal fibre length, tendon slack, max isometric force, and
  the derived PCSA, muscle volume and muscle mass.
- **117 ligament force elements** — attachment points, slack length,
  cross-section, `E·A` stiffness, damping, volumes and areas.

**Deliberately does not scale, and why:**

| held | reason |
|---|---|
| all 17 fibre-class conduction velocities | cellular dimensions. This is the invariant the headline result *is* |
| specific tension, tissue density, elastic moduli | material and molecular properties |
| every strain, angle, fraction and normalised length | dimensionless |
| ligament `peak_strain_over_declared_range` | a strain — so a taller body strains its ligaments exactly as much over the same joint range, and the elements flagged `kinematically_admissible: false` stay flagged at every stature. Scaling neither rescues one nor breaks one |
| **patch count**, held at 1,326 | a decision, defended in 5.4 |
| central and synaptic delays, membrane time constants | not peripheral route lengths |
| receptor gains (Hz/Pa, Hz/°C) | transduction properties of the ending |
| `provenance.source_to_canonical` | provenance of a *file*, not a body quantity. Folding the stature scale into it would make the record of where the bytes came from untrue |
| the passive joint stops, Hunt–Crossley contact parameters | the same gaps `model_scaling.py` already records; inherited unchanged |

**Is anything's absolute size fixed by chemistry?** Measured rather than assumed.
The smallest extent in the atlas is **0.685 mm** (nuchal ligament), **34×** the
largest size-invariant cellular dimension (myelinated axon 1–20 µm, erythrocyte
~8 µm, capillary lumen ~4 µm, sarcomere ~2.7 µm). **Zero of 4,000** entities are
within an order of magnitude of it. So the chemistry-fixed constraints act
*below* this atlas' resolution and appear in this body as invariant **properties**
— conduction velocity, specific tension, density, modulus — which is exactly
where the exceptions in 5.1 live. An entity that did need the exception would be
expressible: placement and size exponents are carried separately per entity for
that reason. None needs it.

## 5.4 The patch count is a choice, and it is made explicitly

`build_dermatome_patches.py` bisects until each patch is under an area budget.
Re-run that on a bigger body with the budget unchanged and you get about `s²`
as many patches — **1,692** at 2.03 m. Scale the budget as `s²` and you get the
same **1,326**, larger. Both are readings of the same script and they are
different bodies.

**Chosen: count fixed, area `s²`, receptor density `s⁻²` — 744.8 → 583.8
patches/m², a 21.6% fall.**

1. A patch is an **afferent channel**, not a square centimetre. Each carries a
   nerve, a relay and a cortical target, and IBM-1 consumes them as a fixed input
   set. A count that followed body size would change the dimension of the brain's
   sensory input with stature — a 1.4 m and a 2.0 m body could not run the same
   trained model — and the parameter would stop being a property of the body and
   become one of the interface.
2. The anatomy agrees: a taller person has the same 31 spinal nerve pairs and the
   same 32 dermatomes.
3. Real receptor counts do not follow body size — tactile acuity is better on
   smaller fingers because Merkel-cell density is higher on them. Innervation is
   laid down as a roughly fixed count and spread over whatever surface grows.

`--count-follows-area` reports the other answer in the same run, so this reads as
a decision rather than as the only thing the code could do.

**What it is not a claim about.** `dermatomes.json` already records that this
density is two-tier by *area*, not by measured receptor density, and that real
fingertip innervation exceeds trunk innervation by an order of magnitude.
Holding the count preserves that limitation exactly. **The density that falls as
`s⁻²` was never a measured density.**

## 5.5 The finding: the engine divides the scaling straight back out

`scripts/native_mechanical_stream.cpp` line 69:

    double mass_scale = target_mass / original_mass;

and every body mass is multiplied by it. So the plant weighs **whatever
`target_mass_kg` it was handed**, and the `s³` that
`materialize_stature_variant.py` applies to the 22 body masses is annihilated.

Ask for 2.03 m and pass the standing `target_mass_kg = 77.6122029` — the literal
in 58 files — and you get a 2.03 m skeleton weighing 77.61 kg. **A person 13%
taller and not one gram heavier.** Every gate on the mechanical side passes while
it happens, because the mass scaling was applied correctly and the annihilation
is downstream in the engine.

The scaled model carries 122.869 kg before the engine touches it. A caller who
wants the geometric-similarity mass must pass **111.835 kg**; the NHANES
population says **99.424 kg** (5.6). Neither is the default.

Measured statically, from the scaled `.osim` and the engine's own formula.
**The standing-weight-equals-m·g and 1e-14 momentum-residual gates need a running
plant, and no native run has ever been accepted on a scaled body on either side**
— `native_acceptance_complete` is `false` in every variant this repository has
written — so they are named open here rather than reported as passed.

## 5.6 Geometric similarity is false, and by how much

`scripts/measure_stature_allometry.py`. The isotropic assumption makes a testable
prediction about a population, and `BMX_J.parquet` — already on disk, already
read by `index_anthropometry.py` for percentiles — can answer it.

Survey-weighted over 4,822 adults aged 20–79:

| | measured exponent | isometry says | z |
|---|---:|---:|---:|
| mass ~ stature^b, men | 2.372 | 3.00 | −4.4 |
| mass ~ stature^b, women | 1.721 | 3.00 | −7.5 |
| **mass ~ stature^b, pooled** | **2.034 ± 0.083** | **3.00** | **−11.6** |
| BMI ~ stature^b, pooled | 0.034 | 1.00 | −11.7 |

The pooled slope lies *between* the two within-sex slopes, so it is not the
artefact that pooling two groups with different means can manufacture; and the
weighted mean heights land on the published CDC values to 0.1 cm, so it is not a
weighting error. BMI being flat in stature agrees with this repository's own
already-measured male/female BMI difference of d = −0.03.

**Reported, not applied.** `mass_kg` is an independent knob and this scaler does
not move it; the number exists so a caller who takes the `s³` default takes it
knowingly. At 1.40 m the default is 21% too light; at 2.05 m, 13.6% too heavy.

`ihm/body_scaling.py ALLOMETRY` carries six such entries, each with the size of
the error and a disposition:

| where isotropy is wrong | isotropic | truth | applied? |
|---|---:|---|---|
| mass vs stature | `s³` | `s^2.034`, measured here | no — reported |
| BMI vs stature | `s¹` | `s^0.034`, measured here | no — the falsification |
| vascular flow | `s³` | Kleiber `s^2.25`; **9.6% over-perfused at 2.03 m**, because Poiseuille's `r⁴` beats the length and resistance falls as `s⁻³` | no — no vessel radii are rewritten and no flow model reparameterised |
| brain volume | `s³` | far shallower; **no catalogued source here measures it**, so no exponent is asserted. At 2.03 m the brain is 1.44×, which is very probably too big | no — flagged as known-wrong, not unexamined |
| characteristic time | `s⁰` | `s^0.5` under equal Froude number — why a taller person walks at a lower cadence. 6.3% at 2.03 m | no — rescaling time would change the meaning of every time constant in the body at once |
| passive joint stops | `s³` | unknown; authored range-of-motion limits, not measured tissue | no — the same gap `model_scaling.py` already records |

Two consequences that are *not* errors but fall out of the algebra and are worth
saying out loud: **self-weight stress rises as `s`** (+12.9% at 2.03 m) and
**strength-to-weight falls as `s⁻¹`** (−11.5%). That is what geometric similarity
means, it is the classical reason large animals are not scaled-up small ones, and
it bears directly on whether a tall body can pick itself up.

## 5.7 The gates, and the ones that would not have caught it

Every stage carries an **independent arm** — a check that reads a different
artifact through a different code path, because a gate that only checks your own
writes cannot see a consistent mistake.

| stage | the independent arm |
|---|---|
| anatomy | scale the **vertices** of the gzipped triangle meshes `anatomy.json` only points at, and re-integrate. A sum of cross products and a divergence-theorem volume over scaled coordinates are not "the record times a power of `s`" |
| nerves | rebuild all 411 route lengths from geometry by the method each route declares — and for the 249 muscle bindings the endpoint is an **entity centroid out of `anatomy.json`**, a different artifact from a different script. All 411 reproduce their declared length exactly at the base, which is what entitles the recomputation to be believed at the scaled one |
| skin | the 1,326 patch areas must sum to the exterior area integrated over 109,183 **scaled** triangles selected by the recorded exterior id list — 109,183 triangles that have never heard of a patch |
| muscle | mechanical muscle volume (PCSA × fibre length, mechanics catalog) over anatomical muscle volume (mesh-integrated, `anatomy.json`). They disagree 5.56× at the base — 98 actuators are not 568 muscles — and the disagreement must be invariant |
| ligaments | stiffness ÷ cross-section must recover **332.2 MPa** from a *third* file. A modulus is a material property and must come back at the same number on a body of any size |

**Seven sabotage arms, all caught, and what did not catch them is the point.**

| sabotage | caught by | passed anyway |
|---|---|---|
| area given the length exponent | mesh integration, `V = A·t` | centroids, distances, the two-body stature agreement, composed mass |
| volume given the area exponent | mesh integration, `V = A·t` | the same four |
| skin thickness left unscaled | `V = A·t`, leaf completeness | everything else |
| patch area given the length exponent | patch-sum-vs-mesh, coverage, budget fraction | the count gate and the route gate |
| conduction velocities scaled too | delay arithmetic, the table | **every geometry arm** |
| routes scaled, body left alone | **only** the geometry arm, by 55.8% | **every delay gate, perfectly** |
| max isometric force given the volume exponent | **only** the check against `KNOWN_ANSWERS` | PCSA, muscle volume and specific tension all still agree with each other — they were derived from the one mistaken formula |
| Blankevoort stiffness read as N/m | **only** the 332.2 MPa known answer | every other ligament gate |

The last three are the lesson. A scaling table can be **self-consistently
wrong**, and every check derived from it will agree with every other one. Only a
comparison against something fixed elsewhere — another module's literal, another
file's published constant, another artifact's geometry — can see it.

And the `max(0.03, …)` floor in the route method is not homogeneous, so whether
it ever bites is measured rather than hoped: the shortest route is 58.1 mm, the
floor activates below `s = 0.517`, and the bottom of the declared stature range
is `s = 0.779`.

## 5.8 What this still is not

Isotropic. One factor everywhere. It changes how **big** the subject is and not
what **shape** it is: every segment length ratio, mass fraction and radius of
gyration is exactly the source subject's at every stature, and section 3's
account of why anisotropic proportions need the fitted paths *refitted* rather
than rescaled is unchanged. The meshes are carried by a transform rather than
rewritten. And no native run has been accepted on a scaled body.

What has changed is that the parametrization no longer stops at the scaffold.
