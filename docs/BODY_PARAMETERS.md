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
    scripts/measure_anisotropic_pelvis_error.py   why a coefficient scale cannot carry them
    scripts/native_polynomial_path_fit.cpp        drive OpenSim's PolynomialPathFitter
    scripts/build_native_path_fitter.py           build it against the local install
    ihm/native/path_refitting.py                  run it, read it back, compare
    ihm/native/anisotropic_scaling.py             per-body, per-axis geometry
    ihm/native/body_measurements.py               measure shape, not only size
    scripts/verify_path_refitting.py              five gates on the refitter
    scripts/materialize_proportional_variant.py   build a re-proportioned body
    scripts/audit_hormone_axis_anatomy.py         which hormone fields have an organ
    scripts/audit_hair_sex_dependence.py          what a female hair field would need
    scripts/verify_proportional_variant_native.py  does the ENGINE accept it

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
| `pelvis_breadth_ratio` | – | 0.85 – 1.20 | **implemented**; hip separation over stature |
| `shoulder_breadth_ratio` | – | 0.85 – 1.20 | **implemented**; acromial separation over stature. **No measured sex value exists** |
| `leg_length_ratio` | – | 0.85 – 1.20 | **implemented**; (femur+tibia) over stature |
| `sex` | – | **{male, female}** | **implemented — proportions only, anatomy unchanged** |
| `age_years` | y | 18 – 90 | declared; reaches nothing |
| `body_fat_fraction` | – | 0.05 – 0.50 | declared; reaches nothing |

*surfaced* means the knob existed and the schema only names, bounds and
attributes it. *implemented* means this work made it one. *declared* means it is
written down and reaches the running body through nothing at all, and the
schema says which of the three every parameter is rather than presenting a flat
list in which they look alike.

`resolve({'sex': 'female'})` used to raise, and the raise was right: accepting
the value would have changed a string in a JSON file and nothing else. It stopped
raising when it stopped meaning nothing. It now sets four measured proportions --
stature, mass, pelvic breadth over stature, leg length over stature -- and each
one changes the running mechanical body by an amount that is gated. Its status
string carries the parenthesis in full: **implemented (proportions only; anatomy
unchanged)**, and every resolved record carries a `sex_realisation` block naming
what moved and what did not, so the label cannot be quoted without the caveat.

What did not move: the anatomical entity set. Still 0 female-specific entities,
still 0 mammary glands in either sex, and the male genital tract still bound to
the pelvis.

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

**Amended 2026-09-10: one female source does exist, for the lower body.** The survey above stands
for whole bodies and is wrong for the pelvis-and-below. `data/sources/visible-human-lower-extremity.json`:
the Visible Human **Female** and **Male** cryosections segmented into 260 geometries each -- 76
muscles, 28 bones, **16 cartilages, 8 ligaments, and subcutaneous and intramuscular FAT** --
CC BY 4.0, public download, no agreement (Andreassen et al., Scientific Data 2023;10:34). The
female cryosections are 0.33 mm apart against the male's 1.0 mm, and NLM has required no licence
for the underlying images since July 2019.

It is a female LOWER body, not a female body: no thorax, no viscera, no head. But it answers three
things this repository has open, and answers them from ONE specimen each:

* **articular cartilage**, which this body has none of and whose absence made the knee pilot fail
  0 of 6 (`docs/TISSUE_MECHANICS.md`) -- here it comes with the femur and tibia it belongs to;
* **fat as geometry**, subcutaneous and intramuscular, against the 4.887 kg of declared fat this
  body carries as denser bone and muscle instead;
* **the foot**, where the skin-on-bone registration fails the stance gate -- a foot whose bones,
  muscles and fat are one person's needs no cross-specimen warp.

What it is not: MakeHuman and the parametric surface models (SMPL-X and its skeleton-inferring
descendants) are skins. MakeHuman's assets are CC0 and its surface is better than anything here,
but the mesh is an envelope with helper cubes for animation bones -- no femur, no pectoralis, no
uterus. They are legitimate as a shape TARGET to warp measured anatomy toward, or for display and
clothing, and they are not anatomy. The whole-body female options that are anatomy remain: the
ICRP-110/145 reference female (a real 43-year-old's CT, 141 tissues, but the phantom data ships
with the publication and its terms are unread here), and segmenting the Visible Human Female
above the pelvis ourselves.

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

### The pelvic term is a knob now, and here is what it cost

`measure_anisotropic_pelvis_error.py` measured the reason it was not one. Widen
the pelvis medio-laterally by the measured female/male hip-over-stature ratio
**1.1265**, hold everything else:

- **56 of 98** paths change; **42 do not** (50 of 80 on the example model);
- among those that change the span is **0.01% to 6.23%** — `addmagProx` 6.2%,
  piriformis 4.8%, `addbrev` 4.0%, `glmax3` 2.0%, and nothing below the knee.

A polynomial coefficient scale admits **one** factor, so any single factor is
wrong for every muscle but at most one. That conclusion stands. What was wrong
was the next sentence: *"refitting means OpenSim's `PolynomialPathFitter` … the
OpenSim Python bindings are not installed here."*

**They are not, and they were never the only route.**
`data/runtime/opensim/install/opensim/lib/libosimActuators.so` — built by
`scripts/build_native_opensim.py`, present since the engine was built — exports
**97 `PolynomialPathFitter` symbols**. The fitter was not missing. It was
undriven. `scripts/native_polynomial_path_fit.cpp` is 300 lines and links against
the install that was already there.

#### What the refit is worth, against the file the engine actually runs

`verify_path_refitting.py`, five arms, all passing. Measured on **27 held-out
frames** of the reference trajectory against the model's own `GeometryPath`s
with their wrap objects, never against the fitter's own printed RMS — that is a
training error, measured on the samples it fitted.

| | path length RMS | moment arm RMS |
|---|---|---|
| refit of the unchanged model | **3.01e-4 m** | 3.55e-3 m |
| the **shipped** `FunctionBasedPathSet` | **5.02e-3 m** | 3.64e-3 m |

The engine runs the shipped set. A refit of the same model is **16.7× closer** to
the muscle paths that model declares.

#### Two things the gate design had to stop assuming

**The fitter is not deterministic.**
`LatinHypercubeDesign::computeRandomHypercube` seeds `std::mt19937` from
`std::random_device` on every call. Two refits of the same model differ by
**3.12e-5 m** RMS. So the shipped coefficients cannot be reproduced bit for bit
by anyone, including whoever produced them, and a gate demanding it would be
testing the wrong property. That run-to-run spread is the **noise floor**, and
nothing measured on a refit means anything unless it clears it.

**The shipped set is not the reference.** The polyline arm — fitted polynomial
length against a straight-line walk of the model's own `PathPoint` locations,
computed in Python by forward kinematics over the `.osim` — was first written to
demand that a refit not move that ratio away from the *shipped* set's. That
quietly made the shipped set correct by definition. Scored against the
`GeometryPath` truth instead: the refit's worst drift is **0.81%** (`gaslat_l`),
the shipped set's is **5.07%** (`ehl_r`).

#### The arm whose answer was fixed in advance

Scale the model **1.10× isotropically**. Path length is homogeneous of degree one
in the geometry at fixed pose, so the true lengths must be exactly 1.10×. They
are, to **1.4e-7 m** — and that residual is the constraint assembler, not the
scaling: this model carries four `CoordinateCouplerConstraint`s and the worst
muscle is `gaslat_l`, which crosses the coupled walker knee. (That `scale_model`
is exact is established without a solver by `verify_stature_scaling.py`, whose
identity arm reproduces every path length to relative error 0.0 from the XML.)
Then refit the scaled body from scratch: it lands on 1.10× the base refit at
**1.74× the noise floor**. Analytic route and numerical route share no code.

#### The arm that had to fail

A refit of the **wide-pelvis** body, scored against the **original** body's
truth: **72× the noise floor**. The same refit against its own body: 3.05e-4 m.
A gate that has never been seen to fail is not a gate.

### Proportions, and what they measurably do

`materialize_proportional_variant.py --sex female`. Nine gates, all passing,
measured on the artifact that was written:

| quantity | change | target |
|---|---|---|
| hip separation / stature | **+12.6461%** | +12.6461% (error 1e-13) |
| leg length / stature | **−2.5526%** | −2.5526% (error 1e-13) |
| shoulder breadth / stature | 0.0000% | 1.0 — **no measured value exists** |
| stature | −8.0058% | exact by construction |
| total mass | −19.97% | — |
| pelvis segment mass fraction | +10.5% | — |
| hip adduction to keep the stance | **+0.2253°** | — |

**Order of operations, and it is not cosmetic.** Every proportional parameter is
defined *relative to stature*, because that is how the anthropometry measures it.
So the shape factors go on first, the shaped body is **measured**, the fixed
point is solved on that measurement, and one isotropic factor then lands stature
exactly. The legs are 48% of stature, so shortening them shortens the very
denominator the request was written against: assuming the raw factor equals the
ratio is wrong by **1.2%**, which is half the size of the effect being asked for.
The solver converges to 1e-13 in 38 iterations and 3.4 s.

**The Q-angle claim that was not made.** A wider pelvis does **not** tilt the
femur at the default pose — it translates the whole leg laterally, and femoral
obliquity is unchanged at 0.264°. What moves is where the foot lands. The honest
consequence is the hip adduction that would put it back: **+0.225°**. The model
carries no ASIS marker and no tibial tuberosity, so it cannot measure a clinical
Q-angle, and `body_measurements.py` names the quantity `femoral_obliquity_deg`
for that reason.

#### Two measurements that say the refit was not optional

The variant's true muscle paths sit **7.78e-3 m RMS** from the isotropically
scaled source subject's — **250× the noise floor** — and that scaled subject is
exactly what a coefficient scale of the shipped polynomials would have produced.

And the per-muscle fibre-length factor, taken by OpenSim's own
`Muscle::extendPostScale` rule from each muscle's own `GeometryPath`, spans
**0.8888 to 0.9737** across the 98 muscles. One number is wrong for 97 of them.

#### A side effect worth stating rather than discovering

The refit covers **all 98** `PathActuator`s; the shipped set covers 80. The 18 it
adds are the `arm26` and `gait2392` trunk muscles, which have never had fitted
moment arms — so `JointPosturalController` can now actuate them. It also removes
their explicit `GeometryPath` in the engine, which is the geometry a retinaculum
or a fascia would need to constrain. Both are in the registration.

#### What is still an assumption

- **Wrap objects.** Every wrap in this model is a rotated `WrapCylinder`, and a
  rotated cylinder under an anisotropic scale is an elliptic cylinder, which
  OpenSim cannot express. The default policy scales length by the stretch along
  the cylinder's own axis and radius by the geometric mean of the two principal
  stretches perpendicular to it — cross-sectional area preserved, exact under
  isotropy. `translate_only` is the alternative and the per-wrap factors are
  written into every variant.
- **`leg_length_ratio` is isotropic per segment**, so it slims the leg as it
  shortens it. A long-axis-only femur scale is *refused*: `walker_knee_r` takes
  its roll-glide translations in a frame rotated by (−1.64, 1.45, 1.57) rad
  inside `femur_r`, where componentwise scaling has no meaning.
- **The correspondence problem is not solved, it is declared.** NHANES measures
  hip *circumference* at the buttocks; `pelvis_breadth_ratio` is a *skeletal*
  medio-lateral factor. Using one as the other assumes soft tissue and bone scale
  together in that direction. Bi-iliac breadth is measured by no catalogued
  source, in either sex. And `BMXLEG` is inguinal crease to proximal tibia, not
  the OpenSim hip-centre-to-knee-centre length.
- **The correspondence problem is not solved, it is declared** (above).
- **Nothing here walks, stands, or is driven by anything.** The native check
  below is an *acceptance*, not behaviour.

#### The engine accepts it, and finding that out found a defect

`verify_proportional_variant_native.py` loads the variant and the base body
through `NativeMechanicalStream` under the same requested mass and prints both,
because one body's numbers alone say nothing. Five 2 ms steps, then close.

| | variant | base |
|---|---|---|
| muscles / coordinates / bodies | 98 / 33 / 22 | 98 / 33 / 22 |
| plant mass | 66.0928 kg | 66.0928 kg |
| normalised fibre length | **0.569 – 1.278** | 0.573 – 1.281 |
| integrated 5 × 2 ms | yes | yes |
| worst step | 0.228 s | 0.294 s |

`native_acceptance_complete` is `true` on this variant, and it means exactly what
the registration says it means: *the plant instantiated this geometry with these
refitted muscle paths and integrated*.

**The first run of this check failed, and the failure was real.** The variant's
normalised fibre length came out at **0.381** against the base's 0.573, on
`fhl_l` and `fhl_r` — **below the Millard active force-length floor of about
0.47**, so flexor hallucis longus was producing no active force at all.

The cause: `scale_model` puts the *global* isotropic factor on every
`optimal_fiber_length` and `tendon_slack_length`, and the per-muscle correction —
0.8965 for the leg against the global 0.9447 — went only into `catalog.json`.
**The engine reads fibre and tendon lengths from the `.osim` and never from the
catalog.** A 5.4% error in a musculotendon that is 90% tendon lands almost
entirely on the fibre: 0.38 m × 0.054 = 20.5 mm of extra tendon against a 43 mm
optimal fibre, which is the 0.594 that was observed.

That is this repository's recurring shape — *a quantity computed correctly, then
applied to the wrong object* — and no XML gate had caught it. Two gates now do,
and both would have caught it without launching the engine:

- **muscle operating point does not move**: `(path length at the default pose −
  tendon slack) / optimal fibre length`, dimensionless, which is what the
  equilibrium solve lands on. It moved by 0.594 under the defect; it now moves by
  **1.1e-15**.
- **per-muscle fibre factor reached the model, not only the catalog**: reads the
  `.osim` back and requires each muscle to carry its own measured path ratio.
  **2.2e-16**.

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

#### The absence is measured now, not described

`audit_hormone_axis_anatomy.py` counts the entities that would carry each
declared field, over the 4,000 segment-bound set:

| field group | supporting entities | with a computed volume |
|---|---|---|
| `reproductive.*` | **11**, all bound to `pelvis`, all male | 3 |
| `uterine.*` | **0** | — |
| `placental.*` | **0** | — |

Two of those rows belong on `docs/DISCONNECTS.md`: they are declared models with
no running object behind them, and until this audit they looked in the inventory
exactly like the rows that do have one.

**The audit's known-answer check failed, and the failure is the useful part.**
The gonad meshes measure **6.7 and 6.8 mL** each against a 12–30 mL adult
clinical reference range, and the prostate mesh is **11.6 mL** against a typical
adult 15–25 mL. Both measurable reproductive organs in this atlas are small by
about the same factor — a property of a single-donor atlas whose skin extent is
1.7195 m and whose organ volumes have never been calibrated against a population.
Any hormone axis scaled by gonadal volume would inherit it, which is why the
number is reported *before* anything is scaled by it.

It also means the one coupling this body could support — a testosterone axis
attached to a real testis — would be attached to a testis half the expected size.
That is a data problem, not a modelling one, and it is now on the record instead
of waiting to be discovered downstream.

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
fourteen that would not.

**But the missing piece is not the androgen variable.** `audit_hair_sex_dependence.py`
reads the built evidence and finds the parameter that sexual dimorphism in body
hair actually runs through already present on every field:
`shaft_bearing_fraction`, the fraction of counted follicles that produce a
visible shaft. The evidence file's own note says why — follicle *number* is
roughly sex-invariant (seago1985 measured no sex difference on thigh or upper
arm); what differs is how many of those follicles make a terminal hair.

Of the seven androgen-dependent fields:

- **six** (axillary, chest, abdomen, back, pubic, perineal) set that fraction to
  **1.0** at tier `assumed`;
- **one** (beard) declares it **unmeasured** and leaves the mass absent rather
  than inventing it;
- **none** have a measured value.

A fraction of 1.0 — every counted follicle bears a shaft — is specifically an
**adult male** assumption, and it is not labelled as one anywhere in the field
record. The beard's refusal is not an inconsistency but the correct call: its
density is transferred from *female cheek vellus* hair, where assuming a terminal
shaft per follicle would be badly wrong.

So wiring an androgen signal into these fields today would drive a quantity whose
male value is itself an assumption. What would close it is **seven numbers**, and
the standard female body-hair instrument, Ferriman–Gallwey scoring, is *ordinal*
and does not supply a fraction. Until they exist, `sex` reaching hair would move
a rendered surface and no measured quantity — which is the state `sex` was
raising over in the first place.

## 4. What was delivered, and what a female body would take

**Delivered.** One declared schema with provenance and consumers for ten
parameters. Stature as a working, gated, geometrically consistent knob, with the
mechanism that makes a mis-scaled body a hard failure instead of a silent one. A
measured account of what sex anatomy this body has. Measured sex-stratified
proportions from data that was already on disk.

And, since: **a working path refitter with five gates**, which was the item this
document called the smallest on the blocker list and mislocated — it needed a
driver, not an installation. Three anisotropic proportional parameters that
change the body by amounts measured against a refit noise floor. `sex` as a
parameter that stopped raising because it stopped meaning nothing.

**Still not delivered.** A female body. Nothing here is labelled female, because
nothing here would be entitled to the label. `--sex female` produces **male
anatomy at female-typical proportions** and the registration says so in those
words.

What it would take, in order of what blocks what. Item 2 is struck because it is
done; the others have moved.

1. **A female whole-body mesh source.** Still the gate on internal genitalia,
   external genitalia and breast, all of which are *additional entities* and none
   of which can be produced by transforming male ones. **Still not delivered, and the
   word that matters is WHOLE-body**: the Visible Human Female is now catalogued
   (`data/sources/visible-human-lower-extremity.json`, CC BY 4.0, 260 segmented
   geometries per sex including cartilage and fat) and it stops at the pelvis, so it
   reaches none of the three entities named in this item. See the amendment under
   "Does any catalogued source ship a female variant". **Moved, as of 2026-09-10.**
   zenodo.org answers from this machine (the HTTP 504 was the old sandbox's
   egress proxy), so both catalogued cards were checked against the archive and
   both were wrong in places:
   * **TotalSegmentator** (record 10047292, v2.0.1, **CC-BY-4.0** verified, one
     23.58 GB zip of 147,361 members). Its own `meta.csv`, read out of the zip by
     HTTP Range without downloading it: 1,228 subjects, **510 female** (not the
     relayed 503), 716 male, 2 blank. The breast is **not in the dataset**: the
     archive holds the 117 total-task labels, and a breast has to be produced by
     running the Apache-2.0 `breasts` subtask model on a female CT. The same CT's
     Apache-2.0 `body` subtask gives a female **skin** envelope, which matters
     because this body's skin is male too.
   * **UT-EndoMRI** (record 13749613). The relayed "CC-BY" is **wrong**: the record
     carries no licence field, and its User Agreement reads *"available for free
     use exclusively in non-commercial scientific research"*. Whether this
     programme qualifies is the **owner's** decision; until it is made nothing is
     downloaded and no derived uterus or ovary mesh goes into a public repository.
     It remains an endometriosis cohort, and the record says ovaries may be
     deformed or surgically absent.

   **Selection had to be by what the scan contains, not by the exam's name.** The
   archive's `study_type` names the clinical exam and the member is a CROP: s1218,
   "neck-thorax-abdomen-pelvis", is a 68 mm slab. `select_totalsegmentator_by_coverage.py`
   uses the zip's central directory -- an empty label mask compresses to ~250 B --
   to find subjects whose 24 ribs, sternum and clavicles are all in the scan (219
   of 510 female; 71 of those without pathology), gated on a subject whose true
   voxel counts are known. `screen_totalsegmentator_fov.py` then reads each CT's
   field of view from its NIfTI header, because three of the first four failed on
   it.

   **One female torso passes every gate: s0790**, age 40, an 807 mm scan, 39 of 39
   registration bones, breasts **610 / 551 mL**, chest whole, breast whole, its own
   bones inside its own body mask, the breast inside the body (to one voxel -- the
   gate was amended after s0790 failed the strict version at 0.9825, and every
   subject's excess is a one-voxel rim) and off the ribs, laterality decided from
   the data. The original gate set had PASSED a 68 mm slab with a 50 mL
   cross-section for a breast; the coverage gates were added because of it, and a
   negative control shows they refuse it. What this is and is not: one clinical
   subject's anatomy as a segmentation model drew it, and "breast" is a single
   soft-tissue label -- no gland, duct, nipple or areola.

   **Registered, and one gate still fails, by millimetres.**
   `scripts/register_female_torso.py` fits one similarity from s0790's thorax to this
   body on 39 bone correspondences (41 BodyParts3D entities; four bones exist twice,
   once from Z-Anatomy -- T5, T11, T12 and the manubrium -- and only one copy may
   enter), refined by ICP (scale 1.064, trimmed surface RMS 5.0 mm).
   * a: recovers a known similarity, 6.4e-16. **PASS**
   * b: leave-one-bone-out median 10.5 mm against a 29.7 mm neighbouring-rib null. **PASS**
   * c: laterality after mapping. **PASS**
   * d: this body's ribs are inside the mapped female trunk (0.997 / 0.999), but its
     sternum is only 0.767 inside, and 4-5% of its rib vertices sit inside the mapped
     breasts. **FAIL**

   The failure has been taken apart:
   * **Not the xiphoid.** TotalSegmentator's sternum label stops short of it (mapped,
     it spans -66.5..+80.4 mm along this body's sternal axis; the xiphoid lies at
     -91..-71), which biased the correspondence -- but refitting without it moves the
     sternum residual 23.8 -> 20.4 mm and gate d **not at all** (0.767 both ways)
     (`scripts/measure_sternum_correspondence.py`).
   * **A symmetric forward offset of this body's chest wall.**
     `scripts/measure_female_chest_wall_offset.py`, classifying every point first:
     the sternum points outside the female trunk are 95-97% ANTERIOR of it, by a
     median 2-4 mm; the rib points inside the breasts are ribs 2-7 on both sides,
     deepest at rib 5 (median 8.0 / 8.3 mm, max 13.9 / 13.6 mm), tapering to 2-3 mm at
     ribs 2 and 7. That is this body's chest wall bulging about a centimetre further
     forward across the breast base -- a chest-SHAPE difference at the scale of the
     fit's own residual, not a breast on the wrong ribs.
     **Caveat added 2026-09-10, and it applies wherever this figure is quoted:** these are
     nearest-point distances, the estimator later measured to carry a d^2/R bias and to
     select for already-aligned regions. An independent reader along the surface normal
     agrees within 5% on s0790 and sits 36-45% below on the other three. Quote it as
     **6-11 mm by a biased estimator, 5.5-6.6 mm by another**, and see the distance-field
     section below before relying on either.

   **Stray fragments, and a third subject that passes everything.** TotalSegmentator
   sometimes labels a speck far from the bone it names. s0970's "cut" ribs were specks
   of 1-25 voxels 138-372 mm from their bones, scattered through the lower scan. The
   extraction now drops a label's piece only if its closest approach to the label's
   largest piece exceeds 50 mm. Across all eight subjects every kept piece lies within
   11.2 mm and every dropped one at least 73.5 mm away; "keep the largest piece" would
   have been wrong, because s0790's and s1159's vertebrae carry real bone split off by
   a sub-voxel gap, 7.5-11.2 mm away. Regression: s0790 and s1159 drop nothing and
   their gates and volumes are identical. **3 of 8 subjects now pass every gate:
   s0790, s1159, s0970.**

   **Three subjects, one answer: it is systematic.** Two more women were registered the
   same way. s1067 (38) on the 34 labels wholly inside her scan
   (`--whole-labels-only`; her scan cuts ribs 10-11 and right 12, keeps ribs 2-7 and the
   breast whole), and s1159 (47), who passes every extraction gate, on all 39. Zenodo's
   file backend was down (HTTP 504 at 30.5 s, every request); both came through a
   HuggingFace copy verified byte-identical to Zenodo's central directory -- all
   147,361 members -- with every member still checked against the zip's own CRC32.

   | | s0790 | s1067 | s1159 | s0970 |
   |---|---:|---:|---:|---:|
   | gates a-c | pass | pass | pass | pass |
   | this body's sternum inside the mapped trunk | 0.767 | 0.635 | 0.315 | **0.090** |
   | ribs inside the mapped breasts | 5.1% / 3.9% | 4.2% / 3.0% | 3.4% / 4.4% | 5.2% / 5.0% |
   | out-of-trunk sternum points that are ANTERIOR | 95-97% | 95-98% | 92-99% | 99-100% |
   | ...by a median of | 2-4 mm | 3-5 mm | 5-6 mm | **11-19 mm** |
   | deepest rib inside the breast | rib 5 | rib 5 | rib 5 | rib 6 (5 close) |
   | its median / max depth, L; R (mm) | 8.0/13.9; 8.3/13.6 | 10.5/17.2; 7.3/14.9 | 6.7/12.6; 7.8/12.7 | **14.1/21.1; 13.4/20.9** |

   In all four the DIRECTION is the same -- this body's anterior chest wall in front of
   hers -- and so is the PROFILE: peaking mid-chest (ribs 5-6) and tapering above and
   below. The MAGNITUDE is not: from a few millimetres to about two centimetres, with
   s0970 roughly twice the others. (A first version of this paragraph, written on three
   subjects, said "about a centimetre at worst"; the fourth overturned it.) s0970 is also
   where the fit is least certain -- its sternum carries the largest residual of any
   subject (34-36 mm) and its fitted scale is the lowest (1.014 against 1.06-1.10) -- so
   one subject cannot say how much of her larger offset is anatomy and how much is
   registration. The direction is a chest-shape difference between this body and female
   chests; the size varies between women.

   **What would move it next.** One systematic correction: seat each breast on this
   body's own chest wall instead of refitting the thorax per subject. It is not built
   yet, and when it is, gate d2 cannot judge it -- a correction built to clear the ribs
   passes that gate by construction. Its judge has to be something it does not
   optimise: that the breast base still CONTACTS the chest wall (does not float), that
   breast volume is conserved, and that the same correction, fitted on two subjects,
   clears the third.

   **Tried, and no whole-breast correction works** (`scripts/fit_breast_seating_correction.py`;
   leave-one-subject-out, the judge above, fixed before fitting):

   | model | ribs inside the breast (<= 1%) | contact median, posterior 10% (<= 3 mm) | volume |
   |---|---|---|---|
   | one rigid forward shift per side (8-16 mm) | fails with s0790 held out (1.2-1.6%) | **fails every fold** (3.3-7.4 mm) | exact |
   | height-dependent push along the chest-wall normal | 0.00-1.51% | **fails every fold** (3.4-6.9 mm) | **fails every fold** (+7.8-9.2%) |

   Rib clearance and chest-wall contact pull against each other under any movement of
   the WHOLE breast: clearing the rib-5 overlap (up to ~16 mm) lifts the base 3-7 mm
   off the chest wall everywhere else, and the height profile buys that back with 8-9%
   more breast. Had gate d2 been the judge, the shift would have passed. Nothing seated
   was written. What is needed is a local deformation of the breast's BASE, conformed to
   this body's own chest wall with its volume compensated -- the breast-sized instance of
   the deformable soft-tissue decision this programme already has open. The fourth subject
   makes the case harder to avoid: a correction of one fixed size cannot cover an offset
   that varies twofold between women, so the base has to conform to each chest wall.

   **The base deformation: approved 2026-09-10, judged by rules fixed BEFORE it is built.**
   The judge written above for a fitted, systematic correction does not fit this. A
   deformation that places the breast's base on the chest wall clears the ribs and makes
   contact by construction, so neither can judge it. The judge, committed before any code:
   * **The solver first, on a case with a known answer.** A near-incompressible block
     (nu = 0.49) under prescribed compression conserves volume within 1%, and on the SAME
     mesh and boundary conditions the in-repo solver (`DeformableRegion`) and FEBio 4.13 --
     independent codes -- agree within 5% of the maximum displacement (RMS over nodes).
   * **Per breast, all four subjects (s0790, s1067, s1159, s0970):** (a) volume within 1%
     of the undeformed breast -- the solver only penalises volume change, so this is a real
     test; (b) no collapsed or inverted element, every tet's J > 0.2, the backend's own
     rejection threshold; (c) the in-repo solver and FEBio agree on the deformed breast
     within 5% of its maximum displacement.
   * **Reported, NOT judged, because the boundary conditions enforce them:** base contact
     (median <= 3 mm) and this body's rib points inside the breast (<= 1%).
   * **Stated up front.** The bed is the anterior surface of this body's pectoralis major,
     the breast's anatomical bed, not the ribs. With only prescribed displacements and no
     body force, Young's modulus cancels from the resting shape; only nu matters, taken as
     0.49 (adipose is nearly incompressible), an assumption. Gravity, and recovering the
     unloaded shape of a breast imaged supine under gravity, are out of scope here and named.

   **Result (673696f): the solver passes; closest-point seating is ill-posed and stopped at the
   gate on the first breast.** Solver: block at nu = 0.49 under 10% compression, volume -0.221%;
   `DeformableRegion` vs FEBio 0.00001% (compression) and 0.00009% (base on a cylinder) of max
   displacement; a mismatched-case control reads 124.8%, so the metric discriminates. Per breast,
   s1159-left: prescribing each of 2,449 penetrating base nodes to its closest point on the
   pectoralis bed is not one-to-one -- 867 of 4,241 base triangles flip and 2,945 are squeezed
   below 0.2 of their area BEFORE any solve, 7,904 of 12,621 base tets start at J <= 0, and FEBio
   stops at t = 0.227 after 104 retries. Gate (b) cannot be met under that boundary condition by
   any solver, so no breast was seated and the modulus-cancels check was not reached. The other
   seven breasts were meshed and not solved (s0970 would need up to 72 mm).

   **The sliding base: its judge, fixed 2026-09-10 before it is built.** The boundary condition
   the result points to is the anatomical one -- the breast is attached to the pectoral fascia
   along the normal and slides on it. So: every base node that penetrates the bed in the
   registered position is held ON the bed along the bed's normal and is free tangentially
   (frictionless, bilateral); every other base node may not enter the bed (unilateral); the rest
   of the surface is free. FEBio's counterpart is sliding contact against a rigid bed with
   tension allowed on the held nodes.
   * **The solver first, on the new boundary condition, with a known answer.** A block at
     nu = 0.49 compressed 10% between two frictionless platens deforms HOMOGENEOUSLY: every
     node's displacement is linear in its position and the lateral stretch is fixed by volume.
     Both solvers must reproduce the homogeneous field within 1% of the maximum displacement,
     and agree with each other within 5% on a block pressed onto a frictionless cylinder.
   * **Per breast, the same four subjects:** (a) volume within 1%; (b) every tet J > 0.2;
     (c) the two solvers agree within 5% of the maximum displacement; **(d) no base triangle
     flipped in the solved state** -- the failure the closest-point rule produced, judged
     because nothing in the new boundary condition prevents it by construction.
   * **Reported, not judged:** normal gap of the held nodes (enforced by the in-repo bound; a
     real number in FEBio's penalty contact), tangential slide per node, this body's rib points
     inside the breast, and whether Young's modulus still cancels (it should; checked by solving
     at two moduli).
   * **Stop rule:** a breast that fails (a)-(d) is a result and is recorded; the boundary
     condition is not changed and re-run on the same subjects.

   **Built: the solver gate passes in both codes; no breast is seated.**
   (`scripts/seat_breast_sliding.py`, `ihm/assembly/sliding_contact.py`,
   `ihm/assembly/febio_sliding.py`.) The in-repo side is projected Newton on
   `DeformableRegion`'s energy with the analytic neo-Hookean Hessian projected positive
   semi-definite per element, each base node's frame rotated to the bed normal so the normal is a
   bound, re-linearised as the node slides.

   | solver known answer | in-repo | FEBio | gate |
   |---|---:|---:|---|
   | nu=0.49 block, 10% between frictionless platens, vs the exact homogeneous field (lateral stretch 1.052925) | 0.00004% | 0.00002% | 1% -- **pass** |
   | block pressed onto a frictionless cylinder, the two codes against each other | RMS 1.58% of max displacement | 5% -- **pass** ||

   Supporting known answers: the Hessian against finite differences of the parent's own gradient
   (1.5e-10 relative), the rotation map against a per-node product (2e-16), and a mismatched-case
   control reading 124.8%, so the agreement metric can tell cases apart. FEBio terminates normally
   on the cylinder with no penetration left and a median held gap of 0.077 mm.

   | breast | (a) volume | (b) J > 0.2 | (c) two solvers | (d) base flips |
   |---|---|---|---|---|
   | s1159 left | not evaluated | not evaluated | not evaluated | not evaluated |
   | s0790, s1067, s0970 (both sides), s1159 right | not attempted | | | |

   **Neither solver completes a breast, so nothing is judged.** A breast starts up to 45 mm inside
   the bed. FEBio will not resolve that by contact: in place, 582 s without finishing one time step;
   with augmented Lagrangian (which closed the cylinder's gap and read 1.47%), the step retries are
   exhausted; with the bed retracted and advanced back, contact never engages. The in-repo solve
   advances about 0.1% of the seating per accepted step -- two elements invert whenever the
   increment grows -- so a full seating would need ~1000 steps of ~30-60 s. The other seven breasts
   were not attempted. Whether the sliding base is the right boundary condition is therefore still
   **untested on a breast**; what is tested is the solver, on the new boundary condition, in both
   codes.

   Four findings worth keeping. **The bed's three parts overlap:** taking each node's NEAREST bed
   point puts 34% of s1159-left's candidates on a sheet on the wrong side of them (median 24.6 mm
   away), so their measured gap flips sign as they slide -- tens of millimetres of apparent
   penetration through a constraint that forbids any. A capped ray along each node's own outward
   normal cannot pick the wrong sheet. **A re-linearisation cannot be judged by the held gap while
   nodes are deep:** at 41.4 mm depth a 1.36 deg turn of the interpolated normal moves the measured
   gap by 0.98 mm while the solve honours its bound to 0 um; the test is whether the linearisation
   settles, plus the gap itself at the end, where the lever arm is zero. **This FEBio build has only
   the skyline solver** (MKL, HYPRE and SUPERLU all OFF), so sliding contact runs with a symmetric
   approximation of its stiffness -- that changes the Newton matrix, not the equilibrium. **numpy 2
   scalars must never reach FEBio's XML:** `{x!r}` writes `np.float64(-0.01)`, which FEBio reads as
   no displacement at all, silently.

   **The bed misses a fifth of the breast base, and two attempts to overturn the 45 mm failed
   their own known answers (2026-09-10).** The judge names the bed as the anterior surface of
   pectoralis major. A breast's base is not confined to that muscle: inferolaterally it lies on
   serratus anterior and the external oblique aponeurosis. Measured, and needing no depth
   instrument -- the share of breast vertices with ANY pectoralis major within 15 mm of them in the
   plane perpendicular to the anterior axis:

   | breast | on pectoralis major | with serratus, obliques, rectus, ribs and sternum added |
   |---|---:|---:|
   | s1159 left | 87% | 99% |
   | s1159 right | 82% | 98% |
   | s0790 left | 78% | 90% |
   | s0970 left | 71% | 90% |

   So 13-29% of each breast has no bed under it at all, and a node there is measured against the
   muscle's EDGE. That is a candidate explanation for penetrations of 45 mm where the chest-wall
   offset this file already measures is about 10 mm (ribs 2-7, deepest at rib 5, max 13.9 mm).

   **It is only a candidate, and my own attempts to settle it failed.** I built two instruments for
   "how far behind the chest wall is this vertex" and both were refused by their own known answer:
   a local-median rule read p99 +190.9 mm testing the wall against itself once ribs were included
   (the wall wraps the body, so the median sits between front and back), and restricting to
   anterior-FACING vertices read median +39.9, p90 +173.4 mm, because the inner surface of a
   posterior rib faces anteriorly too. What passed its known answers is narrower: the pectoralis's
   own interior points read 100% inside and a point 500 mm lateral 0%, and by that test 27.5% of
   s1159-left lies inside the muscle's volume, at most 12.5 mm deep -- which is bounded by the
   sheet's thickness and therefore says nothing about tissue BEHIND it. **The subagent's 45 mm
   stands unrefuted**; what is established is the coverage gap above.

   **The bed is corrected before the next run, and the gates are not.** The judge above named
   pectoralis major because that is the breast's bed over most of its area; the measurement says
   it is not the bed over 13-29% of it. So the bed becomes **pectoralis major, pectoralis minor,
   serratus anterior, the external oblique and rectus abdominis** -- the muscular chest wall the
   breast actually rests on -- and still NOT the ribs, because a breast rests on muscle, and the
   rib overlap is a separate gate this file already reports. This is a modelling correction made
   after a failure and labelled as one: it changes what the tissue is seated ON, not what counts
   as seated. Gates (a)-(d) and the two-step placement rule are unchanged.

   **Predicted before it runs:** the maximum penetration falls from 45 mm to the order of the
   measured chest-wall offset, ~10-15 mm, and the count of base nodes with no bed on their ray
   falls from 672 of 2,850 to under 100. If the penetration stays above 40 mm with the full
   muscular wall under it, the depth is not a coverage artefact and the registration itself is
   what puts this tissue inside the chest.

   **Run on s1159 left. The coverage prediction holds; the depth prediction does not.** With the
   muscular chest wall as the bed, 4,097 of 4,301 posterior nodes have muscle on their ray -- only
   **204** do not, against 1,451 under pectoralis major alone. But the maximum penetration is
   **41.8 mm**, against 45.5 mm before: essentially unchanged. By the rule fixed above, the depth
   is **not a coverage artefact**, and the registration is putting this tissue inside the chest.

   **What the max hid, and what my own reporting hid with it.** The distribution at the registered
   pose, full muscular wall, 4,097 base nodes of which 3,123 penetrate:

   | percentile | 50 | 75 | 90 | 99 | max |
   |---|---:|---:|---:|---:|---:|
   | penetration | 7.4 mm | 12.5 mm | 17.4 mm | 24.7 mm | 41.8 mm |

   Deeper than 30 mm: **13 nodes (0.3%)**. Deeper than 40 mm: **one**. So the 45 mm that drove this
   whole line of work was always a handful of nodes, and the bulk of the base sits at 7-17 mm --
   the order of the chest-wall offset this file measures, which is **6-11 mm by a biased
   estimator and 5.5-6.6 mm by an independent one** (the caveat above; it was written here as
   "~10 mm" before either was checked). Every "45 mm" statement here, mine included, described one
   node as though it described the base.

   **Step 1 with both faults fixed still fails the 25 mm rule**, and the escape is not fully cured.
   The translation magnitude is now capped smoothly rather than boxed per component, and a node
   that loses its bed is counted at its pre-placement penetration. The minimiser still wants
   **73.1 mm** from the registered pose (32-80 mm across starts), and at that optimum the worst
   penetration is WORSE (41.8 -> 54.6 mm) with 858 nodes having lost their bed. Counting a lost node
   at its pre-placement value only charges escape up to what it already cost: a node whose
   penetration would otherwise have grown still profits by leaving. s1159 left is reported
   **unseated: registration failure**.

   **Step 2 was tried anyway from the registered pose, as a labelled diagnostic**, since the
   corrected bed leaves a median of only 7.4 mm. It still stalls: load fractions of ~0.002, minimum
   J down to 0.739, cut back at every attempt to grow the step. The sliding condition is not
   reachable for this breast by either route, so no gate table exists and the other seven were not
   attempted.

   One tool worth keeping: ray casting against the 125k-face muscular wall is spatially indexed in
   two tiers (the bed's face radii run 1.6 mm median but 28.9 mm max, and one oversized triangle
   would set the query radius for every ray), verified to reproduce brute force exactly and 21x
   faster.

   **Fourth attempt: fix the LABEL, not the threshold. Fixed 2026-09-10, before it runs.** The
   question the distribution raises is whether 13 nodes of 4,097 should stop a seating, and the
   answer is not to move a threshold after three runs have failed it -- it is to ask what those 13
   nodes are. **Breast tissue does not lie behind pectoralis major.** A CT label that does is a
   segmentation error, and this cohort's labels come from a model run on a clinical scan, not from
   a dissection. So the correction is to the DATA and it is gated:
   * **The trim.** Remove from the breast label every vertex lying behind the muscular chest wall
     by more than **20 mm** along its own outward ray, together with the tetrahedra they belong to,
     and re-mesh. 20 mm is fixed now, sits below the p99 of 24.7 mm so it is a real test rather
     than one satisfied by construction, and is above the ~10-17 mm the chest-wall offset explains.
   * **Gate on the trim itself, before any seating:** the removed volume is **<= 1%** of the
     breast. If it exceeds 1%, the label is not locally wrong -- the breast is in the wrong place,
     the subject is reported unseated as a registration failure, and no seating is attempted.
   * **Then drop step 1 entirely.** The rigid placement has failed twice on an objective that
     rewards carrying tissue off the muscle, and its escape is not cured by charging a departing
     node its pre-placement cost. With the corrected bed the median penetration is 7.4 mm and 95%
     of nodes have muscle on their ray, so there is nothing left for a rigid move to buy. Seat from
     the REGISTERED pose, which is also where the registration's own evidence puts the tissue.
   * **The four per-breast gates do not move**: volume within 1% of the *trimmed* undeformed
     breast, every tet J > 0.2, the two solvers within 5%, no flipped base triangle. The trimmed
     volume is reported next to them so nobody reads a 1% volume gate over a breast that lost 5%
     to the trim.

   **Predicted, before it runs:** the trim removes well under 1% of each breast (the 13 nodes are
   0.3% of the base and lie at its medial edge), the maximum penetration falls below 20 mm by
   construction, and the solve that stalled at load fractions of 0.002 completes. If the trim
   exceeds 1% on a subject, that subject's breast is misplaced rather than mislabelled, and saying
   so is the result.

   **Run on all eight breasts: every one fails the trim gate. All four subjects are registration
   failures, and no seating was attempted on any of them.** The prediction is refused -- the trim
   takes 2.0% to 7.7%, not "well under 1%" -- and by the rule fixed above that verdict is not about
   the labels but about where the registration puts them.

   | breast | vertices > 20 mm behind the wall | deepest | tets removed | volume removed |
   |---|---:|---:|---:|---:|
   | s1159 left | 301 | 41.8 mm | 1,603 | **1.99%** |
   | s1159 right | 314 | 33.6 mm | 1,832 | **2.26%** |
   | s1067 left | 643 | 46.4 mm | 3,469 | **3.70%** |
   | s0790 right | 833 | 45.3 mm | 4,526 | **4.73%** |
   | s0790 left | 1,042 | 47.7 mm | 5,636 | **4.86%** |
   | s1067 right | 1,218 | 45.1 mm | 6,545 | **5.25%** |
   | s0970 left | 1,645 | 42.6 mm | 8,895 | **7.56%** |
   | s0970 right | 1,605 | 40.8 mm | 8,416 | **7.73%** |

   The amount scales by subject -- s1159 ~2%, s1067 and s0790 ~4-5%, s0970 ~7.6% -- in the same
   order as the twofold spread in chest-wall offset this file already measured between these women.
   A segmentation error would not be systematic across four subjects, eight breasts and both sides,
   nor would it scale with the subject's own offset. **What is wrong is upstream of every seating
   attempt: one similarity transform per subject, fitted to the whole torso, cannot place the breast
   against this body's chest wall.** Four seating instruments have now been built and judged against
   it -- prescribed closest-point, the sliding base, rigid placement, and this trim -- and the
   solver gate passed in two independent codes while no breast was ever seated. The next instrument
   belongs on the REGISTRATION, not on the tissue: a fit that is local to the chest wall, or a
   per-subject correction with its own known answer, before any seating is judged again.

   **The chest-wall registration, and its gates, fixed 2026-09-10 before it is built.** Accepted:
   the seating line stops here, and the next instrument is the fit. Her CT carries the bones to do
   it -- ribs 1-12 per side, sternum, clavicles and T1-T12 as separate labels -- so the fit does not
   need her soft tissue, which is what the whole-torso similarity was leaning on.

   * **The instrument.** One similarity per subject fitted to the CHEST WALL alone: her ribs 2-7,
     both sides, onto this body's, by the one-way trimmed ICP the pelvic registration already uses
     (partial scan surfaces onto complete body bones, never the reverse). Ribs 2-7 because that is
     the breast's own base, measured: the rib points inside the mapped breasts are ribs 2-7,
     deepest at rib 5.
   * **Gate 1, known answer.** This body's own chest bones, truncated the way her scan's field of
     view cuts them, moved by a known similarity, recovered within **2 mm, 2 deg and 1%** -- the
     pelvic registration's gate, which it passes at 0.033 mm.
   * **Gate 2, laterality.** Her left ribs map nearer this body's left ribs than its right.
   * **Gate 3, held out.** Ribs 2-7 are fitted; **sternum, clavicles, rib 1, ribs 8-12 and T1-T12
     are held out**, and on those the median nearest-surface distance must be **no worse than the
     whole-torso similarity achieves on the same structures**, subject by subject. A fit judged on
     what it was fitted to would only report its own objective; this asks whether a chest-local fit
     costs anything elsewhere.
   * **Gate 4, the consequence.** With the new transform, the breast volume more than 20 mm behind
     the muscular chest wall is **<= 1%** on every subject -- the trim gate above, now testing the
     registration that failed it rather than the label that did not.
   * **Reported, not judged:** how far each subject's chest-local transform differs from its
     whole-torso one, and the chest-wall offset table recomputed under it.

   **Predicted before it runs:** gate 4 passes for s1159 and s1067, and s0970 is the one at risk --
   its 7.6% is twice s1159's and it carries the largest offset of the four. If s0970 still exceeds
   1% under a chest-local similarity, then its chest is a different SHAPE and not a different size,
   no similarity will place it, and the honest next step is a deformable chest-wall fit or dropping
   that subject with the reason stated.

   **Built and run on all four (`scripts/register_female_chest_wall.py`). The instrument works; the
   chest-local similarity improves every breast and still fails, on the held-out gate, for every
   subject.** The prediction is half refused: gate 4 passes for s1159 alone, not s1067.

   | subject | 1 known answer | 2 laterality | 3 held out (chest-local vs whole-torso) | 4 behind the wall, left / right (was) | overall |
   |---|---|---|---|---|---|
   | s0790 | pass, 0.218 mm | pass | **FAIL** 5.81 vs 4.55 mm | 3.28% / 4.01% (4.86 / 4.73) | **FAIL** |
   | s1067 | pass, 0.143 mm | pass | **FAIL** 6.47 vs 5.18 mm | 1.76% / 3.24% (3.70 / 5.25) | **FAIL** |
   | s1159 | pass, 0.023 mm | pass | **FAIL** 8.24 vs 5.52 mm | **0.75% / 0.06%** (1.99 / 2.26) pass | **FAIL** |
   | s0970 | pass, 0.051 mm | pass | **FAIL** 8.70 vs 5.60 mm | 2.75% / 3.67% (7.56 / 7.73) | **FAIL** |

   Gate 1 holds to the pelvic standard -- 0.023-0.218 mm, 0.029-0.345 deg, <= 0.012% scale -- with
   this body's ribs truncated to the 0.76-0.89 of each rib her labels actually carry. Laterality is
   unambiguous: her ribs land 3-4 mm from this body's ribs on their own side and 123-127 mm from the
   other. The fit itself is good where it was fitted: trimmed residual 3.6-6.2 mm on ribs 2-7.

   **The chest-local fit helps the breast everywhere.** Volume more than 20 mm behind the muscular
   wall falls on all eight breasts, most of all where it was worst: s0970 7.6/7.7% -> 2.8/3.7%,
   s1159 2.0/2.3% -> 0.8/0.1%. So the diagnosis was right -- the breast sat inside the chest because
   of the fit, not because of the label.

   **And it costs the rest of the torso, which is what the held-out gate was for.** On the 27 held-out
   structures the median nearest-surface distance is WORSE than the whole-torso similarity for every
   subject, and the damage is concentrated in the lower ribcage: ribs 9-12 move 17-27 mm further from
   this body's own. The chest-local transform also scales her systematically less than the whole-torso
   one (s1067 1.032 vs 1.098; s1159 1.009 vs 1.060) and moves her breast centroid 8.6-19.2 mm.

   **So the conclusion the prediction reserved for s0970 holds for all four: the difference is SHAPE,
   not size or pose.** One similarity cannot serve the chest wall and the rest of the torso at once --
   what it buys at the breast it loses at the lower ribs -- and that is now measured rather than
   argued. A deformable chest-wall fit, or per-region fits with an explicit blend and its own
   held-out gate, is the next instrument; a fifth similarity is not.

   One note on the known answer, because two versions of it were wrong before this one. Cutting this
   body's ribs to her field of view's superior-inferior band removed NOTHING (her 369 mm band spans
   ribs 2-7 whole, and it is her segmentation, not her field of view, that leaves them partial).
   Cutting them to what her mapped label lies within removed far too much (a median 0.27-0.46 of each
   rib) because it measures the whole-torso misalignment under investigation. The cut that is neither
   is calibrated by AREA RATIO, which depends on the scale and not the pose, taken from the anterior
   costal-cartilage end where a CT rib label stops. **If her labels are instead thinned all over
   rather than cut at that end, gate 1 is easier than it looks** -- said here because the agent that
   built it said it, and nothing downstream should quote 0.023 mm without it.

   **The deformable chest-wall fit, pre-registered 2026-09-10 before it is built -- and it is the
   instrument the OTHER line already built.** `docs/SEGMENT_CONTACT_SURFACES.md` reached the same
   wall from the other side: one global similarity cannot carry this specimen's skin onto the
   scaffold either, and a smooth warp (`scripts/skin_warp.py`, thin-plate spline with a fixed
   cross-validation rule, a zero-warp control and a matrix-free solver verified against the direct
   one) is already built and gated there. The same instrument answers this, so it is reused rather
   than reinvented.

   * **The instrument.** `W(x) = Gx + d(Gx)`, G the whole-torso similarity this file already
     records, d a regularised thin-plate spline fitted to correspondences from HER chest bones --
     ribs 2-7 both sides -- onto this body's, exactly the set the chest-local similarity used.
   * **Gate 1, known answer, and not the one that needed three tries.** A warp is recovered from a
     warp: take this body's own chest bones, displace them by a KNOWN smooth field of the same
     family, and the fit must recover it to **1 mm RMS** over the fitted bones. This tests the
     instrument without depending on any assumption about where her labels stop.
   * **Gate 2, laterality.** Unchanged.
   * **Gate 3, held out, unchanged and unweakened.** Fit ribs 2-7; sternum, clavicles, rib 1, ribs
     8-12 and T1-T12 held out, and the median nearest-surface distance there **no worse than the
     whole-torso similarity's** -- 4.55, 5.18, 5.52 and 5.60 mm for s0790, s1067, s1159 and s0970.
     The chest-local similarity failed this at 5.81, 6.47, 8.24 and 8.70; a warp that trades the
     lower ribs for the breast fails it too.
   * **Gate 4, the consequence.** Breast volume more than 20 mm behind the muscular wall **<= 1%**
     on every subject, as before.
   * **Reported:** the warp's bending energy, its displacement at the breast base and at ribs 9-12,
     and the chest-wall offset table recomputed.

   **Predicted before it runs:** gate 3 passes, because a warp has the freedom a similarity lacks
   and the held-out ribs only need to stay where they already are; gate 4 passes on s1159 and s1067
   and is uncertain on s0790 and s0970. If gate 3 fails again, the constraint is not the
   transform's freedom but the correspondences -- her ribs are segmented from a CT and this body's
   are an atlas, and two rib surfaces that do not mean the same thing cannot be made to agree by
   any map.

   **Gate 1 fails, no subject is fitted, and the reason is the correspondences -- shown on this
   body's OWN bones, with no segmentation mismatch anywhere near it**
   (`scripts/fit_chest_wall_warp.py`, importing `scripts/skin_warp.py` unchanged).

   | gate 1, a warp recovered from a warp | result | bar |
   |---|---:|---:|
   | surface recovery, random smooth field (3.1 mm: 1.4 normal, 2.6 tangential) | 1.194 mm RMS | 1 mm -- **FAIL** |
   | pointwise recovery, same field | 3.062 mm RMS | -- |
   | surface recovery, purely NORMAL control field (1.5 mm mean, 4.0 max) | 1.274 mm RMS | 1 mm -- **FAIL** |
   | the warp's own fit residual at its correspondences (normal field) | **0.01 mm** | -- |

   Two readings of "recover it to 1 mm RMS" were measured rather than one chosen: the field
   POINTWISE, and the SURFACES. They differ because nearest-point correspondences cannot see motion
   ALONG a surface -- 2.6 of the 3.1 mm here is tangential and unidentifiable by any surface-matching
   transform. That alone would have been a weak excuse, so it was tested: a field applied purely
   along the surface normals, identifiable by construction, recovers no better (1.274 mm).

   **What it is.** The warp reproduces its targets to 0.01 mm and still lands 1.27 mm from the true
   surface, so the targets are what is wrong. Measured directly against the known true preimage, the
   nearest-point target is off by a mean of 0.671 mm, p90 2.281 mm and max 7.152 mm for displacements
   averaging 1.53 mm -- the bias of order d^2/R that nearest-point matching carries on a curved
   surface, here ribs a few millimetres thick. Its RMS is the 1.27 mm the gate measured.

   **So the constraint named as the fallback is real, and it is not about her data.** It was reached
   with this body's own bones on both sides of the fit and a field of the fitted family; no CT, no
   atlas mismatch, no segmentation. Nearest-point correspondence cannot support a 1 mm bar on these
   surfaces, and no increase in the transform's freedom changes that -- a warp fitted to biased
   targets reproduces the bias. The next instrument is a CORRESPONDENCE, not a transform: arc-length
   or landmark parameterisation along each rib, normal shooting, or a symmetric matching that is
   unbiased on curvature, with this same known-answer gate to prove it before any subject is fitted.

   **Accepted, and the correspondence is pre-registered here before it is built (2026-09-10).**
   Failing at gate 1 is worth more than failing at gate 4 would have been: it was reached with this
   body's own bones on both sides and a field of the fitted family, so nothing about her CT, her
   segmentation or the atlas is implicated, and no better warp built on these targets can be
   expected to pass. Both readings of the gate were measured rather than the convenient one, and
   the purely normal control rules out tangential blindness as the explanation. That is the
   standard the rest of this line is held to.

   * **The instrument: symmetric normal shooting.** For each source sample, shoot along the SOURCE
     surface's own outward normal to the target surface and take the first intersection within a
     cap; shoot back from that hit along the TARGET's normal; keep the pair only where the return
     lands within 1 mm of the source. Direction then comes from the surfaces rather than from
     proximity, which is what carries the d²/R bias.
   * **Gate 0, the correspondence itself, before any warp is fitted.** Against the known preimage
     of a known smooth field, the target error must be **<= 0.2 mm mean and <= 0.5 mm p90** -- an
     order below the 0.671 mm and 2.281 mm the nearest-point rule gives on the same test. A
     correspondence is judged as an instrument before anything is fitted with it.
   * **Gate 1 unchanged:** the warp recovers a known field to **1 mm RMS**, reported both pointwise
     and to the surface, exactly as it was failed.
   * **Gates 2-4 unchanged**, and only reached if gate 0 and gate 1 pass.
   * **Reported:** how many samples are dropped for failing the return test, and where.

   **Predicted before it runs:** normal shooting brings the target error under 0.2 mm and the
   recovery under 1 mm. If it does not, the limit is the surfaces themselves -- a rib is a few
   millimetres thick and these meshes are not dense -- and the 1 mm bar would have to be re-derived
   from what the geometry can support, which is a new pre-registration and not a loosening of this
   one.

   **Gate 0 fails as written, on the mean alone, and one line fixes it -- proposed, not adopted**
   (`ihm/anatomy/normal_shooting.py`, written as a reusable instrument because the skin line needs
   it too; `scripts/fit_chest_wall_warp.py` uses it).

   | target error vs the known preimage | mean | p90 | max | kept |
   |---|---:|---:|---:|---:|
   | nearest point, the rule that failed gate 1 | 0.671 mm | 2.281 mm | 7.152 mm | -- |
   | **normal shooting, exactly as pre-registered** | **0.299 mm** | 0.194 mm | 7.142 mm | 4,701 / 4,800 |
   | normal shooting + normal agreement (NOT the gate) | 0.078 mm | 0.157 mm | 5.074 mm | 4,382 / 4,800 |
   | the bars | <= 0.2 | <= 0.5 | | |

   The instrument itself is sound: a closed-form check on a subdivided sphere, where the true partner
   of a point on a 2%-inflated copy is its radial point, recovers it to a mean 4.9e-4 of the radius
   with nothing dropped. On the ribs, p90 beats its bar by 2.6x and beats nearest point by 12x. The
   MEAN fails because of a tail, and the tail is one thing: a normal shot at a rib a few millimetres
   thick CROSSES it and lands on the far cortical wall, where the surface faces the other way and the
   return test still passes because the two walls are roughly parallel. On right rib 2 that is 99 of
   362 kept pairs, normal agreement -0.850, median error 4.52 mm; the failure is concentrated in
   right ribs 2 and 3 (means 1.548 and 1.396 mm) while the other ten sit near 0.12 mm.

   **The fix is to require the target's normal to agree with the source's** -- a far-side hit is
   exactly a disagreement -- which drops 319 of 4,800 samples and gives 0.078 mm mean, 0.157 mm p90,
   inside both bars. It is implemented and OFF by default, because it ADDS a filter to an instrument
   whose pre-registration specifies shooting and the return test, and a gate is not something to
   quietly re-specify after seeing its number. It needs a pre-registration of its own; the numbers
   above are what it would be judged against.

   So no warp was fitted and gates 1-4 were not reached. What is established: normal shooting removes
   the d^2/R bias that made nearest-point targets unusable (0.671 -> 0.299 mm mean, 2.281 -> 0.194 mm
   p90, both without any filter), and what remains is a thin-sheet failure with a known cause, a
   measured size, and a one-line remedy.

   **With agreement adopted and turned on, gate 0 is clean and GATE 1 STILL FAILS. The correspondence
   is at its discretisation floor; the warp is not what the correspondence limits.**

   | | mean | p90 | worst rib |
   |---|---:|---:|---|
   | gate 0, agreement on (REPORTED, not claimed -- the rule was adopted on its argument after these were seen) | 0.078 mm | 0.157 mm | right rib 3, 0.133 mm |
   | gate 0, agreement off, as pre-registered | 0.299 mm | 0.194 mm | right rib 2, 1.548 mm |
   | nearest point | 0.671 mm | 2.281 mm | -- |

   The far-side hits are gone: the ten good ribs were near 0.12 mm before and the two bad ones at
   1.548 and 1.396; now every rib is between 0.10 and 0.133 mm. Of 4,800 samples, 4,382 are kept --
   319 dropped for disagreeing normals, 98 for a return further than 1 mm, 1 for no hit.

   | gate 1, a warp recovered from a known field (3.1 mm: 1.4 normal, 2.6 tangential) | surface | pointwise |
   |---|---:|---:|
   | with the unbiased correspondence (agreement on) | **1.512 mm RMS** | 3.147 mm |
   | with the biased nearest-point targets, for comparison | 1.194 mm RMS | 3.062 mm |
   | the bar | 1 mm | -- |

   **An 8.6x more accurate correspondence made the recovery slightly WORSE, which is the finding.**
   The warp's cross-validated error at held-out correspondences is 0.434 mm (0.735 in the first pass)
   and is FLAT from lambda 1e-6 to 1e-4 -- an almost-interpolating fit predicts an unseen point no
   better than a smoothed one. So the limit is not regularisation and not correspondence accuracy: it
   is that the field between samples is not predictable from the samples. In the fit, 3,437 of 3,600
   samples were kept (0 disagreeing normals, 90 returns too far, 63 no hit), spaced about 3.6 mm.

   **What the geometry can support, measured.** On spheres at three resolutions the correspondence
   error tracks the faceting exactly -- 0.135 x edge^2, ratio 0.134 / 0.135 / 0.137 as the edge
   halves -- so it is discretisation, not the instrument. This body's ribs 2-7 have a mean edge of
   1.94 mm (median 1.44, p90 4.05), and at rib-scale curvature that law predicts about 0.07 mm
   against the 0.078 mm measured. The correspondence is therefore AT its floor and cannot improve
   without denser bone meshes, while the warp sits 20x above it.

   **One caveat on gate 1 that belongs with any re-derivation of its bar:** the known field's spatial
   scale is an implementer's choice that was never pre-registered -- 720 anchors over ribs 2-7 with
   ~4 mm random weights, which varies at roughly a centimetre. A field that varied at anatomical
   scales would be easier to recover from 3.6 mm samples, and a finer one harder. Whatever bar
   replaces 1 mm should fix that scale explicitly, or gate 1 measures the test field as much as the
   instrument.

   **The bar is replaced by a curve, and the question it answers is fixed first (2026-09-10).**
   The caveat is correct and it is mine: 1 mm was a number, not a requirement derived from anything.
   A single bar against a field of unstated scale cannot say whether this pipeline can carry a
   chest wall, and an 8.6x better correspondence scoring slightly worse is the proof that the bar
   was measuring the wrong thing. So gate 1 is replaced, not loosened:

   * **Gate 1a, the resolution curve.** Recover known smooth fields of the fitted family at spatial
     scales **L = 5, 10, 20 and 40 mm**, fixed here, five seeds each, amplitude held at the 3.1 mm
     already used, reported as surface RMS against L. The instrument's resolution **L\*** is the
     smallest L recovered to <= 1 mm. Run at the present correspondence density and at **4x** it,
     because sampling and the warp family are different limits and the curve separates them.
   * **Gate 1b, the question that actually matters.** Measure the spatial scale of the REAL
     deformation -- the per-subject offset field between her registered ribs and this body's, whose
     magnitude this file already records at about 10 mm -- as the distance over which that field
     decorrelates to 1/e. The pipeline can carry this deformation if that scale **exceeds L\***.
   * **Reported beside both:** the correspondence's discretisation floor (0.078 mm, 0.135 x edge^2
     at a 1.94 mm mean edge) and the warp's held-out correspondence error (0.434 mm, flat from
     lambda 1e-6 to 1e-4), so a reader can see which of the three limits binds.
   * **Gates 2-4 unchanged**, and reached only if 1b says the deformation is representable.

   **Predicted before it runs:** L\* lands between 10 and 20 mm at the present density and improves
   less than a factor of two at 4x, because a 0.434 mm held-out error that is flat in lambda is a
   sampling limit rather than a smoothing one; and the real chest-wall field decorrelates over
   several centimetres, so 1b passes and the line continues. If instead L\* is above the real
   field's scale, no warp on these meshes can carry a chest wall, and the honest end of this line
   is to say so rather than to fit one.

   **The curve is FLAT, L\* does not exist below 40 mm, and 4x the sampling changes nothing. This is
   the end of the line, and not for the reason predicted.**

   | surface RMS, mean of 5 seeds, amplitude fixed at 3.1 mm | L = 5 mm | 10 mm | 20 mm | 40 mm |
   |---|---:|---:|---:|---:|
   | present density (~3,300 correspondences, ~3.6 mm apart) | 1.53 | 1.58 | 1.47 | 1.52 |
   | 4x density, paired seed 0 at L = 5 (13,496 correspondences, ~1.8 mm apart) | 1.595 vs 1.600 | | | |

   **L\* is above 40 mm at both densities.** A field varying over 40 mm recovers no better than one
   varying over 5 mm, and quadrupling the correspondences improves the paired seed by 0.3%. So the
   limit is neither the field's scale nor the sampling -- the two things the curve was built to
   separate -- and the prediction of L\* between 10 and 20 mm is refused in a way that makes the
   question moot.

   **What binds, measured.** The warp fits its correspondences to 0.012 mm and is already 0.818 mm
   off 1-2 mm away from one, 1.309 mm at 2-4 mm, 1.652 mm at 4-8 mm, against a 0.288 mm floor for
   the measure itself (an identity warp on coinciding surfaces). The reason is the surface: on this
   body's ribs the normal turns a median 4.6 deg between samples 1-2 mm apart and 7.9 deg at 2-4 mm,
   and a 3.1 mm displacement along normals 7.9 deg apart differs by 0.43 mm in direction. Normal
   shooting returns correct target POINTS, but the displacement VECTORS it returns inherit that
   turning, so they are not the restriction of any smooth 3D field, and a thin-plate spline chasing
   them oscillates between its own data. The three limits reported beside it -- the correspondence's
   0.078 mm discretisation floor, the warp's 0.434 mm held-out error, the 0.288 mm measure floor --
   are all an order below what actually binds.

   **Gate 1b, measured without the correspondence's filters** (the return and agreement tests keep
   only pairs that already nearly coincide -- 3-7% here -- which biases both the magnitude and the
   scale): her ribs stand 4.1-5.4 mm mean (2.7-3.7 median) from this body's under the whole-torso
   similarity, and that offset field decorrelates to 1/e within **2-6 mm** on all four subjects.
   Against an L\* above 40 mm, the deformation is an order finer than anything this instrument can
   represent. **Caveat, and it cuts both ways:** the offsets are measured along each sample's own
   normal, so this number inherits the same normal-turning contamination that limits the warp. What
   is certain is that the two are not separated by the margin the gate needed.

   **So: no warp of this family, fitted on these meshes, can carry this chest wall**, and the
   honest end is to say so. One avenue remains and is deliberately NOT taken here, because it would
   be a new instrument and a new pre-registration: the obstacle is the roughness of the target
   normals, so fitting the normal component as a SCALAR field on the surface, or smoothing the
   target surface before shooting, attacks the thing that binds rather than the two that do not.

   **The scalar offset field: pre-registered 2026-09-10, before it is built.** Taken, because the
   diagnosis names what to change and the change is a different REPRESENTATION rather than another
   transform. A vector field in 3D has to reproduce the normal's turning to describe a deformation
   that is mostly a difference in how far out the surface sits; a scalar does not.

   * **The instrument.** The deformation is `d(x) = s(x) n_smooth(x)` on this body's own rib
     surface: `s` a smooth scalar field fitted over the surface, `n_smooth` the normal of that
     surface after **20 Taubin iterations** (volume-preserving, fixed here, reported), which is the
     roughness that binds. A displacement is then one number per point, and the turning cannot
     enter it.
   * **Gate A, known answer.** Displace this body's own ribs by a KNOWN scalar field along those
     same smoothed normals, at the same amplitude (3.1 mm) and the same four scales
     (L = 5, 10, 20, 40 mm), five seeds. Recover to **<= 0.5 mm surface RMS** at every L. Half the
     old bar, because a representation that matches the deformation should beat a general one, and
     a scalar field that cannot do this is not worth fitting to a subject.
   * **Gate B, what the representation gives up, measured before it is used.** The real deformation
     is not purely normal. Decompose the unfiltered offset field between her ribs and this body's
     into normal and tangential parts and report the tangential fraction. **If tangential exceeds
     50% of the total, a normal-only model is the wrong representation and that is the result** --
     stated now so it cannot be discovered later and explained away.
   * **Gates 2-4 unchanged**, reached only if A passes and B permits.
   * **Reported:** how much the Taubin smoothing moves the surface it smooths, since that motion is
     an error the fit will not see.

   **Predicted:** gate A passes at every L, because the field being fitted is exactly the field
   being represented; gate B comes back tangential-minority at the breast base, where a chest wall
   differs mostly in how far forward it sits. If gate B comes back tangential-majority, the line
   ends for good and the female chest wall needs a different body, not a different fit.

   **Gate B passes decisively; gate A fails at every L. The representation is right and the
   instrument still misses by about a factor of two.**

   | gate A, surface RMS, 5 seeds, amplitude 3.1 mm | L = 5 | 10 | 20 | 40 mm |
   |---|---:|---:|---:|---:|
   | scalar field along smoothed normals | 1.046 | 1.064 | 0.942 | **0.897** |
   | the 3D vector warp, for comparison | 1.53 | 1.58 | 1.47 | 1.52 |
   | the bar | 0.5 | 0.5 | 0.5 | 0.5 |

   Worst over all L and seeds: 1.430 mm. So gate A **fails**, and gates 2-4 are not reached. But the
   scalar field is the better representation by the margin its argument predicted: it takes 1.5 mm to
   about 1.0, and unlike the vector warp it now IMPROVES with scale (0.897 at L = 40 against 1.046 at
   L = 5), which is what a representation that matches the deformation should do.

   | gate B, the real offset field decomposed against this body's smoothed normals | normal | tangential | tangential fraction |
   |---|---:|---:|---:|
   | s0790 | 6.49 mm | 0.54 mm | 7.7% (energy 3.1%) |
   | s1067 | 10.87 mm | 1.12 mm | 9.3% (energy 2.9%) |
   | s1159 | 10.06 mm | 1.02 mm | 9.2% (energy 2.2%) |
   | s0970 | 9.27 mm | 0.62 mm | 6.3% (energy 0.7%) |

   **Gates 2-4 on all four subjects: all fail, and a 10% representation error is not what stood in
   the way. The correspondence collapses on real data.**

   | subject | 2 laterality | 3 held out (bar) | 4 behind the wall, left / right (whole-torso similarity) | overall |
   |---|---|---:|---:|---|
   | s0790 | pass | **4.59** (4.55) FAIL | 6.66% / 5.51% (4.86 / 4.73) | **FAIL** |
   | s1067 | pass | **5.59** (5.18) FAIL | 3.46% / 3.73% (3.70 / 5.25) | **FAIL** |
   | s1159 | pass | 4.77 (5.52) pass | 2.31% / **0.00%** (1.99 / 2.26) | **FAIL** |
   | s0970 | pass | 5.16 (5.60) pass | 7.81% / 7.60% (7.56 / 7.73) | **FAIL** |

   **Why, in one row.** The fit kept 134 to 311 of 3,600 samples per subject -- 3.7% to 8.6% -- with
   2,100 to 2,700 of them finding NO target within the 20 mm cap. Ribs are thin and curved and stand
   5-10 mm apart, so a ray along her rib's own normal passes BESIDE this body's rib rather than
   hitting it. The pairs that do survive are the already-aligned minority: the offsets they see have
   a median of 0.37-0.45 mm where the real deformation is 6.5-10.9 mm. The scalar field is therefore
   fitted to about a twentieth of the deformation it was meant to carry, and it behaves like it --
   gate 3 barely moves in either direction (the largest damage anywhere is +0.6 mm, on T7), and
   gate 4 improves on some breasts (s1159 right 2.26% -> 0.00%, s1067 5.25% -> 3.73%) while making
   others worse (s0790 4.86% -> 6.66%).

   This is the same selection bias that had to be corrected in gate 1b, found this time INSIDE the
   fit: a correspondence with a return test keeps what already agrees. It passed gate 0 at 0.078 mm
   because there the two surfaces were the same surface displaced; it does not survive two rib cages
   that stand a centimetre apart.

   **So the standing 10% criterion is not what to test next, and neither is the representation.**
   Gate B settled the representation (91-94% normal) and gate A's bar was unreachable by
   construction; what blocks the line now is getting a correspondence at all across a centimetre of
   separation on thin bones. Two carried forward as asked: the Taubin smoothing's **5.501 mm max**
   displacement is the kind of local error that would show in held-out ribs, and it does not appear
   there because the field is too small to expose it; and `tps_fit_direct` silently broadcast a
   one-column scalar into three identical columns until it was caught -- a fault that returns a
   plausible wrong answer rather than an error, which is what makes the known-answer gates worth
   their cost.

   **Stop corresponding ribs. Pre-registered 2026-09-10, before it is built.** The whole line has
   been trying to match two rib cages bone to bone, and a rib is a few millimetres thick with 5-10 mm
   of air beside it: a ray that misses by 3 mm finds nothing, which is why 2,100 of 3,600 samples
   found no target at all. But a breast does not rest on a rib. It rests on the CHEST WALL -- a
   continuous sheet of muscle over bone over intercostal -- and that surface has no gaps to miss.

   * **The instrument.** For each body, build the chest-wall ENVELOPE over the breast base: the
     outward surface of the union of ribs 2-7, their costal cartilage, the sternum and the muscular
     wall already used as the seating bed, closed by a morphological closing with a **10 mm** ball
     (fixed here) so that inter-rib gaps are spanned rather than entered, and smoothed by the same
     20 Taubin iterations. Correspondence is normal shooting between the two envelopes; the scalar
     offset field of `de88ef4` is then fitted on the envelope, where it was always the right object.
   * **Gate 0', the correspondence at a REAL separation -- the lesson this line paid for.** The old
     gate 0 scored a surface against itself displaced, which is why 0.078 mm meant nothing at a
     centimetre. The new one scores it between THIS body's envelope and each subject's, and it is
     a coverage gate first: **>= 95% of samples must find a target within 20 mm**, against the
     3.7-8.6% the rib correspondence managed. Only then is accuracy read, on the synthetic case, at
     a separation of 10 mm rather than zero.
   * **Gate A', accuracy.** The standing relative criterion: recovery residual <= 10% of the normal
     deformation carried.
   * **Gates 2-4 unchanged.**
   * **Reported:** how far the closing moves the surface, since a 10 mm ball spanning a 10 mm gap
     invents surface where the body has none, and that invention is exactly where the breast sits.

   **Gate 0' FAILS: coverage 13.5-22.1%, against >= 95%.** Better than the rib correspondence's
   3.7-8.6%, and nowhere near enough, so no field was fitted and gates A', 2-4 were not reached.

   | subject | coverage | no hit | no return | return too far | normals disagreed | offsets seen |
   |---|---:|---:|---:|---:|---:|---:|
   | s0790 | 13.5% | 1189 | 161 | 1456 | 308 | 1.23 mm median |
   | s1067 | 16.1% | 983 | 139 | 1585 | 315 | 1.01 mm |
   | s1159 | 22.1% | 848 | 125 | 1551 | 282 | 1.14 mm |
   | s0970 | 15.4% | 1069 | 165 | 1565 | 247 | 1.11 mm |

   **Two reasons, and the second is about the instrument rather than the bodies.** A quarter of rays
   miss the other envelope entirely and a wider cap does not help (75.8% hit at 20 mm, 77.0% at
   40 mm): the two envelopes span different extents, so rays near their edges leave the far surface
   altogether. And the round-trip test rejects 40-44% more; loosening it tenfold, to 10 mm, raises
   coverage only to 52.1% (diagnostic, not the gate). At a real separation the return has a lever
   arm -- two surfaces whose normals differ by 10 deg over an 8 mm gap miss the round trip by
   1.4 mm -- so a 1 mm tolerance set at zero separation is doing something different here.

   **And the envelope shrinks what it was built to carry.** The two envelopes stand a median 3.0 mm
   apart where the rib-to-rib deformation is 6.5-10.9 mm, and the offsets the surviving pairs see are
   1.0-1.2 mm. A 10 mm ball closing both walls into their own inter-rib gaps inflates each toward the
   other, so the very offset the scalar field is meant to represent is largely absorbed by the
   closing. **Reported as required:** the closing moves this body's wall 1.85 mm median, 2.14 mean,
   8.45 mm max, and **6.7% of the envelope stands more than 5 mm from any real surface** -- invented
   surface, in the place the breast sits.

   **Deviation, recorded.** The envelope is built from ribs 2-7 and the sternum on BOTH sides, not
   from the muscular wall: neither body has a costal-cartilage label, and her CT has no muscle labels
   at all. Including this body's muscle on one side only would bias every offset outward by the
   muscle's thickness, flattering the breast gate and corrupting the held-out one.

   **A known-answer lesson, paid for again.** `igl.marching_cubes` enumerates its grid with x
   FASTEST; C-ordered values scramble against their positions and shatter the surface -- 9,743
   components here. A sphere could not catch it, because a sphere is invariant under axis
   permutation and the wrong ordering still returns a sphere of the right radius. The check that
   caught it uses an ellipsoid with three distinct semi-axes. **A known answer must break the
   symmetry it is testing**, the same way a correspondence's must be scored at the separation it
   will be used at.

   **Stop looking for a correspondence at all. Pre-registered 2026-09-10, before it is built.**
   Every instrument in this line has failed at the same place: pairing a point on one body with a
   point on the other. Nearest point is biased, normal shooting misses thin bone, the envelope that
   removes the misses also removes the deformation, and the return test that guards any of them
   rejects honest pairs once the surfaces are a centimetre apart. The closing radius was mine and
   the trade-off is structural: anything that makes two chest walls matchable erases what
   distinguishes them.

   A scalar offset field does not need a correspondence. For a point `x` on THIS body's chest wall
   with smoothed normal `n(x)`, the offset is **where her wall is along that ray**, read from a
   SIGNED DISTANCE FIELD of her chest wall rather than from a partner point: `s(x)` is the first
   zero crossing of her SDF along `n(x)` within a 25 mm cap. An SDF is defined everywhere, so
   nothing misses, nothing needs a round trip, and no geometry is closed or smoothed to make the
   two surfaces comparable -- the raw ribs, cartilage and sternum stay exactly as they are.

   * **Gate 0'', the field carries the deformation.** The median `|s|` over the breast base must
     match the rib-to-rib offset this file already measures, **6.5-10.9 mm, within 20%** -- and
     NOT the 1.0-1.2 mm the filtered correspondences saw, nor the 3.0 mm the envelopes stood apart.
     This is the gate the last three instruments would have failed, and it is first because it is
     the one that matters.
   * **Gate 0''', coverage.** >= 95% of base points must find a crossing within the cap. An SDF has
     no gaps, so a failure here is a statement about the ribs' extent rather than the instrument.
   * **Known answer, with the symmetry broken as this line has now twice had to learn.** An
     ellipsoid with three distinct semi-axes, displaced by a known scalar field along its normals
     at a **10 mm** separation, recovered to <= 10% of the displacement. Not a sphere.
   * **Then:** the field is smoothed as a scalar ON the surface -- never the geometry -- by the
     committed CV rule, and **the smoothing must be reported as a fraction of the deformation it
     removes**. Gates A' and 2-4 unchanged.

   **Predicted:** gate 0'' passes, because reading a distance field along a normal cannot select
   for agreement the way a matched pair does; coverage clears at the cap; and the smoothing removes
   under 20% of the deformation. If gate 0'' fails -- if even an SDF read along the normal sees
   1 mm where the bones stand 8 mm apart -- then the 6.5-10.9 mm figure is itself wrong, and what
   needs re-examining is the offset measurement rather than any instrument built on it.

   **Predicted:** coverage clears 95% easily, because a closed envelope has no gaps to miss;
   the offsets the correspondence sees then match the 6.5-10.9 mm the decomposition measured,
   rather than the 0.37-0.45 mm the rib pairs saw. If coverage clears and gate 4 still fails, the
   obstacle is not correspondence either, and the remaining explanation is that this body's chest
   is the wrong shape to carry these breasts at all -- which would be a statement about the body,
   and the end of this line.

   **And a rule that outlives this line:** a correspondence's known answer must be scored at the
   separation it will be used at. A surface against itself displaced is a test of arithmetic; two
   different bodies a centimetre apart is the test of the instrument.

   **A normal-only model gives up 6-9% of the deformation**, far inside the 50% that would have
   ended the line, and the normal part -- 6.5 to 10.9 mm -- is the chest-wall offset this file has
   recorded all along. A chest wall differs from another chest wall mostly in how far forward it
   sits, as predicted, and that is now measured rather than assumed. (The pairing here is
   NEAREST POINT, because a decomposition cannot use a correspondence that imposes a direction; it
   carries its own d^2/R bias, which inflates neither fraction by anything like the margin.)

   **Reported, and it matters for the bar:** 20 Taubin iterations move the surface they smooth by
   0.416 mm mean and 5.501 mm max. The gate A bar of 0.5 mm therefore sits BETWEEN the measure's own
   floor (0.288 mm, an identity warp on coinciding surfaces) and the smoothing displacement, so at
   0.9 mm the instrument is roughly twice a bar that is itself at the geometry's noise level. Said
   here, as with the 1 mm bar before it, so nobody quotes the gap without it.

   **Gate A stands FAILED, and the bar was mine and unachievable (2026-09-10).** 0.5 mm was set by
   halving a number that was itself arbitrary. Against the floors now measured -- 0.288 mm for the
   measure and 0.416 mm mean for the smoothing the instrument requires -- a recovery cannot do
   better than about **0.51 mm** on these meshes by construction, so the gate could not have been
   passed by any implementation. That is a defect in the gate, and the result is recorded as a
   failure against it rather than quietly rescored. What the run does establish stands on its own:
   the scalar representation beats the vector warp everywhere (1.5 mm to about 1.0) and, unlike it,
   **improves with scale** -- 0.897 at L = 40 against 1.046 at L = 5 -- which is the signature of a
   representation that matches its deformation. The flat curve is gone.

   **What replaces it, and what it does NOT certify.** An absolute bar at the noise floor asks the
   wrong question; what matters is whether the error left behind is small against the correction
   being made. So: **the recovery residual must be at most 10% of the normal deformation it
   carries** -- 6.5 to 10.9 mm here, so 0.65 to 1.09 mm. That is derived from the deformation, not
   from any score. **But the current numbers were seen before this was written, and 0.9 mm is about
   10% of 9 mm, so this criterion certifies nothing by itself.** It is recorded as the standing
   criterion for future fits, and **the verdict on this line rests on gates 2-4**, which are unseen:
   laterality, the held-out ribs at the whole-torso fit's own 4.55-5.60 mm, and breast volume more
   than 20 mm behind the muscular wall at 1% per subject. If those pass, the chest wall is carried;
   if they do not, a 10% representation error was not the thing standing in the way.

   **Normal agreement is adopted, on the argument and not on the score (2026-09-10).** A hit whose
   surface faces away from the source is on the far wall of the rib: it is not the partner of the
   source point, it is a different part of the bone that happens to lie along the ray, and the
   return test cannot see it because the two cortical walls are parallel. That is a correctness
   condition on what a correspondence IS, and it would have been right to specify it before the
   numbers existed -- the pre-registration simply did not anticipate a surface thin enough to shoot
   through. It is adopted for that reason, with the reason stated, and **its gate-0 numbers are
   reported rather than claimed as a pass**: 0.078 mm mean and 0.157 mm p90 were seen before the
   rule was adopted, so they are evidence about the instrument and not a verdict on it. The agent
   that found this implemented it OFF by default and refused to claim the pass, which is why the
   distinction is available to make at all.

   **What the instrument is now judged on is gate 1, which has not been run with it**, at the same
   1 mm RMS and with both readings reported; then gates 2-4 unchanged. Gate 0 stands as written for
   any future correspondence -- <= 0.2 mm mean, <= 0.5 mm p90 against a known preimage -- and a
   correspondence that needs a filter to reach it declares the filter.

   **Also reported, not judged:** 319 of 4,800 samples are dropped by the agreement test, and where
   they fall is the finding rather than a nuisance -- right ribs 2 and 3 carry means of 1.548 and
   1.396 mm against about 0.12 mm on the other ten, because those are the thinnest sheets in the
   set. **The skin line inherits this**: `ihm/anatomy/normal_shooting.py` is written as a reusable
   instrument with its own sphere known answer (recovering a 2%-inflated copy's radial partner to
   4.9e-4 of the radius), and any line that shoots normals at ribs or scapulae meets the same
   thin-sheet case.

   **Seat it in two steps: place, then conform. Fixed 2026-09-10, before it runs.** 45 mm of
   overlap is not tissue deformation, it is placement -- the breast is another woman's tissue where
   a similarity registration put it, and a quasi-static solve pushing 45 mm of interpenetration out
   is the wrong instrument for that. A surgeon seats an organ and then it conforms; so:
   * **Step 1, rigid.** The rigid motion (translation and rotation, NO scale, so volume cannot
     change) that minimises the sum of squared penetration depths of the base nodes measured along
     their own outward normals, by L-BFGS over the six degrees of freedom from the registered pose.
   * **Step 2, conform.** The sliding boundary condition exactly as judged above, from the placed
     state.
   * **The four per-breast gates are unchanged** -- (a) volume within 1% of the undeformed breast,
     (b) every tet J > 0.2, (c) the two solvers within 5% of maximum displacement, (d) no flipped
     base triangle -- and (a) stays meaningful because a rigid motion conserves volume exactly.
   * **Reported per breast:** the rigid translation and rotation, and the penetration left after
     step 1.
   * **The check that keeps step 1 honest:** a breast needing more than **25 mm** of translation is
     recorded as a REGISTRATION failure for that subject, not a seating -- at that point the
     placement, not the tissue, is what is wrong, and the breast is reported unseated rather than
     moved until it fits. Whether the sliding condition is right for a breast is what step 2 tests;
     step 1 exists so that test is reachable, not to make it pass.

   **Run on s1159 left: unseated, and the objective is why.** Step 1 as written minimises
   penetration ONLY, and that has a degenerate minimum -- carry the breast away and every ray misses
   the muscle. Unbounded, L-BFGS took it in one step (999.7 mm). Bounded to the region where
   placement means anything, the minimiser from the registered pose ran to the CORNER of the box
   (translation [25.0, -25.0, 21.7] mm, magnitude 41.5 mm, rotation 18.9 deg): the objective fell
   15x while the WORST penetration got worse, 45.5 -> 56.6 mm, and 672 of 2,850 base nodes lost the
   muscle on their ray. By the 25 mm rule that is a registration failure and s1159 left is reported
   unseated -- but the verdict is a property of the objective, not evidence about the registration,
   and it is recorded as such. The other seven were not attempted.

   **What a rigid placement can and cannot do**, swept directly with no optimiser: a pure anterior
   translation of the registered breast (mean base normal [0.308, 0.160, 0.938]).

   | anterior shift | penetration max | median | nodes penetrating | nodes with muscle on their ray | objective |
   |---:|---:|---:|---:|---:|---:|
   | 0 mm | 45.5 mm | 18.5 mm | 2597 | 2850 | 1.157 |
   | 10 mm | 58.4 mm | 12.1 mm | 2721 | 3082 | 0.612 |
   | 20 mm | 46.6 mm | 8.5 mm | 2489 | 3291 | 0.460 |
   | 25 mm | 49.4 mm | 8.0 mm | 1838 | 3264 | 0.470 |

   The median halves and MORE nodes come over the muscle -- the opposite of escape -- but the
   deepest nodes stay 45-58 mm inside at every shift. So a rigid placement improves the typical node
   and does not remove the extreme overlap, which is the part neither solver can push out. Step 1
   does not make step 2 reachable for this breast.

   **Two implementation faults to fix before this is re-run, neither applied here.** The 25 mm rule
   is on translation MAGNITUDE; the search bounded each COMPONENT, so the magnitude could reach 43 mm
   and the rotation vector 26 deg -- a magnitude constraint is needed, not a box. And the escape has
   a cheap fix: count a node that loses its bed at its pre-placement penetration, so carrying tissue
   off the muscle cannot pay. Both change a pre-registered step and are argued, not silently made.

   **Uterus and ovaries: the pelvic registration, and its gates, fixed before any fitting.**
   UT-EndoMRI (owner-approved 2026-09-10; an endometriosis cohort, NOT a typical-anatomy
   reference) gives uterus and ovary labels on pelvic T2 MRI. TotalSegmentator's
   `total_mr` labels both hip bones, the sacrum and both femurs on these scans (tested on
   D1-000: 24 classes, the rater's uterus inside the box the bones span). The T2 is 5 mm
   slices over 160 mm, so the femurs and possibly the iliac crests are CUT: bone centroids
   would be dragged toward the scan, as s1067's cut ribs were. So the registration pulls
   the MRI's partial bone surfaces one way onto this body's complete bones. Gates:
   * **known answer:** truncate one of this body's own bones, move it by a known
     similarity, and the one-way fit must recover it (translation within 2 mm, rotation
     within 2 deg, scale within 1%);
   * **laterality:** the MRI's left hip bone maps nearer this body's left hip bone;
   * **containment:** the mapped uterus lies inside this body's pelvic ring (the convex
     hull of the hip bones and sacrum), and overlaps no bone (<= 1% of its volume).
   Overlap with this male body's prostate, seminal vesicles and bladder is EXPECTED -- a
   female pelvis variant does not exist yet -- and is reported, not judged.

   **Registered on the first two subjects; every gate passes** (`scripts/register_pelvic_organs.py`).
   Known answer: this body's own bones, cut as a 160 mm scan would cut them (kept surface
   hips 0.71-0.74, sacrum 0.91, femurs 0.22), moved by a known similarity, are recovered to
   0.033 mm, 0.018 deg and 0.004% scale.

   | subject | scale | one-way residual | laterality | uterus in pelvic ring | uterus in bone | reported only |
   |---|---:|---:|---|---:|---:|---|
   | D1-000 | 1.003 | 4.46 mm | pass | 1.000 | 0.000 | 1.4% inside this body's bladder |
   | D2-000 | 0.929 | 3.36 mm | pass | 1.000 | 0.000 | none |

   **The ovaries, measured rather than reported (2026-09-10).** The gates judged the uterus and
   reported the ovaries as a share inside the ring, which says nothing about whether they land
   where ovaries land. `scripts/measure_registered_ovaries.py` over the completed batch: 66 of 87
   subjects carry an ovary mesh, 75 pieces in all.

   | | median | 10-90% |
   |---|---:|---:|
   | lateral offset from the midline | 35.0 mm (unsigned) | 12 to 41 mm |
   | distance to the nearest hip bone | 24.2 mm | 15.9 to 38.1 |
   | height above the ring floor | 100.2 mm | 85.6 to 118.5 |

   74 of 75 pieces sit at >= 0.99 inside the ring, and the sides are even: 36 on this body's left,
   39 on its right. Known answers first, as always. **(1)** The lateral axis is fixed by this body's
   own hips, 138 mm apart. **(2)** The uterus is a midline organ and the ovaries are not: median
   |lateral offset| 8.2 mm over 87 uteri against 35.0 mm over 75 ovary pieces -- if that had come
   out the other way the frame would have been wrong and nothing else readable. **(3)** A woman has
   one ovary per side, so the two pieces of a two-piece subject should fall on opposite sides:
   **8 of 9 do, and D1-026 does not.**

   **The failure is in my assumption, not in the registration.** D1-026's two pieces sit 1.0 mm
   apart laterally and 10.5 mm vertically, both on the left, 42-45 mm from the hip bone: that is one
   structure segmented into two fragments, not two ovaries on one side -- and the extraction keeps
   any piece within 50 mm of the label's largest, so a fragmented ovary arrives as two pieces by
   design. "Two pieces means one per side" is therefore void as a rule for this archive, and in an
   endometriosis cohort it would also fail honestly on kissing ovaries, where both are drawn
   adherent on one side. The number that survives is the other eight, and the caveat on every
   artefact stands: these ovaries are displaced and enlarged by disease, so this measures where the
   labels land, not where a healthy ovary sits.

   Both subjects in the first pair carry one ovary piece, inside the ring. **One threshold was not
   fixed in advance:**
   the gate said "inside the pelvic ring" without a number; the implementer chose >= 0.99 of
   the uterus surface. Both subjects read 1.000, so it decided nothing here, and it is
   recorded as the implementer's choice. What the gates cannot say: where the uterus sits
   WITHIN the ring, its orientation, or its relation to bladder and rectum -- this body has
   no female soft-tissue ground truth -- and one similarity onto a male bony pelvis absorbs
   real size differences into its scale (D2-000: 0.929).

   **The full extraction, and a claim it overturns.** `scripts/extract_ut_endomri.py` over all
   124 subjects (strict-grid route): **87 with a uterus mesh**, 78 with an ovary label kept.
   Uterus volume 29-766 mL, median 104 -- above a typical adult uterus, as a pathology cohort
   predicts. Ovary pieces per subject: one in 64, two in 12 (surgical absence or one-sided
   labelling; the data cannot say which). Uterus inter-rater Dice over the 12 surviving D1
   pairs: median 0.785.
   The earlier statement that "every label is on the T2 grid" came from two subjects and is
   WRONG for the archive: across D1, 35 of 83 uterus rater labels (24 ovary, 21 endometrioma)
   were drawn on MRI series that are not in the release, and some subjects' raters sit on two
   different released grids. No label is the same shape as an image with a small offset, so
   the tolerance was never the problem and the strict gate's exclusions are correct under
   its rule. Whether those labels can be recovered in scanner world space was tested
   against a criterion fixed before any cross-grid pair was computed -- the median cross-grid
   Dice at least half the median same-grid Dice:

   | organ | same grid | different released grids | one off the release |
   |---|---:|---:|---:|
   | uterus | 0.812 (19 pairs) | 0.798 (6) | 0.708 (17), some near 0 |
   | ovary | 0.602 (20) | 0.626 (5) | **0.162 (5)** |

   Both organs pass AS WRITTEN, and the ovary pass needs saying plainly: the pooled number is
   carried by the released-grid pairs, and off-release ovary labels do not agree. Off-release
   uterus labels mostly agree, with gross failures (D1-019 0.000, D1-027 0.031, D1-004 0.246 --
   most likely the patient moved between series). So the recovery route, added after the
   result and recorded as such: merge raters across RELEASED grids in world space; accept an
   off-release uterus label only where an on-grid uterus label of the same subject confirms it
   at Dice >= 0.406 (the committed criterion, applied per pair); never recover an off-release
   ovary. It was first predicted to recover five uteri; it recovers **four** (D1-001, D1-008,
   D1-009, D1-017), **87 -> 91**. The fifth, D1-034, was my counting error: the prediction took
   a pair where BOTH labels lie off the release, which anchors nothing to a released image, and
   the extractor applies the rule correctly. Built as `extract_ut_endomri.py --route world`
   (`data/derived/ut-endomri-organs-v2-world`). Gates: the strict route through the patched code
   reproduces the committed strict manifest exactly; subjects whose accepted raters are unchanged
   keep their volume to within 5.9e-8 relative (the exact-equality version of that gate FAILED on
   8 subjects, all float noise from two header sources agreeing only to the affine tolerance).
   **D1-004 is ambiguous:** its raters outline 116 mL against 649 and 680 mL; the strict route's
   two-rater intersection gave 89 mL and the world route's majority of three gives 609 -- with two
   raters, one outlier collapses the intersection to the smaller outline.

   **The four recovered uteri, registered.** D1-001, D1-008 and D1-009 pass every gate. **D1-017
   fails containment as written**: 0.688 of its uterus inside this body's pelvic ring. The
   diagnosis, from a check defined AFTER the failure and therefore evidence, not a gate: it is a
   726.5 mL uterus whose lowest 25% (and lowest 10%) lies entirely inside the ring, while the
   part outside bulges 44 mm up and 32 mm forward out of the hull's front face above the pubis,
   below the iliac crests. The fit is ordinary (residual 4.9 mm, laterality correct, 0.37% into
   the sacrum). That is an enlarged uterus seated correctly, and it says the containment gate
   cannot judge an enlarged uterus: it fails correctly placed pathological uteri, which this
   cohort is full of. A lowest-quartile containment test would separate misplaced from enlarged;
   it is recorded as a proposal, to be applied to every subject and labelled post-hoc if
   adopted. The gate is unchanged, and the main batch still holds D1-041 (766 mL) and D1-045
   (618 mL), predicted to fail it the same way.

   **The proposal, applied to every registered uterus (post-hoc, 2026-09-10).**
   `scripts/score_uterus_lowest_quartile.py` imports the registration's own bones, sampler and
   ring. Known answer: it reproduces every manifest's `uterus_in_ring` exactly (18 of 18).
   "Lowest" is along this body's own vertical (calcanei below the frontal bone: axis 1, +).
   The test is lowest 25% of the surface >= 0.99 inside the ring, with bone <= 0.01. Volumes
   below are in the registered frame, extracted x scale^3 (D1-017: 726.5 x 0.941^3 = 605 mL).

   | subject | registered volume | as written | lowest 25% | lowest 10% | gate | proposal |
   |---|---:|---:|---:|---:|---|---|
   | 16 others | 46-219 mL | 1.0000 | 1.0000 | 1.0000 | PASS | PASS |
   | D1-017 | 605.0 mL | 0.6881 | 1.0000 | 1.0000 | FAIL | PASS |
   | **D1-035** | **142.5 mL** | **0.9830** | 1.0000 | 1.0000 | FAIL | PASS |

   D1-017 is the case the proposal was built for. **D1-035 is the reason it is not adopted.** It
   is an ordinary-sized uterus (about 133 mL extracted, near the cohort median of 104). It
   misses the as-written threshold by 1.7% of its surface, which is the implementer's 0.99
   rather than a number fixed in advance, and the proposal passes it too. So the proposal
   separates "seated at the bottom" from "not seated", not enlarged from misplaced. It is more
   permissive than its motivation, and a test like that is how a misregistration gets waved
   through. **The gate as written stays the verdict for every subject.**

   Where D1-035 leaves the ring (`scripts/locate_uterus_outside_ring.py`, evidence after the
   verdict): the uterus sits LOW, 50-92 mm above the ring's floor in a ring 206 mm tall, and the
   1.7% outside is at mid-height, about 35 mm anterior of its centroid and barely above it (+2
   mm), 23-28 mm beyond the nearest ring bone -- out through the convex hull's FRONT face. D1-017's
   bulge instead rises (+39 mm up, +34 mm anterior, 59-205 mm, up to 93 mm from the bone). So
   D1-035 reads as a normally seated uterus whose body leans forward (the usual anteverted
   position) past a hull face that is not an anatomical wall. Neither a volume rule nor a
   lowest-quartile rule describes that; the honest limit is the one already named -- a convex
   hull of three bones is not the pelvic cavity.

   **The batch, complete: 87 registered, 76 pass, 11 fail, in two disjoint modes.** No subject
   fails both, and **no subject fails laterality** -- every one of the 87 maps its left hip to this
   body's left hip. 37 of the archive's 124 subjects have no uterus mesh to register at all. The
   manifest was completed to 87 by re-registering D1-000, D1-002, D1-003 and D1-004, whose records
   (not their meshes) the pre-merge overwrite had destroyed; all four pass.

   | mode | subjects | what fails |
   |---|---|---|
   | containment (6) | D1-035, D1-039, D1-041, D2-061, D2-062, D2-072 | uterus surface outside the ring, 0.766-0.983 against >= 0.99 |
   | bone overlap (5) | D2-026, D2-042, D2-044, D2-048, D2-073 | 1.4-23.4% of the uterus interior inside a bone -- **the sacrum in every one** -- against <= 1%, while the ring share is 0.998-1.000 |

   By sub-dataset: D1 3 of 27, D2 8 of 60 -- alike as rates, but every bone overlap is D2.

   The two modes point in OPPOSITE directions. Every containment failure leans forward: the part
   outside sits 29-45 mm anterior of its own centroid (D1-035 +35, D1-039 +38, D2-061 +29,
   D2-062 +45), out through the hull's front face. The sacrum cases lean back -- D2-044's small
   outside fraction is 25 mm POSTERIOR of its centroid.

   **It is not fit quality, and it is not size.** The five sacrum cases fit like the cohort:
   scale 0.957-1.059, trimmed residual 4.0-4.8 mm, against cohort means of 0.985 +/- 0.063 (D1)
   and 0.962 +/- 0.049 (D2) at ~4.3 mm; their uteri are 52.7-179.4 mL, ordinary FOR THIS COHORT
   (whose IQR is 68.8-154.8 mL) though the top of that span is 1.5x an adult non-gravid
   reference -- see the cohort section below; the argument here is that the five are not
   outliers among their peers, which stands either way. Failure
   rates by sub-dataset are alike -- D1 3 of 27, D2 8 of 60 -- but **every bone overlap is D2**,
   and D1_MHS and D2_TCPW are different sites with different rating protocols (three raters
   against one).

   Scored against the lowest-quartile proposal (`score_uterus_lowest_quartile.py`, known answer:
   the full-surface ring share reproduced exactly for all 91 registered uteri, the 87 here plus
   the 4 world-route ones): of the six containment failures it would pass four (D1-035, D1-041,
   D2-061, D2-072) and still fail two (D1-039 at 0.9892, D2-062 at 0.9262), and it fails all five
   sacrum cases through its bone criterion. It remains unadopted, for the reason D1-035 gave.

   Two candidate causes, and this batch does not separate them:
   * **the target pelvis is MALE** -- a male sacrum is narrower with a more prominent promontory,
     so a uterus mapped correctly in its own pelvis can land inside this body's sacrum. The
     section above already names the male pelvis as what the similarity's scale absorbs;
   * **the cohort is endometriosis**, where cul-de-sac disease retroverts and fixes the uterus
     against the sacrum; the source card records "altered uterine position" as expected here.
   A uterus is never inside a sacrum in a living person, so the overlap is an artefact of mapping
   whichever way it arose. The gate is unchanged and these subjects stay failed.

   **Where in the sacrum, and what it rules out.** This body's sacrum spans 144 mm. The overlapping
   interior points sit at the MIDDLE of it -- 49% of its height for D2-044 and D2-026 (10-90%:
   37-61%), 65% and 64% for D2-048 and D2-042, 41% for D2-073 (39-42%) -- where 0% is the coccyx
   end and 100% the promontory. **All five, 41-65%, and not one at the promontory.**
   (Resampled independently of the registration's own count: 24.1% against its 23.4% for D2-044.)
   A male promontory in the way would put the overlap HIGH, at 90-100%. It is not there. Mid-sacrum,
   S2-S3, is exactly where the rectum and the cul-de-sac lie between uterus and bone, so these uteri
   are mapped 2-4 cm too far back, through a space that is not empty in either sex. That leaves the
   endometriosis reading -- a retroverted uterus adherent to the cul-de-sac, which this cohort is
   selected for -- as the likelier of the two, and it is still a reading, not a measurement: a
   posterior bias in the fit would look the same from here. What would separate them is the rectum,
   which this body has and which nothing has yet been asked about.

   **Prediction, fixed now, before either is registered:** D1-041 (766 mL) and D1-045 (618 mL)
   FAIL containment as written and PASS the lowest-quartile test, as D1-017 does. If either
   fails the lowest-quartile test, it is misplaced, not merely enlarged, and it is reported as
   such.

   **The prediction was half wrong** (both registered by the batch the same day; 37 uteri scored,
   the known answer still reproducing every manifest's ring share exactly):

   | subject | registered volume | as written | lowest 25% | predicted | got |
   |---|---:|---:|---:|---|---|
   | D1-041 | 631.3 mL | 0.8753 FAIL | 1.0000 PASS | FAIL / PASS | **as predicted** |
   | D1-045 | 461.4 mL | **1.0000 PASS** | 1.0000 PASS | FAIL / PASS | **wrong: it passes** |
   | D1-039 | 319.9 mL | 0.9604 FAIL | **0.9892 FAIL** | not predicted | fails both |

   D1-045 is 461 mL registered and sits entirely inside the ring, so **volume does not predict a
   containment failure**: what matters is where the uterus points, not how big it is. That is the
   same lesson D1-035 gave from the other end (142 mL, fails as written), and it further weakens
   the case for any volume-based amendment.

   D1-039 is the first uterus to fail BOTH, and by the rule stated above that reads as misplaced
   rather than enlarged -- but it fails the lowest-quartile test at 0.9892 against 0.99, which is
   0.1% of its surface and no distance at all. Its outside part is 38 mm anterior of its centroid
   and 31-59 mm past the nearest ring bone, at mid-height: the same forward lean as D1-035, larger.
   Calling that "misplaced" would be reading a threshold to three decimals it was never worth to;
   it is recorded as marginal and nothing is inferred from it. The sternum
   (gate d1) is untouched by any breast correction; that is the skin-envelope problem.
2. ~~**A path refitting tool.**~~ **Done.** `libosimActuators.so` already
   exported 97 `PolynomialPathFitter` symbols; what was missing was 300 lines of
   driver. A refit of the unchanged model is 16.7× closer to the model's own
   `GeometryPath`s than the file the engine currently runs.
3. **An anthropometric correspondence.** Unchanged as a *measurement* gap, and
   now explicit as a *declaration*: `pelvis_breadth_ratio` carries, in its own
   schema entry, that NHANES measures hip circumference at the buttocks while the
   parameter is a skeletal factor. Bi-iliac breadth, pelvic inlet shape and
   biacromial breadth are measured by **no** catalogued source, in either sex —
   which is why `shoulder_breadth_ratio` is 1.0 in every sex preset. ANSUR II has
   biacromial breadth by sex and is not catalogued here.
4. **A gonadal axis attached to a gonad.** The absence is now measured rather
   than described (11 supporting entities, all male; 0 uterine; 0 placental), and
   the audit turned up a second problem in front of the first: the testis meshes
   are 6.7 mL against a 12–30 mL reference.

   **Audited across ten solid organs (2026-09-10, `scripts/audit_organ_volumes.py`), and the
   headline is not the testis.** **Seven of the ten are not watertight**, so no volume can be
   computed for them at all: left kidney 10 boundary edges, right kidney 22, spleen 58, pancreas 52,
   gallbladder 6, urinary bladder 24, stomach 22. A divergence-theorem volume on an open surface
   returns a number rather than an error, which is exactly how this would have been missed.

   | organ | volume | reference (stated, not catalogued) | ratio to low end |
   |---|---:|---:|---:|
   | left testis | 6.72 mL | 12–30 mL | **0.560x** |
   | right testis | 6.78 mL | 12–30 mL | **0.565x** |
   | prostate | 11.59 mL | 15–30 mL | **0.772x** |

   **The testis is not alone** -- every organ that can be measured reads low -- but three organs
   cannot distinguish a scale fault from three mesh faults, and the spread (0.560 to 0.772) is wide
   enough that they are not obviously one factor. The reference ranges are typical adult male values
   stated in the script so they can be argued with; they are **not** from any source catalogued in
   `data/sources/`, and no gate is attached to them.

   **A tempting link that does NOT hold, recorded so it is not made later.** The mass ledger records
   tissue volumes at sourced densities weighing 60.00 kg against a profile demanding 70.77, absorbed
   by a uniform 1.179x inflation. Undersized organs look like an explanation and are not: the ten
   audited organs total roughly 1.1 L at the low end of their references, so even at 0.56x they
   account for **under half a kilogram of the 10.77 kg shortfall, under 5%**. The organ defect and
   the mass shortfall are separate problems.

   **The known answer failed first, and the bar was not moved.** The gate initially asked an
   inscribed icosphere to match its *sphere's* closed form within 2e-3 and read 2.17e-3 at
   subdivision 4. That is not a fault in the volume reader -- an inscribed polyhedron understates by
   construction -- so the bar was testing tessellation fineness while claiming to test the reader.
   The bar stands as failed and the **instrument** was replaced, which is the permitted move: the
   cube now gates exactness at machine precision (rel 2.12e-16), and the sphere's deficit must
   converge at **second order**, falling by a factor near 4 per subdivision. It measures **3.981x
   and 3.995x**, with the deficit one-signed as an inscribed polyhedron requires. That is a sharper
   test than the absolute bar could ever have been: a sign error, a factor of two, or a winding bug
   does not produce clean second-order convergence toward the analytic value.

   **The census across all 4,000 meshes, and the correction it forced on my own framing**
   (`scripts/audit_mesh_watertightness.py`). **1,311 of 4,000 (32.8%) are closed**; 2,651 are open
   and 37 non-manifold. By class: lymphatic 98.2% closed, integumentary 81.4%, connective 77.6%,
   skeletal **31.8%**, muscular 22.6%, arterial 7.8%, venous 6.7%, respiratory 2.9%. The largest
   holes are sheet-like structures where openness is the shape rather than a defect -- external
   intercostal 62,748 boundary edges, mesentery of small intestine 54,232.

   **I wrote the census up as showing that mass properties were reading meaningless volumes. That
   was wrong and is withdrawn.** `mechanics.json` already carries a `volume_basis` on every entity
   and already distinguishes these cases -- 1,068 are labelled "unclosed-surface signed integral
   volume prior; open boundaries and orientation are unresolved", 1,425 use "surface area times
   assumed wall thickness", and so on. The ledger knows. A census alone was never going to find an
   unhandled defect, because the handling is declared per entity.

   **What the census IS good for is a cross-check, and that found something specific.** Two
   independent determinations of the same property -- the ledger's declared basis against edge
   incidence measured from the mesh:

   | | count |
   |---|---:|
   | agree, both say closed | 1,197 |
   | agree, both say open | 1,163 |
   | basis makes no closedness claim | 1,577 |
   | **DISAGREE: ledger says closed, the edges say open** | **62** |

   **All 62 sit in one basis class**, and it is the one whose claim is load-bearing: *"a BodyParts3D
   viscus is one closed surface enclosing wall AND lumen, so the enclosed volume carries gas or
   chyme at tissue density. Tissue volume is surface area times an assumed wall thickness; the
   remainder is lumen."* They are the hollow viscera -- **trachea, oesophagus, stomach, duodenum,
   ileum, ascending / transverse / descending colon, rectum**. Of the 65 entities carrying this
   basis, **62 are not closed surfaces**.

   **Why it matters, and how much.** The wall-versus-lumen split is derived from an enclosed volume,
   and an open surface does not enclose one -- the signed integral returns a number governed by
   where the holes are. The basis rests on a premise its meshes do not satisfy. The scale is bounded
   and should not be overstated: this class carries **1.581 L** of allocated tissue volume, roughly
   **2%** of the 70.77 kg target. A real defect in a declared basis, not a crisis in the ledger --
   and exactly the kind that survives because the sentence asserting it reads like a measurement.

   **This also explains a row in the organ audit above**: the stomach read UNMEASURABLE there at 22
   boundary edges while the ledger treats it as a closed viscus. The same fact from two directions.

   **A known answer failed here too, and the expectation was the wrong side.** The detector was asked
   for "2 open, 1 non-manifold" from a cube with one duplicated face. Duplicating a triangle raises
   each of its three edges from two uses to three, so the right answer is **3 non-manifold and 0
   open** -- which is what the detector printed. The expectation was corrected, the detector was
   not, and the census had not yet run, so nothing was tuned to a result. A failed known answer is
   only as useful as the care taken to find out which side of it was wrong.

   **Measured, and my own finding mostly dissolves** (`scripts/measure_viscus_hole_error.py`). The
   size of an error decides whether it is worth repairing, so the 62 viscera were capped -- boundary
   edges chained into loops, each loop fanned to its centroid -- and the capped volume compared
   against the raw one. Known answer first: a unit cube with one and with two faces deleted must cap
   back to **exactly 1.0**, which tests boundary detection, loop chaining, fan orientation and the
   volume reader together, and catches a backwards fan because that gives 1 − 2×(cap) rather than 1.
   Both pass at machine precision.

   **The enclosed-volume shortfall across the 62 is median 0.0%, min −11.9%, max +0.5%.** The holes
   are real and they are tiny -- one to five loops of a few edges on meshes of tens of thousands of
   triangles. **So the premise violation recorded above has essentially no numerical consequence for
   the ledger.** The basis sentence is still wrong -- those surfaces do not enclose anything, and the
   claim should be corrected -- but nothing downstream of it needs re-deriving, and the 1.581 L is
   not in doubt. A defect in a stated premise is not automatically a defect in a number, and I had
   written the first as though it implied the second.

   **Two residues that do not dissolve:**

   1. **`mucosa of stomach` reads −412.46 mL, capped −368.58, a −11.9% shortfall.** The sign is the
      real finding: a **negative** signed volume means inside-out winding, which is a different
      defect from a hole and one the capping cannot fix.
   2. **`jejunum` at −4.5%** is the only other viscus above half a percent.

   **A census of winding over every closed mesh, since the sign is checkable everywhere.** Of the
   **1,311** closed surfaces, exactly **4** have negative signed volume: **left and right acetabular
   labrum (−1.91 mL each) and left and right glenoid labrum (−0.76, −0.75 mL)**, 5.3 mL in total.

   **The ledger is not affected and already says why**: their basis is "absolute signed surface
   integral", and the stored `volume_m3` for each is the positive magnitude (1.9083, 1.9058, 0.7577,
   0.7541 mL). Taking the absolute value is the right thing for a mass and it makes the sign
   invisible to the ledger -- which is precisely why nobody had looked at it.

   **Where it is NOT invisible is anything that uses the surface's normal direction** -- contact,
   signed distance, inside/outside tests, normal shooting. An inside-out labrum has its normals
   pointing into the joint. All four are fibrocartilage joint structures promoted as candidate
   tissue elements in `scripts/build_joint_promotion_candidate.py`. **Whether any of them currently
   carries a normal-dependent element is not established here and is not claimed** -- it is flagged
   as the question worth asking, and it is cheap to answer for whoever next touches the joint
   element set.
5. **Seven measured female `shaft_bearing_fraction` values.** Sharper than
   "measured female values for the hair fields": the parameter exists on every
   field already, six of the seven androgen-dependent ones assume 1.0, and 1.0 is
   an adult male assumption that is not labelled as one.

An honest summary of the current state, for anything that quotes it: *male
anatomy including a complete genital tract and no breast of either sex; one
implemented isotropic size parameter and three implemented anisotropic
proportional parameters, gated against a measured refit noise floor and accepted
by the plant; a path refitter that beats the shipped path set by 16.7× on the
model's own muscles; a `sex` parameter that moves four measured proportions and
zero anatomical entities; and no female-specific geometry anywhere.*

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
3. Real receptor counts appear not to follow body size — tactile acuity is
   reported to be better on smaller fingers because Merkel-cell density is higher
   on them, innervation being laid down as a roughly fixed count and spread over
   whatever surface grows. **A published prior, not measured here and not
   catalogued in this repository**, and the weaker of the two substantive reasons.

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


   **The reader works and both subject gates fail. The premise that an SDF has no gaps holds for the
   FIELD and not for the ANATOMY.**

   The known answer passes: an ellipsoid with 90/60/40 mm semi-axes, displaced by a known scalar field
   along its own normals at 10 mm separation plus 2 mm of structure, is recovered at 100% coverage
   with 0.00% error. (Stated plainly: that error is near-exact partly by construction, because the
   target was built by displacing the same vertices the rays start from, so each crossing passes
   through a vertex. It tests the reader's mechanics, not its behaviour between vertices.)

   | subject | reader, mean | reference, mean | reader, median | reference, median | coverage |
   |---|---:|---:|---:|---:|---:|
   | s0790 | 6.20 mm | 6.55 mm | 4.00 | 5.38 | 40.5% |
   | s1067 | 6.63 | 11.02 | 5.19 | 10.12 | 41.9% |
   | s1159 | 5.55 | 10.15 | 3.25 | 6.12 | 33.3% |
   | s0970 | 5.98 | 9.31 | 3.38 | 7.49 | 30.8% |

   **Gate 0'' fails** (38-68% away as written, and 5-45% like-for-like), and the gate as written
   compared unlike statistics -- the reference 6.5-10.9 mm was a MEAN over nearest-point pairs, the
   reader reports a MEDIAN. Compared mean to mean, s0790 agrees within 5% and the other three are
   36-45% apart, so the 6.5-10.9 figure is corroborated on one subject and not on three.

   **Gate 0''' fails: coverage 30.8-41.9%.** An SDF cannot miss a thin surface, which is what it was
   brought in for, but a ray along this body's wall normal still finds NO BONE within 25 mm about
   two thirds of the time -- her ribs are separate bones with air between them, and the field is only
   crossed where bone lies along that particular ray. So the gap problem was never about the
   representation of the surface; it is the anatomy, and it has now defeated a pairing, an envelope
   and a field in turn.

   **And the selection survives too.** The rays that do find bone are the ones where the two walls
   are aligned, which is why the reader's median (3.3-5.2 mm) sits so far below its own p90
   (13.8-16.5 mm). Every instrument this line has built reports the aligned minority unless
   something forces it to account for the rest.

   ### The registration-correction branch ends here, and what it leaves behind (2026-09-10)

   Six instruments, in order: prescribed closest-point seating, the sliding base, rigid placement,
   the label trim, a chest-local similarity, a vector warp, a scalar field on smoothed normals, an
   envelope, and a distance-field reader. **Every one failed at the same place** -- establishing
   what this body's chest wall is, relative to hers, densely enough to correct it -- and the last
   three failed for the same reason in three different representations. A rib cage is struts with
   air between them. Pair against it and the rays miss; close it and the closing eats the
   difference; read a field along it and two thirds of the rays cross no bone at all. That is not
   an instrument fault and no further instrument of this kind is pre-registered.

   **What is now known, and it is not nothing:**
   * the deformation is **91-94% normal** (`de88ef4`), so the representation question is settled;
   * the correspondence floor on these meshes is **0.078 mm**, set by faceting as 0.135 x edge^2;
   * a scalar field on the surface beats a vector warp everywhere and improves with scale, where the
     vector warp's curve was flat -- the first instrument in the line that behaved like a matched one;
   * and three rules that outlive it: score a correspondence at the separation it will be used at,
     break the symmetry a known answer is testing, and a tolerance calibrated at zero separation is
     a different test at eight millimetres.

   **A figure this file must stop quoting without its caveat.** The "~10 mm chest-wall offset"
   (line 971, and the 8.0 / 8.3 mm at rib 5 above) is a **mean over nearest-point pairs**, which is
   the estimator this line measured to be biased by d^2/R and to select for already-aligned regions.
   The distance-field reader, which selects differently, agrees with it on s0790 within 5% and is
   36-45% below it on the other three. **It is corroborated on one subject of four and should be
   written as "6-11 mm by a biased estimator, 5.5-6.6 mm by another" wherever it is used**, until
   something measures it without selecting. Nothing downstream of it is safe to quote more precisely
   than that.

   **Where the breast work goes instead, pre-registered before it is built.** The target was always
   to put a female breast on this body, and correcting the registration was one route to it, chosen
   because seating failed first. The other route uses only surfaces that ARE continuous: this body's
   own muscular chest wall, which has no gaps, and her breast's own shape. **Seat the breast on this
   body's wall and let the deformation absorb the registration error**, with the result labelled for
   what it is -- her breast adapted to this chest, not a claim about her anatomy.
   * the sliding boundary condition of `5699cb2`, whose solver gate passed in two independent codes;
   * from the registered pose, on the muscular bed of `88178df`, with the median penetration of
     7.4 mm that the corrected bed leaves -- not the 45 mm that stopped the first attempt, which was
     13 nodes;
   * gates (a)-(d) unchanged: volume within 1%, every tet J > 0.2, the two solvers within 5%, no
     flipped base triangle;
   * **and one addition, because the registration error is now being absorbed rather than fixed:**
     report the deformation's own magnitude per breast, since a seating that moves tissue by more
     than the breast's own dimension is no longer that breast.

   **Predicted:** it seats, because 7.4 mm of median penetration on a continuous bed is what the
   solver was gated for; and the deformation magnitude lands in single-digit millimetres, which is
   a breast adapted to a chest rather than a different breast.

   **The return to seating: it does not seat, and the sliding condition is not what fails.**
   On the continuous muscular bed, from the registered pose, with median penetration 7.4 mm: the
   in-repo solver accepts about 0.1% of the seating per step and inverts elements whenever the step
   grows, and FEBio fails inside its first time step with 42 retries and negative Jacobians. Neither
   is the culprit -- their shared solver gate passed in two independent codes, and the boundary
   condition is the one whose block cases agree to 1.58% of maximum displacement.

   **What fails is the field the seating is asked to apply: it is not a deformation.** The depth each
   base node must travel varies between NEIGHBOURS faster than the tissue can follow.

   | between neighbouring base nodes (11,963 edges, median length 3.02 mm) | |
   |---|---:|
   | change in required depth: median / p90 / max | 0.72 / 3.32 / 34.44 mm |
   | that change as a gradient (per unit edge length): median / p90 / max | 0.24 / 1.08 / 74.7 |
   | edges whose ends must slide past each other by MORE than their own length | **1,311 (11.0%)** |
   | by more than twice their length | 603 (5.0%) |

   A gradient above 1 means the two ends of an edge pass through one another. No boundary condition
   and no solver preserves an element through that, so 11% of this base cannot be seated by any
   instrument of this kind.

   **And the obstruction is roughness, not size.** The motion required is small against the breast
   itself -- 143 x 187 x 116 mm, 263 mm across, 562 mL -- at a median 7.4 mm (2.8% of its diagonal)
   and a maximum 41.8 mm (15.9%). It is not that the breast must move far; it is that neighbouring
   points of its base must move by very different amounts over 3 mm. That is the same roughness that
   collapsed 69% of the base triangles under the first, closest-point rule, seen now through a
   boundary condition built to avoid exactly that.

   **So the seating line and the registration line meet at one statement.** The registered breast's
   base and this body's chest wall differ by an amount that varies faster than the tissue's own mesh
   resolution. Removing that requires either a locally accurate registration -- which the four
   instruments above established cannot be had from a rib cage of struts and air -- or smoothing the
   depth field before applying it, which is a modelling choice that has to be declared: the breast
   would then follow the chest wall rather than the registration, and what is delivered is a breast
   shaped to this body, not this subject's breast placed on it.

   ### Smooth the field until it is a deformation, and declare what that costs. Fixed 2026-09-10

   Option 1 is closed by four instruments. Option 2 is taken, with the smoothing bound by the
   requirement rather than by taste, and with the thing it removes reported as the headline rather
   than a footnote -- because smoothing is exactly what the envelope did when it ate the offset it
   existed to carry, and the only difference here is that this time it is declared.

   * **The rule.** Smooth the depth field ON the base surface (never the geometry), increasing the
     bandwidth monotonically from zero, and **stop at the first bandwidth where the field is a
     deformation**: max |gradient| between neighbouring base nodes **<= 0.5**, half the value at
     which an edge's ends pass through each other. The stopping rule is the invertibility the
     solver needs, so nothing is chosen to make a result come out.
   * **Reported as the headline:** the bandwidth required, the fraction of the field's magnitude it
     removes, and the median and p90 depth before and after. A reader must be able to see how much
     of this subject survived.
   * **The stopping condition, fixed now:** if reaching a gradient of 0.5 removes **more than half**
     the field's magnitude, the breast is being reshaped by the chest rather than seated on it,
     and that is the result -- reported, not smoothed further.
   * **Gates (a)-(d) unchanged** -- volume within 1%, every tet J > 0.2, the two solvers within 5%,
     no flipped base triangle -- plus the seating actually completing, and rib points inside the
     breast reported as before.
   * **The label travels with the artefact:** what this produces is *a female breast adapted to
     this body's chest wall, derived from subject s1159 and not a model of her*. Any manifest,
     figure or card carrying it says so, the way the endometriosis caveat travels with the pelvic
     organs.

   **Predicted:** a bandwidth of 10-20 mm suffices, it removes 20-40% of the magnitude, the seat
   completes, and the deformation stays in single-digit millimetres. If it removes more than half,
   the honest reading is that this chest cannot carry this breast without becoming its author, and
   the female-torso work needs a female body rather than a male one corrected -- which is the
   conclusion `data/sources/visible-human-lower-extremity.json` was catalogued against.

   **The smoothing costs almost nothing and the seating still does not complete.**

   | bandwidth (s1159 left, 4,097 base nodes) | worst gradient | edges over 1.0 | magnitude removed | depth median | depth p90 |
   |---:|---:|---:|---:|---:|---:|
   | 0 mm (raw) | **74.71** | 1,311 | 0.0% | 7.38 mm | 16.06 mm |
   | 5 mm | 1.47 | 71 | -0.4% | 7.14 | -- |
   | 10 mm | 0.72 | 0 | -0.4% | 6.74 | -- |
   | 12 mm | 0.56 | 0 | -0.2% | 6.81 | -- |
   | **20 mm (first at the bar)** | **0.39** | **0** | **0.3%** | **7.48** | **11.87** |

   **The bandwidth required is 20 mm and it removes 0.3% of the field's magnitude** -- not the 20-40%
   predicted, and nowhere near the half that would have stopped the run. The roughness that made the
   field impossible to apply carries almost none of its size: smoothing takes the worst neighbour
   gradient from 74.71 to 0.39 and the edges that had to fold through themselves from 1,311 to zero,
   while the median depth moves 7.38 -> 7.48 mm and the p90 16.06 -> 11.87 mm. (The worst gradient is
   not monotonic in bandwidth -- 0.56 at 12 mm, 0.74 at 15, 0.39 at 20 -- so "the first bandwidth at
   the bar" is 20 mm by the rule as written.)

   **And it still does not seat.** With the smoothed field the solve reaches a fraction of **0.0156**
   -- 1.6% of the seating -- and then stops: the same **13 elements** invert at every step size tried,
   down to 9e-5 of the travel, so the increment is no longer what limits it. Those elements are not
   slivers: their reference quality is a median 0.478 against the mesh's own 0.487. Gates (a)-(d)
   were therefore not reached and no artefact was written.

   Two convergence tolerances were declared at the bed's facet scale (1 mm) on the way, because
   0.05 mm asks the solver to resolve the muscular wall finer than the wall is defined: its normals
   turn several degrees over 2-4 mm, so a sliding node's hit point moves ~0.5 mm per pass and its
   measured gap wanders ~0.2 mm. They are convergence criteria, not gates.

   **What this does and does not establish.** It does not establish the fallback reading -- this chest
   is NOT reshaping the breast, because the smoothing that makes the field appliable costs 0.3% of it.
   The obstruction is no longer the registration's roughness, which is now measured and removed at
   negligible cost; it is that a base of 3,123 held nodes driven along their own normals into a
   faceted anatomical bed drives a handful of ordinary elements to degeneracy, and the solver's own
   J > 0.2 floor then blocks every further step. That is the fifth suspect this line has isolated by
   measurement, and the first that lies inside the solver rather than in the data.

   **My fallback is withdrawn, not merely unsupported.** I had held in reserve the reading that this
   line fails because a male chest cannot host a female breast -- the same reading that would say the
   female work needs a female body before it can proceed. The smoothing result kills it: the field
   that could not be applied loses **0.3%** of its magnitude on becoming appliable, so the geometry
   the registration asked for was very nearly the geometry it gets. A host that was reshaping its
   guest would not survive that test. Anything quoting the male-chest explanation for this line
   should stop.

   **Pre-registered before the solver is touched, because the fix and the diagnosis are the same
   people's work.** Two stepping changes are proposed -- remeshing the base layer, and relaxing the
   interior between increments instead of clipping constrained nodes first. Both are plausible and
   neither is evidence. The question they assume the answer to is whether the solver can complete a
   motion that is KNOWN to be completable, and that has a known answer available before any fix:

   * **Control R, the rigid control.** Drive the same 3,123 held nodes through the same bed, same
     solver, same J > 0.2 floor, with a **rigid translation** of the whole breast base equal to the
     smoothed field's median displacement. A rigid motion preserves every element's Jacobian
     **exactly**, so no element can invert under it for any reason of physics or mesh quality. It
     breaks the symmetry the sphere could not: it separates "this configuration is infeasible" from
     "the stepping is wrong", which the current failure cannot distinguish.
   * **Gate R, fixed now:** Control R completes to fraction **1.0 with zero inversions**.
   * **If R passes:** the solver handles feasible motions through this bed, the smoothed field is
     genuinely infeasible for this mesh, and a remesh is attacking the right thing. The 13 elements
     are then a property of the mesh under that specific field and are reported per element.
   * **If R fails:** the fault is inside the stepping and NOT in the field, the mesh, or the breast,
     and every one of this line's five suspects was diagnosed by a solver that cannot complete a
     translation. That outcome invalidates more than it fixes, so it is the one to want.
   * **Predicted: R passes.** The 13 elements sit at reference quality 0.478 against a mesh median of
     0.487, which is ordinary, and ordinary elements do not invert under a rigid motion. The
     prediction is recorded so that a pass is not read as confirmation of anything beyond the gate.

   No remesh and no relaxation change is made until R has a number. A stepping fix tried first would
   be judged by whether the seat completes, and a seat that completes is exactly what everyone
   involved wants to see -- which is the condition under which this programme has been wrong 30
   times.

   **One honest weakness in the bandwidth rule, recorded and NOT acted on.** The worst gradient is
   not monotonic in bandwidth (0.56 at 12 mm, 0.74 at 15, 0.39 at 20), so "the first bandwidth at the
   bar" depends on the sweep's spacing: a finer sweep between 12 and 20 mm could return a smaller
   bandwidth and a different field. The rule was written assuming monotonicity and that assumption is
   false. It is left as written because re-running the sweep after seeing the result would be
   choosing a bandwidth with the answer in view, and because at 0.3% removed the choice cannot matter
   to any conclusion here. If a later result turns on the bandwidth, this paragraph is the reason to
   re-derive it properly first.

   **CONTROL R FAILS, and it invalidates this line's solver-side diagnoses.** Driven by a rigid
   translation of 7.48 mm -- the smoothed field's median, along the held nodes' mean normal -- the
   same 3,123 held nodes through the same bed, solver and J > 0.2 floor stall at fraction **0.0312**,
   with 27 elements inverting at every step size tried down to 9e-5. A rigid motion preserves every
   element's Jacobian exactly and satisfies every held constraint at zero strain energy, so nothing
   can invert under it for any reason of physics or mesh quality. **The fault is inside the stepping.**

   The mechanism is visible in the failure message, which is the same one every seating attempt gave:
   "bounded start inverts N elements". Each increment begins by CLIPPING the constrained degrees of
   freedom to their new bounds while every other coordinate stays at the previous solution. For a
   motion whose feasible path also moves the tangential and interior coordinates -- which a rigid
   translation certainly does -- that start is infeasible by construction, however small the step.

   **What this withdraws.** Every statement of the form "the seating stalls because the field is
   infeasible" is now unsupported, including the reading that a handful of ordinary elements
   obstruct it: the solver stalls the same way on a motion that is feasible by construction. Those
   thirteen elements are not established as a property of the mesh under that field.

   **What survives, because it was measured without the solver.** The edge-gradient census (11.0% of
   base edges requiring their ends to slide past each other by more than the edge's own length,
   worst 74.7) is a property of the field and the mesh, computed directly. So is the smoothing sweep
   (74.71 -> 0.39 at 20 mm bandwidth, 0.3% of magnitude removed), gate B's 91-94% normal
   decomposition, and the whole registration line's coverage, envelope and distance-field results,
   none of which ran a solver. The block-case solver gate also stands as what it always was: two
   codes agreeing on prescribed displacements, where the clipped start is feasible.

   So the seating question is reopened rather than answered: the sliding condition has still not
   been tested on this breast by a stepping that can carry a translation.

   **Blast radius, checked rather than assumed.** Three files import `sliding_contact`:
   `seat_breast_sliding.py`, `verify_contact_broadphase.py` and `assess_soft_body_coupling.py`. Only
   the first calls `solve_sliding`; the other two import the module and never reach it, and no
   document cites either. **The contamination is confined to the breast seating line.** The skin
   warp, the knee cartilage registration and the bed work use different solvers and are untouched by
   this. That was worth ten seconds of grep before anyone worried about it.

   **The diagnostic named the wrong cause, and that is the transferable lesson.**
   `sliding_contact.py:73` raises `"bounded start inverts N elements; reduce the load step"`. The
   second clause is false -- the load step is not what limits it, and shrinking the increment shrinks
   the bound change without touching the inconsistency. That message was read about fifteen times
   across this line and it steered every reading toward the data, because **it asserted a cause with
   the authority of a measurement**. A diagnostic that names a wrong cause is worse than no
   diagnostic, since a silent failure gets investigated and a confident one gets believed. The
   message should state what was observed -- the clipped start inverts -- and stop there. Fixing the
   message is part of fixing the stepping.

   **Why the existing validation could not catch it.** The block case compared two codes on
   *prescribed displacements*, where the clipped start happens to be feasible, and passed at 1.58%.
   It never exercised the sliding path at all. That is `docs/LOG.md` row 30 word for word -- *a
   control only tests the path it exercises* -- and the second time today that sentence has been the
   answer.

   **Why my prediction was wrong, precisely.** I argued that ordinary elements, at reference quality
   0.478 against a mesh median of 0.487, do not invert under a rigid motion. True of the *motion*,
   irrelevant to the *algorithm*: the rigid motion is never applied, because the solver builds a
   non-rigid intermediate state and inverts there. I reasoned about the physics of the target while
   the failure lived in an intermediate the physics never sees. Worth keeping, because designing a
   known answer around the target's properties is exactly how a control gets built that the code can
   still fail for reasons the control does not name.

   **Control R is the standing gate on the stepping, fixed before any fix is written.** A feasible
   warm start is the proposed repair -- extrapolate the previous increment, or apply the bound change
   inside the Newton loop instead of clipping before it. Whichever is built, **R must pass at
   fraction 1.0 with zero inversions before any breast is judged again**, and R is re-run on every
   later change to the stepping. No seating result from this line is readable until it does.

   **The stepping fix is authorised, with a second gate and three named outcomes (2026-09-10).**
   Either repair may be built; the second is preferred as removing the bug rather than routing
   around it, since the active set at line 80 exists to handle exactly this. The choice is the
   implementer's and must be stated.

   * **Gate R** (unchanged): the rigid translation completes to fraction 1.0 with zero inversions.
   * **Gate R2, new and equally binding:** the block case must still pass. It agrees between two
     codes at 1.58% on prescribed displacements, and it passes *because* the clipped start happens
     to be feasible there. Changing the start can break the one path that already worked, and a fix
     that trades a broken case for a working one is not a fix. Both numbers reported, before and
     after.
   * The message is fixed with the code: `"bounded start inverts N elements; reduce the load step"`
     asserts a false cause and cost this line five suspects' worth of investigation. It states what
     was observed and stops.

   **Three outcomes, named before the work so none can be read as a surprise:**

   1. **R still fails** -- the repair did not work. Iterate; the breast is not touched.
   2. **R passes and the seat completes** -- judge by gates (a)-(d) as written, nothing loosened.
   3. **R passes and the seat still stalls** -- the important one. Only then is "the field is
      infeasible for this mesh" an *earned* reading rather than an artefact of the start, and only
      then do the thirteen elements become interpretable as a property of the mesh. That would
      restore, on sound footing, the conclusion withdrawn above, and would be the first thing this
      line has been able to say about the breast itself.

   **No prediction is offered on which outcome obtains**, and the omission is deliberate. Today's
   record on the neighbouring skin line is five refused predictions against one directional hit,
   while fixing criteria before the numbers exist has worked without exception. Guessing the number
   adds nothing here and would only give a result something to be read against.

   **One rule carried over from the skin line, applying the moment a seat completes.** There the fit
   residual proved flat while the truth error doubled, then fell while the error rose -- not a
   degraded measure of accuracy but not a measure of it at all -- and a correspondence with 4.4x
   better targets produced a 26% *worse* warp, three quarters of that 4.4x being the filter
   declining to answer for the hard points. **So agreement with the prescribed displacements is not
   to be reported as quality.** Report a quantity with a known truth behind it; where none is
   available, say so rather than substituting a residual.

   **Control R's premise is wrong, and I am recording that BEFORE its outcome is known.** At the
   time of writing R is mid-run at fraction 0.6250 -- already twenty times past the 0.0312 where the
   clipped start stalled, so the stepping fix plainly works -- with `min J` falling steadily,
   0.596 -> 0.407 -> 0.307 -> 0.275 -> **0.226** against a 0.2 floor. That fall is the tell.

   I justified the gate with: *"A rigid motion preserves every element's Jacobian exactly and
   satisfies every held constraint at zero strain energy. Nothing can invert under it."* **That is
   true only if every node moves rigidly.** Here only the **3,123 held nodes** are driven; the
   interior and the bed-contact nodes are free, and the bed is enforcing unilateral contact
   throughout (0.0556, 0.0725, 0.0236 mm of penetration in the log). A body whose held surface
   translates into a fixed bed **cannot** translate rigidly as a whole, so the rigid translation is
   not a zero-energy admissible solution and elements can legitimately deform -- which is exactly
   what `min J` shows.

   **So "R completes at fraction 1.0 with zero inversions" is not achievable by construction, and
   it is withdrawn as a pass/fail gate.** The flaw is in the premise, not the threshold: I asserted
   feasibility that the bed contradicts. What R remains, and remains good for, is a **comparative**
   diagnostic on a fixed target: the same drive, the same bed, the same floor, before and after the
   stepping change -- **0.0312 against whatever it reaches.** That comparison is valid precisely
   because both arms face the same possibly-infeasible target.

   **Consequences fixed now, before the number:**

   * A stall at, say, fraction 0.8 must **not** be read as the stepping still being broken. Against
     0.0312 it would be a twenty-five-fold improvement on the identical target.
   * If R *does* reach 1.0, that is a **stronger** result than the flawed gate ever required, and it
     will be reported as such rather than as "the gate passed".
   * **Gate R2 is untouched** and remains a real pass/fail: the block case must still agree at
     1.58%, because there the prescribed displacements are feasible and the premise holds.
   * Outcome 3 of the three named above -- "R passes and the seat still stalls" -- loses its
     sharpness with R weakened, so the earned reading that "the field is infeasible for this mesh"
     now needs a control whose feasibility is guaranteed. **The obvious one: drive the held nodes
     with the bed's constraint removed**, where a rigid translation genuinely is zero-energy and
     nothing may invert. That is the control R should have been, and it is pre-registered here as
     R′ before either result is read.

   ### Outcome 1: the repair works, R2 passes, and R′ still fails — on a solution it has already found

   **Gate R2 passes and the converged block answers are identical: the repair changed the path, not
   the result.** Platens against the exact homogeneous field go 0.00004% → **0.00000%** against a 1%
   bar; the cylinder holds at 6.382 µm. **R improved 24-fold** on an identical target, 0.0312 →
   **0.7500**, reported as that and not as a gate, since `73baa16` withdrew its premise before the
   number landed. That withdrawal was right for the reason given: the agent's own log showed the bed
   penetration that makes a whole-body translation inadmissible.

   The second half of the repair was earned by measurement rather than assumed — the bound change
   enters through the **stiffness**, and the bounds that response violates are carried by an
   **active set inside that solve** rather than clipped onto the answer afterwards:

   | start of an increment | distance from the rigid translation | inverted elements |
   |---|---:|---:|
   | clip alone (the bug) | median 7.48 mm | **9,865** |
   | linear response alone | **0.0000 mm** median and max | **0** |
   | response, then clipped | max 7.34 mm on 79 one-sided nodes | 185 |

   **R′ is the one that counts, and it fails at 0.6250 on a single element.** Its first increment is
   *exactly* right — one Newton iteration, `min J` **1.000**, the translation reproduced — then
   `min J` falls 1.000 → 0.710 → the 0.2 floor. **The solver finds an exact zero-energy solution and
   then drifts off it.** With the bed removed that solution provably exists, so this is not a
   feasibility question dressed as a solver question, which is the mistake that cost Control R its
   premise. Here the premise is **demonstrated by increment 1** rather than asserted.

   **A guess tested and refuted, worth recording because it inverts the obvious reading.** The
   association pinned at 0.4994 mm every pass against a 0.5 mm jump limit while nodes moved 0.94 mm
   per increment, so updates were being refused and the constraint normals went stale. Raising the
   limit to 5 mm made it **worse**: R′ then stalled at fraction **0.0010**, association jumping
   3.9 mm per pass, 16.4 mm gap error. **The limit was masking an unstable association, not causing
   it.**

   **Pre-registered before it is built: freeze the frames within a load step**, re-associating only
   *between* steps, which makes each step a fixed-constraint problem whose zero-energy answer under
   R′ is the translation itself.

   * **Gate R″:** R′ with frames frozen within each step completes at fraction **1.0**, zero
     inversions, `min J ≥ 0.99` at **every** increment. Unlike R, this premise is not asserted —
     increment 1 already reaches `min J` 1.000, so a findable zero-energy solution is demonstrated.
   * **Localise the drift in time, which a pass/fail cannot:** report `min J` immediately **before
     and after each re-association**. Falling at the re-association boundaries implicates the
     association; falling *within* a step whose constraints are fixed puts the cause elsewhere and
     makes the association a red herring twice over.
   * **No prediction offered.** My last four calls across these two lines were one hit, two refusals
     and one premise found false before its own result. The design is worth more than the guess, and
     both outcomes of the before/after `min J` reading are informative.

   **What this does not license.** The breast stays untouched. Outcome 3 — "the field is infeasible
   for this mesh" — requires a stepping that can carry a motion known to be feasible, and R′ says we
   do not have one yet.

   **GATE R'' FAILS at the same fraction as R', and its reading exonerates the association.**
   With the frames frozen inside each load step and re-association only between steps, the rigid
   drive with the bed removed still stalls at **0.6250** on one element. But freezing changed the
   TRACE, which is what the gate was built to read. Re-association moves no nodes, so each step's
   entry min J is the previous step's exit:

   | step | fraction | min J at entry | min J at exit | change inside the step (constraints fixed) |
   |---:|---:|---:|---:|---:|
   | 1 | 0.125 | 1.000 | 1.000 | 0.000 |
   | 2 | 0.250 | 1.000 | **0.710** | **-0.290** |
   | 3 | 0.375 | 0.710 | 0.712 | +0.002 |
   | 4 | 0.500 | 0.712 | 0.710 | -0.002 |
   | 5 | 0.625 | 0.710 | 0.712 | +0.002 |

   **Every bit of the loss happens inside one step whose constraints were fixed, and none of it across
   a re-association.** Under R' the same drive fell 1.000 -> 0.710 -> 0.307 -> 0.275 -> 0.226 to the
   floor; frozen, it falls once and then holds to within 0.002 for three more steps. So freezing the
   frames stopped the progressive drift and did not stop the stall, and the association -- suspected
   through the jump-limit test and again here -- is not the cause of either. It is exonerated twice.

   What remains is sharper than "the stepping is broken". Step 1 solves the increment exactly: one
   Newton iteration, min J 1.000, the translation reproduced to 0.0000 mm. Step 2 is the same problem
   with the same fixed constraints and an exact zero-energy solution available -- and the solver takes
   142 Newton iterations to arrive somewhere else, at min J 0.710, and never recovers. A single
   discrete event, not an accumulation: the iteration leaves the zero-energy branch in one step and
   the subsequent steps merely carry that state along. The breast stays untouched.

   **GATE S fails at min J 0.7098, and the evidence supports NONE of the three named causes. The
   cause is the re-association, which R'' appeared to clear.**

   | cause | test | result |
   |---|---|---|
   | 1, the Hessian projection | step 2 with `project=False` | **cleared**: the same converged state, E 3.664278e-04 and min J 0.7098, in 13 Newton iterations instead of 142. The projection costs speed, not correctness. |
   | 2, the line search | total energy at every accepted step | **cleared**: monotone within the solve (4.82e-04 -> 3.66e-04). The one energy RISE is across the step boundary -- the start construction -- not an accepted step. |
   | 3, a second branch in the formulation | step 2 with the association untouched | **cleared**: the same step, same solver, same bounds arithmetic converges to **min J 1.0000** and the exact translation, error **0.0000 mm**. No second branch: the zero-energy solution is found when the constraints are consistent. |

   **What it is.** Step 2's first Newton iteration already begins at min J 0.4245 -- the state is
   deformed before the iteration runs -- and Newton then IMPROVES it to 0.7098. The damage is in the
   start, and the start is built from the re-associated bounds. Measured directly, after one 0.935 mm
   rigid step: 5 of 3,123 held nodes (0.2%) lose the bed entirely, their ray no longer meeting it
   within the 60 mm reach; for the rest the closest point moves a median 0.4499 mm but up to
   **46.6 mm**; and the persistence rule then freezes **1,455 of 3,123 (46.6%)** of the associations
   while the others update. Half the constraint set then refers to new bed points and half to old, and
   the difference enters as a bound change of up to **5.80 mm** for a step whose true motion is
   0.935 mm. The solve lands 4.79 mm from the exact answer.

   **This corrects the previous entry.** R'''s reading -- all the loss inside a step, none across a
   re-association -- was measured correctly and interpreted too coarsely. The re-association moves no
   nodes, so nothing changes AT the boundary; its effect materialises in the next step's start, which
   is inside the step. The cause sits at the boundary and the symptom appears after it. The
   association is not exonerated: it is the cause, and the jump limit that masked it is part of the
   mechanism rather than a protection against it.

   **The coarse interpretation was mine, and it was made twice.** `04f32a0` took "all loss inside a
   step, none across a re-association" and concluded the association was "cleared twice, once by a
   refuted guess and once by design". The measurement was right and the inference was wrong, for a
   reason I should have seen when designing the diagnostic: **a re-association moves no nodes, so it
   cannot show up at the boundary — only in the next step's start.** I built a before/after reading
   around a boundary and then read the absence of a jump at that boundary as absence of a cause.
   And `73baa16` compounded it, recording the jump limit as "masking an unstable association, not
   causing it"; freezing 47% of constraints while updating the rest **is** the inconsistency, so the
   limit is a mechanism, not a bystander. Both of my three named causes for gate S, and my reading
   of R″, were wrong; the agent's direct measurement of what the bounds actually do is what settled
   it.

   **Pre-registered before the fix: consistency is all-or-nothing per step.** Two arms, both
   internally consistent, differing only in whether staleness accumulates:

   * **T-none** — no association update for the whole drive.
   * **T-all** — every association updated at every step boundary.
   * **Gate T:** R′ completes at fraction **1.0**, zero inversions, `min J ≥ 0.99` at every
     increment.
   * **The 5 nodes that lose the bed need a declared rule** — released, or held at their last valid
     association — stated in the artefact with the count reported per step. Chosen silently, it
     would be the next borrowed constant, which this programme has just paid for once in `CORR_TRIM`.
   * **A structured expectation rather than a point prediction:** T-none should pass, a rigid
     translation requiring no re-association at all. T-all may fail, some associations genuinely
     teleporting 46.6 mm. **If T-none passes and T-all fails, "update everything" is not the fix for
     real seating** — where the motion is not rigid and the bed is present — and a teleport-rejection
     rule is needed, which is exactly where the 46.6 mm outliers live. That is the outcome worth
     preparing for, because it is the one that does not transfer from this control to the breast.

   **GATE T: T-none PASSES, T-all FAILS, and the solver is sound.**

   | arm | association | result |
   |---|---|---|
   | **T-none** | never updated after the first | **PASS**: fraction 1.0, min J **1.0000**, zero inversions, the exact rigid translation -- median and max \|u - d\| **0.000 mm**, volume ratio 1.000000 -- in 36 s, **one Newton iteration per step**, gap error 0.0000 mm at every step |
   | **T-all** | every association updated at every step | **FAIL**: stalls at fraction 0.1250 with **548** inverted elements, worse than the mixed rule's 0.6250 |

   **The stepping and the solver are sound.** Given a constraint set that stays internally consistent,
   R' completes exactly: nothing is lost anywhere, every increment is solved in one Newton iteration,
   and the answer is the known one to 0.000 mm. That is the claim gate R was built for and could not
   make, and it is now measured rather than assumed.

   **And "update everything" is not the fix.** T-all is worse than the per-node jump limit it was
   meant to replace -- 0.1250 against 0.6250 -- because some associations genuinely teleport 46.6 mm
   after a sub-millimetre step, and taking those updates wholesale wrecks the constraint set in one
   move rather than in halves. Both extremes fail for real seating: 'none' cannot follow a bed that a
   non-rigid motion slides along, and 'all' imports the teleports.

   **Declared rule for the nodes that lose the bed:** HOLD at the last valid association, with the
   count reported per association (T-none: [0], since it associates once). The alternative --
   RELEASE, dropping the constraint for that step -- is implemented and named in the code, and was not
   used here: releasing would change WHICH nodes are driven midway, so the two arms would have
   differed in more than staleness.

   **What follows for the breast.** Gate S showed the damage comes from MIXING updated and stale
   constraints; gate T shows that neither pure extreme works for a motion that needs the bed followed.
   The rule that satisfies both readings is all-or-nothing PER STEP: if any association would move
   further than a declared threshold, refuse the whole update for that step and keep every previous
   association; otherwise take them all. Every step then has an internally consistent constraint set,
   which is the property the solver needs, and no step mixes epochs. That threshold is a number that
   must be pre-registered rather than chosen by analogy.

   ### The threshold, derived rather than chosen — no free parameter

   > **threshold = |the step's node motion| + the bed's facet scale**

   Both terms are already measured on this line and neither is picked. The node motion is known per
   step — **0.935 mm** in this drive. The facet scale is **1 mm**, declared here earlier as a
   convergence tolerance because the muscular wall's normals turn several degrees over 2–4 mm. The
   physical content: a closest point on a locally smooth bed cannot outrun the node that owns it by
   more than the surface's own resolution, so anything further has **changed feature**, which is
   what a teleport is.

   On this drive the bar is **1.935 mm**. The median association move of **0.4499 mm** sits well
   inside it; the **46.6 mm** moves exceed it **twenty-four-fold**. It separates the two populations
   by a factor of 24, not by a judgement call — which is the property `CORR_TRIM = 0.10` never had.

   **Gate U, fixed before implementation:** under this rule R′ passes exactly as T-none does —
   fraction 1.0, zero inversions, `min J ≥ 0.99` — because some associations exceed the bar at every
   step, the whole update is refused, and T-none's behaviour is recovered. **That is a consistency
   check on the rule, not a new result, and must be reported as such.** A rule that failed it would
   be wrong; a rule that passes it has only shown it did not break what already worked.

   **The real test is the refusal rate, and it is where this rule can fail.** Report the fraction of
   steps whose update is refused.

   * **Refusal rare** — the rule follows the bed while staying consistent, which is what real
     seating needs.
   * **Refusal near 100%** — the rule has degenerated to T-none, which cannot follow a sliding bed.
     A single all-or-nothing switch then cannot deliver both consistency and bed-following, and
     saying so is more useful than tuning the bar.
   * **In that case the answer is the step size, not the rule:** shrink the increment until no
     association exceeds the threshold, making a consistent update always available. That is
     adaptive step control driven by **association motion** — the quantity that actually breaks —
     rather than by `min J`, and it follows from the same measurement rather than introducing a new
     constant.

   **No prediction on the refusal rate.** The breast's motion is non-rigid and slides along the bed,
   which is exactly the regime neither this control nor its predecessors exercise, and I have been
   wrong twice today on this line by reasoning past the end of what the controls actually cover.

   **GATE U: the consistency check passes and the REFUSAL RATE is 93.8%. The single switch has
   degenerated to T-none, and the answer is the step size.**

   Threshold, derived and not chosen: the step's own node motion plus the bed's 1 mm facet scale --
   1.935 mm on this drive, against a median association move of 0.4499 mm and teleports of 46.6 mm.
   It separates the two populations by a factor of 24.

   **Consistency check (not a new result):** under the rule R' completes to fraction 1.0 with min J
   1.0000, zero inversions and the exact rigid translation, 0.000 mm -- exactly as T-none does, which
   is what a correct rule must do and shows only that it did not break what already worked.

   **The real number: 15 of 16 updates refused (93.8%).** Some association exceeds the bar at nearly
   every step, the whole update is refused, and the rule reduces to never updating. So a single
   all-or-nothing switch cannot give both consistency and bed-following, stated plainly rather than
   tuned: lowering the bar refuses more, raising it re-admits the 46.6 mm teleports that deformed the
   body in the first place.

   **But the step size can, and this is measured rather than argued.** For the same rigid drive, with
   the association compared against the same threshold:

   | rigid step | threshold | association move, median | over the bar | lost the bed |
   |---:|---:|---:|---:|---:|
   | 0.007 mm | 1.007 mm | 0.0036 mm | **0** | 0 |
   | 0.075 mm | 1.075 mm | 0.0361 mm | 6 | 0 |
   | 0.374 mm | 1.374 mm | 0.1797 mm | 25 | 3 |
   | 0.935 mm | 1.935 mm | 0.4499 mm | 98 | 5 |

   There is **no floor**: the count over the bar falls to zero as the increment shrinks, and the nodes
   that lose the bed go with it. A consistent update is therefore always available at a small enough
   step, so the next instrument is adaptive step control driven by ASSOCIATION MOTION rather than by
   min J -- the quantity that actually breaks -- and it needs no constant that is not already
   measured here. It is not built: it needs its gate first.

   **Gate V, fixed before implementation.** Shrink the increment until no association exceeds the
   threshold, then take the whole update. R′ must complete at fraction **1.0**, zero inversions,
   `min J ≥ 0.99` — **and the refusal rate must be ≈ 0, with associations actually updating.** That
   conjunction is the gate. Completion alone would be passed by T-none with extra steps, which is
   not the point: the aim is bed-following *with* consistency, and a run that completes while
   refusing every update has achieved only the half already in hand.

   ### The regime nothing has tested: every control so far has been a RIGID motion

   The agent's own caveat, and the sharpest thing in this round. R, R′, R″, S, T, U are all rigid
   drives. **The breast's motion is non-rigid and slides along the bed** — which is exactly where
   large *legitimate* association motion lives. In a rigid drive there is almost none, so **no
   control run so far can distinguish "following a bed" from "teleporting off it"**. Both my last
   two inferences and the agent's were made about sliding, from controls containing none.

   **Gate W, pre-registered: an analytic PLANAR bed, driven by a rigid translation TANGENT to it.**
   A plane makes the answer exact — a tangential rigid translation is zero-energy, the body does not
   deform, and every association slides by **exactly the translation distance** with **zero**
   teleports. It is the minimal control that exercises bed-following, and its answer is known in
   closed form rather than by comparison with another arm.

   * **Pass:** fraction 1.0, `min J ≥ 0.99`, refusal rate 0, and association motion matching the
     translation distance to the solver's tolerance.
   * **Fail:** the rule cannot follow even a perfectly smooth bed — decisive, and it would move the
     problem from the threshold to the association's construction.
   * **No result about the breast is readable until a control exercises the regime the breast is
     in.** That is the standing condition on this line now, and it outranks the step-control work:
     adaptive stepping that passes on rigid drives would still be untested where it matters.

   **GATE W PASSES: the rule follows a bed. The regime every earlier control missed is now covered.**
   (`scripts/gate_w_planar_bed.py`.) An analytic plane, driven by a rigid translation TANGENT to it,
   with the same held set and the all-or-nothing rule under test:

   | | result |
   |---|---|
   | completion | fraction 1.0, **one Newton iteration per step** |
   | min J | **1.0000**, zero inversions, volume ratio 1.00000000 |
   | **refusal rate** | **0.0%**, 0 of 16 |
   | association motion per step | **0.9350 mm at every step**, against a translation distance of 0.9350 mm |
   | solution vs the exact tangential translation | median 0.000000 mm, max 0.000001 mm |

   One design point, because a normal-only constraint cannot express this test: under the sliding
   boundary condition a TANGENTIAL translation is driven by nothing -- u = 0 already satisfies
   n.u = 0 -- so the base is driven in all three components here. That is what makes the association
   follow the bed rather than sit still, and it is why every previous control, all of them normal-only
   rigid drives, contained almost no legitimate association motion to test.

   **What it separates.** The rule now has both readings on the record: on a smooth bed it follows
   sliding perfectly, refusing nothing while associations slide by exactly the distance the body
   moves; on the anatomical bed it refuses 93.8%. So the refusals are a property of THAT BED's
   features -- the sheets a closest point jumps between -- and not of the rule's ability to follow
   motion. Before gate W neither this line nor its coordinator could have told those two apart,
   and three wrong readings came from reasoning about sliding using controls that contained none.

   ### A weakness in my derived threshold that gate W could not expose

   In W **every association moved by the same 0.9350 mm**, because the drive is rigid. So the
   threshold's first term — *"the step's node motion"* — was a single well-defined scalar. **Under a
   non-rigid motion there is no such scalar:** different nodes move different amounts within one
   step, and a global bar is then either too tight for the fast nodes or too loose for the slow ones.
   My derivation quietly assumed the regime it happened to be tested in.

   **The generalisation, pre-registered, and forced rather than chosen:**

   > **threshold_i = |node i's own motion this step| + the bed's facet scale**

   The refuse-or-take decision stays **global** — that is what keeps the constraint set in one epoch —
   while the bar it is measured against becomes per-node. **It reduces exactly to the current rule
   when all nodes move equally**, so it is a generalisation, not a loosening. Gates T, U and W must
   be **re-run under it and must return their existing numbers unchanged**: a change that alters a
   passed control's result is a different rule, not a generalisation of one.

   ### Gate X: non-rigid, sliding, and still exact

   **A homogeneous deformation of the held base against an analytic plane.** Non-rigid and sliding,
   yet closed-form: a homogeneous field is the exact solution for a homogeneous neo-Hookean material
   under compatible boundary conditions — the same class gate R2 already validates against platens.
   It is the first control in which **association motion varies from node to node**, which is the
   only untested property left in the rule.

   * **Pass:** fraction 1.0, `min J ≥ 0.99`, refusal rate 0, solution matching the exact homogeneous
     field to solver tolerance.
   * **Fail:** the per-node threshold does not survive a spatially varying field, and the
     all-or-nothing decision needs rethinking rather than regeneralising.
   * **Ordering, proposed by the agent and adopted:** run gate V on the **plane** first, where any
     refusal is a defect in the stepping, before the anatomical bed, where refusals are the bed's
     doing. That stops the same confusion recurring one level up.
   * **The standing condition is unchanged and X does not satisfy it.** X is on a plane; the breast
     is non-rigid motion against the *anatomical* bed. X is necessary, not sufficient.

   **The threshold is per-node now, T/U/W return unchanged, and GATE X PASSES.**
   (`scripts/gate_x_homogeneous.py`, standalone like W so both re-run as regressions.)

   The bar is `|node i's own motion this step| + the bed's facet scale`, with the refuse-or-take
   decision still GLOBAL -- that is what keeps every constraint in one epoch. Under a rigid drive all
   nodes move equally and it reduces to the previous rule, which is what the regressions verify:

   | regression | before | after |
   |---|---|---|
   | T-none | PASS, fraction 1.0, min J 1.0000, 0 inversions | **identical** |
   | T-all | FAIL at 0.1250, 548 inverted | **identical** |
   | U | 93.8% refused (15 of 16), completes, min J 1.0000 | **identical** |
   | W | PASS, 0.0% refused, min J 1.0000 | **identical** |

   | gate X, homogeneous shear tangent to an analytic plane | result |
   |---|---|
   | completion | fraction 1.0, one Newton iteration per step |
   | min J | **1.0000**, exact for an isochoric shear; zero inversions; volume ratio 1.00000000 |
   | refusal rate | **0.0%**, 0 of 16 |
   | solution vs the exact homogeneous field | median **0.000000 mm**, max **0.000000 mm** |

   Association motion varies from 0 to 7.5 mm across the base here -- the property no earlier control
   contained -- and nothing is refused, so a per-node bar survives a spatially varying field.

   Two notes on how it had to be built. The WHOLE boundary is driven, not the base alone: a
   homogeneous field is exact only if the entire boundary is compatible with it, and driving the base
   alone would let the free surface relax, leaving nothing closed-form to compare against. And X's own
   log prints a "held gap error" rising to 46.55 mm, which is an INAPPLICABLE metric rather than a
   defect: it measures the distance from a normal-gap-closing drive that X does not apply, since X
   prescribes displacement directly. The number is real and means nothing here.

   **X IS NECESSARY AND NOT SUFFICIENT.** It runs on a plane. The breast is non-rigid motion against
   the ANATOMICAL bed, where a closest point jumps between sheets, and the standing condition on this
   line is unchanged: no result about the breast is readable until a control exercises that regime.
   A passing X is not clearance.

   ### The standing-condition control has no analytic answer. Two known answers that need none.

   On the anatomical bed with a non-rigid motion there is no closed form to compare against, and the
   obvious repair — prescribe the whole boundary, as X does — **removes the bed's role entirely**,
   since a fully prescribed boundary leaves the constraint nothing to do. That is the same bind the
   skin line hit: the regime that matters is the one where no truth exists. Two known answers do not
   require one.

   **Gate Y — reversibility.** Drive the base forward through the anatomical bed by a non-rigid
   motion, then drive it back. For frictionless unilateral contact on a hyperelastic body the
   quasi-static solution is a property of the configuration, not of the path taken to it, so the
   return must reproduce the start.

   **Gate Z — path independence.** Reach the *same* final configuration by two different step
   sequences (16 increments against 64). The two answers must agree.

   **Z is the sharper of the two**, because the failure mechanism this line has identified is
   *association staleness*, and staleness is exactly what a change of step count varies while holding
   everything else fixed. A stale-association error is a path-dependent error by construction.

   **Tolerance, derived not chosen: the bed's facet scale, 1 mm** — the same measured quantity the
   threshold already uses. The association is discrete, so two step counts may legitimately select
   different bed facets and differ at that scale; anything larger is path dependence the
   discretisation cannot excuse.

   * **Both pass:** the stepping is path-independent on the real bed under a real motion, which is
     the property every earlier control was built to approach and none could reach.
   * **Either fails:** the failure is quantified in millimetres on the geometry that matters, and it
     is the first such number this line has had.
   * **Honest limit, recorded now so a pass is not over-read:** path independence and reversibility
     are **necessary, not sufficient.** A consistently wrong answer is reversible and path
     independent too. They bound the *solver's* contribution to the error; they say nothing about
     whether the seated breast is anatomically right, which remains judged by gates (a)–(d) and by
     nothing else.
   * **Expectation, and it is the agent's rather than mine:** refusals reappear, since the anatomical
     bed produced 93.8% under a *rigid* drive while X's 0% came from a surface with no features to
     jump between. If they do, gate V's adaptive stepping is not an optimisation but the mechanism
     that makes this regime tractable at all.

   **Gates Z and Y are built and running; no numbers yet, and the runtime is itself the first
   measurement of this regime.** (`scripts/gate_zy_anatomical.py`, standalone like W and X.)

   Z drives the same non-rigid motion on the ANATOMICAL bed in 16 increments and in 64 and requires
   the final configurations to agree; Y drives forward and back in one schedule and requires the
   return to reproduce the start. Tolerance for both is the bed's facet scale, 1 mm -- the same
   measured quantity the refusal threshold uses -- because the association is discrete and two step
   counts may legitimately select different facets at that scale.

   **What is already measured: cost.** The 16-increment leg alone has run **37 minutes without
   completing**, against **36 s** for the rigid drive on the same mesh (T-none), seconds for the
   planar gate W and 30 s for the homogeneous gate X. The regime the breast is actually in is at least
   60x more expensive per leg than every control that has been used to reason about it, and Z needs
   three such legs. That is a property of a non-rigid motion against a bed whose association is
   contested at nearly every step, and it is worth knowing before anything is planned around a full
   seating run.

   **Stated in advance, as required, so a pass cannot be over-read:** reversibility and path
   independence are **NECESSARY, NOT SUFFICIENT**. A consistently wrong answer is reversible and
   path-independent too. They bound the SOLVER's contribution to the error and say nothing about
   whether the seated breast is anatomically right, which is judged by gates (a)-(d) and by nothing
   else.

   ### The cost of the real regime is itself the first measurement of it

   | control | regime | time |
   |---|---|---:|
   | T-none | rigid, no bed | **36 s** |
   | W | rigid tangential, analytic plane | seconds |
   | X | non-rigid, analytic plane | **30 s** |
   | **Z, first leg** | **non-rigid, anatomical bed** | **>37 min, unfinished** |

   **The regime the breast is actually in is at least 60× more expensive per leg than every control
   used to reason about it**, and Z needs three legs. Nothing in the control series gave any warning,
   because each control removed exactly the thing that costs: W and X replace the bed with a plane,
   which has no features for an association to contest, and T-none removes the bed entirely.

   **The cost is not incidental — it is the same finding in another currency.** The anatomical bed
   refused **93.8%** of association updates under a *rigid* drive. A non-rigid motion contests them
   at nearly every step, and a solver fighting a constraint set that will not settle spends its
   iterations there. **So gate V's adaptive stepping is not an optimisation and not merely a
   correctness fix: it is plausibly what makes this regime affordable at all.** That reverses the
   priority set two entries ago, where V was placed behind the sliding controls.

   **Ordered now, because it costs nothing to find out:** run **gate V on the analytic plane** — where
   X completes in 30 s, the answer is closed-form, and any refusal is a defect in the stepping rather
   than the bed's doing — **in parallel with Z** rather than after it. If adaptive stepping works on
   the plane it is a candidate for making Z affordable; if it does not, we learn that in seconds
   instead of after hours.

   **Z still runs on the current stepping**, because it measures the path-dependence of what exists
   today, and because a V that changed the stepping would require Z to be re-run regardless. They are
   not competing for the same answer.

   **A silent long job cost a scheduling decision for the second time today.** The Z/Y driver passed
   no log, so 37 minutes of compute produced no visible progress; it logs per increment now. The first
   instance was the v3 anchored fit running 4h34m without output, which forced a keep-or-kill call to
   be made blind. Same class as the wording bugs this line has removed — **not a wrong number, but no
   number where one was needed** — and on today's evidence it is the more expensive of the two.

   **GATE V PASSES ON THE PLANE, and the harness nearly reported a false pass.**

   | gate V, adaptive stepping on the plane, closed-form field | result |
   |---|---|
   | completion | fraction 1.0, min J **1.0000**, zero inversions |
   | steps rejected | **0 of 8** |
   | associations actually updated | **9** -- the conjunction, since a rule that never updates shows 0 |
   | solution vs the exact homogeneous field | median **0.000000 mm**, max **0.000000 mm** |

   Adaptive stepping rejects the STEP and shrinks it when an update would exceed any node's own bar,
   rather than refusing the update and carrying stale constraints forward. The step size is then
   driven by association motion -- the quantity that actually breaks -- instead of by min J.

   **The near-miss, recorded because it is the same class as the silent 37-minute run.** The first
   `--adaptive` run reported a pass while its own header read `stepping: allornothing`: the script
   sets `sys.argv = ["x"]` before loading the seat module, which wiped the flag before the mode was
   read. It was checked rather than assumed, and the flag is now read before anything overwrites
   argv. The adaptive path was then verified to DISCRIMINATE rather than merely run: with the derived
   1 mm facet scale it completes with zero rejections, and with a 0 mm facet -- no slack at all -- it
   rejects and stalls, which is what a live mechanism must do.

   **Still only the plane.** V has not been run on the anatomical bed, which is where it might make
   the regime affordable, and gate Z continues on the current stepping so that today's path
   dependence is measured for what exists rather than for what replaces it.

   **GATE V FAILS ON THE ANATOMICAL BED, and the failure has a floor that the rigid drive did not.**

   | gate V, anatomical bed, the smoothed seating drive | result |
   |---|---|
   | completion | **stalled at fraction 0.0002** after 705 s |
   | furthest fraction reached | 0.0625, then backwards |
   | rejection rate | **48.6%** (17 of 35 proposals) |
   | associations over their own bar at a rejection | **median 1**, min 1, max 79, of 3,123 held |
   | increment at the stall | below 1e-6 of the drive |

   **Adaptive stepping does not make this regime affordable**, and the reason is that median of ONE.
   The bar is |node motion| + 1 mm, so as the increment shrinks the bar tends to the facet scale; an
   association that jumps more than 1 mm under an ARBITRARILY SMALL motion can never be accommodated
   by a smaller step. On the rigid drive the over-bar count fell to zero as the step shrank -- that
   was the "no floor" measurement -- and under the seating drive on this bed it does not. The earlier
   table generalised from a rigid motion, which is the same failure this line has now made at four
   levels: the controls removed the feature that costs.

   So the affordability question is answered and the answer is no: 705 seconds bought 0.02% of one
   breast's drive, against gate Z's current stepping which passed 42 minutes unfinished on the same
   drive. Neither stepping makes a full seating run reachable, let alone eight breasts.

   One log line was corrected rather than reported: in adaptive mode the per-step line printed
   "association moved 0.0000 mm", which is the PRE-step call that by construction returns zero, not
   evidence that associations fail to update. It now says the association is held for the step.

   ### My "there is no floor" was measured on a rigid drive and does not survive the seating drive

   `173bc9e` recorded a step-size table in which the over-the-bar count fell to **zero** as the
   increment shrank, and concluded: *"There is no floor... a consistent update is always available at
   a small enough step."* That table was measured on a **rigid** drive. Under the seating drive on
   the anatomical bed there **is** a floor, and the mechanism is arithmetic:

   > the bar is `|node motion| + 1 mm`, so **as the increment shrinks the bar tends to the facet
   > scale**. An association that jumps more than 1 mm under an *arbitrarily small* motion can never
   > be accommodated by a smaller step.

   Gate V's own numbers show it: at a rejection, the count of associations over their own bar has a
   **median of 1** — one node out of 3,123, jumping more than the facet scale no matter how small the
   step, refuses every proposal. **That is a discontinuity, and no step size fixes a discontinuity.**

   **So the affordability question is answered and the answer is no.** 705 s bought **0.02%** of one
   breast's drive (stalled at fraction 0.0002, having reached 0.0625 and gone backwards, increment
   below 10⁻⁶); Z's current stepping has passed 56 minutes unfinished on the same drive. **Neither
   stepping makes a full seating run reachable, let alone eight breasts.**

   **This is the fourth level at which the same error has been made on this line**, and it is mine
   this time: every control removed the feature that costs, and I generalised from the control to the
   regime. Rigid drives have no feature-jumping, so they have no floor; the anatomical bed under a
   non-rigid motion has both.

   ### The fix is not a threshold, it is the association's construction

   A global closest-point query on a faceted, multi-sheet bed **can** jump between sheets, and every
   rule tried so far — per-node limit, all-or-nothing, derived threshold, adaptive stepping — has
   been an attempt to *detect* or *tolerate* that jump after the query made it. **Pre-registered
   instead: constrain the query so the jump is impossible.**

   > **Gate AA — local-search association.** The new association must lie within a bounded geodesic
   > neighbourhood of the previous one on the bed surface, rather than being the global closest point.
   > A jump between sheets is then impossible by construction, the refusal machinery becomes
   > unnecessary, and step size is governed by the physics rather than by association motion.

   * **Regressions, which must return their existing numbers unchanged:** T-none, W, X, and V on the
     plane. A construction that alters a passed control's result is a different method, not a
     better association.
   * **Discrimination, per `CLAUDE.md`:** force the neighbourhood to zero and the behaviour must
     become T-none's — never updating. A pass that looks identical to the mechanism's absence is not
     evidence.
   * **Affordability, the number that decides it:** wall-clock and fraction reached on the same
     seating drive, against V's 0.02% in 705 s.
   * **The failure mode to report rather than hide:** a local search can hold the *wrong* sheet when
     the true contact genuinely moves to another one. Count how often the local best is worse than
     the global closest point, and by how much. If that count is large, the constraint is buying
     tractability by getting the contact wrong, which is worse than being slow.
   * **No prediction.** Every prediction I have made about this regime ahead of a measurement has
     been wrong, most recently the floor.

   ### Gate AA fails, and the SHAPE of its failure says the stepping was never the problem

   | attempted fraction | step cut back to | inverted elements **at the start of the increment** |
   |---:|---:|---:|
   | 0.1250 | 0.06250 | **44** |
   | 0.0625 | 0.03125 | 38 |
   | 0.0312 | 0.01562 | 39 |
   | 0.0156 | 0.00781 | 41 |
   | 0.0078 | 0.00391 | 43 |
   | 0.0039 | 0.00195 | 43 |
   | 0.0020 | 0.00098 | **42** |
   | 0.0010 | 0.00049 | **42** |
   | 0.0005 | 0.00024 | **42** |
   | 0.0002 | 0.00012 | **42** |
   | 0.0001 | 0.00006 | **42** |

   **The count does not move.** Across a 2,000-fold reduction in step size it sits at 42, and the
   local-search association cannot be the cause because it is constrained not to jump. **An inversion
   count that is invariant to the step size is not caused by the step.**

   **And as the attempted fraction tends to zero, the configuration tends to the starting one.** The
   straightforward reading is that **the starting configuration for this drive already contains
   about 42 inverted elements**, and every stepping, threshold and association scheme tried on this
   line — R, R′, R″, S, T, U, V, AA — has been trying to step away from a state that was already
   broken.

   **That is a hypothesis with a one-command test, pre-registered here before it is run:** evaluate
   `det J` on the seating drive's **initial bound configuration**, before any increment is applied,
   and count the non-positive elements.

   * **~42 inverted at the start** → the stepping was never the problem, and neither was the
     association. Every control on this line was measuring its own ability to escape a broken
     initial state, which no amount of stepping can do. The work moves to how that state is built.
   * **Zero inverted at the start** → the inversions are created by the first increment and then
     persist regardless of its size, which is a different and stranger finding, and the stepping
     question reopens on new terms.
   * **T-none is the reason this was not obvious**: it reached `min J` 1.0000 with zero inversions,
     so the mesh is sound *in that configuration* — rigid drive, no bed. Whatever is wrong appears
     when the bed's constraints are applied to the seating drive's starting pose, not in the mesh
     itself.
   * **No prediction**, and for a specific reason: I have now been wrong about this regime at four
     levels, each time by generalising from a control that had removed the thing that mattered.

   **The agent's own honesty about what AA could NOT measure is worth keeping**: "at the starting
   positions the local association IS the global one, so a comparison there would read 0 by
   construction." The wrong-sheet price I asked for is unmeasurable from a run that never left the
   start — reporting that, rather than a reassuring zero, is what makes the rest of the run readable.

   **THE PRE-REGISTERED TEST COMES BACK ZERO. The starting configuration is sound, and the second
   branch is the one we are in -- say it loudly.**

   | configuration | non-positive elements | min J |
   |---|---:|---:|
   | the undeformed mesh | **0 of 54,704** | 1.000 |
   | the bounds imposed at fraction 0, 1e-6 and 1e-4, initial association | **0, 0, 0** | 1.000, 1.000, 0.997 |
   | the same after one association update | **0, 0, 0** | 0.999 |
   | the clip alone, at fractions 1e-2 down to 1e-6 | **0 at every fraction** | -- |
   | the clip plus the linear elastic response | **0 at every fraction** | -- |

   So the seating drive does NOT begin from a broken state, and neither the clip nor the response
   that replaced it inverts anything from the pristine configuration.

   **AA's 42 inversions were my script's defect, and AA's numbers are withdrawn as evidence about
   local search.** Called twice at IDENTICAL positions, the local association moves 141 of 3,123
   nodes by more than 1 mm and up to 6.5 mm, then stabilises -- a one-time disagreement with the
   global query rather than drift. The cause is the starting face: it is assigned by CENTROID
   PROXIMITY, and on a bed whose faces run 1.6 mm median but 28.9 mm maximum the nearest centroid is
   routinely not the face the association sits on, so the neighbourhood is grown from the wrong
   place. That is the third appearance of the same two-tier hazard -- after the ray index and the
   radius test -- and it means the step-size-invariant count of ~42 measured the script, not the
   method.

   Incidentally that disagreement IS the price metric gate AA was asked to report: 141 nodes, up to
   6.5 mm. It is unusable as a measurement of local search until the starting face is taken from the
   ray query that found it rather than inferred from a centroid.

   **What remains open.** The earlier seating stalls under the GLOBAL association -- 13 and 2 inverted
   elements, from states reached after accepted steps -- are not explained by this defect and are
   still unaccounted for. What is now established is only that they do not originate in the initial
   configuration.

   **THE SEED IS FIXED, AA IS HONEST, AND THE OPEN QUESTION HAS AN ANSWER: A STEP OF SIZE ZERO FAILS
   IDENTICALLY.**

   *1. The idempotence gate, run before any AA number is quoted.* The local search now takes its
   starting face from the face the RAY QUERY returned, reconstructed from the same primitives
   bed_rays uses and checked against it as a known answer -- agreement 5.55e-17 m in the point and
   2.90e-14 in the normal over 4,097 of 4,097 hits. A second defect was found in the same object
   while fixing the first: the seeding call took its normal from bed_rays' area-weighted vertex
   normals and every later call from vertex_normals(), two different fields, so a node could be
   handed a different normal at ZERO motion. Both paths now use one field and one interpolation.

   | two calls at identical positions | faces changed | point moved | normal moved |
   |---|---:|---:|---:|
   | at the initial positions | **0 of 4,097** | 0.000000 mm | 0.00e+00 |
   | at 0.5 mm along each ray | **0 of 4,097** | 0.000000 mm | 0.00e+00 |
   | back at the initial positions | **0 of 4,097** | 0.000000 mm | 0.00e+00 |

   *2. Gate AA re-run.* It reaches **fraction 0.0625** -- 6.25% of the drive, against gate V's 0.02%
   on the same bed -- and then stalls with **48 inverted elements, invariant from ds = 0.023 down to
   ds = 0.00009, a 256x range**. So the step-invariance is REAL and was not only my artefact; what
   was my artefact was the particular number 42 and the claim that it said anything about local
   search. **THE PRICE, now measurable because the stall hook reports at a state that has actually
   moved**: the local best is worse than the global closest point for **124 of 4,096 nodes, by a
   median 0.336 mm and at most 13.831 mm**. A 13.8 mm miss is a held wrong sheet, exactly the
   failure mode the gate was built to expose.

   *3. The open question, answered.* At the stall the drive was asked for a step of size ZERO -- the
   same fraction it had already accepted, so no held node is asked to travel at all.

   | request at fraction 0.0625 | inverted | bounds move >1 mm | worst bound move |
   |---|---:|---:|---:|
   | ds = 0.023 | 4 | 39 nodes | 6.80 mm |
   | ds = 0.0029 | 4 | 25 nodes | 7.14 mm |
   | ds = 0.00073 | 4 | 25 nodes | 7.18 mm |
   | **ds = 0 exactly** | **4** | **26 nodes** | **7.1912 mm** |
   | ds = 0 again, same state | 4 | 26 nodes | 7.1912 mm |

   As the step shrinks 128x the bound displacement does not shrink; it CONVERGES to 7.19 mm. A step
   of size zero reproduces the failure exactly, twice.

   **THE MECHANISM, derived rather than guessed.** A held node's bound is
   `on_plane + gap0 + fraction * travel`, where `on_plane = n . (c - X[base])` is the displacement
   that puts the node on its association plane. At a zero increment `gap0` and `fraction * travel`
   are bit-identical to the accepted step's, so the entire 7.19 mm is `on_plane`: the association
   term. **The load fraction scales one of the two terms in the prescribed displacement and not the
   other.** The association term is recomputed absolutely from the current association and applied
   instantaneously, so once the association has drifted, a demand of several millimetres is waiting
   at the head of the next increment whatever its size. That is why no stepping rule helped, why
   adaptive stepping failed on this bed, and why every test from the pristine state read zero -- at
   the start the association has not drifted, so the uncontrolled term is zero.

   **A misleading sentence of my own, caught and corrected inside the round.** The first version of
   this message reported one number for "the bounds displaced the previous configuration". It read
   5,660 nodes when only 4,097 are bounded at all, which is what exposed it: by that line `v` has
   been overwritten by the active-set loop, so the quantity was the WHOLE start displacement
   including the elastic response on free nodes, not the bound displacement. The message now reports
   the two separately. Fifth instance of a sentence that named more than it measured.

   *Added to the solver:* `seat_on_bed(..., on_stall=)`, a hook called at the stall with the drive's
   own `advance` and accepted displacement. It is what makes a zero-sized step askable at all, since
   every cut-back re-associates and so never tested the step alone.

   ### The root cause: the load fraction scales one term of the bound and not the other

   **A step of size zero fails identically.** Requested at fraction 0.0625:

   | requested ds | inverted | bounds moving >1 mm | worst |
   |---:|---:|---:|---:|
   | 0.023 | 4 | 39 | 6.80 mm |
   | 0.0029 | 4 | 25 | 7.14 mm |
   | 0.00073 | 4 | 25 | 7.18 mm |
   | **0 exactly** | **4** | **26** | **7.1912 mm** |
   | 0 again, same state | 4 | 26 | 7.1912 mm |

   Shrink the step 128× and the bound displacement does not shrink — it **converges** to 7.19 mm.

   **The mechanism, which follows without guessing.** A held node's bound is
   `on_plane + gap0 + fraction*travel`, with `on_plane = n·(c − X[base])`. At a zero increment `gap0`
   and `fraction*travel` are bit-identical to the accepted step's, so **the entire 7.19 mm is
   `on_plane` — the association term.** The load fraction scales `travel` and does **not** scale
   `on_plane`, which is recomputed absolutely and applied instantaneously. Once the association has
   drifted, several millimetres are waiting at the head of the next increment **whatever its size**.

   **That accounts for the whole line.** Why no stepping rule ever helped; why adaptive stepping
   failed on this bed; why every test from the pristine state read zero — at the start the
   association has not drifted, so the uncontrolled term is exactly zero. **Every gate from R to AA
   was adjusting the term that was already small.**

   ### Un-withdrawing my inference, and the part of it that stays wrong

   `5131bd8` withdrew `a9d5fe5` entirely on the grounds that AA's instrument was broken. That was too
   broad, and the agent says so against its own retraction: **the step-invariance is real.** It
   survives the seeding fix and holds across a **256× range** of ds. So:

   * **"An inversion count invariant to the step size is not caused by the step" — stands.** The
     logic was sound and the phenomenon is genuine.
   * **"Therefore the starting configuration is already broken" — remains wrong.** The start is
     provably clean, 0 of 54,704.
   * **What I missed was a third branch**: *not caused by the step* admits *caused by a term the step
     does not scale*, which is neither "at the start" nor "in the first increment". My two branches
     were both about **where in the increment** the inversions were born; the answer was a term that
     is not under the increment parameter at all.
   * **The generalisable move, which I did not make and should have:** when a quantity is invariant
     under a parameter, ask **what in the model is not a function of that parameter** — rather than
     asking where in the parameter's range it lives.

   The agent argues the third branch was not available to write down in advance, and I half agree:
   the frame was productive — it is what made "ask for an increment of size zero" the obvious next
   question — but "or something the fraction does not scale" was writable, and I did not write it.

   ### The repair and its gate, fixed before it is built

   > **Scale the association term with the load fraction, or apply it incrementally rather than
   > absolutely, so that every term of the prescribed displacement is under the increment parameter.**

   * **Gate BB, and it is a known answer rather than a threshold: a step of size zero must move
     nothing.** `ds = 0` must produce zero bound displacement and zero inversions, from any state.
     That the current code returns 7.19 mm is the defect stated as a test.
   * **Regressions unchanged:** T-none, W, X, V on the plane, and R2's block cases.
   * **Then** the affordability question re-opens on honest terms: fraction reached and wall-clock on
     the seating drive, against AA's 6.25% and V's 0.02%.
   * **No prediction.** I have been wrong about this regime at five levels, and right once by an
     argument I then withdrew.

   **AA's honest numbers, now that its instrument is a function**: it reaches **fraction 0.0625**,
   6.25% of the drive against gate V's 0.02%, and the wrong-sheet price it was built to report is
   **measurable and real** — the local best is worse than the global closest point for **124 of 4,096
   nodes, median 0.336 mm, worst 13.831 mm.** A 13.8 mm miss is a held wrong sheet. Local search buys
   a 300× improvement in reach and does sometimes hold the wrong surface; both belong in any decision
   about it.

   **A second way the association was not a function**, found while fixing the first: the seeding call
   took its normal from `bed_rays`' area-weighted vertex normals and every later call from
   `vertex_normals()` — two different fields, so a node could be handed a different normal at zero
   motion. One field and one interpolation now, with idempotence gated at **0 of 4,097 faces changed**
   and 0.000000 mm of point motion across repeated and returning calls.

   ### Gate BB passes on the bounds, and the agent refuses to call it clean

   **The bound displacement at a step of size zero is exactly 0.000e+00 mm over 0 nodes**, with no
   inversions, at both accepted steps probed (fractions 0.0625 and 0.1094). The same request before
   the repair returned **7.19 mm over 26 nodes and 4 inversions**.

   **But the solve still moves at a zero-sized step** — 0.335 mm, then 0.477 mm. The bounds impose
   nothing, yet re-association changes the feasible *set*: where the plane has moved away from a
   node, its unilateral lower bound **loosens**, and the body relaxes into it. That motion is carried
   by the minimisation rather than imposed as a jump, which is exactly the difference that stops
   inversions — but "a step of size zero must move nothing" is not literally satisfied, and the
   agent declined to file it as a pass.

   **It should not be satisfied literally, and the gate as I wrote it was too strong.** A loosened
   unilateral bound *ought* to let the body relax; forbidding that would forbid contact from ever
   releasing. What must be true is not that the motion is zero but that it **converges**.

   > **Gate BB′, replacing the literal form:** apply zero-sized steps **repeatedly at the same
   > fraction**. The motion per step must fall below the solver's tolerance within a small number of
   > repeats — the zero-step must reach a **fixed point**.

   * **Motion decays to tolerance** → the residual is one-time relaxation into a legitimately
     loosened feasible set, and gate BB is satisfied in the only sense that was ever meaningful.
   * **Motion persists across repeats** → the association and the solve are chasing each other at a
     *fixed* load, which is the same non-convergence in a new place and is not fixed by the repair.
   * The two numbers already in hand (0.335 and 0.477 mm) are at **different fractions**, so they are
     not a convergence sequence and must not be read as one — a rise from 0.335 to 0.477 across
     different states says nothing about whether either converges.

   **R2 is the open question and must not be guessed at.** With persistent association its `min J`
   falls **1.000 → 0.636 by fraction 0.75**, where a rigid-translation control should hold 1.000. The
   agent is running the identical stage on pre-repair code (`50e7583`) in a separate worktree to
   separate *the repair* from *pre-existing*. **The distinction that decides it:** R2's earlier pass
   (platens 0.00000%, cylinder 6.382 µm) was measured **without** persistent association. If the
   pre-repair baseline also falls, this is a property of persistent association and not a regression
   of the repair — and R2's earlier pass stands for the configuration it was run in, which is not
   this one. Running the baseline rather than reasoning about it is the right call, and the baseline
   being stuck on its second increment at 25 minutes while the repaired run reached 0.75 is
   suggestive and not yet an answer.

   **Affordability, interim:** the drive is at **fraction 0.1094**, already past the 0.0625 where AA
   stalled, and still advancing. **Regressions so far: W, X, V-on-plane and T-none all PASS** with the
   repair in, T-none reproducing the rigid translation at median and max 0.000 mm, `min J` 1.0000.

   ### R2: the repair did not cause it. Persistent association degrades the control in both codes.

   The pre-repair baseline (`50e7583`, separate worktree) reaches fraction 0.25 at **`min J` 0.710**;
   the repaired run at the same fraction reads **0.920**. **Persistent association degrades the
   rigid-translation control in both codes**, so the repair is not the cause — which is all the
   comparison is asked to establish.

   **The agent's caveat against its own favourable number is right and is adopted:** 0.710 against
   0.920 is one fraction in two runs that took different paths there, not a controlled comparison. It
   suffices to refute "the repair caused it", which needs only the baseline to fall at all. It is
   **not** a measurement of how much better the repair is and must not be quoted as one.

   **So a new question is open: why does persistent association cost anything on a control with a
   known answer?** T-none holds `min J` 1.0000 under maximal staleness — but T-none has the bed
   *removed*, so there is no association to go stale. R2's block cases have platens and a cylinder.

   **The decomposition, pre-registered, and it needs no chosen tolerance.** Stale normals cost
   nothing on a surface whose normal does not rotate, and cost a computable amount on one that does.
   Translating a distance `d` along a cylinder of radius `R` rotates the closest point's normal by
   about `d/R`, so the constraint error goes as **d²/R** — the same quantity the skin line spent a
   day on.

   * **Gate CC-flat, a known answer with no tolerance at all:** the **platens** under persistent
     association must still read **0.00000%** and hold `min J` **1.000**. A plane's normal never
     rotates, so staleness costs exactly zero there. **Any deviation is a defect, not a limitation** —
     and it would mean the degradation is a bug rather than the price of persistence.
   * **Gate CC-curved:** the **cylinder's** degradation must track the analytic `d²/R` staleness
     bound. Report measured against predicted. Much larger means something beyond staleness; matching
     means persistence has a known, quantified price on curved geometry, which is a limitation to
     carry rather than a fault to fix.
   * **This separates a bug from a cost**, which the single R2 number cannot do, and it is why the
     block cases under persistent association are a new gate rather than a re-run of the old one.

   **R2's earlier pass is re-scoped, not retired.** Platens 0.00000% and cylinder 6.382 µm stand for
   the configuration they were measured in, which is **association-free**. Quoting them for the
   persistent-association configuration would be the false-pass shape at the level of what a gate
   covers rather than what a flag set.

   **BB′ is built with its threshold declared before the run** — 6 zero-sized steps at one fixed
   fraction, each from the last one's answer, fixed point at **1 µm**, a thousandth of the bed's 1 mm
   facet scale, and the whole sequence printed either way.

   **THE REPAIR LANDS, BB' FAILS HONESTLY, AND A JUMP LIMIT WAS SILENTLY TURNING PERSISTENT
   ASSOCIATION INTO NO ASSOCIATION.**

   *The repair.* A held node's bound was `on_plane + gap0 + fraction*travel`, and the fraction scales
   `travel` and not `on_plane`. It is now interpolated from where the node IS to where the constraint
   wants it, by `theta` -- the share of the REMAINING drive the increment consumes -- so every term
   is under the increment parameter. At theta = 0 the bound is the node's own coordinate, formed with
   the SAME einsum over the SAME frames the solver uses for `v_prev`, so the two agree bit for bit.

   *Gate BB, the bounds half.* A step of size zero displaces the bounds by **exactly 0.000e+00 mm
   over 0 nodes, with no inversions**, against **7.19 mm over 26 nodes and 4 inverted elements**
   before the repair. The defect is closed as a test.

   *GATE BB' FAILS, and it is the branch that was pre-registered as the bad one.* Six zero-sized
   steps at ONE fraction (0.0625), each from the last one's answer, in mm:

   | repeat | 1 | 2 | 3 | 4 | 5 | 6 |
   |---|---:|---:|---:|---:|---:|---:|
   | the solve moved | 0.335152 | 0.122123 | 0.034554 | 0.036361 | 0.037410 | 0.036119 |

   It decays for three repeats and then **plateaus at about 36 um and stays there**. The declared
   fixed point was 1 um, written into the script before the run. This is not a one-time relaxation
   into a loosened feasible set: **at a FIXED load the association and the solve are chasing each
   other**, and the repair does not address it. The motion is almost entirely on the base (36.1 um
   against 6.0 um in the interior), so it is the contact that is cycling. Reported as a failure.

   *Gate R2 / R-prime under persistent association -- the repair turns a FAILING control into a
   passing one.* Same stage, same data, the only difference being the code:

   | | outcome | min J | against the ideal rigid translation |
   |---|---|---:|---|
   | pre-repair (50e7583, a separate worktree) | **STALLS at fraction 0.6250**, 1 inverted, bounds moved 4 nodes up to 7.3841 mm | -- | -- |
   | with the repair | **completes to 1.0**, 0 inversions, 1866 s | 0.4174 | median 0.256 mm, max 6.325 mm |

   That is the first evidence that the repair is good rather than merely neutral. The earlier
   comparison at one fraction (0.710 against 0.920) is superseded by this and should not be quoted.

   *Regressions, all with the repair in:* T-none PASS (rigid translation reproduced to median and max
   0.000 mm, min J 1.0000), W PASS, X PASS, V-on-the-plane PASS.

   **A JUMP LIMIT THAT SILENTLY DISABLES WHAT IT PROTECTS.** Building gate CC exposed this, and it is
   the most transferable thing in the round. With `jump_limit_m` at 0.5 mm and a per-step slide of
   0.935 mm, the association is refused EVERY step and the run reports:

   | jump limit | association moved per step, mm |
   |---|---|
   | 0.5 mm | 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 |
   | raised clear | 0.9233 0.9317 0.9382 0.9428 0.9453 0.9457 0.9442 0.9406 (slide 0.9350) |

   **"association moved 0.0000 mm" reads as SETTLED and means REFUSED.** Persistent association was
   not being exercised at all, and the frozen/re-linearised arms were bit-identical for that reason.
   The anatomical runs report 0.4999 mm against the same 0.5 mm limit, which is the same clamp one
   step from binding. Same shape as the four earlier catches: a reassuring number meaning its
   opposite.

   *Gate CC, with the limit raised so persistence is real.* CC-flat reads **max |u - d| 0.000002 mm,
   min J 1.000000** -- two nanometres, the same noise floor gate W reports at 0.000001 mm, and
   2.7e-7 of the 7.48 mm motion, at the solver's own rtol of 1e-7. **My declared criterion was
   "exactly zero", and 2 nm is not exactly zero, so the script prints FAIL.** I am not rewriting the
   criterion after seeing the number; I am reporting that "exactly zero" is stricter than a
   double-precision solve can be and leaving the judgement where it belongs.

   **FOUR CONSTRUCTIONS OF GATE CC WERE WRONG BEFORE ONE WAS RIGHT, ALL MINE.** Recorded because
   three of the four returned a plausible number rather than an error:
   1. a tangential drive on the contact face -- **nothing moved**, because u = 0 already satisfies
      n.u = 0. Every case returned max |u - d| = 7.480000 mm, min J exactly 1, volume ratio exactly
      1, in 0 s. The trap gate W exists to expose, and it is quoted in a comment I had read aloud
      earlier the same session.
   2. tilting the drive out of the tangent -- the block moved only along the normal, and |u - d| came
      back as **exactly the tangential component of d**. The tangential direction is undetermined.
   3. dragging the contact face through the block, but leaving `travel` at its default `-gap0`, so
      the face was ALSO being seated 1 mm. The plane read 1.402927 mm. Running it against the
      pre-repair code gave **bit-identical numbers**, which is what proved the reference wrong rather
      than the solver.
   4. `travel = 0` with the drive still tilted -- now the face is pinned normally while the far face
      pulls it up, so the block stretches. The plane read 2.602647 mm.

   A guard that says **"NOTHING MOVED -- the drive is not driving; this case is void"** is now in the
   script, and it fired on all five cases of construction 1. The lesson is the standing one in a new
   place: ask whether the mechanism can do anything at all, not only whether it can fail.

   ### BB′ fails into my bad branch: the contact is cycling at a fixed load

   Six zero-sized steps at one fixed fraction, each from the last one's answer, in mm:

   > **0.335152 → 0.122123 → 0.034554 → 0.036361 → 0.037410 → 0.036119**

   It decays for three repeats and then **plateaus at ~36 µm**. That is not one-time relaxation into
   a loosened feasible set; it is the branch `f3b7c08` named as the bad one — **at a fixed load the
   association and the solve are chasing each other.** The motion is almost entirely on the base
   (36.1 µm against 6.0 µm interior), so it is the contact cycling, and **the repair does not touch
   it.** BB's bounds half is closed (0.000e+00 mm over 0 nodes, against 7.19 mm over 26 nodes and 4
   inversions); BB′ is a separate, open defect.

   ### R2: the repair is good, not merely neutral

   | same stage, same data, only the code differs | outcome |
   |---|---|
   | pre-repair (`50e7583`, worktree) | **stalls at fraction 0.6250**, 1 inverted |
   | with the repair | **completes to 1.0**, 0 inversions, `min J` 0.4174, median \|u−d\| 0.256 mm |

   **This supersedes the 0.710-vs-0.920 comparison, which must not be quoted** — it was one fraction
   in two runs that took different paths. A stall-versus-completion on the identical stage is a real
   comparison. Regressions all pass: T-none (0.000 mm median and max, `min J` 1.0000), W, X,
   V-on-plane.

   ### A jump limit that silently disabled the thing it protects

   At `jump_limit = 0.5 mm` against a per-step slide of **0.935 mm**, the association is refused
   **every step**:

   | limit | association moved per step |
   |---|---|
   | 0.5 mm | **0.0000 mm ×8** |
   | raised clear | 0.9233 … 0.9406 (slide 0.9350) |

   **"association moved 0.0000 mm" reads as settled and means refused.** Persistent association was
   not being exercised at all — **which is why frozen and re-linearised came back bit-identical.**
   That comparison is void: it compared two runs that were both frozen. And the anatomical runs print
   **0.4999 mm against the same 0.5 mm limit** — one step from binding, which is not a margin.

   ### CC-flat FAILS as written, and the bar is not being moved

   CC-flat reads **0.000002 mm, `min J` 1.000000** — 2 nm, 2.7e-7 of the motion, against a declared
   criterion of *exactly zero*. **The script prints FAIL and that stands.** This programme records a
   bar below the noise floor as FAILED and does not rescore it: gate A did exactly this in `de88ef4`.
   My criterion was unachievable in floating point for any nontrivial computation, which is my error
   in writing it, not a reason to rewrite it now.

   **And the replacement is not a new bar — it is a convergence test**, which needs no chosen number:

   > **CC-flat′: re-run CC-flat at successively tighter solver tolerances.** If the deviation falls
   > with `rtol`, it is the convergence floor and the flat case is clean. If it plateaus, it is a
   > defect and the 2 nm is real.

   That is the same move as BB′ — turn a threshold question into a convergence question — and it is
   the only form that can distinguish "at the solver's floor" from "small but wrong". Noting against
   it: 2.7e-7 is about **3× the declared rtol of 1e-7**, so a principled bar at rtol would also have
   failed. The measurement is the honest instrument here, not any threshold.

   ### Four wrong constructions, three of which returned a plausible number

   Building CC went wrong four times: a tangential drive that moved nothing (`u = 0` already
   satisfies `n·u = 0` — **the exact trap gate W exists for, quoted in a comment the agent had read
   aloud earlier in the same session**); a tilted drive where `|u−d|` came back as exactly the
   tangential component; a default `travel = -gap0` that seated the face as well as dragging it; and
   `travel = 0` still tilted, which stretched the block. **Only one of the four produced an error;
   three produced numbers that looked like results.**

   The fix is a guard, not the knowledge: a check that prints **"NOTHING MOVED — the drive is not
   driving; this case is void"**, which fired on all five cases of the first construction. **Knowing
   a trap does not protect against it** — the agent had the tangential-drive trap in front of it and
   walked into it anyway. Only the assertion does.


   ### CC-flat′ passes: staleness costs exactly zero on a plane

   | rtol | max \|u − d\| as a fraction of \|d\| |
   |---|---:|
   | 1e-4 … 1e-8 | 2.511e-07, **identical at all five** |
   | 1e-9 | **1.241e-14** |
   | 1e-10 | **6.408e-15** |

   Tightening `rtol` removes the residual entirely, to double precision. **The 2 nm was the
   convergence floor, and the flat case is clean** — the cost of staleness on a surface whose normal
   never rotates is zero, exactly as the decomposition requires. CC-flat stays recorded as FAILED
   against its as-written criterion; CC-flat′ is what establishes the fact.

   **So the decomposition has done its job on the half that separates a bug from a cost.** The
   platens are clean, therefore R2's degradation under persistent association is **not** a defect in
   the association — it is the price of staleness on curved geometry, which CC-curved is measuring.

   **A statistic fitted to a shape it does not have**, caught by the agent against its own script.
   The first version fitted a power law and reported *"falls as rtol^1.335"*. That number is
   meaningless: the measured shape is a **plateau across five decades and then a collapse**, because
   the loose tolerances all stop at the same Newton iterate. A slope through a plateau-plus-cliff is
   a number with no referent. The slope is gone and the verdict now turns on what actually decides
   it — whether the deviation reaches machine precision at all. It does, by seven orders.

   **And it vindicates running the sweep rather than reasoning about it.** 2.7e-7 is ~3× the declared
   rtol, so a principled bar at rtol would also have failed; the plateau is why. No threshold could
   have separated "at the floor" from "small but wrong" here — only varying the tolerance could.

   **CC-curved, interim and not to be read as a result:** staleness **0.001856 mm at R = 400 mm** and
   **0.004481 mm at R = 200 mm** — a factor **2.41** for a halving of R against a predicted 2, and
   sitting at **0.21–0.26 of the `d²/2RN` prediction**, so the sagitta behaves as an upper bound with
   the solve compensating part of it rather than as an equality. **Two points are not an exponent**;
   R = 100 and 50 mm are still running.

   **CC-FLAT' PASSES, CC-CURVED CONFIRMS THE SCALING AND FAILS MY BAND, AND THE BAND WAS THE PART
   THAT WAS WRONG.**

   *CC-flat', the convergence test rather than a new number.* The identical flat case at successively
   tighter solver tolerances, as a fraction of the 7.48 mm drive:

   | rtol | 1e-4 | 1e-5 | 1e-6 | 1e-7 | 1e-8 | 1e-9 | 1e-10 |
   |---|---:|---:|---:|---:|---:|---:|---:|
   | max \|u - d\| / \|d\| | 2.511e-07 | 2.511e-07 | 2.511e-07 | 2.511e-07 | 2.511e-07 | **1.241e-14** | **6.408e-15** |

   Tightening rtol removes the residual entirely, to double precision -- a fall of 3.9e7. **The 2 nm
   was the convergence floor, not a defect, and the flat case is clean: the cost of staleness on a
   plane is zero, as the decomposition requires.** CC-flat itself stays FAILED under its declared
   "exactly zero" and is not rescored.

   A NUMBER OF MINE WITHDRAWN INSIDE THE ROUND: the first version of this sweep fitted a power law
   and printed "falls as rtol^1.335". The measured shape is a PLATEAU across five decades and then a
   collapse -- the loose tolerances all stop at the same Newton iterate -- so a fitted exponent is a
   number without a meaning. It is removed from the script; the verdict now turns on whether the
   deviation reaches machine precision at all.

   *CC-curved, the full sweep.* Frozen frames against re-linearised ones, the difference being the
   staleness and nothing else:

   | R | staleness, mm | predicted d^2/2RN, mm | ratio |
   |---:|---:|---:|---:|
   | 400 mm | 0.001856 | 0.008742 | 0.212 |
   | 200 mm | 0.004481 | 0.017484 | 0.256 |
   | 100 mm | 0.007927 | 0.034969 | 0.227 |
   | 50 mm | 0.018077 | 0.069938 | 0.258 |

   **The staleness scales as (1/R)^1.068 against the sagitta argument's 1.000.** The exponent is the
   part of the prediction that has a derivation, and it matches. The prefactor is a constant 0.24 --
   the measured cost is about four times SMALLER than the sagitta bound, consistently at every
   radius, which is what an upper bound behaves like when the solve carries part of the offset and
   the per-step errors partly cancel rather than accumulating as N times one of them.

   **MY DECLARED BAND WAS ratio IN [0.5, 2.0], SO THE GATE PRINTS FAIL, AND I AM NOT REWRITING IT.**
   The band was mis-designed rather than mis-set: it folded two claims into one test. The derivation
   supports the EXPONENT; nothing in it fixes the CONSTANT to within a factor of four, and I had no
   basis for predicting one. A future run should test the two separately. This round's verdict stands
   as printed, and the finding stands as stated: **on curved geometry the cost of persistent
   association is staleness, it tracks 1/R, and the sagitta is an upper bound on it by about 4x.**
   Combined with CC-flat' reading zero on a plane, the decomposition is complete: the degradation is
   the price of curvature, not a bug.


   ### CC-curved: the exponent confirms, the prefactor was never derivable, and the gate fails

   | R | staleness | predicted `d²/2RN` | ratio |
   |---:|---:|---:|---:|
   | 400 mm | 0.001856 mm | 0.008742 mm | 0.212 |
   | 200 mm | 0.004481 mm | 0.017484 mm | 0.256 |
   | 100 mm | 0.007927 mm | 0.034969 mm | 0.227 |
   | 50 mm | 0.018077 mm | 0.069938 mm | 0.258 |

   **Staleness scales as `(1/R)^1.068` against the sagitta argument's 1.000.** The prefactor is a
   constant **0.24** — four times *smaller* than the bound, at every radius, which is how an upper
   bound behaves when the solve carries part of the offset and per-step errors partly cancel rather
   than accumulating as N copies of one.

   **The declared band was ratio ∈ [0.5, 2.0], so the gate prints FAIL, and it is not rewritten.**

   **But the band was mis-DESIGNED, not mis-set, and that distinction is the lesson.** It folded two
   claims into one number. The sagitta derivation fixes the **exponent** — and the exponent matched,
   1.068 against 1.000. **Nothing in that derivation fixes the constant**, and there was no basis for
   predicting one to within a factor of four. A single ratio test cannot pass the half that was
   derived and fail the half that was guessed; it can only report the guess.

   > **Gate only the part of a prediction that has a derivation.** Where a scaling argument gives an
   > exponent and no prefactor, the gate is on the exponent; the prefactor is *measured and
   > reported*, not predicted. Recorded in `CLAUDE.md`.

   **The decomposition closes, and it answers what CC was built for.** CC-flat′ puts the cost of
   staleness at **zero on a plane** (2 nm falling to 6.4e-15 as `rtol` tightens, a fall of 3.9e7).
   CC-curved puts it at **0.24 × d²/2RN, tracking 1/R**, on curved geometry. **So R2's degradation
   under persistent association is the price of curvature, not a defect in the association** — which
   is the bug-versus-cost question the decomposition existed to settle.

   **Affordability, interim and moving:** fraction **0.2324** after 2,977 s, `min J` 0.406, having
   earlier accepted 0.2852 before a cut-back — roughly **4× the reach** of the 0.0625 where it
   stalled before the repair. Not quoted as final, and when it is, it must say what it was
   exercising: an association **refused** at 0.4999 mm against a 0.5 mm jump limit, not a converged
   one.

   ### Acceptance does not test convergence, and that gap reaches past the affordability number

   The drive reached fraction **0.3115** and stopped on the **wall-clock cap, not a stall** — but the
   last three accepted steps each ran **`Newton 300`, which is `max_newton`.** They exhausted the
   iteration budget rather than converging, and **`settled` checks only association motion and
   penetration**. So a step can be *accepted while unconverged*, and part of the extra reach may have
   been bought that way.

   **0.3115 is therefore not quotable.** The honest statement is that the repair took the drive past
   the 0.0625 where it stalled, and it was still advancing when the clock ran out. `UNCONVERGED` is
   now printed on the step line and the run relaunched.

   **The gap is in the acceptance rule, so it is not confined to the affordability number.** Every
   verdict on this line that rests on *accepted steps* inherits the question, and they do not all
   inherit it equally:

   * **T-none, W, X, V-on-plane are safe by inspection** — they report **one Newton iteration per
     step**, which cannot be an exhausted budget.
   * **R′ completing to 1.0 and the anatomical runs are not**, and R′-completes-where-pre-repair-
     stalls is the headline evidence that the repair is good.
   * **Pre-registered:** re-check every recorded verdict for accepted-but-unconverged steps and
     report which verdicts rest on any. A verdict with none stands as recorded; a verdict with some
     is **re-run**, not re-read.
   * **Predicted: R′ is clean and the anatomical runs are not.** R′ has no bed, so its steps are
     easy; the anatomical bed is where the iteration budget would be exhausted. Stated with the usual
     caveat about my record here.

   **This is the same shape as the jump limit printing `0.0000 mm` for "refused".** An acceptance
   criterion that omits the solve's own convergence reports success for a step that did not finish —
   a green result from a check that never looked at the thing it is named for.

   **THE ACCEPTED-BUT-UNCONVERGED AUDIT: EVERY RECORDED VERDICT SURVIVES, AND THE PROXY IS A PROOF.**

   `settled` checks the association's motion and the unilateral penetration and NEVER looks at the
   solve's own convergence, so a step can be accepted while the Newton loop merely ran out of budget.
   Every verdict on this line that rests on accepted steps inherits the question, so every one was
   re-checked rather than re-read.

   **First, the iteration count was established as a PROOF and not a proxy.** The Newton loop has
   exactly one `break` (`res <= tol`, convergence), one nested `break` that leaves only the Armijo
   line search, and one `raise` (line-search failure, which propagates and triggers a cut-back rather
   than returning a result). There is no other way out. **Therefore `iterations < max_newton` implies
   the solve converged**, and a step at the 300 ceiling is the only kind that can be accepted
   unconverged.

   | verdict | solves | max Newton | at the 300 ceiling |
   |---|---:|---:|---:|
   | T-none | 8 | 1 | **0** |
   | W | 8 | 1 | **0** |
   | X | 8 | 1 | **0** |
   | V on the plane | 8 | 1 | **0** |
   | **R2 / R-prime, repaired, completes to 1.0** | 8 | **148** | **0** |
   | pre-repair R2, the stall at 0.6250 | 5 | 142 | **0** |
   | gate AA | 1 | 30 | **0** |
   | gate BB-prime | 7 | 30 | **0** |

   **Every recorded verdict stands as recorded. None rests on an accepted-but-unconverged step.**
   R-prime completing to 1.0 where the pre-repair code stalls is clean on BOTH sides -- the stall is
   a real stall and the completion is a real completion -- so the headline evidence that the repair
   is good is undamaged. The prediction that R-prime would be clean holds, though not for the reason
   given: its steps are not easy, running to 148 iterations. Only the anatomical seating drive
   reached the ceiling, which is the half of the prediction that was right.

   **Fraction 0.3115 is withdrawn from the record.** The honest statement is that the repair took the
   drive past the 0.0625 where it stalled, and it was still advancing when the clock ran out.

   **A PAIR WORTH NAMING TOGETHER.** `0.0000 mm` printed for REFUSED, and `accepted` printed for RAN
   OUT OF ITERATIONS. Both are green results from a check that never looked at the thing it is named
   for. The first voided a comparison; the second could have voided a verdict and did not. The step
   line now prints `UNCONVERGED`, and gate CC now reports its maximum Newton count and its count of
   unconverged accepted steps, because a gate that cannot report this cannot be audited later.


   ### The convergence audit: every verdict survives, and the check was made a proof first

   **Before trusting the iteration count, its meaning was established by enumerating the loop's
   exits.** The Newton loop has exactly one `break` (`res <= tol`), one nested `break` that leaves
   only the Armijo line search, and one `raise` (line-search failure, which propagates into a
   cut-back rather than returning a result). There is no other exit. **Therefore
   `iterations < max_newton` implies the solve converged**, and a step at the 300 ceiling is the only
   kind that can be accepted unconverged.

   That is the difference between a proxy and a proof, and it is why the audit below is worth
   anything: the iteration count is not *evidence about* convergence here, it is *equivalent to* it.

   | verdict | solves | max Newton | at the ceiling |
   |---|---:|---:|---:|
   | T-none / W / X / V-on-plane | 8 each | 1 | **0** |
   | **R2 / R′ repaired, completes to 1.0** | 8 | **148** | **0** |
   | pre-repair R2, the stall at 0.6250 | 5 | 142 | **0** |
   | AA | 1 | 30 | **0** |
   | BB′ | 7 | 30 | **0** |

   **Every recorded verdict survives. R′ is clean on both sides** — the stall is a real stall and the
   completion is a real completion — so the headline evidence that the repair is good is undamaged.

   **My prediction was right about which half was at risk and wrong about why.** I said R′ would be
   safe because "R′ has no bed, so its steps are easy". Its steps are **not** easy: they run to 148
   iterations. It is clean because it **converged**, not because it had little to do. Right
   conclusion, wrong mechanism — which is the third time this week I have reached a correct verdict
   through an argument that did not hold, and the reason to report mechanisms rather than just
   calls.

   **Both gaps are closed at the source rather than worked around.** The step line prints
   `UNCONVERGED`, and gate CC now reports its maximum Newton count and its count of unconverged
   accepted steps. **CC is re-running under that reporting**, because its verdict was recorded with
   `log=None` and carried no counts at all — the agent refused to leave CC as the one gate exempt
   from its own audit, which is the right instinct and not one I asked for.

   **A candid note on what the audit actually relied on.** This verdict survived only because
   convergence could be *reconstructed* from the iteration count after the fact. That was a property
   of this particular loop, not instrumentation — if the loop had had a second silent exit, the
   reconstruction would have been impossible and the verdicts unrecoverable. The lesson is the
   printed flag, not the successful reconstruction.

   ### The affordability number with the convergence flag on: 0.1094, not 0.3115

   The re-run prints `UNCONVERGED` on the step line, and it changes the answer:

   | accepted at fraction | Newton | held gap to the aim | `min J` |
   |---:|---|---:|---:|
   | 0.0625 | 30 | 14.80 mm | 0.800 |
   | **0.1094** | **77** | 14.78 mm | 0.667 |
   | 0.1797 | **300 UNCONVERGED** | 13.72 mm | 0.509 |
   | 0.2324 | **300 UNCONVERGED** | 14.43 mm | 0.406 |
   | 0.3115 | **300 UNCONVERGED** | 13.22 mm | 0.397 |

   **Three of the five accepted steps exhausted the iteration budget, and every step past 0.1094 is
   one of them.** So the honest reach is **fraction 0.1094**, the last converged step — not the
   0.3115 that was already withdrawn, and not the 0.2852 quoted as interim before that. Against the
   0.0625 where the drive stalled before the repair, the repair is worth **about 1.75×**, not 4×
   and not 5×.

   **`min J` falls monotonically across the unconverged steps** — 0.667 → 0.509 → 0.406 → 0.397 —
   which is what accepting an unfinished solve looks like: each one starts further from a converged
   state than the last.

   **A second thing in the same table that nobody asked about: the held gap to the aim is not
   closing.** 14.80 → 14.78 → 13.72 → 14.43 → 13.22 mm while the load fraction nearly quintuples.
   **Whether that is a lag behind a moving aim or a failure to approach a fixed one decides whether
   the fraction measures progress at all**, and it is not answerable from the log as written. If the
   aim advances with the fraction then a roughly constant gap is a constant lag; if the aim is the
   final target then the drive is barely approaching it. Asked of the agent rather than assumed,
   because a reach number that does not correspond to getting closer is another quantity that does
   not measure what its name says.

   **What stands unchanged:** the repair itself (bounds at a zero step exactly 0.000e+00 over 0
   nodes), R′ completing where pre-repair stalls with **zero** steps at the ceiling on either side,
   and CC's decomposition. Those were audited for exactly this and came back clean. Only the reach
   number moves.

   **`aim` IS THE FINAL TARGET, SO THE LOAD FRACTION MEASURES LOAD APPLIED AND NOT DISTANCE CLOSED --
   AND BY THE LOGGED MEASURE THE DRIVE HAS CLOSED NOTHING.**

   Answering it from the code rather than the log: `held_err = |gap[held] - (gap0[held] + aim*step)|`
   and `aim` is the PHASE TARGET, passed as `advance(trial, u, theta, target)`. It is 1.0 for the
   forward phase. I changed it to that in the repair, because the increment now closes theta of what
   remains and a nominal mid-drive schedule is a number the drive is not trying to hit. So the logged
   "held gap to the aim" is the distance to FULL SEATING, and it is the coordinator's second branch.

   | fraction | 0 (derived) | 0.0625 | 0.1094 | 0.1797 | 0.2324 | 0.3115 |
   |---|---:|---:|---:|---:|---:|---:|
   | distance to full seating, max over held nodes | **14.18 mm** | 14.80 | 14.78 | 13.72 | 14.43 | 13.22 |
   | Newton | -- | 30 | 77 | **300 UNCONVERGED** | **300 UNCONVERGED** | **300 UNCONVERGED** |

   **At the last CONVERGED step the worst held node is 14.78 mm from full seating, against 14.18 mm
   before the drive started. It has gone backwards.** The load fraction quintupled and the distance
   closed is negative on that measure. A reach number that does not correspond to getting closer is
   the fifth quantity on this line that does not measure what its name says.

   *The reach, stated honestly.* **Fraction 0.1094**, the last converged step, against the 0.0625
   where the pre-repair drive stalled: about **1.75x**, not the 4x I wrote or the 5x the raw fraction
   suggests. The three steps past it exhausted the 300-iteration budget, confirmed by the new
   `UNCONVERGED` flag rather than inferred from the count. min J falling monotonically across them --
   0.667, 0.509, 0.406, 0.397 -- is what accepting an unfinished solve looks like: each starts
   further from converged than the last.

   *What the max cannot tell me, and what is now logged.* The step line reported only the maximum, so
   there was no way to know whether one node is stuck or the whole sheet is. It now reports the
   MEDIAN alongside it. Until that run returns, the distance closed is reported as unknown for the
   sheet and negative for the worst node, and the fraction is quoted only with both.

   *Also: what the anatomical runs were exercising.* Their association was REFUSED, pinned at
   0.4999-0.5000 mm against the 0.5 mm jump limit, not converged. Every number from them carries
   that.

   *Stale logs marked.* `gate_aa.log`, `gate_aa2.log` and `gate_aa3.log` in the subject directory
   reported 42 and 100 inverted elements from the pre-seed-fix runs, with a "THE PRICE: not
   measurable" line that read as current. Each now carries a SUPERSEDED header naming the defect, the
   replacement numbers and the commits. A stale log with a confident conclusion in it is what this
   line has been bitten by four times.


   ### `aim` is the final target, and the distance to it went BACKWARDS

   `held_err = |gap[held] − (gap0[held] + aim*step)|` with `aim` the **phase target**, 1.0 for the
   forward phase. So the logged "held gap to the aim" is the distance to **full seating**, not to a
   moving mid-drive waypoint. The fraction measures **load applied**, not distance closed.

   | fraction | 0 (derived) | 0.0625 | **0.1094** | 0.1797 | 0.2324 | 0.3115 |
   |---|---:|---:|---:|---:|---:|---:|
   | distance to full seating, **max** | **14.18 mm** | 14.80 | **14.78** | 13.72 | 14.43 | 13.22 |
   | Newton | — | 30 | 77 | 300 **U** | 300 **U** | 300 **U** |

   **At the last converged step the worst held node is 14.78 mm from full seating, against 14.18 mm
   before the drive started.** The load fraction quintupled and, on this measure, the distance closed
   is **negative**.

   **The scope is not established and must not be overstated.** The step line reported only the
   maximum, so whether one node is stuck or the whole sheet is cannot be told from it — a third
   instance on this line of a summary that cannot answer the question its name implies. The median is
   now reported alongside the max and the drive is re-running; **from here the fraction is quoted
   only with both.**

   **What this does to the framing.** "Affordability" asked how far the drive could get. The prior
   question — whether it was getting anywhere — had not been asked, and on the only measure of it
   available the answer so far is no. The reach of **fraction 0.1094**, about **1.75×** the 0.0625
   stall, stands as a statement about *load carried*, and carries no implication about seating.

   **What still stands:** the repair (bounds at a zero step exactly 0.000e+00 over 0 nodes), R′
   completing where pre-repair stalls with zero ceiling steps either side, CC's decomposition, and
   the three stale `gate_aa*.log` files now carrying SUPERSEDED headers naming the seeding defect.
   The anatomical runs, whenever quoted, were exercising an association **refused** at 0.4999–0.5000
   mm against a 0.5 mm jump limit.

   **The lesson the agent drew against itself is the transferable one.** The convergence audit
   checked *verdicts against the runs that produced them*; it did not check *whether the headline
   metric meant what its name said*. **Those are two different audits**, and passing the first says
   nothing about the second — which is how a number that had gone backwards survived a clean audit
   and three separate quotations.

   **IT IS 39 NODES OF 3,123, NOT THE SHEET. THE DRIVE IS SEATING, AT THE MEDIAN, AT THE SCHEDULED
   RATE.**

   The step line now reports the median beside the maximum, the median distance CLOSED, and the count
   of held nodes further from the aim than they started. Each held node begins exactly |travel| away,
   so that count is what separates one stuck node from a stuck sheet -- which neither a maximum nor a
   median can do alone, since a flat median is equally consistent with nothing moving and with half
   moving each way. At the first accepted step:

   | fraction 0.0625 | value |
   |---|---:|
   | distance to full seating, max | 14.8039 mm |
   | distance to full seating, median | **7.0196 mm** |
   | distance CLOSED, median | **+0.4540 mm** |
   | held nodes further away than they started | **39 of 3,123 (1.25%)** |

   Median travel is 7.484 mm, so the nominal closure at fraction 0.0625 is 0.468 mm and the measured
   is 0.454 mm -- **97% of the scheduled rate**. The alarming maximum was one of about 39 outliers.

   **SO "THE DRIVE GOES BACKWARDS" IS REFUTED, AND THE CAUTION THAT PREVENTED IT FROM BEING RECORDED
   WAS WORTH MORE THAN THE OBSERVATION THAT PROMPTED IT.** The scope was declared unestablished
   rather than asserted, and the scope is what the answer turned on. What remains true and unchanged:
   the load fraction still measures load applied rather than distance closed, the reach still stands
   at fraction 0.1094 as a statement about load carried, and about 39 held nodes ARE moving away from
   their target and are now a named, bounded problem instead of an unbounded one.


   ### Scope: 39 nodes of 3,123. "The drive goes backwards" is REFUTED.

   | at fraction 0.0625 | |
   |---|---:|
   | distance to full seating, max | 14.8039 mm |
   | distance to full seating, **median** | **7.0196 mm** |
   | distance **closed**, median | **+0.4540 mm** |
   | held nodes further away than they started | **39 of 3,123 (1.25%)** |

   Median travel is 7.484 mm, so the scheduled closure at fraction 0.0625 is 0.468 mm against
   **0.454 mm measured — 97% of the scheduled rate.** The drive is closing the gap at very nearly
   the rate the schedule asks for. **The alarming maximum was one of about 39 outliers.**

   **The caution was worth more than the observation that prompted it.** The previous entry recorded
   the max going backwards and said in terms that the scope was unestablished and must not be
   overstated — and the scope is exactly what the answer turned on. Had "the drive goes backwards"
   gone in as a finding it would have been the sixth withdrawn claim on this line.

   **And it would have been withdrawn from the harder direction.** The other five came from numbers
   that flattered, which invite re-examination. This one was alarming, and an alarming number is the
   one least likely to be checked again — it already agrees with the pessimism a hard week
   generates. Worth naming as its own hazard: **a bad number gets less scrutiny than a good one, and
   for the wrong reason.**

   **What survives unchanged:**

   * the load fraction measures **load applied**, not distance closed — that finding stands on its
     own and was not what the scope question touched;
   * the reach stands at **fraction 0.1094**, as a statement about load carried;
   * **about 39 held nodes are genuinely moving away from their target.** That is now a named,
     bounded problem rather than an unbounded one, and **1.25% sits close enough to this line's own
     1% gates that it should be sized properly rather than waved through.**

   **Next, and it separates the two behaviours:** whether the median keeps closing at ~97% of
   schedule at 0.1094 and beyond while those 39 diverge. If it does, the outliers can be
   characterised on their own without contaminating the closure measurement.

   ### One jump limit, three ways of being invisible

   The 0.5 mm association jump limit has now produced three distinct failures of observability, and
   they are worth listing together because each was found separately and read as its own problem:

   1. **It disabled what it protects, silently.** Against a 0.935 mm per-step slide it refused the
      update *every step*, so persistent association was never exercised — which is why frozen and
      re-linearised runs came back bit-identical. That comparison was void.
   2. **It reported refusal as settlement.** `association moved 0.0000 mm` reads as *converged* and
      means *refused*. The anatomical runs print 0.4999–0.5000 mm against the same limit, one step
      from binding.
   3. **It hid which nodes it acted on.** `association moved` is the maximum over **accepted**
      updates, so a refused node contributes nothing to it and vanishes from the log entirely. The
      refused set was **unobservable from outside the solver** until a per-node mask was added to ask
      the current question. A quantity reported only over the accepted population cannot describe the
      rejected one, and nothing in the log said so.

   **The guard built into the characterisation, and its lineage.** The question "are the diverging
   nodes the refused ones?" is an overlap, and an overlap between two sets drawn from 3,123 has a
   chance level: with *b* diverging and *R* refused, agreement by coincidence is `b·R/3123`. The
   instrument reports the overlap **beside that chance level**, and its verdict line says "the
   refusal is implicated" only above **twice chance** — declared in the file before the numbers,
   not chosen after.

   That is `docs/LOG.md` row 5 in a new domain: a cross-modal association read 4.4× baseline because
   the metric averaged over 10⁹ *pairs* when only 187,484 were *edges*, diluted 5000:1. Same shape as
   today's error message reporting 5,660 nodes when only 4,097 are bounded — which was caught only
   because that number was arithmetically impossible. **An overlap without its chance level, and a
   count without its denominator, are the same mistake**, and the second one is caught by luck
   whenever the arithmetic happens not to be absurd.

   **THE 39, SIZED. THEY SHRINK, THEY ARE SCATTERED, AND THEY ARE NOT THE REFUSED NODES -- THE CAUSE
   I PROPOSED IS REFUTED, BELOW CHANCE.**

   | fraction | 0.0625 | 0.1094 | 0.1797 |
   |---|---:|---:|---:|
   | held nodes moving away, of 3,123 | **39** | **34** | 27 |
   | worst, backwards | 1.445 mm | 4.214 mm | 3.677 mm |
   | their own asked travel, mean | 6.951 mm | 6.468 mm | 5.324 mm |
   | spread on the base | 103.6 mm | 103.9 mm | 119.8 mm |
   | same nodes as the previous fraction | -- | **34 of 39 (87%)** | 22 of 34 (65%) |
   | association REFUSED this step, of 3,123 | 963 | 982 | 1,144 |
   | diverging nodes among the refused | **0** | **1** | 3 |
   | expected by chance | **12.0** | **10.7** | 9.9 |

   Convergence, from the identical drive logged with the flag: 0.0625 (Newton 30) and 0.1094 (77) are
   CONVERGED; 0.1797 exhausted the 300 budget and is shown for shape only.

   **WHERE:** scattered, not clustered. They span 103.6 mm on a base whose whole extent is of that
   order, so this is not one bad patch of bed.

   **THE SAME?** largely yes -- 87% carry over between the two converged fractions, so there is a
   persistent core rather than a rotating cast.

   **HOW MANY:** the set SHRINKS, 39 -> 34 -> 27. Whatever they are, the drive is not recruiting more
   of them as it proceeds.

   **REFUSED? NO, AND EMPHATICALLY SO.** I proposed the jump-limit refusal as their cause because it
   would have connected them to something already on the record. The overlap is 0 of 39 and 1 of 34
   against 12.0 and 10.7 expected at random -- **far BELOW chance, not above it**. The diverging nodes
   are drawn almost exclusively from those whose association updates are being ACCEPTED. The refusal
   is not their cause, and a tidy connection to an existing finding is exactly the kind I should want
   refuted rather than confirmed.

   The chance column is why this reads as a refutation rather than a weak confirmation. With 963 of
   3,123 nodes refused, roughly a third of the base, ANY set of 39 would overlap about 12 by
   construction. An overlap quoted without its chance level would have been the same error as an
   inversion count without a denominator -- the one that produced "5,660 nodes" against a base of
   4,097 earlier today, and which was caught then only because that particular number happened to be
   arithmetically impossible. Here it would not have been.

   **A NUMBER THAT WAS INVISIBLE UNTIL NOW:** 963 to 1,144 of 3,123 held nodes -- **31% to 37% --
   have their association update REFUSED at every step.** The log reported only "association moved
   0.4999 mm", which is the maximum over ACCEPTED updates, so a third of the base was being frozen
   every step with no trace in the output. This is the jump limit's third distinct form of blindness:
   it disabled persistence entirely in gate CC, it reads as "settled" rather than "refused", and it
   hid which and how many nodes it was acting on.


   ### The 39 characterised: persistent, scattered, shrinking — and NOT the refused ones

   | fraction | 0.0625 | 0.1094 | 0.1797 (unconverged) |
   |---|---:|---:|---:|
   | moving away, of 3,123 | **39** | **34** | 27 |
   | worst, backwards | 1.445 mm | 4.214 mm | 3.677 mm |
   | spread on the base | 103.6 mm | 103.9 mm | 119.8 mm |
   | same nodes as previous | — | **34 of 39 (87%)** | 22 of 34 |
   | association **refused** this step | 963 | 982 | 1,144 |
   | diverging nodes among the refused | **0** | **1** | 3 |
   | **expected by chance** | **12.0** | **10.7** | 9.9 |

   Scattered across 103.6 mm rather than one bad patch of bed; **87% persist** between the two
   converged fractions; and the count **shrinks**, 39 → 34 → 27. The drive is not recruiting more.

   **My proposed cause is refuted, and below chance.** I suggested the jump-limit refusal because it
   would have tied these nodes to something already on the record — which is precisely the kind of
   tidy connection that deserves suspicion. The overlap is **0 of 39 and 1 of 34 against 12.0 and
   10.7 expected at random**: not merely unsupported but *anti*-correlated. The diverging nodes come
   almost exclusively from those whose association updates are being **accepted**.

   **The chance column is what made that readable**, and it is the reason the guard was built before
   the numbers. With 963 of 3,123 refused — a third of the base — *any* set of 39 overlaps about 12
   by construction. The raw overlap would have been the denominator error again, and unlike this
   morning's "5,660 of 4,097" it would have looked entirely plausible.

   ### The larger finding: a third of the contact constraints are frozen every step

   **963–1,144 of 3,123 held nodes — 31% to 37% — have their association refused at every single
   step.** This was invisible because the log printed only `association moved 0.4999 mm`, the maximum
   over *accepted* updates. **Every anatomical seating result to date was computed with roughly a
   third of the contact constraints held stale by a limit nobody was reading.** That is a live
   setting on the production drive, not a diagnostic artefact.

   **And it leaves the line in a genuine bind rather than with an obvious fix:**

   * the per-node 0.5 mm limit freezes **31–37% of nodes** every step;
   * the all-or-nothing rule with the derived threshold refuses **93.8% of whole steps** (gate U);
   * gate T showed the pure extremes fail in both directions — `none` cannot follow a sliding bed,
     `all` imports the 46.6 mm teleports.

   **Pre-registered, and it needs no curvature model:** the refused association displacement is
   *directly* the staleness — it is exactly how far each association would have moved and did not.
   **Accumulate it per node across the drive and report the distribution.** CC's `0.24 × d²/2RN`
   gives the price of staleness on curved geometry, but it needs a local radius; the accumulated
   refused displacement measures the same thing in millimetres with nothing assumed.

   * **If the accumulated staleness is small** — a few tenths of a millimetre over the drive — the
     31–37% is a bookkeeping fact and the seating results survive with a caveat.
   * **If it is millimetres**, every anatomical result on this line was computed against constraints
     pointing somewhere the bed no longer is, and they are re-run rather than re-read.
   * **No prediction.** The refused nodes are refused *because* their associations want to move
     furthest, so they are the worst-case population by construction — which is an argument that this
     could be large, not a forecast that it is.

   ### The ratchet: a prediction with a mechanism, recorded before the distribution lands

   The agent's argument, which cuts against the "worst-case population" framing by making it
   sharper: **a node is refused when its association wants to move more than 0.5 mm in one step — and
   the refusal freezes the association, so the next step measures the wanted displacement from the
   same frozen point.** If the node keeps sliding, the wanted displacement keeps growing, and it is
   refused again. **Permanently.**

   That predicts a **heavy tail** rather than a modest median: a small set of sliding nodes whose
   staleness grows without bound, not a population all slightly behind. It is a mechanism rather than
   a guess, and it is testable from `refused_steps`, which records how many times each node was
   refused — a node refused at every step is distinguishable from one refused once.

   **Weak evidence already on the record, recorded now so it is not claimed as confirmation later:**
   the refusal counts run **963 → 982 → 1,144** across fractions 0.0625, 0.1094 and 0.1797. The set
   is growing. But the first two are only 2% apart and the third is an **unconverged** step, so this
   is consistent with the ratchet without being evidence for it. The distribution decides.

   **Three figures, not one**, because any single one misrepresents the others: the median over *all*
   held nodes (flattered by the ~65% never refused), the median over the *ever-refused* nodes (what
   the affected population suffers), and the maximum (one node).

   **The boundary is 1 mm, declared in the file before the run** — where the two branches meet, and
   the bed's own facet scale, the length below which this surface cannot distinguish anything. A
   staleness under it cannot point a constraint at a different feature; one over it can. The
   derivation is the reason for the number, not a justification found afterwards.

   **Two accounting decisions, stated so they can be disputed rather than buried:**

   * **Refusals inside a step that was later cut back are not counted.** The accumulator is saved and
     restored alongside the association, so only refusals on accepted steps contribute. Without it,
     every discarded attempt in the cut-back loop would inflate the total — and with four cut-backs
     before the first accepted step, that inflation would have been large and invisible.
   * **Nodes whose ray finds no bed at all are counted as a count, never as a distance.** "No bed to
     move to" is not a displacement, and folding it in as one would put a fabricated number into the
     very distribution the decision rests on.

   **If the ratchet holds, the earlier statement understates the problem.** "31–37% of constraints are
   stale each step" describes a fixed lag. A ratchet means a *growing* set of *permanently* frozen
   constraints, which is a different and worse fact, and it would make re-running the anatomical
   results the only option rather than one of two.

   ### The pre-registered measurement confirmed the pre-registered hypothesis — through its own artefact

   | | 0.0625 | 0.1094 |
   |---|---:|---:|
   | refused ≥ once | 963 (30.8%), max 1 time | 982 (31.4%), max **3** times |
   | median over affected | 0.8174 mm | **3.0842 mm** |
   | max | 47.8232 mm | **143.1322 mm** |

   On its face: the "millimetres" branch, catastrophically. **But 143.1322 ≈ 47.82 × 3** — the same
   ~47 mm teleport, refused three times, summed as though it were accumulated drift.

   **The accumulator cannot distinguish an association falling slightly further behind each step from
   one large jump rejected repeatedly**, and only the first is staleness. A 46.6 mm teleport is
   precisely what the jump limit exists to reject. **Counting a successful rejection as accumulated
   error credits the mechanism's correct behaviour to its failure mode.**

   **This is the day's sharpest failure, and none of the standing rules would have caught it.** The
   measurement was pre-registered with its boundary, its three figures and its accounting decisions
   declared in the file before the run. The hypothesis it appeared to confirm — the ratchet — was
   also pre-registered, with a mechanism, one message earlier. **Everything was in the right order and
   the answer was still an artefact of the instrument's own summation.** Pre-registration protects
   against choosing a threshold to fit a result; it does nothing about a quantity that is
   ill-defined for the question.

   **It was caught by a ratio being too clean**: 143.13 / 47.82 = 2.99. The median alone — 0.82 →
   3.08 across two steps — looked like an ordinary ratchet and would have passed. Had the numbers
   been reported, they would have triggered a re-run of every anatomical result on this line, and the
   re-runs would have been chasing a figure that partly measures the jump limit working as intended.

   > **The rule, recorded in `CLAUDE.md`: a sum over repeated events cannot tell "accumulating" from
   > "the same event, repeated". And confirmation is the moment to check the instrument, not the
   > moment to stop** — a result that agrees with the hypothesis you pre-registered is the one you
   > are least likely to take apart.

   **The repair, and it keeps both quantities.** The verdict moves to the **instantaneous** wanted
   displacement — how far each refused constraint is from the bed *right now* — which is what the
   question actually asks. The accumulated figures stay in the output, because **the ratio between
   them is itself the diagnostic**: accumulation ≫ instantaneous means the same teleport rejected
   repeatedly; accumulation ≈ instantaneous × refusal count means genuine drift.

   **The ratchet is neither confirmed nor refuted.** It may still hold — the instantaneous numbers
   will say — but the confirmation that appeared to arrive was produced by the summation, not by the
   body.

   ## THE SLIDING-SEATING LINE, STOPPED AND LEFT RESUMABLE

   The owner stopped the association-refusal investigation. The instantaneous-displacement
   measurement was never taken. This section is what someone picking the line up in a month needs.

   **THE ONE SENTENCE THAT MATTERS MOST: the accumulated refused-staleness figure is an ARTEFACT and
   must not be quoted.** It read median 3.0842 mm over affected nodes and max 143.1322 mm at fraction
   0.1094, and 143.13 is 47.7 x 3 -- the same ~47 mm teleport refused three times, summed as though
   it were accumulated drift. A teleport of that size is exactly what the jump limit exists to
   reject, so the figure partly credits the mechanism's CORRECT behaviour to its failure mode.
   Quoting it would trigger a re-run of every anatomical result on this line, chasing a number that
   partly measures the limit working.

   **ESTABLISHED, and none of it depends on the refusal question.**
   * *The repair.* A held node's bound was `on_plane + gap0 + fraction*travel` and the fraction scaled
     `travel` and not `on_plane`, so a drifted association demanded its whole drift instantly however
     small the step -- 7.19 mm over 26 nodes at a step of size ZERO. Every term is now under the
     increment parameter via `theta`. A zero-sized step displaces the bounds by exactly 0.000e+00 mm
     over 0 nodes. That one fact explains gates R, R', R'', S, T, U, V and AA: all were adjusting the
     term that was already small.
   * *The repair is good, not merely neutral.* R' completes to fraction 1.0 where the pre-repair code
     STALLS at 0.6250, same stage, same data, zero ceiling steps on either side.
   * *The solver is sound.* T-none, W, X and V-on-the-plane all pass, at one Newton iteration a step.
   * *CC's decomposition.* On a plane the cost of persistent association is zero -- CC-flat' falls to
     6.408e-15 of the drive as rtol tightens. On a cylinder it scales as (1/R)^1.068 against the
     sagitta's 1.000, at a constant prefactor of 0.24. The degradation is the price of curvature, not
     a bug.
   * *The load fraction never measured seating.* It measures load APPLIED. `aim` is the final target,
     so the logged gap is the distance to full seating, and the reach of 0.1094 is a statement about
     load carried and nothing else.
   * *A third of the contact constraints are frozen every step.* 963 to 1,144 of 3,123 held nodes,
     31-37%, have their association update REFUSED at the 0.5 mm jump limit, at every step.
   * *The 39, characterised.* 39 -> 34 -> 27 held nodes move AWAY from their target: scattered over
     103.6 mm rather than one patch, 87% the same nodes between converged fractions, and SHRINKING
     rather than recruiting. They are NOT the refused nodes -- 0 of 39 and 1 of 34 against 12.0 and
     10.7 expected by chance, anti-correlated.

   **OPEN, and honestly so.**
   * The INSTANTANEOUS staleness was measured at ONE fraction only, 0.0625, before the line was
     stopped, and the single reading is recorded here rather than discarded: over the 963 nodes
     refused at that step, **median 0.8174 mm, p90 1.4379 mm, max 47.8232 mm**. At one accepted step
     the accumulated and instantaneous figures necessarily coincide, so this reading carries none of
     the artefact -- and equally it CANNOT distinguish the ratchet, since separating them needs two
     or more steps. One point is not a trend: it says a typical refused constraint was about 0.8 mm
     from the bed after the first step, and says nothing about whether that grows. The measurement
     that would have decided the question was not taken.
   * The RATCHET is neither confirmed nor refuted: a refused node measures its next wanted
     displacement from the same frozen point, so refusal may compound. Its apparent confirmation in
     the accumulated figures was the artefact above.
   * GATE BB' FAILS: zero-sized steps at one fraction give 0.335 -> 0.122 -> 0.0346 -> 0.0364 ->
     0.0374 -> 0.0361 mm, plateauing at ~36 um instead of reaching a fixed point. At a fixed load the
     association and the solve chase each other. The repair does not address it.
   * *The bind.* The per-node 0.5 mm limit freezes 31-37% of nodes per step; the all-or-nothing rule
     refuses 93.8% of whole steps; gate T showed both pure extremes fail. No setting among the three
     is simply correct.

   **THE CAVEAT THAT TRAVELS WITH EVERY SEATING NUMBER FROM THIS LINE**, whether or not anyone chases
   it: the production drive (`stage_dr`) runs `jump_limit_m = 5e-4`, and under it about a third of
   held nodes have their association refused at every step and held stale. `association moved 0.4999
   mm` in those logs is the maximum over ACCEPTED updates and does not show it.


   ## THE UT-ENDOMRI COHORT, IN NUMBERS RATHER THAN AS A SENTENCE

   Measured 2026-09-11 by `scripts/characterise_ut_endomri_cohort.py` over all 124 extracted
   subjects. The endometriosis caveat has travelled with every artefact from this dataset as a
   sentence — *"an endometriosis cohort, NOT a typical-anatomy reference"* — which is true and
   unfalsifiable as written. These are the figures behind it.

   Reference ranges were stated in the script before the data were read, so they could not be
   fitted to it: uterus ~30–120 mL, ovary ~3–10 mL per side, an ovary above ~20 mL in this
   cohort being endometrioma-scale. **Nothing passes or fails on them.** They are the axis the
   cohort is described against, not a gate.

   | | |
   |---|---|
   | uterus labelled | 91 of 124 (absent in 33) |
   | uterus median / IQR | **108.5 mL** / 68.8–154.8 |
   | uterus range | 23.4 – **765.7 mL** |
   | above the 120 mL reference | 38 (**42%**) |
   | above 240 mL | 9 (10%) |
   | largest | **6.4×** the top of the reference range |
   | ovary labelled | 82 subjects (absent in 42), 95 pieces |
   | subjects with a SINGLE ovary | **68 of 82 (83%)** |
   | ovary piece median / range | 7.76 mL / 1.45 – **115.02 mL** |
   | pieces above 20 mL (endometrioma-scale) | 11 (**12%**) |
   | uterus inter-rater Dice (D1, 35 pairs) | median **0.774**, range 0.224–0.901, **26% below 0.70** |
   | subjects carrying extraction flags | 40 of 124 (32%) |

   **What this means for anything derived from here.** 42% of labelled uteri are above an adult
   non-gravid reference; 83% of subjects with an ovary label have only one, which the dataset's
   own README attributes partly to surgical resection; and the raters disagree about where the
   uterus is badly enough on a quarter of D1 pairs that the worst pair overlaps at Dice 0.224.
   **A mean over this cohort is not a typical anatomy**, and no organ taken from it should be
   presented as one.

   ### Two corrections to how this was recorded

   **The manifest asserted a verification it had not performed.** `extract_ut_endomri.py` wrote
   `md5=MD5` — the expected constant — into every manifest unconditionally, including runs under
   `--skip-md5`. The v2-world extract did exactly that: its manifest carries
   `7ace6e1b08efa10d0a1967073b0ba41c` and its own `extract.log` has no `GATE md5` line, so the
   claim was not earned when it was written and a reader could not tell. **The bytes are fine** —
   the archive was hashed independently on 2026-09-11 and matches — but the manifest now records
   `md5_computed` and `md5_verified` separately from `md5_expected`, so a skipped check says so.

   **A gate labelled honestly in a script header was not labelled at all in the artefact.**
   `volume_in_bound` is a coarse typo catcher at 10–1000 mL and the script's header says so; the
   manifest just said `volume_in_bound: true`, on a cohort where it passes 91 of 91 including a
   765.7 mL uterus. Likewise `count_expected`, which asks only `pieces <= 2` — correct, because
   unilateral absence is real here, but not a check on anything. Both now carry their semantics
   in the manifest under `gate_semantics`. Neither gate was wrong; neither was readable from the
   thing it was stamped on. *A caution that lives only in a script header is one revision from
   being lost.*

   ## THE PRODUCTION DRIVE HAS NEVER RUN SINCE THE REPAIR — prediction fixed before it does

   Found 2026-09-11 while resuming this line. `data/derived/female-breast-sliding-v1/s1159/left/dr.log`
   is timestamped **Sep 10 21:50**; the `theta` bound repair landed in
   `ihm/assembly/sliding_contact.py` on **Sep 11 ~03:31** (`6d1fb89`…`025b09a`). `stage_dr` calls
   the repaired `seat_on_bed`. **So the only recorded production failure on the anatomical bed —
   `load stepping stalled at fraction 0.0005: bounded start inverts 2 elements` — is a
   PRE-REPAIR number, and it has been reasoned from since as though it were current.**

   That matters because the repair addresses exactly this mechanism. The bound was
   `on_plane + gap0 + fraction*travel`, which scaled `travel` but not `on_plane`, so a drifted
   association demanded its whole drift **instantly, however small the step** — 7.19 mm over 26
   nodes at a step of size zero. "Bounded start inverts elements at fraction 0.0005" is what that
   bug looks like on a real bed: the tiny step is not the problem, the bound is.

   **PREDICTED, before the re-run:** the production drive gets **past fraction 0.0005**. If it
   stalls there again with the same message, the repair does not reach this failure and the 1.75×
   it is worth on control R′ does not transfer to the anatomical bed — which would be worth more
   than a seat, because every subsequent plan on this line assumes it does.

   **Not predicted, and not to be read into a pass:** how far it gets. Affordability is a separate
   question and the last converged anatomical step is fraction 0.1094. A drive that passes 0.0005
   and then stalls at 0.002 has still refuted the pre-repair reading and still not seated a breast.

   **Also recorded, because it was being planned around:** no breast has ever been seated. There
   is no `judge.json` anywhere in `data/derived/`, so per-breast gates (a)–(d) have **never once
   been evaluated on a real breast** in either the closest-point or the sliding line.

   **And the judge in force is the one at lines 862–871**, not the earlier one: (a) volume within
   1%, (b) every tet J > 0.2, (c) the two solvers within 5% of maximum displacement, (d) no
   flipped base triangle, plus the seat completing. **Base contact (median ≤ 3 mm) and rib points
   inside are REPORTED, not judged** — the sliding boundary condition enforces them by
   construction — and the *fitted-on-some, tested-on-held-out* fold belongs to the rigid-correction
   branch that ended at line 2384. Judging the base deformation by a held-out gate would apply a
   superseded judge.

   **The caveat that travels with whatever this produces:** `stage_dr` runs `jump_limit_m = 5e-4`,
   under which 31–37% of held nodes (963–1,144 of 3,123) have their association refused every step
   and held stale. `association moved 0.4999 mm` is a maximum over ACCEPTED updates and does not
   show it.
