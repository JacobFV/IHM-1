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
   and 0.962 +/- 0.049 (D2) at ~4.3 mm; their uteri are 52.7-179.4 mL, ordinary ones. Failure
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

