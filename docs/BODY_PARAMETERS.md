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
   of which can be produced by transforming male ones. **Moved, as of 2026-09-10.**
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

   Both carry one ovary piece, inside the ring. **One threshold was not fixed in advance:**
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
   (618 mL), predicted to fail it the same way. The sternum
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
