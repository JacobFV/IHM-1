# The wardrobe's rigid-envelope assumption

Handoff for the simulation overhaul. Every number here is measured on the
current tree, and the commands that produce each one are given. Nothing in this
document is a plan that has been agreed; it states a gap and the decision the
overhaul has to make.

## What the wardrobe does today

`scripts/fit_garments_to_envelope.py` registers each of the 33 acquired garments
onto `data/derived/outer-envelope/outer-envelope.npz` — a frozen, watertight
**surface**, 54,360 vertices — by a similarity transform, a piecewise-rigid limb
pose correction and a Laplacian-regularised shrinkwrap with a per-garment
standoff of 2 to 10 mm. `scripts/build_garment_wardrobe.py` then runs each
garment, **alone**, as elastic cloth against that same static surface.

The pipeline's own verification states the assumption plainly:

> `'body': 'static prescribed surface; no whole-body reaction applied'`

That is the only line in either script that mentions the body as a mechanical
participant. Neither script references `step_coupled`, `SymmetricContact`, the
whole-body material domain, or any tetrahedral tissue:

```sh
grep -n "step_coupled\|SymmetricContact\|whole-body\|material-domains" \
  scripts/build_garment_wardrobe.py scripts/fit_garments_to_envelope.py
```

Three consequences follow directly, and all three are now measured rather than
suspected:

* the body cannot yield, so cloth must be stretched over every feature;
* garments never meet each other — **150 of 402** wearable combinations
  interpenetrate (`post-simulation-layering.json`);
* garments never meet themselves — **44,766** self-intersecting face pairs
  across **33 of 33** garments (`post-simulation-penetration.json`).

## The general machinery already exists and is unused

This is the substance of the gap. The repository already contains the general
form of what the wardrobe special-cases:

| capability | where | generality |
| --- | --- | --- |
| N-body elastic step | `ihm/assembly/contact_dynamics.py:164` `step_coupled` | arbitrary dict of mass owners on one clock, unique material ownership, shared explicit stability limit |
| interface declaration | `ihm/assembly/contact_broadphase.py:458` `SymmetricContact` + `expand_symmetric` | any two owners; one declaration expanded to both directional passes so they cannot disagree |
| whole-body soft continuum | `data/derived/material-domains/whole-body-0.005m` | 1,829,976 tets, 501,046 nodes, 2,404 anatomical sources (a 10 mm variant has 229,692 tets) |
| calibrated genital tissue | `ihm/assembly/penile_volume.py`, Khorshidi 2024 ex vivo | per-tissue constitutive laws on a 2 mm pelvic domain |

`ihm/assembly/garment_contact.py` is **not** the counter-example. Its own
docstring says it is "a constrained front-panel experiment, not whole-garment
donning or calibrated genital containment": a hand-carved anterior panel of the
*old generated shorts* against a hand-carved pelvic tet subset. It is a second
special case, not the general path.

The assessment records how far the general path has actually been taken:

> "Two deformable bodies is the largest body-body contact ever exercised in this
> repo." — `soft-body-coupling-assessment-v1/existing_machinery.json`

That exercise was 135 panel nodes against pelvic tissue, 0.24 s simulated in
359 s of wall time. The wardrobe after resampling is **291,751 nodes and
576,136 faces** across 33 garments, and the body it would have to meet is a
1.8 M-tet continuum. The distance between those two numbers is the honest size
of this problem, and it is why nobody has done it.

## The sharpest instance: the genital region

Reported as bad fitting around the penis, and the measurement bears it out. In
the band y ∈ [−0.11, 0.00], anterior, near the midline, essentially *all*
remaining body penetration in lower-body garments is concentrated:

| garment | penetrating samples in region | of total | concentration |
| --- | --- | --- | --- |
| briefs | 15 | 15 (100%) | 6.0× |
| casual-trousers | 91 | 97 (93.8%) | 19.6× |
| suit-trousers | 82 | 86 (95.3%) | 21.2× |
| denim-shorts | 2 | 2 (100%) | 8.0× |

Two causes, and they need separating because only one is a fitting defect.

**The asset cause.** The envelope reaches up to 80 mm further anterior than the
MakeHuman authoring body at the same height; 43% of envelope points in the
region lie outside the authoring body, median 15.5 mm, max 37.0 mm. The garments
were drawn on that authoring body, so they carry **no material** for the
feature. Several are explicitly female patterns (`f_panties_01.obj`,
`wolgade_female_top_01`). No amount of resolution creates cloth that was never
drawn.

**The rigid-body cause.** The tissue there is soft and is being treated as an
immovable obstacle. Measured, the region is a vertically elongated mass — 103 mm
tall, 84 mm wide, 35 mm deep off the local body wall — not a horizontal
projection. Under a garment, real tissue displaces; here the cloth must instead
stretch over a frozen surface.

A tempting third reading is a mis-calibrated constant: `SAMPLE_MAX_DEPTH_M =
0.010` means the shrinkwrap only corrects surface samples between 0 and −10 mm
deep, and every penetration in this region is 13–17 mm, so the fitter is
currently ignoring exactly these failures. Raising it *would* clean the receipt.
**It should not be raised.** Forcing cloth to conform to a rigid obstacle
produces a clean number for the wrong physics. The constant is a symptom of the
frame, not the defect.

For the same reason, a tucked rest configuration must not be authored into the
envelope. Posing it by hand repeats the original error — asserting a
configuration instead of deriving it. Under coupling, the resting configuration
is whatever tissue-and-cloth equilibrium produces.

## What the shrinkwrap actually is

Registration exists because the cloth solve cannot don a garment. The solve has
no bending law, no self-contact, no continuous collision detection, and runs
5.7 to 30.1 ms of simulated time per garment — 19 of 33 garments hit the
1400-step cap before reaching the requested 30 ms. Given that, a garment cannot
be placed loosely and settled, so it is instead *fitted* onto a frozen body and
the result declared.

Under general dynamics there is no fitting step. A garment is placed around the
body as an ordinary mass owner and the coupled solve settles it. The whole
registration stage — similarity, pose correction, shrinkwrap, standoff, fit
modes, and the resolution-scaled regularisation that makes them behave — is
scaffolding for the absence of that solve.

The served geometry is worth noting here: `build_garment_wardrobe.py:596` writes
`geometry/*.json.gz` from the **fitted** positions, not the simulated ones. What
the viewer shows has never passed through a cloth solve at all.

## What holds regardless of architecture

Two pieces of today's work are frame-independent and should survive the
overhaul.

**The receipts.** The penetration check tested vertices only and reported ~0 mm
for 29 of 33 garments; restated on face centroids and edge midpoints, 0 of 33
were clear. Self-intersection and a between-garment layering audit are new.
These are the only instruments that can score a coupled solve — in particular,
the layering audit is what would say whether a tuck *emerged* correctly rather
than being asserted. `--stage fitted|simulated` selects which state is audited.

```sh
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/fit_garments_to_envelope.py --layering <ids>
.venv/bin/python scripts/verify_garment_remesh.py
```

**The resampling.** `ihm/assembly/garment_remesh.py` rebuilds each garment at a
uniform edge length, per connected component, projected back onto its own
surface and its own boundary polyline. This matters *more* under coupling than
under shrinkwrap: well-shaped elements are a contact-stability requirement, not
a cosmetic one. Across 33 garments the median body-penetration rate fell 0.439%
→ 0.048% and median edge-ratio p99 fell 2.71 → 1.91, with source-surface
deviation of 0.0000 mm.

Two garments are un-improved and carried as authored, with the refusal recorded:
`cloche-hat` (resampling produced 43 non-manifold edges) and `wool-trousers`
(1). Both had non-manifold input already.

The target is a ceiling on edge length, not a value to drive every garment to.
The gloves are authored at 6.4 and 7.4 mm, finer than the 10 mm target, and an
earlier revision coarsened them — costing `gloves` faces and raising its body
penetration from 0.285% to 0.479%. Resampling now refines a garment or leaves
its resolution alone, and the receipt records which was used. After that fix
`long-gloves` reaches 0.462%, better than its authored 0.524%; `gloves` reaches
0.311% against 0.285% authored and is the one garment of 33 still marginally
worse than it was, though its worst corner angle improved from 5.39° to 28.94°.

## The decision this leaves

`data/derived/soft-body-coupling-assessment-v1/recommendation.json` ranks five
architectures for coupling the 2,404 per-entity bodies, recommends per-entity
bodies with a classified interface graph (4,934 of 19,440 adjacent pairs
licensed to slide, tie-by-default), and names the real risk as interface
classification error rather than solver cost.

**Garments appear nowhere in it.** All five architectures concern body-internal
interfaces. A garment is a different kind of participant — an open shell, not a
closed volume; it slides against skin over nearly its whole area; it arrives and
leaves; and it must contact *other garments* and itself, which no body-internal
interface has to do. The assessment's contact-patch budget of 6.842 m² is
body-internal only; the wardrobe adds its own surface on top.

So the overhaul has to settle three things this document deliberately does not:

1. whether a garment is an owner in the same graph as the 2,404 bodies, or a
   separate layer coupled at the skin surface;
2. what sets the common clock, given `step_coupled` enforces one explicit
   stability limit across all bodies and cloth edges are far stiffer per unit
   mass than tissue tets;
3. whether the acquired garment assets survive at all, given they were drawn for
   a body this one is not — which the genital region makes unavoidable, and
   which no solver fixes.

## Not claimed

No coupled garment/body solve has been run. No garment has been simulated
against tissue. The resampled wardrobe has not been re-simulated even against
the static surface, so `post-simulation-*.json` still describes the pre-resample
cloth states. The region band used for the genital measurement is an
axis-aligned box in the display frame, not a segmented anatomical boundary.
