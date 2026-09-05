# Clothing components and mechanics

## Design and implementation plan

The assembled generic body remains the reference. Clothing is an engineered
addition; it does not require a complete pre-existing garment or anatomy atlas.
The default view will enable an opaque shirt and shorts independently. Clothing
is separate from the anatomical layer filters, opacity, and native thermal
`clo` input. Toggling a rendered garment must not silently change physiological
insulation or erase underlying anatomy.

The display materialization samples transverse sections of the actual canonical
skin triangles, separates disconnected limbs, and fits garment rings around
their convex envelopes with explicitly assumed ease and pattern boundaries.
A sleeveless shirt has neck and arm openings. Shorts have a sewn crotch and two
leg openings. Each materialization records its source geometry identity and
construction assumptions. The body’s recorded thoracic displacement may deform
the shirt kinematically; that is not a solved cloth/body contact trajectory.

An independent Python mechanical component will provide elastic edge forces,
lumped cloth mass, unilateral point contact, and static/kinetic Coulomb friction
against prescribed body surfaces. It will return equal-and-opposite body
reactions and contact dissipation. Its coefficients are declared engineering
priors, not measured fabric/skin calibration. This first mechanics implementation
does not replace the body’s mechanical state owner and does not move genital
tissue. Physical tucking requires deformable tissue, attachment constraints,
two-way contact and verified stresses; geometric coverage is not that result.

Implementation sequence:

1. Add failing geometry/component tests and mechanical force/friction tests.
2. Implement source-conditioned garment meshes and separate app controls.
3. Implement and validate elastic/contact mechanics with generated audit data.
4. Run unit, native Python, build and real-browser clothing checks; retain results
   and list remaining coupling limitations here.

## Primary methodological sources

- Baraff and Witkin, [Large Steps in Cloth Simulation](https://publications.ri.cmu.edu/large-steps-in-cloth-simulation),
  1998: triangular cloth, internal elastic forces and implicit integration. The
  initial component here uses explicit substeps and edge springs, not their
  continuum formulation or implicit solver.
- Bridson, Fedkiw and Anderson,
  [Robust Treatment of Collisions, Contact and Friction for Cloth Animation](https://lightfield.stanford.edu/papers/cloth-sig02/),
  2002: contact thickness, repulsion/collision handling and Coulomb friction. The
  initial point contact implementation is narrower than their robust continuous
  collision and self-collision algorithm.

These papers support algorithmic choices; they do not supply calibrated fabric
coefficients or establish clinical validity.

## Implemented and retained

`app/src/clothing.js` builds the same garment mesh in the browser and in the
retained materialization. The held skin contains 102,467 vertices and 203,382
triangles; the output shirt has 2,038 positive-mass nodes and 3,840 triangles,
and the shorts have 2,527 nodes and 4,864 triangles. Boundary-loop checks require
four shirt openings and three shorts openings. The crotch seam shares vertex
indices between the two legs. Opaque garments use their own switches and remain
independent of anatomy opacity. The initial anatomy view includes the skin.

The source's joined axillary sections and nearby hand sections can contaminate
torso/hip fits. The construction therefore records a transverse regional mask
of ±190 mm, reduced to ±175 mm above the axilla, scaled with body height. This
is a declared pattern prior, not an inferred biological boundary. A held-source
regression excludes the initial arm/hand-induced flares. There is no anatomical
deletion or genital displacement in the fitting code.

`ihm/assembly/clothing.py` exposes finite cloth mass and nodal state. Internal
forces derive from the edge energy `k (length - rest_length)^2 / 2`, and central
edge forces conserve net force and torque. The integration step is transactional,
has a conservative linearized timestep bound, and returns a numerical energy
defect instead of assigning it to a physiological reservoir. The planar contact
operator uses zero restitution, a static Coulomb cone and kinetic friction;
it returns per-node opposite body impulses and prescribed-support work. These
body reactions are **not yet applied to the canonical mechanical owner**.

Reproduce the mechanical checks and source-pinned garment artifacts with:

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_clothing.py
```

Retained outputs are `data/derived/clothing/garments.json` and
`data/derived/clothing/verification.json`, with source and constructor SHA-256
identities. The current independent checks report:

- Elastic force versus energy-gradient error: `6.46e-11 N`; net torque zero.
- Moving-support contact energy residual: `1.39e-17 J`.
- Supported-strip momentum residual: at most `5.09e-21 N s`.
- Elastic time-refinement maximum energy errors over 0.1 s: 1.55%, 0.77%, and
  0.38% at 100, 50 and 25 microseconds. This is a first-order explicit solver.

The supported contact test reports its finite-step projection/impact energy
defect explicitly; the force and work checks do not establish exact time-domain
energy conservation. Point contact has no triangle-level continuous collision
detection, self-contact, anisotropic fabric constitutive law, bending stiffness,
or measured friction calibration. Garment mass estimates use an assumed
0.18 kg/m², yielding approximately 82 g per garment; these are engineering
estimates, not clothing measurements.

## Next coupling work

The garment node API is compatible with finite-mass tissue contact: positions,
velocities, nodal masses, triangle connectivity, and external nodal forces. A
whole-body extension must apply the opposite impulse to the actual tissue owner,
respect attachments and rest lengths, and validate penetration, pressure,
momentum, energy and timestep convergence. In particular, a compliant genital
domain must deform under garment contact before any result can be called
physical tucking. The current displayed shirt follows the recorded thoracic
field kinematically and does not claim that coupled result.

## App acceptance, 2026-09-05

`npm test` passed all 23 unit tests, `npm run build` passed, and the clothing,
body, and workspace Playwright suites passed all nine checks against the real
loopback app. The acceptance covers default opaque dress, independent garment
switches, preserved anatomy, model-switch isolation, real body playback and
responsive panel behavior. A reviewed screenshot is retained at
`artifacts/verification/clothing/clothed-human.png`. This visual acceptance does
not substitute for coupled mechanical validation.
