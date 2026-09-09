# what the body can feel, per contact element

Measured with `scripts/audit_native_contact_inventory.py`; artifact at
`data/derived/contact-inventory/inventory.json` (gitignored, regenerable in
about 40 s).

## The claim this corrects

> The body's ENTIRE contact set: `foot_load_fraction {r,l}`,
> `foot_contact_force_n {r,l}`, `foot_centre_m {r,l}`, `fall_support_force_n`
> scalar. Two feet and a fall-catch plane. Prone locomotion is not merely
> untuned, it is unrepresentable: a crawling body would pass through the floor
> everywhere except its soles.

Those four fields are what **`scripts/walk_gait.py` chooses to record**. They are
a projection, not the contact set. The engine emits a full per-element `contacts`
array — name, body frame, force vector, moment, centre, radius — and
`walk_gait.fall_load()` sums every non-foot element into the one scalar.

## What is actually there

In the `upright` environment the engine installs, unconditionally
(`scripts/native_mechanical_stream.cpp`, the `environment=="upright"` branch):

* the 12 source foot contact spheres, and
* one `SmoothSphereHalfSpaceForce` against the floor for **every** non-foot body,
  with a radius inscribed in that segment's inertia ellipsoid at its centre of
  mass.

**28 contact elements over 20 bodies**: `calcn_{l,r}`, `toes_{l,r}`, `femur_{l,r}`,
`patella_{l,r}`, `tibia_{l,r}`, `pelvis`, `torso`, `humerus_{l,r}`, `ulna_{l,r}`,
`radius_{l,r}`, `hand_{l,r}`.

Dropped face down and left to settle for 1.5 s, the body comes to rest and does
not pass through the floor:

| element | force |
|---|---|
| `fall_support_torso` | 326.0 N |
| `fall_support_femur_r` | 126.1 N |
| `fall_support_femur_l` | 120.6 N |
| `fall_support_pelvis` | 85.1 N |
| `contactMedialToe_l` | 32.3 N |
| `contactMedialToe_r` | 31.5 N |
| `fall_support_hand_l` | 20.2 N |
| `fall_support_hand_r` | 20.1 N |

Total contact force **761.5 N** against a body weight of **761.4 N**. The floor
carries the whole body.

`fall_support_force_n` is not a catch surface that ends a run by construction.
`walk_gait.diverged()` *chooses* to end a run when it exceeds 5 N, because for
bipedal walking any trunk contact is a fall. That is a controller convention.
`scripts/crawl.py` does not use it and runs an 8 s prone horizon to completion.

## So: configuration change, or build?

**Neither, for contact to exist — it already does.** For contact to be
*anatomical*, a build.

Three tiers exist, and the gap is fidelity, not existence:

1. **`upright`, on now, no configuration.** 28 elements. Bears full body weight
   prone. But every non-foot element is *one* sphere at the segment centre of
   mass: a 0.258 m ball at the chest, 0.086 m balls at mid-thigh, 0.023 m balls
   at the hands. A prone body rests on four balls, not on a chest, two knees and
   two forearms. Knee contact in particular is not at the knee — the nearest
   elements are the mid-femur and mid-tibia spheres.

2. **`supine` + `surface_contact_manifest`.** A 21,381-point posterior skin
   quadrature taken from the canonical surface, with layered confined
   neo-Hookean skin over an optional measured mattress. Real skin geometry. But
   it is **posterior-only** by construction (`posterior_selection`: minimum
   source-X triangle per 5 mm raster cell) and `scripts/native_surface_foundation.h`
   hardcodes the support normal to the source X axis — `penetration = plane -
   location[0]`, `force = (normal, 0, 0)`. `NativeMechanicalStream.__init__`
   enforces `environment == 'supine'` to match. It models lying on your back and
   nothing else.

3. **Anatomical ventral and limb contact, for prone locomotion.** Needs (a) a
   ventral/limb quadrature selection alongside the posterior one, and (b) the
   foundation's support normal made a parameter instead of axis 0, so a
   quadrature can bear against a floor under `upright` gravity. Both are C++ and
   a new manifest builder. **That is a build**, not a configuration change.

## What a prone result may and may not claim

It may claim the body is supported by, and pushes against, real contact
elements on the trunk, thighs, hands and feet. It may not claim the body feels a
chest, knees or forearms: those surfaces are inertia-inscribed engineering
proxies, and any report of prone locomotion has to say so. `scripts/crawl.py`
carries that sentence in every report it writes.
