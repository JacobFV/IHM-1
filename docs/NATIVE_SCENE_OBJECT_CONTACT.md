# the body can feel the ball

Measured with `scripts/scene_object_contact.py`; artifact at
`data/derived/scene-object-contact/report.json` (gitignored, regenerable in about
five minutes).

## The question

The environment catalogue's own `not_selectable` list says:

> **body-object contact** — "No engine solves it. The interactive scene declares
> `body_object_contact=false` and objects contact only the environment plane.
> Objects are placed and simulated as free bodies; the body cannot yet rest on,
> push or grasp them."

Is that a missing solver, or a connection nobody made?

## A connection nobody made

It is accurate about `ihm/assembly/interactive_scene.py`. That file integrates
its 0.065 m / 0.4 kg ball in **Python**, with `sphere_plane_step` against a half
space, and never against the body.

It is not accurate about what is available. Three things were already here:

* `SimTK::CollisionDetectionAlgorithm::SphereSphere` **and** `SphereTriangleMesh`
  are compiled into the installed Simbody
  (`data/runtime/opensim/install/simbody/include/simbody/simmath/internal/CollisionDetectionAlgorithm.h`).
* `OpenSim::HuntCrossleyForce` exposes them over named `ContactGeometry` —
  `ContactParameters` takes a list of geometry names, not a sphere/half-space
  socket pair the way `SmoothSphereHalfSpaceForce` does.
* The native plant already carries **28 contact elements over 20 bodies** in the
  `upright` environment (`docs/NATIVE_CONTACT_INVENTORY.md`).

So: **no new solver.** The answer to "what does it take" is a contact set that
names both sides.

## What was built

`NativeMechanicalStream` takes `scene_objects=` — free rigid spheres given a
`Body`, a `FreeJoint` to ground, a `ContactSphere`, and:

* ground contact through the **same** transferred source foot law the rest of the
  environment uses, so the object is not standing on a different floor;
* body contact through **one `HuntCrossleyForce` per (object, body element)
  pair**.

One force per pair, never one force over all the geometries at once. A single
contact set would also test the body's own elements against each other, and
those interpenetrate permanently by construction — the torso proxy has a 0.258 m
radius and its centre sits about 0.3 m from the pelvis proxy's 0.095 m — so one
shared set would invent a large permanent self-contact force. Pairing confines
every test to object-against-body, which also keeps object-object contact out of
scope, as instructed.

The pairing includes `ContactMesh` as well as `ContactSphere`, so an object
touches whatever the body is currently wearing: the inertia-inscribed proxies
today, and the real segment meshes wherever those are loaded.

Contact material is transferred from the source foot contact set — stiffness
1e6 Pa, dissipation 2 s/m, friction 0.8/0.8/0.5, transition velocity 0.2 m/s —
not invented.

## What it measures

**Drop.** The catalogue's own ball released 0.9 m above a prone body:

| body element | frame | first contact | peak force | contact frames |
|---|---|---|---|---|
| `fall_proxy_torso` | `torso` | 0.27 s | 47.2 N | 45 |
| `fall_proxy_pelvis` | `pelvis` | 0.91 s | **218.1 N** | 54 |
| `fall_proxy_femur_l` | `femur_l` | 1.20 s | 16.3 N | 37 |
| `fall_proxy_femur_r` | `femur_r` | 1.20 s | 15.2 N | 37 |

The ball lands on the chest, rolls down the back, and comes to rest in the notch
between the two thighs, touching both. Four named body elements, in the order a
ball rolling down a body would touch them.

**Push.** The same ball resting on the floor 0.75 m ahead of a crawling body:

* still at exactly its start position through 2.47 s while the body crawls toward it
* contact on `fall_proxy_torso` from **2.48 s to 2.77 s**, peak **18.3 N**
* contact impulse on the object **(+1.277, −2.388, −0.028) N·s**
* ball ends at **0.761 m** from where it started, having been struck by a chest

**Control.** The identical run with the ball placed at x = 3.0 m, beyond reach:
**0 contact frames, 0.000 m horizontal displacement.** Both arms show the same
3.76 mm of vertical settling into the floor's compliance, which is why the
discriminator reported is horizontal displacement and not peak speed — the
settling drop gives *both* arms a vertical transient.

## Two things this result is not

**The impulse is a vector.** A running sum of |F|·dt over this rolling contact is
2.709 N·s, which is 36× the ball's actual momentum change (0.074 kg·m/s peak),
because most of the contact force is normal to the motion and goes straight into
the floor. The report carries the vector integral and the object's peak momentum
side by side so the two cannot be confused.

**The surfaces are not skin.** The ball is touching `fall_proxy_torso`: a 0.258 m
sphere inscribed in the torso's inertia ellipsoid at its centre of mass. A ball
rolls off a body-sized ball much as it would roll off a chest, and the force is a
real contact force between two real geometries — but nothing here may be reported
as the body feeling a chest. That gap closes when segment contact meshes replace
the proxies, and the pairing above already accepts them.

## Still out of scope, and why

Object-object contact (each object pairs only with body geometry), non-sphere
colliders, walls, terrain. `not_selectable` documents each; none of them is
touched here.
