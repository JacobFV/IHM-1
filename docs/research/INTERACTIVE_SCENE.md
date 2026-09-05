# Interactive scene mechanics

The workbench provides Select, Gimbal and Force modes. The latter two apply a
bounded cursor spring to the selected mechanical owner; the gimbal does not
teleport anatomy. Studio, floor and bed environments use separate, recorded
live sessions. Historical physiological replay pauses while a live scene owns
the display, and Reset releases that ownership.

The body uses 2,408 canonical linked translations and affine tissue regions.
Reference orientations remain constrained: off-centroid body moments are
constraint reactions, not free articulated rotations. Floor and bed use named
ideal supports about an **assumed balanced gravity preload**. Gravitational
prestress, support distribution and body-surface contact are not solved. This
is an incremental model, not validated standing, lying or walking dynamics.

The free sphere has finite mass and rotational inertia, ground contact,
inelastic normal impulse and Coulomb friction. Its force point must lie inside
or on the sphere and follows material rotation across solver substeps and
cursor commands. Body-object contact, clothing contact and feedback into
physiology are not yet coupled to this scene; the monitor states these limits.

Each session retains the exact consumed canonical mechanics bytes, import-time
source bytes and normalized loaded-code fingerprints. Changed source files
require restarting the server. Immutable compressed events are atomically
published with predecessor hashes. Failed publication rolls back mechanical
state; failed close preserves the active session. Strict integer sequences
reject duplicated advances. Following an uncertain HTTP response the client
reads current state and pauses, without replaying a potentially committed
force.

Endpoints are local-workbench-only: GET `/api/scene/catalog`, POST
`/api/scene/sessions`, GET `/api/scene/sessions/{id}`, and POST
`/api/scene/sessions/{id}/step` or `/close`. Events are retained under
`data/derived/interactive-scenes/{id}/`. A server restart does not resume an
in-memory session from those records.

The actual-body force check retained at
`data/derived/audits/interactive-scene-wdowc0_3/verification.json` verifies a
0.02 N s total impulse and sphere energy/friction accounting. Later lifecycle
hardening uses a one-owner fixture rather than repeating whole-body runs:

```sh
nice -n 10 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python scripts/verify_interactive_scene.py --light
nice -n 10 node --test app/test/scene-interaction.test.js
```

To share the machine with IBM-1 training, rendering is capped at 30 FPS with
pixel ratio 1, geometry loads use two concurrent requests, and hidden tabs
perform no render or simulation updates. Hair strand dynamics are opt-in;
visible synthesized strands remain available without continuous rod solves.
Source meshes and full follicle populations retain their original precision.
No background scene stepping occurs without the browser's explicit Start.
