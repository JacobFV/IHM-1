# Physical environments

Environment selections now cross the live body API boundary. `environment_selection`
contains catalogue scene/object/component IDs, which the server validates against the
base environment and dependency rules. Browser geometry is not a physics owner.

The server instantiates fixed compound furniture/wall colliders, spring cloth,
a three-dimensional pillow spring lattice, and six-degree-of-freedom movable balls
and blocks. It samples the canonical skin envelope at 45 mm voxel spacing and
attaches the samples to the existing named native bone owners. There is no second
human mass. Contact forces act at these material stations through
`ArticulatedBodyPlant.advance`, including the force moment about the native segment.
The same forces enter the embodied runtime's existing respiratory load projection;
mechanical observations continue into its neural and metabolic exchange.

Environment integration uses 1 ms substeps and a 20 ms explicit exchange with the
body. Cloth/body, prop/body, cloth/prop and prop/prop interactions exchange reciprocal
impulses. Fixed furniture and the anchored pillow foundation supply external support.
Environment checkpoints are restored alongside body checkpoints on a pre-native
failure. An uncertain native commit retains the existing fail-stop behavior.

The accepted `environment_state` reports the clock, force ports, impulse, meshes,
rigid transforms, selection, coefficients and asset hashes. The browser consumes
these positions; it performs no independent cloth or rigid-body integration.
Factory outputs retain source and selected asset copies. Mattress material and skin
foundation choices reach the existing native factory; mattress support is not
counted again by this environment solver.

## Appearance

Paint relief, wood grain, textile weave, stone, rubber and metal use deterministic
albedo/bump textures and physical materials. Furniture has beveled edges; trees use
instanced leaves within the declared canopy envelopes. Shadows follow the gravity
frame. Camera framing also follows that frame, including supine rooms. Cutaway walls
remain colliders. Bare native floors are drawn at their declared level, without a
second floor overlay on scene ground. Cloth is preconditioned for 1.5 s against held
canonical skin before time zero, with its foot hem and two side corners tucked at the mattress; the preview
and runtime share this preparation. A sampled-skin collision barrier converts
cloth position corrections to reciprocal body impulses.
This initialization is separate from continuing body dynamics and is not equilibrium
certification. Initial pillow geometry has a rounded surface and a volumetric lattice.

## Verified and unverified

- Nine causal solver/rollback tests pass, including native force-port delivery through
  an embodied runtime fixture and reciprocal prop/cloth impulses.
- JavaScript tests verify server mesh updates, stable duplicate instance identities,
  rigid rotation about the correct center, and textured material families.
- All six catalogue scenes complete 3 s against a stationary canonical body fixture.
  This tests environment integration, not whole-body feedback stability.
- A real native OpenSim comparison with and without bedroom contacts changes body
  position by approximately 32.9 mm after the first 20 ms. The 0.2 s native run
  finishes with 82 contact force ports. Retained result:
  `data/derived/environment-native-1snh2tfi/report.json`. This demonstrates a causal
  load path, not realistic transient magnitude or calibrated biomechanical behavior.
- All six scene previews were inspected in Chrome. Screenshots and browser errors
  are retained under `test-results/environment-upgrade/`.
- **The full physiology/IBM/body run is not passing.** The adapter first required a
  rebuild because its retained `native_tissue_ports.h` hash was stale. After rebuild,
  the native signed-muscle budget rejects a decrement of approximately -41.0 W
  against an allowed 8.69 W. The pre-existing test without environment coupling also
  fails at approximately -36.5 W with the same budget guard. No energy check was
  disabled or clamped. `verify_environment_embodied.py` records this failure.

Contact stiffnesses, damping, friction and pillow/textile coefficients are engineering
values, not measured or clinically validated properties. The body surface attachment
is approximate. Fixed cylinders and canopies use bounding boxes; rigid manifolds use
sampled vertices, so edge-edge and high-speed collisions are not comprehensive.
There is no cloth self-collision, resolved cutaneous pressure, fabric thermal
insulation, or separate leaf/blade mechanics. The native support law remains the
owner of bed/floor support. Live temperature changes are unavailable in this adapter
and are disabled in the UI rather than silently accepted.

The cloth implementation is a damped spring model, **not XPBD**. The
[XPBD paper](https://mmacklin.com/xpbd.pdf) is a useful reference for a future compliant
constraint solver with more robust stiffness control; this implementation does not
claim its iteration/timestep properties.

## Reproduce

```sh
PYTHONPATH=. .venv/bin/python scripts/verify_environment_dynamics.py
PYTHONPATH=. .venv/bin/python scripts/verify_environment_scenes.py
PYTHONPATH=. .venv/bin/python scripts/verify_environment_native.py
PYTHONPATH=. .venv/bin/python scripts/verify_environment_embodied.py
npm --prefix app test
npm --prefix app run build
node app/test/environment-preview.mjs
```

The native and embodied scripts create fresh retained run directories and close
only their own owners. The preview script starts and closes its own local server on
port 8766. The embodied command currently reproduces the budget failure above.
