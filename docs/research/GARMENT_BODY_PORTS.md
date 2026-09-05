# Source-bound garment/body mechanical ports

Implement objective barycentric spring attachments between a garment node and a body surface triangle, with an explicit scalar rest length and stiffness. Derive both endpoint forces from one stored energy, distribute the body force with the exact barycentric weights, and report body force and torque about a stated origin. A support port is an engineered elastic tether; it is not an adhesive or friction law and is not enabled implicitly in the viewer.

Pair selected full-shirt shoulder/neck support nodes and full-shorts waistband nodes to the closest point on the retained canonical skin triangles. Retain exact source triangle/node indices, barycentric weights, initial distance and regional pattern label. Retain additional torso and thigh interface sample pairings as contact query metadata, not tethers. Check source hashes and geometric reconstruction. Nearest atlas surface pairing does not identify anatomical attachment biology.

For the next baseline, prescribe this source body surface at rest and run actual whole-cloth inertia/elasticity, gravity and the existing face-interior contact kernel. Its prescribed-target inverse mass is exactly zero; any positive placeholder masses required by the existing API must be proven irrelevant, never interpreted as body mass. Return every contact and attachment reaction to a body force/wrench ledger, with the equal/opposite support taking it in this prescribed experiment. Record unresolved nearest-edge contacts and repeated discrete contacts, numerical projection work, input/output positions and coefficients. Compare contact/friction-disabled and timestep-refined runs before exporting a shape; do not label a finite transient equilibrium or full containment.

## Verified force and source interfaces

`BarycentricAttachment(garment_node, body_nodes, barycentric, rest_length_m, stiffness_n_m)` and `attachment_forces(garment_positions_m, body_positions_m, attachments, wrench_origin_m=...)` now return both owners' nodal forces, one stored elastic energy, body force/torque and paired conservation residuals. The central spring uses a scalar rest length, so both energy and force transform objectively under a rigid rotation/translation. It rejects the undefined direction of a collapsed nonzero-rest spring. No body mass is invented and no live scene is mutated.

`build_source_pairings(root)` produces 16 optional shirt support pairings near the upper arm-opening/shoulder pattern and 16 shorts waistband pairings, plus 12 torso and 24 hip/thigh contact samples. Every pairing has its original skin triangle ID, three original body node IDs, convex barycentric coordinates, reconstructed source point, distance, closest feature and region label. Source face winding is retained without claiming that the atlas shell's exterior/cavity convention is resolved.

The fresh receipt `data/derived/garment-body-ports-aydnbsru/report.json` passes analytic force, both-owner energy-gradient, equal/opposite force/torque, rigid-objectivity, invalid-weight, face/edge/vertex nearest-point and actual source-reconstruction tests. Its inputs retain the exact source skin, garment mesh, constructor and mechanical code. An independent existing-kernel contact fixture also confirms that a fully prescribed target's response is unchanged for placeholder masses from 1e-9 to 1e9 kg: its inverse mass is zero, so these API placeholders are not physical body masses. Contact reaction impulses are equal/opposite and kinetic loss equals contact dissipation.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 10 .venv/bin/python scripts/verify_garment_body_ports.py
```

The full gravity/contact drape trajectory is **not run or exported yet**. It is deferred while the separate IBM-1 job loads the user's machine; these small checks completed in under one second with one low-priority CPU thread. Existing whole-garment gravity fixtures and the default viewer are unchanged. The force ports and geometric pairings are ready for the planned explicit drape experiment, without treating an elastic support as skin friction or hiding unresolved contacts.
