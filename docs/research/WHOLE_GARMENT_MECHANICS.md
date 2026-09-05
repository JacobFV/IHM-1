# Whole-garment elastic materialization

## Bounded implementation plan

Retain every node and face of the existing source-conditioned shirt and shorts. Validate the retained garment bytes against their canonical skin geometry and constructor identities. Audit connected components, oriented edge incidence, degenerate faces, and ordered boundary loops before creating a mechanical owner. Classify the shirt's lower hem, neck and two arm openings and the shorts' waistband and two leg hems using their retained generic frame. These classifications are pattern priors, not newly measured anatomical boundaries.

Materialize each entire garment with the existing SI `Cloth` energy and area-derived nodal masses using explicit caller-supplied density and spring stiffness. Preserve exact source indices, with no geometry decimation or duplicate material owners. Existing shared indices are perfectly connected seams; no added stitch stiffness or seam strength is implied.

Verify topology with synthetic open/invalid examples and actual source meshes. Verify full-garment force gradients, net force/torque and mass. Run a short gravity-loaded hanging fixture with the highest boundary constrained, compared with free flight, retaining all positions, velocities, support impulses, energy defects and exact inputs. This is a support/load-transmission fixture, not a worn-body contact experiment. Whole-garment body/self-contact and nonlinear fabric bending remain subsequent gates; the existing viewer and prior artifacts are untouched.

## Implemented and verified

`ihm/assembly/whole_garment.py` exposes `materialize_garments(root, areal_density_kg_m2=..., edge_stiffness_n_m=...) -> (bodies, identity)`. Each returned `Cloth` owns the complete source mesh, original node/face index mapping, material owner ID and pattern audit. Source and constructor hash mismatches reject materialization. The identity record hashes the garment artifact itself; this is an exact retained materialization with verified upstream references, not an independent reconstruction of every constructor operation.

The source shirt has 2,038 vertices, 3,840 triangles, one connected component and four boundary loops; the shorts have 2,527 vertices, 4,864 triangles, one connected component and three loops. Every interior edge has exactly two consistently oriented incident faces. No source faces or coordinates are changed. Topological checks do not certify absence of self-intersection.

The fresh run `data/derived/whole-garments-uc8tp3zc/report.json` retains both 20 ms gravity experiments per garment, complete sampled state and detached exact input bytes including the skin and constructor. At the explicitly assumed 0.18 kg/m² density, mesh masses are 82.1916 g and 83.1802 g. Virtual-work force-gradient errors are 5.51e-13 N and 1.08e-13 N. Free cases accelerate together under gravity without elastic strain; constrained upper boundaries transmit load and develop positive spring energy. Upward support impulses over 20 ms are 0.001308 N·s and 0.003215 N·s. No settling, convergence of drape, cloth-body friction or body containment is asserted from this short transient.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_whole_garments.py
```

A subsequent worn-drape solve must replace the demonstration support constraints with stated torso/shoulder/hip/thigh contact or force ports. Source-bound finite-mass reactions, edge/vertex/self-contact and seam/bending constitutive evidence remain separate requirements; this module does not silently attach garments to skin.
