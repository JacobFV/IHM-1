# Exact-source collision and energy-observer preparation

The bounded neutral-pose source scan found no distinct eligible face-interior candidate within 4 mm of a free hair node with absolute face-normal/shaft-direction cosine below 0.8 among the 16 nearest face centers. This limited search is negative evidence, not a global collision-clearance certificate. It loaded the hash-bound 5 MB skin source, used 113,688 KiB peak RSS, and ran in 3.05 s. Original roots, strand mass and geometry are unchanged.

`data/research/hair_collision_face_scan.json` also retains eight exact exterior hand faces. Each face's three nodes has the same existing inferred nearest-bone-envelope owner, `hand_r`. No synthetic face or nearest-face root replacement is introduced.

The proposed isolated pose query targets the center of one such face 10 micrometers from a scalp free node, with its normal perpendicular to the shaft. Five right-arm native coordinates are searched using copied native States. Source shoulder bounds are broad (±10 rad), so the explicitly engineered search restricts flexion/addition/rotation to ±2.5, ±2.5, ±3.1 rad; elbow and pronation retain source bounds. These are search bounds, not calibrated anatomical ROM. There is one starting seed, at most 120 native evaluations, and a 20 s total timeout. The native continuing body's coordinates, transforms and time must remain unchanged. Kinematic feasibility would not certify other-body containment, muscle strain, or an actual collision.

`scripts/prepare_hair_energy_observer.py` prepares a fresh isolated native source copy. Its two narrow changes expose raw `getMusclePotentialEnergy` and a copied-state bounded `hair_pose` transform query. No force law is modified. Compilation is separately gated and queued, with 45 s wall time, 3 GiB address space, one CPU and reduced priority. No default runtime pointer is changed.

Support's `scripts/effective_passive_energy_observation.py` supplies the source-bound Thelen correction. Native `model.calcPotentialEnergy` already contains raw muscle passive potential: add **only** `source_passive_correction_j` to the native total. Do not add the complete corrected muscle total again. Raw native energy remains retained. The known getter defect and the still-investigated wrapped-path virtual-work mismatch are distinct; correcting one does not close the entire energy ledger.

Prepared scripts pass parsing and the existing one-owner source fixture. Native compilation, pose feasibility and subsequent collision acceptance remain separately gated steps. None of this preparation enables a live model or establishes a 20 ms interval domain.

## Actual bounded pose feasibility

The isolated observer compiled in 12.97 s with 1,228,880 KiB peak RSS. The copied-state native pose query then passed in 1.82 s with 94,844 KiB peak RSS. `data/research/hair_hand_pose_acceptance.json` retains the exact source face, scalp sample, solved coordinate values, all queries, build identities and limits. Final center error was 8.48e-11 m and normal error 1.76e-10. A final native observation verified unchanged continuing body transforms, coordinates and time. Both children exited and the shared slot was released. No collision or trajectory was advanced. The solved pose remains an engineered test initial condition, not a whole-body clearance or muscle-strain certificate.
