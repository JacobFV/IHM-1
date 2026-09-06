# Anatomical registration residuals

The existing single proper rigid world transform is physically consistent but anatomically approximate. The retained initial/step audit gives **62.356 mm RMS / 132.032 mm maximum** COM-to-bone-envelope-center residual. Those points are not homologous landmarks, so this is a fit diagnostic, not anatomical accuracy.

## Reproduce without native execution

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_anatomical_joint_registration.py
```

Optional `--output /tmp/fresh-registration-audit.json` creates a fresh JSON file exclusively. The script reads canonical mechanics metadata, the retained supine smoke state, its retained source XML, and execution manifest. Each input is capped at 16 MiB; no geometry arrays, native process, builds, or full anatomy assets are loaded. JSON includes SHA-256 identities, input sizes, all 21 inter-body joint records and six endpoint proxies. It verifies the retained XML against its execution manifest. This is a historical state audit, not a new simulation or a validation of current native binaries.

## Findings

The common map is identical across all 22 bodies, its rotation is proper, and joint-origin distances survive transformation to within 1e-12 m. Pin/weld frame origins agree in both retained states to at most 2.78e-17 m. Source hip and walker-knee frame origins also coincide initially. Patellofemoral offset origins are separated by 56.023 mm: the CustomJoint prescribes relative translation, so treating that separation as a solver closure failure would be incorrect.

Distances below are from the mapped source joint's corresponding parent/child offset origin to the **owning canonical bone group's AABB**. They are lower bounds to bone-surface distance, not joint-center calibration errors. A zero value only means the point lies inside the box.

| Joint | Parent / child outside-box distance, right (mm) | Left (mm) |
|---|---:|---:|
| Hip | 0.000 / 5.249 | 0.000 / 5.886 |
| Walker knee | 0.000 / 25.160 | 0.000 / 24.594 |
| Ankle | 55.531 / 30.372 | 55.927 / 30.659 |
| Shoulder (`acromial`) | 60.213 / 32.443 | 64.453 / 32.080 |
| Elbow | 4.403 / 8.347 | 6.022 / 7.768 |
| Radius–hand | 8.129 / 44.865 | 10.231 / 42.944 |
| Subtalar | 108.776 / 64.646 | 108.466 / 64.724 |
| MTP | 84.196 / 57.646 | 85.342 / 58.282 |

A second deliberately coarse screen pairs the distal Y-face center of the parent envelope with the proximal Y-face center of the child envelope:

| Envelope pair | Endpoint separation R / L (mm) | Signed Y gap R / L (mm) |
|---|---:|---:|
| Femur–tibia | 19.769 / 19.872 | -5.770 / -5.806 |
| Tibia–talus | 11.363 / 11.226 | -9.148 / -9.193 |
| Humerus–ulna | 40.118 / 40.092 | -23.250 / -23.260 |

Negative Y gap means overlap of projected intervals, not bone penetration. Envelope face centers are not articular points; these numbers cannot establish physical surface clearance or soft-tissue continuity. The ankle native origins are 70–75 mm from these envelope endpoint proxies despite the smaller envelope-to-envelope separation. Canonical geometry may look adjacent while its native pivot is misplaced.

The worst center correspondences are torso (132.032 mm), hands (99.847/98.394 mm), and toes (93.822/86.848 mm). Canonical `source_files` names `opensim__Rajagopal__Rajagopal2016.json` and `native_corrected/baseline/mechanics.json`, whereas the audited plant consumes `example3DWalking/subject_walk_scaled.osim`. Reusing the retained per-segment scale/translation is therefore neither a justified model calibration nor a valid replacement for the common world map. For example, the native tibia's local ankle offset has a 465.174 mm Y component, versus approximately 379.957 mm canonical tibia envelope height; these are different measurements, but their mismatch motivates homologous landmark extraction.

## Concrete correction proposal

1. Retain exactly one world transform `C` for native bodies, gravity, forces, contacts and environment. Keep native segments as the sole inertial owners; canonical mass remains a material-partition/reference quantity. Do not restore segment-specific world maps or integrate canonical mesh mass on top of native mass.
2. Build a versioned landmark record tied to **this exact source XML** and canonical geometry hashes. Retained source joint frames give kinematic pivot/axis candidates. Use matching anatomical features on source and canonical bone surfaces, with explicit uncertainty: femoral head sphere centers and acetabular centers for hips; femoral condylar axis plus tibial plateau for knees; malleolar axis/talar dome for ankles; humeral head/glenoid for shoulders; trochlear axis/ulnar notch for elbows; distal radial articular surface and proximal carpals for wrists. Do not use COM, whole torso envelopes, or whole toe groups as substitute joint centers. The source shoulder label alone is not proof that its pivot matches a glenohumeral center.
3. Fit one weighted proper rigid transform using those homologous points/axes, report per-landmark and held-out errors, and compare with the current center fit. Rigid fitting cannot remove body-proportion or rest-pose mismatch. First distinguish pose mismatch from segment-length mismatch; where needed calibrate a single consistent native reference pose and joint geometry using shared adjacent-segment landmarks. Updating native lengths requires consistent muscle/path, contact, COM and inertia treatment with explicit provenance and preserved total target mass, not merely visual bone resizing.
4. Recompute each canonical material embedding as `inverse(T_body_reference) * inverse(C)` after the consistent native reference is set. Segment-local embeddings are legitimate; segment-dependent definitions of world are not. If canonical meshes are adjusted, constrain opposing joint surfaces together and preserve their rest geometry provenance. Whole-surface soft tissues crossing joints need a documented multi-support/deformable treatment; assigning one rigid owner cannot guarantee continuity.
5. Gate a future implementation on common-frame/virtual-power invariants already exercised by `verify_articulated_registration.py`, anatomical landmark residuals, shared-pivot motion tests, and selected mesh-to-mesh joint clearances under small prescribed joint rotations. Include knee translation laws explicitly. Set anatomical tolerances from the landmark uncertainty and intended task, not from the existing 62 mm fit. Add no extra inertial mass.

No registration correction is applied by this audit. True articular landmarks and surface clearances remain unmeasured; extracting a small selected set of retained bone surfaces is the next bounded geometry task.
