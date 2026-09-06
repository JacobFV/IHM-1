# Retained lumbar and shoulder muscle coverage

The current 92-muscle source model has no lumbar-crossing muscle path. Its 80 fitted paths reference only lower-limb coordinates. The other 12 are the bilateral Arm26 extension: **six cross the shoulder**, namely TRIlong, BIClong and BICshort on each side; TRIlat, TRImed and BRA are elbow-only. Shoulder muscle action is therefore incomplete, not absent. A path crossing does not certify useful torque in every shoulder direction, particularly with unilateral tension and wrapping. This is a source audit, not a supported-pose result.

Reproduce without OpenSim or geometry decoding:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/audit_lumbar_shoulder_coverage.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_lumbar_shoulder_coverage.py
```

The generated `data/research/lumbar_shoulder_coverage/audit.json` hashes the exact immutable bytes parsed from five retained inputs, lists all 80 fitted coordinate dependencies, extracts all 12 arm paths and wrapping references, and retains 38 relevant canonical muscle entities with their geometry pointers, anchors and parameter evidence. `lumbar_candidate_forces.xml` is a six-muscle ForceSet insertion fragment, not an activated model. Neither script imports OpenSim, loads a mesh, changes a native model, or assigns excitation.

## Concrete lumbar candidate

The retained Gait2392 Thelen2003 model contains bilateral `ercspn`, `intobl`, and `extobl`, each a two-point pelvis→torso path with no wraps. Its embedded credits identify CC BY 3.0, Delp and colleagues, and the Anderson/Pandy trunk model lineage. Exact source bytes and model parameter values are retained by hash; this audit does not claim those population model values were measured on the current subject.

| Donor muscle, each side | Fmax (N) | Optimal fiber (m) | Slack tendon (m) | Neutral path (m) | Right neutral moment arms, extension/bending/rotation (mm) |
|---|---:|---:|---:|---:|---|
| Erector spinae | 2500 | .12 | .03 | .149708 | 42.690, 38.089, 9.691 |
| Internal oblique | 900 | .10 | .10 | .199096 | −52.824, 93.847, −36.106 |
| External oblique | 900 | .12 | .14 | .274059 | −62.791, 62.253, 26.005 |

Left bending/rotation signs reverse; extension signs remain. The neutral 3×6 moment-arm matrix has rank three. This is local geometric rank, not proof of a physiologically feasible torque cone, force-length capacity at a displaced posture, or equilibrium. For example, an isometric Fmax multiplication omits fiber/tendon state and activation and must not be used as an available support force.

Registration uses the homologous `back` joint, whose donor and target rotation sequences/axes match and whose body offset orientations are zero. For each attached body, `target_station = target_back_offset + donor_station − donor_back_offset`. Unit scale is an explicit conservative transfer assumption. This preserves donor path lengths and moment arms at neutral while positioning the donor joint at the existing target joint. Pelvis offset changes from (−.1007,.0815,0) to (−.1158377177,.0921303039,0) m; torso offset remains zero. It is a body-local attachment transfer, not a new segment-specific world transform. Target morphology and insertions still need independent anatomical review; copying target scale factors from an unrelated donor is not justified. The fragment preserves the complete Thelen parameter elements and renames GeometryPath to the modern `path` property as required by the current target document version.

The verifier checks deterministic extraction, exact retained fragment bytes, donor-relative station preservation, and each neutral moment arm against an independent central finite difference of path length. It also checks equal-and-opposite endpoint force torque cancellation. No loading result is claimed.

## Shoulder and canonical evidence limits

Missing named native actuators include the deltoid parts, pectoralis major parts, supraspinatus, infraspinatus, subscapularis and teres muscles present in the retained canonical inventory. Their canonical rows are **inferred attachments**, generally extrema projected to nearby bones. Their strengths use `PCSA = estimated volume / assumed fiber length`, 0.3 MPa specific tension, assumed zero pennation and a 10% path slack tendon. These are engineering priors, not retained experimental origin/insertion or wrapping measurements. All exact candidate rows are exported so these limitations cannot disappear in a manual transcription.

The current native model has no scapula or clavicle bodies or scapulothoracic coordinates. Attaching scapular canonical points directly to torso would explicitly freeze scapular motion. It can only represent an acknowledged fixed-girdle approximation; pectoralis minor and other girdle-only paths would then produce no useful shoulder-joint torque. Generic straight lines from canonical endpoints can also pass through humerus or joint surfaces. The held Arm26 donor provides biceps/triceps/brachialis paths only and cannot supply missing deltoid/rotator-cuff parameters by renaming its muscles. No source-supported complete shoulder extension is established by these files.

## Safe next integration boundary

1. Insert only the six named Gait2392 lumbar forces into a separately hashed opt-in model variant, keeping all 22 existing bodies and their exclusive inertias. These massless force elements represent tissue already included in body mass; add no tissue mass. Preserve the original 92 muscle order and append six explicit source-bound catalog entries. Do not silently replace all ForceSet contents with the fragment.
2. Retain GeometryPath evaluation for the six additions. The existing 80-path fitted set contains no lumbar functions and must neither overwrite nor disable these paths. Confirm the adapter can accept 98 muscles, including control ordering, metabolic ownership, signed muscle ledgers, checkpoint state, and factory receipts.
3. In an authorized native validation slot, load the separate variant, compare native neutral moment arms against this artifact, scan physiological lumbar configurations for path continuity and force-length/tendon feasibility, and verify action–reaction, virtual work, passive state, and replay. Report reserve torque separately. No excitation chosen just to settle the support fixture is a biological control policy.
4. Extend brain/motor ownership with explicit lumbar drive bindings and provenance only after the force model is verified. Existing lumbar CoordinateActuators remain independently identified reserves; their torque is not muscle excitation. A later support solve must report whether legitimate muscle states, contact and reserves balance, rather than concealing residuals in a new actuator.
5. For shoulders, acquire a licensed complete shoulder donor with origin/insertion and wrapping evidence, register humerus plus scapula/clavicle landmarks to the retained anatomy, and decide explicitly whether to partition actual girdle bodies or declare a fixed-girdle approximation. Validate its tension-feasible torque cone before claiming shoulder coverage. Keep these roles separate from the actor-recovery thoracic/rib/diaphragm lane, which currently contributes geometry and no calibrated trunk muscle forces.
