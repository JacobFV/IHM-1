# Actual E/I cortical stance imitation: research result

This experiment connects a copied 256-site IBM E/I cortex to all 98 native muscle commands through 258 explicitly named native state inputs. It is **not accepted as a stable standing controller** and is not enabled as a runtime mode.

The teacher is the 10 ms sampled native margin-aware LQR artifact at `data/research/locomotion_control/linearization_eyxmbi_z/discrete_margin/linearization.npz`. The 70 kg native model identity is `c3a936233b06e694340504ee14e5b3e2ae8687e093b65cc185aba31d5fa5a29a`. Inputs come from complete retained native snapshots in `lqr98_margin_push/trajectory.json`; labels are recomputed from each snapshot, because recorded commands belong to the preceding mechanical interval. No 86-muscle data is padded into this 98-muscle model.

A fixed random projection maps the complete normalized native state into sensory cortical sites. It does not contain the teacher gain matrix or teacher commands. Motor sites are disjoint. Only the copied association embedding and cortical motor decoder are trained. The decoder adds a learned correction to the explicitly retained equilibrium excitation `u0`; severing the association preserves exactly that same tonic baseline. Excitations respect the native 0.01 activation floor.

## What failed

The first margin-teacher imitation achieved small held-out command errors but introduced a nominal excitation error of 0.000625. In actual matched native trials the trained controller fell by **1.12 s**, earlier than the severed baseline at **1.76 s**. This result is retained at `data/runtime/motor-learning/cortical-stance-native-margin-20260908/report.json`. Imitation MSE was not evidence of balance.

## Explicit neural baseline normalization

The revised experiment maintains two persistent trajectories of the **same IBM E/I equations and kernel**: one receiving the actual native feedback, and a counterfactual receiving zero state error. The decoder reads their motor-activity difference. This is an engineering normalization of intrinsic neural activity, not a biological baseline circuit or an additional mechanical controller. It provides no PD shortcut. At zero state error it produces exactly the tonic command baseline throughout the measured one-second drift check.

Training the normalized version took 120 bounded epochs. On an interleaved holdout from the same native trajectory:

| Command comparison | MSE |
|---|---:|
| Full trained cortex | 6.90e-8 |
| Severed association, same baseline | 1.44e-6 |
| Original embedding, trained decoder retained | 1.17e-6 |

These are correlated within-trajectory examples, not independent-trial generalization. In a separate two-second persistent replay of recorded states, full MSE was 7.41e-7 versus 1.99e-6 for sever. A real recorded perturbation held as a step input produced 10%, 50%, and 90% of its 100 ms cortical response at **18, 37, and 63 ms**. That response time is distinct from the teacher's 10 ms sampling interval.

## Actual native comparison

Both arms restored the same native mechanical checkpoint. A correctly located 5 N forward plus 5 N lateral force acted at the pelvis center of mass from 1.0 to 1.1 s. No stabilizer, coordinate prescription, support, or direct teacher feedback bypass was added.

| Arm | Outcome | Peak COM displacement |
|---|---|---:|
| Trained normalized cortex | Completes 2 s; pelvis height 1.004 m | 5.95 cm |
| Severed cortex, same tonic baseline | Falls by 1.76 s | 37.48 cm |

The full arm already drifted 2.64 mm before the push and ended with COM speed **0.168 m/s**. It improved the short mechanical outcome but did **not** demonstrate stable recovery. Longer horizons, broader perturbation data, and compensation for cortical response time remain necessary before acceptance.

The earlier teacher trajectories used a force-location convention subsequently found to pass local pelvis coordinates to a ground-coordinate force port. They remain real native perturbations, but are not labeled as center-of-mass pushes here. The paired cortical native test above converts the local center-of-mass station to ground coordinates correctly.

Artifacts:

- Trained normalized model: `data/runtime/motor-learning/cortical-stance-reference-20260908/cortical_stance.pt`
- Contract, baseline, persistence and response-time checks: sibling `verification.json`
- Paired native result: `data/runtime/motor-learning/cortical-stance-native-reference-20260908/report.json`
- Negative earlier result: `data/runtime/motor-learning/cortical-stance-native-margin-20260908/report.json`

Entry points are `scripts/train_cortical_stance.py`, `scripts/verify_cortical_stance.py`, and `scripts/evaluate_native_cortical_stance.py`. These are native-mechanics experiments without the full physiological exchange and make no walking claim.

## Patient-specific temporal training and port audit

The later patient model uses 77.6122029 kg, 258 exact native states and 98 muscle commands. The slower margin-aware teacher (R multiplied by 1000) independently recovered a native pelvis-COM push over 10 s even with a four-sample command queue. That is a result for the engineering teacher, not for the cortex.

A 512-site temporal cortical imitation still fell at 2.46 s in the held-out negative diagonal push; its matched severed baseline fell at 1.94 s. The failed trial is retained at `data/runtime/motor-learning/cortical-stance-native-patient-slow-negative-20260908/report.json`.

A subsequent port audit found the fixed random sensory encoder had only 107 independent channels for 258 state variables. Of the normalized teacher gain's Frobenius norm, 76.98% lay in the encoder nullspace. Projecting the exact slow teacher onto those observable directions made its delay-free local closed-loop spectral radius 1.7441, with 10 unstable modes. Projection onto the motor decoder's column space alone remained locally stable. These are linear architecture diagnostics, not a nonlinear cortical controller proof; recurrent history could in principle support an explicitly trained observer.

The next research architecture therefore assigns 516 signed identity sensory ports and 196 disjoint motor ports in an actual 1024-site E/I materialization. This is an explicit engineering port allocation without an anatomical-region claim. Its sensory matrix has rank 258 and normalized teacher-nullspace error 3.33e-16. It removes the measured instantaneous information loss; it does not establish native balance. The exact encoder matrices, materialized kernel, equations, and source provenance are retained with each training artifact. The original shared IBM checkpoint remains unchanged.

## Follow-up with completed native recovery

The full-state small-signal controller subsequently recovered both opposite pelvis-COM pushes for ten seconds, with matched severed arms falling. The complete final architecture, engineering assumptions, failed intermediate trials and measured outcomes are documented in [IBM_CORTICAL_STANCE.md](IBM_CORTICAL_STANCE.md). The native-tested bundle is `data/models/ibm_cortical_stance_v1`; this does not retroactively validate the pilot artifacts above.
