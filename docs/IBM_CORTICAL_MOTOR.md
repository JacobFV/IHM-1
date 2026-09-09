# Trained actual IBM cortical ankle experiment

A copied 128-site IBM `CorticalDynamics` E/I materialization now makes a measurable contribution to an isolated native ankle task. This is a short motor-control experiment, not walking, balance, or a validated human brain.

The fixed sensory encoding receives an explicit target error and native ankle angular velocity. Its sites are disjoint from the precentral decoder sites. Training updates the copied association embedding and a bias-free decoder through the actual E/I equations at 1 ms. Cortical state persists between runtime exchanges. The original 30,000-site fused checkpoint is unchanged.

Training uses synthetic engineering supervision `clip(3 * error - 0.22 * speed, -0.45, 0.45)`, not human motor demonstrations. In 256 held-out synthetic examples:

| Arm | Command MSE |
|---|---:|
| Trained cortex | 0.00075446 |
| Association severed; decoder retained | 0.04474818 |
| Embedding reset to source; trained decoder retained | 0.04300785 |

Signed target/output correlation is 0.9923. Severed command variation is zero. Resetting only the embedding demonstrates that this result depends on training the core, rather than fitting only a readout.

Three matched 0.3-second native articulated-body trials per arm then tested previously held-out target angles. State remains persistent throughout each trial.

| Target (rad) | Trained cortical tracking MSE (rad²) | Severed MSE (rad²) |
|---|---:|---:|
| -0.12 | 0.0232107 | 0.0454583 |
| +0.12 | 0.0103297 | 0.1942319 |
| +0.22 | 0.0263574 | 0.2902209 |

These native trials use a **free-falling articulated body**, privileged joint feedback, and no downstream cord. They demonstrate a contribution to short ankle tracking, not sustained standing or locomotion. Zero-controller controls are retained in the receipt. The full native physiology/world controller mode must be assessed separately.

The runtime adapter `IBMCorticalAnkleController` (`implicit_cortical_ankle`) runs this trained E/I cortex as the sole cortical motor owner. An optional segmental cord follows it. The existing `implicit_ankle_primitive` mode uses a different eight-site algebraic kernel; its results are separate. Raw `implicit` mode remains untrained.

The trained adapter preserves transactional cortical/cord state, explicit muscle mapping, motor/sensory blocks, and per-arc telemetry. MAP/oxygen availability and temperature scale the E/I dynamics through declared IHM engineering priors. Additional regional sensory rates provide a bounded supplemental sensory drive; that transfer is not a learned anatomical encoding.

Reproduction:

```sh
.venv/bin/python scripts/train_cortical_motor.py --epochs 150
.venv/bin/python scripts/verify_cortical_motor_controller.py
.venv/bin/python scripts/evaluate_native_cortical_motor.py --help
```

Artifacts and evidence:

- `data/runtime/motor-learning/cortical-20260908/cortical_motor.pt`
- `data/runtime/motor-learning/cortical-20260908/training.json`
- `data/runtime/motor-learning/cortical-20260908/embedding_ablation.json`
- `data/runtime/motor-learning/cortical-native-20260908/receipt.json`

The trained artifact SHA256 is `49557ac8562db3811e0004d4489900c6301023cb394df89a0ed767912a5e8164`. Retain its accompanying exact `pretrain_video_loop.py`; the loader rejects an equation-source hash mismatch.
