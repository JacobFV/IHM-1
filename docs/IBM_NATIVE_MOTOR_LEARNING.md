# Reduced IBM kernel learns native ankle feedback

`scripts/train_native_motor_learning.py` copies eight uniformly selected embedding
rows from IBM's fused association checkpoint and trains only those rows to imitate
an engineering ankle PD controller. The sensory encoder and signed motor readout
are fixed and have no bias. Position error and velocity must cross the normalized
embedding/tanh association matrix to produce muscle excitation. Severing this
matrix produces exactly zero excitation regardless of target or feedback.

This is a deliberately reduced motor materialization. It does **not** execute
IBM's cortical E/I dynamics, anatomically locate its eight selected sites, train
the full fused kernel, or use the native peripheral/cord pathway. The externally
provided desired ankle angle is a task instruction, and native coordinate/speed
are privileged feedback. The full native articulated plant receives only muscle
excitations; the evaluator never prescribes coordinates or adds external forces.
The free environment has gravity and no contact support; this is isolated joint
control in a falling body, not standing or walking.

Reproduce from IHM root (a fresh output directory is required):

```sh
.venv/bin/python scripts/train_native_motor_learning.py \
  --output data/runtime/motor-learning/ankle-new --seconds .3
.venv/bin/python scripts/test_native_motor_learning.py
```

The evaluation restores an identical complete native checkpoint for every arm:
trained kernel, severed kernel, untrained copied kernel, no controller, and the
engineering teacher. Three held-out desired angles are compared with the same
tracking MSE objective. Training uses synthetic position-error/velocity examples,
not these native trajectories. A separate `motor_kernel.pt` preserves the original
checkpoint hash, selected rows and training definition. `receipt.json` records
trajectories, source hashes and all arms, including failures of the improvement
criteria. The original IBM checkpoint is never written.

The artifact can be loaded by `ihm.native.motor_learning.load_motor_policy`.
It cannot be loaded directly into `IBMImplicitController`: the materializations
have different ports and dynamics. Treat this as a narrow learned motor primitive
and an explicit causal training experiment, not a claim that the preexisting
fused perceptual kernel had motor competence. Matching the teacher demonstrates
successful imitation; it does not establish biological learning or a benefit
from the checkpoint's prior training over random initialization.

## Recorded evaluation

`data/runtime/motor-learning/ankle-20260908/receipt.json` records 15 native
rollouts, each 0.30 seconds, using target offsets -0.12, +0.12 and +0.22 rad.
Mean tracking MSE in rad² is 0.01085397 for the trained kernel, 0.17261620
for severed/no-controller, 0.44361657 for untrained, and 0.01085396 for
the engineering teacher. All three target directions/values improve against
severed and untrained. The 93.7% aggregate improvement establishes a causal role
for the trained primitive in this experiment. It does not establish a useful
role for prior perceptual training; the trained primitive reproduces its teacher.

## Unified runtime adapter

`ihm.native.motor_learning_controller.IBMAnklePrimitiveController` implements
the unified runtime's controller protocol. It advances the persistent IBM E/I
controller, then explicitly selects the trained ankle primitive as motor owner.
Its metadata says `cortical_motor_output_active=false`; inactive cortical
commands and cord arc statistics remain separately named. The primitive's motor
output does not pass through the cord and has zero modeled neural delay.
Blocking either ankle sensory channel zeros this privileged joint feedback;
motor blocks zero the corresponding output. Oxygen/perfusion availability scales
primitive excitation. Checkpoints bind the original controller state, learned
artifact hash, target and primitive commands. Retained source copies include the
separate learned artifact and both Python adapter files.

The default artifact is the recorded local training output above; the default
target is +0.12 rad. Factory/API exposure is a separate integration change. This
adapter makes the primitive usable in the same native world/physiology runtime,
but the standalone 15-arm experiment does not certify that longer unified run.

## Actual cortical E/I follow-up

`scripts/evaluate_native_cortical_motor.py` evaluates the separate 128-site
`cortical_motor.pt` artifact created by `scripts/train_cortical_motor.py`.
This follow-up runs the actual IBM E/I equations with disjoint sensory and motor
ports and persistent cortical state throughout each native rollout. The trained
parameters are copied association embeddings and a bias-free motor decoder.
Severing zeros association edges while retaining the decoder and local E/I
dynamics. Three desired native ankle targets each have full, severed and
zero-controller arms, with identical restored native starting state. Both the
native mechanics and cortex advance 20 ms per exchange, cortex in 1 ms ticks.
This evaluation is distinct from the reduced eight-site experiment above.

```sh
.venv/bin/python scripts/evaluate_native_cortical_motor.py \
  --output data/runtime/motor-learning/cortical-native-new
```

The task instruction and feedback remain privileged desired ankle angle and
native joint angle/speed. Training is engineering PD imitation, not biological
learning. It does not add a gait, balance controller, or anatomical sensory path.

The completed follow-up receipt at
`data/runtime/motor-learning/cortical-native-20260908/receipt.json` records nine
0.30-second rollouts. Mean tracking MSE is **0.01996592 rad² full** versus
**0.17663703 rad² severed/no-controller**, an 88.7% reduction. Each held-out target
improves separately: -0.12 rad (0.02321069 versus 0.04545829), +0.12 rad
(0.01032972 versus 0.19423190), and +0.22 rad (0.02635736 versus 0.29022091).
This establishes a contribution through the trained actual cortical association
pathway on this limited task. It does not establish prior perceptual pretraining
benefit, long-horizon control, or anatomical/biological fidelity. The reduced
experiment samples at 10 ms and this cortical experiment at 20 ms, so their
aggregate MSEs should not be compared as a matched controller benchmark.
