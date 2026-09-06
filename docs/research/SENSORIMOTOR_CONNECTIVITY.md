# Reusable sensorimotor exchange

`SensorimotorController` now closes computed muscle feedback through delayed spinal pathways and an explicit IBM population decoder to named motor excitation ports. This is a reusable controller, not a walking demonstration or a calibrated voluntary motor policy. Native mechanics remains the sole owner of activation, tendon/fiber state, muscle force and mechanical work.

## Exact retained capabilities and evidence

`BodyBrain` verifies and executes the preserved Wilson–Cowan and shunting-inhibition functions in `data/derived/canonical/brain-sources/ibm-neural.py` (SHA256 `f8c7d3172af71b399d9e6301f068d8fde8adf3363941ef0d3d3ec50ef62f663b`). Eighty registered brain regions evolve through the existing population ODE. The underlying cortical anatomy/association weights and physiology bridges have their existing population-transfer limitations. No live IBM checkout is imported or edited.

The retained IBM tactile transfer functions instead return signed donor-response units. They do not identify receptor firing rates or muscle recruitment. This controller therefore does **not** silently relabel their outputs as Hz or muscle force. It uses explicitly parameterized strain/force-to-population-rate bridges. Identifying a tactile observation law remains necessary before the signed tactile response can replace these inputs.

The [Geyer and Herr primary paper](https://www.cs.cmu.edu/~cga/tmp-public/song.pdf), DOI [10.1109/TNSRE.2010.2047592](https://doi.org/10.1109/TNSRE.2010.2047592), was acquired and retained at `data/raw/sensorimotor/geyer_herr_2010.pdf`. Table I and Appendix I supply normalized soleus force gain 1.2, gastrocnemius force gain 1.1, tibialis length gain 1.1 and offset .71, reciprocal soleus-to-tibialis gain .3, basal excitation .01, and the 20 ms distal neural delay. These are published model parameters tuned for that paper's walking plant, **not measurements calibrated to the held OpenSim subject or a generic body**. The receipt is `data/sources/sensorimotor/card.json`; `from_root` verifies the retained PDF hash.

## Coverage and explicit transfer assumptions

`native_muscle_catalog(root)` reads all 80 muscles, original attachment body names, maximum force and optimal fiber length from the retained `example3DWalking/subject_walk_scaled.osim`. These are bilateral **lower-limb** effectors; the source has no active upper-limb muscle set. Attachment names provide hip/pelvis, knee and ankle/foot group labels. Contralateral pre/postcentral assignment is a regional engineering prior, not measured fine somatotopy. Every catalog row retains the source path and SHA.

All 80 muscles provide normalized fiber length and tendon force feedback and accept individually addressed descending drive. Only eight effectors have spinal primitives: bilateral tibialis anterior, soleus and both gastrocnemius heads. The two gastrocnemius heads share a lumped force sensor, computed as summed force divided by summed maximum force, and common source reflex gain. This split-to-lumped transfer is an explicit assumption. Other muscles receive zero excitation unless specifically requested through their descending port.

Actual foot normal force gates stance force feedback and soleus inhibition. The 5 N threshold is an engineering prior. Tibialis length feedback continues in swing, while soleus/gastrocnemius retain source basal stimulation. This is only the distal subset of the published controller; it omits hip/trunk, knee-protection, swing-placement and balance laws. None of these missing laws is inferred from anatomy labels.

Measured muscle strain/force is converted to cortical excess rate through declared coefficients (200 Hz per positive normalized fiber extension, 40 Hz per normalized tendon force). A requested effector drive enters its precentral population (100 Hz per request fraction). Evolved precentral activity changes reflex gain by .005 per Hz relative to its initial rate and provides a .002 per Hz drive, gated by that exact effector's request. These coefficients are **uncalibrated engineered decoder parameters**; biological motor recruitment and learned intentions are unresolved. Sensory effects on the motor population traverse the existing IBM network dynamics, rather than bypassing the brain under a brain label.

## Runtime contract

```python
c = SensorimotorController.from_root(root)
result = c.step(
    0.02,
    {"time_s": t, "muscles": muscles,
     "foot_contact_force_n": {"r": right_normal_force, "l": left_normal_force}},
    descending={"addbrev_r": 0.2},
    sensory_blocks=["soleus_r"], motor_blocks=[],
)
next_commands = result["motor_excitations"]
```

Each muscle observation requires `fiber_length_m` (actual CE) **or** explicitly named `fiber_length_proxy_m`, `optimal_fiber_length_m`, `tendon_force_n`, `max_isometric_force_n`, and `sensor_basis`. Missing channels, unknown targets, nonfinite values and clock mismatch are rejected. Input dictionaries are never mutated. The plant supplies actual geometry/path, force and contact observations; a proxy is never silently promoted to native fiber state.

The controller splits the source lumped 20 ms neural latency into 10 ms afferent and 10 ms efferent delays; that split is unmeasured. At time t, it samples the current plant, advances the brain using previously arrived afferents, delivers queued events, and schedules computed excitation after the motor delay. Returned excitation at t+dt applies to the **next** plant interval. Exchange discretization adds latency beyond transport delay; the returned `exchange_interval_s` exposes it. There is no interpolation of unavailable future observations. The parent runtime must advance mechanics using the previous command, then exchange the new endpoint observation on the next tick.

Sensory blocks clear current and in-flight afferents for selected muscles. Motor blocks immediately clear that output and its queued motor commands; actual activation decays in the plant. A sensory block leaves descending and basal drive intact. Unblocking starts new causal transport. The controller supplies no second activation filter.

`checkpoint/restore` includes all delayed events, delivered sensory features, excitation state and the complete brain population state. Identity binds the source catalog, brain source descriptors, decoder parameters and executed controller/brain module hashes. Restoration validates shape, finiteness, keys, event order, clock and model identity before committing; step failures restore the entire neural state. The parent must include this checkpoint in its plant/native transaction.

## Bounded verification

Run single-threaded:

```
OPENBLAS_NUM_THREADS=1 nice -n 10 .venv/bin/python -m scripts.verify_sensorimotor
OPENBLAS_NUM_THREADS=1 nice -n 10 .venv/bin/python -m scripts.probe_sensorimotor_feedback
```

Four tests cover delayed output, unilateral sensory/motor blocks, eight source-law values, all-80 catalog mapping and effector-specific non-ankle commands, actual IBM state changes, immutable observations, checkpoint replay and rollback after invalid neural input.

The separate 1 kg linear mechanical fixture integrates position, velocity and actuator activation; its computed length and force feed back into the controller. At approximately 0.10 s, active forces are -4.099 N intact, -0.100 N with sensory block and zero with motor block. Releasing the motor block restores force. Halving exchange time from 2 ms to 1 ms changes final position by 0.172 mm. This is a numerical causality/refinement check on an explicitly synthetic test plant, not anatomical, native OpenSim, gait, or human validation. Frames and model identity are retained in `data/derived/audits/sensorimotor_feedback_v1/report.json`. Actual native integration belongs to the coordinated mechanical/runtime acceptance tests.
