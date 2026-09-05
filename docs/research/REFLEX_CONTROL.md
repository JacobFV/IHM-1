# Source-grounded stretch feedback in the canonical mechanics reduction

2026-09-05. `BodyReflex` and the neuromechanical experiment execute a mechanics → length sensor → delayed spinal controller → muscle activation → mechanics loop. The tested effector is the retained right tibialis-anterior actuator, `body-muscle-opensim-tibant_r`, attached to canonical tibia/calcaneus identities. This is a reduced stretch-feedback experiment, not a validated withdrawal reflex, articulated ankle, autonomous gait, or cortical motor controller.

## Controller source and transfer

[Geyer and Herr (2010)](https://dam-prod.media.mit.edu/x/files/wp-content/uploads/sites/3/2013/04/A-Muscle-Reflex-Model-that-Encodes-Princlples-of-Legged-Mechanics-Produces-Human-Walking-Dynamics-and-Muscle-Activities.pdf), Table I and Appendices I–II, supplies the TA positive-length primitive, gain 1.1, length offset 0.71, prestimulation/minimum stimulation 0.01, distal delay 20 ms and activation time constant 10 ms. These are model parameters, not individual measurements. The implementation selects the uninhibited TA length primitive; it does not reproduce their complete state-dependent locomotor controller or soleus inhibition. Normalized length-feedback primitives and a first-order activation equation are also explicit in [Wang et al. (2012), §4.1 and §3.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4523558/).

For normalized sensed length `L`, external descending drive `d`, and external gain multiplier `g`, the implemented excitation is

```text
u(t) = clip(0.01 + d + g × 1.1 × (L(t − 0.020 s) − 0.71), 0.01, 1)
da/dt = (u − a) / 0.010 s
```

The numeric defaults retain the published identity while their application to this mechanical reduction is a transfer assumption. The newly retained SCONE `ControllerGH2010v9.scone` independently lists TA `KL=1.1`, `L0=0.71`, delay 0.020 as paper values. Its separately optimized H0914 results are not substituted into this experiment.

The inspected paper is retained at `data/derived/audits/neuromechanical/sources/geyer_herr_2010.pdf`, SHA-256 `71619c24b0c7f2f4a0c2eca2149e6b1a17542487a82a558f3bff140327e7d381`. The builder fetches it from the cited author-hosted location only when absent and requires that digest.

## State, sensor and ownership

`ihm/assembly/reflexes.py` owns the event queue, arrived sensory feature, activation and clock. Its input is a computed mechanical state with explicit time, muscle tension and every path attachment's translation, rotation and deformation. Missing transforms or forces fail; identity transforms and zero force are not silently substituted. Each attachment is transformed about its reference centroid, exactly as in `BodyMechanics`; the sensor then sums the resulting polyline segments.

The existing reduced mechanics does not solve native contractile-element length. The declared observation law is

```text
fiber_length_proxy = optimal_fiber_length
                   + (current_path_length − reference_path_length) / cos(pennation)
L = fiber_length_proxy / optimal_fiber_length
```

This fixes tendon length and pennation and aligns the reference path with optimal fiber length. It does not use the unrelated affine muscle-belly contraction as another length sensor, infer spindle firing in Hz, or claim physiological spindle calibration. The tension port is observed and recorded but does not drive this selected length-only controller.

At time `t`, the sensor sample is queued for `t + loop_delay_s`. The activation ODE is integrated analytically between sample arrivals. The caller applies the returned endpoint activation during its next mechanics interval. Consequently, the experiment adds one explicit 50 µs exchange interval to the lumped neural loop. Prehistory has no arrived length deviation and activation starts at zero; the specified prestimulation still acts. The 20 ms delay is not divided into invented separate afferent, synaptic and efferent measurements.

Sensory block clears in-flight/arrived feedback and retains the explicit basal/descending drive. Motor block clears the queue and sets excitation to zero; existing activation decays with its own pole. Controller removal suppresses the feedback term while retaining basal/descending inputs. Unblocking requires new sensory propagation. The example uses the retained anatomical-prior nerve label `peripheral-nerve-right-deep_fibular`; that label is not a reconstructed axonal route.

`checkpoint()` includes model digest, clock, activation, arrived feature, serial and every queued sample. `restore()` verifies identity and finite ordered state before mutation. Invalid step inputs also leave the controller unchanged. The experiment separately saves and restores mechanics state and the exchanged activation, so its continuation check covers the entire loop.

## Mechanical experiment

Run from the repository:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.verify_body_reflexes
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_neuromechanical_experiments
```

The builder verifies `mechanics.json`'s retained source-file hashes. It extracts the existing right tibia, talus, calcaneus and TA tissue carrier, the complete registered TA polyline, and the existing talus–calcaneus support link. No new stiffness, mass or insertion is fitted. Proximal bones and the TA carrier have prescribed zero translation; calcaneus translation is solved by the existing mechanics engine. Orientations remain constrained. The source support is an inferred geometric elastic/damping prior, and the muscle law is the existing Gaussian active force-length / linear passive-tension reduction, not OpenSim's equilibrium muscle.

The extracted experiment is an independent mechanical materialization, not additional mass added to a simultaneous whole-body ledger. Initial conditions are the source's unloaded, zero-activation reference without gravity. A 20 N external force along the final TA path segment acts from 0.10 to 0.25 s. It changes free-body motion; no subsequent path or position trajectory is prescribed. The total duration is 0.40 s.

The intact and no-load cases isolate mechanical input to the sensor and its downstream motor consequence. Matched cases apply motor block, sensory block, controller removal or half descending gain at 0.10 s. A recovery case releases motor block at 0.25 s. The delayed case instead changes the loop-delay parameter from 20 to 60 ms from the beginning; it is a parameter comparison, not a matched-history intervention at 0.10 s.

## Verification and numerical evidence

The latest retained run is `data/derived/audits/neuromechanical/verification/run-z8_c__4w/`. It includes complete frame traces, controller checkpoints, source receipts and a report. `reflex_response.png` and `.svg` show the oscillatory path/activation/tension responses.

| Check or comparison | Result |
|---|---:|
| First load-induced activation difference | 0.1282 s, after load at 0.10 s and the 20 ms neural delay |
| Intact versus no load, maximum path difference after load | 7.608 mm |
| Intact versus motor block | 40.239 mm |
| Intact versus sensory block/controller removal | 37.931 mm |
| Intact versus half descending gain | 20.570 mm |
| Intact versus 60 ms loop delay | 34.918 mm |
| 50 versus 25 µs integration, maximum path difference | 0.9853 mm |
| Maximum internal force residual, intact | 1.15e-13 N |
| Maximum internal torque residual, intact | 2.93e-14 N m |
| Full controller + mechanical checkpoint continuation | Exact frame equality |

Unit checks independently verify the delayed response against the analytic activation solution, constant-input time partitioning, block decay and release latency, descending gain, rigid-motion invariance of path sensing, missing-input rejection, cross-model checkpoint rejection and atomic invalid-input behavior. Mechanical checks require finite state, bounded activation, pre-intervention frame equality for action-switch comparisons, restored activity after nerve release, and actual motion differences.

The original coarse refinement failed and remains retained: 1 versus 0.5 ms differed by 33.87 mm over a 1.2 s trial; 100 versus 50 µs differed by 1.938 mm over 0.4 s. The source gain produces strong oscillations when transferred to this small inferred support model. The final 50/25 µs comparison passes the declared 1 mm numerical gate over 0.4 s; it does not establish asymptotic stability or human accuracy. The published gains were not silently retuned to make this gate pass.

At 0.4 s, intact accumulated external work is 0.5106 J and active mechanical work is 6.8824 J. The mechanics energy residual is −0.3644 J, retaining the existing discrete work quadrature and implicit numerical dissipation; exact global energy closure is not claimed. Active mechanical work is derived from the actual mechanics work increment after subtracting measured external and prescribed-boundary work. It is not a metabolic cost estimate.

## Integration boundary and remaining gaps

The component exposes bounded `motor_activations` compatible with the existing mechanics input. It is not yet installed in shared `BodyRuntime`, whose current peripheral path remains separate; an integrating coordinator must choose one motor-delay/activation owner rather than stacking this component on `BodyPeripheral`'s activation pole. The controller's normalized sensory feature is not mapped to IBM regional Hz without a supported observation law. No IBM microcircuit acquisition or donor edit is involved.

Native BioGears brainstem chemoreflex and autonomic control remain independent owners. No mechanical-work-to-native-exercise/metabolic mapping is implemented here: measured mechanical work alone cannot determine activation heat, maintenance cost or efficiency. A later source-supported metabolic law and native demand port are necessary for that additional causal edge. The verified result is a spinal stretch-feedback primitive producing computed reduced muscle-driven motion, with explicit source-transfer and numerical limits.
