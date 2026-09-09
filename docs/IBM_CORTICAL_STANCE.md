# Cortical stance: bounded native recovery with a causal kernel ablation

The promoted `data/models/ibm_cortical_stance_v1` controller uses persistent IBM E/I dynamics to correct all 98 muscle commands in the identified 77.6122029 kg patient model. It recovered two opposite small pelvis pushes for ten seconds. Severing its association kernel retains the same tonic muscle baseline and causes a fall in both matched trials.

| Pelvis COM force, applied for 100 ms at t = 1 s | Full cortex peak COM | Full cortex final COM at 10 s | Severed outcome |
|---|---:|---:|---|
| (+5, 0, +5) N | 1.784 mm | 0.196 mm | Falls at 2.90 s |
| (−5, 0, −5) N | 2.264 mm | 0.172 mm | Falls at 1.94 s |

The full arms finish with COM speeds of 3.61 and 5.20 µm/s respectively. Both arms restore the same native checkpoint within each experiment. Forces are applied at the actual ground-space pelvis center of mass. No coordinate prescription, external support, live LQR, or cord motor bypass stabilizes the full arm. These are native mechanics experiments at baseline neural physiology, not proof of walking or whole-world physiological recovery.

## What runs

- 1,024 materialized sites of the actual IBM `CorticalDynamics` E/I equations, integrating at 1 ms and exchanging commands with native mechanics at 10 ms.
- 258 privileged native state variables, explicitly including joint value/speed and every identified muscle activation/fiber state. A signed identity encoder supplies 516 disjoint sensory sites; 508 other sites supply motor activity to a bias-free decoder.
- A copied association embedding trained by temporal imitation, followed by an offline decoder fit against actual local E/I responses. The engineering slow-LQR teacher provides training targets only. Its gain is absent from the runtime motor path.
- A fixed native equilibrium command baseline plus cortical corrections. A persistent zero-input trajectory of the same kernel provides neural baseline subtraction. This counterfactual reference is an engineering normalization, not an anatomical brain circuit.
- Explicit small-signal scaling: native input scales and decoder gain are each multiplied by 100, preserving local DC feedback gain while reducing nonlinear cortical drive. Float64 computation and exact native equilibrium values avoid amplified subtraction noise. No E/I equation was replaced.

The original shared IBM checkpoint remains unchanged. The promoted artifact SHA-256 is `0747ac07a4baa02d30964933c1f1fd1bf1285d9a30a9c35e9eac7de2e2c98042`. Its manifest retains lineage and both paired native receipts.

## Why the earlier versions failed

The initial 512-site implementation compressed 258 states into a rank-107 random sensory projection. That discarded directions containing 76.98% of the normalized teacher gain. Simply adding full-rank sensory ports did not solve control: its learned neural/native sample eigenvalue was 1.3501 and it diverged before the push. The negative native receipts are retained.

Expanding the motor readout to 508 sites yielded a rank-258 neural transfer with condition number 18.86. An offline local decoder fit reduced static gain error to 4.52e-8 and removed the measured unstable augmented modes. However, the initial finite push drove that cortex outside its fitted local regime; it fell at 1.85 s, before its severed baseline at 1.94 s. The explicit input/decoder scaling and float64 evaluation then produced the paired recovery above. Local spectral stability alone was not accepted as native recovery.

## Runtime contract and limits

`IBMCorticalStanceController` in `ihm/native/cortical_stance_controller.py` implements strict 10 ms current-interval actuation, exact patient/model/registration binding, and an eight-tensor transactional cortical checkpoint. Any sensory block flushes in-flight cortical correction to the tonic baseline; motor blocks apply to the corresponding output channels. The spinal arcs are not motor owners in this stance policy and are reported unavailable. Separate raw IBM modes retain the measured reflex integration.

The adapter applies explicitly engineered MAP/O2 availability to neural drive and temperature Q10 to E/I derivatives. The ten-second receipts above do not validate those altered physiological conditions. Additional regional sensory rates remain visible in telemetry but do not control this narrow stance policy. There is no trained descending target or walking port.

This establishes a bounded causal cortical contribution to standing in a native mechanical body. It does not establish a general human brain, biological correctness of the privileged ports/reference circuit, or walking competence.

## What the severed control does not establish

A later matched experiment rebuilt this entire pipeline on two other kernels —
IBM-1's 16-objective consolidated kernel, and that same kernel with its site rows
permuted so its learned association is destroyed and every marginal statistic
preserved. Both produce a policy that recovers the same push, at peak COM 1.7469
mm and 1.7452 mm against this bundle's 1.7844 mm, and both severed arms fall on
the same trajectory at the same 2.90 s.

That does not withdraw anything above. The persistent cortical dynamics do own
all 98 corrective commands here, and severing them does cause a fall. What it
fixes is a reading the severed arm cannot support: **severing removes the
dynamics, not the learning**, so it shows the cortex is the motor owner and says
nothing about whether the IBM kernel's trained content contributed. On this task
it demonstrably does not — the analytic readout fit succeeds equally on a kernel
with the learning destroyed. See `IBM_CURRICULUM16_KERNEL.md` and the receipt at
`data/runtime/motor-learning/ibm-kernel-comparison-20260908.json`.
