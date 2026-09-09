# Native mechanical stretch-reflex latency

The actual native articulated muscle model now has an in-situ timing measurement. This uses `NativeMechanicalStream` and IBM's `SegmentalCord` at a 1 ms exchange, with a fixed zero descending request to isolate the spinal mechanism. It does not claim a 1 ms full physiology runtime.

A 300 N upward force acts on a forefoot station of `calcn_r` from 5 to 15 ms. Native mechanics computes the resulting motion and soleus fiber length; no deformation gradient or muscle length is fabricated. The evoked spindle channel is `max(0, 200 * (native_length / initial_native_length - 1))` Hz. This explicitly references the initial native length to isolate perturbation-evoked afference, rather than interpreting tonic fiber/optimal-length strain as a new stimulus.

Three arms restore the same native checkpoint: pulse, pulse with stretch afference blocked, and no-pulse sham. Pulse-minus-sham comparisons separate the force response from passive motion of the free body.

| Event | Native clock |
|---|---:|
| First pulse-specific soleus spindle increment | 12 ms |
| First named stretch contribution | 42 ms |
| First alpha excitation difference | 42 ms |
| First mechanics interval receiving that excitation | 43 ms |
| First differing native activation endpoint | 51 ms |

The measured Ia-to-named-stretch delay is **30 ms**. Alpha changes at the same neural endpoint. The explicit next-interval exchange adds 1 ms before application. Native activation changes later: the reported activation remains at 0.01 initially despite small alpha differences. The activation-endpoint timestamp is measured separately rather than inferred from the neural command.

Peak evoked Ia is **8.3708 Hz** and peak named stretch contribution is **0.02511**. The stretch-blocked arm has zero named stretch contribution. Renshaw is retained and reported separately. By 80 ms, native activation peaks at **0.01779** in the pulse/closed-reflex arm versus **0.01** in the blocked arm, and actual fiber-length trajectories differ.

This demonstrates force → native muscle stretch → spindle increment → delayed spinal alpha → applied excitation → native muscle activation. It does not establish anatomical calibration of the chosen spindle gain or arbitrary whole-body reflex circuits. Only the tested soleus afferent channel is active, and the body is free rather than standing.

Reproduce:

```sh
.venv/bin/python scripts/verify_native_reflex_latency.py \
  --output data/runtime/native-reflex-latency-new
```

Full per-millisecond traces, source hashes, and timings are retained at `data/runtime/native-reflex-latency-20260908/receipt.json`.
