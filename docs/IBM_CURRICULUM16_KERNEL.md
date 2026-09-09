# The 16-objective IBM kernel in the embodied runtime

IBM-1 at commit `267ad26` publishes a consolidated association kernel trained
against sixteen objectives at once. This document records what that kernel is,
what changed when it was wired into the IHM-1 embodied runtime, and — the part
that matters — which of the resulting behaviours are caused by it and which are
not.

This is an engineering sensorimotor simulation. Nothing below is a claim about a
validated human brain.

## What the latest kernel actually is

`IBM-1 ckpt/ibm1_curriculum16.pt`, schema `ibm1/implicit-v1`, a `(30000, 128)`
`dyn.embed` association kernel plus sixteen task heads. Retained here as
`data/models/ibm_curriculum16_kernel_v1/`.

Two facts decide how it can be used.

**It is the same interface, not a new one.** Schema string, tensor name, shape
and dtype are identical to the fused `ckpt/ibm1_implicit.pt` the existing IHM
bundles were built from. Its curriculum config records `init:
ckpt/ibm1_implicit.pt`, so it is that kernel trained further, not a different
object. Measured against it: mean per-site cosine 0.9701 (median 0.9721, min
0.8825), relative Frobenius difference 0.2250, mean per-site norm 0.2284 against
0.2243, effective rank 112.5 against 117.25. So: a real, moderate move in the
same space, slightly more concentrated.

**None of its sixteen objectives is a motor objective.** They are
`visual_eeg`, `audio_meg`, `video`, `audio_visual`, `optic_nerve`,
`cochlear_nerve`, and ten per-subject `visual_eeg_s00…s09` terms — all sensory
retrieval. IBM-1's own handoff says there is no motor corpus in that programme,
and its `embody.py` three-arm run found severing the cortical kernel changed the
command variation not at all (sd 0.02648 full against 0.02569 severed, against
0.00003 with the cord removed). Only the sixteen task **heads** carry the
sensory competence, and a head is not transferable to a motor task by
construction. So the kernel arrives here as **an initialisation, not a policy**.

The heads are therefore not retained: the bundle copies only `dyn.embed`
(15,362,274 bytes out of the donor's 1,062,702,035) and records the donor's full
digest so the omission is auditable rather than silent.

### Provenance

| field | value |
|---|---|
| donor repository / commit | IBM-1 `267ad2630e61d833e93ce553848e46b0214e90f1` |
| donor path | `ckpt/ibm1_curriculum16.pt` |
| donor sha256 | `0604a5f423db807795de2e4fc4ce1f7cf7945c6a9268165941de2cfb40ef2a4d` |
| retained kernel sha256 | `d21289cdbd76b7836ff1439a1a6af026d5a11b33ad5f4ef9aa0a944427d81acb` |
| consolidation step | 11,500 |
| IBM equation source sha256 | `9329ce624bd0ec6676fc76dc307210e8a0724ef3f79dece01c3668d6e4beb9e4` |

The donor checkpoint was **not modified**; it was read and copied. It was also
still being written by a live 60,000-step curriculum run at the time, so
`scripts/build_ibm_kernel_bundle.py` hashes and deserialises one in-memory copy
and re-hashes the file afterwards, refusing to record a digest that does not
describe the bytes it took. The manifest says the donor run was in progress:
this pins a snapshot, not a finished run.

## How it is wired

Two controller kinds were added alongside the existing ones. Nothing was
replaced; `implicit_cortical_stance` and its receipts are untouched.

| kind | motor owner | kernel |
|---|---|---|
| `implicit_curriculum16` | raw kernel then segmental cord, untrained motor heads | retained 16-objective kernel |
| `implicit_curriculum16_stance` | persistent 1024-site trained IBM E/I stance policy at 10 ms | retained 16-objective kernel |

- `ihm/assembly/controller_selection.py` now names the kind set once
  (`KINDS`, `STANCE_KINDS`, `CORTICAL_STANCE_KINDS`, `RAW_IMPLICIT_KERNELS`) and
  the other call sites read those instead of repeating literals.
- `ihm/native/cortical_stance_controller.py` gained a `BUNDLES` table. One
  adapter now serves both stance bundles; each entry names its artifact and the
  kernel it descends from, and the artifact's own `provenance` still enforces the
  pairing — `IBMCorticalStanceController` refuses to start if the kernel the
  inner implicit controller loaded does not hash to the artifact's
  `brain_checkpoint_sha256`.
- `IBMImplicitController` takes the selected kind and reports `kernel_name`
  alongside the existing `checkpoint_sha256`, so a receipt says which kernel ran.
- The workbench selects both kinds and shows an **Association kernel** row.

`scripts/train_cortical_stance_sequences.py` gained `--checkpoint`, which is the
only training change: the stance pipeline is otherwise byte-for-byte the one that
produced `ibm_cortical_stance_v1` (same teacher `f8740ca0…`, same four axis
trajectories plus the same declared extra, same 1024 sites, signed-identity
encoder, 100 epochs, lr 0.001, reserve weighting, then the same local readout fit
and the same ×100 small-signal scaling). Only the kernel differs.

## The ablations

Three kernels, one pipeline. Everything downstream of `dyn.embed` is held fixed,
so any difference is attributable to the kernel and nothing else.

| kernel | what it is |
|---|---|
| fused | `ibm1_implicit.pt`, the 34-source kernel the existing bundle was built from |
| curriculum16 | the new 16-objective consolidated kernel |
| site-permuted | the 16-objective kernel with its site rows permuted under seed 20260908 |

The permuted kernel is the control that matters. It has the same per-site
vectors, the same norms and the same marginal distribution as the real kernel —
after the adapter's 30,000 → 1,024 resampling its mean per-site norm is 0.1966
against 0.2083 — but the learned correspondence between an embedding and the
cortical site it occupies is gone. If a policy built on it works, the kernel's
learned content was not what made the policy work.

### Stage 1 — temporal imitation separates nothing

`temporal_training_axis_mse`, full against its own severed arm:

| kernel | full | severed |
|---|---:|---:|
| fused | 2.618e-07 | 2.749e-07 |
| curriculum16 | 2.967e-07 | 2.749e-07 |
| site-permuted | 2.888e-07 | 2.749e-07 |

The severed number is identical across all three because a severed policy emits
exactly the fixed tonic baseline `u0`, which no kernel touches. Two of the three
full arms are **worse than their own severed control**. Nothing here is evidence
of a learned motor contribution from any kernel, including the new one.

### Stage 2 — the analytic readout is indifferent to the kernel

The deployed policy's decoder is solved, not trained: `(I − J)⁻¹G` on the frozen
1 ms E/I Jacobian, fitted so the closed-loop DC gain matches the LQR gain.

| kernel | motor rank | condition number | gain relative error | decoder norm |
|---|---:|---:|---:|---:|
| fused | 258 | 18.857 | 4.516e-08 | 3.603 |
| curriculum16 | 258 | 18.418 | 3.898e-08 | 3.002 |
| site-permuted | 258 | **17.319** | 4.007e-08 | 2.721 |

Full rank, near-identical conditioning, gain error at 1e-8 in every case — and
the *permuted* kernel is the best-conditioned of the three. Whatever this stage
needs from a kernel, a kernel with its learned structure destroyed supplies it.

### Stage 3 — matched ten-second native pushes

100 ms pelvis-COM push at t = 1 s on the 77.6122029 kg model, full against a
severed arm that retains the same tonic baseline.

| kernel | push (N) | arm | outcome | peak COM | final COM | final speed |
|---|---|---|---|---:|---:|---:|
| fused | (+5, 0, +5) | full | recovers 10 s | 1.7844 mm | 0.1956 mm | 3.608 µm/s |
| fused | (+5, 0, +5) | sever | falls 2.90 s | 769.71 mm | 769.71 mm | 2.020 m/s |
| curriculum16 | (+5, 0, +5) | full | recovers 10 s | 1.7469 mm | 0.1923 mm | 3.498 µm/s |
| curriculum16 | (+5, 0, +5) | sever | falls 2.90 s | 769.71 mm | 769.71 mm | 2.020 m/s |
| fused | (−5, 0, −5) | full | recovers 10 s | 2.2640 mm | 0.1717 mm | 5.201 µm/s |
| fused | (−5, 0, −5) | sever | falls 1.94 s | 376.68 mm | 376.68 mm | 2.190 m/s |
| curriculum16 | (−5, 0, −5) | full | recovers 10 s | 2.1287 mm | 0.1827 mm | 6.321 µm/s |
| curriculum16 | (−5, 0, −5) | sever | falls 1.94 s | 376.68 mm | 376.68 mm | 2.190 m/s |

Two things to read carefully.

**The severed rows are numerically identical across kernels**, to every printed
digit and to the same fall time. That is correct rather than suspicious: severing
sets the association weights to zero, the mean-free readout then returns exactly
`u0`, and `u0` is the same fixed native equilibrium in both artifacts, so both
severed arms are the same open-loop trajectory. It also says exactly what the
severed control can and cannot establish. It establishes that the cortical
dynamics own the corrective motor output. **It does not establish anything about
the kernel**, because it removes the dynamics, not the learning.

**The full rows are the same to within a few percent.** 1.75 against 1.78 mm,
2.13 against 2.26 mm, on one trial per direction. The 16-objective kernel is
neither better nor worse than the fused one at this task by any margin these runs
can resolve.

### The whole-runtime arms, now measured

Every number above comes from `evaluate_native_cortical_stance.py`, which drives
native mechanics directly. The equivalent matched arms through the full runtime —
cord, BioGears physiology and the world in the loop — were originally lost four
times to `EmbodiedRuntime`'s source-freeze guard, which correctly refuses to start
when the checkout changes during initialisation; the checkout was being edited
concurrently. They have since been run on a quiescent checkout and are retained at
`data/derived/unified-world-ofb2dp7z/`:

```bash
PYTHONPATH=. .venv/bin/python scripts/verify_unified_world.py --seconds 5 \
  --controller implicit_curriculum16_stance --environment upright --scene '' \
  --arms full sever --pelvis-push 5 5
```

| arm | duration | exchanges | fallen | peak COM | final COM speed |
|---|---|---|---|---|---|
| full | 5.000 s | 250 | no | 1.7344 mm | 1.745e-05 m/s |
| sever | 2.900 s | 145 | **yes** | 499.12 mm | 2.020 m/s |

The full arm's 1.7344 mm agrees with the 1.7469 mm the isolated harness measured
for the same kernel and push, and the severed arm falls at 2.90 s, the same time
the isolated harness reports. So the kernel carries its stance behaviour intact
through the full runtime; the extra owners neither rescue the severed arm nor
degrade the full one.

This closes an *integration* gap only. It is not additional evidence about kernel
content: the site-permuted control in Stage 4 recovers the same push just as well,
and nothing in a whole-runtime arm distinguishes a learned kernel from a permuted
one.

### Stage 4 — the control that settles it

The same pipeline, same push, on the kernel whose learned site correspondence was
destroyed:

| kernel | arm | outcome | peak COM | final COM | final speed |
|---|---|---|---:|---:|---:|
| curriculum16 | full | recovers 10 s | 1.7469 mm | 0.1923 mm | 3.498 µm/s |
| **site-permuted** | full | **recovers 10 s** | **1.7452 mm** | **0.2025 mm** | **3.901 µm/s** |
| site-permuted | sever | falls 2.90 s | 769.71 mm | 769.71 mm | 2.020 m/s |

The control recovers. Peak COM excursion differs from the real kernel's by 1.7 µm
— 0.1% — and its severed arm falls at the same 2.90 s with the same trajectory,
because that arm never depended on the kernel in the first place.

Every row above is retained, with the digest of the report it came from, in
`data/runtime/motor-learning/ibm-kernel-comparison-20260908.json`, regenerated by
`scripts/collect_ibm_kernel_comparison.py`.

**So the honest reading is that the IBM kernel's learned content is not what makes
this stance policy work.** What makes it work is that the E/I network is a
well-conditioned linear filter of privileged full native state, plus a decoder
solved offline to invert that filter onto an engineered LQR gain. Any kernel that
leaves the filter full-rank and well-conditioned supports the same policy, and a
permuted kernel does. Nothing in the sixteen sensory objectives is reaching the
motor path, which is exactly what IBM-1 predicted when it wrote that there is no
motor corpus in that programme.

This does not withdraw the existing `implicit_cortical_stance` result. That
result's claim — the persistent cortical dynamics own all 98 corrective commands,
and severing them causes a fall — still holds, and holds identically for the new
kernel. What Stage 4 adds is the boundary of that claim: **the severed control
establishes that the dynamics are load-bearing; it never established that the
learning was, and it now demonstrably does not.**

## What the new kernel does and does not add

**Adds.**

- A pinned, self-contained copy of IBM-1's current shared kernel inside IHM-1,
  with donor commit and digest, so the embodied runtime no longer depends on the
  state of a live training run in a sibling repository to reproduce a bundle.
- A second stance bundle that is a genuine matched pair with the first: same
  everything, one variable. That pair is what made Stage 4 possible.
- A named ablation control kernel and the pipeline plumbing (`--checkpoint`) to
  build a policy from any kernel, which is reusable for the next kernel IBM-1
  produces.
- Slightly better numbers at every stage of the readout fit (condition number
  18.42 against 18.86, gain error 3.90e-8 against 4.52e-8, decoder norm 3.00
  against 3.60) and at three of four native push measurements. All of these are
  within what one trial per condition can resolve, and none is claimed as an
  improvement.

**Does not add.**

- Any measurable improvement in native stance recovery. 1.75 against 1.78 mm and
  2.13 against 2.26 mm, n = 1 per direction.
- Any motor competence. Its sixteen objectives are sensory; its temporal
  imitation arm is *worse* than its own severed baseline; and its policy is
  matched by a kernel with the learning destroyed.
- Any new capability in the runtime. It does not walk, does not generalise beyond
  this fixed body, mass, model and push magnitude, and is not evidence about a
  biological brain.
- Anything to the raw untrained motor path. That path is what IBM-1's own
  three-arm measurement covers, and its conclusion — the cord produces the
  dynamics, the kernel produces none — is unchanged by a newer kernel.

## What would actually test this kernel

The kernel carries sensory association, so the experiment that could show it
contributing is a sensory one: route real optic or cochlear afference through the
declared nerve path into the materialised sites and ask whether a downstream
discrimination is better with this kernel than with the permuted control. The
sixteen heads that already do that in IBM-1 are the thing to port, and they were
deliberately not copied here because nothing in the embodied runtime consumes
them yet. Until that path exists, wiring this kernel into the motor path is
plumbing, correctly done, that changes no measured behaviour — and this document
exists so that is on the record rather than inferred from a body that balances.

## Reproduce

```bash
# retain the kernel and its matched control
PYTHONPATH=. .venv/bin/python scripts/build_ibm_kernel_bundle.py
PYTHONPATH=. .venv/bin/python scripts/build_ibm_kernel_bundle.py \
  --shuffle-seed 20260908 --output data/models/ibm_curriculum16_shuffled_kernel_v1

# stance pipeline, once per kernel
PYTHONPATH=. .venv/bin/python scripts/train_cortical_stance_sequences.py \
  --sites 1024 --encoder-kind signed_identity --epochs 100 --lr 0.001 --reserve-weight \
  --teacher data/research/locomotion_control/delay_design_e5s_lepu/candidate_state1.0_R1000.0.npz \
  --extra-trajectory data/research/locomotion_control/patient_slow_lqr_delay4_com_push/trajectory.json \
  --checkpoint data/models/ibm_curriculum16_kernel_v1/kernel.pt --out <stage1>
PYTHONPATH=. .venv/bin/python scripts/fit_cortical_stance_local_readout.py \
  --artifact <stage1>/cortical_stance.pt --output <stage2>
PYTHONPATH=. .venv/bin/python scripts/scale_cortical_stance_regime.py \
  --artifact <stage2>/cortical_stance.pt --output data/models/ibm_curriculum16_stance_v1 --scale 100

# matched native push, per direction
PYTHONPATH=. .venv/bin/python scripts/evaluate_native_cortical_stance.py \
  --artifact data/models/ibm_curriculum16_stance_v1/cortical_stance.pt \
  --registration data/models/engineering_stance_v1/registration.json \
  --seconds 10 --perturb-x 5 --perturb-z 5 --output <receipt>

# whole runtime, matched arms
PYTHONPATH=. .venv/bin/python scripts/verify_unified_world.py \
  --controller implicit_curriculum16_stance --environment upright \
  --arms full sever --pelvis-push 5 5 --seconds 3
```

Native runs are far slower than real time and the `verify_unified_world.py` path
rejects a source tree that changes mid-initialisation, so it must be run on a
quiescent checkout.
