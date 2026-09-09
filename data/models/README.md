# Local learned motor materializations

These are separate task artifacts. The original IBM fused checkpoint is unchanged.

* `ibm_curriculum16_kernel_v1/kernel.pt`: the association kernel from IBM-1
  `ckpt/ibm1_curriculum16.pt` at commit `267ad26`, consolidation step 11,500,
  copied out of the donor repository without modifying it. Only `dyn.embed` is
  retained; the donor's sixteen task heads are sensory retrieval readouts and
  none is a motor head, so carrying them would be 1.06 GB that no motor path can
  use. `manifest.json` pins the donor digest, the donor commit, the objective
  list and the fact that the donor run was still training when the snapshot was
  taken.
* `ibm_curriculum16_shuffled_kernel_v1/kernel.pt`: the matched ablation control
  for the above — the same kernel with its site rows permuted under a declared
  seed. Every per-site vector, its norm and the whole marginal distribution are
  preserved exactly; only the learned correspondence between an embedding and the
  cortical site it sits on is destroyed. It exists so that a stance policy which
  works can be asked whether it needed the kernel's content at all.
* `ibm_curriculum16_stance_v1/cortical_stance.pt`: the same 1,024-site stance
  pipeline as `ibm_cortical_stance_v1`, trained from the 16-objective kernel
  instead of the fused one. Same teacher, same trajectories, same epochs, same
  readout fit, same small-signal scaling; only the kernel differs. See
  `docs/IBM_CURRICULUM16_KERNEL.md` for its matched receipts and for what the
  ablations do and do not establish.
* `ibm_curriculum16_shuffled_stance_v1/`: the ablation control for that bundle,
  built from the permuted kernel through the identical pipeline. It is **not**
  wired to any controller kind and is not a deliverable controller. It is
  retained because it recovers the same native push as the real-kernel bundle,
  which is the evidence that this stance task does not depend on what the IBM
  kernel learned.

* `ibm_cortical_stance_v1/cortical_stance.pt`: 1,024-site persistent IBM E/I cortex
  with trained copied embeddings, a locally fitted decoder, 516 signed sensory
  ports and 508 disjoint motor ports. The 258-state / 98-muscle contract binds the
  77.6122029 kg patient model. Float64 and explicit small-signal scaling preserve
  cortical feedback accuracy. The manifest retains lineage, equations and two
  matched ten-second native COM-push receipts. Full cortex recovers both pushes;
  sever retains the same tonic baseline and falls at 1.94 and 2.90 seconds. This
  establishes bounded native cortical motor causality, not walking or a human
  brain. Whole-physiology/world validation is separate from these native receipts.
* `ibm_cortical_ankle_v1/cortical_motor.pt`: 128-site persistent IBM E/I cortex,
  trained copied embedding and decoder. Exact equation source, training summary
  and embedding-reset ablation are alongside it.
* `ibm_ankle_primitive_v1/motor_kernel.pt`: earlier eight-site association kernel
  with fixed ports, without cortical E/I integration.

The ankle artifacts imitate synthetic PD feedback using privileged joint angle
and speed. The stance artifact uses privileged full native state and an explicit
paired zero-input neural reference. Its engineering LQR teacher supplies offline
training targets only; no live teacher controls its runtime motor output. None
is a walking controller or a validated human brain.

Native results and scope are documented in `docs/IBM_CORTICAL_STANCE.md`,
`docs/IBM_CORTICAL_MOTOR.md` and `docs/IBM_NATIVE_MOTOR_LEARNING.md`. Integration
adapters default to these retained artifacts; they also require the paired IBM
checkout for its anatomical/cord model and original implicit source identity.
