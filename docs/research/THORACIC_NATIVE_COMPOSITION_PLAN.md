# Opt-in cervical and thoracic native composition plan

The current native torso cannot be replaced by the standalone thorax subsystem.
`data/research/thoracic_mechanism/native_composition_v1/plan.json` binds the exact
current model, cervical registration/inertia, and thoracic recipes and retains
all nine cervical frame transforms and 48 thoracic material moment rows.
The target native model SHA-256 is
`899148fee8932126fd266a64f38c95a7fbecaa6651df908c8e4e890a8197b9d6`.
No additional patient mass normalization is permitted by this plan.

| Partition | Mass (kg) |
|---|---:|
| Current native torso | 27.654676965260336 |
| Nine cervical donor bodies | 7.418919568221102 |
| Thoracic subsystem, including core | 20.235757397039237 |
| Thoracic residual core | 17.93770438819818 |
| Moving thoracic material shares | 2.298053008841056 |

For each row the retained ledger gives mass, `H=integral x dm`, and
`Q=integral x x^T dm` about the same original torso origin. This Q denotes a
second moment, not generalized force. Cervical body-frame COMs and tensors are
rigidly transformed to that frame before summation. At the declared zero neck
and thorax reference, cervical plus thoracic sums reproduce original torso
mass, COM `[-0.03,0.32,0] m`, and diagonal COM inertia
`[1.520014507406313,0.7788205903211117,1.4755841071015214] kg m²`.
The independent ledger test reconstructs those quantities from H and Q.

## Assembly sequence and required mechanical interface

1. Verify the frozen target identity and all source receipts. Keep the original
   upstream torso joint, one torso pose/twist, one canonical registration and
   the original other body mass records. Create the nine cervical frames with
   the retained transforms and corrected physical mass priors. Debit their
   complete reference moments once; never add a donor whole torso mass.
2. Replace only the cervical-reduced torso ledger by its residual thoracic core
   plus the existing moving-material map. The 48 material partitions are not
   48 independently moving rigid bodies: intercostal and diaphragm fields blend
   moving anchors. Installing their masses as extra rigid bodies or a custom
   force on an unchanged torso would not implement the derived kinetic energy.
3. Assemble one coupled kinetic system. For a cervical body with COM map
   `x_i(q_neck)`, use its translational Jacobian
   `[I,-[x_i]cross,J_neck_i]` and angular Jacobian in the same parent body frame.
   Add `m_i A_i^T A_i + B_i^T I_i B_i` to the shared parent/neck matrix;
   embed the thoracic 32-coordinate metric in the same parent block. Add both
   subsystems' convective/gyroscopic terms and upstream parent-joint velocity
   maps. Apply the cervical coupler N-map with its derivative/convective terms.
   Solving the standalone thorax RHS beside an independently integrated neck
   omits these parent cross-inertias and is not an acceptable composition.
4. A native implementation must support this full configuration-dependent
   kinetic metric, or an equivalent exactly verified constrained material
   representation. Ordinary additive force components cannot substitute for
   the missing inertial operator. Backend implementation choice remains open;
   the current source-only RHS and reference ledger provide validation targets.
   Verify full-rank independent coordinates after couplers and the two rib
   locks, reference mass/H/Q, and perturbed-state material kinetic energy.
5. Preserve the existing 92 muscle paths while auditing the 78 donor paths,
   24 neck coordinates and 18 couplers. Convert old XML with unique component
   names and resolved frame sockets. Verify donor rotation convention and
   muscle moment arms by native perturbations; do not infer calibration from
   the retained centroid fit. Resolve missing C7 display geometry and uncertain
   fixed shoulder anchors without changing the physical mass ledger.
6. Transfer head/neck and thoracic surface/contact ownership exactly once, with
   source-node mappings and common world force transformations. Physical
   locked-rib joint reactions remain unavailable; the free dynamics API's
   formal eliminated-slot balancing loads are not native joint wrenches.
7. Keep native lung/airway/gas ownership explicit. Before physical recoil or
   drive is enabled, transfer overlapping native chest compliance and ideal
   pressure-source ownership and verify signed geometric displacement/work.
   Then run free-trajectory refinement, contact/supine acceptance and synchronized
   mechanical/native work checks in the opt-in model under the native queue.

Reference moment conservation does not require the articulated model's inertia
to stay constant as the neck and thorax move. Its current mass distribution
must follow both material maps. No backend conversion, default model mutation,
recoil law, muscle activation, native job or coupled acceptance is performed by
this plan. Unresolved native requirements remain enumerated in the JSON.
