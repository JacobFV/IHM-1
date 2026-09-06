# Bounded supine support acceptance

The existing `data/derived/articulated-acceptance-noo34jds` receipt establishes short native branching, checkpoint, metabolism, and force behavior for the 92-muscle v2 model. It explicitly does not establish settled support. At 2 ms, its contact resultant remains approximately 97.3% of body weight away from static force balance, despite a near-zero **dynamic** momentum audit. The latter checks Newton's law including acceleration; it does not establish equilibrium.

## Diagnostic protocol

`scripts/verify_supine_support.py` defaults to native-free classifier fixtures and a retained-transient check. `--run-native` is reserved for a coordinated heavy-resource slot. The model, transferred contact coefficients, default excitation, and solver parameters are unchanged. Garments and external forces are disabled. Canonical attachments add no inertia.

A fresh articulated plant uses `whole_body_arm26_v2/registration.json`, 92 muscles, and 77.6122029 kg. Native observations are sampled every 5 ms. Stages end at 0.02, 0.1, 0.3, 0.6, and 1.0 simulated seconds. A 60-second wall deadline terminates the native process on timeout. Divergence or convergence ends the sequence early; an unsettled horizon is a failed equilibrium acceptance, not a request to automatically run longer.

Each stage acquires and restores an in-process checkpoint before advancing. The receipt verifies equal native state observations excluding response kind, records snapshot hashes, stage outcome and wall time, then releases the checkpoint. Full native initial/stage snapshots and every diagnostic sample are retained. These JSON observations cannot restart the process; the complete checkpoint exists only while its native session remains alive. Existing native execution manifests retain model, source, build and augmentation identities. A final receipt binds the harness SHA-256. If a deadline interrupts a stage, its last receipt remains incomplete and the final report fails closed.

## Measured quantities

- Static imbalance is `norm(contact + mass * gravity) / (mass * norm(gravity))`, in native source coordinates. Normal support is projected opposite gravity and divided by weight.
- COM acceleration is reconstructed from the native identity `M(a-g)-contact = residual` with no applied load. Independent interval acceleration uses the difference of mass-weighted COM velocities, including each body's angular contribution to COM velocity. Neither a small audit residual nor a small initial velocity alone is equilibrium evidence.
- Native kinetic energy is normalized by total native body mass. The trailing-window linear slope tracks continued accumulation or decay. Maximum segment angular speed and total COM speed prevent a motionless-COM interpretation of ongoing segment motion.
- Proxy penetration is `max(0, plane_x + sphere_radius - sphere_center_x)` for each native posterior sphere. This measures engineering proxy overlap, not tissue deformation or mattress compression.
- Constraint error is the larger native QErr/UErr norm. These aggregate mixed-coordinate norms are numerical checks, not physical lengths. Momentum residual is normalized by weight.

## Predeclared engineering criteria

Convergence requires at least 0.4 simulated seconds and **every sample** in the final 0.2-second window to satisfy static imbalance ≤2% weight, COM acceleration ≤0.02 g, interval COM acceleration ≤0.02 g, kinetic energy ≤0.00005 J/kg, COM speed ≤0.01 m/s, maximum angular speed ≤0.05 rad/s, penetration ≤0.01 m, aggregate constraint error ≤1e-5, and dynamic audit residual ≤1e-7 weight. Absolute fitted kinetic-energy slope must be ≤0.0001 J/kg/s. These conservative diagnostic tolerances are engineering choices, not validated biological acceptance bands.

A stage stops immediately on a sampled support magnitude above 10 weights, COM acceleration above 20 g, kinetic energy above 2.5 J/kg, penetration above 0.05 m, constraint error above 0.01, or audit residual above 0.01 weight. It also stops when a complete 0.2-second window shows kinetic energy above 0.1 J/kg and greater than four times its starting value (with a 0.001 J/kg denominator floor). Nonfinite values, native failure, checkpoint mismatch, or deadline expiration fail closed. Sampling cannot rule out unobserved peaks between samples.

## Cost and source-only verification

The retained short acceptance used 1.32 wall seconds and approximately 99 MiB peak child RSS for initialization plus four 2 ms advances. This diagnostic makes at most 200 advances and five checkpoint/restore checks. Provisional full-run budget is 15–60 wall seconds; the hard deadline is 60 seconds, and the existing native adapter enforces a 4 GiB address-space ceiling. This is an estimate, not an extrapolated timing guarantee. One authorized native run should provide the next decision; repeated long runs or coefficient sweeps are outside this protocol.

Source-only checks cover exact static balance, free fall with a perfect dynamic audit, insufficient initial observation, sustained synthetic equilibrium, nonzero kinetic energy, penetration divergence, kinetic growth divergence, penetration sign, and nonfinite rejection. The retained 2 ms result is checked only for its failure to meet static balance. No actual support-convergence result is claimed by these fixtures.

## Authorized native result — 2026-09-05

One bounded run retained at `data/derived/supine-support-2qf4zcjh` stopped after **0.6487 wall seconds**, before completing its first 0.02 s stage. The first 5 ms advance succeeded; the attempted second 5 ms advance raised `ValueError: Out-of-domain native muscle metabolism`. The exact traceback and interrupted checkpoint receipt are retained. The stage-start restore matched its native observation. The owned process was closed and no native mechanical stream remained running afterward. This run fails equilibrium acceptance and does not justify extending the horizon.

At the last successful sample (0.005 s), normal support was **0.026236 weights**, static imbalance was **0.973781 weights**, COM acceleration was **0.973781 g**, and independently differenced COM acceleration was **0.973163 g**. Kinetic energy rose from zero to **0.029936 J/kg**, COM speed reached **0.047734 m/s**, and maximum segment angular speed was **9.725908 rad/s**. Sampled proxy penetration was zero and aggregate constraint error was 7.85e-17. The dynamic momentum audit was only 9.09e-17 weights; that excellent force-accounting residual coexists with near-free-fall acceleration and therefore cannot establish support equilibrium.

The stopped native metabolism guard is a separate blocker requiring diagnosis before further support runs. No solver/contact coefficient tuning or additional native run was performed. The original native output remains untouched. Evidence SHA-256 values:

- `report.json`: `8abaa30b71d154ba43adc7cacfdfbc415efa103ca416f8ab878a6261f5be8753`
- `failure.json`: `5ec19ff17ec5a8f37bc53d4025bdedf33ed2fe68b46ab7b051a593c066bdecab`
- `checkpoint_receipts.json`: `99ec22932b7e99b67c53904fd81ad1d545e892d039902f4f7c07135b74bf3078`
- `samples.json`: `5c81bdf52653cd06c760cdd15feb0ca9e83e6f8b984bcfa604cb321ed44f6402`
- `initial_native.json`: `6e2b9bb9b68c4a8b02a6b83a0f3e1137a7989bdb097762326e3109fefa48ae39`
- `plant/native/execution.json`: `858addb1c17250ebeabef7c939dd782930539fbb76ef51a2a7911c6a9dd95a87`
