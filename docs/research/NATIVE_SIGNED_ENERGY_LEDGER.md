# Native signed muscle energy ledger

The mechanical stream exposes three cumulative joule fields:
`muscle_metabolic_energy_j`, `signed_active_fiber_work_j`, and
`muscle_heat_energy_j`. Metabolic energy and signed active-fiber work use the
same before/after native metabolism samples and endpoint trapezoidal weights
for each accepted command interval. Heat is their difference. This matches the
existing metabolic quadrature; it does not claim exact integration through
internal integrator steps.

The separate `positive_active_fiber_work_j` retains its existing positive-only
mechanical meaning. In the articulated projection, `positive_muscle_work_j`
is the latest exchange's positive-only work, while all three ledger fields
remain cumulative. Never substitute positive work into the signed heat ledger.

`metabolic_reference` contains immutable initialization powers `M0_w`, `W0_w`,
and `H0_w=M0_w-W0_w`. It is not automatically a validated resting baseline.
The coupling owner must bind that reference to model/state/control identity.
For two accepted snapshots separated by dt, difference each cumulative ledger
and divide by dt to obtain Mbar, Wbar, Hbar. Subtract their respective frozen
reference powers to obtain DeltaM, DeltaW, DeltaH. The identity
DeltaM=DeltaH+DeltaW holds within scale-aware floating-point roundoff. Signed
work preserves eccentric absorption. No demand increment is clamped here.

Checkpoint/restore includes both new cumulative ledgers. Native command failure
restores both alongside state, excitation, load and existing accumulators;
the reference is immutable. Tendon/passive storage and external boundary work
retain their separate mechanical ownership. This accounting is an incremental
demand boundary, not proof that native chemistry supplied all requested work
or that BioGears has complete absolute chemical-to-thermal closure.

Run `scripts/verify_native_signed_energy_ledger.py --run-native` only in the
coordinated heavy slot after rebuilding the native stream. It checks matching
endpoint quadrature, signed interval incidence, nonzero-ledger restore/replay,
and rejection after a partially parsed excitation command. The articulated
native verifier also checks the projected fields and reference replay.

## Retained acceptance, 2026-09-05

Build `data/runtime/mechanical-stream/build-pq9lgmjn` passed. Actual native
fixture `data/derived/signed-energy-ledger-21z1omhe/report.json` passed two
0.002 s intervals, replay and rejected-command rollback (70,076 KiB maximum
child RSS). The second interval had metabolic energy 0.2665360523 J, signed
work −0.0210903699 J, and heat 0.2876264222 J. Cumulative signed work was
−0.0187190077 J while positive-only work was +0.0259328507 J, demonstrating
that the fixture exercises eccentric absorption. The retained
`projection_report.json` confirms all cumulative/power/reference fields survive
the canonical projection across 2,408 entities, without another native run.
