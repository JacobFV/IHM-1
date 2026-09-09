# Mechanical metabolic reference

The embodied factory now measures its fixed metabolic reference from a separate
zero-excitation mechanical counterfactual, then restores the exact initial
mechanical checkpoint before the common physiology/controller/mechanics clock
starts. It does not change the native substrate budget or suppress rejection.

The previous reference was one instantaneous sample at the model's default
initial muscle activations. In the 77.612 kg native patient instance, it was
130.681 W. A controller sending low activation subsequently produced a sustained
negative increment against that initial condition. The unmodified signed native
port rejected a request of -10.429 W against 10.223 W available decrement near
1.2 simulated seconds, despite the two-second exchange lag.

The new reference commands every native muscle to zero, advances 0.2 s, and
integrates chemical, signed active-fiber work and heat over the next 0.6 s.
Those durations and zero-command condition are engineering choices. The body
can move under gravity during measurement; this is neither a settled standing
measurement nor identified human resting metabolism. In the measured supine
instance the resulting reference is M=72.403 W, W=-11.072 W, H=83.476 W.
Its meaning is the cost of a declared quiescent-command counterfactual. It is
frozen before task execution and does not adapt to controller output or to the
native allowance. Upright and free measurements are explicitly separate
environment-specific counterfactuals, not a claim of resting posture.

The factory retains `metabolic_reference.json`, its hash and contents in the
session manifest, and the reference in each completed frame's coupling record.
The native plant still reports its original instantaneous sample independently.
The reference receipt records the commanded muscle IDs and measured interval
energies. Calibration advances are followed by checkpoint restoration;
physiology is not advanced during calibration. Native physiology commands and
acknowledgments are retained separately in its receipts journal.

The existing exchange lag remains unchanged. Each accepted step now accumulates
the signed chemical, heat and work energy awaiting exchange. For each channel,
integrated raw incremental energy equals integrated native exchanged energy plus
reported pending energy. Pending quantities may be negative for a decrement:
they are accounting balances, not measured ATP or thermal reservoirs. Rejected
steps do not commit these balances. Chemical = heat + signed work is preserved
for the reference, raw rates, sent rates and pending balances.

Reproduce the actual native regression with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=. \
  .venv/bin/python scripts/verify_sustained_metabolic_exchange.py --seconds 5
```

It writes a bounded-duration acceptance or failure report and energy trace under
`data/derived/sustained-metabolic-*`. It checks native/mechanical clock agreement
and raw/sent/pending energy conservation. Successful completion proves only the
requested duration and conditions; the retained physiology model's long-run
glucose and acid-base limitations remain separate issues.

Fast regression:

```bash
.venv/bin/python -m unittest scripts.verify_embodied_runtime \
  scripts.verify_candidate_embodied_factory scripts.verify_regional_embodied_factory
```

Measured on 2026-09-08: the calibrated supine native run completed 10.0 s / 500
signed exchanges, with zero unmet muscle energy and exchanged chemical demand
0.618–23.918 W. The command/acknowledgment-derived receipt is retained at
`data/derived/physiology-stability-0w0st0kr/verification.json`. This run used the
baseline regional controller, not the newly added implicit-kernel adapter.
The three fast regression modules above passed 22 tests.
The latest-source retained conservation regression separately passed 3.0 s / 150
exchanges in `data/derived/sustained-metabolic-m08emxwj/verification.json`:
89.456 J raw chemical increment = 47.744 J sent + 41.711 J pending. Its heat
and signed-work channels satisfy the same per-step accounting checks.
