# Independent cortical stance adapter audit

The promoted `IBMCorticalStanceController` passed eight independent actual-weight
adapter checks plus five scheduler checks on 2026-09-08. These use the trained
1024-site float64 policy and explicit native-port/physiology fixtures. They do
not establish native standing, recovery, walking or biological fidelity; root's
whole-native acceptance runs and the retained mechanical recovery evaluations
are separate evidence.

Audited artifact: `data/models/ibm_cortical_stance_v1/cortical_stance.pt`, SHA-256
`0747ac07a4baa02d30964933c1f1fd1bf1285d9a30a9c35e9eac7de2e2c98042`.

Verified behavior:

- Two current-boundary 10 ms controller evaluations consume fresh accepted
  mechanical observations inside one 20 ms physiology exchange. The attached
  world advances four 5 ms steps, each holding its current controller output.
- All 98 muscle IDs reach each mechanical command. The energy fixture sums all
  98 channels and checks that their full interval energy reaches the coupling
  ledger and one physiology exchange. This is a routing test, not a native
  metabolic calibration.
- Runtime motor commands pass through the persistent trained E/I policy.
  Poisoning the raw implicit controller, spinal cord, and LQR teacher command
  methods does not affect successful execution. The declared fixed equilibrium
  tonic excitation remains part of the motor output, including the severed arm.
- Severed association edges produce exactly that tonic baseline from a fresh
  state. A sensory block flushes existing correction to the paired neural
  reference; it does not substitute an engineering balance controller.
- All eight neural tensors replay exactly. Injected late cortical failure and
  failure on the second controller evaluation restore neural tensors, accepted
  mechanical fixture state, world state, clocks and commands before physiology
  is touched. Nonfinite values in each of the eight checkpoint tensors are
  rejected without mutation.
- Nonbaseline oxygen and temperature inputs change cortical outputs. Unsupported
  descending targets and malformed sensory, motor and physiology inputs reject.
- Model hash, registration bytes, 77.6122029 kg patient mass and 10 ms sample
  bindings reject mismatches. Adversarial NaN/boolean mass and interval values,
  malformed reference declarations, and unsupported four-state runtime artifacts
  reject. The audit found missing finite checks and missing eight-state admission
  checks; the adapter author corrected these before final acceptance.

The policy consumes privileged full native state through engineering sensory
ports. Additional regional sensory inputs are validated but explicitly do not
control this stance decoder. Training used an offline LQR teacher; the audit
does not mistake absence of a runtime teacher for absence of an engineering
prior or for training without supervision.

Reproduce without native jobs:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=.:scripts \
  .venv/bin/python -m unittest scripts.verify_cortical_stance_controller_audit \
  scripts.verify_control_exchange
```

Result: 13 tests passed. No adapter, factory or policy source was changed by
the independent auditor; fixes were coordinated with their owner.
