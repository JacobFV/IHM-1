# IBM body compatibility: executable bounded adapter

The IHM adapter runs actual IBM package bytes, not extracted registry metadata or a
reimplementation of donor transfer functions. It does **not** establish that the
whole IBM materialized graph runs, or that whole-body neural physiology is complete.

## Source identity

`scripts/vendor_ibm_backend.py` copies the actual complete `ibm` source package to
`data/derived/canonical/ibm-backend/source/ibm`, preserving the donor checkout. The
artifact manifest lists SHA-256 for every source file and the donor tracked dirty
status. The recorded snapshot contains 94 files, 1,939,813 bytes, from IBM-1 commit
`63a5668f6e456b8ca5df5848e5d839cf77357bc7`. This commit alone does not identify the
runtime: the checkout contains modified files, including `runtime/step.py` and
`runtime/ensemble.py`, whose actual bytes are included. The package identity is
`0bcb65606a350d586c1c5fc1634913aee7f90c7de4201eef57770c6af715346f`, the SHA-256 of the
sorted JSON per-file digest map. IHM pins that digest in executable code and
verifies every file before importing. No donor dataset, learned weight, or complete
brain artifact is implicitly included. Another IBM version already imported in
the process is rejected. NumPy and SciPy are needed; this route does not need Torch.

Regeneration reads the current donor bytes; if they change, the new artifact is
rejected until a deliberate source pin update is reviewed. The generated artifact
is local derived data; copying this repository without its artifact requires
regenerating exactly the pinned snapshot, including its recorded dirty contents.

## Public execution surface

```python
from ihm.brain.ibm_backend import IBMBackend
backend = IBMBackend()                         # verifies the source artifact
model = backend.materialize_tactile(sites_m)  # actual MaterializationRequest + build
result = backend.run_window(indentation_um, sites_m=sites_m, kind="rapid",
                            dt_s=.002, delay_s=.01)
values, residual, audit = result.values, result.residual, result.audit
```

`sites_m` is a finite nonempty `(sites,3)` array in the canonical body frame.
The geometry adapter converts metres to millimetres and declares an explicit
discrete `body_surface` sampling. It does not relabel body coordinates as
`subject_t1`, use an invented head volume, or claim a receptor census. The actual
request targets `transduction.mechanoreceptor`, selects the body-surface region,
and overrides transduction with `mechanoreceptor_rapid` or `mechanoreceptor_slow`.
Both build paths must allocate receptor sites or raise. Dependency tracing still
reaches many unsupported processes. `strict=False` is used to inspect that bounded
view; `materializer_missing`, `materializer_notes`, and
`full_materialized_graph_executable=False` remain visible in the audit.

The runtime adapter explicitly selects the materialized receptor block and
constructs runtime `Coupling` objects. It cannot directly call donor
`couplings_of(model)`: `MaterializedModel.edges` maps topology names to `EdgeSet`,
whereas that runtime convenience function treats iteration elements as coupling
objects. No generic graph conversion is claimed. The selected receptor transfer
comes from the registered donor implementation; its scalar parameters use donor
prior medians, filtered to the callable signature. For example rapid adaptation
uses the declared 20 ms prior, not the transfer function's 30 ms default. Unknown
caller parameters are rejected; unused priors such as `threshold_um` do not
silently become dynamics. Uncertainty in those priors is not propagated.

`run_window` also supports `twitch` (the actual donor primitive) and `activation`
(the registered activation-plus-twitch transfer). Those two are explicit effector
probes, **not materialized whole-body effector graphs**. Omitted site coordinates
use a declared synthetic test support; application integration should pass actual
canonical coordinates. Inputs accept `(samples,)` or `(sites,samples)` and outputs
always have shape `(sites,samples)`.

## Units and temporal semantics

The receptor transfer gain is specified per micrometre indentation. Body contact
pressure in Pa cannot be substituted as its input. `pressure_to_indentation_um`
requires explicit positive `stiffness_pa_per_m` and computes `d_um = p/k * 1e6`.
This is a caller-chosen local linear foundation assumption, not measured skin
mechanics or a donor calibration. Donor mechanical displacement is mm, requiring
another explicit factor of 1000 to reach micrometres. Donor head-volume mechanical
and thermal carriers have no validated whole-body support mapping here.

The actual runtime solves drift equations. With the selected transfer and an
explicit IHM unit-rate self leak, the equation is

`dy/dt = H*x - y`, hence `y_hat = H*x_hat/(i*omega + 1)`.

The additional one-second leak is an adapter convention, not a physiological fit.
Outputs are labeled donor drift response, not calibrated mV or muscle N. A pure
transfer plot and this runtime output are measurably different; the parity verifier
asserts that distinction. The hard-clamped input is applied through donor `Clamp`,
and output is generated by donor `solve_window`, with convergence checked.

Windows are Fourier-periodic boundary-value solves. Wraparound and pre-stimulus
response are possible. There is no claim of causal streaming, replay equivalence,
or continuity across calls. `run_window` does not persist `Carry`; wrapping repeated
calls with Python history alone does not change these semantics. Donor `advance`
starts a new carry within each invocation and is not used to imply continuity.
Delays supplied to the coupling use the actual runtime spectral phase shift, but
are caller values, not reconstructed tract lengths or validated afferent anatomy.

## Verification and remaining gaps

Run `.venv/bin/python scripts/verify_ibm_body_parity.py`. It checks real nonempty
materialization for both receptor choices, analytical parity with donor source
transfers, actual runtime parity for rapid/slow/twitch/activation with a 10 ms
delay, rejection of empty geometry/invalid stiffness/unused parameters, and source
corruption detection. A two-site, 512-sample pulse at 2 ms yields runtime solution
errors below `3e-17` against an independently evaluated frequency-domain solution.
That is numerical parity, not biological validation. No reconstructed topology
delay, body-to-brain closed loop, autonomic circuit, thermal/nociceptive/chemical
transduction, learned controller, motor-unit recruitment, force calibration,
long-window causal behavior, or uncertainty propagation is validated by this test.
