# Causal selected-transfer materialization

`ihm.brain.causal.CausalIBM` realizes the selected pinned IBM donor LTI
transfer as persistent finite-dimensional state, independently at each supplied
support point. It does **not** make IBM's materialized graph or periodic
`solve_window` runtime stream. The source identity is provided by the validated
`IBMBackend` source snapshot.

```python
from ihm.brain.ibm_backend import IBMBackend
from ihm.brain.causal import CausalIBM

receptor = CausalIBM(IBMBackend(), kind='rapid', sites_m=[[0., 0., 0.]],
                     delay_s=.012, semantics='direct')
response = receptor.advance([5.], dt_s=.002)  # 5 micrometres indentation
saved = receptor.checkpoint()                # JSON-serializable
receptor.restore(saved)                      # rollback or identically configured restart
```

Input is zero before time zero and piecewise constant over each `advance`
interval. Output is the signed value at the interval's right endpoint. A
timestamped input queue implements the explicit caller delay, including delays
that are not integer multiples of the timestep. Integration splits at arrival
times; there is no response before input arrival. The delay does not infer
axonal topology or conduction velocity. Support coordinates must be finite and
nonempty; each site has independent filter state.

## Exact selected equations and provenance

The coefficients follow the actual donor functions in pinned
`ibm/processes/transduction.py` and `ibm/processes/effector.py`, using the same
registered prior medians selected by `IBMBackend`, with function defaults for
parameters absent from the registry. With Laplace variable `s`:

| Kind | Direct transfer | State order |
| --- | --- | --- |
| rapid | `g ta s / ((1+ta s)(1+tm s))` | 2 |
| slow | `g (ta s+f) / ((1+ta s)(1+tm s))` | 2 |
| twitch | `g / (1+tc s)^2` | 2 |
| activation | `1 / ((1+tr s)(1+tf s)(1+tc s)^2)` | 4 |

Rapid has zero DC response and signed negative release. Slow retains the
donor's static fraction. Activation preserves both activation poles and the
double twitch pole: it is the donor's fixed linear cascade, not a nonlinear
switch between separate rise/fall rates. No gains are fitted or inferred.

`semantics='direct'` realizes `H(s)`. `semantics='drift'` adds a pole at
minus one inverse second, realizing `H(s)/(s+1)` (one additional state), the same transfer
semantics as the selected periodic donor drift experiment but with a causal
initial value problem. A finite zero-state transient is not expected to equal
a periodic window with circular prehistory.

SciPy `tf2ss` constructs continuous state space and `cont2discrete(method='zoh')`
provides exact integration for held input, up to floating-point arithmetic.
Every construction independently calls the actual donor transfer at DC and
801 logarithmic frequencies from `1e-4` through `1e5` rad/s and compares it
against the state-space resolvent. A normalized maximum absolute response
error above `1e-9` rejects the materialization. This is coefficient
materialization, not rational fitting or broadband gain calibration. Receipt
`audit` records the donor identity, parameters, poles, order, semantics,
supports, units, and measured error.

Time constants must be in `[1e-5,1e3]` seconds with maximum ratio `1e6`;
static fraction must be in `[0,1]`, and absolute gain at most `1e6`.
Unsupported/nonfinite parameters and failed numerical parity reject execution.
These are numerical applicability bounds, not biological validation ranges.
Checkpoints bind source, kind, parameters, supports, delay, and semantics;
changing any of these requires a new state trajectory. Invalid restores do
not mutate the existing trajectory.

## Units and limits

Rapid/slow input is physical indentation in **micrometres**, not pressure in
pascals. A skin mechanics model must supply the conversion; passing contact
pressure directly changes the model's units. Twitch/activation input is unit
drive. Output remains donor response in its selected transfer convention,
not calibrated receptor millivolts, firing rate, muscle activation fraction,
or newtons. Coupling this output to a generic afferent rate is an additional,
unfitted modeling assumption. These four scalar LTI mechanisms do not establish
whole-body neural, chemical, thermal, visceral or motor circuitry.

## Verification

Run `.venv/bin/python scripts/verify_ibm_causal.py`. Its JSON output includes
all eight construction receipts. Checks cover independent actual-donor vs
discrete ZOH frequency response, analytic twitch step and rectangular pulse,
fractional conduction delay, zero pre-delay output, rapid adaptation/release,
step halving, checkpoint rollback and JSON restart, changed parameter/delay
identity rejection, empty support, and invalid inputs. The frequency comparison
allows the explicitly bounded small phase error introduced by finite input
holds; the analytic held-input tests use `1e-12` absolute tolerance.
