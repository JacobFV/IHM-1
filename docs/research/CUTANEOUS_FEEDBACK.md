# Executable cutaneous feedback

`ihm.assembly.cutaneous_feedback.CutaneousFeedback` converts explicitly registered
material-site contact forces and caller-supplied native skin temperature into
causal pinned IBM receptor responses. It returns `sensory_inputs_hz` for the
existing `BodyBrain` operator. It never advances that shared brain itself.

This is implemented neural stimulation with declared engineering conversions, not
an anatomically identified axon graph, receptor census or calibrated perception.
The tests invoke the actual preserved IBM transfer functions through `CausalIBM`
and show that their output changes the existing `BodyBrain` population activity.

## Interface and ownership

```python
feedback = CutaneousFeedback(
    root,
    sites=[{
        'id': 'caller-registered-material-site',
        'position_m': [x, y, z],
        'normal': [nx, ny, nz],  # unit outward normal
        'contact_area_m2': area,  # None means no pressure transduction
        'stiffness_pa_per_m': foundation_stiffness,
        'sensory_region': 'brain-rh-postcentral',
        'reference_temperature_C': reference_temperature,
        'support_basis': 'caller material registration and provenance',
    }],
    recruitment_hz_per_response=explicit_gain,
    delay_s=explicit_delay,
)
result = feedback.step(dt_s, {
    'time_s': feedback.time_s,
    'contacts': [{'id': site_id, 'force_n': force_vector_n}],
    'skin_temperature_C': native_skin_temperature,
}, sensory_blocks=[])
```

Every force ID must name a constructor site. If a contact also supplies `point_m`,
it must equal that registered site's position within 1 nm; silently retargeting a
material point is forbidden. Moving body support requires an explicit registration
update/materialization; this bounded module does not infer a moving surface normal.
Missing contacts mean released contact. Duplicate contacts are rejected: aggregate
physical forces explicitly before this boundary. The caller owns spatial registration,
area, the outward normal and physical provenance. Synthetic fixture supports must
remain labeled as such; a site is not registered merely because an arbitrary force
payload includes coordinates.

The report preserves physical units, sample and response times, per-site signed
responses, block status, support provenance and cortical assignments. An omitted or
`None` skin temperature disables and resets the thermal channel with an explicit
engineering dropout policy. It is never interpreted as a physical return to the
reference temperature. Scalar native skin temperature is transferred uniformly to
sites; this is not a locally resolved temperature measurement.

`SensorimotorController.step(..., additional_sensory_inputs_hz=previous_inputs)`
now validates exact canonical population IDs and finite rates in [0, 1000] Hz,
then adds these inputs to the existing regional sum before its sole `brain.step`
call. The report includes both `additional_sensory_inputs_hz` and pooled
`brain_sensory_inputs_hz` (each region capped at 1000 Hz).
Use the cutaneous endpoint from the previous exchange. Thus a physical sample held
on `[t, t+dt]` produces receptor output at `t+dt`, available to the next shared brain
interval. Do not run a second `BodyBrain` instance as part of the body session.
Root body-session wiring is a separate integration responsibility.

## Physical and neural conversions

Compression is `max(0, -force_n dot outward_normal)`. Tangential force does not
become normal pressure. Known area permits `pressure_pa = compression_n / area_m2`.
The existing `pressure_to_indentation_um` operator applies the explicit linear
foundation assumption `indentation_m = pressure_pa / stiffness_pa_per_m`, then
converts metres to micrometres. Foundation stiffness has units Pa/m and is not a
Young modulus. Unknown area yields null pressure and indentation, preserves raw
force, and supplies no mechanical receptor stimulus.

Rapid and slow receptor channels use the existing pinned `CausalIBM` exact
zero-order-hold realization and source prior medians. The thermal channel selects
the actual registry implementation `thermoreceptor_static_dynamic`, whose transfer
is the same source static/dynamic operator supported by the slow realization.
Its own source parameter medians are used; the adapter verifies frequency-domain
parity with the actual selected donor function. Thermal stimulus is degrees Celsius
relative to the explicitly supplied reference. This reference subtraction is an
engineering operating-point assumption. Source parameters for nonlinear warm/cold
bells are not used by the source LTI law and are not invented here.

Each source response remains signed, including release undershoot. The explicit
cortical recruitment prior pools the absolute response magnitudes, multiplies by
`recruitment_hz_per_response`, and caps each region at 1000 Hz. This pooling is an
inferred observation law, not identified receptor firing; it can recruit the
caller-assigned canonical population on both onset and release. Gain, spatial
assignment, area, foundation stiffness, reference temperature and delay are all
visible in materialization identity. Zero recruitment gain disables cortical input
while retaining physical and receptor observations.

A sensory block immediately clears the blocked site's queued stimuli and resets
its receptor state. Unblocking starts fresh and observes the configured delay.
This is an explicit engineering ablation, not a detailed local-anesthetic model.
No refractory period, spikes or additional adaptation are invented. Adaptation
comes solely from the selected pinned source transfers.

## Verification and remaining scope

```sh
PYTHONPATH=. .venv/bin/python scripts/verify_cutaneous_feedback.py
```

Eight bounded tests cover dimensional conversion, source parity, actual population
response, conduction delay, queued-signal block purge, release undershoot, thermal
transients, missing temperature, unknown area, normal orientation, point identity,
nonfinite input rejection, checkpoint replay and failed-restore atomicity. They do
not launch native mechanics, browser sessions or heavy jobs.

Remaining limitations include moving/deforming skin registration, calibrated
contact compliance, regional thermal fields, receptor distributions and thresholds,
nonlinear temperature tuning, pain, axon topology, anatomical somatotopy and
calibrated cortical recruitment. No IBM-1 source is edited. No force or heat is
fed back into mechanics or systemic physiology by this receptor adapter.

## Native deformation integration contract

When a native material foundation already solves contact indentation, its actual
normal deformation should be supplied directly to the receptor as
`indentation_um = native_indentation_m * 1e6`. Do not apply the linear Pa/m
foundation conversion again. Preserve native triangle force, known area, pressure,
material/triangle ID, sample time, and deformation provenance alongside it.

The direct-indentation port is implemented. Sites declare
`mechanical_input='native_indentation'`, a manifest SHA256, quadrature and triangle
identity, indentation basis and area basis. Each present contact supplies `id`,
`force_n`, `indentation_m`, matching material identity and indentation basis.
This mode rejects a second stiffness parameter. Native area supports a pressure
audit without becoming an input to another deformation model. Eleven cutaneous
fixtures cover both input modes and malformed or mismatched source identities.
Moving supports and native endpoint field binding remain integration work.

The controller integration is tested by
`PYTHONPATH=. .venv/bin/python scripts/verify_sensorimotor.py`. The matched fixture
holds all mechanical and descending inputs equal, supplies real skin-force-driven
IBM receptor output to one controller, and blocks the other receptor. After forty
10 ms exchanges the shared brain activity and actual controller motor excitation
differ, while both brain and receptor clocks equal 0.4 s. The first exchange is
identical, confirming that future receptor endpoints are not injected early.


## Shared runtime exchange

`EmbodiedRuntime(..., cutaneous=registered_receptor)` now checkpoints receptor
state alongside mechanics and neural state. It requires an explicit
`mechanical_state.cutaneous_contacts` list: absence is an error. Generic static
force sites may omit released contacts, but native material sites must each emit
one explicit observation, including zero indentation when released. Missing or
duplicate selected native rows fail both protocol and runtime coverage checks. It samples the accepted physical state at the start
of the exchange and delivers the previously accepted receptor endpoint to the
single shared brain step. New receptor endpoints become available next exchange.
`skin_sensory_blocks` removes both queued site contributions and receptor state;
pre-native failures restore receptor state with the other reversible owners.
Frames retain the receptor output and source identity. Native skin temperature
is used only when the corresponding observation exists.

The runtime fixture executes the actual pinned receptor model and verifies this
latency, block behavior, synchronized clocks and missing-contact rollback.
The workspace factory now accepts `surface_contact_manifest` and an explicit
`cutaneous_configuration` with `regions` mapping native `skin-contact-N` IDs to
canonical brain populations, `recruitment_hz_per_response`, and
`reference_temperature_C`. It requests those native quadrature sensors, binds
actual material identities and current positions/normals through the same rigid
canonical frame, and retains receptor priors in the session manifest. Direct skin
indentation drives transduction; bed indentation is separate. Selection is bounded
to 64 sensory sites while all source contact quadrature points remain in physics.
No cortical mapping or gain is silently represented as measured connectivity.
The default factory does not enable this optional setup without its explicit
surface and cortical mapping inputs. Supported resting pose and accepted coupled
whole-body contact/reflex trajectory remain outstanding.

Actual native factory acceptance is retained at
`data/derived/audits/cutaneous-factory-ylro2d66/verification.json` (17.55 s,
initialization and confirmed cleanup only). The real 92-muscle/signed-native
session requested quadrature 0, received its canonical current position/normal,
zero indentation and source triangle 143278, and bound exactly the same material
identity into the pinned receptor. All clocks remained zero. The fixture's
right-postcentral assignment and gain 0.1 are explicit engineering inputs.
This verifies actual process/material identity wiring, not touch activation,
supported equilibrium or a coupled advancing whole-body trajectory.


A read-only integration review found that generic missing-contact release semantics
could hide missing native sensor rows. The regression now rejects missing,
duplicate and unrequested native samples, distinguishes control acknowledgments
from physical-state responses, and verifies failure occurs before physiological
mutation. Material receipts are copied so consumers cannot alter later identities.
