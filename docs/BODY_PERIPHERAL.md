# Peripheral somatic bridge

`BodyPeripheral` connects the existing IBM-derived regional brain to canonical
skin and the actual `BodyMechanics` actuator IDs. It is an executable generic
reference, not a calibrated complete peripheral nervous system. It does not add
an autonomic controller or drive respiratory muscles.

Build with `.venv/bin/python scripts/build_body_peripheral.py`; verify with
`.venv/bin/python scripts/verify_body_peripheral.py` and
`.venv/bin/python scripts/verify_body_brain.py`.

The generated `data/derived/canonical/peripheral.json` now declares 144 sided
nerve records and 20 relay groups, including all 71 IBM trunk names bilaterally
and the preserved bilateral `sciatic_fibular` route. New route declarations are
topology-only. See the [v2 IBM join contract](IBM_IHM_NERVE_JOIN.md) for endpoint
semantics, the companion bridge patch and fibre timing. The current build has
16 cutaneous patches, 249 actuator bindings and 715 unsupported actuator routes.
Multiple OpenSim compartments may bind the same canonical muscle. The builder
uses explicit anatomical name rules; it never equates the nearest nerve surface
with measured innervation. The unsupported set includes muscles with mixed supply,
unresolved compartments, and absent naming rules. It is deliberately not filled
by geometric guesses. The counts depend on the existing canonical anatomy and
mechanics artifacts.

`peripheral_display.json` has `lines` with `id`, `nerve_id`, `entity_ids`,
`points_m`, `kind`, and `evidence_kind`; `nodes` contains relays and patches with
`id` and `position_m`. Coordinates use the canonical body's meter frame, where
left has positive x. Centerlines are authored schematic routes through coarse
relays, not nerve fascicles or dissected trajectories. Patch centers are projected
to skin mesh vertices from authored regional landmarks. Radius and coverage are
visual priors. No canonical registry entity is created or replaced.

## Runtime contract

```python
peripheral = BodyPeripheral.from_dict(peripheral_data)
p = peripheral.step(
    0.02,
    stimuli={'peripheral-skin-left-palm': {'pressure_pa': 20000.}},
    mechanical_state=previous_mechanics_state,
    brain_state=previous_brain_state,
)
b = brain.step(0.02, physiology, sensory_inputs_hz=p['brain_inputs_hz'])
m = mechanics.step(0.02, {'activation': p['motor_activations']})
```

For sequential co-simulation, use outputs in the next physical interval; the
integration owner must retain and timestamp the previous states. This example
shows port wiring, not permission to erase an additional interval of latency.
The model returns its own `time_s` and never invokes or advances another solver.

`dt_s` must be 0.0001–1 second. Each call holds its supplied stimulus for the
interval; omitted stimuli return to zero pressure/stretch and 32 C skin baseline.
Patch keys are `peripheral-skin-{left,right}-{palm,forearm,upper_arm,thigh,calf,foot,face,trunk}`.
Supported stimulus units/ranges are pressure 0–1,000,000 Pa, temperature -20–80 C,
and fractional stretch 0–1. These are numerical interface bounds, not exposure
safety limits. Unknown IDs, nonfinite numbers, and unsupported modalities raise
before runtime mutation. Excess firing is capped at 200 Hz/channel. Zero input
means zero *excess* rate; biological spontaneous firing is omitted.

Each patch has pressure, stretch, warm and cold channels, with exponential
20 ms rate relaxation. Gains (0.005 Hz/Pa, 200 Hz/strain, and 8 Hz/C), saturation,
and 32 C reference are explicitly uncalibrated. Warm and cold travel separately
at 0.5 and 2.1 m/s; touch uses an assumed 50 m/s. A 12 ms central relay allowance
is added. Existing binding distances are 1.15 times the straight-line endpoint-to-relay
distance, with a 30 mm floor; central allowance is temporal and not a
measurement of cortical fiber length. These delays are not patient latencies.
Pressure channels are sustained low-pass responses; receptor subtypes, vibration,
rapid adaptation, nociception, axon spikes and injury are not represented.

Rate changes enter a timestamped delay queue. Sensory rates are sampled at the
end of each interval, then delayed, so coarser steps add up to one sampling
interval. There is no delivery before the physical delay. Motor commands are
sampled at interval start and activation is integrated exactly between motor
arrival events with a 30 ms time constant. Queued changes are suppressed below
1e-9 Hz/activation; no random generator, wall clock, or external process affects
stepping. Runtime work and memory are bounded by the finite channels, permitted
minimum time step, maximum rate and route delays.

Returned ports include `brain_inputs_hz` (contralateral postcentral regional IDs),
`motor_activations` (actual mechanics actuator IDs), `receptor_rates_hz`,
`proprioceptor_rates_hz`, `relay_activity_hz`, `nerve_activity_hz`, and
`pending_events`. Relay/nerve activity is an aggregate visualization of arrived
rates, not an action-potential recording. The coarse pathway does not resolve
DCML versus spinothalamic tracts, thalamic nuclei, or decussation locations.

Brain sensory input adds 0.02 nS per excess afferent Hz to the IBM-derived neural
conductance drive outside the preserved functions. Regional cortical activity
has no identified muscle recruitment law. The peripheral runtime recruits
muscles only through explicit motor commands; sensory input alone does not
invent a motor policy.

Explicit `brain_state.motor_commands={actuator_id: fraction}` can command any
supported binding at 0–1. Motor conduction reads the binding alpha-class delay plus a separate 12 ms
central prior. Proprioceptive transport reads the Ia-class delay. Velocities
come from a preserved IBM table snapshot; scalar delays are compatibility aliases
for these specific classes. Gamma timing is separate, without a gamma controller. Commands
omitted in a later interval return toward zero after conduction and activation
relaxation. A normal resting brain state cannot spontaneously command muscles.

Proprioception combines positive axial affine strain, positive change in the
mechanical attachment path, and tendon force normalized to maximum isometric
force. The reduction is `min(200, 200*strain + 40*force/Fmax)` Hz, then 20 ms
relaxation and 60 m/s afferent transport plus 12 ms central allowance. This lumps
spindle and tendon-organ excess signals; it omits fusimotor control, shortening
responses and motor-unit recruitment. Returned mechanics transforms and forces
are used, so the feedback is caused by actual mechanical output. The attachment
path reconstruction uses translations and does not resolve local soft-anchor
rotations. It is not a detailed proprioceptive organ simulation.

## What is reused from IBM-1

The existing `BodyBrain` still executes the exact ASTs of IBM's `_sigmoid`,
`wilson_cowan_excitatory`, and `shunting_inhibition_rate` from hash-verified
snapshots. Its 80 regional populations and transferred 68 DK cortical labels
remain unchanged. Sensory conductance and somatic readout are new IHM priors.
Original IBM neural functions were not modified.

IBM's `ibm/processes/transduction.py`, `ibm/processes/effector.py`, and
`ibm/runtime/step.py` were inspected and preserved byte-for-byte under
`data/derived/canonical/peripheral-sources`, with hashes in `source_receipts`.
Their role is `inspected_reference_not_executed`. IBM declares receptor filters,
activation dynamics and descending transport, but its current runtime solves
spectral trajectory windows with overlap-based causality checks. This bridge
uses a different time-domain delayed-rate reduction. In particular, its single
30 ms activation pole is not the original IBM two-pole rise/fall transfer.

The whole IBM multimodal materialization, spectral solver, uncertainty
representations, data fitting, EEG/MEG forward models, sensory epithelia and
learned motor control are **not** copied or run by this peripheral addition.
Calling this integration “the whole IBM brain” would be inaccurate.

## Evidence receipts and limits

The artifact retains source URLs and scoped findings; the preserved IBM files
retain local SHA-256 receipts. Literature was checked September 5, 2026.

- [Human fingertip afferent recordings](https://pmc.ncbi.nlm.nih.gov/articles/PMC6666033/)
  support force-related sensory firing and timing. The implemented scalar gain
  and regional patch layout are not fitted to those recordings.
- [Human warming/cooling conduction study](https://pubmed.ncbi.nlm.nih.gov/3225599/)
  reports mean estimates 0.5 and 2.1 m/s. Applying these means to every authored
  regional route is a population transfer assumption.
- [Human muscle spindle stretch responses](https://pubmed.ncbi.nlm.nih.gov/2141632/)
  support stretch sensitivity and show heterogeneous response classes. The
  implemented positive-strain/tendon-force sum is a simplification.
- [Median motor branches in cadaver forearms](https://pmc.ncbi.nlm.nih.gov/articles/PMC7580294/)
  and [anterior arm innervation variants](https://pmc.ncbi.nlm.nih.gov/articles/PMC4079987/)
  support named anatomical relationships and document exceptions.
- [Human lower-limb segment stimulation](https://pubmed.ncbi.nlm.nih.gov/1766452/)
  demonstrates broad overlapping roots and frequent variants; one regional relay
  per nerve is not a validated segmental assignment.
- [Vastus medialis nerve dissections](https://pubmed.ncbi.nlm.nih.gov/2254163/)
  support femoral/lumbar supply and more complex compartment patterns than the
  present body reduction.

Other named anatomical rules, 50/60 m/s velocities, central delay, geometric
path lengths and cortical gain/readout are explicit generic priors, not values
measured by these sources. There is no reported biological validation or
calibrated confidence percentage. Verification checks software causality and
resolved identities rather than asserting physiological accuracy.
