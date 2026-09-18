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

## Route-length audit, 18 September 2026: the frame is right; the relays were not on the cord (FIXED the same day, below)

**Why it was run.** IBM-1 withdrew a claim on 17 September that the two body atlases share a
frame. The Z-Anatomy display export is a normalised ±1 box and the canonical body is metres,
15.7% apart. Every IBM conduction delay is `path_length_m / velocity` over the lengths in
`peripheral.json` and `dermatomes.json`, so a route measured on the normalised atlas would
make every delay on it 1.157x wrong. The audit script is IBM-1
`scripts/audit_nerve_route_lengths.py` and its output is IBM-1 `out/nerve_route_audit.json`.
The full argument is in IBM-1 `docs/LOG.md`, 18 September.

**No length is in the wrong frame.** The evidence:
- `peripheral.json`, `dermatomes.json` and `anatomy.json` declare `bodyparts3d-display-m`.
- The builders read positions only from `anatomy.json` entity centroids, the FJ2810 skin
  mesh, and authored coordinates. None of them reads the extended manifest.
- The frame is adult-sized in metres: skin 1,719 mm tall, femur 467 mm, tibia 380 mm,
  humerus 310 mm.
- Z-Anatomy entities used as endpoints entered `anatomy.json` through the landmark
  registration, held-out RMS 4.7 mm, from blend world coordinates.
- A re-implementation of this repo's length method reproduces all 146 declared nerve lengths
  exactly.
- No common 0.865 factor appears anywhere.

**What is wrong is the relays.** They are authored at `[±0.012, y, −0.04]`
(`build_body_peripheral.py`, `ys`) as `regional_group_prior`. Each spinal relay sits below
the cord segment it stands for. The medulla-level relays are the exception, and they are
close.

| relay | y | anatomical level | y there | off |
|---|---|---|---|---|
| cranial / solitary | 0.650 / 0.655 | medulla oblongata (FJ1769) | 0.698 | −48 / −43 mm |
| cervical | 0.530 | C5–T1 segments (C5 vertebra) | 0.619 | −89 mm |
| thoracic | 0.300 | mid-thoracic segments (T6) | 0.452 | −152 mm |
| lumbar | 0.040 | lumbar segments (T11) | 0.316 | −276 mm |
| sacral | −0.080 | sacral segments, conus (L1) | 0.257 | −337 mm |

Two authored endpoints are also off the structures they name:
- The vagus "gastric wall" endpoint `[.06, .18, .03]` in `peripheral_route_catalog.py` is
  132 mm below the stomach centroid and 73 mm below its lowest vertex.
- The phrenic "diaphragm" endpoint is 72 mm below the diaphragm centroid.

The same `[.06, .18, .03]` is used for the greater splanchnic.

**Known answers.** The references are Z-Anatomy nerve centrelines taken from the staged
`Startup.blend` and carried to metres. That transform lands the blend's femur on the canonical
femur.
- Medial plantar, conus → sole: 1,207 mm along the centrelines, against 879 mm here.
- Tibial, from the tibialis posterior to the root entry: 876 mm, against 649 mm here.
- Femoral: 515 mm, against 333 mm here.
- Median, from the pronator quadratus: 660 mm, against 610 mm here.
- Right vagus, brainstem → stomach: 415 mm (first contact) to 593 mm (lowest), against
  508 mm here.
- Phrenic (literature): right 246–300 mm, left 306–330 mm, against 335 mm here. That match
  is two errors cancelling (relay 89 mm low, endpoint 72 mm low), not a measurement.

Keeping this repo's method (`1.15 x` straight line, authored polylines) and moving only the
relays onto the cord reproduces the centrelines to within 5–9%. **The `1.15` factor is
adequate; the relay positions are the error.**

Moving the relays changes each route by x0.53 to x3.49, median x1.27:
- legs x1.4–2.9;
- splanchnics x1.6–3.5 (greater 168 → 314 mm, pelvic 128 → 446 mm);
- vagus x1.08 (508 → 551 mm);
- arms x1.1–1.45;
- cranial motor routes x0.5–0.9;
- special senses unchanged.

`dermatomes.json` routes every patch through the same relays, so all 1,326 patch lengths
inherit the bias.

**Not changed by the audit** (superseded by the fix below). Moving the relays alters every route, every patch, the stature-scaling
recomputation and every downstream delay at once. It is a contract decision for this repo,
not a correction to slip in. Until it is made, spinal-route lengths are schematic lower
bounds, and the lumbosacral and splanchnic ones are short by 1.4–3.5x. The records already
say `geometry_kind: schematic_route` and `measured_axon_geometry: false`. The finding is that
the schematic is biased, not merely unvalidated.

**A related seam in stature scaling** (fixed below). `scale_nerve_conduction.py` scaled route lengths by
`stature / MECHANICAL_STATURE_M` (1.7973 m). The lengths were measured on the anatomical body,
`ANATOMICAL_STATURE_M` = 1.7195 m. A "2.03 m" variant therefore carries the routes of a
1.942 m body, 4.5% short. `STATURE_DISAGREEMENT` records the two statures, but the nerve
scaling does not cite it.

### FIXED, 18 September 2026: the relays are on the cord, four organ endpoints are on their organs

**What moved, and where every coordinate now comes from.**
`scripts/build_spinal_cord_levels.py` writes `data/derived/canonical/spinal_cord_levels.json`;
`build_body_peripheral.py` and `enrich_peripheral_routes.py` read the relays from it and refuse to
build without it. Each relay record carries `position_source`.

- **The cord's path** is the centreline of Z-Anatomy's `Spinal dura` mesh (the object IBM-1
  `scripts/blender_zanatomy_nerves.py` extracts from the staged `Startup.blend`), carried into
  `bodyparts3d-display-m` by the same landmark registration that placed every `body-za-*` entity
  in `anatomy.json` (99 bones, affine + thin-plate spline, held-out RMS 4.7 mm). It is not the site
  exporter's bounding-box chain; the two differ by 9.8 mm at C5 and under 4 mm below it.
- **Segment levels** are BodyParts3D vertebral bodies: the midpoint of the discs above and below.
- **Which vertebra a segment sits behind** is an authored clinical prior, stated per relay.
- **The brainstem relays** are the BodyParts3D pons (`cranial`) and medulla oblongata (`solitary`)
  halves, symmetrised.

| relay | stands for | old position | new position (m) | level | moved up by | to dura centreline | vertebral canal at level |
|---|---|---|---|---|---|---|---|
| cervical | C1–T1 (brachial C5–T1, cervical plexus, phrenic) | `[±.012, .530, −.040]` | `[±.004, .6192, −.0289]` | C6 segment, C5 body | 89 mm | 4.0 mm | inside: 24/24 directions, 4.0 mm clear |
| thoracic | T1–T12 | `[±.012, .300, −.040]` | `[±.004, .4855, −.0685]` | T7 segment, T5 body | 186 mm | 4.0 mm | inside: 24/24, 6.4 mm |
| lumbar | L1–L4 | `[±.012, .040, −.040]` | `[±.004, .3270, −.0578]` | L3 segment, T11 body | 287 mm | 4.0 mm | inside: 23/24, 3.8 mm |
| sacral | L4–S4 | `[±.012, −.080, −.040]` | `[±.004, .2686, −.0408]` | S2 segment (conus), L1 body | 349 mm | 4.0 mm | inside: 18/24, 4.7 mm |
| cranial | III–VII, XII nuclei | `[±.012, .650, −.040]` | `[±.0105, .7296, −.0210]` | pons | 80 mm | — | — |
| solitary | nucleus of the solitary tract | `[±.012, .655, −.035]` | `[±.0067, .6984, −.0329]` | medulla oblongata | 43 mm | — | — |

"Moved up by" is the change in y. The thoracic relay sits 33 mm higher than the audit's T6 figure
because a relay standing for T1–T12 belongs at its middle segment, T7, which lies behind the T5
body. Each left/right pair is an exact mirror image. The centreline's own x (under 0.3 mm) is dropped,
and the brainstem halves' centroids are averaged in y and z, because IBM-1's visceral join collapses
the two sides and raises if their lengths differ by 1e-9 m.

**Endpoints re-anchored** (`ANCHORED` in `peripheral_route_catalog.py`). Every authored endpoint
whose label names an atlas structure was put beside that structure. The trunk-and-viscera points sat
10–17 cm low, like the relays did. The ones re-anchored onto a vertex of the named mesh:

| route | label | was | now |
|---|---|---|---|
| vagus | gastric wall | 87 mm from the stomach surface | anterior stomach wall in front of its centroid |
| greater splanchnic | foregut | the same point as the vagus | the vagus's wall point |
| least splanchnic | renal | 101 mm from the left kidney | left kidney hilum (surface vertex nearest the midline) |
| pelvic splanchnic | pelvic viscera | 171 mm below the bladder, in the thigh | bladder wall vertex nearest the authored point |
| phrenic | diaphragm | 29 mm off the diaphragm, 72 mm below its centroid | the dome's highest vertex at the authored lateral offset, per side |

**The visceral right side is a mirror, not anatomy.** The right vagus really reaches the
*posterior* gastric wall, and the right kidney sits lower than the left. Anchored per side, the vagus
is 448 mm on the left and 423 mm on the right. IBM-1 `ihm_bridge.visceral_routes` then raised
("the sides can no longer be collapsed"), and four interoception tests and one cortical-sheet test
failed. So each visceral right route is declared as the mirror of the left. Its record says
`endpoint_source.rule: mirror_of_left` and `on_named_structure: false`, and it carries no
`endpoint_id`. Lifting this needs the brain-side join to accept per-side lengths.

Other authored endpoints were checked and left alone. Muscle-region labels fell 8–41 mm from the
named muscle; `lower_subscapular` is 84 mm from subscapularis, but that nerve also supplies teres
major. The generic visceral labels ("upper abdominal", "lower abdominal", "thoracic visceral")
name no single structure.

**Instrument.** `scripts/verify_relay_placement.py` re-implements this repo's length method.
- It reproduces all **146/146** new declared lengths (worst 1.1e-16 m).
- With the old relays and authored endpoints put back, it reproduces **146/146** lengths of the
  pre-change file (worst 5.6e-17 m).
- IBM-1's audit script, run against the new file, also reproduces 146/146.

**Known answers.** Z-Anatomy centreline measured from the point nearest the route's own endpoint to
the cord at the relay's level, in this frame. The target was within 10%.

| route | before (mm) | after (mm) | Z-Anatomy (mm) | before/Z | after/Z | verdict |
|---|---|---|---|---|---|---|
| tibial (L / R) | 649 | 1,049 | 958 | 0.68 | 1.09 / 1.10 | within |
| deep fibular | 790 | 1,189 | 1,100 | 0.72 | 1.08 | within |
| superficial fibular | 743 | 1,143 | 1,050 | 0.71 | 1.09 | within |
| median | 610 | 705 | 679 | 0.90 | 1.04 | within |
| medial plantar (sciatic–tibial, to the sole) | 879 | 1,279 | 1,211 | 0.73 | 1.06 | within |
| femoral | 333 | 656 | 597 | 0.56 | 1.10 | Z stops 57 mm short of the endpoint: lower bound |
| radial | 562 | 658 | 482 | 1.16 | 1.36 | Z ends at the elbow, 156 mm short: lower bound |
| vagus, left (declared) | 508 | 448 | 252 | — | — | Z left vagus ends mid-thorax (y 0.484): no reference |
| vagus, right, anchored on the posterior wall (not declared) | 508 | 423 | 491 | 1.03 | 0.86 | **OUTSIDE** |

Of the 11 same-scope comparisons, **10 are within 10% now, against 1 before**. The one outside is
the vagus. It is a polyline route, and the polyline method has no `1.15` tortuosity factor
(`1.15 × 423 = 486` against 491). The method was not changed to close it. The "Z stops short"
rule (the nearest Z point is the curve's own distal end) was written after the first table showed
radial and left vagus. Both rows are still printed.

**Over a metre.** Conus to sole along Z is 1,219 mm. IHM medial plantar is **1,279 mm** (was 879),
the sole receptor patch 1,231 mm, and the longest dermatome patch (L5, foot) 1,202 mm. **PASS.**

**Every route, old → new.** The ratio runs from ×0.54 to ×3.18, median ×1.33 (left side, full list
in `relay_placement_report.json`).
- Legs ×1.4–3.0: sciatic_tibial 201 → 600 mm, femoral 333 → 656, tibial 649 → 1,049.
- Splanchnics: greater 168 → 279, lesser 217 → 411, least 232 → 258, lumbar 129 → 404,
  pelvic 128 → 301.
- Arms ×1.1–1.5.
- Vagus 508 → 448: the endpoint moved up onto the stomach.
- Phrenic 335 → 284 (left) and 256 (right), against the literature's right 246–300 and left
  306–330 mm.
- Cranial motor routes ×0.5–0.8, because the pons is nearer the orbit than the old relay was.
- Special senses unchanged.

The 1,326 dermatome patches changed only `path_length_m`: median ×1.29, foot 819 → 1,165 mm,
thigh ×2.05, head ×0.71.

**Gates re-run, and every verdict that changed.**

| gate | before | after |
|---|---|---|
| `build_spinal_cord_levels.py`: frame check, extraction = registration's blend world | — | PASS (1.7e-8 m) |
| same: canal sweep C3–L1, bar fixed before first run | — | **FAILED 17/18**: C3, 1.5 mm clear against a 2 mm bar. Relays are gated on their own four levels, which pass. Under the 15.7% frame error the sweep passes 8/18, so it can fail, but it is weak in the thoracolumbar canal |
| `scale_nerve_conduction.py --scale 1.0`: 30 mm floor inactive | PASS | **FAILED**. The right sternocleidomastoid (accessory) binding is now 29.6 mm by straight line from the C6 relay. The real accessory nerve goes up through the foramen magnum and out the jugular foramen, and a straight line cannot represent that |
| `scale_nerve_conduction.py --stature-m 2.03`: independent arm, bindings from scaled geometry | PASS | **FAILED** (1.3%): the same floored binding, since `max(0.03, ·)` is not homogeneous |
| `materialize_body_variant.py --stature-m 2.03` (into scratch) | 48/48 | **48/49**. The variant was NOT written over `data/derived/body-variants/stature_1p129489`, which still carries the pre-fix nerves |
| `materialize_body_variant.py --scale 1.0` (into scratch) | 52/52 | **50/51**. The identity gate no longer applies, because the route scale there is 1.0452 (below) |
| `verify_peripheral_join.py` | FAILED (expects 144 nerves; 146 since `dorsal_ramus`) | FAILED, unchanged |
| `verify_provenance_pins.py` | FAILED (`body.json` pinned peripheral.json bytes that were already stale) | FAILED, unchanged |
| `verify_body_peripheral`, `verify_dermatome_patches`, `verify_peripheral_coverage`, `verify_sensorimotor_loop`, `verify_body_integration`; `scale_skin_patches` at 1.0 and 2.03 | — | PASS |
| IBM-1 `tests.test_interoception`, `test_cortical_sheet`, `test_innervation_coverage`, `test_cord` | — | PASS; they failed on the per-side visceral version, above |

The floor failure is not rescored. The fix it points at is to route the accessory nerve through the
jugular foramen as a polyline, which is a change to the method, and it is left to this repo.

**Stature seam, fixed.** Nerve and skin-patch routes were measured on the anatomical body, whose
skin is 1.7195 m. They now scale by `stature / ANATOMICAL_STATURE_M`:
- `resolve()` gains `derived['anatomical_stature_scale']`.
- `scale_nerve_conduction.py` and `scale_skin_patches.py` use it.
- `materialize_body_variant.py` passes it to both stages.

A new nerve-stage gate re-measures the skin mesh's extent. It requires that extent to equal the
reference, which the old 1.7973 m reference fails by 4.5%. At 2.03 m every delay is now **+18.06%**,
not +12.95%.

The quoted "vagal C 507.84 → 573.60 ms" becomes:
- **599.5 ms**, if only the reference changed;
- **447.91 → 528.80 ms** with the relays and endpoint fixed too.

The anatomy stage still scales the displayed entities by the mechanical factor, because it
preserves the declared registration between the two bodies. So a 2.03 m variant displays a
1.942 m anatomy while its routes belong to a 2.03 m body. `registration.json` records this as
`known_seam`. The unscaled canonical `peripheral.json` is the 1.7195 m body's; the default mechanical
stature, 1.7973 m, now gets routes ×1.0452.

**Still open.** See the limitations in `spinal_cord_levels.json`:
- One relay stands for a group, so C2 scalp patches are routed to the C6 level.
- The lesser splanchnic (T10–T11) is routed to T7.
- The segment-to-vertebra relation is a population rule, ±1 segment.
- The dura is Z-Anatomy's, not this body's.
