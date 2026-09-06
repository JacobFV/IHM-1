# Observing an external intake mass boundary

`ihm.assembly.intake_mass.consumed_intake_delta` compares two authoritative
`ihm.native-consumed-intake.v1` envelopes. It emits a boundary only when the
native consumed count advances by exactly one. A queued meal, repeated snapshot
or internal GI transfer produces no new mass. Epoch changes, skipped counts,
changed counters without consumption, inconsistent components, and incompatible
clocks or command ordering reject.

The payload mass is the actual native `SENutrition::GetWeight` result, checked
against its retained source convention: nutrient grams plus calcium milligrams
and water milliliters converted to kilograms. The water conversion is the
native model's convention; it is not a new density measurement. Source audit:
`NATIVE_INTAKE_MASS_BOUNDARY_AUDIT.md`.

The output identity is `[native owner epoch, consumed count, meal sequence]`.
This function is a pure observer. It does not own mass, apply an impulse or
provide exactly-once mechanical delivery by itself. The mechanical owner must
bind that identity to its own acknowledged transaction, preserve momentum and
separate mass-flux/capture energy, and treat uncertain delivery as terminal.
Historical consumed quantities cannot be silently replayed into a newly attached
mechanical owner. Internal absorption and excretion require their own ownership
and boundary accounting.

Five focused tests cover consumption, queued/no-op observations, cumulative
meals, malformed identity/count/component increments and clock rejection. These
are contract tests; actual adapter envelope acceptance is recorded separately
when the native receipt implementation runs. No mechanical/native integration
claim follows from the parser tests.

## Endpoint mechanical bridge

`IntakeMassBridge` binds a fresh zero-count native intake owner to a fresh,
explicitly enabled mechanical mass owner. Its target torso station and
registration identity must be supplied. It queries that exact station's native
source-frame velocity and records the explicit co-moving-at-ingestion assumption.
This is not measured swallowing momentum. Transfer occurs at the matching
mechanical interval endpoint without advancing time or revising muscle energy.

It validates mechanical epoch, sequence, enabled state, individual owner/site
inventory, cumulative mass and returned boundary mass before recording success.
The validated native port range is at most 0.5 kg of added payload in total;
the bridge does not split larger intakes or infer excretion. A zero-mass native
consumption advances boundary accounting without issuing a mechanical impulse.

An already-consumed boundary is retained before mechanical preflight. Any
subsequent failure is terminal for this bridge, including capacity, clock or
point-query failure. Uncertain mass commands are never retried. Snapshots retain
pending and acknowledged boundary identities but are not native checkpoints.
The surrounding whole-body owner must abort if physiology has consumed intake
that mechanics cannot account for.

Contract tests use the actual production port shape with a finite fake owner.
Review counterexamples cover wrong point/site velocity and missing/disabled
owner acknowledgment; both reject. Actual native coupled acceptance and embodied
factory activation remain required before this bridge can be claimed as live
whole-body intake feedback.

### Actual endpoint acceptance

`data/derived/audits/intake-mass-bridge-ho3xdh6h/verification.json` records two
actual native water consumptions (10 and 20 mL) projected into the actual
92-muscle mechanical model at the hash-bound canonical stomach centroid on the
torso. Total mechanical mass increased by 0.03 kg. Both common endpoints, exact
coordinate values, frozen muscle reference/accumulated ledgers, constraint and
momentum residuals, explicit owner inventory, and no-op receipt reobservation
passed. Both processes were reaped; wall time was 1.90 seconds.

Run: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 prlimit
--as=2147483648:4294967296 nice -n 10 .venv/bin/python -m
scripts.verify_intake_mass_native --run` with the shared heavy slot reserved.

This is a free-body, two-owner mass-boundary acceptance. Physiology used its
ordinary native advancement; it does not validate the whole-body signed
metabolic loop, supported equilibrium, swallowing mechanics, regional organ
deformation or excretion. The receipt preserves that distinction.
