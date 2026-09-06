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
