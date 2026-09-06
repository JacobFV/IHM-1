# Signed cardiovascular reader and separate resistance experiment

## Exact source finding

The retained signed Cardiovascular.cpp is SHA-256
`ad93e3930f307f487e166dac1795234074f73062fc11b44eead288fe5d305f6d`, inherited
from `whole_body_integrity_signed_muscle_v2`. At lines 2083–2096,
`AdjustVascularTone()` calls `MetabolicToneResponse()` only inside
`abs(ResistanceChange)>ZERO_APPROX`. ResistanceChange is the drug-induced MAP
change in mmHg divided by cardiac output in mL/s, then multiplied by the native
tuning factor. It has units mmHg s/mL. Ordinary non-drug intervals therefore
skip the signed effective reader in `MetabolicToneResponse()` (reader index 0,
line 1656).

A separate gate inside the tone method, lines 1682–1686, enables the resistance
modifier only for generic `HasExercise()`. Signed muscle coupling explicitly
rejects generic Exercise. Moving the tone call outside the drug branch alone
would consequently expose the reader and update an Exercise MAP target, but
would not cause direct muscle/peripheral resistance feedback. That target is
read in `Nervous::BaroreceptorFeedback` (Nervous.cpp line 553). It is initialized
in Energy.cpp line 124 and is not routinely reset there. The tone method's early
return at metabolicFraction==1 can leave an old target in place. A naive
nonzero-only call would thus have a release latch and could jump from an old
drug-gated target when current background metabolism differs from basal.

HeartDriver calls AdjustVascularTone every active-heart step (line 1456), not
only at heartbeat boundaries. Cardiac arrest skips this call. Native controller
order is Cardiovascular PreProcess before Energy PreProcess, Endocrine, Drugs,
Tissue and Nervous; the cardiovascular reader therefore sees the preceding
native background TotalMetabolicRate plus the current boundary's signed delta.
This phase difference must not be mistaken for an inconsistent unit conversion.
The reader adds delta watts multiplied by `Convert(1 W, kcal/day)` to native
kcal/day and records its result back in watts. MetabolicFraction divides by
patient basal kcal/day; it is dimensionless.

## Observation-only correction

`patch_signed_cardiovascular_reader.py` produces a `reader_only` overlay adding
an else branch to the existing nonzero-drug-resistance conditional. With an
active signed record and available native total metabolic rate it calls the
same effective reader once and discards the returned value. The drug branch
retains its existing call. No native MAP target, resistance, chemical store or
heat source is written. Inactive physical behavior and signed zero-demand
physical behavior are preserved on valid native input. Reader counters now
record ordinary active signed intervals, including zero; this is instrumentation,
not completion of physiological feedback. Invalid negative effective demand
can be rejected at this earlier reader, without any rollback claim.

## Separate transferred resistance prior

The distinct `vascular_prior` overlay additionally applies a transparent
normalization of the original native exercise rule:

- Native exercise q(F) = (1.5 F + 5.5) / 7.
- B is current native background total metabolic rate, P is patient basal rate,
  and D is the signed muscle chemical increment; all three are in watts.
- Relative ratio = q((B+D)/P) / q(B/P) = 1 + 1.5 D / (1.5 B + 5.5 P).
- Multiply muscle NextResistance by ratio^-2 and extrasplanchnic/splanchnic
  NextResistance by ratio^-1 after native autonomic reset and drug handling.

D remains the existing command's increment relative to its frozen mechanical
muscle reference. B is the engine's current background, not another cumulative
chemical owner or a new frozen basal demand. Normalization makes D=0 exactly
one even when B differs from P. Positive D lowers resistance. Negative D raises
it without rectifying or clipping D, with the explicit mathematical domain
B>=0, P>0, B+D>=0 and finite inputs. This is an extrapolation of an exercise
prior to negative muscle demand, not a validated physiological response.
The existing native muscle substrate/decrement budget still applies independently; this ratio grants no additional chemical availability. The inherited 0.1 mmHg s/mL resistance floor remains in force and is reported
as part of the inherited model domain, not a newly fitted limit.

The experiment requires all three native regional autonomic resistance scales
so the preceding native code resets each affected path from its baseline each
active-heart step. Cardiac arrest skips the inherited tone dispatcher and is outside this release experiment. Otherwise it rejects, avoiding cumulative multiplication. At D=0 it makes
no resistance write; the existing native reset removes the local modifier.
A release test should show local before/after equality, not equality to the
whole-body trajectory before a perturbation: physical pressure, flow and neural
states may retain causal consequences. No Exercise MAP target is changed, and
no energy/heat amount is reassigned. A separate diagnostic trace records local
pre/post resistance sums and factors; this trace adds no physiological store.

The numbers 1.5, 7, powers -2/-1 and resistance floor are transferred from the
pinned native exercise law. The anatomical distribution and magnitude are not
calibrated for signed active-fiber demand. This experiment does not resolve
cardiovascular conditioning, local oxygen control, disease domains, load-state
nervous defects, mechanical energy-supply limitations or supported equilibrium.

## Reproducible isolated validation

Five source/algebra checks run with `verify_signed_cardiovascular_source.py`.
The preparer emits separate overlays with pinned input/output/header hashes;
`build_signed_cardiovascular_experiment.py --build` compiles one isolated
Cardiovascular translation unit and links a fresh explicit variant. It preserves
all inherited object identities and changes no default library, adapter or pointer.
Native compilation requires the coordinated heavy slot.

`verify_signed_cardiovascular_native.py --run --reader <variant> --prior <variant>`
compiles a separate probe and runs parent/reader/prior libraries serially under
resource caps, using detached retained runtime inputs. The planned checks are
inactive/zero physical parity, all-mode reader-only physical parity, one
cardiovascular reader and one tissue/heat owner, positive and negative direct
vascular response, and zero local modifier release. Native acceptance has not
been run at source preparation time; any actual receipt must be reported
separately. No production promotion is implied.

## Native falsification (2026-09-06)

Both isolated libraries compiled: reader-only `whole_body_integrity_cardio_reader_v1`
in 3.126 s / 571,988 KiB maximum compiler child RSS, and prior
`whole_body_integrity_cardio_vascular_prior_v1` in 3.179 s / 574,260 KiB.
The 12 serialized three-step native sessions completed at
`data/derived/audits/signed-cardiovascular-native-oevybgqr/frames.json`, but the
intended vascular response assertion **failed**. The failure is retained in
`falsification.json`; there is no accepted vascular response receipt.

The reader-only correction recorded one cardiovascular read per signed step
and preserved all observed physical values for all four modes. Disabled/zero
physical values matched the parent exactly. At the first endpoint the native
background and basal rate were both 82.18220139584905 W; reader values changed
by the requested +/-0.01 W. The experimental ratio moved to
1.0000260744675424 or 0.9999739255324577. Nevertheless the muscle resistance
sum in its cached path vector remained zero, and Aorta1ToMuscle1 resistance
remained 3.92438997821351 mmHg s/mL. Only the manually added portal-vein
splanchnic entry existed; its 0.1 floor suppressed positive-demand dilation.

This exposed a second, separate defect: `Cardiovascular::SetUp` fills its
regional resistance vectors only from paths with `HasCardiovascularRegion()`.
The existing Circuit.xsd already defines that optional field, but
`io/cdm/Circuit.cpp` neither writes nor reads it. SEFluidCircuitPath initialization
sets it Invalid, and the legacy retained state contains no such labels.
Consequently routine regional autonomic feedback is also lost on this LoadState
path. Repair requires correct IO and an explicit source-defined topology
reconstruction for legacy snapshots, not additional metabolic scaling.

The observer initially assumed a venous path had resistance; retained XML
showed Muscle2ToVenaCava has flow only. The final raw probe keeps its actual
flow and omits an invented resistor; `falsification.json` explicitly records
resistance as unavailable. Early compile/resource/probe assumptions are retained
in q81sp5u6, pm4ongn3 and gqnu4i5n audit directories. These were fixture failures,
not evidence of a physical response.
