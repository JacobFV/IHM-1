# Native respiratory control feasibility audit

Audited 2026-09-05. **A causal blood-gas-controlled respiratory actuator is feasible in the retained native engine. An anatomically resolved brainstem–phrenic–diaphragm pathway is not already present.** First expose and perturb the existing native controller; then introduce an explicitly owned neural/actuator replacement through a native pre-solve port. Do not label a pressure multiplier or a replayed diaphragm animation as that neural pathway.

## Retained implementation and build

The donor is `data/raw/physiology/biogears`, repository `https://github.com/BioGearsEngine/core.git`, clean at revision `3f16a5fa1dade9c511b88d923606fa51cc35e95d`. The acquisition receipt is `data/derived/physiology/collection.json`. All source locations below are relative to `projects/biogears/libBiogears/` within that donor.

`scripts/build_native_biogears.py` builds an ARM aarch64 native library and C++20 adapter using the local extracted Eigen/Xerces/log4cpp/CodeSynthesis dependencies. A temporary C++ probe compiled and linked `SEApnea::GetSeverity`, `PhysiologyEngine::ProcessAction`, respiratory driver getters and the arterial CO2 getter against the existing library; its version-only execution printed `8.0.0+Source-Unknown`. This verifies API availability, not physiological intervention behavior. Generated version metadata cannot identify the donor accurately; preserve the revision and byte hashes separately. No native intervention simulation was run for this audit.

For a source extension, follow the isolated translation-unit replacement and link-response-file pattern in `scripts/build_biogears_saturation_variant.py`: retain original donor bytes, exact patch and source hashes, and create a separately named library manifest. Compose with the already selected saturation/thermal fixes rather than silently reverting them. A new public class layout/API needs consistent headers and recompilation of dependent translation units; replacing only `Respiratory.cpp.o` is insufficient if member layout changes.

## Existing causal chain and timing

| Source anchor | Actual behavior |
| --- | --- |
| `src/engine/Systems/Nervous.cpp:239,274,803` | `PreProcess → AfferentResponse → ChemoreceptorFeedback` reads instantaneous native arterial PaO2 and PaCO2. Peripheral afferent firing has a 2 s first-order response; central frequency/pressure contributions have 180 s constants and peripheral contributions 100 s constants. Central input is arterial CO2 deviation, with lower gain below baseline; peripheral input combines O2 and CO2. Drug and metabolic modifiers also apply. |
| `Nervous.cpp:888–919` | Publishes baseline plus central/peripheral frequency changes to `RespirationDriverFrequency` (/min), and baseline minus central/peripheral pressure changes to `RespirationDriverPressure` (cmH2O). This pressure is a requested negative amplitude, not the instantaneous circuit source. |
| `src/engine/Systems/Respiratory.cpp:423,817,927` | Preprocessing generates a phase-dependent muscle-pressure waveform. At cycle boundaries it latches requested amplitude/frequency and applies driver actions, limits and drug effects. Apnea, ventilator state and no-breathing state subsequently modify the waveform. |
| `Respiratory.cpp:454` | `m_Calculator.Process(RespirationCircuit, m_dt_s)` solves native pressures, flows and volumes, followed by gas/aerosol transport and calculated respiratory outputs. |
| `src/engine/Systems/Tissue.cpp:525,838,1942` | Native alveolar–pulmonary-capillary substance transfer and tissue diffusion continue; blood chemistry then updates blood-gas outputs. Gas exchange is not solely owned by `Respiratory`. |
| `src/engine/Controller/BioGears.cpp:971–1004` | Respiratory preprocessing precedes nervous preprocessing. Thus the waveform consumes the prior published nervous request, and its amplitude/frequency also wait for a cycle boundary. The next nervous update sees the previously completed blood-gas state. Preserve and record this discrete timing. |

The donor [nervous methodology](https://www.biogearsengine.com/documentation/_nervous_methodology.html) describes aggregate independent central/peripheral respiratory feedback and assumes rapid blood–brain CO2 equilibration. Its [respiratory methodology](https://www.biogearsengine.com/documentation/_respiratory_methodology.html) describes the muscle pressure source, circuit and gas transport. The retained source is authoritative where public documentation names older target-ventilation interfaces.

There is no spiking respiratory rhythm generator, identified pre-Bötzinger network, phrenic conduction, motor-unit recruitment, neuromuscular junction or diaphragm force–length law in this chain. `Nervous::CentralSignalProcess` also controls autonomic cardiovascular signals; it is not an anatomical respiratory central pattern generator. Pulmonary stretch feeds autonomic processing here; its existence does not establish a native Hering–Breuer respiratory phase controller.

`ihm/assembly/body_runtime.py:53–61` prescribes thoracic/diaphragm motion from native gas-volume telemetry. Its brain input includes MAP, saturation and temperature, not PaCO2; its neural output does not drive the native engine. The generic medulla activity and uncalibrated autonomic readouts in `ihm/assembly/brain.py` do not supply the missing pathway.

## Supported interventions and exact extension

**Immediate experiment:** add an adapter timeline action using `SEApnea` (`include/biogears/cdm/patient/actions/SEApnea.h`), severity 1 for drive removal and severity 0 for restoration, through `bg->ProcessAction(action)`. `Respiratory::Apnea` at line 1754 scales

`P_applied = P_default + (P_generated - P_default) * (1 - severity)`.

Severity zero removes the action in `src/cdm/scenario/SEPatientActionCollection.cpp:441`. Nervous feedback continues accumulating while the actuator is suppressed. This is a lumped actuator intervention, not selective phrenic blockade. `P_default` initializes to **−5 cmH2O** and is serialized; the circuit path's baseline is **0**, so neither zero absolute source nor `GetPressureSourceBaseline()` is the correct drive-removal reference. Use the actual loaded default. Apnea acts in spontaneous breathing, not the separate conscious-breathing branch.

Other actions are unsuitable as the main neural port: conscious respiration prescribes breath commands; mechanical ventilation imposes airway pressure/flow and suppresses spontaneous muscle drive; drug neuromuscular block has additional pharmacological effects. Most `SEOverride` respiratory fields only overwrite reported outputs in `PostProcess`, after the solve. Respiration-rate override also affects respiratory event/output processing; it does not replace the native muscle waveform generation.

**Recommended native extension:** an engine-instance-owned `RespiratoryDrivePort`, applied once in `Respiratory::Process` immediately before `m_Calculator.Process`. Inputs carry an explicit mode and timestamp; outputs expose the requested, latched, generated and applied values separately. Apply the port to both `m_DriverPressure_cmH2O` and `m_DriverPressurePath->GetNextPressureSource()` so internal telemetry/state agrees with the solved source. Validate mode compatibility with conscious breathing, apnea and ventilation; do not inadvertently restore spontaneous drive under a ventilator. Writing the path before `AdvanceModelTime` fails because `RespiratoryDriver()` overwrites it.

Use mutually exclusive ownership modes:

1. `native`: exact original waveform, including native chemical control; the disabled port must preserve trajectory parity.
2. `native_actuation`: preserve native chemoreception/rhythm and pass its active waveform through an explicitly modeled transmission/actuator or test block. A unity passthrough is a numerical port check; gains and delays require evidence and must not be presented as calibrated phrenic physiology.
3. `external_neural`: an independently justified brainstem/phrenic/actuator model supplies active muscle pressure about the native default. Replace the native respiratory waveform as the applied source. Prevent native chemoreceptor driver publication from acting as a second controller; retain its sensing/autonomic functions or mark its respiratory output as diagnostic only. Never add both pressure waveforms.

For an actual phrenic/diaphragm implementation, the model must own rhythm and chemosensory dynamics, conduction state, activation state and a pressure/force conversion with provenance. Choose whether native chemical feedback supplies that model's demand or is replaced; do not have both estimate corrective ventilation. The BioGears circuit has one `RespiratoryMuscle` node feeding both pleural sides (`BioGears.cpp:4722,4780–4791`): selective left/right phrenic deficits and separate intercostal drive require additional actuator/circuit mechanics. Removing diaphragm drive alone cannot be asserted equivalent to removing all respiratory-muscle drive. Preserve IBM donor and pinned import unchanged until a separately versioned model supplies those missing mechanisms.

## Port signals, state and units

| Signals/state | Owner and contract |
| --- | --- |
| PaCO2, PaO2 (mmHg), arterial pH (dimensionless), saturation [0,1], native metabolic demand (W) | Native read-only feedback with native step timestamp; never substitute saturation for PaO2/PaCO2. |
| Peripheral afferent firing (Hz), central/peripheral pressure deltas (cmH2O), frequency deltas (/min), baselines | Native nervous state in native modes. Already serialized by `src/io/biogears/BioGearsPhysiology.cpp:489–580`; add read-only diagnostic exposure, not duplicate Python integration. |
| Requested and cycle-latched amplitude (cmH2O), requested and latched frequency (/min), phase/time (s), generated and applied source (cmH2O), default source (cmH2O) | Native respiratory state. Label pressure amplitude versus total waveform distinctly; negative active pressure drives inspiration. |
| Neural firing (Hz), conduction delay/event queues (s), activation [0,1], muscle force (N), actuator pressure (Pa or cmH2O) | New neural/actuator owner only when defined. Never relabel native pressure as measured phrenic firing. Use an explicit conversion contract: 1 cmH2O = 98.0665 Pa; 1 L = 10⁻³ m³. |
| Flow (L/s), gas volumes (L or mL), pleural/alveolar pressure (cmH2O), O2/CO2 transfer (mL/s), respiratory rate (/min), tidal volume (mL) | Native circuit, gas and blood systems; body geometry follows solved state. Full 3D mechanical reaction remains a separate unresolved constitutive coupling. |
| Port mode, action state, event cursor, actuator filters, delay queues and timestamp | One live runtime checkpoint with native state and external state; callbacks or process globals alone are not serializable state. |

The current immutable telemetry replay cannot run this closed loop. Introduce persistent stepwise native execution with aligned native timesteps and coordinated rollback/checkpoints. Do not run the complete native trajectory first and then call replayed neural outputs feedback.

### Exact first streaming export surface

For the next implementation, keep **native chemoreflex as the sole respiratory controller**. Read completed native telemetry after each tick, deliver scheduled actions before the next tick, and let body geometry consume the resulting native volumes. The optional external-controller mode above is a later extension, not a requirement for this first streaming port. No current cortical IBM output should be mapped to medullary drive.

Public typed readers, using `const SERespiratorySystem* r = bg->GetRespiratorySystem()`, `const SEBloodChemistrySystem* b = bg->GetBloodChemistrySystem()` and `const SENervousSystem* n = bg->GetNervousSystem()`. Every system getter listed below also passed a separate C++20 syntax-check probe against the retained headers:

| Export | Exact getter |
| --- | --- |
| Chemical controller requested amplitude/frequency | `r->GetRespirationDriverPressure(PressureUnit::cmH2O)`, `r->GetRespirationDriverFrequency(FrequencyUnit::Per_min)` |
| Native cycle fraction, not directly measured neural phase | `r->GetRespirationCyclePercentComplete()` (dimensionless [0,1]) |
| Applied driver source and mouth/trachea flow | Existing `EnvironmentToRespiratoryMuscle` path `GetPressureSource(PressureUnit::cmH2O)` and `MouthToTrachea` path `GetFlow(VolumePerTimeUnit::L_Per_s)` after the tick |
| Aggregate flow and pressures | `r->GetInspiratoryFlow(VolumePerTimeUnit::L_Per_s)`, `GetExpiratoryFlow(...)`, `GetMeanPleuralPressure(PressureUnit::cmH2O)`, `GetTranspulmonaryPressure(PressureUnit::cmH2O)` |
| Blood-gas feedback | `b->GetArterialCarbonDioxidePressure(PressureUnit::mmHg)`, `GetArterialOxygenPressure(...)`, `GetArterialBloodPH()`, `GetOxygenSaturation()` |
| Native autonomic effector factors | `n->GetHeartRateScale()`, `GetHeartElastanceScale()`, `GetComplianceScale()`, `GetResistanceScaleExtrasplanchnic()`, `GetResistanceScaleMuscle()`, `GetResistanceScaleMyocardium()`, `GetResistanceScaleSplanchnic()`; all dimensionless factors, not nerve firing or probabilities |

Check each `Has...()` before serializing; represent missing values explicitly. Public `GetRespirationMusclePressure(cmH2O)` is available but `Respiratory.cpp:1805` subtracts a fixed 1033.23 cmH2O from the node pressure, so the actual circuit pressure-source getter is the reliable actuator export, especially outside reference ambient pressure. Inspiratory and expiratory flow outputs are signed opposites of the same tracheal flow; they are not separate nonnegative flow magnitudes.

Live internal nerve diagnostics need an explicit read-only source extension: `Nervous.h` protects `m_AfferentChemoreceptor_Hz`, `m_AfferentPulmonaryStretchReceptor_Hz`, `m_AfferentBaroreceptorAortic_Hz`, `m_AfferentBaroreceptorCarotid_Hz`, `m_AfferentCardiopulmonary_Hz`, `m_SympatheticSinoatrialSignal_Hz`, `m_SympatheticPeripheralSignal_Hz`, `m_VagalSignal_Hz`, central/peripheral frequency and pressure deltas, and the baroreceptor operating point (mmHg). Provide an immutable diagnostic snapshot getter if these signals are needed; none has a public `SENervousSystem` scalar getter. They are aggregate native model signals, not phrenic output. Serialized XML can audit chemical internal state but omits several instantaneous signals and should not be used as a high-rate control API. Native effector factors already influence the cardiovascular model; exporting them must not apply them again through the body bridge.

## Validation experiment to implement

Start matched branches from the same stabilized native state and selected library variant. Record full-rate circuit/driver telemetry (the adapter supports 50 Hz; verify actual native `GetTimeStep`) plus slower blood/chemical state. Proposed durations below are experiment settings, not established physiological guarantees.

1. Baseline 60 s; matched control continues untouched. Verify disabled-port equality and unity-actuation parity against native output at native precision.
2. Remove all active respiratory drive for 30–60 s, retaining passive default pressure and open airway. After passive recoil/transients, pressure excursion and ventilation should collapse; native PaCO2 should rise relative to the matched control and chemical demand should respond. O2 and pH changes must be measured rather than hard-coded. Nonzero residual passive flow is not automatically a failure.
3. Restore drive and run 300–600 s to observe native response and recovery. Verify applied pressure and flow resume, and PaCO2 trends toward control without resetting chemical state. Use actual circuit flow/volume, not only breath-derived rate/volume outputs, which can lag or go stale during apnea.
4. Independently perturb inspired CO2 with a native environment gas-composition action, balancing fractions to one, then wash out. Run for several central-controller time constants where practical. Compare feedback-enabled response with a held pre-perturbation drive-demand branch to establish that blood gas changes alter control, rather than merely showing breathing accompanies metabolism.
5. Save/reload during baseline, suppression, and immediately before action removal; compare against uninterrupted branches. `m_hadApnea` is reset in `Respiratory::SetUp` and is absent from the inspected serializer, so checkpoint equivalence at removal is specifically unresolved and must be tested. Check controller phase, action activity, event cursor and actuator history, not just final gas values.

Success establishes causal native actuation and feedback under the tested numerical model. It does not validate anatomical brainstem dynamics, phrenic laterality, diaphragm recruitment, patient-specific control, or conservation of muscle metabolic work through the separately replayed 3D thorax. Preserve those limitations in the implementation status and UI.
