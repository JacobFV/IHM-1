# Whole-engine respiratory observer acceptance

This isolated probe composes the [native work observer](NATIVE_RESPIRATORY_WORK_PORT.md) with the held BioGears advance method. It does not modify a production adapter, any shared library, chest compliance, native controller, muscle recruitment, or thermal/metabolic source.

`verify_native_respiratory_engine.py` defaults to source preparation. It requires the pinned GI variant and held engine source; explicit `--run` additionally checks the entire variant object inventory and compiles an isolated probe. The generated method preserves original readiness/event guards, PreProcess, Process, PostProcess, events, clocks and tracking byte-for-byte after removing the observation hooks. Duplicated lifecycle anchors are rejected by source-only tests.

The probe takes a pre-PreProcess native lung-volume and requested-amplitude/frequency snapshot. After PreProcess it captures the actual generated pressure waveform, records the later requested amplitude/frequency and cycle phase, optionally adds a held 10 Pa external pressure, and calls `WorkObserver.begin`. After Process it calls `solved` using native next circuit values and the respiratory system's volume output. Original PostProcess, events, clocks and return follow. `commit` occurs only after successful native return and a completed zero signed-muscle record. All modes use that same zero record. No respiratory source work is inserted into the signed-muscle boundary; a failure poisons the probe and work observer terminally.

Four independent loads use the same frozen native human state (`cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3`) and eight 20 ms intervals:

- Baseline executes the original native advance method.
- Observer executes the generated method with a zero external boundary.
- Pulse control executes the generated method without work observation, applying +10 Pa in intervals 3–4 and zero thereafter.
- Pulse observer uses the same pulse plus work observation.

The verifier requires exact equality of all retained native body ports and serialized state content for baseline/observer and pulse-control/pulse-observer. Raw byte equality is also recorded. Native serialization iterates an unordered set of active-substance pointers (`SESubstanceManager.h:57`, `io/biogears/BioGears.cpp:462`), so the content comparator sorts only root `ActiveSubstance` siblings by their unique Name; all tags, attribute values and non-formatting text remain exact. No numerical tolerance or omitted state field is permitted. It checks source/chest flow and constitutive residuals, generated/external work partition, clocks, cumulative committed work, actual airflow and lung-volume perturbation, and zero external work after pressure release. Source-versus-lung stroke discrepancy is retained as evidence, not forced to zero. Serialized-state equality is an observation-parity check; it does not prove exact restore/rollback, and the probe performs no retry or restore.

This short acceptance does not establish complete-breath energetic balance, long-run stability, a mechanical diaphragm/thorax, IBM phrenic recruitment, respiratory muscle chemistry or distributed chest contact. The existing native pressure source still supplies inspiration and native chest/lung elements still own recoil. Production integration must retain these ownership labels and cannot count observed ideal-source work as additional muscle heat or chemical consumption.


## Retained result: full-state acceptance failed

`data/derived/audits/native-respiratory-engine-1s_h5ouy/verification.json` records `passed=false` and `respiratory_observation_checks_passed=true`. All four loads and 32 total native intervals completed and were reaped. Every retained native body-port value matches exactly between its observed/control pair. Source KCL residual is at most `3.401006689815467e-17 m³`, chest constitutive residual `1.525506338004995e-17 m³`, generated/external work-partition residual `1.951563910473908e-18 J`, and source-versus-lung stroke difference `2.029152128442402e-17 m³`. The +10 Pa pulse changed lung volume by `−2.7040860894867365 mL` versus its zero-load control at interval 4 and airway flow by `−0.07262684774126504 L/s` at interval 3. After release the external-pressure work component is exactly zero; the native source and gas circuit continue advancing.

Raw saved-state bytes differ in active-substance ordering. After sorting only that source-defined unordered collection, the sole state-content mismatch is `System[37]/ReactionTime[17]`:

| Pair | Control value (s) | Observer value (s) |
|---|---:|---:|
| Baseline / observer | `0` | `9.65804430779721e-313` |
| Pulse control / observer | `9.46096256785514e-313` | `9.54707288710903e-313` |

The source exposes a pre-existing continuation defect: `Nervous::Initialize()` assigns `m_ReactionTime_s` and `m_TiredTime_hr`, but the constructor/Invalidate/SetUp path does not, and `BiogearsPhysiology::UnMarshall(BioGearsNervousSystemData)` (`491–547`) restores neither private variable. `CalculateSleepEffects()` (`1313–1331`) assigns reaction time only in conditional branches and otherwise reads its uninitialized private value, writing it to the public scalar. The private tired-time accumulator is likewise read conditionally without continuation restoration. The public reaction-time scalar is serialized, but that does not restore its private writer. This source defect prevents an exact whole-state observer-parity claim. The probe does not initialize those fields by guess, filter the differing value, round it away, or change any native library.

The initial failed strict byte comparison and all native binary/log/state artifacts remain retained. `--analyze-retained <directory>` reruns only source binding and offline receipt checks; it performs no compilation or native advance and intentionally returns failure while state content differs. Source preparation, generated-method restoration and duplicate-anchor checks pass. Production adapter integration remains pending a separately reviewed, source-aware native continuation correction and renewed parity acceptance. No respiratory mechanical actuator or new recoil model is enabled by these results.

Resources remained bounded: compile 1.95 s / 543476 KiB peak RSS; each native load/run 0.27–0.29 s / at most 80272 KiB RSS. The offline XML comparator additionally tests that unordered substance reordering is permitted, a `0` → `1e-313` scalar change fails, duplicate substance names fail, ordinary child reordering fails, and non-whitespace tail text cannot disappear from comparison. Comparison and missing-continuation-field source bytes are checked against the held Git revision.
