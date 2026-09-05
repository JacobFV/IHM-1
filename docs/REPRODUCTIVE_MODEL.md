# Reproductive endocrine source model

Run `.venv/bin/python scripts/collect_reproductive.py` to acquire and execute the official Physiome Schlosser–Selgrade model, then `.venv/bin/python scripts/verify_reproductive.py`. Raw CellML, official generated Python, code-generation page, license page and provenance are retained in `data/raw/reproductive/schlosser_selgrade_2000/`. Acquisition is bounded to 250 KB per file and no downloaded source file is overwritten on repeat collection. The generated outputs are `data/derived/reproductive/index.json` and `trajectory.json`.

## Source and attribution

The model is [Schlosser and Selgrade, 2000, Physiome exposure](https://models.physiomeproject.org/exposure/38fcbccf8eacd4c93e5aa531625a8950/schlosser_selgrade_2000.cellml), derived from [permanent workspace commit 73464388213e0eeb5ac101b80c211544db027123](https://models.physiomeproject.org/workspace/schlosser_selgrade_2000/rawfile/73464388213e0eeb5ac101b80c211544db027123/schlosser_selgrade_2000.cellml). Its generated-code page explicitly identifies that same commit. Original article: Paul M. Schlosser and James F. Selgrade, *A model of gonadotropin regulation during the menstrual cycle in women: qualitative features*, Environmental Health Perspectives 108:873–881 (2000), PMID 11035997. CellML author: Catherine Lloyd, Auckland Bioengineering Institute, University of Auckland. Repository citation: Lloyd et al., *The CellML Model Repository*, Bioinformatics 24(18):2122–2123 (2008).

The model's [license and citation page](https://models.physiomeproject.org/exposure/38fcbccf8eacd4c93e5aa531625a8950/schlosser_selgrade_2000.cellml/license_citation) specifies Creative Commons Attribution 3.0 Unported. Attribution, exposure URL, workspace revision, acquisition UTC timestamp, license URL and each file's SHA256 are retained. The executable source is unchanged; the local adaptation replaces only integration/output orchestration, not model equations or physiological coefficients.

Reviewed CellML SHA256: `46a4cd42dda1299862ffba0086cbb910710dd9a55b1212a4b97d17341e9cb573`. Reviewed official generated Python SHA256: `910a322c4ac848f4467915a86858d56c3e0bf8da62f6faaf40929d851370c195`. Both are checked before import/execution. The inspected generated module contains math/NumPy imports, equation definitions and a guarded plotting entrypoint; the latter is not executed.

## Executed mechanisms and scope

The source contains four differential states: releasable LH pool, circulating LH concentration, releasable FSH pool and circulating FSH concentration. Synthesis fills pools, regulated release transfers hormone into blood and clearance removes it. Twenty original parameter constants and four initial values are checked against the same CellML variable definitions and units. All MathML numeric literals, including prescribed input-function coefficients, are separately inventoried so those values do not disappear behind the generated code.

Estradiol E2, progesterone P4 and inhibin Ih are prescribed nonperiodic functions of time. Their shifted forms use fixed source parameters dE=0.42 days, dP=2.9 days and dIh=2 days. They are not ovarian state variables responding to the simulated LH/FSH. This is therefore a forced gonadotropin response model, **not an autonomous menstrual-cycle oscillator**, and no emergent 28-day period is claimed or tested. Ovulation, follicle maturation, uterine tissue, pregnancy, individual endocrine calibration and patient coupling are outside its executed state space.

The source documentation includes an older statement that delays were not represented and a later modification note that fixed delays were handled. The actual pinned equations contain shifted analytic input functions. That executed meaning is retained transparently; the code does not invent a delay differential equation or silently resolve all historical documentation differences.

## Units and integration

The source defines one day as 86,400 seconds. Pool states use micrograms and circulating LH/FSH use micrograms/litre. E2 uses ng/L, P4 nmol/L and inhibin uses its source assay unit U/L. Source unit definitions accompany the artifact. Channel metadata provide SI conversion scales where physically specified: microgram is 1e-9 kg, microgram/litre is 1e-6 kg/m³, ng/L is 1e-9 kg/m³ and nmol/L is 1e-6 mol/m³. The source supplies no inhibin assay-to-mass/molar conversion, so its SI unit and scale remain null.

The default run preserves the official code-generation example's 0–10-day window, sampled at 501 times. The wrapper uses SciPy BDF with rtol 1e-9, atol 1e-10 and maximum step 0.1 day. It independently repeats the integration with DOP853. Numerical tolerances and timestep limits are solver settings, not new physiological coefficients. The public wrapper permits at most 30 days and 10,001 output samples; extending the input functions does not make them periodic or establish biological validity beyond the example.

## Results and verification

All four states remain nonnegative. Maximum normalized BDF/DOP853 disagreement is approximately 7.11e-9. Maximum absolute differences per state are below 2.46e-6 in the respective source units. An independent conservation identity checks that internal release cancels between pool mass and circulating concentration times the source distribution volume. Residuals are below 7e-13 micrograms/day, leaving only external synthesis and clearance.

On the saved sample grid, LH peaks near 109.98 micrograms/litre at day 0.76 and FSH near 206.76 micrograms/litre at day 1.04. These are outputs of this particular source initialization/input history, not clinical reference ranges or an independently validated patient fit. Progesterone and inhibin prescribed inputs peak at day seven. The cycle-period target is explicitly null because this source has no autonomous cycle-generation target.

The tests check pinned source identity, initial values and units against CellML, release mass cancellation, positivity, two-solver agreement, conversion of the ten-day time horizon to 864,000 seconds, meaningful SI metadata and bounded invalid-input rejection. No measured holdout is available or fabricated.

## API and artifact schema

`load_source(root)` returns the pinned official equation module and source audit. `run_reproductive(root, duration_days=10, samples=501, method='BDF')` executes the source and returns a JSON-compatible result. `build_reproductive(root)` runs both solvers, checks agreement and writes artifacts.

`trajectory.json` contains `id`, `system`, `source_kind='source_model_simulation'`, source/license provenance, audit, `time_days`, `time_s`, state-major `state_values`, and 16 `channels`. Each channel includes source id/component/unit, kind (`state`, `algebraic` or `prescribed_time_input`), values, source-unit extrema, peak time and SI conversion metadata. Parameters, solver settings, validation and limitations remain alongside the trace. The compact index uses `models` with source metadata, parameter evidence, channel summaries and `trajectory_path`. This separate source family must not be labeled as measured human evidence, illustrative invented coefficients or a coupled BioGears patient.
