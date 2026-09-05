# Acquired integrated physiology source and coefficient evidence

This collection contains actual upstream equations, configuration constants, benchmark assets, and literature reference targets. It does **not** make all extracted constants experimentally calibrated, nor attach these engines to the IHM runtime.

Rebuild the index with `.venv/bin/python scripts/collect_physiology.py`. Add `--download` to clone missing repositories. Existing checkouts are not updated automatically. Exact revisions, licenses, byte counts and collection time are in `data/derived/physiology/collection.json`; every tracked file has its SHA-256 in `assets.jsonl`. Source downloads occupy about 657 MB with Git metadata (527,932,492 tracked source bytes).

| Source | Local checkout | Material actually acquired |
|---|---|---|
| [BioGears](https://github.com/BioGearsEngine/core) | `data/raw/physiology/biogears` | C++ integrated engine; patient/substance/environment XML; cardiovascular, respiratory, renal, endocrine, energy, nervous, GI, hepatic and tissue methodology; bibliography; validation reference CSVs, prior engine output CSVs, validation XLSX and scenarios |
| [HumMod standalone](https://github.com/HumMod/hummod-standalone) | `data/raw/physiology/hummod` | DES/XML equations, response curves, subsystem references, benchmarks and distributed standalone binaries |
| [Physiomodel](https://github.com/physiology/Physiomodel) | `data/raw/physiology/physiomodel` | Whole-body Modelica model derived from HumMod 1.6.1, annotated types/equations, default/steady-state assets and literature references |
| [Physiolibrary](https://github.com/MarekMatejak/Physiolibrary) | `data/raw/physiology/physiolibrary` | Modelica physiological component and unit definitions; **current** checkout, not the historical dependency revision required by Physiomodel |
| [BETSE](https://github.com/betsee/betse) | `data/raw/physiology/betse` | Electrodiffusion engine, ion channels/pumps, biochemical networks, example configurations and published 2016 Frontiers attractor example configs |

Public-use license routes: BioGears Apache-2.0 (third-party components separate); BETSE BSD-2-Clause; current Physiolibrary BSD-3-Clause; Physiomodel public GPL-3.0 route under its Physiomodel License 1.0. HumMod XML descriptions are GPL-2.0, while its executables and other material have distinct restrictive terms in its README. The presence of binaries in the upstream archive is not permission to embed, modify or commercialize them. No access agreement or account was used; no HumMod executable was run.

## Derived records and evidence strength

`coefficients.jsonl` contains 14,672 records: HumMod 7,842; Physiomodel 251; BioGears 3,370; BETSE 3,209. This is a discovery index, **not 14,672 independent fitted physiological coefficients**. It includes equation definitions, curve coordinates, model parameters, configuration values and C++ numeric assignment candidates (including initializers/algorithm constants). Physiologically meaningful numbers embedded in C++ expressions, Modelica component modifications and other syntactic forms are not exhaustively extracted.

Each record includes source, relative source file, line, file SHA-256, parameter name, literal numeric value if available, full expression, units if explicit, units status and nearby equation context. Source paths are relative to `data/raw/physiology`. No expression is evaluated by the collector. Missing units stay null. Typed Modelica declarations retain their type and require version-specific unit resolution; no display unit is mistaken for an SI storage unit. HumMod response-curve x/y/slope coordinates remain separate and require their original interpolator and units. Many HumMod constants cannot be interpreted independently of the original equation graph.

`reported_fitted_coefficients.jsonl` is a separate, manually scoped evidence route: three BioGears epinephrine exercise-response logistic parameters (`e50_W=190 W`, `eta=0.035 1/W`, `maxMultiplier=18.75`, yielding an asymptote of 19.75 after the added baseline). The source code states that the shape was adjusted to Tidgren 1991 data and the maximum was adapted from Tidgren 1991 and Stratton 1985. The methodology describes fitting the logistic curve. Eta's units follow dimensional analysis of the source logistic equation and are labeled as an inference. These are **source-authors-reported fits**, not an independently reproduced calibration. Original observation points, objective function, fit residuals and uncertainty were not recovered. The code's basal epinephrine scaling also differs from the approximate value in its methodology, so the archived code revision is authoritative for execution.

`calibration_evidence.jsonl` indexes upstream methodology statements about fitting/calibration/tuning, preserving their source and citation keys. These statements require manual mapping to exact parameters; they do not automatically promote the generic coefficient index.

`validation_targets.jsonl` contains 475 BioGears literature target rows, with value/range text, units, citation keys, bibliography file, notes and system. These are reference values as curated by upstream, not raw participant observations. Forty-four rows have at least one citation key not found in the included bibliography, explicitly flagged. Results CSVs shipped by upstream are retained as simulated outputs, **not measurements**, and are not labeled locally reproduced.

## Cross-system coverage and compatibility gaps

The acquired equations cover cardiac output and pressure/flow, gas exchange and oxygen delivery, renal filtration/reabsorption, water/electrolyte regulation, endocrine control, metabolism/heat and neural feedback. HumMod/Physiomodel supply numerous slower fluid/endocrine couplings absent from simple cardiopulmonary scaffolds. BioGears methodology explicitly discusses its feedback loops and known omitted endocrine interactions; acquisition does not remove those omissions. BETSE contributes cellular electrodiffusion, gap junction and biochemical coupling; its generic tissue model is not a subject-specific human skin model.

These engines are not plug-compatible coefficient libraries. Their compartment volumes, patient baselines, species assumptions, solvers, state definitions, response-curve semantics and time/unit conventions must be preserved. BioGears tunes some circuit properties during stabilization against the configured patient. Flattening those tuned quantities into fixed universal constants loses the original behavior. Physiomodel requires compatible Modelica tooling and its historical Physiolibrary version; no Modelica compiler was available, so it was not compiled. The BioGears C++ engine was acquired but not built or executed in this collection task.

## Actual BETSE execution

An isolated `.venv-physiology` was created and the acquired BETSE revision installed editable. Its complete dependency freeze is in `data/derived/physiology/betse_run/requirements-lock.txt`. Generated default config, geometry, state pickles, logs and outputs are under `betse_run/`. The upstream tracked files remain unchanged.

Commands used sequentially:

```sh
python3 -m venv .venv-physiology
.venv-physiology/bin/pip install -e data/raw/physiology/betse
.venv-physiology/bin/betse --headless config data/derived/physiology/betse_run/sim_config.yaml
.venv-physiology/bin/betse --headless -v seed data/derived/physiology/betse_run/sim_config.yaml
.venv-physiology/bin/betse --headless -v init data/derived/physiology/betse_run/sim_config.yaml
.venv-physiology/bin/betse --headless -v sim data/derived/physiology/betse_run/sim_config.yaml
```

The generated default profile uses the full BETSE solver. Seeding completed, and initialization ran 500 time steps over 5 seconds, finishing at average membrane voltage -44.884 mV, intracellular sodium 11.9376 mmol/L and potassium 139.058 mmol/L. This is execution evidence and numerical smoke verification; there was no comparison to experimental observations and no claim of wound-healing validation. The historical published attractor configuration was retained but not replicated in this run. The main simulation also completed successfully: 350 steps over 0.035 seconds on 212 cells, engine-reported final average cell membrane voltage -43.09 mV (the unweighted mean over 1,221 membrane edges is -43.6258 mV). All stored membrane voltages and intracellular concentration arrays were finite. `mean_membrane_voltage.csv` exports the 34 sampled average voltages; `run_summary.json` records shapes, extrema and final results. Random geometry means those summaries can differ on reseeding.

A first initialization command accidentally started while the seed gzip was still being written and failed with EOF; after seed completion, initialization was rerun successfully. No code or numerical setting was changed to suppress a failure. Final completed-command logs are retained; `run_summary.json` records the final outcome and this execution history.

Integrity verification: `.venv/bin/python scripts/collect_physiology.py --verify` checked all 8,850 tracked file hashes, all 14,672 coefficient provenance rows and source anchors for the BioGears epinephrine coefficient, Physiomodel ventricular EDV and HumMod ADH/kidney/skin coverage. Three empty HumMod template declarations were detected and excluded. Results are saved in `verification.json`. This check verifies extraction integrity, not physiological validity.
