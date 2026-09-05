# IHM-1 — implicit human body model

A body-wide probabilistic substrate, modeled after `../IBM-1`. Explicit predictors
are dependency-traced views of shared physical state, evidence and processes.

The primary 3D view is now **one canonical generic body** with 2,408 anatomical
entities, including registered lung lobes, thyroid and lymphatic structures,
explicit skin-layer priors, source-informed mechanics and an 80-region IBM-derived
brain. Source viewers remain available as evidence inspection. Acquired,
registered and directly synthesized structures retain separate provenance and
uncertainty; an unmeasured property never acquires a confidence percentage merely
by being rendered.

The [canonical body runtime](docs/CANONICAL_BODY_RUNTIME.md) materializes actual
native physiological samples into synchronized tissue and neural states. The
30-second generic baseline has 1,500 native samples and 301 display frames, with
finite Laplace spectra computed from the native clock. Source hashes, synthesized
transfer assumptions and numerical approximation remain attached. **This is an
executable research assembly, not a complete or independently validated digital
twin.** Mechanical and neural feedback into native physiology, articulated motion,
whole-body voltage coupling and calibrated tissue-resolved drainage remain open.

The anatomy viewer preserves full source surface topology: 6,681,030 triangles
from BodyParts3D and 8,111,194 from Z-Anatomy, with source frames and overlapping
ancestry retained separately. [Anatomical grouping](docs/ANATOMY_COVERAGE.md) now
uses recorded ontology paths rather than name substrings.

The acquired evidence now materializes into an executable native multisystem
engine, measured population beliefs, human wound-field predictors, source-derived
skin–lymph circuits and temporal predictors. A local 3D workbench displays the
actual anatomy, wrapped muscle paths, vascular flow states and physiological
trajectories. **Independent whole-human calibration is not established.**

```bash
# Dependencies and acquired data are already installed in this workspace.
.venv/bin/python -m ihm integrated --output artifacts/integrated-human.json
.venv/bin/python -m ihm serve --port 8765
# Open http://127.0.0.1:8765
```

For a fresh checkout, see [workbench setup](docs/APP.md) and the source-specific
collection/build scripts. Build the frontend with `npm ci && npm run build`
inside `app/`. Raw data and native binaries are intentionally not in Git.

The [implicit human API](docs/IMPLICIT_HUMAN.md) selects a predictor appropriate
to each evidence family. [Temporal models](docs/TEMPORAL_MODELS.md) distinguish
finite Laplace transforms of observations from causal circuit transfer functions;
[the native circuit predictor](docs/CIRCUIT_PREDICTOR.md) exposes actual fast and
slow vascular/interstitial/lymphatic modes. Published CSF, reproductive endocrine, and supine thermal models also execute
through the same API. The original illustrative scaffold remains available separately.

## Acquired real data

| Collection | Local contents | Details |
|---|---|---|
| BodyParts3D | 2,234 OBJ meshes and anatomical ontology tables | [Anatomy](docs/ANATOMY_DATA.md) |
| Z-Anatomy | 2,581 full-resolution evaluated surfaces; 495 muscles, 158 lymph-node groups, 256 body-surface regions | [Extended anatomy](docs/EXTENDED_ANATOMY.md) |
| Published lymphatic graph | 996 vertices, 1,117 structural edges and source lengths; model-derived, not measured | [Lymphatic network](docs/LYMPH_NETWORK.md) |
| OpenSim | 34 model variants, 343 meshes, 4,480 muscle path/attachment points across variants | [Anatomy](docs/ANATOMY_DATA.md) |
| Vascular Model Repository | Three human meshes, nine flow waveforms, BCs, 201 cerebral pressure/velocity time states | [Vascular](docs/VASCULAR_DATA.md) |
| Integrated physiology | HumMod, Physiomodel, Physiolibrary, BioGears and BETSE sources; 14,672 parameter/equation candidates and 475 literature targets | [Physiology](docs/PHYSIOLOGY_DATA.md) |
| NHANES 2017–2018 | 129 data tables + 129 codebooks, 9,254 participants, 4,292 variable definitions | [Population](docs/POPULATION_DATA.md) |
| Human Reference Atlas | 1,759,409 structure/cell/biomarker rows; 1,082 vessel records; 193 geometric summaries | [Semantics](docs/SEMANTIC_DATA.md) |
| Supine thermoregulation | Pinned JOS-3, 85 thermal nodes, 17 regions, four one-hour cases and measured bedding resistance | [Thermal model](docs/THERMAL_MODEL.md) |
| Human skin wound fields | 40-person experiment, 155 numeric table cells, 16 group means/SEMs | [Skin measurements](docs/INTEGUMENTARY_DATA.md) |

The local raw corpus is approximately 6.00 GB including retained archives and
extracted assets. The searchable catalog contains 36,656 parameter/constraint
records, **not 36,656 experimentally calibrated coefficients**. All source models
and overlapping variants retain their identities. BETSE's actual full solver has
also run locally on 212 cells; this was a generic-tissue execution check, not human
wound validation.

The original **BioGears integrated engine now also runs natively on ARM64**:
60 seconds after stabilization, 3,000 samples at 50 Hz, and 22 requested
physiological outputs. Its original cardiovascular, respiratory, renal, endocrine,
energy and blood-chemistry mechanisms remain together. See the
[native backend and reproduction steps](docs/NATIVE_BACKEND.md) and
[actual multisystem traces](artifacts/native-multisystem-physiology.png).
This is a resting upstream patient; native supine posture and independent calibration are not established. The
3D viewer distinguishes a rigid display pose from the separate physical
hydrostatic adjunct; the canonical body adds explicit registered geometry and reduced one-way tissue drivers, without claiming a volumetric whole-body contact/flow solution.

The [supine thermal model](docs/THERMAL_MODEL.md) uses the authors' JOS-3
implementation and separately identified bedding boundaries. It retains 85 node
temperatures, explicit heat balances and timestep refinement. This is a separate
source subject. BioGears hour-long rest runs expose upstream cooling; the
[native audit](docs/NATIVE_BACKEND.md) preserves that behavior and the reviewed
bounds, evaporation telemetry and experimental dimensional corrections.

To rebuild the views from the acquired corpus and verify the numerical models,
API, production frontend and an actual native scenario through the browser:

```bash
.venv/bin/python scripts/build_workbench.py
RUN_NATIVE_BROWSER=1 .venv/bin/python scripts/verify_all.py --native --app
```

The verification suite starts its own local server. Its command outputs, timings
and return codes are retained in `artifacts/verification/report.json`.

```bash
uv pip install --python .venv/bin/python -e '.[data,plot]'
.venv/bin/ihm data summary
.venv/bin/ihm data search --term insulin --limit 10
.venv/bin/ihm data search --term e50_W --source biogears
.venv/bin/ihm build whole-body --subject adult-population-reference \
  --population-prior data/derived/population/nhanes-2017-2018/joint-population-prior.json \
  --output artifacts/empirical-body-prior.npz
```

The last command initializes 23 physical state components using survey-weighted
measured covariance. It does not alter or validate the illustrative dynamic
coefficients. The nonlinear cardiopulmonary example does not yet accept this
Gaussian population initialization.

Source/derived assets are local under `data/raw/` and `data/derived/`; they are
excluded from source distributions. Reproduction scripts retain download URLs,
revisions, hashes, codebooks and source-specific licensing/access limitations.

## Original executable modeling scaffold

```bash
uv venv .venv
uv pip install --python .venv/bin/python -e .
.venv/bin/ihm inventory
.venv/bin/ihm models
.venv/bin/ihm demo --output artifacts
```

The package uses Python 3.11+, NumPy and SciPy and has no dependency on IBM-1.

## One ontology; heterogeneous supports and graphs

The four primitives remain **fields, anatomy, topologies, processes**. A field
component is a physical quantity with units, support and region. Components can
coexist at one location. Anatomical partitions can overlap across systems.
Topologies identify permitted interactions; process parameters determine strength
and dynamics. No body-wide adjacency matrix is imposed on unrelated tissues.

The inventory includes cardiovascular/vascular, blood, interstitial, lymphatic,
immune, hematopoietic/splenic, respiratory, renal/urinary, digestive/hepatic/pancreatic,
endocrine, adipose/metabolic, musculoskeletal, nervous/sensory/CSF, integumentary,
and reproductive/uterine/placental state. See [system coverage](docs/SYSTEMS.md).
Reproductive and placental entries describe possible anatomy; materialize only
applicable targets for a subject. The `whole-body` view is an inventory reference,
not a claim that every person has every registered structure.

Four executable representations share the same registry and runtime:

- A nonlinear, supine cardiopulmonary loop that generates actual beats and tidal breaths.
- Organ/compartment Gaussian state with illustrative affine coupling.
- Skin patches with individually parameterized cells, ion reservoirs, membrane
  capacitances, ion conductances, gap junctions, epithelial barriers and a separate
  extracellular conduction graph.
- Conservative vascular–interstitial–lymphatic fluid networks, with distinct
  filtration, uptake, collecting, lymph-node and lymphovenous interfaces.

## Does it breathe and does its heart beat?

The `cardiopulmonary` view does: four contracting chambers, pressure-operated
valves, a closed pulmonary/systemic blood loop, and respiratory effort driving
tidal airflow. Pleural pressure carries one-way coupling from breathing into circulation. The baseline is
lying supine at rest; movement is deferred.

```bash
.venv/bin/ihm build cardiopulmonary --subject supine-example \
  --time 30 --output artifacts/cardiopulmonary.npz
.venv/bin/python scripts/plot_cardiopulmonary.py
```

It uses prescribed pacing and respiratory rhythms, not full electrical conduction,
reflex regulation or gas chemistry. This nonlinear view is deterministic and has
no measurement assimilation yet. [Mechanisms and limits](docs/CARDIOPULMONARY.md).

## Integumentary bioelectricity

Membrane voltage, transepithelial voltage and extracellular electric fields are
different quantities. `skin` keeps the membrane and extracellular networks
separate, with explicit tissue geometry and cell phenotypes. The supplied example
is a synthetic one-dimensional cross-section; custom 3D node/edge geometries are
accepted without inventing a universal nearest-neighbor graph.

```bash
.venv/bin/ihm build skin --subject synthetic-skin \
  --source data/examples/synthetic-skin-voltage.card.json data/examples/skin.json \
  --time 5 --output artifacts/skin.npz
.venv/bin/ihm forecast artifacts/skin.npz --times 5 10 30
```

```python
from ihm.materialize.skin import SkinPatch, skin_model, electric_field

patch = SkinPatch.line(n=9, wound=True)
model = skin_model(patch, subject="example", view="combined")
model.advance(5)
field = electric_field(model, patch)  # V/m and propagated standard deviations
```

`view="surface"` and `view="membrane"` materialize different dependencies.
`SkinPatch(cells, surface, gap_edges, surface_edges, ...)` accepts heterogeneous
cell parameters and independent contact graphs. `--geometry patch.json` accepts
the JSON representation of `dataclasses.asdict(patch)`.

The membrane mechanism uses Ohmic ion currents and a Nernst linearization around
declared concentrations. The extracellular sheet uses current sources and
barrier conductances. A wound creates a low-resistance shunt. Basal extracellular
potential and intracellular ion reservoirs are prescribed; this is not a full
bidomain or electrodiffusion model. There is no implemented cell migration,
proliferation, pattern-memory or wound-closure predictor. See [bioelectric evidence
and limitations](docs/BIOELECTRICITY.md).

## Vascular state, blood state and lymph drainage

Vascular mechanics are distinct from blood composition and interstitial fluid.
The specialized fluid materialization exchanges volume across explicitly named
interfaces and conserves total volume in a closed network.

```bash
.venv/bin/ihm build fluids --subject example --time 10 --output artifacts/fluids.npz
.venv/bin/ihm build blood --subject example --output artifacts/blood.npz
.venv/bin/ihm build lymphatic --subject example --output artifacts/lymphatic.npz
```

`FluidNetwork` allows arbitrary compartment sizes, compliances, pressure offsets,
conductances and pump heads. Its valve regime is fixed; it reports reverse flow
rather than clipping it. It does not yet transport solutes or implement active
lymphatic valve gating. The organ-scale `lymphatic` view remains illustrative;
use `fluids` to exercise conservative drainage physics.

## Heterogeneous evidence

CSV, JSON and NPZ loaders consume explicit measurement records:

```
id,subject,time,component,value,unit,variance
lab-0,synthetic-person,0,metabolic.glucose,108,mg/dL,81
```

Time is seconds relative to a shared subject-specific origin; variance is in
squared source units. Source cards declare permitted bindings, evidence kind,
species and reference. Unit conversion scales variance as well as values.
Missing values supply no evidence. Unknown bindings/units, nonfinite data,
duplicate records, nonhuman transfer and mixed subjects fail explicitly.

```bash
.venv/bin/ihm build glucose --subject synthetic-person \
  --source data/examples/synthetic-lab.card.json data/examples/labs.csv \
  --source data/examples/synthetic-cgm.card.json data/examples/cgm.json \
  --output artifacts/glucose.npz
.venv/bin/ihm forecast artifacts/glucose.npz --times 600 900 1200
.venv/bin/ihm inspect artifacts/glucose.npz
```

The CLI assumes independent measurement errors within and between batches.
For simultaneous correlated errors, use
`model.assimilate(evidence, error_covariance=R)`. Cross-batch error dependence is
not modeled. Source references, input SHA-256 hashes, evidence identities, units,
selected dynamics and interventions survive export. Raw images, omics and
uncalibrated voltage-sensitive dyes need measurement operators; they are not
silently converted into direct physical-state observations. The active acquisition manifests under `data/derived` identify what has actually
been downloaded and indexed; remaining catalog-only sources stay explicitly marked.

## Fitting and materialization

```python
from ihm import body, materialize, Request
from ihm.forge.fit import fit_affine, evaluate

model = materialize(body(), Request(("blood.pressure", "renal.filtration"), "person-1"))
# X contains physical process inputs; dx_dt is that process's pressure on output.
# fit = fit_affine(X, dx_dt, input_ids, output_id, source_id, training_subjects)
# process = fit.process(existing_process_id, topology_id, noise=diffusion_intensity)
# fitted_model = materialize(body(overrides=(process,)), request)
```

Fitted forms preserve a process's input/output/topology contract. Multiple
processes may write one output, so fitting its total derivative to an individual
process requires subtracting the other contributions or isolating the process
experimentally. Held-out evaluation checks training-subject overlap; source
labels and data partitioning remain the caller's responsibility.

Forecasting is exact for the selected affine stochastic model. It propagates
state/process-noise covariance, not parameter-posterior uncertainty. Gaussian
support is unbounded, so biologically impossible states can occur outside the
local operating regime. Dense covariance limits the current runtime to small
materializations; cell-scale whole-body simulation is not implemented.

Artifacts are portable compressed NPZ containing JSON metadata and numeric
arrays, loaded without pickle. Forecasts do not mutate the saved posterior.
State clamps are hard interventions over the requested interval.

## Verification

```bash
.venv/bin/python scripts/verify.py
.venv/bin/python scripts/verify_evidence.py
.venv/bin/python scripts/verify_skin.py
.venv/bin/python scripts/verify_fluids.py
.venv/bin/python scripts/verify_cardiopulmonary.py
.venv/bin/python scripts/verify_population.py
.venv/bin/python scripts/verify_collected_data.py
```

Checks cover analytic Gaussian conditioning, correlated observations, source
provenance, subject isolation, units, duplicate rejection, dependency closure,
budgets, covariance scale, fitting, export, skin RC relaxation, assembled current
conservation, wound-field symmetry and closed-network fluid conservation.
These establish software/numerical behavior, not biological predictive validity.

Optional plot: install `matplotlib` and run `scripts/plot_skin.py` to export the
intact-versus-wound comparison. [Architecture](docs/ARCHITECTURE.md),
[evidence contract](data/README.md), [design](docs/superpowers/specs/2026-09-04-human-body-design.md).
