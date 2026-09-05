# Evidence contract

`examples/` contains synthetic records only. Actual public source data now live
under `raw/`, with parsed and indexed artifacts under `derived/`. Acquisition
manifests retain URLs/revisions, hashes and source-specific terms. The consolidated
catalog is `derived/evidence-catalog.sqlite`; `derived/collection-summary.json`
reports actual coverage. Some legacy `sources/` entries remain catalog-only.

See [source-grounded design](../docs/DIGITAL_HUMAN_DESIGN.md),
[anatomy](../docs/ANATOMY_DATA.md), [vascular](../docs/VASCULAR_DATA.md),
[physiology](../docs/PHYSIOLOGY_DATA.md), [population](../docs/POPULATION_DATA.md)
and [semantics](../docs/SEMANTIC_DATA.md). Public data and upstream model parameters
are not automatically a calibrated individual or jointly validated whole body.

A local `SourceCard` specifies `id`, `kind`, `components`, `reference`, `species`
and `status`. CSV and JSON use one row per measurement. NPZ uses equal-length
one-dimensional arrays with these exact column names:

| Column | Meaning |
|---|---|
| id | stable record identity within source |
| subject | exact subject identifier; already reconciled upstream |
| time | seconds relative to an explicitly shared subject time origin |
| component | registered physical component identifier |
| value | numeric measured/derived value; blank/null means no evidence |
| unit | exact unit token with an explicitly supported conversion |
| variance | strictly positive error variance in squared source units |

`load(path, card, registry)` returns validated `Evidence` objects and hashes input
bytes. It rejects unknown bindings, units, nonfinite values, invalid variance,
duplicate rows and nonhuman/catalog-only conditioning. Hashes establish identity,
not truthfulness or permission. Source references and hashes persist in artifacts.

Waveform samples may be supplied in NPZ, but arbitrary waveform amplitude is not
a measurement variance. No implicit resampling, synchronization or averaging is
performed. All clocks and identifiers must be reconciled before loading. A stable
record identity must survive format changes to prevent double counting.

Correlated same-time evidence must carry its full error covariance through the
Python API. The CLI assumes independent errors. Separate times/batches are assumed
independent conditional on state. Subject cohorts must be explicitly split; a
population average is not an observation of a particular person.

Images constrain anatomy only after segmentation/registration with a declared
frame. Omics constrain identified physical quantities or process parameters
through a justified operator. Fluorescence needs dye calibration before it is
voltage evidence. Species transfer needs an explicit model and uncertainty.
None of these operations is silently guessed by the generic loaders.
