# Measured population physiology

The CDC/NCHS 2017–2018 NHANES public release was collected from the official
[category listings](https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory&Cycle=2017-2018):
129 SAS XPORT tables and 129 codebooks across demographics, examination,
laboratory, dietary and questionnaire components. Original files and codebooks
are retained with HTTP source URLs and SHA-256 receipts. This is real public,
deidentified participant data, not the synthetic examples from the initial build.

`index_nhanes.py` converts all tables to Parquet without pretending they are all
one-row-per-person: dietary and other repeated records retain their multiplicity.
There are 9,254 unique demographic participant IDs and 4,292 table-specific
variable definitions; repeated-table row counts are not people counts.

## Population priors that can be materialized

Twenty-three physical measurements are joined by exact participant identity with
age, sex, survey strata/PSUs and survey weights. The adult subset is age 20+ and
includes illness and treatment; it is not labeled a healthy reference population.

- Hemoglobin/hematocrit/platelets and neutrophil/monocyte counts.
- Blood albumin, urea, creatinine, bilirubin and chemistry/electrolytes.
- Fasting glucose/insulin, CRP, peripheral pulse and urine collection flow.

The released canonical/SI columns are preferred; explicit conversions are applied
where needed. Plasma/serum chemistry is not silently mapped to interstitial
measurements, nor voided urine to nephron output. Fasting analytes and joint
complete-case estimates use positive `WTSAF2YR`; other marginal estimates use
positive `WTMEC2YR`. Sampling design columns remain in the linked output.

The stable participant hash split reserves an internal holdout. Joint covariance
is computed only among complete cases with positive fasting weights; selection
bias remains possible. Marginal means/dispersion are also retained by sex and
combined adults. Dispersion is a population distribution, not an assay error or
a survey-design standard error. Detection-limit flags are preserved; CDC-provided
substitutions are retained and reported rather than treated as uncensored truths.

```bash
.venv/bin/python scripts/collect_nhanes.py
.venv/bin/python scripts/index_nhanes.py
.venv/bin/ihm build whole-body --subject adult-population-reference \
  --population-prior data/derived/population/nhanes-2017-2018/joint-population-prior.json \
  --output artifacts/empirical-body-prior.npz
```

Only initial state means/covariance and their provenance change. Dynamic process
coefficients and forward-validation status do not change. This is a source-backed
population prior, not a calibrated digital human.

## Decoder verification

Review reproduced an installed pandas XPORT decoder defect: an all-zero SAS
floating-point byte sequence was returned as `2**-260`, including zero survey
weights. The pipeline now uses the independent `pyreadstat` decoder. Codebook
checks confirm 325 zero fasting weights and 550 zero MEC weights, so nominal
zero-weight participants are excluded. Original XPORT bytes remain unchanged.
