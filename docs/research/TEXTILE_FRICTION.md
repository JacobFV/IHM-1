# Human skin–cotton friction evidence

The calibration API exposes two measured cohort estimates from [Temel, Johnson and Lloyd (2022)](https://link.springer.com/article/10.1007/s11249-021-01560-5), rather than a body-wide friction constant.

| Region | Static mean ± SD | Dynamic mean ± SD |
|---|---:|---:|
| Chest | 0.90 ± 0.19 | 0.79 ± 0.20 |
| Dorsal forearm | 0.31 ± 0.07 | 0.25 ± 0.05 |

These dimensionless values are numeric text in §3.2 and the abstract. Ten participants (3 male, 7 female; age 22.3 ± 2.9 years) contacted pre-shrunk single-jersey 100% cotton, 0.62 mm thick and 140 g/m². Conditions were 25°C, 50% relative humidity, 2 ± 0.5 N, 0.04 m/s, 40 mm strokes, and superior motion. Subjects were supine after ten minutes' acclimation; skin received no special cleaning or treatment. Table 1's regional temperature and hydration measurements accompany each estimate. Table 2 reports repeatability ICC, not friction coefficients. Other regions' plotted coefficients have not been digitized.

## API and evidence boundaries

```python
from ihm.calibration.textile import measurement_conditions, textile_friction
estimate = textile_friction(
    'chest', fabric='temel_2022_cotton_single_jersey',
    conditions=measurement_conditions())
```

The complete protocol is a required evidence selector. It does not certify that a calling simulation matches the experiment. Returned SD describes reported cohort dispersion; it is not a posterior, confidence interval, or a proposed sampling law. Static/dynamic covariance is unavailable. Unsupported regions, fabrics, modified or incomplete protocols raise `ValueError`; there is no wetness scaling, genital extrapolation, or fallback coefficient. Values are returned as fresh copies.

## Acquisition and verification

`data/sources/human-skin-textile-friction.json` records the DOI, German National Library PDF acquisition URL, UTC receipt, license, extraction commands, byte counts and SHA-256 hashes. The acquired PDF and `pdftotext -layout` output remain under ignored `data/raw/textile/temel_2022/`; the numeric transcription is committed in `data/measurements/textile/temel_2022.json`. Reacquisition uses the source card command; extract text with its recorded command. Different PDF or extraction bytes should trigger review rather than silently updating the receipt.

Run `OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_textile.py`. It verifies acquisition hashes, exact numeric values, unsupported requests and defensive copies. It also calls the existing `PlaneContact.resolve` with 0.05 kg, tangential speed 0.02 m/s and inward normal speed 0.04 m/s: the requested impulse ratio is 0.5. Chest ends at zero tangential velocity; dorsal forearm remains at 0.01 m/s. Equal/opposite reaction impulses are checked. Solver source coefficients and behavior are unchanged. Results are retained in `data/derived/audits/textile/verification.json`.

That impact calculation transfers measured sliding coefficients into an explicit Coulomb model. It is a deterministic contact demonstration, not a replication of the tribometer's sustained loading or a validated garment/body interaction. Cloth stiffness, weave mechanics, sweat effects, pressure dependence and other regions remain outside this dataset.
