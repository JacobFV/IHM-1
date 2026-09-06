# Measured mattress compression candidate

A usable finite-thickness engineering candidate is retained in `data/research/bed_material/`: Hong et al. (2022), *Biology* 11, 1030, [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9311775/), DOI 10.3390/biology11071030. Article and Figure 1 are CC BY 4.0; the manifest records author attribution, source URLs, bytes and SHA-256 identities.

The study measured three Sinomax foam mattresses, labeled SM/MM/HM, using an Instron 5569. Each mattress was 1.9 × 1.2 × 0.2 m; reported 25% indentation loads were 20/42/120 lbf. Figure 1b provides compression curves, and Table 1 uses Poisson ratio 0.01. Coupon size, test rate, density, temperature, conditioning and stress-measure definition are unreported in the inspected methods. No unloading curve or damping fit is supplied. Its friction coefficient 0.6 describes pillow–mattress contact, not skin–bed friction.

## Retained candidate

`hong2022-compression.csv` contains ten digitized points per curve, including a declared zero-load origin. `digitize_hong2022.py` reproduces the table from explicit manual pixel picks and verifies the source figure hash. The curves remain separate; no selection is made to force body equilibrium.

| Curve | Stress at final digitized strain 0.66449 (kPa) |
|---|---:|
| SM | 2.964 |
| MM | 8.183 |
| HM | 24.549 |

Manual line placement uncertainty of roughly ±3 pixels corresponds to ±0.387 kPa and ±0.0076 strain; these are screening estimates, not confidence intervals. The JPEG and medium/pillow curve overlap add uncertainty. This is a digitized figure, not raw instrument output.

Use monotone piecewise-linear interpolation on the retained domain only. For bed thickness `h=0.2 m`, strain is `delta_bed/h`, giving a maximum tabulated travel of 0.132899 m. Assuming the plotted stress is nominal, per-reference-area energy is `h * integral(p(epsilon), 0, delta_bed/h)`. The analytic integral of each linear segment makes the force–energy derivative check explicit. Reject out-of-domain compression; do not silently clip pressure or invent a densification branch.

Place this bed response **in series** with the existing bounded skin law: equal traction, `delta_total=delta_skin+delta_bed`. Bed displacement must not enter skin compression receptors. Preserve finite bed thickness, one physical world, and equal/opposite external wrenches. The local independent-column approximation still lacks lateral/shear coupling validation.

For conservative static feasibility, the loading curve can supply a declared elastic candidate; reusing it on unloading is an engineering simplification with zero hysteresis, not a measured unloading law. Do not infer damping, creep, friction or recovery times. Dynamic settling requires separately supported dissipative parameters. No supine equilibrium or clinical calibration follows from this material acquisition.

A second primary cyclic-foam paper, DOI 10.3390/psf2022004020, was identified, but publisher acquisition returned 403/429; no cyclic-law bytes or coefficients were retained from it.

## Verification

Run `.venv/bin/python data/research/bed_material/digitize_hong2022.py`. The retained run verified the source hash, 30 rows, finite compression domain, nonnegative monotone stress, and deterministic CSV identity. No native process, full anatomy asset, or simulation model was modified.
