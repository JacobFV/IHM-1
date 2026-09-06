# Conditional muscle microvascular materialization

The new `ihm.assembly.vascular_priors` module materializes a bounded, seeded human-muscle-prior graph and exposes an optional passive network diagnostic. It uses the previously held human vastus lateralis Table 3 measurements. It adds no native blood volume, perfusion, transport law or feedback. Other organs remain explicitly unsupported for conditional capillary materialization.

## Registry and evidence

`data/sources/organ_microvascular_priors.json` preserves all four age/sex cohort summaries: capillary density, capillary/fiber ratio, capillary domain area and logarithmic domain-area SD, with species, region, biopsy scope, n, units, mean and between-participant SD. Actual observed ranges and distribution families are null because the source table does not supply them. Age comes from Table 1, with its different young-male sample count explicitly recorded. The source XML hash is verified every time the registry loads. The registry and both generator modules are hashed in each returned materialization.

These are population-transfer constraints, not individual human anatomy. The registry's skin/lung/kidney/liver rows state the evidence gap and `materialization_supported=false`; they do not borrow muscle coefficients. The existing measured kidney arterial materializer remains separate.

## Concrete callable path

```python
from ihm.assembly.vascular_priors import materialize_muscle, solve_passive

graph = materialize_muscle(
    '.', cohort='young_men', section_mm=(0.2, 0.2),
    extent_mm=0.5, seed=19,
)
# 331 capillaries/mm² × 0.04 mm² -> 13 transverse crossings.
# Achieved density: 325/mm²; residual: -6/mm².
assert graph['radius_m'] is None
assert graph['native_owner'] == 'Muscle'
assert not graph['independent_store']

# Optional diagnostic only: every coefficient below is a caller scenario.
flow = solve_passive(
    graph['edges'], [1e15] * len(graph['edges']),
    pressure_boundaries_pa={0: 4000., 1: 1000.},
)
```

The target defaults to the selected cohort mean. A caller may provide a positive `density_per_mm2` override; it is labeled as a scenario target. There is no invented distribution family or automatic Gaussian sampling from summary mean/SD. The seed determines transverse geometry only.

The module reuses the existing `microvascular_unit` paired supply/return topology, but **discards its arbitrary radius coefficients**. Its transverse coordinates are scaled to the caller's section rectangle. It checks actual edge intersections of the x=0 plane, confirming the count and that crossings are capillary edges. Count is nearest-integer target density × area; the reported residual is bounded by half a capillary divided by section area. Positive finite sizes, densities and seeds are required. Requests outside 1–128 crossings fail instead of truncating an intended whole-body graph.

The graph is donor-local synthetic geometry. Its paired binary branching, longitudinal extent and placement are engineering assumptions; the mean-density constraint does not validate those dimensions. `radius_m` and `flow_solution` remain null in the materialization. Anatomical registration and correspondence to native terminals are unresolved.

## What is and is not closed

Every generated capillary connects the inferred supply/return topology. The optional `solve_passive` accepts explicit positive edge resistances in Pa·s/m³ and finite pressure boundaries in Pa. It rejects floating components and computes nodal pressures and edge flows using Kirchhoff conservation. It reports internal flow residual and signed boundary outflows. The analytic test uses two parallel routes to confirm pressure, flow signs and total boundary cancellation; a generated-graph test checks internal conservation.

That mathematical closure is conditional on caller-supplied resistances and boundary pressures. It is not evidence of actual muscle resistance, radius, hematocrit, oxygen delivery or flow. The diagnostic has no native feedback and does not create a native perfusion owner. No library coefficient supplies an unstated viscosity or pressure assumption.

Domain-area distributions, capillary/fiber ratio and spacing heterogeneity are retained in conditioning and explicitly listed in `unmet_constraints`. Matching density alone cannot satisfy their joint distribution; this increment does not claim otherwise. The next geometry improvement requires fiber cross sections and a stated spatial-placement objective, followed by independent validation of domain statistics. Longitudinal topology, radii and exchange laws still require additional human evidence.

A runtime consumer may attach the result to native anatomical supports only with an explicit registration and conservative quantity-allocation contract. This module provides geometry and diagnostics; it does not call `body_exchange`, change its partitions or double-count native volume/mass/flow.

## Verification

`.venv/bin/python scripts/verify_vascular_priors.py` passes eight bounded tests: source values/units, source-tampering rejection, deterministic density and missing-field semantics, invalid/oversized requests, numeric overflow rejection, analytic passive closure, generated-graph closure and unanchored-component rejection. Tests were first observed failing for the absent module; the overflow check also failed before validation was added. No native engine, browser, large graph or new dataset download was used.

## Explicit radius and flow realization

`materialize_muscle(..., hydraulic_scenario=...)` now optionally produces definite positive radii, geometric lumen volumes, cylindrical resistance and a solved pressure/flow field. The option requires every hydraulic assumption; omission still returns the earlier geometry-only result. There is no hidden measured-radius default.

[Fine structure of the human skeletal muscle capillary, Gidlöf et al. 1988](https://pubmed.ncbi.nlm.nih.gov/3350624/) provides primary quadriceps-femoris biopsy diameter measurements perpendicular to muscle fibers. Major/minor lumen diameter means are **3.82/2.61 µm before estimated preparative-shrinkage correction**, and **5.31/3.62 µm after correction**. Reported ± quantities are 0.82/0.72 µm and 1.14/1.00 µm respectively. The held abstract does not identify SD versus SE, sample n or major/minor covariance; these remain unknown. Its network volume/surface estimates use assumed section density, so they were not promoted to independent measured constraints. No human length distribution was acquired. Rat length/corrosion-cast diameter reports located during search were excluded from human calibration.

The option `human_quadriceps_1988_equal_area` transforms the selected pair of diameter means into `r = sqrt(D_major_mean * D_minor_mean)/2`. This is an area-equivalent circle formed from two means, **not** an observed mean radius or a hydraulic equivalence to a flattened capillary. Corrected means give approximately 2.1922 µm. The transfer combines different human quadriceps specimens with the separate vastus-lateralis density cohorts; it is conditional synthesis, not same-donor anatomy.

[Structural Microangiopathies in Skeletal Muscle Related to Systemic Vascular Pathologies in Humans, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7013089/) was also acquired. It explicitly warns that glutaraldehyde-fixation shrinkage cannot be excluded and reports a factor-two correction of prior radius calculations. Its human morphometry and mouse ultrastructural validation are not pooled. The 2020 figure's radius values were not digitized. This prevents an unjustified import of earlier radius numbers merely because a paper labels them “radius.” The primary payload hashes and byte sizes are in `data/research/organ_microvascular/radius_acquisition.json`; 1988 bibliographic/abstract access is not treated as permission for unrestricted full-article reuse. The 2020 article is CC-BY.

```python
scenario = {
    'radius_mode': 'human_quadriceps_1988_equal_area',
    'diameter_state': 'shrinkage_corrected',
    'capillary_radius_cv': 0.2,
    'supply_radius_multiplier': 2.0,
    'return_radius_multiplier': 2.4,
    'viscosity_pa_s': 0.003,
    'pressure_boundaries_pa': {0: 4000., 1: 1000.},
}
graph = materialize_muscle('.', seed=13, hydraulic_scenario=scenario)
```

The CV, distribution family, supply/return multipliers, viscosity and pressure boundaries in this example are **engineering scenario choices**. None is attributed to the diameter study. `engineering_lognormal` is an alternative explicit mode requiring `capillary_radius_mean_m`; it permits a fully declared synthesis when source transfer is unsuitable. The implementation uses seeded lognormal shape then rescales the finite capillary sample to the chosen mean. It records realized mean, population-form SD, min/max and seed. Thus realized CV need not equal requested lognormal shape CV, especially for small fixtures; n=1 necessarily has zero realized SD. CV is bounded to [0,2] for this small numerical fixture. Constant supply/return radii are explicit multiplier assumptions with no inferred Murray-law claim.

Actual edge lengths come from the generated graph. Poiseuille resistance uses declared viscosity and realized radius/length, and `solve_passive` solves those resistances rather than a separately supplied resistance vector. Returned geometric volumes and lumen surface are geometric demand, **not independent blood stores**; they are neither allocated to nor added on top of the native Muscle owner. Native feedback remains disabled. The Newtonian circular-cylinder scenario omits RBC effects, phase separation, compliance and exchange; its numerical solution does not establish physiological accuracy.

The verification suite now has eleven passing tests, including deterministic source-conditioned moments, positive geometry, pressure bounds/forward monotonicity, Kirchhoff conservation, and the analytic scaling that doubling all radii multiplies flow by 16 and geometric volume by 4. An incomplete scenario is rejected. The retained bounded example `data/research/organ_microvascular/muscle_hydraulic_fixture.json` contains actual radii and solved flows with input/source/generator provenance.
