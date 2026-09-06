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
