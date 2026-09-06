# Conditional kidney microvascular beds

`ihm.assembly.kidney_microvascular_patch.materialize_kidney_patch` creates a bounded kidney-specific graph spanning afferent supply, vascular chambers, conduits, anastomotic glomerular loops, efferent drainage and a downstream peritubular mesh. It uses new human primary evidence and does not reuse the muscle graph generator. It is an explicit local hydraulic prototype, not a recovered donor glomerulus or native renal solver.

## Held evidence

`data/research/kidney_microvascular/acquisition.json` records five acquired primary/bibliographic payloads totaling **13,913,917 bytes**, with URLs, UTC timestamps and SHA-256. No raw human glomerular or PTC 3D graph was acquired. The 2018 tables are also extracted to JSON, and the machine-readable registry is `data/sources/kidney_microvascular_priors.json`.

[Neal et al., Novel hemodynamic structures in the human glomerulus](https://pmc.ncbi.nlm.nih.gov/articles/PMC6293306/) studied unused human transplant kidneys. Fourteen resin-reconstructed glomeruli from four kidneys supply these measurements (mean ± **SE**, not individual-vessel SD): afferent arteriole diameter **21.5±1.2 µm**, afferent vascular chamber **43.2±2.8 µm**, afferent conduit **15.9±0.7 µm**, efferent first-order vessel **9.9±0.4 µm**, efferent chamber **38.4±4.9 µm**, and efferent arteriole **15.9±1.2 µm**. Afferent conduit counts average 6.6±0.6 (observed 2–11); efferent first-order counts 12.6±1.4 (3–22). These are actual human measurements under the study's perfusion/fixation conditions. Chamber diameters average three ellipsoidal axes; circular tubes in the prototype are a geometric approximation.

The study used nine unused transplant kidneys overall, with different imaging subsets. It distinguished fresh and fixed aqueous imaging from resin reconstruction. In-plane glomerular max/min measurements pooled from 14 glomeruli were 194.4±5.1 µm (n=28 measurements). The prototype does not claim to reproduce that diameter or packing. The paper's **3.5 µm filtration-capillary radius appears in a hoop-stress calculation (Table A1)**; the registry retains it as a published model assumption, not a measured radius distribution. Article XML is 172,769 bytes, CC-BY-4.0.

[Shimizu et al. 1988](https://pubmed.ncbi.nlm.nih.gov/3185306/) reconstructed portions of five human glomeruli from each of three cortical zones and analyzed cycle rank. Its held abstract supports cyclic topology and regional caution but provides no transferable numeric cycle-rank distribution or downloadable graph. [Bohle et al. 1998](https://pubmed.ncbi.nlm.nih.gov/9736285/) describes seven lobule-like structures with anastomoses and 0.95 cm total capillary length per glomerulus. Its abstract gives no variance or specimen-level graph. Applying that aggregate length to the prototype's filtration-edge budget is explicitly a transfer assumption because the abstract does not resolve inclusion of conduit segments. Bibliographic payloads are 3,787 and 5,670 bytes; full-article reuse licenses were not established.

[Reconfiguration and loss of peritubular capillaries in chronic kidney disease, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10640592/) includes human CD31-stained biopsies (10 allograft and 30 IgA-nephropathy cases). Preserved and injured regions are not separate healthy donors. Human measurements are 2D shape statistics; the article's reconstructed 3D networks are **mouse**, not human. The article XML is 144,131 bytes and its held supplementary PDF is 13,587,560 bytes, CC-BY-4.0. Table S6 supplies actual human 2D measurements: preserved-region luminal area **71.0±63.3 µm²**, ellipse major/minor diameters **12.8±4.1 / 5.0±1.8 µm**, and area fraction **10.3±3.2%**; corresponding injured-region values are **45.0±46.9 µm²**, **9.2±3.8 / 4.8±2.1 µm**, and **7.4±2.5%**. These are group mean±SD of biopsy/region mean measurements, not individual-vessel radius distributions. PDF text extraction is held separately. No mouse topology or radius is relabeled as human.

## Callable materialization

```python
from ihm.assembly.kidney_microvascular_patch import materialize_kidney_patch
scenario = {
    'viscosity_pa_s': .003,
    'peritubular_radius_mode': 'human_control_area_equivalent',
    'arterial_pressure_pa': 12000.,
    'venous_pressure_pa': 1000.,
    'bowman_pressure_pa': 2500.,
    'interstitial_pressure_pa': 500.,
    'glomerular_hydraulic_conductance_m3_s_pa': 1e-17,
    'peritubular_hydraulic_conductance_m3_s_pa': 1e-18,
}
patch = materialize_kidney_patch('.', owner='RightKidney', seed=11, scenario=scenario)
```

The PTC mode maps the selected human Table S6 mean area to `sqrt(mean_area/pi)`, giving a 4.7539459 µm equal-area circular radius for preserved regions. This is not measured mean radius or hydraulic equivalence to elongated PTCs. `human_injury_area_equivalent` selects injured regions; an explicit `peritubular_radius_m` instead permits a declared engineering radius. All pressure, viscosity and exchange-conductance numbers in the example are engineering choices. The constructor requires the complete scenario and rejects unsupported owners or invalid values. `LeftKidney` and `RightKidney` select native quantity ownership only; no anatomical side registration is fabricated.

The graph uses seven lobular loop/chord networks, seven afferent conduits and fourteen first-order efferent connections. Loops and interlobular anastomoses provide positive cycle rank; they are newly synthesized, not measured individual connectivity. Arteriolar/conduit/chamber radius means come from the human registry. Peritubular geometry is an explicitly inferred 3×4 lattice downstream of the efferent arteriole. This preserves the serial vascular beds, unlike an afferent-to-venous paired muscle tree.

Each filtration edge has an actual 3D sinusoidal polyline, with amplitude solved to meet its share of the 9.5 mm total length target. Length is computed from polyline increments and used in resistance, not from endpoint chords. These winding paths are engineering geometry, not measured tortuosity; no collision-free packing, diameter-envelope or individual topology validation is claimed. Straight segments describe other vessel classes. Returned `positions_m`, `edges`, `edge_polylines_m`, `radius_m`, `edge_length_m` and edge kinds provide explicit geometry for later rendering.

## Pressure and exchange boundaries

Poiseuille resistance uses the realized radius and polyline length with declared Newtonian viscosity. Blood passes through afferent arteriole → afferent chamber → conduits → glomerular loops → efferent first-order vessels → efferent chamber → efferent arteriole → peritubular mesh → venous exit. Source and sink nodes of that sequence remain explicit.

Two additional external pressure ports preserve Bowman's-space and interstitial boundaries. Positive total conductance adds equally allocated hydraulic shunts from glomerular or PTC nodes to the relevant pressure port; zero conductance leaves that named port disconnected while retaining its pressure metadata. Internal Kirchhoff conservation is solved across all blood and exchange edges. `hydraulics` reports the full diagnostic pressure/flow solution; `physical_node_pressure_pa` and `physical_edge_flow_m3_per_s` select the actual blood-geometry arrays. `exchange_flows_m3_s` separately accounts for hydraulic transfer to each external pressure port.

This is a purely hydraulic demonstration. There is no oncotic model, membrane selectivity, charge or solute flux, and the result is not calibrated GFR or tubular reabsorption. `solute_fluxes`, `charge_fluxes` and `oncotic_model` remain null. The exchange conductances are not inferred from geometric surface area. Filtration fraction is explicitly uncalibrated.

All volumes remain geometric demand. `independent_store=false`, `native_volume_allocated_ml=null` and `native_feedback_enabled=false` prevent interpreting the prototype as another native fluid store or feedback path. The 2018 specimens are not the same donor as the separately held LADAF-2021-17 arterial example; no arterial-to-glomerular registration is claimed. App/API integration remains deferred until an organ-generic contract is agreed.

## Verification

`OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_kidney_microvascular_patch.py` passes six bounded tests: source-conditioned deterministic topology and target total length; four-boundary pressure bounds, exchange direction and Kirchhoff conservation; polyline-versus-reported physical lengths and positive radii; invalid scenario/owner rejection; human PTC mean-area conditioning; zero-exchange and serial-bed disconnection checks. Tests first failed for the missing module. The retained example `data/research/kidney_microvascular/conditional_patch.json` contains definite geometry and a solved hydraulic state, with registry and generator hashes.
