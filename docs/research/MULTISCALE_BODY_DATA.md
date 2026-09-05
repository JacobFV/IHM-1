# Organ-specific conditional microstructure

The user's September5 steering extends the implicit/explicit architecture to unresolved body microstructure. Brain vascular and cortical acquisition remains owned by the concurrent IBM-1 work. IHM keeps its present imported package pinned; future IBM revisions must produce new byte receipts, compatibility/parity checks and explicit migration records. Live donor edits cannot silently change an existing body materialization.

Each body microstructure component distinguishes measured geometry in the source specimen, population/other-specimen evidence transferred to the generic body, and conditional synthesis. All three can participate in one substrate. A materialization records which observations constrained its structure, unsupported dimensions, deterministic seed, scale, registration and resolution censoring. Biological variation, inference uncertainty and display precision remain separate.

Matching radius/length histograms does not establish equivalent physiology. Conditioning and validation must include the observable of interest, correlations and topology: connectivity, permeability, transit-time distributions, exchange area, oxygen extraction, fluid/protein return or mechanical response. Murray scaling is a candidate statistical constraint, not an exact universal law imposed on every vessel bed. Transferring a cortex-specific generator to skin, glomeruli, muscle or liver is not justified.

| Tissue/domain | Microstructure and conditioning | Required functional constraints |
|---|---|---|
| Skin/subcutis | Dermal plexuses/loops, layer depth, body region, appendages, endothelial/lymphatic phenotype | Perfusion, exchange, tissue water/protein, lymph return, heat transport and electrical barriers |
| Kidney | Renal/interlobular arterial tree, cortex/medulla, glomerular afferent/efferent pathways, vasa recta/tubules | Filtration, oxygen delivery, pressure losses and serial/parallel exchange ownership |
| Skeletal muscle/tendon | Fiber/fascicle orientation and size, capillary alignment, motor-unit territories, tendon entheses | Force/strain, oxygen extraction, fatigue and conduction; no independent random attachment topology |
| Liver | Portal/arterial supply, sinusoidal network and central drainage, lobular zones | Dual-inlet transport, zonal oxygen/metabolism and hepatic venous return |
| Lung | Airway/acinar/alveolar organization, septal capillary sheet and recruitment | Regional ventilation/perfusion, gas exchange and pressure-volume behavior |
| Heart | Myocardial fiber/sheet orientation, coronary beds, conduction territories | Active work, coronary oxygen delivery, activation timing and ventricular coupling |
| Gut/pancreas | Crypt/villus architecture, enteric/muscular layers, islets/ducts | Absorption, secretion, motility and portal transport |
| Lymph/immune | Initial capillary connectivity, collectors/valves and nodes with compartment types | Uptake, pumping, one-way flow and fluid/protein/cell balance |
| Peripheral nerves/sensors | Fascicles, axon classes, terminal density and target assignments | Modality-specific conduction, selective blocks and effector response |

Initial concrete acquisitions/constraints are the regional human hair-density table, retained native skin/lymph coefficients, a reduced public human kidney arterial graph example, and a paired image/dense-label slab. The published final graph remains unavailable at its cited DOI. `scripts/collect_kidney_microstructure.py` records actual acquisition status; a source card alone is not held data. The [2025 kidney study](https://www.nature.com/articles/s44303-025-00090-2) identifies one donor and arterial segments down to interlobular scale, with substantial finer topology unresolved and explicit departures from simple vascular scaling. These observations must not be labeled a complete capillary or venous tree.

The [kidney acquisition record](KIDNEY_MICROSTRUCTURE.md) retains 148 nodes, 147 arterial polylines and 19,924 measured points, plus 32 consecutive image/label pairs at a documented 50 µm voxel spacing. Same-donor identity is established; graph-to-image and generic-body registration remain unresolved. The implicit body exposes this evidence without silently transferring it into a population prior:

```python
from ihm.human import ImplicitHuman
human = ImplicitHuman.open()
evidence = human.microstructure_evidence()
graph = human.materialize('kidney-arterial-geometry')
# graph['nodes_m'], graph['points_m']: original donor-local coordinates in SI
# graph['body_registered'] is False; radius_m and flow_solution are None
```

Evidence manifests pin both the binary graph and its interpretive source card. The raw thickness field remains uninterpreted because radius/diameter semantics are unverified. Geometry materialization neither generates missing vessels nor treats an ex vivo arterial example as generic-body capillary ground truth.

The [Human Organ Atlas](https://human-organ-atlas.esrf.eu/) offers multiscale ex vivo organ images; its [repository paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12324500/) describes nested whole-organ and local resolutions and multi-organ donors. Preserve donor/condition identities to support later joint inference; do not combine unmatched specimens as a measured single body. Large acquisitions are allowed by the user's accuracy-first scope, with checksummed raw inputs and derived products retained separately.

Pathological topology, unusual boundaries, local lesions and unobserved variation require scenario-specific validation. A synthesized network is one conditional realization, not recovered individual anatomy or a calibrated patient twin.
