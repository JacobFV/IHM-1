# Conditional human cutaneous microvascular patch

`ihm.assembly.skin_microvascular_patch.materialize_skin_patch` creates a bounded
arteriole → superficial plexus → papillary loops → venular plexus flow graph on
one exact retained skin triangle. Every response includes original face index,
source coordinates/hash, engineered territory ID/hash, seed, declared scenario,
SI centerlines/radii/pressures/flows, and the complete source-conditioning card.
It is not a recovered individual vessel graph.

## Held primary evidence

The new acquisition retains 165,587 bytes: complete confocal article HTML
154,092 bytes, its bibliographic record 6,837 bytes, and a primary dermal-plexus
record 4,658 bytes. `data/research/skin_microvascular/acquisition.json` records URLs,
byte lengths, SHA-256, retrieval date and reuse status. No raw vascular images,
subject graphs or source videos were acquired. The full-text XML endpoint returned
404; the complete public HTML succeeded. No open reuse license was identified;
free reading does not imply a permissive redistribution license.

[Altintas et al., DOI 10.1007/s10278-009-9219-3](https://pmc.ncbi.nlm.nih.gov/articles/PMC3046666/)
measured living healthy middle-volar-forearm skin by reflectance confocal imaging.
The selected baseline cohort had 11 adults (5 female, 6 male), age 22.4 ± 3.1 years.
Perfused-loop density was 7.04 ± 0.62/mm² and lumen diameter 9.59 ± 0.25 µm.
Reported variation is SD; individual vessel distributions are unavailable. Diameter
used two papillae in four fields; density used eight fields at the epidermal–dermal
junction. Instrument resolution was 0.4 µm lateral, 1.9 µm axial, maximum depth
350 µm. No fixation/shrinkage correction applies to these in vivo measurements.
The reported 61.09 ± 3.21 blood cells/min is not volumetric blood flow and is unused
in the hydraulic calculation. This is functional observed-loop density, not total
vascular length density. The lower-leg application is an explicit site transfer.

[Yen and Braverman, PMID 1249441](https://pubmed.ncbi.nlm.nih.gov/1249441/)
reconstructed normal human forearm papillary horizontal vascular connections using
1 µm plastic sections and ultrathin sampling at 10–20 µm intervals across 450 µm.
The retained abstract supports arteriole/capillary/venule identity, not this
particular graph or deep-plexus coefficients.

Previously held DERMA-OCTA remains metadata/workbook-only; no source image bytes
were silently converted into a measured graph. Its noncommercial license and
subject/scan accounting discrepancies remain in `ORGAN_MICROVASCULAR_EVIDENCE.md`.

## Exact geometry and limits

The territory manifest is pinned to SHA-256
`de669f96f6db5b150dd467f149464677e19f22ec2e7be5a77d895d873b3fef1a`
from `95d8dc4`. Only a triangle belonging to both the selected named lower-leg mask
and its original exterior-eligible inventory is accepted. Inner/seam faces and
unresolved territory IDs are rejected. The full original mesh is read and hashed;
no source simplification, deletion, retessellation or canonical mutation occurs.

A centroid-scaled triangle gives the requested physical footprint, 0.15–8 mm²,
with a finite-radius margin to each original face edge. All graph surface
projections lie in that convex footprint. If no explicit face is requested, the
first original face in the selected territory satisfying the size/margin checks
is chosen deterministically. Requested capillary count is round(density × area),
limited to 1–56; integer rounding error and realized density remain explicit.
Seed controls loop positions; exact source bytes and all inputs identify the patch.

Depth is an inward source-face-normal extrusion, not a closed-surface certificate.
The exterior proxy is open; physical dermal containment, folds, self-intersections,
vessel packing and patient-specific tissue depth are unverified. Nominal epidermis
ends at 0.1 mm, loop apex at 0.15 mm, superficial plexus at 0.3 mm, deep plexus at
1.4 mm, dermis ends at 1.6 mm. These are scenario inputs, not new measured priors.
Every full centerline plus radius is checked against that nominal dermal slab.
No vessel is placed in the nominal epidermis. Original source geometry and
engineered territory identity are preserved even though microscopic strata are
conditional.

## Physics and ownership

`default_scenario()` returns all engineering assumptions explicitly: 20 µm arterial
and 25 µm venular radii, viscosity 0.003 Pa·s, 4,500/1,500 Pa boundaries, the above
depths, and affirmative forearm-to-lower-leg transfer. Callers must supply this
scenario or an explicitly edited bounded copy. Capillary radius is half the source
mean lumen diameter; no SD is misused as an individual radius distribution.

Two deep arterial/venous branches connect separate superficial arterial/venous
star plexuses through ascending/descending vessels. Each intervening papillary
loop retains a 33-point curved centerline. Template topology and all noncapillary
radii are engineering choices. Segment length is the sum over the full retained
polyline; Poiseuille resistance uses that length and the realized radius. The
Kirchhoff solve returns every node pressure and signed flow, with boundary and
internal conservation residuals. This Newtonian scenario does not model red-cell
partition, autoregulation, exchange, charge, solutes, filtration or lymphatics.

`prospective_native_owner` is Skin. Macro vessel endpoints and native boundary
binding are null and `binding_ready` is false. The geometry lumen-volume sum is
only a diagnostic: additional native volume is exactly zero, no independent blood
store exists, no command is emitted. Regional lymph terminology in the engineered
territory does not establish vascular territories or create lymphatic vessels.
Actual native binding requires receipt-backed macro endpoints and an explicit
conservative allocation of the existing Skin owner before any feedback.

## Use and verification

```python
from pathlib import Path
from ihm.assembly.skin_microvascular_patch import materialize_skin_patch, default_scenario
patch = materialize_skin_patch(Path('.'),
    territory_id='left_lower_leg_anteromedial',
    scenario=default_scenario(), patch_area_mm2=1, seed=3)
```

Run `.venv/bin/python scripts/verify_skin_microvascular.py` for actual retained-face
membership and coordinate projection, measured density/diameter conditioning,
conservation, pressure bounds/dissipation, ownership, reproducibility, pressure
scaling, malformed/bounded requests and nominal slab failures. These are light
Python checks, with no native engine, browser, download or whole-body generation.
