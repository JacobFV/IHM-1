# Engineered source-skin territories and exterior eligibility

`materialize_engineered_skin_territories.py` produces eight exclusive lower-leg
quadrant masks on the retained canonical skin, using the retained left/right tibia
and fibula surfaces. The result is an **inferred geometry partition**, not measured
lymph drainage, a venous watershed, or individual anatomy. No native region is
rebound and no native command is emitted.

## Double-shell discovery and reusable membership

The original skin's 3.502598974493317 m² surface measure combines two large shells
and small fragments. All coordinates remain in the source-declared canonical metre
frame; its overall height is about 1.71947 m. No clinical body-surface-area target
or rescaling was applied.

| Index-connected component | Faces | Area (m²) | Boundary edges | Signed-volume integral (m³) |
| --- | ---: | ---: | ---: | ---: |
| 0, inferred exterior proxy | 109,183 | 1.7804602548390722 | 827 | +0.06971094613047656 |
| 1, inferred inner shell | 94,055 | 1.722073687808178 | 345 | −0.06627586508909047 |
| 2–99, small fragments | 144 | remainder | 340 combined | not interpreted as volume |

Component 1 has nested axis-aligned bounds and opposite orientation. No coincident
faces of either winding were found after coordinate rounding to 1e-8 or 1e-6 m;
776 coincident vertex positions exist. There are 1,512 boundary edges, zero edges
with more than two incident faces, and zero inconsistent orientations among paired
interior edges. **The surfaces are open:** the signed integral is not a closed
volume and cannot by itself prove exterior identity. General triangle intersections,
near-coincident patches and physical exclusivity were not tested. Selecting the
largest positive-integral component is an explicit engineering eligibility prior.

The frozen artifact is
`data/research/engineered_skin_territories/materialization.json`:

- SHA-256: `de669f96f6db5b150dd467f149464677e19f22ec2e7be5a77d895d873b3fef1a`.
- Bytes: 1,204,791.
- `source_receipts[0].sha256`: `e724b62e6fc21764a0cf79c1874e96b9852755dc420711360a02d0b7afd93d05`.
- `surface_diagnostic.face_component_ids`: component ID for each of the 203,382
  original source-face indices, preserving original order.
- `surface_diagnostic.selected_component_id`: 0.
- `contact_eligible_triangle_ids`: the 109,183 faces of that component.
- `surface_diagnostic.components`: per-component face count, area, bounds, open
  edges, integral and eligibility basis; all 100 components retained.

This manifest is reusable by contact and hair-root audits through exact face IDs.
An inner/seam attachment is ineligible under this prior; it is not automatically
moved to the exterior. Neither canonical skin bytes nor canonical synthesized layer
volumes were changed by this generator. Multiplying whole-mesh area by tissue
thickness would count both shells and is not justified by this audit.

## Reproducible mask recipe

The five geometry receipts (skin, bilateral tibia/fibula) retain original source
identity, transformation, byte length, SHA and license metadata. Each tibia's
vertex-PCA major axis is oriented superiorly. The projected fibula centroid defines
the lateral direction; a perpendicular axis oriented toward canonical +z defines
anterior. These are geometric reference frames, not traced lymph channels.

The default recipe uses the 5th–95th percentile interval of tibia vertices along
that axis, a 0.13 m maximum radius, same-side canonical x, and a 3 mm exclusion band
around quadrant planes. Every vertex of an assigned triangle must meet every
condition. Faces crossing any boundary remain unresolved rather than duplicated.
The percentile, radius and exclusion band are bounded, explicit parameters. Changes
produce a new materialization rather than silently retaining old registrations.

Only exterior-proxy faces are assignable. The eight masks contain **7,448 faces**
with combined area 0.1881951269286595 m². All other 195,934 source faces remain in
the complete unresolved material inventory; 94,199 of those are additionally
excluded from contact as inner/seam geometry. Exterior unresolved area is stored
separately so inner-shell area never dilutes the contact pressure denominator.

The four-pathway terminology comes from the retained primary study
[Shinaoka et al., Radiology 2020](https://pubmed.ncbi.nlm.nih.gov/31746690/), documented
in `data/research/lymph_registration/source_cards.json`. That study supports distinct
lower-leg drainage pathways, but does not establish these bone-axis planes,
percentiles, radii, masks or area weights. No node-chain identity is inferred.

## Contact reduction and prospective native configuration

`reduce_contacts(materialization, contacts)` rejects inner/seam source-face IDs,
then uses the existing exclusive inventory validator to preserve contact area,
normal force and full supplied force vector by material owner. The unresolved
pressure denominator uses only exterior unresolved area. One average surface
pressure per region remains an engineering reduction: it does not preserve local
load variance/moments or identify interstitial pressure. No native pressure port is
commanded, and supplied force/normal-pressure consistency is not inferred.

`three_region_proposal(materialization, first_ids, second_ids)` requires two
explicit, nonempty, nonoverlapping unions of named masks. The third region is the
entire remaining exterior proxy, including unassigned exterior faces. It returns
positive area fractions and explicit selected IDs, with **no installed native
binding**. For a caller explicitly choosing all four left masks and all four right
masks, the illustrative fractions are:

```
left selected lower-leg union   0.052809544820190284
right selected lower-leg union  0.052890723645319306
remaining exterior proxy        0.8942997315344904
```

These are source-geometry fractions under uniform thickness/extracellular-fraction
assumptions, not measured extracellular volumes. Selecting all four quadrants in
one union also loses their separate drainage distinctions. The generator's saved
artifact selects **no** such unions by default.

Three existing native circuit branches can in principle represent three explicitly
chosen unions while conserving the aggregate. Required work is a **new immutable
variant**, not editing or relabeling the current installation:

1. Replace `NativeRegionalSkin::fractions`' static .2/.3/.5 with validated positive
   configuration weights and receipt-bearing semantic IDs. `partition()` and
   `copy_laws()` must use that same immutable configuration throughout installation
   and runtime refreshes. Keep third-allocation residual correction.
2. Continue partitioning baseline/current/next volume, compliance and flow sources;
   divide resistance by each weight and preserve shared pressure sources. For
   positive weights summing to one, parallel conductance and summed compliance
   retain the aggregate source values algebraically. Partition every species mass
   exactly once; the original Skin owner remains an aggregate. This establishes
   an allocation identity, not full nonlinear engine parity.
3. Update the separate adapter's exported fraction/ID receipt, command validation,
   and Python observer's fixed `FRACTIONS` check. A config digest must bind source
   face masks, eligible surface, selected unions and engineering volume assumptions.
   Mapping named unions to engine indices must be an explicit configuration entry.
4. Repeat actual native installation conservation, zero-load matched-engine parity,
   load/unload, sweat owner/waste closure and graph/cache lifecycle tests for the new
   weights. Smaller branches increase resistance and can affect conditioning; the
   existing .2/.3/.5 acceptance does not validate them. Skin intracellular/vascular
   remain shared; native Lymph is still lumped and no node chain is created.

Verification: `.venv/bin/python scripts/verify_engineered_skin_territories.py`
passes four groups in about four seconds, including byte-reproducible actual-source
masks, exclusive inventory, bounded parameter sensitivity, explicit conservative
union proposals, exact material contact incidence and inner-face rejection.
No downloads, canonical rebuild or native execution are part of these checks.
