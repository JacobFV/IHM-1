# Skin shell eligibility and layer volume

The retained skin mesh has two dominant opposing components, plus small seam
fragments. Its raw triangle area is 3.5025989745 m². The selected largest
positive-orientation component has 109,183 triangles and area 1.7804602548 m².
It has open boundaries; selecting it as the exterior is an explicit engineering
prior, not a closed-volume proof, clinical body-surface measurement or verified
self-intersection-free surface. The original 203,382 triangles remain intact.

The earlier canonical layer generator multiplied all source triangle area by
uniform layer thickness. This counted both dominant shells as independent layer
support. The generator now requires the hash-bound component evidence in
`data/research/engineered_skin_territories/materialization.json`, verifies its
source geometry, exact component mask and selected triangle area, and uses that
area for synthesized epidermis, dermis and hypodermis volume. The original raw
surface area remains descriptive source metadata. Mechanical generation carries
explicit shell thickness and physical-support provenance into its output.

`ihm.assembly.skin_layers.physical_skin_support` performs these checks against
the exact geometry and evidence bytes. Five small tests cover a two-shell
counterexample, changed source, stale receipts, duplicated/misassigned face IDs
inconsistent area, and area overflow in excluded faces. The component assignments
and exterior-selection rationale are trusted engineering evidence; this checker
does not independently reconstruct topology or authenticate human review. The actual held-source calculation is retained in
`data/research/engineered_skin_territories/layer_volume_candidate.json` (0.29 s,
1 GiB address-space cap). At the unchanged 0.1/1.5/5 mm thickness priors, proposed
volumes are 0.000178046, 0.002670690 and 0.008902301 m³ respectively. These are
geometric priors, not measured tissue inventories.

The current canonical anatomy/mechanics artifacts have **not** been regenerated.
Migration must preserve the source skin, update dependent layer volumes and
mechanical mass-normalization receipts together, and separately validate each
consumer. It must not remove native patient mass: the native 22-body model is
scaled independently to patient mass, and canonical tissue proxies are not
independent native inventory stores. Existing frozen acceptance fixtures remain
historical evidence.

The supine quadrature uses minimum-X ray intersections, rather than summing all
faces; the existence of two shells alone does not prove doubled bed contact.
Its selected face membership is audited separately. Its legacy layer thickness
calculation divided old layer volume by raw source area, cancelling the old area
inflation. A corrected future artifact must use explicit shell thickness so this
ratio does not accidentally halve the contact material thickness. Hair roots
also require exact exterior-component eligibility; silently reattaching roots
would erase provenance.
