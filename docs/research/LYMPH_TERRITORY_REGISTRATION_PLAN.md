# Skin drainage territory registration: source evidence and executable boundary

The retained canonical assembly cannot yet identify anatomical drainage territories
for the installed native Skin regions. The new source-only inventory keeps every
skin triangle unresolved. It makes future annotation testable without changing the
fixed native .2/.3/.5 engineering installation or claiming anatomical identity.

## Retained geometry audit

`data/research/lymph_registration/canonical_audit.json` records byte lengths,
SHA-256 identities, licenses, transforms and source lineage. Reproduce it with
`.venv/bin/python scripts/audit_lymph_registration.py` (about 0.3 seconds, no network
or native execution).

| Retained structure | Audited inventory | Meaning |
| --- | --- | --- |
| BodyParts3D skin | One surface, 203,382 triangles | Material-index inventory; no lymph watershed masks |
| Z-Anatomy node groups | 158 named geometry objects; all geometry hashes verified | Authored group geometry, not 158 measured nodes |
| Published lymph graph | 996 vertices, 272 tagged lymph nodes, 1,117 edges | Undirected authored computational topology; zero named graph nodes |
| Nearest-group associations | 242 accepted spatial candidates out of 272 | Distance acceptance is not drainage or node identity |

The skin frame is metres, +x left, +y superior, +z anterior. The full geometry
receipt is `e724b62e6fc21764a0cf79c1874e96b9852755dc420711360a02d0b7afd93d05`.
Its 5,069,453 retained compressed bytes were read locally, not downloaded again.
The graph pose registration includes transverse/anterior scale and arm-pose priors;
no spatial threshold promotes a candidate to a sentinel-node connection.

BodyParts3D skin provenance reports CC BY 4.0. The named Z-Anatomy node objects
report CC BY-SA 4.0, with BodyParts3D-derived ancestry and third-party attribution
retained separately. They are not independent subjects. The published graph's
article declares CC BY 4.0; its supplement has no separate license file, and that
article declaration is not asserted for the original commercial PlasticBoy assets.
All already-retained publisher files and archive members have been hash checked.
See the [publisher model](https://doi.org/10.3390/math8122236) and the exact local
license evidence in `data/raw/lymphatic/provenance.json`.

## Bounded primary-source findings

Four source cards and retrieval receipts are in
`data/research/lymph_registration/`. Direct PubMed retrieval returned cookie-challenge
HTML with HTTP203; these failed bytes remain identified as failures. A bounded
Europe PMC metadata request retained the primary-paper abstracts successfully.
The four records report no open-access license; public availability does not license
reuse of figures. No full articles, figures, CT scans or watershed masks were
acquired, and no diagrams were digitized.

| Primary study | Registration implication |
| --- | --- |
| [Suami et al., upper limb, 2007](https://pubmed.ncbi.nlm.nih.gov/17440362/) | Preserve bypass alternatives to the dominant axillary destination; do not force all channels into one chain. |
| [Suami et al., upper torso, 2008](https://pubmed.ncbi.nlm.nih.gov/18349641/) | Midaxillary landmark is a candidate annotation constraint; anterior and posterior upper-torso destinations differ. |
| [Shinaoka et al., lower limb, 2020](https://pubmed.ncbi.nlm.nih.gov/31746690/) | Four pathway groups motivate a lower-leg template. Posterolateral drainage reaches popliteal nodes; inguinal destinations must retain shared versus distinct groups. |
| [Granoff et al., upper limb, 2022](https://pubmed.ncbi.nlm.nih.gov/35939638/) | Forearm-to-upper-arm connections vary. Pathway occurrence frequencies cannot be used as fluid allocation fractions. |

These studies constrain semantic destination alternatives. They do not establish
canonical triangle masks or subsequent node-to-node chains. Source-backed chains
currently stop at the reported first basin/group; continuation to iliac, lumbar,
thoracic or venous endpoints remains an explicitly unresolved edge until separately
supported. In particular, a popliteal endpoint is not matched automatically to the
canonical **deep** popliteal object by its nearest position or similar name.

## Parameterizable partition design

1. **Material ownership:** choose one pinned whole-skin mesh. Each triangle belongs
   to exactly one territory ID or `unresolved`; reject duplicates, holes disguised
   as exclusion, and source-hash changes. Boundary triangles require deterministic
   refinement with parent-triangle barycentric provenance before sharing. Never
   assign one full triangle to two territories. Initial candidate templates are
   side-specific lower-leg AM/AL/PM/PL, upper-torso anterior/posterior, and named
   forearm pathway catchments. These are proposed annotation slots, not implemented
   masks. Head/neck, thigh, genital and remaining skin stay unresolved.
2. **Registration evidence:** every submitted mask carries source DOI, observation
   type, subject/side, coordinate transforms, landmark correspondences, uncertainty,
   annotation author/revision and geometry hash. A curated node crosswalk stores
   source node/group identity, canonical candidate IDs, laterality and superficial/
   deep distinctions. Name or proximity alone cannot validate it. Multiple drainage
   alternatives may belong to one material territory; they do not duplicate skin.
3. **Exclusive native stores:** a future generalized installer takes arbitrary
   territory IDs and explicit installation weights. A proposed engineering prior
   is `w_i = sum(A_triangle * thickness * extracellular_fraction)` normalized across
   **all** territories including unresolved. Thickness and extracellular fraction
   are parameters with evidence/uncertainty, not measured values supplied by this
   audit. Area-only allocation is the explicitly uniform-parameter special case.
   Partition baseline/current/next volume and compliance with residual correction;
   preserve aggregate species current mass exactly as the existing installer does.
   Skin vascular/intracellular stay shared. A receiving node ID has one fluid/protein
   owner regardless of how many skin territories drain into it. No node stores are
   added here: current native Lymph remains the single lymph owner.
4. **Flow topology:** store validated directed source edges separately from the
   undirected display graph. Each edge records endpoint identity evidence; unknown
   direction, radius, resistance, valve law, compliance and protein transport remain
   unknown. Only a new tested native node-chain extension can create these stores
   and laws. Chemical feedback remains aggregate until independently implemented.
5. **Force/contact mapping:** map contact patches through pinned material triangles,
   clip patch areas, and preserve the summed force vector and normal load. Aggregate
   normal load divided by total territory area yields an engineering whole-territory
   surface-pressure average with unloaded facets included. Preserve patch-load
   variance and residual moments in a future mechanics adapter; one scalar loses
   spatial gradients and does not identify interstitial pressure. Do not send this
   reduction to existing engineering native IDs. Full force/moment and virtual-work
   parity, pressure-transmission evidence and generalized native lifecycle acceptance
   are separate activation gates.

`partition_template.json` makes the first eight side-specific lower-leg annotation
slots machine-readable. Each has null triangle membership, native binding, volume
weight and hydraulic parameters. Named canonical basin candidates are explicitly
unaccepted; AM/PM shared-destination constraints are preserved without creating
duplicate node owners. The template cannot activate a native variant.

## Implemented source-only validation

`SkinTerritoryInventory` in `ihm/assembly/lymph_registration.py` verifies the mesh
hash, positive finite triangle areas, submitted mask exclusivity and explicit
unresolved complement. It rejects native binding fields, so it cannot relabel the
installed engineering compartments. `report()` always states that anatomical
validity is unestablished; passing geometry accounting is not an anatomy review.

`reduce_contacts()` preserves supplied vector force, area and supplied compressive
normal load by triangle owner. It fails closed for multiple patches on one triangle
until a geometric disjointness proof is implemented, excess area, negative pressure,
invalid vectors or duplicate triangle IDs. The vector force and normal pressure are
separate supplied observations: this utility does not infer their normal orientation
or establish their physical consistency. Output contains no native command and
claims no interstitial-pressure identification.

`.venv/bin/python scripts/verify_lymph_registration.py` passes four bounded groups:
actual canonical/retained-source receipt checks; complete inventory and force/area
conservation; corrupt/overlapping masks and native rebinding rejection; invalid
contact handling. Synthetic masks exercise the validator only. No actual drainage
mask, clinical calibration, new native store or whole-body chain was fabricated.
