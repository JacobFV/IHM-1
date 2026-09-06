# Native regional lymph and interstitial coverage

`ihm.assembly.lymph_coverage.NativeRegionalLymphCoverage` provides a read-only,
conservative anatomical view of the existing native tissue ports. It creates no
solver, fluid reservoir, protein reservoir, transport law, or extra feedback loop.

## Integration contract

```python
from ihm.assembly.lymph_coverage import NativeRegionalLymphCoverage

coverage = NativeRegionalLymphCoverage.from_workspace(root, native_identity)
report = coverage.observe(native_snapshot)
```

The snapshot is the existing `{time_s, values}` native acknowledgment. Optional
`organs` selects the standard native organ names. Construction requires an explicit
native identity in the workspace path; the direct constructor accepts tiny fixture
anatomy and optional identity, hashes and topology. The returned schema is
`native_regional_lymph_coverage_v1`. Inputs are not mutated. Workspace construction
checks and retains hashes of anatomy, support geometry, native port definitions,
the mapping implementation, and available canonical/source lymph graph receipts.

| Report field | Meaning |
| --- | --- |
| `native_owners` | One authoritative record per selected organ extracellular compartment plus the single `Lymph` compartment; native volume, Albumin mass, pressure and Albumin concentration |
| `regional_views` | Non-owning allocations of extracellular volume and Albumin mass to existing anatomical surface supports |
| `allocation_conservation` | Per-owner volume and Albumin allocation residuals; null when the corresponding native observation is unavailable |
| `pathways` | Organ E3-to-L1 drainage, L1-to-L2 and return-valve flows, plus `LymphToVenaCava`; absent ports remain unknown |
| `topology` | Source graph availability distinguished from absent native regional mapping, regional lymph stores and regional flow solution |
| `unlocalized_native_owners` | Lymph plus selected organs without a registered anatomical support |
| `unobserved_source_keys` | Exact missing or NaN source ports; never substituted with zero |

Every quantity retains its native source key and unit (`mL`, `g`, `mmHg`, `g/L`,
or `mL/s`). Regional quantities are explicitly classified as allocated outputs;
they are not independently measured or integrated states. Flow signs are retained.
Infinity, negative mass/volume, booleans and nonnumeric quantities are rejected.
Native null/NaN availability becomes JSON null with `unobserved` status.

## Ownership and anatomical assumptions

`body_microstructure.anatomical_supports` supplies the existing organ-to-anatomy
rules. Normalized source surface areas distribute each extracellular owner's
observed volume and Albumin mass. The final support receives the floating-point
remainder, preserving the native total. This is a display allocation prior, not a
claim of measured drainage territories, perfused tissue volumes, disjoint atlas
meshes, or regional pressure gradients. No matches means one explicitly unlocalized
view retaining the whole owner's quantity. Duplicate anatomy IDs and duplicate or
unknown native organ selections fail closed.

These views overlap the existing `body_exchange` extracellular partitions. A
consumer must choose a representation for an owner; it must never sum the native
owner with its regional views, or combine these views with the existing exchange
partitions as extra stores. `accounting_owner` and `independent_store` make this
contract machine-readable. Lymph itself remains one unlocalized native compartment
and receives no anatomical allocations. VenaCava is a destination of the observed
return path, not an additional reservoir copied into this report.

Albumin is the available protein species, not total protein. The mapper reports no
protein flux, reconstructed concentration-times-flow transport, or whole-body mass
closure. Serial lymph pathway flows must not be added as separate drainage volume.
A snapshot flow is not integrated transported volume.

## Held coverage and explicit gaps

The held anatomy currently gives localized supports to Bone (257), Brain (24), Gut
(60), each Kidney (1), Liver (10), LeftLung (2), RightLung (3), Muscle (425),
Myocardium (3), Skin (1), and Spleen (1). Fat has one unlocalized fallback. Counts are
source-dependent mapping facts, not physiological validation.

The held canonical lymph graph contains 996 nodes and 1117 edges. Its retained
source graph is a published computational topology and is undirected. Its presence
does not establish correspondence between native owners and graph nodes, regional
fluid/protein inventories, flow directions, radii, nodal pressures or regional
lymph transport. The report says these mappings and solutions are absent even when
the display graph exists. Without a supplied graph, its source-topology status is
`absent_source_graph`; no graph is fabricated.

## Bounded verification

Run `.venv/bin/python scripts/verify_lymph_coverage.py`. Six small tests reuse the
existing `verify_body_exchange.snapshot` fixture and check conservative allocations,
no second Lymph store, input immutability, unit/source identity, unknown observations,
unlocalized fallback, source-topology distinction, serial-flow reporting, invalid
owners/quantities, and geometry-receipt rejection. No native simulation is launched.
A workspace mapping-only construction also verified the held support/source graph
receipts and the counts above. This checks mapping integrity, not physiological
calibration or native integrated mass closure.
