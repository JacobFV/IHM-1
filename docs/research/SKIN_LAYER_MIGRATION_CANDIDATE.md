# Isolated skin-layer migration candidate

`prepare_skin_layer_migration.py` prepares copies of current `anatomy.json` and
`mechanics.json` in a fresh directory. It never installs them or runs either
full builder. Geometry references and their source bytes remain in the existing
repository; only the exact skin geometry is read. The native 22-body model,
its mass parameters, muscle paths and native checkpoints are not modified.

Run:

```sh
.venv/bin/python scripts/verify_skin_layer_migration.py
.venv/bin/python scripts/prepare_skin_layer_migration.py --output /tmp/new-skin-candidate
```

The candidate requires matching current anatomy receipts in retained mechanics
and the component evidence. It validates entity identity, original mechanics
normalization, inertia and every link damping law before changing anything.
Three layer volumes become reviewed exterior-proxy area times their explicit
shell thickness. Depth intervals, thickness and geometry references survive.
All proxy masses and inertias are then recomputed under the existing builder's
fixed total-mass allocation, along with every link's reduced-mass damping.
Support stiffness, topology, registration and actuator parameters are retained.
The mechanics anatomy hash refers to the candidate anatomy bytes, under the
logical canonical path; it does not claim installation at that path.

Output files are `anatomy.json`, `mechanics.json`, and `manifest.json`. The manifest
records input/output hashes, helper/builder/migration source hashes, changed
volumes and old/new normalization. The existing evidence receipt is retained as
historical input to the migration. It must not be relabeled as newly generated
against candidate anatomy: doing so also changes its hash embedded in anatomy
and creates a circular identity dependency. Later territory regeneration should
create a separate receipt against candidate anatomy, leaving this historical
component-mask evidence intact.

The bounded actual-input candidate at `/tmp/ihm-skin-layer-candidate-20260906`
computed 1.7804602548390722 m² exterior-proxy support versus 3.502598974493317 m²
raw shell area. Total shell thickness remains 0.0066 m. Unscaled proxy mass
changes from 64.4865523954032 to 53.12043684568518 kg; the uniform normalization
changes from 1.195763985446006 to 1.4516201574924263. Thus **all non-carrier proxy
masses change**, including unchanged non-skin geometry. Their final sum plus
six 1 mg carriers remains the existing 77.1107029 kg constraint. This is a
correction to the existing overlapping atlas allocation, not measured weight
loss, native body mass loss or a physiological material transfer.

## Dependency and acceptance order

| Consumer/source | Affected output or identity | Required next action before adoption |
| --- | --- | --- |
| `build_body_mechanics.py` | `canonical/mechanics.json`, `source_files[canonical/anatomy.json]`, mass allocation, inertias, link damping | Candidate regenerates these laws from retained topology; compare preserved non-mass fields and normalized totals. |
| `build_supine_surface_contact.py` | Fresh support `materialization.json` and contact material thickness | Use explicit carried `shell.thickness_m`; volume divided by raw two-shell area would incorrectly reduce 6.6 mm to 3.355 mm. Support owner manages this migration. |
| `ihm/assembly/articulated.py:115`, `embodied.py:61` | Per-run copied `canonical_mechanics.json`, canonical reference mass/registration hashes | New owners must freeze the new mechanics identity. Existing owners and checkpoints remain tied to old frozen input. |
| `materialize_engineered_skin_territories.py:129` | Fresh territory `materialization.json` contains `anatomy_sha256` | Re-materialize/revalidate a new artifact epoch; preserve old component-mask source evidence used by candidate. |
| `ihm/assembly/respiration.py`, `build_body_respiration.py` | `canonical/respiration.json` contains anatomy receipt | Regenerate in candidate environment; systemic projection rejects stale `anatomy_sha256` at `systemic_projection.py:108`. |
| `build_hair_strands.py:55` | Hair `manifest_fragment.json` anatomy/source hashes | Geometry/root positions are unchanged; perform dependency-equivalence validation or fresh generation. Do not silently rewrite historical receipt hashes. |
| `build_body_material_domains.py:21` | Fresh material-domain manifests/ownership, anatomy and mechanics inputs | Revalidate allocation-dependent metadata; local geometric partitions may remain equal but changed input hashes require a new receipt. |
| `body_runtime.py:27`, `regional_touch.py:34` | Body state registry identity; touch source hashes | New registry/touch outputs must name new anatomy epoch. |
| `body.py:72`, `body.py:107` | Canonical body assets/trajectory consume anatomy and related profile | Reassemble and verify dependent outputs before publishing a coherent canonical set. Existing research trajectories remain historical. |

These are traced direct consumers, not an exhaustive automatic dependency
resolver. Whole-file anatomy hashes can invalidate unrelated geometry-only
research receipts; equality of geometry does not justify asserting that those
runs used the new anatomy. No shared canonical regeneration is authorized by
running this tool.

Bounded acceptance completed: fixture rejection of stale source receipts,
mass normalization and damping; exact retained skin evidence evaluation; actual
retained metadata baseline law checks; candidate mass/inertia/damping recomputation.
No complete geometry scan, native acceptance, shared canonical overwrite or
supine/native mass change was performed. Adoption still requires the downstream
owner checks above and a coordinated, reviewable canonical publication.
