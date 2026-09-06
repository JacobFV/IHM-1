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

## Persisted candidate-root acceptance (2026-09-06)

The accepted bounded staging epoch is
`data/derived/skin-layer-epoch-20260906-v2`. It contains an isolated `root/`,
`registration.json`, `differences.json`, historical `migration_manifest.json`,
and `acceptance.json` with SHA-256 hashes for 28 output files. Derived payloads
remain local retained artifacts; acceptance/difference receipts are committed.
The earlier v1 stage remains historical and is not the selected candidate.

```sh
.venv/bin/python scripts/stage_skin_layer_candidate.py \
  --candidate /tmp/ihm-skin-layer-candidate-final-20260906 \
  --output data/derived/a-new-skin-epoch
.venv/bin/python scripts/verify_skin_candidate_stage.py \
  data/derived/skin-layer-epoch-20260906-v2
```

The stager copies bounded canonical metadata, the profile's exact patient file,
and the body manifest's 14 runtime source files into its private root. It
regenerates respiration from staged anatomy and calls the actual `body.build`
and `CanonicalBody.from_workspace` interfaces. All 2,408 entities load with
matching staged source/executor/patient receipts. This is loader acceptance,
not a running physiology/mechanics/body-runtime session.

Registration is recomputed twice with `CanonicalRegistration`: original and
candidate mechanics use the same retained
`data/derived/supine-support-5ma720yd/initial_native.json`. The resulting full
registration manifests are exactly equal. No native executable is instantiated;
this registration comparison provides no new native equilibrium evidence.

Exact old/new differences, enumerated in `differences.json`:

- Anatomy: four changed entity rows, consisting of the parent support receipt
  and three layer support/volume/method updates. Assumption ledger also changes.
  Every entity's geometry reference, bounds, centroid and principal axis is equal.
- Mechanics: 2,399 rows change only proxy mass/inertia; three layer rows also
  change material volume, volume basis, shell and support; parent skin gains
  support only. All 3,341 support-link damping values change coherently.
- Respiration: only `anatomy_sha256` changes; reconstructed physical parameters
  and bindings are equal to the retained artifact.
- Body manifest: only source and executor receipts change. Registration is equal.
  Native mass migration is explicitly false.

**Source-mask reuse can use geometry identity.** `physical_skin_support` already
checks exact skin bytes and the matching source geometry receipt independently
of the evidence's whole-anatomy epoch. The stager additionally proves equality
of every source-mask input entity's geometry reference, bounds and centroid,
plus equality of all candidate anatomical geometry metadata. It retains the
historical evidence bytes and old `anatomy_sha256`; an explicit new equivalence
receipt establishes the restricted reuse. This does not relabel the original
territory materialization as a new run or establish anatomical drainage.

Direct-consumer evaluation is recorded under `acceptance`:

| Consumer | Bounded result | Remaining publication requirement |
| --- | --- | --- |
| Body manifest/loader | Actual build and loader pass, all hashes checked | Start a new owner only after coordinated publication. |
| Articulated registration | Actual old/new reconstruction exactly equal | Keep existing native owners pinned; create new provenance for future owners. |
| Respiration/systemic projection input | Actual respiration rebuilt; only anatomy receipt differs | No systemic/native trajectory was rerun. |
| Engineered territory component mask | Exact skin bytes and referenced geometry metadata agree | Keep source artifact historical; separate equivalence receipt is sufficient for mask support, not a new territory claim. |
| Supine contact | New shell metadata is consumable by the migrated builder; source geometry equal | Support owner must build a separate contact/statics epoch; old accepted equilibrium remains pinned. |
| Hair fragment | Referenced skin geometry and all anatomical geometric metadata equal | Existing full-anatomy receipt remains historical; hair owner needs explicit equivalence or fresh fragment generation. |
| Material-domain ownership | Geometry equal, canonical allocation changed | Dedicated domain metadata regeneration/acceptance still required; no hash-only relabeling. |
| Body-runtime registry/local touch | Candidate body assets load; anatomy epoch necessarily changes | Registry/touch outputs must be created under the new owner; no runtime stepping performed. |
| Viewer asset serving | Source references preserved, only skin geometry physically copied | Stage is not a full viewer root; other geometry remains in original source storage. Explicit source-root asset serving or complete asset materialization is required before viewer adoption. |

No heavy build is needed to repeat the completed checks. Full material-domain,
contact or native acceptance must be separately scoped; the candidate is not
ready for wholesale canonical promotion until those required consumers are
resolved. The root stage has no symlinks into writable live canonical assets.

## Remaining consumer candidates and complete geometry snapshot

The next selected metadata epoch is
`data/derived/skin-consumers-20260906-v2`; v1 is historical. Run
`verify_skin_consumer_candidates.py data/derived/skin-consumers-20260906-v2`
to recheck retained outputs. `complete_skin_candidate_consumers.py` reproduces
this epoch from the staged root into a fresh output directory.

Both 4 mm and 2 mm retained pelvic domains were rebound without voxelization.
Their original geometry source hashes, tetrahedra, vertices, owner indices,
constitutive arrays, boundary triangles and fixed nodes are checked and retained.
Per-owner tetra volume is checked against the parent materialization. Candidate
allocation updates registry masses, per-owner inertial densities and per-element
density arrays; `DynamicTetrahedra` reconstructs nodal mass and the explicit
stability limit. Exclusive domain mass is checked through `MaterialOwnership`.
No domain is activated, and no native mass is debited.

| Domain | Tetrahedra | Old/new allocated kg | Old/new explicit limit s |
| --- | ---: | ---: | ---: |
| Pelvis 4 mm | 5,796 | 0.0762281941 / 0.0925386485 | 0.0000892907 / 0.0000983807 |
| Pelvis 2 mm | 46,890 | 0.0762281941 / 0.0925386485 | 0.0000459624 / 0.0000506415 |

All three owner densities scale by 1.2139687891260489, matching the existing
canonical reallocation. Geometry and elastic material priors do not change.
Parent manifests remain historical; new manifests bind the candidate anatomy
and mechanics and their exact new arrays/registry files.

`hair-equivalence.json` records the retained elastic_v3 fragment, exact skin
geometry, three population/render asset hashes and unchanged source evidence.
It does not regenerate roots, claim a new historical hair build, migrate
residual native inventory, or certify evolving hair runtime code. Explicit
old/new anatomy hashes establish only the geometric dependency equivalence.

The separate `support/` candidate inherits the exact current
`lumbar-supine-reference-1qex2x9i/contact` quadrature and native input bytes.
Only layer thickness provenance and candidate source receipts change. The
original omitted index 15948 / source face 165251 remains omitted; no deeper
ray hit replaces it. Material thickness remains 0.0066 m. `accepted_support`
and `native_integration` remain false: old equilibrium evidence cannot be
relabelled as a solve under a new canonical source epoch. A new native solve
is not required to prove byte equality, but any claim of a new accepted running
owner belongs to its separate native acceptance workflow.

The authorized geometry snapshot is
`data/derived/skin-viewer-epoch-20260906-v1/root`. Its `acceptance.json` records
2,405 unique geometry files totaling 192,730,414 bytes. A low-priority sequential
pass hashes each source, makes an independent reflink-or-copy, checks destination
hash, and rechecks source hash, inode, size and modification time. Every
source modification time remains unchanged; no writable hardlink aliases exist.
The plan in consumer v1 and selected consumer v2 is byte-identical and shares
the exact plan hash recorded by this geometry snapshot.

This resolves the incomplete canonical geometry serving inventory: the private
root now contains every anatomical geometry reference, and the real body loader
accepts all 2,408 entities. The original 193 MB dataset was not mesh-decoded,
remeshed or sent to a native process. `materialize_skin_viewer_root.py` rejects
existing destinations and enforces the explicit byte budget before copying.
No browser, app server or viewer deployment was performed. Root promotion still
needs its coordinated application/owner switch; all live canonical and frozen
native references remain unchanged.
