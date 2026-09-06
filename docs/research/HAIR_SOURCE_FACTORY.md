# Exact hair source factory preparation

`data/research/hair_source_factory.json` is a source-only preparation receipt.
It does not enable a native owner or change source geometry, strand counts,
native mass, or viewer dynamics.

## Exact retained correspondence

The factory joins `attachment.sample_ids` to the original population `ids` and
its retained `face_index`. It then selects those exact indexed source faces,
retains their global source node/face indices, and verifies the attachment
triangle coordinates, barycentric weights, roots, face-normal directions,
strand radii, lengths and complete initial centerlines. No nearest-face query,
resampling, root snapping or guide-population weighting is used.

| Retained group | Guides | Exact faces / nodes | Guide nodal mass | Full population prior |
| --- | ---: | ---: | ---: | ---: |
| Scalp | 64 | 64 / 189 | 15.6347 mg | 92,934 strands; 22.8101 g |
| Body | 32 | 32 / 95 | 0.0250671 mg | 201,528 strands; 0.155011 g |

Guide moments use the exact lumped-node masses of `ElasticHairState`. Population
moments integrate independent straight circular cylinders using retained roots,
normals, radii, lengths and the stated 1,312 kg/m³ density reference. The latter
are full-population material/geometric priors, not measured subject masses and
not weights assigned to a guide. The artifact retains first moments and full
second-moment/inertia tensors in the canonical reference frame.

Inputs are captured and hashed once, then parsed from the same immutable bytes.
Large sibling render arrays and complete source meshes are not decoded: only
small hair members and requested indexed source coordinates are materialized.
The body population digest is pinned by the prior elastic manifest. That
manifest omitted a scalp-population digest; this receipt pins the current
retained scalp bytes and records that provenance limitation explicitly.

## Reusable source skin owner protocol

`ihm/assembly/source_skin.py` provides `extract_source_patch` and the immutable
`SourceSkinPatch` interface. Extraction establishes geometry correspondence.
`SourceSkinPatch.bind` additionally requires an explicit binding containing:

- Matching source entity and geometry SHA-256.
- Matching frozen registration SHA-256 and a stated assignment basis.
- Exactly one supplied registered owner for each retained global source node.

Missing, duplicate, unknown or mismatched owners are rejected. No owner is
inferred by this module. Existing nearest-bone-envelope ownership, if supplied
by a caller, must remain labeled as a fixed registration prior rather than an
anatomically calibrated compartment assignment. A follicle spanning different
owners is rejected because its clamp orientation derivative is unresolved.

`prepare_hair_factory` validates the complete retained guide definition and
returns the local follicle faces, global source correspondence, barycentric
data, guide owners and dependencies. It always returns `native_enabled: false`.
The anatomical audit currently has no supplied exact node-owner binding, so it
records validated geometry correspondence with `owner_binding: null`.

The extracted faces support registration only. They are sparse follicle faces,
not a complete neighborhood for body/garment collision coverage.

## Native inertia partition proposal

Canonical metadata gives the head-hair proxy 0.472760339 kg and pubic-hair proxy
0.018625841 kg. These proxies are **not separately integrated native masses**.
In particular, the retained body-hair guide regions are not a justified mapping
onto the pubic proxy. Neither proxy is an automatic debit budget.

The source model totals 85.269848542 kg. The actual retained cutaneous factory's
`assembled_model.osim` totals 77.6122029 kg, matching its mechanical identity.
The audit stores all 22 actual local mass/COM/full inertia records plus their
native reference transforms and canonical-reference-to-body-local embeddings.
These are read from retained production artifacts, not reconstructed by blindly
scaling guessed values. No running native process is launched.

A conservative **inferred partition**, once source ownership is supplied, is:

1. Transform each selected independent guide's actual nodal mass moments into
   its exact native owner's reference frame. `transform_mass_properties`
   accepts proper rigid transforms only.
2. Subtract only those selected moments from that actual owner's mass, first
   moment and second moment. Leave unrepresented population in the residual
   native budget. Use full-population priors to assess plausibility without
   treating them as calibrated partitions or bundle multipliers.
3. Recompute the residual COM and full inertia tensor. Reject exhausted mass or
   a nonpositive centered second-moment tensor. Check that residual native owner
   plus explicit guides recover the original mass and moments.
4. Freeze the allocation basis and source dependencies. Install nothing until
   the native mass/COM/inertia path and one shared-clock coupled interval are
   verified. Hair and garment feedback must not advance the same body twice.

`partition_native_inertia` implements only this algebraic admissibility check.
It requires explicit same-owner/common-frame inputs and an allocation basis,
returns `native_enabled: false`, and mutates neither input nor native state.
No anatomical debit has been selected or applied.

## Bounded checks

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_hair_source_factory.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/audit_hair_source_factory.py
```

The tiny verifier checks exact joins, immutable retained bytes, source/hash/
owner rejection, cross-owner clamp rejection, coordinate mismatch rejection,
rigid moment transforms, mass/moment recombination, invalid residual tensors,
and persistent non-enablement. The anatomical audit verifies correspondence and
retains the source/inertia receipt. Neither establishes native dynamics,
complete contact coverage or authoritative viewer synchronization.

The bounded anatomical audit completed in 0.68 s at 160,308 KiB peak RSS, with
one BLAS thread. All 17 retained dependency hashes were checked after generation.
The small synthetic verifier completed in 0.12 s at 57,116 KiB peak RSS.
