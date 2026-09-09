# Skin rendering follows the contact material support

`SegmentSurfaceBinding` exposes the same hard material attachment used by the
world's sampled skin contacts. The whole skin no longer needs to inherit one
segment transform merely because it is one source mesh.

For every **canonical rest vertex**, choose the nearest named bone group's
axis-aligned rest envelope by Euclidean distance; ties choose the lexicographically
first native segment name. The first retained canonical bone in that group owns
the attachment. Subtract its reference centroid, rotate the resulting offset by
that bone's current rotation, then add its current centroid. Apply the common
canonical-to-world/view transform after this operation, once.

The mechanical frame carries frozen `surface_binding` metadata (including source
hashes, reference envelopes/centroids, segment and bone IDs, a stable
`binding_identity`, and `surface_entity_ids`) plus small per-frame
`surface_transforms` containing the 22 current segment centroid/rotation pairs.
The binding is fixed in rest coordinates and does not get recalculated from a
moving posed skin. There is no display smoothing or independent pose generator.
Respiration's display-only displacement cannot be added while retaining exact
agreement with the physical skin contact support.

Other extensive attached surfaces can use this same declared prior. A garment
with its own simulated positions must use those positions instead. At a
mixed-owner rendered triangle, nearest bound vertex selection is an approximate
force-picking rule; sending the selected bone ID preserves the intended segment
and material-point attachment. The frame adds neither forces nor body mass.

This is **not deformable skin FEM**. Hard assignment can expose seams at joints;
canonical anatomical registration is still approximate. The contract fixes the
render/contact disagreement without pretending to solve tissue continuity.

Verification: `PYTHONPATH=. .venv/bin/python scripts/verify_surface_binding.py`.
The retained receipt `data/derived/surface-binding-verification.json` checks all
102,467 full-resolution skin vertex owners against the original registration
ranking, then compares 1,193 contact samples in actual recorded native poses at
20 and 80 ms. Both canonical and world coordinates are bit-for-bit identical
(maximum error 0 m). The fixture `data/derived/surface-binding-fixture.json`
contains frozen metadata and 24 posed vertices for cross-language renderer tests.
These two short poses establish attachment parity, not long-duration stability or
anatomical validation.

## Continuous shared candidate

The hard attachment is retained as a comparison, **not a satisfactory skin
embedding**: an actual native pose at 80 ms stretches a 0.779 mm edge to 325 mm
where adjacent vertices chose hand and pelvis supports. The renderer parity test
correctly exposed this physical prior defect rather than concealing it.

`ContinuousSurfaceBinding` replaces the hard assignment with weights derived
from the retained skin topology. Seeds choose the nearest retained bone mesh
vertex (an explicit vertex-distance approximation to bone surface distance).
A screened weighted graph Laplacian diffuses those seed weights along skin
triangle edges, using an uncalibrated 80 mm length and 1 mm minimum edge length.
Exactly coincident vertices share supports; no distance-based weld joins nearby
body parts. The source has 100 raw connected components and 776 exact duplicate
vertices; the exact weld yields one connected graph with 101,691 vertices.
Negative solver roundoff is removed, tails below 1e-8 are discarded, and weights
are quantized to float32 then normalized once. The frozen sidecar is consumed by
both physics and display; poses never alter its weights.

The position of each skin point is the weighted sum of its segment-specific
posed material stations. A force **F** at that point scatters **wᵢ F** to each
segment at its own station **xᵢ**. This is the transpose of the velocity map. It
preserves total force, total moment about every origin, and virtual power;
applying those distributed forces at the blended point instead would generally
violate the rotational work relationship.

The sidecar `data/derived/canonical/continuous_surface_binding.json.gz` contains
source hashes, algorithm parameters, reference vertices, segment IDs and weights.
`asset_records()` supplies exact immutable JSON bytes keyed by SHA256;
`manifest()` declares `weights_url` and `weights_sha256`. Per-frame transforms
remain the small 22-segment pose set. Other attached surfaces transfer weights
from the nearest canonical rest skin vertex, breaking exact-distance ties by
source vertex index. Independently simulated cloth overrides this attachment.

Reproduce the artifact with:

```
PYTHONPATH=. .venv/bin/python scripts/build_continuous_surface_binding.py \
  --registration <retained-runtime>/mechanics/registration.json
PYTHONPATH=. .venv/bin/python scripts/verify_continuous_surface_binding.py
```

The receipt `data/derived/continuous-surface-binding-verification.json` measures
maximum edge stretch at 80 ms falling from 417.53 to 3.485 and the 99th percentile
from 12.82 to 1.474. Force, moment and virtual-power residuals are below 2e-15;
a shared rigid motion is reproduced within 3e-16 m. These are geometric and
coupling checks. The minimum edge ratio is still 0.0375 in that pose: substantial
local compression remains, and the method neither preserves tissue volume nor
prevents self-intersection. Native body/world regression is required before
promoting the changed contact embedding. This is not skin FEM or a calibrated
material law.
