# Kidney microstructure acquisition

The final arterial graph from [Rahmani et al., npj Imaging 2025](https://doi.org/10.1038/s44303-025-00090-2) is **not acquired**. Supplementary Note 7 names DOI `10.5281/zenodo.15275990` and the cc1, cc4 and cc9 Amira graphs. On 2026-09-05 that record's API returned HTTP 404, “The persistent identifier is not registered.” The older [7708966 record](https://zenodo.org/records/7708966) exposes metadata but restricts its files. The collector retains the real response bodies and the publisher's supplement, with sizes and SHA-256 digests in `data/raw/microstructure/kidney/acquisition.json`.

## Actual graph bytes

The paper's [author-maintained analysis repository](https://github.com/HiPCTProject/Skeleton_analysis/tree/e35024b33018e010a057a3e69147459e7cadcc4e) contains `Initial_Ordering/Test.am` and `Test2.am`. Their embedded original path identifies `LADAF-2021-17-right/Artery1_only.Spatial-Graph_Reduced.am`. These are a reduced example from the donor, with **148 nodes, 147 edges and 19,924 polyline points**, not the paper's 10,190-vessel network. Test2 supplies Strahler and topological labels. The collector verifies that its coordinates, connectivity, point offsets and thickness values exactly match Test.am before accepting those labels. Its damaged coordinate-unit text is not independently interpreted.

`example_graph.npz` contains `nodes_mm[N,3]`, `edges[E,2]` (zero-based), `point_offsets[E+1]`, `points_mm[P,3]`, `thickness_native[P]`, `strahler[E]`, `topo[E]` and `source_transform[4,4]`. Coordinates are converted from explicitly declared micrometers to millimeters and remain in the original local coordinate frame. The source matrix is retained without applying it or registering the graph into the body.

The point attribute is named `thickness`; its unit metadata is unset. Upstream code itself includes a `radius?` comment. We preserve those numbers and **do not infer radius or diameter**. Radius distributions and Murray residual statistics are null. No cubic Murray constraint is imposed. Repository-wide or graph-specific licensing was not found; do not assume the article's CC-BY license grants reuse of these GitHub assets. Correction status is also unknown for the example. Consequently it is marked `measured_other_donor`, with `usable_for_population_priors=false`.

## Reproducible geometry measurements

`statistics.json` reports global and supplied-order-conditioned segment length (sum of Euclidean polyline increments), tortuosity (path length / endpoint chord), and mean native thickness. The latter remains an uninterpreted quantity. It also reports the node-degree histogram: 75 terminal nodes and 73 degree-three nodes. All 148 nodes belong to a single tree. Mean segment length is 1.103236 mm; mean path/chord tortuosity is 1.416294. These describe this reduced, potentially unsmoothed example only, and do not establish organ-wide physiology or population variance.

The published donor was a 63-year-old male who died of pancreatic cancer. Formalin embalming, formaldehyde postfixation and ethanol dehydration preceded ex vivo imaging. Acquisition used 25, 6.5 and 2.6 µm voxels; reduced analysis scales were 50, 13 and 5.2 µm. Published collapse corrections include manually checked perimeter-equivalent radii and local outlier replacement. None are newly applied here. The 97% detection statement concerns vessels with radius >50 µm in the paper's overview validation, not this example's completeness. Neither the published arterial tree nor this example is a capillary-complete, glomerular or microvenous graph.

## Paired raw image and label slab

The viable next source is [Zenodo 19685382](https://zenodo.org/records/19685382), associated with the [2026 kidney vessel segmentation paper](https://doi.org/10.1038/s41467-026-74050-8). Its public CC-BY-4.0 archive is 51,177,854,763 bytes and contains 33,706 ZIP members. The server supports byte ranges. The collector retrieves the ZIP64 central directory and a bounded slab of **32 consecutive TIFF image slices and matching dense vessel masks** from `competition-data/train/kidney_1_dense`, starting at sorted index 1100. This is a real contiguous 3D slab, not isolated illustrative screenshots; retained files preserve the publisher's TIFF bytes.

The slab contains 114,093,568 retained TIFF bytes, with array shape 32 × 1303 × 912, uint16 images and masks containing 0/255. The slab manifest records each archive member, source DOI, original size, ZIP CRC32 and SHA-256. Member CRCs are verified; the whole archive MD5 is explicitly not verified because the complete archive was not downloaded. Original images and gold labels remain separate. Predictions from winning models are not substituted for gold labels. This dataset is not silently identified with the 2025 donor, and physical voxel spacing must be verified before deriving morphometric priors. No segmentation-derived graph has yet been created from these images.

## Commands and verification

```bash
.venv/bin/python scripts/collect_kidney_microstructure.py
.venv/bin/python scripts/collect_kidney_microstructure.py --offline --slab-slices 32 --slab-start 1100
.venv/bin/python scripts/verify_kidney_microstructure.py
```

Graph verification checks source digests, finite values, declared array counts, valid node references, polyline endpoints, connected-tree topology, supplied-order partition counts, and an analytic 3-4-5 bent-polyline example that detects unit/path/chord mistakes. Invalid connectivity is rejected. Slab verification additionally uses Pillow (available in this workspace) to decode every TIFF, check shape, data type, mask values, consecutive slice indices and image/mask pairing, and recheck both SHA-256 and CRC32. `--offline` applies to the graph acquisition only; requesting a slab may access the network for missing members. Existing slab members are checked by ZIP CRC before reuse. Larger slabs can be retained with `--slab-slices`; they are never called whole-kidney coverage.
