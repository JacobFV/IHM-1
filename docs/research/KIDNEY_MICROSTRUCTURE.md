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

The slab contains 114,093,568 retained TIFF bytes, with array shape 32 × 1303 × 912, uint16 images and masks containing 0/255. The slab manifest records each archive member, source DOI, original size, ZIP CRC32 and SHA-256. Member CRCs are verified; the whole archive MD5 is explicitly not verified because the complete archive was not downloaded. Original images and gold labels remain separate. Predictions from winning models are not substituted for gold labels. No segmentation-derived graph has yet been created from these images.

## Verified physical spacing and donor mapping

[Table 1 in the 2026 paper](https://www.nature.com/articles/s41467-026-74050-8/tables/1) explicitly identifies the donor-1 densely labeled whole-kidney dataset as the **LADAF-2021-17 right kidney, male, age 63, BM05, 50 µm voxels**. Its high-resolution VOI is a separate 5.2 µm row. Supplementary Table 1 maps donor 1's whole right kidney to `10.15151/ESRF-DC-1773966439`, the same scan cited by the 2025 arterial-network paper. Thus the slab and arterial example originate from the same donor and organ; they are not independent donors. This establishes nominal isotropic voxel spacing `[50, 50, 50]` µm, without guessing from the approximate overview description. Image/slab-to-example coordinate registration and their exact spatial overlap remain unknown.

The raw 2026 article, full Table 1 HTML, and supplement PDF are retained under `segmentation_19685382/metadata/`. Their evidence receipt records URLs, byte sizes and SHA-256 hashes. Table 1's retained hash is `7f523bc8d48840c96b85daa3d46312c344c1ae82a68f6c93f7ff71f8ddd1bca2`; the supplement hash is `d4286901266a9b68101dec6ce2407559da8e57265fbb1bcde5d2248b2214451f`. Verification rechecks those receipt digests and the exact nine-column Table 1 row. The publisher lists 2279 slices, whereas the archive contains 2280 nonempty image TIFFs indexed 0000–2279. This discrepancy is retained explicitly; no slice is silently removed or renumbered.

At the published spacing, the slab contains **490,087 labeled vessel voxels, equivalent to 61.260875 mm³ of segmented volume**. Its rectangular image box is 4753.344 mm³, giving a labeled fraction of 0.01288795 of that box, including background. This is not vessel volume fraction within kidney tissue because no kidney-tissue denominator mask was acquired. The labels intersect both slab cuts: 17,420 foreground voxels at the first z-face and 12,084 at the last. The measurement is therefore boundary-censored. It describes segmented volume inside the acquired slab only; it does not estimate whole vessels, in-vivo lumen volume, capillary completeness, organ-wide density, or population variation. Nominal voxel calibration does not validate radii or Murray-law priors.

## Commands and verification

```bash
.venv/bin/python scripts/collect_kidney_microstructure.py
.venv/bin/python scripts/collect_kidney_microstructure.py --offline --slab-slices 32 --slab-start 1100
.venv/bin/python scripts/collect_kidney_microstructure.py --offline --refresh-slab-metadata
.venv/bin/python scripts/verify_kidney_microstructure.py
```

Graph verification checks source digests, finite values, declared array counts, valid node references, polyline endpoints, connected-tree topology, supplied-order partition counts, and an analytic 3-4-5 bent-polyline example that detects unit/path/chord mistakes. Invalid connectivity is rejected. Slab verification additionally uses Pillow (available in this workspace) to decode every TIFF, check shape, data type, mask values, consecutive slice indices and image/mask pairing, and recheck both SHA-256 and CRC32. `--offline` applies to the graph acquisition only; requesting a slab may access the network for missing members. Existing slab members are checked by ZIP CRC before reuse. Larger slabs can be retained with `--slab-slices`; they are never called whole-kidney coverage.
