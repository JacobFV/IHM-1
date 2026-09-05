# Published lymphatic structural network

The local model preserves **996 vertices, 1,117 undirected structural edges and 272 lymph-node flags** from the supplementary files of Savinkov et al. (2020), [Graph Theory for Modeling and Analysis of the Human Lymphatic System](https://doi.org/10.3390/math8122236). Coordinates, node flags, endpoint indices and source edge lengths are copied without downsampling. The app fragment uses all 1,117 source edges as line segments in a separate model/frame.

This is **published model geometry, not measured human anatomy**. The antecedent paper, Tretyakova et al. (2017), [Developing Computational Geometry and Network Graph Models of Human Lymphatic System](https://doi.org/10.3390/computation6010001), explicitly explains that the PlasticBoy basis is an anatomical estimate rather than a measured dataset. The local acquisition contains the authors' graph, not the original PlasticBoy polygon assets.

## Reproduce

```sh
python3 scripts/collect_lymph_network.py
python3 scripts/build_lymph_network.py
python3 scripts/verify_lymph_network.py
# Explicit app integration, after the base app manifest exists:
python3 scripts/build_lymph_network.py --append
```

`build()` only writes `data/derived/lymphatic`; it leaves the shared app manifest alone. Integration is also available as `append_to_manifest(manifest_path)` from the builder. The fragment's geometry `kind` is `lines`, with a `line_segments` primitive payload compatible with the app. The append hook copies geometry into the app geometry directory. It does not merge this graph with another atlas or claim cross-atlas registration.

## Files and provenance

The collector pins SHA-256 identities for both papers, publisher XML and [the publisher supplement ZIP](https://mdpi-res.com/d_attachment/mathematics/mathematics-08-02236/article_deploy/mathematics-08-02236-s001.zip). ZIP hash:

`7c7b05eb6ae5dd7682ed0ed1fc89e8f08eb0c33d8bc3ac8711f97a5c07a0d390`

The publisher XML links this ZIP as supplementary material. All five files are retained, including the rule-based graph generator and supplementary metrics figure; **the generator is not executed**. The geometry uses only `Data-based HLS graph/graph_vertices.txt` and `graph_edges.txt`.

- `data/raw/lymphatic/provenance.json`: URLs, hashes, byte counts, archive-member paths, version evidence, attribution, license caveat and units evidence.
- `data/derived/lymphatic/graph.json`: lossless parsed source vertices, topology and source lengths, plus explicitly derived chord lengths and unknown physical quantities.
- `data/derived/lymphatic/manifest_fragment.json`: isolated model/structure descriptions with provenance and display transform.
- `data/derived/lymphatic/geometry/published-lymphatic-network.json.gz`: all source edge endpoint pairs in display meters, with edge IDs.
- `data/derived/lymphatic/verification.json`: independent comparison against archived source rows, topology checks, rendered coordinate checks and hashes.

Raw/derived data follow the repository's existing ignored-data convention and are reproducible with the checked-in collector/build scripts.

The publisher article/XML declares CC BY 4.0. The supplement has no separate license file; provenance records this distinction rather than inventing an independent license grant. Credit the 2020 authors, the 2017 antecedent and PlasticBoy as the basis. The article license is not asserted for PlasticBoy's original commercial assets.

The requested [Wolfram resource](https://datarepository.wolframcloud.com/resources/63e705c2-bcc6-4458-959d-bd5f8baafd30/) was discoverable in search but returned HTTP 404; its cloud endpoint showed a scheduled-upgrade page. The acquired dataset therefore comes directly from the primary publisher supplement, not an asserted Wolfram export. Wolfram's published 996/1,117 counts match this source graph.

## Graph schema and interpretation

Vertex IDs are the **one-based source data-row order**, excluding comments. `nodes` contain `id`, `position_mm`, `is_lymph_node` and `anatomical_name: null`; the source supplies no per-node names. `edges` contain `id`, `source_from`, `source_to`, `source_length_mm`, independently computed `endpoint_chord_length_mm`, and null `radius_mm`, `flow_ml_s` and `flow_direction`.

`directed: false` is intentional. The source headers say `From To Length(mm)`, but the paper calls the original graph undirected and uses a separate Poiseuille-based simulation to propose orientation. Endpoint ordering alone is not physiological direction evidence; no such flow solve is reproduced here.

The paper's section 2.2 says coordinates were scaled to a baseline 1,750 mm body height. It also documents author-added connections/output vertices (987–993, 993–995, 994–995, 995–996) and removal of some degree-two non-lymph-node vertices. These are model construction choices already present in the source; no local missing vessel is invented. The 272 source flags are not interpreted as 272 independently observed lymph nodes in a measured person.

Display applies the documented rigid rotation `[x,y,z] → [x,z,-y]`, converts mm to m, and recenters the bounding box. No fitting, stretching, smoothing, nearest-neighbor connections, radii, tissue surfaces or registration are introduced. Displayed edges are straight endpoint chords. Source length rounding differs from independently computed chords by at most 0.000506 mm. The graph is connected and has no duplicate undirected edges.

Unknown radius, pressure, flow, valves and anatomy beyond this source remain unknown. Passing extraction verification establishes data/geometry fidelity, not physiological or subject validation.
