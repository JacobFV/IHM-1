# Extended source anatomy: Z-Anatomy

This collection adds genuine authored geometry from the [official Z-Anatomy repository](https://github.com/Z-Anatomy/Models-of-human-anatomy/tree/702db27cf1395bf067158e86c7659b6e37a7be2c), pinned to revision `702db27cf1395bf067158e86c7659b6e37a7be2c`. The 86,734,957-byte source archive is retained, alongside its README, license, extracted 306,838,281-byte Blender file, immutable checksums, and download URLs. It partly derives from BodyParts3D. It is an authored atlas, **not an additional measured human, population sample, or independent validation**.

The app family is `z-anatomy`. Its geometry shares the source atlas frame internally; no registration to BodyParts3D or other families is assumed. Classification uses source collection membership plus documented name rules. Anatomical object counts are not counts of distinct whole organs or individual nodes: bilateral parts and multi-node group objects retain their original names.

## Fidelity and evaluation

`source_index.json` records each exported object's collections, source/evaluated matrices, base and evaluated triangle counts, bounds, geometry hashes, and authored modifier settings. `raw_object_inventory.json` inventories every Blender object, including excluded objects, labels, collections, base counts and modifiers. The unmodified archive remains the authoritative complete source.

Both original and evaluated geometry are retained in separate NPZ files:

- `base_geometry/`: original mesh datablocks, or tessellated source curves before object modifiers.
- `source_geometry/`: Blender dependency-graph evaluation using the source's existing viewport modifier settings. This includes authored Mirror, Solidify, Subdivision, and Geometry Nodes. Original subdivision levels and curve tessellation settings are not raised.
- Both preserve local vertices, triangulated local faces, world vertices, and world faces. The source/evaluated `matrix_world` maps local into atlas coordinates. Face winding is reversed only when the matrix has a negative determinant, with local face ordering retained separately.
- `geometry/`: app JSON geometry containing every evaluated triangle, **without decimation or upsampling**. Triangle indices remain exactly equal to the evaluated source arrays. Display positions are rounded to eight decimal places in normalized units and normals to six; the authoritative NPZ arrays are unrounded float64 exports of Blender coordinates.

The source scene declares `METRIC` and `scale_length=1`; those settings and the Blender/source-file versions are recorded. This is not independently measured physical calibration. For display only, a fixed Blender-Z-up rotation and one uniform normalization scale center this atlas and make its longest extent two units. The transformation is recorded and verified independently from the source world geometry. No per-organ rescaling or rearrangement is performed.

The evaluator runs with `--disable-autoexec`: embedded application Python does not execute. Blender's complete stdout/stderr is retained at `blender_evaluation.log`. The source produces an oesophagus/profile dependency-cycle warning. Its evaluated output is retained with that warning; it is not certified to be an unambiguous independent reconstruction. The log also reports a removed obsolete text-editor UI region and Blender shutdown memory bookkeeping. These are retained, not silently suppressed.

## Coverage and limitations

The verified export contains **2,581 anatomical surface objects, 8,111,194 evaluated triangles, and 4,091,749 evaluated vertices**. Original pre-modifier/tessellated base surfaces have 6,928,996 triangles. All 8,111,194 evaluated triangles reach the display output. The full inventory covers 7,184 source objects; 4,603 are excluded from display with reasons. `coverage.json` preserves these counts, and the verifier checks them against every source and display array.

| System | Surface objects |
|---|---:|
| Lymphatic | 163 (158 node-group objects, 2 thymic lobes, 2 tonsils, spleen) |
| Muscular | 495 |
| External body regions | 256 |
| Connective | 601 |
| Skeletal | 277 |
| Arterial / venous / cardiac | 425 / 230 / 17 |
| Digestive / respiratory | 48 / 39 |
| Endocrine / urinary / reproductive | 10 / 6 / 14 |

 The collection includes regional lymph-node groups across the head, neck, trunk, pelvis and limbs, spleen, thymic lobes, tonsils, muscles, connective structures, skeletal structures, arteries, veins, small visceral organs/ducts and external body regions.

No named thoracic duct, cisterna chyli, or lymphatic-vessel geometry exists in this pinned source inventory. A collection heading or anatomical definition is not counted as geometry. External body-region surfaces are useful external anatomy, **not** measured skin thickness, histological layers, follicles or sweat glands. Mesh names and left/right suffixes remain source-authored. No measurements, functional lymph flow, or physiological calibration are invented.

Annotations are excluded by their `.j`, `.g`, `.t`, or `.i` source suffix; text/cameras/lights and insertion-annotation-only or unlinked objects are excluded. Empty evaluated surfaces are explicitly reported. Source geometry and modifier-derived surfaces are distinguished, rather than treating authored procedural node instances as independent source measurements.

## Licensing and attribution

The [source README](https://github.com/Z-Anatomy/Models-of-human-anatomy/blob/702db27cf1395bf067158e86c7659b6e37a7be2c/Readme.md) declares **CC BY-SA 4.0** for the project. Derived exported atlas geometry remains CC BY-SA 4.0 with attribution; its license is recorded per structure and family. Attribute:

> Z-Anatomy — The libre 3D atlas of anatomy — Gauthier Kervyn and contributors. BodyParts3D — The Database Center for Life Science — CC-BY-SA 2.1 Japan (as credited in this Z-Anatomy source).

The README also credits an inner-ear adaptation from the University of Dundee under CC BY-NC-SA 4.0 and a kidney adaptation from Lissie Cowley under CC BY-NC 4.0. The blanket project license does not resolve these component restrictions. The complete raw archive retains those third-party assets and attribution notices; the display export conservatively excludes the nervous/sense-organ source collection and all kidney-named objects. This is a collection/name-based exclusion, not a forensic per-vertex authorship determination. The source's additional credits remain in the acquired README. No claim that all raw archive contents are unrestricted is made.

## Reproduction and integration

Run from the repository root:

```bash
python3 scripts/collect_extended_anatomy.py
uv pip install --python /usr/bin/python3 --target .cache/blender-python numpy==2.5.2
.venv/bin/python scripts/build_extended_anatomy.py --force-extract
.venv/bin/python scripts/verify_extended_anatomy.py
```

This environment uses system Blender 4.0.2 and its Python 3.12; the source file records Blender version 3.5.10. Dependency-graph mode is `VIEWPORT`. The subprocess defaults `PYTHONHOME=/usr` to prevent an unrelated uv interpreter's standard library from being selected. Set `EXTENDED_BLENDER_PYTHONHOME` when using a different Blender Python installation. NumPy is installed in a project-local cache for Blender, leaving system packages unchanged. Application display serialization uses the project environment's NumPy and trimesh.

The builder's default output is an isolated manifest fragment, leaving the shared app manifest untouched. Integration is explicit and idempotent:

```bash
.venv/bin/python scripts/build_extended_anatomy.py --integrate-only --append-manifest data/derived/app/manifest.json
```

The append operation replaces only `z-anatomy` records, copies this family's geometry to the app geometry directory, and preserves unrelated families. The verification script tests that behavior in a temporary directory, validates every stored checksum and face index, verifies local/world/display transformations, checks source/display triangle equality, and requires evaluated procedural nodes and mirrored anatomy to be present.
