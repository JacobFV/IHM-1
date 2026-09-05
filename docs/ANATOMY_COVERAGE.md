# BodyParts3D anatomy coverage

Regenerate with `.venv/bin/python scripts/audit_anatomy_coverage.py`; verify with `.venv/bin/python scripts/verify_anatomy_coverage.py`.

All 2,234 acquired source meshes have an explicit classification status. 0 are unresolved; 230 display groupings changed from the historical substring classifier.

The crosswalk uses exact recorded FMA concept anchors and directed child-to-parent paths in the downloaded is-a and part-of tables. Is-a identity takes precedence over part-of context. All supported system tags remain available; a primary display group does not imply an exclusive biological function. Each mesh record retains source concept IDs, paths with labels, six source-table SHA-256 hashes, a crosswalk hash, and its source OBJ hash.

| Primary system | Historical | Source crosswalk |
|---|---:|---:|
| arterial | 621 | 639 |
| cardiac | 62 | 23 |
| connective | 27 | 48 |
| digestive | 70 | 110 |
| endocrine | 5 | 4 |
| integumentary | 4 | 4 |
| lymphatic | 3 | 3 |
| muscular | 403 | 425 |
| nervous | 107 | 146 |
| reproductive | 10 | 12 |
| respiratory | 103 | 101 |
| sensory | 0 | 33 |
| skeletal | 283 | 285 |
| urinary | 5 | 6 |
| venous | 360 | 395 |

## Corrected examples

| Source mesh | Recorded name | Historical | Source crosswalk |
|---|---|---|---|
| FJ1961 | posterior intercostal arteries | other | arterial |
| FJ1975 | posterior intercostal arteries | other | arterial |
| FJ2090 | set of perforating arteries | other | arterial |
| FJ2418 | anterolateral head of lateral papillary muscle of left ventricle | other | cardiac |
| FJ2419 | anterior papillary muscle of right ventricle | other | cardiac |
| FJ2429 | lateral papillary muscle of left ventricle | other | cardiac |
| FJ1291 | left common tendinous ring | other | connective |
| FJ1324 | tarsal plate of left upper eyelid | other | connective |
| FJ1328 | tarsal plate of left lower eyelid | other | connective |
| FJ1252 | gingiva of upper jaw | other | digestive |
| FJ1253 | gingiva of lower jaw | other | digestive |
| FJ1858 | hepatovenous segment ix | venous | digestive |
| FJ1462 | set of right levatores costarum breves | other | muscular |
| FJ1462M | set of left levatores costarum breves | other | muscular |
| FJ1463 | set of right levatores costarum longi | other | muscular |
| FJ1734 | anterior commissure | other | nervous |
| FJ1735 | brachium of left superior colliculus | other | nervous |
| FJ1736 | brachium of right superior colliculus | other | nervous |
| FJ3135 | left deferent duct | other | reproductive |
| FJ3140 | right deferent duct | other | reproductive |
| FJ1282 | anterior chamber of left eyeball | other | sensory |
| FJ1285 | left choroid | other | sensory |
| FJ1286 | left choroid | other | sensory |
| FJ3153 | xiphoid process | other | skeletal |
| FJ3178 | body of sternum | other | skeletal |
| FJ3148 | urethra | other | urinary |
| FJ2656 | great cardiac vein | cardiac | venous |
| FJ2657 | anterior interventricular vein | cardiac | venous |
| FJ2658 | anterior interventricular vein | cardiac | venous |

## Ambiguity and contextual membership

2008 meshes have is-a identity evidence; 226 rely on part-of context. 4 primary assignments have equally near competing system anchors.

Pancreatic structures retain digestive and endocrine membership. Source part-of aggregates also associate some cardiac structures with the systemic arterial system and liver with venous structures; explicit is-a anchors for valve cusps, chamber cavities, papillary muscle regions and liver segments take precedence. Aggregate tags describe recorded source associations, not an assertion that a valve is an artery or a liver segment is a vein.

## Skin and remaining anatomical limits

Source skin `FJ2810` has 102,467 vertices and 203,382 triangular faces; bounds in the official millimeter frame are `[[-334.119, -246.783, -78.1112], [332.825, 45.2476, 1641.36]]`, with axis spans `[666.94, 292.03, 1719.47]` mm. This is a source surface shell, not layered epidermis/dermis or a skin physiology model. Display decimation must be reported separately from these source counts.

BP3D lymphatic coverage here consists of spleen and two thymus lobes. No lymph-node or lymph-vessel meshes are present in these downloaded element tables. Three classified lymphatic surfaces do not constitute a lymphatic network. Additional sources must remain separate coordinate frames until registration is validated.

The sensory display group covers eye surfaces, chambers, lacrimal structures and external ear. Nerves retain nervous identity; associated muscles retain muscular identity. Thyroid cartilage is skeletal, sternothyroid is muscular, and inferior thyroid artery is arterial. Papillary cardiac muscles retain cardiac identity. Generic anatomical terms without a recorded path remain explicitly unresolved rather than being inferred from substrings.

This classification audit establishes source inventory and traceability. It does not establish tissue completeness, watertightness, subject calibration, or physiological validity. The downloaded BP3D ontology is a selected representation, not the complete FMA release.
