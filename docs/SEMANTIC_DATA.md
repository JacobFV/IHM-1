# Anatomical and cellular evidence

The imported [HRA ASCT+B record graph](https://github.com/x-atlas-consortia/hra-pop/blob/main/output-data/v0.11.1/reports/hra/asctb-records.csv.zip)
contains 1,759,409 structure–cell–biomarker rows. It is indexed locally in
`data/derived/semantics/human-reference-atlas.sqlite` with original labels and
normalized FMA/UBERON identifiers. A cell-marker association is not a kinetic
coefficient or an interaction edge.

The [HRA VCCF repository](https://github.com/hubmapconsortium/hra-vccf) supplies
1,082 named vessel records, organ crosswalks and 193 literature-derived geometric
measurements with population, sex, sample-size and method fields. Data are CC BY
4.0; source SQL is MIT. Revision and file hashes are retained.

The collector found 191 upstream `Geometry.csv` rows with extra CSV columns in
the citation suffix. The well-formed identity/numeric/unit/reference-URL prefix
is retained after validation, while the remaining citation tail is preserved
verbatim as parsed raw fields and marked malformed. It is not silently repaired
into a fabricated citation. Primary studies are not independently verified merely
because the atlas cites them.

BodyParts3D FMA annotations provide 461 vessel-to-mesh annotation candidates.
These include ancestor memberships and are explicitly **not coordinate
registrations or guaranteed anatomical equivalences**. Branching relations also
remain anatomical relations; venous flow is not assigned from tree direction.

```bash
.venv/bin/python scripts/index_semantics.py
```
