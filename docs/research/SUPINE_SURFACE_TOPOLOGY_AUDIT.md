# Frozen supine contact topology and thickness migration

The retained materialization `data/research/engineered_skin_territories/materialization.json` identifies component 0 as an inferred exterior proxy (109,183 faces, 1.780460255 m²), component 1 as the oppositely oriented nested inner shell (94,055 faces), and 98 tiny seam components. This is an engineering eligibility selection, not a closed, intersection-free surface certificate; the exterior has 827 boundary edges.

The read-only `scripts/audit_supine_surface_components.py` verifies geometry, topology and frozen quadrature hashes and maps all 21,381 sampled source face IDs. **21,380 samples belong to exterior component 0, one to seam component 45, and none to the nested inner shell.** Minimum-X raster selection therefore prevents broad double-shell contact duplication in this artifact, but does not exclude every seam.

The exception is quadrature index 15,948, source face 165,251, assigned to torso, with projected area 0.000025 m² (about 0.00468% of total quadrature area). Its reference source point is [−0.184211453, 1.192523741, −0.146564569] m. Restricting that same ray to eligible exterior faces selects face 164,276 at X=−0.065089816 m, **119.121636 mm deeper**. The full and exterior-only minimum-X plane is unchanged at −0.268803031 m. Because the exterior proxy is open, a much deeper next intersection cannot be assumed to be an anatomically repaired posterior surface. Excluding the seam and choosing replacement versus omitting the unresolved cell require an explicit geometry policy.

This audit does not emit individual posed contact forces; it does not claim the exceptional point has any particular force in the last static candidate. Replacing or removing it changes the contact geometry and input identity. Existing frozen receipts remain historical evidence; their exact evaluation cache must not be promoted across that change. No canonical artifact, quadrature, native input or reference was regenerated during the audit.

## Future layer thickness

`build_supine_surface_contact.py` now reads `entity['shell']['thickness_m']` directly and retains its shell/prior/surface-support provenance in each output layer. A layer carrying `physical_surface_support` but missing explicit thickness fails even if a legacy option is supplied. Legacy artifacts without carried shell metadata can use volume/raw-area only through the explicit `--legacy-thickness-basis legacy_raw_full_skin_area_volume_ratio` option, and the output records that legacy basis and limitation.

The fixture uses the actual raw/exterior areas 3.502598974493317/1.7804602548390722 m². Explicit epidermis/dermis/hypodermis thicknesses 0.1/1.5/5 mm continue to total **6.6 mm** after layer volumes use exterior area. It catches the erroneous 3.354948 mm result from dividing corrected volume by raw double-shell area, rejects missing new thickness, and requires an explicit legacy basis. Existing envelope, force/energy, action/reaction and compression-domain fixtures also pass. No existing thickness or constitutive parameter changed.

Evidence SHA-256:

- `data/derived/supine-surface-contact-exmzq9pq/component_audit.json`: `6187701b8b437706a0382eb2e2b165ff77d8a0ad39c94043b546390677278d0c`
