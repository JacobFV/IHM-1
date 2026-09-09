# Current surface support rematerialization

The fresh alternative support artifact is
`data/derived/supine-surface-contact-5jqy1juo/manifest.json`. It is experimental;
production support defaults were not changed.

The builder now accepts explicit retained native reference, contact XML and native
execution receipts. This artifact uses the current 98-muscle source reference,
current canonical metadata, verified source skin geometry and explicit skin-layer
thicknesses. It also honors the current `physical_surface_support` domain and its
component evidence hash.

The freshly rasterized first-hit seam cell is index **15942**, source triangle
**165251**. It is excluded by the current exterior-component domain. Its 0.000025 m²
cell is omitted without substituting the much deeper opposite ray hit; the
missing support remains explicit. Copying the old index 15948 would have removed
the wrong current cell. The artifact retains 21,382 quadrature points and
0.53455 m² projected area. Native forces replace the old sphere owner; they are
not added to it.

The initial geometry still has substantial gaps: pelvis 39 mm, femurs 88–89 mm,
tibias 112 mm above the reference plane. A rigid 6.6 mm skin layer cannot by itself
provide distributed support across those gaps.

At 77.6122 kg (761.38 N weight), actual native runs with identical source default
coordinates and no world/blanket forces produced:

| Support | Outcome | Support at 0.1 s | Additional observation |
|---|---|---:|---|
| Default spheres | Completed 0.1 s | 624.49 N | Still moving |
| Rigid skin foundation | Native integrator stopped at 0.03807 s | — | No domain guards relaxed |
| Skin + retained MM mattress | Completed 0.1 s | 269.88 N | 41.02 mm bed compression; 1.249 mm skin compression |

MM remains in a falling transient at 0.1 s, with support only 35.4% of weight.
These runs do not establish equilibrium or a better initial posture. Report:
`data/derived/surface-support-current-xhc7djfd/report.json`.

Controlled native/Python MM static force and elastic-energy parity, reciprocal
resultants, dissipative power and a 2 ms native step pass. Report:
`data/derived/native-surface-foundation-e55d0sak/report.json`.

Reproduce source-domain fixtures and the native comparison:

```sh
PYTHONPATH=. .venv/bin/python scripts/verify_supine_surface_support.py
PYTHONPATH=. .venv/bin/python scripts/verify_current_surface_support.py \
  --manifest data/derived/supine-surface-contact-5jqy1juo/manifest.json
```

`build_supine_surface_contact.py --reference-native ... --contact-force-set ...
--reference-execution ...` creates a new artifact instead of retagging old hashes.
The current retained reference is
`data/derived/surface-reference-current-3xvkopcl/initial_native.json`.
