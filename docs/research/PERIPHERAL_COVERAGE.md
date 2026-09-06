# Peripheral coverage of the registered body

`ihm.assembly.peripheral_coverage` provides a static, machine-readable audit of
`data/derived/mechanics/whole_body_arm26_v2/registration.json`. This is the actual
92-muscle XML/catalog assembly: 80 retained walking-model muscles plus 12 bilateral
Arm26 shoulder/elbow actuators. Coverage describes implemented ports; it does not
establish that a session is running or that the pathways are biologically validated.
The audit neither changes IBM-1 nor launches native mechanics or neural integration.

Generate an artifact using the existing Python environment:

```sh
PYTHONPATH=. .venv/bin/python -m ihm.assembly.peripheral_coverage \
  --output data/derived/mechanics/whole_body_arm26_v2/peripheral_coverage.json
PYTHONPATH=. .venv/bin/python scripts/verify_peripheral_coverage.py
```

Runtime/UI consumers can call `build_peripheral_coverage(root,
controller_catalog=controller.catalog)` and cache its result by `audit_sha256`.
`write_peripheral_coverage(root, output_path, registration=...)` validates before
writing. The optional controller catalog must equal the verified registered catalog,
including provenance and assignments. An ankle fixture or the default 80-muscle
controller is deliberately rejected for this 92-muscle plant. Construct the controller
with `whole_body_effector_catalog(root, manifest)` when targeting the registered body.
No integration changes to existing neural files are required by this module.

The `ihm.peripheral-coverage.v1` JSON contract contains:

| Field | Meaning |
| --- | --- |
| `effectors` | Exact native muscle IDs, source hashes, attachment bodies, receptor IDs, afferent/efferent mappings and separate spinal-reflex support |
| `receptors` | 184 muscle length/force observation ports and two foot-load stance gates; these IDs identify software ports, not anatomical receptors |
| `brain_region_ids` | Canonical population IDs available to the bridge |
| `unsupported_paths` | Stable gap IDs with explanatory reasons |
| `source_hashes` | Stream-verified preserved source bytes, model/catalog and current bridge/audit implementation identity |
| `summary` | 92 engineered afferent/efferent targets, eight source reflex targets, zero anatomically measured connections |
| `audit_sha256` | Deterministic hash of the complete report before this field is added |
| `runtime_active` | Always false for this static audit |
| `registration_native_verified` | The manifest's claim, carried verbatim; not a new native verification |

The source-derived reflex primitives are the bilateral tibialis anterior, soleus,
medial gastrocnemius and lateral gastrocnemius routes in `sensorimotor.py`. Tibialis
anterior uses its length and ipsilateral soleus force; gastrocnemius heads use the
explicitly lumped two-head force transfer. Foot load determines stance. The pinned
Geyer/Herr PDF hash and implementation hash are recorded. These primitives are
source transfers with engineering assumptions, not measured subject-specific spinal
connections.

All 92 cortical afferent routes aggregate normalized mechanical length/force into
contralateral postcentral Hz. Efferent routes use contralateral precentral rates for
an engineered gain and requested-effector-gated drive. The preserved IBM laws own
regional neural dynamics; the bridge owns these engineering conversions. Native
mechanics owns activation and muscle force. No pathway is labeled as identified
motor recruitment. Fiber-length proxy observations must remain explicitly named and
must retain `sensor_basis`, as required by the controller.

The audit rejects stale model/catalog/source hashes, root-escaping sources,
nonexistent or duplicated effectors/brain IDs, invalid laterality, mismatched
cortical assignments, and attachment IDs absent from the actual XML. All canonical
brain-source files are hashed in a stream; donor meshes are not parsed or loaded as
geometry. It checks structural/provenance compatibility, not native force validity,
neural stability, receptor calibration or runtime latency.

Remaining systems include explicit spindle/Golgi/Ia/II/Ib populations, axons and
neuromuscular junctions; nerve roots/trunks and measured somatotopy; reflexes outside
the eight ankle targets; a whole-body cutaneous receptor census and calibrated
cutaneous-to-cortex conversion; cranial/special senses; IBM visceral/autonomic
pathways; and an identified autonomous motor policy. The separate regional rapid/
slow IBM touch materialization produces signed donor responses without an identified
conversion to cortical Hz, so it is recorded as an unsupported end-to-end route.
Native systemic autonomic control remains a separate owner. Hand/finger, neck,
facial, respiratory and independent scapular muscle effectors are not supplied by
this 92-muscle catalog.
