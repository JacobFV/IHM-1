# Pressure command topology

Every new embodied frame reports `input_capabilities.skin_pressure` with `whole_skin` and exact `regional_ids`. The runtime derives these from the current native owner's callable commands, explicit whole-skin support flag and regional membership. Unknown owners report no pressure capability. Initial snapshots, accepted steps and subsequent snapshots carry the same owner contract.

The frontend passes the actual accepted frame into `bodyCommand`. `skinPressureCapability(frame)` exposes whole, regional or unavailable status and an actionable message. Missing capabilities are not inferred from tissue observation availability, region label patterns, UI defaults or selected configuration.

For a whole-skin owner, commands retain whole pressure including zero release. For a regional owner, the panel's neutral whole-pressure zero is omitted; a nonzero whole pressure fails before any request, with instructions to reset it and use regional boundaries. Explicit `regional_skin_pressures` are validated against exact owner IDs and the backend's finite 0–5000 Pa domain. Serialization never includes both pressure topologies. Invalid or missing capabilities cannot accept nonzero pressure. Other mechanical/neural commands remain usable when pressure is zero and its capability is unavailable.

Validation: 15 focused JS command/live tests and 28 Python capability/runtime/actor/HTTP tests pass. No native process, build, browser trajectory or physiological advance was required. Existing server processes without the new capability field remain pressure-unavailable until restarted; the client does not guess their topology.
