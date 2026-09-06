# Regional native Skin in signed physiological sessions

`ihm.native.regional_session.RegionalSignedNativeSession` uses the separately
built `native_biogears_regional_signed` adapter and the verified
`whole_body_integrity_regional_skin_graph_v2` library. Existing adapters remain
available. Construction requires a retained native state; fresh patient
stabilization with the added topology is not accepted.

The owning regions are `region_a`, `region_b`, and `residual`, with fixed
engineering volume/compliance fractions 0.2/0.3/0.5. These names do not identify
anatomical territories. The original Skin extracellular compartment becomes an
aggregate view, with no independent duplicate inventory.

`regional_skin_pressure(region, pressure_pa)` accepts 0–5000 Pa, acknowledges
without advancing time, and supplies the next native solve's external pressure.
` signed_step(reference_id, delta_m_w, delta_h_w, delta_w_w)` retains the signed
muscle/heat/work protocol. Ordinary stepping rejects after reference binding.
Whole-Skin compression and state serialization reject because their original
contracts do not describe the new topology.

Values under `tissue.regional_skin.<region>.*` report actual native volumes,
pressures, requested/applied external pressures, path flows and species masses.
Aggregate ownership, per-step fluid incidence, and native sweat donor/waste
ledgers are separate observations. Sixteen detached legacy Skin circuit fields
remain unavailable; no old path values or guessed zeros replace them.

The native wire fixture at
`data/runtime/physiology/regional-signed-adapter-axpxzdgk/protocol/verification.json`
passes five steps and 15 commands, including load/unload, signed and respiratory
coupling, invalid commands, and ownership checks. The actual Python session
fixture at `data/derived/audits/regional-python-vvdedt24/verification.json`
passes two signed steps and regional load/release in 0.718 s wall time.

This does not yet enable regional Skin in the default embodied factory. Its
exchange views must represent aggregate versus owning compartments explicitly,
and spatial registration must identify any physical tissue territory before
local body contact can supply these compartment pressures. Long-run and
whole-body mechanical equivalence remain unverified.

The shared runtime accepts `regional_skin_pressures` when its native owner
supports those boundaries. It validates all region IDs/values and incompatible
whole-Skin topology inputs before advancing any owner. Commands are delivered
before the signed physiological exchange; an uncertain command terminates the
session under the existing native-outcome policy. Unspecified pressures retain
native boundary state; release is an explicit zero command. This input interface
adds no anatomical region assignment or reconstructed contact pressure.
