# IBM-1 ↔ IHM-1 nerve join, contract v2

The generated `data/derived/canonical/peripheral.json` declares 144 sided nerve
records: all 71 IBM trunk names bilaterally, plus the preserved bilateral
`sciatic_fibular` identity. Strip `peripheral-nerve-{side}-` for the join and
retain `sciatic_tibial → sciatic`. No IDs were renamed. `sciatic_fibular` needs
an explicit IBM composition decision; it is not silently aliased.

Read **`nerves[].path_length_m`** for a representative route per nerve and side.
The old bridge only collected muscle and skin bindings and would substitute IBM
trunk lengths for the new routes. Apply [the tested bridge patch](patches/ibm-ihm-bridge-v2.patch)
in IBM-1 with `git apply /path/to/IHM-1/docs/patches/ibm-ihm-bridge-v2.patch`.
It consumes the nerve-level field, propagates limitations and targets, rejects
missing/nonpositive/nonfinite lengths and unknown names, and explicitly warns
for the two preserved unmatched sciatic-fibular routes. The patch was exercised
against the local IBM table: 142 joined routes, no trunk-length fallback.

## What the length means

`path_length_scope: representative_endpoint_to_relay` means a full **inferred
path between the declared endpoint and relay**, not the length of the named
trunk alone, and not receptor-to-cortex. Positions and lengths are in metres;
delays are in seconds; velocities are in metres/second. Synaptic delays and
central transport beyond the relay are excluded. Muscle-specific lengths remain
on muscle bindings; a nerve-level representative must not replace them when the
actual endpoint is known.

The earlier description of IHM lengths as measured mesh routes was incorrect.
Existing binding estimates are `max(0.03, 1.15 * Euclidean(endpoint, relay))`.
Skin endpoints are selected on a mesh, but the intervening nerve is not traced
on that mesh. Existing nerve representatives select the upper median binding.
New routes sum segments between authored landmarks. Their coordinates and lengths
are modelling hypotheses, not quantities taken from dissections or literature.
Every nerve and binding carries evidence/geometry flags and limitations.

A branching nerve, plexus or sympathetic chain has no unique whole-route length.
The catalog names one representative distal endpoint for each new route. This
is coverage of route **declarations**, not validated full embodiment. In
particular, a visceral afferent path does not specify the lengths of the separate
preganglionic and postganglionic efferent stages. Do not apply one class velocity
to a path spanning those two neurons and interpret it as an autonomic loop.

## Special senses and viscera

Optic ends at an authored lateral-geniculate relay; cochlear and vestibular end
at separate brainstem relay groups; olfactory ends at the bulb. Downstream
cortical IDs appear in `downstream_brain_target_ids`, separately from the length:
pericalcarine, transverse temporal, insula and entorhinal, respectively. Both
hemispheres are listed as coarse downstream annotations, with no resolved
projection weights, decussation or central delay. New routes have
`runtime_support: topology_only`; no special-sense transducers are implied.

The [CN VIII reference](https://www.ncbi.nlm.nih.gov/books/NBK537359/)
supports distinct brainstem relays and downstream auditory cortex. The
[pelvic dissection study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9720663/)
supports the distinction of pelvic splanchnic pathways. Neither supplies the
coordinates or lengths authored here. HumMod's local `VagusNerve.DES` was
inspected: it supplies physiological activity/control declarations, not a
geometric nerve trace. No HumMod geometry was inferred from those rates.

## Fibre timing

IBM remains the source of the velocity table. The build preserves its exact
`nerve.py` bytes and SHA-256 under `peripheral-sources`, extracts the velocity
constants without executing IBM code, and records them in `fibre_velocity_m_s`.
The table is an IBM prior snapshot, not an independently validated IHM source.
IBM continues to own trunk composition and the full per-class split.

For IHM's muscle runtime, bindings additionally emit pure conduction
`delays_s` for Ia, Ib, II, alpha and gamma. `BodyPeripheral` reads Ia and alpha
explicitly. The retained scalar fields are compatibility aliases for those
specific classes plus the separate 12 ms central prior. They do not describe
all classes in a trunk. Cutaneous warm/cold/pressure timing remains separately
parameterized. The proprioceptive response remains an illustrative Ia-timed
stretch/force proxy, not independent spindle, tendon-organ or gamma dynamics.

`ReflexParameters.from_binding` computes Ia + alpha + an explicit synaptic
prior (default 1 ms), without adding cortical transport. The canonical TA
experiment uses that factory; its delayed contrast adds 40 ms. This changes
experiment timing from the old fixed 20 ms loop. Gamma has its own conduction
delay but no fusimotor controller is claimed.

## Rebuild and checks

Run from IHM-1, with the sibling IBM-1 checkout available:

```sh
.venv/bin/python scripts/build_body_peripheral.py
.venv/bin/python scripts/verify_peripheral_join.py
.venv/bin/python scripts/verify_body_peripheral.py
PYTHONPATH=. .venv/bin/python scripts/verify_body_reflexes.py
```

Derived data is ignored by git in this project. The catalog, enrichment code,
contract, bridge patch and checks are source files; rebuilding writes the
canonical artifact consumed by IBM.
