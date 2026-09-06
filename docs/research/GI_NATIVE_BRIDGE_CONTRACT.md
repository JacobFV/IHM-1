# GI ownership bridge toward native integration

`ihm/assembly/gi_native_bridge.py` provides an executable, stateless transaction plan. It consumes authoritative snapshots, returns paired water/Na/K/Cl deltas and fecal receipts, and holds no persistent physiological inventory. It composes the measured conditional colonic response with the conservative lumen-transit foundation. This is useful implementation of the donor arbitration and transfer interface; it does not silently install incomplete chemistry into native physiology.

## Exact ownership and scheduling boundary

The inspected native `Gastrointestinal::PreProcess` performs digestion, `ChymeSecretion`, and `AbsorbNutrients`. The latter manually mutates solute masses and the SI circuit source node's `NextVolume`, while the fluid circuit later resolves volume transport. Therefore the bridge snapshot must be taken **after the existing SI uptake calculation**, using consistent prospective volume and already-updated mass, before downstream circuit commit. Reading current volume beside next-step solute mass would produce an inconsistent concentration. The eventual adapter must establish the exact current/next scalar convention with a tiny native fixture before activation.

Five authoritative owners are required in this order:

| Owner | Native binding or proposed addition | Paired boundary |
|---|---|---|
| SI lumen | Existing `SmallIntestineChyme`, mapped `SmallIntestineC1` | Debit SI→colon; existing wall uptake is never repeated |
| Colon lumen | New native fluid node/compartment, serialized once | Credit SI inflow; signed wall exchange; debit rectal transit |
| Rectum lumen | New native fluid node/compartment, serialized once | Credit colon inflow; debit fecal boundary |
| Colon blood | Existing `LargeIntestineVasculature` / `LargeIntestine1` | Equal opposite signed wall exchange; downstream portal flow stays native |
| Fecal ledger | New explicit cumulative external water/substance output | Credit exactly the rectal debit; never erase at Ground |

Distinct owner tokens are mandatory; a single blood/chyme pointer cannot masquerade as multiple reservoirs. The bridge has no default initial colon contents. Initial contents must be measured or labeled as engineered priors and added to the global initial mass/volume ledger once.

## Executable conservative transaction

The planner applies colon wall exchange before downstream transit. This priority is a declared engineered operator-splitting choice, not a measured timing law. The empirical signed water/Na/K/Cl vector is capped by one common available-donor time fraction, including vascular donors for secretion. Transit then draws from the remaining aqueous stores using donor composition. All three lumen edges use the same phase snapshot, so a newly arriving packet cannot pass several regions in one step. A future joint solve can replace the splitting choice without changing ownership.

Snapshots carry a native timestep epoch and a hash of owner bindings, phase, inventory values and current-ledger owner. `validate_commit` rejects a stale or changed snapshot. This is a read-only precondition; a future native writer must perform the check and mark the transaction consumed inside one mutation boundary. The Python function alone cannot guarantee atomicity or prevent replay against an unchanged snapshot. Native integration must never both set a circuit flow that transports a quantity and also manually debit/credit that quantity again.

## Explicit incomplete pathways

The implemented payload is limited to Na/K/Cl and water. A native transit packet must enumerate **all** dissolved substance quantities, including glucose, amino acids, TAG representation, calcium, urea and applicable drug quantities; generic transport cannot drop unrecognized cargo. Per-drug CAT arrays must either participate through their existing owner or be explicitly excluded from the same physical-lumen interpretation. They are not an extra aqueous-volume source. Solid residue requires separate handling.

The empirical colonic vector generally has nonzero `Na + K − Cl` equivalents. A named current ledger owns that residual for diagnosis; naming it does not complete a physical current path. H/bicarbonate, membrane charge/current and coupled exchange mechanisms remain missing. The planner therefore returns `native_commit_ready=False` and a concrete missing-interface list. No arbitrary chloride compensation or extra ATP/heat debit is generated. Native tissue/energy remains the eventual sole metabolic owner.

An engineered prior can already be used through the conditional response's explicit transfer-prior flag outside its measured perfusion protocol. That enables bounded exploration without pretending universal calibration. It does not bypass inventory ownership, complete species transfer or electrical conservation requirements.

## Verification and next native increment

`.venv/bin/python scripts/verify_gi_native_bridge.py` checks finite water/species closure across all five owners, competing wall/transit donor use, no repeated SI wall uptake, no same-tick rectal pass-through, stale snapshot/phase rejection, duplicate-owner rejection and explicit charge-gap ownership. It does not execute native physiology.

The next native increment is a tiny isolated fixture exposing consistent `NextVolume`/solute snapshots and atomic paired delta application to actual compartments, with zero-baseline colon/rectum additions and fecal serialization. First activate complete-species **lumen transit only** with prescribed finite volume requests; this transfer does not need an invented epithelial current. Then implement measured signed wall exchange only after explicit charge/bicarbonate/current ownership is available. This sequence produces useful native transit independently of universal absorption calibration, while preserving one store per physical compartment.
