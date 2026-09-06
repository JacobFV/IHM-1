# Installed regional extracellular exchange observations

`ihm/assembly/regional_exchange.py::observe_regional_skin(snapshot, *,
native_identity=None, source_hashes=None)` reads the actual regional signed adapter
ports. It returns `available=False` when the entire regional prefix is absent.
A partial installation with missing volume, fraction, pressure, or circuit flow
fails closed. This helper executes no native process and integrates no state.

The three accounting keys are `Skin.region_a.extracellular`,
`Skin.region_b.extracellular`, and `Skin.residual.extracellular`. Their native
compartment names are `IHM_<region>_SkinTissueExtracellular`. The original
`Skin.extracellular` / `SkinTissueExtracellular` is returned in `aggregate_views`
with `accounting_owner=False`, `independent_store=False`, and the three children.
Skin vascular and intracellular owners, Lymph, and the remaining 12 organs remain
shared/original owners. The .2/.3/.5 installation fractions are engineering volume
and compliance assumptions; current volume and species quantities are read directly
and need not retain those proportions. No anatomical attachment is inferred.

For integration into `NativeTissueExchange.observe`, detect the regional prefix
before constructing required original Skin paths. Merge the three owning records,
retain the Skin aggregate only as a non-owning view, and omit its incidence entry.
Replace the original Skin extracellular spatial partitions with views of these
actual leaves; do not also allocate an independent Skin parent. Replace the old
`SkinE1ToSkinE2`, `SkinE3ToSkinI`, and `SkinE3ToSkinL1` transfers with the helper's
12 transfers. Add the helper's shared-owner rates to other selected-organ rates;
never count them twice. Original Skin path exports may be missing or null and
must not be fabricated to pass the legacy observer.

Each region preserves all nine path observations under actual native path names
`IHM_<region>_<original_path>` and exact snapshot source keys. The selected
extracellular fluid incidence is:

```
Q(SkinE2ToSkinE3) - Q(SkinE3ToSkinI)
  - Q(SkinE3ToSkinL1) - Q(SkinSweating)
```

E2ToE3 inflow is attributed to the shared Skin vascular boundary. E3ToSkinI reaches
the shared intracellular owner; E3ToSkinL1 reaches the shared Lymph boundary.
Serial E1/E2 and L1/L2/valve observations are retained without additional transfer
entries. E3ToGround is a compliance boundary, not an additional fluid sink.
`external.sweat` is a transfer endpoint and a separate external rate, never an
owning compartment. Signed flows are preserved, including reverse intracellular
flow. These are selected circuit incidence rates, not a complete vascular or
whole-body volume derivative.

Every exported species mass is retained with its source key and explicit ug-to-g
conversion. Missing masses remain `None`; seven core species remain represented
even if all their regional exports are absent. Concentration, molarity, and gas
partial-pressure fields remain `None`. No solute flux is reconstructed. Regional
stores continue to share the native aggregate Tissue/Diffusion chemical laws.
The helper retains all raw regional ports, including load/work and sweat ledger
observations, without treating those diagnostics as extra stores.

Independent child-sum volume/species residuals are checked against the aggregate;
legacy aggregate observations, when present, must agree. Native exported ownership
and installation residuals are checked too. Each check includes units and an
absolute tolerance `max(1e-8, abs(aggregate quantity)*1e-10)`. Unknown mass comparisons
are explicitly `unobserved`. Native last-step fluid residuals are required after a
completed step, optional at installation, and checked to 1e-8 mL; they do not
establish conservation across arbitrary snapshot intervals. No whole-body mass
closure is claimed.

Bounded verification: `.venv/bin/python scripts/verify_regional_exchange.py`.
Five test groups use the actual retained adapter lifecycle snapshots in
`data/derived/audits/regional-python-vvdedt24/`, plus controlled corruption and
missing-port cases. No new native run is required. Fixture SHA-256 identities:

| File | SHA-256 |
| --- | --- |
| initial.json | cb37d4e4f5460d48e2938a374ab5a92148890559f7070482926732b2903e0b4d |
| loaded.json | fea8f66846d3a12f1d7cea48cecd1a6ff935f3220248c6528887b5d4efccd922 |
| released.json | a80f276241564280e431cde0436ba71a8794683e08a17486740f734f38e80337 |
