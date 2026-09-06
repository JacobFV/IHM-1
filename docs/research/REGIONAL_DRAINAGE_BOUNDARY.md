# Regional traction to native drainage boundary contract

`RegionalDrainageBoundary` in `ihm/assembly/regional_drainage.py` maps a registered
traction field to an available native command only when the field is representable
by that boundary. It reports unsupported local loading without averaging it into
lumped Skin or pretending to predict regional drainage.

## Actual native ownership

In held BioGears source `projects/biogears/libBiogears/src/engine/Controller/BioGears.cpp`
(lines 4479–4510), `SkinE3ToGround` is the extracellular compliance reference.
`SkinE3ToSkinL1` supplies lymph drive pressure; `SkinL1ToSkinL2` supplies lymph
resistance; `SkinToLymphValve` joins SkinL2 to Lymph. Skin extracellular maps E1,
E2, E3, L1 and L2 into **one liquid compartment**. These node names do not create
independently owned regional lymph volumes. `SkinIToGround` is a separate
intracellular compliance boundary.

The audited `scripts/native_tissue_compression.h` accepts only Skin. It replaces
the active E3-to-Ground compliance path with `IHMSkinE3ToExternalPressure`, creates
`IHMSkinExternalPressure` and drives it through `IHMGroundToSkinExternalPressure`.
The detached native compliance object remains available to native laws, and the
adapter copies its scalar state before and after the native solve. The existing
`CoupledNativeSession.skin_compression(pressure_pa)` command supports 0–5000 Pa.
Oncotic paths and intracellular references are not regional pressure controls.

## Mapper input and output

Construct `RegionalDrainageBoundary(cells, owner_scope_receipts)` and call
`.map(tractions, native_snapshot)`. Each registered cell supplies `id`,
`native_owner` (e.g. `Skin.extracellular`), `area_m2`, unit `outward_normal`, and
`source_identity={mesh_sha256, triangle}`. Mesh/triangle pairs must be unique
across all owners and display IDs. Use non-overlapping triangles of one registered
surface discretization; do not register overlapping source atlases or mesh LODs as
additional area. Geometry registration and source-hash verification remain the
caller's responsibility; this mapper validates the supplied identity/area contract,
not mesh coordinates or actual physiological territory.

Every traction supplies `cell_id`, matching `area_m2`, and a three-component
`traction_pa`. Positive compression is `-traction · outward_normal`. The report
retains normal and tangential traction, each vector force `area × traction`,
normal load, summed owner area and resultant force. Missing cells remain unknown;
there is no implicit zero loading. Duplicate identities and changed areas fail.

A command candidate requires all registered cells for that owner, an explicit
`owner_scope_receipts[owner]` with `complete_native_owner_surface=True` and
`source_identity`, uniform nonnegative pressure, negligible tangential traction,
and observed finite `tissue.compression.Skin.requested_pa` / `applied_pa` ports in
the supplied native snapshot. The scope receipt is a caller assertion grounded in
its registration; normalized atlas surface weights alone cannot establish it.
Uniform pressure is recovered by area quadrature and its normal-load residual is
reported. This is not a homogenization rule for nonuniform contact.

The output `native_command_candidates` names the session method and actual native
node/path names. The mapper never sends commands; root integration must apply any
accepted candidate to the same native session whose snapshot supplied capability.
Commands are whole-owner boundaries and create no regional native state. Partial,
nonuniform, tensile, shear, excessive-pressure, missing-capability and other-organ
loads retain their area/force audit under `unresolved`. Merely constructing an
anatomical regional lymph coverage report does not enable a pressure command.

## Precise native extension needed for local drainage

A local Skin patch with a different pressure cannot be represented by the existing
single E3 pressure/compliance and single extracellular substance inventory. The
required extension is `native_regional_interstitial_pressure_boundaries`:

1. Register stable region IDs with source support receipts, disjoint volume/mass
   allocation fractions and an explicit residual region. Surface-area fractions
   alone are insufficient physiological compartment calibration.
2. Create regional E3 fluid nodes and per-region external reference pressure nodes;
   preserve the selected parent's total baseline/current/next volume and compliance.
   A new endpoint would accept `{region_id, external_pressure_pa}` and return the
   exact region-to-node/compliance-path receipt. It must reject unmapped IDs.
3. Partition the native liquid ownership and Albumin/other transported masses once,
   including the residual. The old parent must become an aggregate view or be
   replaced; it must not retain an independent duplicate inventory. E1/E2 vascular
   exchange, E3-to-I, E3-to-L1, lymph resistance/valve and sweating paths all require
   explicit native mapping. L1/L2 are currently mapped into Skin extracellular;
   copying their names is not a regional ownership solution.
4. Preserve Tissue/Diffusion references, oncotic/osmotic updates and the active
   cardiovascular solve lifecycle. Add region observation ports for external and
   interstitial pressure, compliance, volume, species mass and drainage flow, with
   native-owned protein transport receipts rather than Python concentration × flow.
5. Before enabling regional commands, verify zero-pressure parity, summed native
   fluid/species ownership and native integration closure, locality under a patch
   load, unload recovery, and native law pointer validity. No such regional solver
   or physiological calibration is claimed by this mapper.

No shared native header is changed by this increment. The extension is a concrete
missing interface, not an implemented or validated regional drainage solver.

## Bounded verification

`.venv/bin/python scripts/verify_regional_drainage.py` runs six small tests with the
existing tissue snapshot fixture. They check uniform force/area preservation,
nonuniform rejection, incomplete coverage and missing scope/ports, duplicate source
ownership, area/normal validation, unsupported shear/tension/range and other native
owners. No native compilation or simulation is run.
