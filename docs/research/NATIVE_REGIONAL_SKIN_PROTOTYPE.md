# Opt-in native regional Skin prototype

This increment implements a native circuit/compartment split and verifies it on a
tiny executable native fixture. It does not activate a patient adapter or claim
validated regional chemical transport.

## Implemented topology and ownership

`scripts/native_regional_skin.h` defines `ihm::NativeRegionalSkin`, accepting an
existing native fluid circuit, compartment manager and Skin extracellular owner.
It preflights the actual five-node ownership (SkinE1, E2, E3, L1, L2) and nine
attached paths, rejecting unregistered connections or unsupported path states.
It creates region_a, region_b and residual with engineering fractions 0.2, 0.3 and
0.5. These fractions are explicit assumptions, not measured tissue territory.

Each region gets native E1/E2/E3/L1/L2 nodes and all corresponding capillary,
oncotic, intracellular-flow, drainage, lymph-resistance/valve, compliance and
sweating paths. Vascular supply, SkinI, Ground and Lymph remain shared endpoints.
Regional E3 compliance references a regional external-pressure node. Pressure
sources retain their values; resistance is divided by the region fraction;
compliance, fluid volume and flow sources are partitioned. Baseline/current/next
volume, compliance and law scalars are copied when available. The residual gets
remaining extensive quantity, preventing intentional mass/volume duplication.
Floating-point parity is verified with explicit tolerances, not bitwise identity.

The original extracellular species masses move once to three native liquid child
compartments. The original parent's node mappings are removed before `AddChild`;
its native mass and volume getters then produce child aggregates. Parent mass is
read-only, and no parent node remains active. Native species have one current mass
scalar; they have no separate baseline/next mass slots to fabricate. Unknown mass
states remain unknown. Water fraction and pH are retained when available.

Detached original paths stay in the circuit manager so native Tissue,
Cardiovascular and Energy law pointers remain valid. `after_preprocess()` fans
out their current laws and applies region pressures. `after_postprocess()` commits
detached law scalars and synchronizes original path flow caches with the summed
regional native flows. These caches are observations/law objects, not active
reservoirs. The host must install the split after native state load/setup and call
both hooks around every native solve; this prototype is not wired into the current
patient adapter or exposed as a user-load command.

## Actual architectural failure and correction

Two native Diffusion reads construct Skin path names and query the **active**
cardiovascular circuit: capillary Albumin flow through SkinE2ToSkinE3 and passive
lymph flow through SkinE3ToSkinL1. Removing the originals from the active solve
would make those reads dereference null. The isolated source transformation
changes exactly those two lookups to `GetCircuits().GetFluidPath(name)`, retrieving
the retained aggregate-law flow objects. The native transport algebra and corrected
shared-donor mass allocation remain unchanged. Two source tests reverse the exact
substitutions and recover the inherited source byte-for-byte.

The second native invariant is that `SELiquidCompartment::AddChild` rejects a
node-mapped parent. The header fixes this by transferring node mappings into
children before aggregation. The executable fixture verifies the parent has three
leaves, has no mapping, and returns read-only summed species masses.

Preserving aggregate source transport laws initially means region-specific solute
transport is **not** established by this variant. Native Diffusion distributes
aggregate changes using its existing mass/volume-weighted child operations; local
pressure is resolved in the fluid circuit while chemical constitutive feedback
remains based on the parent tissue law. Region-specific Albumin drain and oncotic
feedback need a subsequent native-law extension and native mass-ledger tests.

## Immutable builder and lineage

`scripts/build_biogears_regional_skin_variant.py` inherits
`whole_body_integrity_gi_absorption` (the signed-muscle lineage plus the native GI correction), verifies its pinned manifest/library and
all 360 object hashes, preserves its signed header receipt, and copies the
regional header into a fresh output directory. The inherited Diffusion source must
match the pinned raw source plus the existing shared-extracellular-donor fix. Only
the Diffusion object is scheduled for replacement; every other object is inherited.
Original sources, libraries and variant manifests are unchanged.

Default invocation prepares source, patch, object response and compile/link
commands only. `--build-library` explicitly compiles and links a new immutable
output (4 GiB address-space / 180 s CPU per child, one thread). Native compilation
must use the root's scheduled slot. The directory must not already exist. The
prepared `whole_body_integrity_regional_skin_v1` and `whole_body_integrity_regional_skin_prepared_v2` artifacts on the earlier signed parent were source preparation
only; it is not a loadable variant. Production integration must require a completed
manifest and verified executable fixture, and must preserve the parent's signed
header identity when extending the lineage further.

## Executable native evidence

Run `.venv/bin/python scripts/verify_native_regional_skin.py` for source/parent pin
checks. With a root-granted compile slot, add `--compile-and-run`. The verifier
links the immutable signed parent, checks actual `ldd` resolution, caps the fixture
at 1 GiB address space / 120 s CPU, and never initializes or advances a patient.
It uses source-shipped substance definitions through isolated runtime symlinks.

Retained successful evidence:
`data/derived/audits/native-regional-skin-2__axp91/verification.json`.

| Check | Native result |
| --- | --- |
| Regional ownership | Three native leaves, parent aggregate-only, Albumin/Glucose/Sodium initial masses conserved |
| Fluid allocation | Baseline/current/next volume sums 100 mL; baseline compliance sum 2 mL/mmHg |
| Zero-load parity | 10 native circuit steps of 0.02 s with a changing oncotic pressure source |
| Maximum aggregate volume difference | 1.4211e-14 mL |
| Maximum aggregate drainage-flow difference | 9.0414e-14 mL/s |
| Local 133.3224 Pa load | Regional pressure difference 0.9998002 mmHg |
| Local drainage response, normalized by allocation fraction | 0.0099980 mL/s regional difference |
| Unload | Requested region's applied pressure returned to zero |
| Compile resource | 5.70 s, 545908 KiB maximum RSS |
| Native fixture resource | 0.13 s, 43364 KiB maximum RSS |

The initial fixture construction failure exposed a missing substance runtime
working directory: BioGears setup returns before creating its circuit manager if
substances cannot load. The verifier now supplies the held runtime asset symlinks
before constructing BioGears. No source-law or solver adjustment concealed that
failure. Native protein transport and full patient zero-boundary parity remain
separate required gates before activation.
