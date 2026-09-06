# Explicit named native Skin configuration: prepared, not activated

`regional_skin_configuration.py` binds two caller-selected unions of engineered
exterior skin masks and the entire exterior complement to three semantic region
names. It does not relabel the accepted `.2/.3/.5` installation. No shared native
header, adapter, Python observer, or canonical geometry is modified.

The retained example selection uses all four engineered left lower-leg masks,
all four right lower-leg masks, and the remaining exterior proxy. The names are
`engineered_left_lower_leg`, `engineered_right_lower_leg`, and
`engineered_exterior_complement`. Its three explicit, uniform priors are 0.001 m
thickness and 0.2 extracellular fraction. These numbers are **engineering weighting
parameters**, not measured tissue properties. Their product cancels during this
uniform normalization; callers can provide separately evidenced nonuniform priors.

For each region, the weight is proportional to `exterior proxy area * thickness *
extracellular fraction`. Weights are positive and normalized, with the third obtained
by residual correction. This example produces approximately
`0.052809544820190284 / 0.0528907236453193 / 0.8942997315344904`.
Every existing native Skin extensive quantity is partitioned once across the three
stores. The proxy-volume sum does **not** replace native Skin volume. Excluding
inner/seam faces from contact excludes duplicate surface representations, not native
fluid, protein, unlocalized tissue or the remainder's inventory.

## Receipts and runtime mapping

`data/research/configured_regional_skin/example_selection.json` is the explicit
selection/prior input. `prepared_v1/configuration.json` contains source geometry
receipts, the materialization identity, all three named prior records, normalized
weights and `face_to_native_index` for all 203,382 source faces. Values 0/1/2 identify
the named configured owner; −1 excludes inner/seam contact representations. Exactly
109,183 exterior-proxy faces are eligible. The open-shell exterior assumption from
the materialization remains explicit; no physical surface exclusivity or lymph
watershed validation is claimed.

The configuration's canonical-content SHA-256 is
`ad93aa3dcc4e9bcfeb54f076c7bbe76d316f0b2e6d52ec4865536a20e511aefb`.
It covers the mapping, source receipts and all weighting assumptions. The staged
configured header SHA-256 is
`7401a91d4d1c705334db6a65a4bd9f7ff1bc3bb8d0b0ecebe8b0a99096e50413`.

`map_contacts(config, contacts, root=...)` resolves each exact source-face ID to the
configured semantic name. It rereads and verifies the pinned geometry, checks each
patch against its source triangle area, rejects repeated faces and inner/seam faces,
and preserves supplied vector force and normal load. Whole-region surface-pressure
proposals use that region's entire eligible area, including unloaded portions.
Those proposals are not identified interstitial pressures and produce no native
commands. A future opt-in adapter must explicitly choose a pressure-transmission
prior and use the same configuration digest; no anonymous `region_a` assignment is
performed here.

## Exact patch boundary and scaling laws

`prepare_configured_regional_skin.py` reads the immutable accepted
`whole_body_integrity_regional_skin_graph_v2` header and produces a separate
`NativeConfiguredRegionalSkin` class. The guarded patch changes only class name,
three fraction constants, three semantic names and a configuration digest member.
Reversing those edits must restore the accepted header byte-for-byte. The runtime
configuration is immutable and compiled into this opt-in header; the matching JSON
supplies face-to-owner mapping. No library build or production variant is claimed
by this source preparation.

The inherited algorithm preserves these laws for fraction `f_i`:

| Native quantity/law | Configured branch |
| --- | --- |
| Node baseline/current/next volume | Extensive partition with last-branch residual |
| Species current mass | Same extensive partition; native read-only aggregate parent |
| Baseline/current/next compliance | Extensive partition |
| Baseline/current/next resistance | `R / f_i`, preserving summed parallel conductance |
| Baseline/current/next flow sources, including sweating | Extensive partition |
| Current/next path flow snapshots | Extensive partition |
| Current/next node pressure and pressure-source states | Copy unchanged |
| Valve state | Copy native state |
| Regional mechanical drive | Zero initially; independent explicit local pressure |
| Aggregate native Tissue/Diffusion chemistry and lymph destination | Retained shared laws/owners |

Global leaf caches, transporter graph replacement, detached native law-cache updates,
parent read-only aggregation and inherited owning-leaf sweat/waste handling remain
unchanged. Source-law preservation is not sufficient to claim nonlinear engine parity
for different weights. It only establishes the algebraic total inventory, compliance,
conductance and flow-source identities.

## Native fixture preparation and gates

The parent manifest is pinned to
`52de403e3dab2682f774755b7bf1370711bc5ab218244c581460b84a82aad130`,
with the accepted graph-v2 library
`bc91cbab829c04bfa2df7490433af5bee752e7514332ac715e9f0bd11a17650a`.
`preparation.json` embeds its full parent manifest/ancestry and all 360 object hashes.
No unaccepted cardiovascular, sleep, or other experimental descendant is combined.

The staged fixture reuses the actual native circuit solver acceptance source,
changes only its header/class reference and replaces the two hardcoded `.2/.3`
specific-flow denominators with configuration weights. It checks installation
baseline/current/next volume, compliance, species ownership and read-only aggregation;
10 matched zero-pressure steps with changing native source laws; local pressure
separation/drainage response; and pressure release. It does not initialize a patient
or execute regional protein transport.

- `.venv/bin/python scripts/verify_regional_skin_configuration.py`: four lightweight
  groups pass, including actual-face contacts, exact preparation reproduction,
  algebraic conservation, invalid priors/overlap/tamper rejection and reversible patch.
- `.venv/bin/python scripts/verify_configured_regional_native.py`: source receipts
  verified, no native process executed.
- `--compile-and-run`: explicit queued native gate; one compiler process with 1 GiB
  address space, 60 s CPU and 75 s wall limits, then one ≤10 s circuit fixture.
  Outputs, failures, timing, runtime linkage and receipts are retained in a new audit
  directory. Do not invoke until the root grants the shared heavy slot.

The circuit gate has now passed as recorded below. It must be followed by matched full-patient zero-load and load/release lifecycle checks,
fluid/species/sweat ledgers and configuration-aware observer/adapter integration in
new files before activation. The accepted existing installation remains unchanged.


## Completed isolated native circuit acceptance

The granted serial run completed on 2026-09-06 with no patient initialization or
advancement. `data/derived/audits/configured-regional-native-k0i3zjpb/` retains
compile/run stdout, stderr, resources, linkage, staged sources and verification.
The checked-in `native_circuit_acceptance.json` records their byte hashes and the
same configuration digest `ad93aa3d...`.

| Check | Observed result |
| --- | --- |
| Matched zero-load steps | 10 passed |
| Maximum aggregate volume error | 4.263256414560601e-14 mL |
| Maximum aggregate selected-flow error | 1.7436052601738083e-13 mL/s |
| Local regional pressure difference | 0.9998001855873788 mmHg |
| Local specific-drainage difference | 0.009998001855567137 mL/s |
| Pressure release | Passed native pressure-source check |
| Compiler | 2.04 s, peak RSS 549,240 KiB |
| Native fixture | 0.10 s, peak RSS 43,724 KiB |

Installation inventory and source-law updates passed the inherited fixture checks.
These are floating-point tolerance parity results, not bitwise equality. Regional
protein transport, full-patient initialization/stepping and coupled mechanical
pressure transmission were not exercised. No production activation is enabled.
The process was reaped and the shared native slot released immediately after the run.
