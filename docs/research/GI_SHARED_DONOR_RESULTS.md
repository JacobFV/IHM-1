# Shared-extracellular-donor transport probe

2026-09-05. The held facilitated-diffusion source reproduces a local availability defect in a three-compartment native fixture. Two independently capped outflows can consume the same extracellular starting pool twice. Native leaf and parent compartment layouts fail differently. This establishes neither occurrence nor causal contribution during the retained exertion_v3 trajectory.

## Native execution and retained evidence

Run `python3 scripts/verify_gi_shared_donor.py` only while holding the shared heavy execution slot. The script compiles one small executable and runs eight independent one-call fixtures. It loads substance definitions, creates native compartments, and calls facilitated diffusion once per process. It never loads, initializes or advances a whole-body patient state.

The completed receipt is `data/derived/audits/gi-shared-donor-gepygiov/report.json`: **21/21 diagnostic checks passed**, including checks that the defects are reproduced. This is a successful reproduction gate, not a physiological acceptance gate. The receipt retains compiler command, exact source and executable hashes, resolved dependency hashes, each fixture's stdout/stderr, and `/usr/bin/time -v` resource reports. Two earlier fixture compile errors are retained in `gi-shared-donor-8oq2035c` and `gi-shared-donor-4yaqlke8`; they concerned protected construction APIs, not numerical source behavior. The final fixture registers glucose through the native public compartment-manager interface.

The probe extracts these three complete method bodies verbatim from held `Diffusion.cpp` and compiles them into the executable:

- `CalculateFacilitatedDiffusion`
- `DistributeMassbyMassWeighted`
- `DistributeMassbyVolumeWeighted`

The entire source must match SHA-256 `9f90e429a87ae1ba3b7071fc14525905e768d23a27a0f21d15a2581a49b7d48b` before compilation. This explicitly binds the methods executed to the inspected source without assuming a historical shared library contains those bodies. Native BioGears compartment classes and native CDM `Balance` remain linked from the dependency inventory. No production source or header is edited. The selected supporting BioGears library is the substrate-availability variant; its and every resolved dependency's hashes are in the receipt. This is source-executed native method evidence, not execution of a historical full-body binary.

Ending masses and concentrations are printed after the actual native call. Raw/capped VE/EI amounts in the JSON are **independently mirrored algebra**, not injected native counters. The report labels them `mirrored_algebra_not_native_counters`. Their incidence ledger is `[-VE, VE−EI, EI]`; final native mass changes are compared against that ledger, with the held parent distribution helper's extra debit cap accounted for separately. The native end states directly establish the reported failures regardless of that mirrored calculation.

## Inputs and outcomes

All three compartments have 1 mL total volume. For the parent case, E has two native leaf children, each with 0.5 mL and half the starting E mass; parent and children are not double-counted. The held Glucose.xml is loaded and its hash checked against the source definition. Native Km is 0.8, an untyped scalar used with concentrations in **g/mL** by this method. Native maximum diffusion flux is 1 g/cm²/s. Every fixture explicitly supplies a combined coefficient of 1 g/s; no tissue-area calibration or tissue mass is inferred.

| Fixture | Initial V/E/I, µg | dt, s | Native final V/E/I, µg | Total residual, µg |
| --- | --- | --- | --- | --- |
| Zero gradient | 1000 / 1000 / 1000 | 0.02 | 1000 / 1000 / 1000 | 0 |
| Small positive gradients | 1100 / 1000 / 900 | 0.02 | 1097.500312461 / 1000 / 902.499687539 | 0 |
| Matched negative gradients | 900 / 1000 / 1100 | 0.02 | 902.500312539 / 1000 / 1097.499687461 | 0 |
| E above both neighbors, ample | 0 / 1000 / 0 | 0.02 | 25.031289111 / 949.999921875 / 24.968789014 | 0 |
| Shared scarce E, leaf | 0 / 1000 / 0 | 1 | 1000 / **−1000** / 1000 | 0 |
| Shared scarce E, parent | 0 / 1000 / 0 | 1 | 1000 / 0 / 1000 | **+1000** |
| Adversarial gradient exactly −Km | 0 / 800000 / 800000 | 0.02 | 800000 / 0 / 800000 | 0 |
| Adversarial gradient below −Km | 1000 / 900000 / 900000 | 0.02 | 0 / 901000 / 900000 | 0 |

The scarce cases deliberately increase dt to 1 s while keeping the same masses and coefficient as the ample case. This is an adversarial availability fixture, not evidence that the retained engine's timestep/area combination produced it. The negative-gradient singularity fixtures use extraordinary concentrations and are explicitly outside the retained physiological range.

## Mechanism and interpretation

In the scarce cases, mirrored uncapped VE is −1251.564455569 µg and EI is +1248.439450687 µg. Each independent cap permits 1000 µg from the same initial E donor. The recipients therefore receive 2000 µg from an available 1000 µg.

For a leaf E, the native method writes the combined −2000 µg debit directly, leaving −1000 µg. Its native final concentration is **zero**, because CDM concentration calculation locally substitutes zero for a negative input mass. It does not repair the negative mass scalar. A nonnegative concentration snapshot can therefore conceal this negative store.

For a parent E, `DistributeMassbyMassWeighted` caps the combined debit to the parent's available 1000 µg. Both recipient credits remain unchanged. The donor stays nonnegative while total glucose increases by 1000 µg. Merely clamping donor mass to zero cannot solve this defect.

Matched small positive/negative gradients conserve mass but have slightly different flux magnitudes, consistent with the signed denominator in the held law. At exactly −Km, the mirrored raw result is negative infinity; the native call returns finite final masses because its availability cap limits the transfer. This fixture does **not** establish finite-input rejection. Below −Km, native V loses glucose even though its concentration is below E: the denominator changes sign and produces transport against the initial gradient. These adversarial results do not diagnose physiological-range delivery.

## Proposed correction and next evidence gate

The proposed shared-donor correction is a single conservative transfer allocation before any compartment mutation. After existing individual donor caps, collect all simultaneous outflows from each starting donor. For shared E, compute `W = max(-VE, 0) + max(EI, 0)`. If `W > E_available`, scale those E outflows together by `E_available/W`, then derive every recipient credit and donor debit from the same final transfer pair. With proportional allocation, the scarce fixture would end at V/E/I = **500/0/500 µg** in both layouts. This policy preserves the relative requested allocation and must be reviewed as a model choice; sequential allocation would be order-dependent and different. No production correction is implemented in this task.

Do not count simultaneous incoming mass as starting availability without an explicit integration policy. Do not repair conservation with independent recipient credits and a later donor clamp. Before a fix ships, require finite/nonnegative pools, total and per-owner transfer closure, unchanged zero/ample controls, and identical leaf/parent aggregate outcomes. Then add event counters for actual native muscle V/E/I gradients, timestep, coefficient and competing E withdrawals around glucose collapse. Only those observations can establish whether this mechanism contributed to exertion_v3. A signed-law or Km-unit correction is a separate design requiring parameter evidence.

## Resource bounds

Each compiler/fixture child had a 1 GiB address-space limit, 120 s CPU limit, 150 s wall timeout, and single-thread BLAS/OpenMP settings. The successful compiler took 1.88 s wall time with peak RSS **516,812 KiB**. Each fixture took 0.05–0.08 s with peak RSS at most **40,500 KiB**. No heavy execution remained active when the slot was released. These tiny fixture measurements do not justify a larger whole-body job.
