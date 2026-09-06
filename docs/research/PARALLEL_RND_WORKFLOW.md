# Parallel R&D workflow

The scope remains [the whole-body objective register](WHOLE_BODY_REMAINING_CONCERNS.md).
This workflow schedules work; it does not relax any acceptance requirement.

## Capacity and delegation

After the user-requested restart, this session exposes 17 active slots total: 16 workers plus root.
Children can delegate, but descendants share that limit. Do not create a tree
of waiting coordinators. Prefer independent workers with bounded ownership and a root that also
performs integration. Reassign a completed worker immediately when an independent
ready task exists. Coordinate slot assignments with root before nested delegation.

Configuration follow-up: the installed Codex CLI 0.153.4 loads the documented
`agents.max_concurrent_threads_per_session = 16` setting successfully (`codex
doctor --json`, config.load=ok). The user config was backed up to
`/home/brandonin/.codex/config.toml.bak`; a high-effort bounded `chore` role was
added at `/home/brandonin/.codex/agents/chore.toml`. This does not demonstrate
hot-reconfiguration of this already-running session: its exposed capacity
remains four total slots until the host supplies a different limit. Do not
bypass that limit with extra CLI sessions. Current official guidance is
https://learn.chatgpt.com/docs/agent-configuration/subagents.

The increase from 8 to 16 was requested explicitly; its pre-change backup is
`/home/brandonin/.codex/config.toml.bak.20260905-194117`.

## Task contract

Each assignment specifies the objective, exact owned files, read-only dependencies,
input/interface contract, smallest meaningful acceptance check, resource budget,
and expected handoff. Shared files have one writer. A worker proposes changes to
another owner's interface before editing it. Commit only owned files; preserve
other agents' and the user's uncommitted work.

Separate acquired data, implemented code, isolated validation, integrated
validation, and viewer availability in every report. A successful component
test does not close a whole-body objective.

## Two scheduling queues

Research, source inspection, implementation and small bounded checks run in
parallel. Native builds/runs, bulk processing and full-asset browser work enter
one heavy-job queue. Root grants one explicit heavy-job slot at a time. A grant
is conditional on adequate current machine headroom; defer rather than risk
IBM-1 or the desktop. A per-process cap is not an aggregate machine budget.

For current thin native checks: single numerical threads, nice 10, native/build
address-space cap 4 GiB, and caller cap 1 GiB where compatible. Record commands,
input hashes, PID, output receipts, measured peak memory when available, and
whether validation actually completed. Larger jobs require a fresh measured
budget; do not silently raise a failed cap. No concurrent heavy browser check.

Freeze native source/header inputs from compilation through verification. Before
restarting after interruption, check surviving owned processes and existing
receipts. See [the incident note](INCIDENT_2026-09-05_RESOURCE_CONTENTION.md).

## Current wave

| Owner | Bounded deliverable | Writes | Resource lane |
| --- | --- | --- | --- |
| Root | Coordination, runtime integration and acceptance decisions | Workflow document; root runtime/API files | Light until mechanical slot returns |
| clothing | Compile and verify 92-muscle mechanics/metabolic adapter | Existing owned mechanical adapter files and receipts | Exclusive heavy slot; sequential build then tiny test |
| gi_integrity | Discriminating next experiments for retained GI failures | GI_FAILURE_NEXT_EXPERIMENT.md | Read-only evidence review |
| locomotion_implementation | Adversarial runtime/API/UI integration review | EMBODIED_ACCEPTANCE_REVIEW.md | Read-only review and bounded small checks |

After each handoff, root reads the evidence, resolves interface issues, and queues
the smallest integrated check. The next ready tasks should come from mechanics,
regional blood/lymph/skin, digestion, neural coverage and viewer acceptance;
do not keep polishing a finished local example while a higher-priority lane is
unstaffed. Research can proceed while a native verification occupies the heavy
slot. A production UI build follows actual native acceptance.

## Completion and recovery

Each handoff names the commit or changed files, exact checks and receipts,
remaining limitations, and the next dependency it unblocks. Root owns the
integrated status and updates the objective register only with matching evidence.
Checkpoint progress frequently. Never mark uncertain commands as unapplied or
rerun them automatically. Do not alter concurrent IBM-1 work.
