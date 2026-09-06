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
added at `/home/brandonin/.codex/agents/chore.toml`. At configuration time the old session still exposed four total slots.
After the requested restart, the host explicitly exposed 17 total slots;
that capacity is now confirmed. Do not
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

## Continuing implementation wave

| Owner | Current bounded responsibility | Resource lane |
| --- | --- | --- |
| Root | Shared runtime, signed adapter, live intake/API integration, integration acceptance | Heavy queue coordinator |
| GI transport | Corrected native shared-donor variant and linked conservation fixtures (completed) | Released |
| Signed native port | Muscle-only chemical/heat/work variant and native local ledgers (completed) | Released |
| Mechanical energy ledger | Signed work/heat/reference checkpoints and projection (completed) | Released |
| Support equilibrium | Anatomical surface foundation and constrained supported initialization | Source work; request native slot |
| Intake scheduler / live wiring | Scheduler and live monitor/controller integration (completed) | No production viewer promotion yet |
| Neural coverage | Actual cutaneous operator to shared controller integration | Bounded operator fixtures |
| Vascular evidence | Explicit human-muscle geometry/radius/flow materializations | Bounded graph fixtures |

The same workers can receive subsequent bounded implementation work. Completed
local checks do not close whole-body objectives. The signed adapter's small
boundary smoke passed, but actual mechanical reference coupling and supine
support remain unresolved. The prior static optimizer found an unsupported
free-fall candidate; it was rejected, not promoted.

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
