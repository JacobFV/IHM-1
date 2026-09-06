# 2026-09-05 freezes, interruption and parallel-work conflicts

Purpose: retain causes and operating lessons so repeated failures do not consume
R&D time. This separates host failures from ordinary failed numerical probes and
source/build coordination mistakes.

## Confirmed host evidence

At18:39 PDT, kernel journal inspection found repeated **global out-of-memory
events**, not merely a slow simulation:

| Time PDT | Killed Python PID | Anonymous resident memory, KiB |
| --- | --- | --- |
| 15:48:10 | 3220274 | 114626652 |
| 16:09:04 | 3245778 | 116463080 |
| 16:15:07 | 3253089 | 109716320 |
| 16:33:49 | 3315825 | 99855196 |
| 16:42:04 | 3322806 | 107208768 |

These processes individually held approximately95–111 GiB of anonymous resident
memory. The kernel marked the events `global_oom`. System service watchdog
timeouts also appeared around16:15 and16:42, consistent with severe host-wide
resource pressure. The journal names the allocation requester separately from
the killed process; `nxrunner`, `systemd` or NetworkManager invoking the OOM
killer does not mean they consumed the missing100 GiB.

The filtered log is retained at
`artifacts/incidents/2026-09-05-resource-contention/kernel-excerpt.txt`, SHA256
`59c7de029450a6a71ec444e918284a10a6f7ff80dec66400153cbe53c797ef9a`.

## Attribution limits

The exited PIDs' exact commands and working directories have not been recovered.
They shared the logged application scope, but that alone does not establish
whether a particular process belonged to IBM-1, IHM-1 or another task. Do not
attribute them to a project without evidence.

The user later reported another crash/interruption. These earlier OOM events
explain documented host freezes, but no event has yet been time-matched to that
specific later interruption. At inspection, uptime was1 day13 minutes: the host
had not rebooted during this work period. A browser/app/session crash and a host
reboot are different events.

Earlier IHM work also contributed avoidable concurrent load: full-asset browser
checks used substantial CPU, and elastic hair guide updates were too expensive
under contention. These are observed performance concerns, not proven causes of
the logged OOM kills. Reducing frame rate and using nice does not cap memory.

## Separate development conflicts

1. **Source/build race.** An agent changed `native_tissue_ports.h` after the thin
   adapter compiled. The source hash guard correctly rejected the executable,
   costing a rebuild. Fix: explicitly freeze shared headers with their hashes,
   serialize the build, and prohibit edits until its verification finishes.
2. **Native scalar ownership violation.** The first compression probe aborted
   with `Scalar is marked read-only`. The helper copied native current scalar
   state through ordinary setters. Fix: use the same per-scalar `Override`
   protocol as the native circuit solver; do not disable global read-only rules.
   The corrected0.4-second compression probe passed. This isolated process
   exception is not evidence of a host OOM or desktop crash.
3. **Zero-boundary numerical perturbation.** Rewriting an unchanged pressure in
   Pa caused a tiny unit-round-trip difference and failed exact source parity.
   Fix: a zero external boundary leaves the original native scalar untouched.

## Operating rules for the continuing work

- Keep independent research/code work parallel. Serialize memory-heavy native
  runs, compilation, bulk geometry conversion and full-asset browser work on
  this host; avoid launching another job merely because an agent is free.
- Record expected and peak resident memory for substantial jobs, their PID,
  process tree, project and receipt path. Preserve enough machine headroom for
  the desktop and concurrent IBM work; a successful small test is not a memory
  budget for its scaled-up version.
- Before increasing problem size, use measured memory scaling and an explicit
  bound. Use per-job memory/CPU limits where compatible. **A hard memory limit
  has not yet been installed by this incident note**; nice and single-thread
  BLAS settings alone do not prevent OOM.
- Keep IHM rendering capped, hidden-tab rendering/simulation paused, geometry
  loading bounded and hair dynamics opt-in until a measured affordable path
  exists. Preserve source resolution independently of visualization cost.
- Do not restart uncertain jobs after a crash until surviving owned processes
  and retained receipts are checked. The latest compression verification had
  already completed successfully, so restarting it would have wasted work.
- Keep immutable input/build/command receipts and frequent source checkpoints.
  Close only owned processes; do not stop or reconfigure IBM-1 training.
- Share this incident evidence across project scheduling decisions. Limiting
  IHM alone cannot prevent a different unbounded process from exhausting RAM.

The immediate lesson is to budget **memory and simultaneous workload**, not just
CPU thread count, while separately enforcing shared-file freeze/build windows.

## Later protocol interruption: preserve the physical error

The static-pose stream wrote its JSON marker before evaluating the requested
pose. A native diagnostic could therefore interrupt the protocol record. The
original malformed line was not retained, so its exact contents are unknown.
A replay of the retained pending coordinate request after the framing fix
returned a clean material-domain rejection; it did not reproduce a host crash.
See `SUPINE_INITIALIZATION.md` and receipt
`data/derived/pending-static-pose-rssrue1z` for the bounded replay.

The fix builds the complete response before emitting its marker (`ac4ddcd`).
Pending commands, successful responses and rejected trials are now retained,
and continuing native state is checked after rejected evaluations. Reusing
cached responses across that rebuild required an explicit framing-only source
comparison; physical headers, libraries and inputs still had to match.

Operational lesson: do not rerun a long optimization merely because transport
failed. Retain the exact pending request, distinguish protocol failure from a
physical model rejection, replay the smallest failing case, and preserve source
identity before resuming cached work. This incident is separate from the earlier
resource-contention evidence.
