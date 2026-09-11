# working notes for agents on IHM-1

IHM-1 is the body. IBM-1 is the brain. **Read
`/home/brandonin/Documents/IBM-1/CLAUDE.md` as well** — its metrics discipline,
its randomness traps and its corrections ledger apply here without restatement.

## Read these before moving the body

- **`docs/ACTUATION_STAGES.md`** — how the body comes to move, and which parts
  are scaffold. Read it before touching anything that actuates, contacts or
  poses. It is the single most misreadable thing in this repo.
- `../IBM-1/docs/DIRECTION.md` — the programme target and the standing
  corrections that produced it.
- `../IBM-1/docs/DISCONNECTS.md` — where a declared model and the running model
  are different objects.

## The one thing that is easiest to get wrong

**There are two bodies and the one that runs is the scaffold.**

The engine integrates 22 rigid bodies, 80 muscles, 33 coordinates. That is the
crude body, and it is the only thing that touches the world. The real body —
3,816 anatomical entities — is posed KINEMATICALLY from those segments through
`data/derived/anatomy-segment-binding/` and feels nothing itself.

The crude body has always been a scaffold with a planned disposal. Making it
better *as a scaffold* is worth doing; deepening dependence on it is not. Never
report what the crude body did as what the body did.

## Traps measured here, each of which cost real time

- **Nothing enforces the joint ranges the model declares.** Every rotational
  coordinate carries `<range>` and `<clamped>true</clamped>`, the model holds
  **zero** `CoordinateLimitForce`, and OpenSim does not clamp during forward
  dynamics. A prone search found 973 mm of travel with the ankles folded to 145°.
  Anything free to search will find this. Opt-in stops now exist — use them, and
  check any trajectory against the declared ranges before believing it.
- **`advance()` cost is unbounded by `dt`.** It is an error-controlled Simbody
  integration. Cost ran 0.112 s → 24.3 s per 10 ms step when 26 unconnected
  `CoordinateActuator` ports left the arms a rag doll and a muscle went outside
  the force–velocity domain. Hold every declared actuator port.
- **A reported summary field is not the system.** The gait trajectory's
  `contacts` reports foot load; the plant carries 28 contact elements over 20
  bodies. Concluding "the body has two feet of contact" from the former was wrong
  and cost a redirection.
- **`SegmentalCord` without `muscle_bindings`** matches only bare OpenSim ids and
  silently gives arcs to 66 of 249 channels. Pass the bindings.
- **Files with coordinate-named columns are not all poses.** `_Kinematics_dudt`
  is acceleration and `_Actuation_power` is watts; both parse cleanly and produce
  muscle excursions of 1e22. Gate on physiology (path/optimal in 0.2–20), not on
  filenames alone.
- **Two bodies of different mass.** Mechanical 77.6122029 kg (a bare literal in
  67 places across 58 files, no derivation found) against anatomical 70.7713 kg.

## A pre-registration can carry an unmeasured assumption

**Fixing a threshold before the result protects against one failure. It does nothing about a
parameter adopted by ANALOGY.** `CORR_TRIM = 0.10` -- discard the largest-separation 10% of
correspondences per segment -- entered the v1 skin-warp pre-registration as "the same 10% trim the
per-segment fit itself uses", reasoned from ICP's outlier rejection. It was never measured. It then
rode into v1 spline, v2 flow and v3 anchored as production behaviour.

Measured against a known truth a month later, in a 2x2 that varied it against the correspondence
rule: **the trim costs 22% of recovery accuracy and the correspondence rule is worth ±0.02 mm.**
A day had been spent comparing correspondence rules -- nearest point, normal shooting, a hybrid --
while the term that carried all of the effect sat in a constant nobody was varying.

The guard is not more pre-registration. It is: **when a pre-registration imports a parameter by
analogy, write down that it is unmeasured and what would measure it.** A borrowed constant is an
assumption wearing a pre-registration's clothes, and the discipline that catches post-hoc threshold
moves is blind to it by construction.

## A control that cannot fail has not passed

**Two separate checks, and the second is the one everyone skips.**

*Did the run actually use the configuration it claims?* A gate-V run reported **PASS** while its own
header read `stepping: allornothing` — the driver set `sys.argv = ["x"]` before importing the seat
module, wiping `--adaptive` before the mode was read. The verdict was true of the run and false of
the thing the run was supposed to test. **Read the mode back out of the run's own output; never
infer it from the flag you passed.**

*Can the mechanism fail at all?* On the plane the adaptive rule rejected **0 of 8** steps. A rule
that never engages reports 0 rejections too, so 0 is indistinguishable from a no-op. Forcing the
facet slack to 0 mm made it reject and stall, which is what a live mechanism must do. **A control
whose pass looks identical to its absence is not evidence until you have made it fail on purpose.**

This is adjacent to *a known answer must break the symmetry it is testing*, but distinct: there the
test case was too weak, here the test never ran. Both produce a green result from an instrument that
measured nothing.

## Jobs

- Never commit weights, caches, meshes or large payloads. `data/derived/`,
  `artifacts/` and `logs/` are gitignored; commit scripts, schemas and reports.
- `pkill -f <pattern>` matches your own shell's command line. Kill by PID.
  Orphaned pool workers have been a recurring problem here — clean up.
- The machine is shared with several agents and with IBM-1 training. Be frugal
  with worker counts.
- Branch is `feat/integrated-human`. Commit and push as you go.
- **`git add <file> && git commit` commits the WHOLE INDEX, not your file.** With
  several agents working in this one checkout, another agent's staged files ride
  along silently and land under your commit message. It happened in `4386509`,
  which claims to be a solver finding and carries 53 lines of `TISSUE_MECHANICS.md`
  and 103 lines of `register_knee_cartilage.py` that its author never touched,
  while the actual author's own commit `5c35a13` reads as a 5-line change. Nothing
  was lost and the branch is correct; the attribution is not, and a log that
  misattributes work is a log you cannot use to find out when something changed.
  **Always name the paths on the commit itself** — `git commit -- <paths>` (or
  `--only`) — so the index you did not build cannot follow you in. Do not leave
  files staged between steps either: stage and commit in one action.
