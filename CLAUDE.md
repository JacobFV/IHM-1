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

## Call it twice at the same input

**A query that is meant to be a pure function of position must return the same answer when called
twice at the same position.** `LocalAssociation` did not: called a second time at *identical* node
positions it moved **141 of 3,123** associations by more than 1 mm, up to **6.5 mm**, then stabilised
from the third call on. That one-time disagreement, re-imposed by a cut-back loop's discarded
attempts, produced a step-size-invariant count of ~42 inverted elements that was read as a property
of the mechanics — of the mesh, the stepping, the association concept, the starting state — for a
full round of work before anyone called the function twice.

The check costs one line and it is not the same as a known answer: a known answer tests the value,
idempotence tests whether the instrument is a function at all. Run it on anything stateful, cached,
or seeded from a previous result.

**The cause was the third instance of one hazard on this line: a coarse proxy standing in for an
exact query on a mesh with a long facet tail.** The local search seeded its starting face by
*centroid proximity*, and this bed's faces run 1.6 mm median against a 28.9 mm maximum, so the
nearest centroid is routinely not the face the association sits on. Before it, the unsigned ray
index and the radius test failed the same way. **When a mesh has a long facet tail, a proxy that is
right for the median face is wrong for the tail, and the tail is where contact lives.**

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

## Gate the part you derived, measure the part you guessed

A scaling argument usually fixes an **exponent** and not a **prefactor**. The sagitta bound
`d²/2RN` for association staleness predicted a 1/R scaling; the measurement gave `(1/R)^1.068`,
which matches, at a constant **0.24 × the bound** — four times smaller, at every radius, because
the solve carries part of the offset.

The gate had been written as *ratio ∈ [0.5, 2.0]*, which folds both claims into one number, so it
failed on the prefactor while the derived half was confirmed. **A single ratio test cannot pass the
derived half and fail the guessed half; it only ever reports the guess.**

So: gate the exponent, and *measure and report* the prefactor rather than predicting it. If you
cannot say where a constant in your bar comes from, it does not belong in the bar.

## A sum over repeated events cannot tell drift from a rejected jump

An accumulator measuring "how stale is this constraint" added up every *refused* association update.
It returned a median of 3.08 mm and a maximum of **143.13 mm** — the "catastrophic" branch of a
pre-registered decision. Then: **143.1322 / 47.8232 = 2.99.** It was one ~47 mm teleport, refused
three times, summed as though the constraint had drifted 143 mm.

A jump limit exists to reject teleports. **Counting a successful rejection as accumulated error
credits the mechanism's correct behaviour to its failure mode.** The fix was the *instantaneous*
wanted displacement, with the accumulated figure kept beside it because the RATIO is the diagnostic:
accumulation ≫ instantaneous means the same jump rejected repeatedly; accumulation ≈ instantaneous ×
count means real drift.

**Two things about how it nearly got through, and neither is covered by the other rules here.**

The measurement was pre-registered — boundary, three figures, accounting decisions, all declared in
the file before the run. The hypothesis it appeared to confirm was pre-registered too, with a
mechanism. Everything was in the right order. **Pre-registration protects against choosing a
threshold to fit a result; it does nothing about a quantity that is ill-defined for the question.**

And it was caught only because a ratio was suspiciously close to an integer. The median alone, 0.82
→ 3.08 over two steps, looked like exactly the ratchet that had been predicted. **A result that
agrees with the hypothesis you pre-registered is the one you are least likely to take apart** — so
confirmation is the moment to check the instrument, not the moment to stop.
