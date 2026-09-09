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

## Jobs

- Never commit weights, caches, meshes or large payloads. `data/derived/`,
  `artifacts/` and `logs/` are gitignored; commit scripts, schemas and reports.
- `pkill -f <pattern>` matches your own shell's command line. Kill by PID.
  Orphaned pool workers have been a recurring problem here — clean up.
- The machine is shared with several agents and with IBM-1 training. Be frugal
  with worker counts.
- Branch is `feat/integrated-human`. Commit and push as you go.
