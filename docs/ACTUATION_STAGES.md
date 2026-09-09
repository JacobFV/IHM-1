# how the body comes to move: the staged plan, and what is scaffold

Read this before working on anything that moves the body. It exists because the
low-resolution body is easy to mistake for a component of the design when it is
a **scaffold with a planned disposal**, and because this programme's recorded
failure mode is drifting toward whatever is measurable on the proxy.

The end state, stated first so nothing below is mistaken for it:

> **The brain directly controls individual muscles. Those muscles' contractions
> control the skeleton. The whole thing is a highly complex rigid AND soft body
> system, actuated by the dynamic elasticity of the muscles, which is determined
> by the brain's outputs.**

**Contact with the world is never bone against world.** It is bone, mediated
through the interlayering of soft bodies — fatty tissue, muscle, skin and the
rest — and it is the skin that meets the floor. A bone mesh is a collider against
*other bones and its own soft tissue*, not against the ground.

**And the end state is not a later stage. It is to be fully implemented now.**
The body must be selectable between two modes:

| mode | meaning |
|---|---|
| **driven** | the body is posed from the crude scaffold. Fast, and what generates training corpora. |
| **fully present participant** | the body is dynamically simulated in its own right: muscles actuate, soft tissue deforms and mediates contact, the skeleton moves because muscles pulled it. |

Both exist, and the caller chooses. The stages below describe how work is
sequenced and what the scaffold is *for* — they do not license shipping only the
scaffold and calling the participant mode future work.

---

## The two bodies, and which one is real

| | crude body | real body |
|---|---|---|
| what | 22 rigid bodies, 80 muscles, 33 coordinates | 3,816 anatomical entities |
| driven how | dynamically integrated by the engine | **kinematically posed** from the crude body through the segment binding |
| touches the world | yes — all contact happens here | no |
| status | **scaffold** | the thing being built |

The binding is `data/derived/anatomy-segment-binding/` (99.88% of entities
assigned, gates in `report.json`). `ihm/assembly/body.py` states the real body's
transforms come from the run and never hold a last value.

**The crude body has always been a scaffold.** It is not a simplification we
have settled for; it is a stage. Work that makes the scaffold better *as a
scaffold* is worth doing. Work that deepens our dependence on it is not.

---

## The stages

### 1. Now — the scaffold actuates, and generates training data

The crude body actuates the real body's muscles and bones. It determines rough
contact with the floor and with objects. That motion, and the afference it
produces, is **collected as training data for teaching the brain to actuate the
real body**.

This is why forced motion is legitimate here and only here: a pose trajectory
imposed on the body produces what real muscles and skin experience, and that
corpus bootstraps circuits that would otherwise have nothing to condition on.
`scripts/collect_pose_corpus.py` (68 motions, 17,622 frames) and
`scripts/collect_forced_gait_corpus.py` are that step.

Contact is rough on purpose at this stage. It only has to be good enough to
generate plausible load and afference.

### 2. Possibly — the brain learns to actuate the scaffold too

Letting the brain drive the low-res body may earn its place as a curriculum
stage: a lower-dimensional control problem to solve before the real one. Optional,
and justified only by whether it helps the transfer.

### 3. The scaffold's jobs move to the real body

Not because the participant mode arrives later — it is to be built now — but
because each job the scaffold still holds is a place the real body is not yet
sovereign. Contact first, then actuation. When a job has moved, the scaffold
keeps it only as the `driven` mode's fast path.

### 4. The end state

Brain output → individual muscle activation → muscle contraction → skeleton
motion, in a coupled rigid and soft body system where the muscles' dynamic
elasticity is what does the work.

---

## What this means for decisions you are about to make

**The participant mode needs the real body to be dynamic, and that is the point.**
Earlier guidance here said not to give the anatomical body its own dynamics; that
was written when the end state was mistaken for a later stage, and it is
withdrawn. What remains true is that the `driven` mode must stay cheap — it is
the corpus generator — so the two modes are different code paths over the same
anatomy, not one compromise between them.

**Contact belongs to skin, not bone.** A design that puts the ground collider on
a bone mesh has skipped the fat, muscle and skin that actually meet the floor.
The 21,381-point skin contact quadrature and `surface_contact_manifest` already
exist and are the right place to start.

**Do not optimise the scaffold's fidelity as an end.** Making the 22-segment
body's gait beautiful is not progress toward a brain controlling real muscles.
Making its contact good enough to generate honest afference is.

**Do not report scaffold behaviour as body behaviour.** A trajectory the crude
body executed is a trajectory the crude body executed. The real body was posed
from it. Both are true; only one is what the programme claims to be building.

**Forced motion is never a result.** It is stage-1 data generation. This is
recorded in `docs/DIRECTION.md` in the IBM-1 repo and repeated here because it is
the single easiest thing to misreport.

---

## Where the pieces are

- segment binding — `data/derived/anatomy-segment-binding/`, `docs/ANATOMY_SEGMENT_BINDING.md`
- pose corpora — `scripts/collect_pose_corpus.py`, `scripts/collect_forced_gait_corpus.py`
- the plant and its contact geometry — `scripts/native_mechanical_stream.cpp`
- innervation, so the brain reaches muscle and skin only through nerves —
  `scripts/measure_innervation_coverage.py` in IBM-1
- the standing direction and its corrections — `docs/DIRECTION.md` in IBM-1
- where declared and running models differ — `docs/DISCONNECTS.md` in IBM-1
