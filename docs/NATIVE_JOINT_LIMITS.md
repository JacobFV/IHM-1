# the model declares joint ranges and the plant does not enforce them

Every rotational coordinate in `data/models/engineering_stance_v1/model.osim`
carries a `<range>` and `<clamped>true</clamped>`. The model holds **zero**
`CoordinateLimitForce`, and OpenSim does not apply coordinate clamping during
forward dynamics. So the ranges are declared and nothing keeps the body inside
them.

The only joint stops in the plant are the exponential terms of the source
`ExpressionBasedCoordinateForceSet`, and they do not cover the ranges:

| coordinate | passive force | stop |
|---|---|---|
| `ankle_angle_{l,r}` | `-0.1*qdot` | **none — damping only** |
| `pro_sup_{l,r}` | *absent from the file* | **none at all** |
| `hip_rotation_{l,r}` | `∓0.03*exp(14.94*(±q-0.92))` | at ±0.92 rad, against a declared ±0.698 |
| `knee_angle_{l,r}` | `6.09*exp(33.94*(-q-0.13))…` | at −0.13 rad, against a declared 0 |

Two consequences, and the second is the one that cost a result.

**The declared range and the model's own passive stop disagree.** A body lying
prone and doing nothing — 0.02 tonic excitation, 2 s — already sits **0.240 rad**
outside its declared range at `hip_rotation_l`, because the passive stop is
centred 0.222 rad wider than the range. No guard set at the declared range can
be used, because a body at rest violates it.

**A search will buy travel by leaving the ranges, because nothing stops it.**
The first prone crawl parameter search produced 16.0 s of sustained locomotion
and 973 mm of forward travel, completing its horizon with no stall. Measured
against the declared ranges afterwards:

| coordinate | spanned | declared | past |
|---|---|---|---|
| `ankle_angle_r` | −2.524 … 0.000 | ±0.873 | **1.651 rad (95°)** |
| `ankle_angle_l` | −2.514 … 0.000 | ±0.873 | 1.641 rad |
| `pro_sup_r` | −1.590 … 1.281 | 0 … 2.090 | 1.590 rad |
| `pro_sup_l` | −1.474 … 1.349 | 0 … 2.090 | 1.474 rad |
| `hip_rotation_l` | −1.186 … 0.515 | ±0.698 | 0.488 rad |

145° of ankle plantarflexion: the feet folded back on themselves, paddling. **That
result is withdrawn.** It was not a human body crawling; it was a body with no
ankles crawling.

## What was done

`NativeMechanicalStream` takes an optional `coordinate_limits=` list and the
engine installs a `CoordinateLimitForce` per entry. It is **opt-in**, because a
stop changes the plant and the identified stance linearization
(`linearization.npz`, and therefore everything in `scripts/walk_gait.py`) was
solved without one. `scripts/crawl.py` passes the model's own declared ranges;
`walk_gait.py` does not, and its plant is unchanged.

The *limit* is the model's own. The stiffness, damping and transition width are
explicit engineering constants stated by the caller, not ligament measurements,
and they are written into `execution.json` on every run.

Swept over 3 s of the seed pattern at identical port gains:

| stiffness | s per 10 ms advance | worst excursion past the declared range |
|---|---|---|
| none | 0.515 | 1.450 rad (`ankle_angle_r`) |
| 300 N·m/rad | 0.937 | 0.139 rad |
| 100 N·m/rad | 0.505 | 0.192 rad |
| **30 N·m/rad** | **0.264** | **0.220 rad** |

A stop at 30 N·m/rad is **half the wall clock of no stop at all** and holds the
plant 6.6× closer to its declared range. That is not a trade: a plant kept out of
absurd configurations is a plant the error controller can integrate. Only 300,
which is genuinely stiff, costs anything.

The first version of this sweep compared a no-stop row against stopped rows at
*different port gains* and read as though the stops cost wall clock. Same shape
as every entry in the IBM-1 log: a quantity computed correctly and compared
against the wrong thing.

## What is still not enforced

The stop tolerance in `crawl.diverged()` is 0.35 rad, not zero, and that is not
a safety margin — it is the measured amount by which the source model's own
passive stops disagree with its declared ranges. Until those two agree, "inside
the declared range" is not a condition this body can satisfy, and every crawl
report carries `worst_excursion_past_declared_range_rad` so the number is never
implicit.

## Recorded human walking also violates the declared knee range, and by 5°

Measured 2026-09-11 by `scripts/measure_recorded_motion_admissibility.py`, against
gates fixed in IBM-1 `docs/LOG.md` before the script was written. 68 recorded
motions in `data/derived/pose-corpus/` — walking, running, jumping, crouching —
already in OpenSim coordinate space with sha256 provenance to their `.mot`.

**The result: 20 of 42 scorable motions are admissible, 47.6%, against a
pre-registered 90% bar. FAIL.** Real humans produced these trajectories.

Two confounds had to be removed first, and each alone would have manufactured a
total violation out of nothing:

- **Units.** `runningModel_Kinematics_q` holds *degrees* in the JSON — knee down
  to −114.0 — while its siblings hold radians. Each file's `in_degrees` flag
  describes the source `.mot`, not the JSON. Detected by "a radian hinge angle
  cannot exceed 2π", which does not consult the ranges under test.
- **Sign.** The gait2392 family declares `knee_angle_r ∈ [−120°, +10°]`; this body
  declares `[0°, +140°]`. Mirror images — the two disagree on which direction of
  knee rotation is positive. Uncorrected this reads as a 100% violation, and it
  did: the first run scored 4.8%.

That second one is worth keeping, because the gate designed to catch exactly this
class **could not see it**. G1 checks each motion against the model it was
generated from — but that model shares the motion's convention, so the flip
cancels on both sides. *A gate that compares like with like is blind to a
difference between the two likes.* The conventions are now mapped from the two
models' declarations (the sign of the larger-magnitude bound), never from what
makes the data fit.

### What survives: one bound, off by 10°

After both corrections, **every residual violation is below a lower bound and not
one is above an upper bound.**

| coordinate | motions violating | worst excursion |
|---|---|---|
| `knee_angle_r` | 25/42 | **−4.92°** |
| `knee_angle_l` | 24/42 | **−5.19°** |
| `hip_flexion_r` | 8/42 | −0.56° |

The flexion side is clean — this body allows 140° where gait2392 allows 120, and
there are zero violations above any upper bound. The whole failure is that **this
body declares the knee's lower bound at exactly 0, permitting no hyperextension at
all**, and normal human walking reaches about 5° of it at terminal stance.

Three independent sources disagree with that 0:

| source | permitted knee hyperextension |
|---|---|
| this model's own passive stop, `ExpressionBasedCoordinateForceSet` | **7.4°** (−0.13 rad, the table above) |
| `gait2392`, which produced these motions | **10°** |
| recorded human walking, measured here | **uses 5.2°** |
| `engineering_stance_v1` declared `<range>` | **0°** |

The declared range is the only one of the four that says zero, and the model's own
plant already contradicts it — which is the same disagreement this document opens
with, now with a second witness and a number attached.

**Not claimed.** That moving the bound to −10° fixes the corpus. A diagnostic sweep
of the lower bound is indicative only (it is computed with different accounting
than the gate and reproduces 38% where the gate reports 47.6%); it puts
admissibility at roughly 95% by −10°, so nearly all of the failure does rest on
this one number. But the pre-registered G2 threshold is **not** moved, and the
recorded verdict stands at FAIL until the bound is re-derived against measurement
rather than adjusted to pass.
