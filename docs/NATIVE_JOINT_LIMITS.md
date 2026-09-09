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
