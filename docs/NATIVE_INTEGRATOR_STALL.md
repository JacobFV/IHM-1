# why one 10 ms advance stops returning

`TimeoutError: Native mechanical stream response timed out` at 1.83 s of a 5 s
gait horizon. Reproducible to the instant across machine loads, which is what
made it tractable.

## It is not a timeout

The message comes from a 120 s `select` deadline in
`ihm/native/mechanical_stream.py`. It says the engine did not answer; it does not
say why. One `advance` builds a fresh OpenSim `Manager`
(`setIntegratorAccuracy(1e-7)`, max step 5e-4, internal step limit 100000) and
runs an error-controlled Simbody integration over the interval, so **the cost of
an advance is not bounded by `dt`**.

Timed per advance with the deadline raised to 3600 s
(`scripts/diagnose_gait_timeout.py`):

| sim time | wall for one 10 ms advance |
|---|---|
| ≤ 1.75 s | 0.112 s (median over 175 advances) |
| 1.79 s | 0.88 s |
| 1.80 s | 1.24 s |
| 1.81 s | 2.65 s |
| 1.82 s | 24.3 s |
| 1.83 s | had not returned after **22 minutes** |

A factor of ~10 per 10 ms of simulated time. Raising the deadline converts a
fast crash into a slow one and nothing else.

## What collapses the step

The state entering the first slow advance, measured:

```
pro_sup_r        11.03 rad/s
elbow_flex_l      4.97 rad/s
arm_flex_r        4.93 rad/s
arm26_BIClong_l   3.05 m/s of fibre velocity
                  = 23.9 optimal fibre lengths per second,
                  against the muscle's maximum contraction velocity of 10
```

Every one of those is in the arm chain. Past 10 optimal fibre lengths per second
the force–velocity curve is being extrapolated, the equilibrium muscle's
inversion is ill-conditioned, and the error controller has to shrink the step
without bound.

(An earlier version of this note said `pro_sup_r` had bounced off its own q=0
clamp. That was an inference from `<clamped>true</clamped>`, and it is wrong:
nothing enforces a coordinate clamp in this plant. See
`docs/NATIVE_JOINT_LIMITS.md` — measured, `ankle_angle_r` reaches 2.52 rad
against a declared range of ±0.873.)

And the arm chain had no motor drive at all. The source model declares 26
`CoordinateActuator` torque ports — lumbar ×3, and per side shoulder
flexion/adduction/rotation, elbow, forearm pronation — and **no controller was
ever connected to them**, so their controls were identically zero for the whole
run. The arms were a rag doll riding a walking body. `pro_sup` is worse than the
others: the source `ExpressionBasedCoordinateForceSet` damps shoulder, elbow,
hip, knee, ankle and toe, and has no term for forearm pronation at all, so that
coordinate has neither drive nor damping nor any limit whatsoever.

## The fix

The engine grows an explicit port over the actuators the model already declares
(`ExcitationPorts` instance `declared_coordinate_actuator_ports`), reported in
its own `coordinate_actuators` block so nothing can count a torque port as a
muscle. `NativeMechanicalStream.advance` takes `coordinate_actuation=`. With
every command at zero the plant is exactly what it was before the port existed.

`scripts/walk_gait.py` then holds the arms at the registered initial pose with a
PD sized to each segment's own inertia.

**A/B, same binary, search-J best parameters, 5 s horizon:**

| | simulated | wall | stop |
|---|---|---|---|
| ports at zero | 1.83 s | 180 s | `native_failure` — the timeout |
| ports held | 1.49 s | 18.5 s | `pelvis_tilt_beyond_0.55rad` |

With the ports held there is no step-size collapse anywhere: the whole run costs
0.124 s per advance, the same as the baseline's healthy region, and it ends on a
measured mechanical event instead of an engine failure. Travel and the single
genuine step are preserved (160 mm, foot advance 163 mm).

The arm hold changes the plant, so 1.49 s and 1.83 s are not the same
experiment. The claim is only the narrow one: **the integrator no longer
collapses, and the run ends for a reason the mechanics can state.**

## The same failure is not specific to gait

The first prone crawl attempt stalled at 3.88 s with a 110 s advance. Same
signature: `elbow_flex_r` at 10.6 rad/s dragging `arm26_BIClong_r` to **38.1**
optimal fibre lengths per second. A PD port is not a guarantee — it is possible
to command targets the port then chases too fast.

So `scripts/crawl.py` carries a guard: stop when any fibre exceeds 15 optimal
fibre lengths per second, and report in every run how long it spent past the
model's own maximum of 10. That turns an unbounded wall-clock stall into a fast
termination that names the muscle. A run that ends normally still says how far
outside the declared muscle model it went.

**Over-damping the ports is the other cliff, and it is closer than it looks.**
The commands are held at 100 Hz, so a zero-order-held damper is stable only
while `kd·dt/I < 2`. A sweep at 0.6 s of settling, counting commands driven to
their ±1 ceiling out of 780:

| ω (rad/s) | ζ | peak coordinate speed | saturated commands |
|---|---|---|---|
| 8 | 0.7 | 18.7 rad/s | 1 |
| 8 | 1.0 | 18.7 | 2 |
| 14 | 1.0 | 18.6 | 8 |
| 20 | 1.0 | 18.8 | 178 |
| 14 | 2.1 | **47.9** | 392 |
| 20 | 2.1 | **110.3** | 500 |

(The 18.7 rad/s floor is the initial drop transient, not the ports.) `crawl.py`
uses ω = 10 rad/s, ζ = 0.9, which puts every port at `kd·dt/I = 0.18`.

## Reproducing

```
.venv/bin/python scripts/walk_gait.py --mode single --horizon 5.0 \
    --params data/derived/gait-search/searchJ_best_params.json            # holds the arms
.venv/bin/python scripts/walk_gait.py --mode single --horizon 5.0 --no-arm-hold \
    --params data/derived/gait-search/searchJ_best_params.json            # reproduces the stall
IHM_NATIVE_RESPONSE_TIMEOUT_S=3600 .venv/bin/python scripts/diagnose_gait_timeout.py
```

`IHM_NATIVE_RESPONSE_TIMEOUT_S` exists for that third line only. It is a
diagnostic knob, not a correctness parameter.
