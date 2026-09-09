# sustained prone locomotion, and what it is allowed to claim

`scripts/crawl.py`, search B (`data/derived/crawl-search/searchB.json`), best
parameters replayed at a 16 s horizon into `data/derived/crawl-best/`.

## The number

| | |
|---|---|
| simulated | **16.0 s**, full horizon, `stop_reason: horizon` |
| pelvis forward travel | **0.920 m** |
| mean forward velocity | 0.0575 m/s |
| lateral excursion | 0.244 m |
| peak fibre velocity | **9.69 of 10** optimal fibre lengths per second |
| time outside the muscle force–velocity domain | **0.0%** |
| worst excursion past a declared joint range | **0.200 rad** (`ankle_angle_r`) |
| frames | 1600 at 10 ms |
| wall clock | 306 s; slowest single advance 0.68 s |

Fifteen contact elements carry load over the run — torso, pelvis, both femurs,
both patellae, both tibias, both hands, one radius, and the foot contacts — for a
mean total contact force of 854 N against a body weight of 761 N.

## What drives it

Muscle excitation in [0,1] on the 98 source muscles, plus the source model's
declared `CoordinateActuator` torque ports for the lumbar and both arms. **The
torque ports are not muscles.** The model carries no shoulder musculature, so
this motion is not wholly muscle-driven and must not be reported as such.

No coordinate is prescribed, no external force is applied to any body, no motion
constraint exists. Joint values are integrated native OpenSim/Simbody output.

The controller is open loop — a fixed-period pattern, no feedback. A prone body
cannot fall, so the balance regulator that bipedal gait needs is simply absent,
and whatever travel appears is produced by the muscles and the floor rather than
by a controller correcting the plant.

## What it may not claim

**The surfaces are not skin.** The body slides on the engine's
inertia-inscribed proxies: one sphere per segment at its centre of mass, 0.258 m
at the chest, 0.086 m at mid-thigh. A crawl on four balls is supported and
propelled by real contact forces, but it is not a chest, two knees and two
forearms on a floor.

**0.200 rad is not zero.** The plant is still 11.5° outside its declared ankle
range, held there by the opt-in joint stops
(`docs/NATIVE_JOINT_LIMITS.md`). For scale: a body lying still, driven at 0.02
tonic excitation, already sits 0.240 rad outside its declared range at
`hip_rotation_l`, because the source model's own passive stop is centred 0.222 rad
wider than the range it declares. This run is *inside* that resting figure. It is
not inside the declared ranges, and no result here should say it is.

**It is a crawl, not a gait.** Bipedal walking on this plant still gets one
genuine step and then falls at 1.49 s (`docs/NATIVE_INTEGRATOR_STALL.md`).

## Against the withdrawn result

The first prone search returned 973 mm over 16 s and was withdrawn: it had folded
the ankles to 145° of plantarflexion, 1.651 rad past their declared range, and
paddled. With the stops in place the same search returns 920 mm — **95% of the
travel, at 12% of the range violation** — and, unlike the withdrawn run, spends no
time at all outside the muscle model's force–velocity domain.

The stops did not cost the result. They cost 5%.
