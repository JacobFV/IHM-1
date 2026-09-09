# the anatomy and the segments the controller drives

`docs/DISCONNECTS.md` item 1, in IBM-1:

| | |
|---|---|
| declared / rendered | 3,816 BodyParts3D entities in the simulated body |
| actually driven | 22 OpenSim bodies, 80 muscles |
| mapping between them | **none** |

There is now a mapping, and everything below is measured on it.

    scripts/bind_anatomy_to_segments.py   -> data/derived/anatomy-segment-binding/binding.json
                                          -> data/derived/anatomy-segment-binding/report.json
    scripts/render_anatomical_motion.py   -> a render, an app-format trajectory, a motion report

`data/derived/` is gitignored. The scripts are the artifact; regenerating the
binding takes about a minute.

## What the binding is

Each of the 4,000 atlas entities is assigned one of the 22 OpenSim bodies.

334 of them are bones and are assigned by the atlas' own label — "right femur"
is `femur_r`, every rib and vertebra and the whole skull is `torso`, because
`torso` is the model's HAT body and carries `hat_skull`, `hat_jaw`,
`hat_spine` and `hat_ribs_scap`. Those 334 bones then *are* the reference
shape of each segment.

The other 3,666 are assigned geometrically: each entity's surface is sampled,
every sample takes the segment whose bone group is nearest, and the entity
takes the majority. The vote share is kept as `coherence`. An entity that
straddles a joint has a low one, and a rigid binding will tear it — 220
entities are below 0.6 and the renderer does not pretend otherwise.

Five entities are **not** bound and are excluded by name and reason: the skin,
its three layers, and the lymphatic network graph. No single rigid segment
carries a whole-body surface. The skin already has a continuous
graph-regularised linear-blend attachment over the same 22 segments
(`data/derived/canonical/continuous_surface_binding.json.gz`) and the renderer
uses it, so the integument stretches across a joint instead of tearing at it.

**3,995 of 4,000 bound, 99.88%.**

## The registration, and what it costs

The OpenSim model is a scaled Rajagopal subject; the atlas is the BodyParts3D
adult male; and the model's default pose is not the atlas pose. So the
reference pose (33 coordinates) and one similarity (scale, rotation,
translation) are fitted together by least squares against the 22 bone groups.
Scale 0.963. What is left over:

| | |
|---|---|
| bone-group centroid residual | **24.7 mm RMS**, 51.8 mm max (humerus_r) |
| joint centre vs the anatomical joint | **38.4 mm median**, 101.4 mm max (shoulder) |

Neither is measurement error and neither is subject variation. It is the
residual of putting one specimen's skeleton onto another's.

The joint-centre offset is the number that matters at playback, because the
anatomy is rotated about those points. Its consequence is measured directly:
at the reference pose the 64 closest surface point pairs across a joint sit a
fixed distance apart; bind the two sides to different rigid segments and that
distance changes. The largest change over a clip is exactly the gap or
interpenetration the binding introduces.

## Gates

None of these is the assignment rule restated.

| gate | result |
|---|---|
| held-out bones — drop a bone from its own segment's reference shape and re-assign it | **317 / 322, 1.6% misassigned** |
| native muscle attachments — the bound segment must lie on the kinematic chain between the muscle's declared attachment bodies | **78 / 80, 2.5% misassigned** |
| hand-checked cases with an answer anatomy fixes | **21 / 21** |
| second opinion — same rule, run against the registered OpenSim bone meshes instead of the atlas' own bones | **93.5% agreement** |

Only 322 of the 334 bones can be held out; the other twelve are the only bone
in their segment, so removing one removes the segment. All five held-out
failures are bones that are more than 30% of their own segment's surface —
calcaneus, both hip bones, right tibia — where the same thing happens in
miniature.

Endpoint-only scoring gives the muscle gate 64/80, but that is the wrong test:
a two-joint muscle's belly lies on the bone *between* its attachments, and
`attachment_bodies` names the endpoints. Hamstrings bind to the femur,
gastrocnemius to the tibia; those are right, not wrong.

Both remaining misses are **psoas major**, and they are a real disagreement
rather than a binding error: OpenSim attaches psoas to the pelvis, and the
atlas originates it on the lumbar spine, which is `torso`. The anatomy is
right and the model's attachment list is the approximation.

## How far the anatomy actually moves

The comparison that matters is the body's own canonical run, which moves
**6.35 mm** at most across 30 s. It is perfusing and breathing, not moving.

| trajectory | duration | max entity travel | median | worst joint opening | median opening |
|---|---|---|---|---|---|
| `unified-world-ofb2dp7z/full` — trained cortical kernel, in the loop | 4.98 s | 262 mm | 100 mm | 44 mm | 5 mm |
| `gait-best` — the gait controller: one genuine step, then a fall | 1.49 s | 1,028 mm | 397 mm | 73 mm | 6 mm |
| `forced-gait-corpus` — measured human gait, played in place | 2.39 s | 1,193 mm | 1,001 mm | 69 mm | 15 mm |
| `unified-world-ofb2dp7z/sever` — cortex severed, falls | 2.88 s | 1,835 mm | 949 mm | 142 mm | 22 mm |

**Travel is not skill.** The severed run travels furthest of all because it
collapses, and the falling gait controller beats the standing cortical one.
The only thing this column establishes is that the anatomy is now driven at
all, which is what item 1 says it was not. Skill belongs to
`data/derived/gait-best/report.json` and to the ablations, not here.

## What this does not do

- One rigid segment per entity. There is no soft-tissue deformation, no
  volume preservation, no sliding, and the joint openings above are the price.
- The model has no neck, no shoulder girdle and no spine joint. Skull, jaw,
  teeth, every vertebra, every rib, both scapulae and both clavicles ride
  `torso` as one body, because that is what the model is. The head cannot nod.
- A replay carries no physiology and no internal state. `/api/body/segment-bound`
  lists these separately from the canonical trajectory and says so on every
  entry, so a kinematic replay can never read as the body's own physics.
