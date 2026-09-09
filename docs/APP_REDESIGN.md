# The app: one live brain, one prompt bar, one way to look inside

Three problems, and the third is the one that needs thought.

## 1. "Start simulation" vs "play" is two verbs for two different things

It reads as one choice with two buttons. It is not. Today:

- **play** replays a *recorded* trajectory. Its own tooltip says so: "This does
  not run the simulation."
- the live path posts to `/api/embodied/sessions` and steps the real body.

So the app offers a recording and a simulation through affordances that look
interchangeable, and the tooltip is doing work the interface should do.

**The fix is to stop treating them as the same control.** The simulation is the
app. It is always live, and the transport is a single play/pause that runs or
holds *it*. A recording is a different noun — you *open* one, and when you do,
the transport scrubs the recording and the header says which one, in the way a
video player is obviously not a camera. No mode toggle, no shared button.

## 2. The prompt bar

The brain takes images, audio and text, and returns images, audio and text. No
integrations, no tools, no retrieval from outside the corpus. Six paths, and each
is a materialization that already exists or is a trivial composition of two:

| in | out | how | measured |
|---|---|---|---|
| image | image | drive the optic nerve, embed, retrieve nearest in the EEG bank | **63.5%** top-1, 200-way |
| image | text | the same retrieval, read the concept name off the matched entry | identical to the above, by construction |
| text | image | look the concept up in the vocabulary, drive with its image | exact |
| audio | audio | cochleagram in, embed through the cochlear nerve, retrieve nearest cochleagram, resynthesise | audio_meg head, ~6.5% vs a 7.12% ceiling |
| audio | text | same embedding, nearest labelled entry | as above |
| text | audio | vocabulary lookup to a cochleagram, resynthesise | exact |

**What must be said in the UI, not buried:** every "output" is a *retrieval from
the corpus*, not a generation. The brain does not draw an image; it finds the one
its cortical state is nearest to. Showing a returned image without that label
would be the single most misleading thing this app could do. So the panel shows
the top-5 with cosine similarities and the chance level, the way
`scripts/introspect.py` already does — and it shows them for the same reason:
the top-1 similarity is **not** a confidence. It falls only 0.706 → 0.686 under a
scramble that costs 63 points of accuracy. A wrong answer looks exactly as
certain as a right one, and the interface has to make that visible rather than
hide it behind a single confident-looking result.

## 3. Looking inside: the ring, and what a click should mean

The IBM-1 site rings the brain with its materializations and draws leader lines
to the structures each one touches. The instinct to reuse it here is right, but
the app has something the site does not: the signals are *live*, and there is a
body on the other end. So the ring should not be a menu of models. It should be
**a ring of the things currently exchanging signal**, and clicking one should
answer a question you cannot otherwise ask: *what is this contributing, right
now, to the thing I am looking at?*

That question has a real answer here because the routes are declared, not drawn:

- `ibm/topologies/ihm_bridge.py` joins 142 routes (71 nerves x 2 sides), each
  carrying its measured path length from the body mesh and a **delay per fibre
  class** — optic 2.5/4.2/8.3 ms across three retinal populations, vagus
  6.4-350 ms across five.
- `ibm/embodiment.py` declares 483 ports.
- `ibm/processes/cord.py` maps 81 of 96 muscles onto 20 of 31 spinal segments.

So a click can trace an actual path, and the annotation can carry a real number.

**The interaction.** One thing is *highlighted* — a system: the cortical sheet, a
lobe, the cord, a limb, the viscera. The ring shows every materialization and
route that touches it. Selecting one dims the rest and draws its contribution:
the ports it writes, the fibre classes it writes them on, the delay on each, and
a live magnitude. Selecting the *system* instead of a contributor shows the
converse — everything feeding it, ranked by how much it is currently carrying.

**The honest part, which is also the interesting part.** Most of those arrows
are currently near zero, and the interface should show that rather than animate
a signal that is not flowing. Measured: only **0.03-0.07%** of a sensory drive
reaches the disjoint motor region, and severing the cortical kernel changes the
motor command by nothing. An interface that drew a confident glowing arrow from
sensory cortex to motor cortex would be asserting something this project has
measured to be false. The width of the arrow should *be* the measured magnitude,
so a dead path looks dead — and the one genuinely load-bearing result, that
severing the cortical dynamics drops the body at 2.90 s while the intact one
holds at 0.17 mm, should be visible as the thick arrow it actually is.

That is the design's whole argument: the ring is not decoration, it is the
ablation made continuous. You are looking at which parts are carrying the
system, live, with the arrows sized by measurement.

## What that needs from the backend

- `GET /api/brain/graph` — the static structure: materializations, ports, routes,
  fibre classes, delays. Derived from the declarations, cached.
- `GET /api/brain/state` — per-tick magnitudes on each edge, so the arrows are
  live rather than illustrative.
- `POST /api/brain/prompt` — {modality, payload} in, {modality, results[], chance}
  out, where results carry similarity and the corpus entry they came from.
