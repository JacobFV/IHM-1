# Unified live body controls and monitors

The source UI's default **Start Body** now uses `/api/embodied/sessions`. There is no automatic fallback to `/api/scene`: the older constrained mechanics experiment is available only through Advanced execution. Environment choices map Free space→`free`, Bed→`supine`, Floor→`upright`. Native integration acceptance and production build remain separate from these interface tests.

Live body frames drive canonical `entities`, including their current rotations/deformation gradients. Respiratory skin displacement remains applied in canonical reference coordinates before each entity transform; diagnostic `respiration.skin_field.entity_transforms` is not applied again. Historical physiology and regional replay remain separate owners and are labeled as recorded materializations.

## Real controls and display ownership

- **Live body** displays actual native values and units supplied in `physiology.signal_metadata`.
- **Live signal 1–3** each select an independent native variable. A bounded 600-sample history retains actual accepted timestamps and sequence numbers. Null/nonfinite observations are displayed as unavailable; no synthetic interpolation or invented unit is used. Histories clear when the actor identity changes.
- **Motor & skin inputs** discovers actual `mechanics.muscles` keys. Requested descending drive, sensory block and motor block remain individually addressed. Observed native activation/excitation/fiber length/tendon force are displayed separately from requested drive. Native muscle activation dynamics remain the mechanical plant's responsibility.
- **Whole-Skin compression** is explicitly a uniform native compartment input in Pa, not pressure beneath a cursor. Release clears requested drives, blocks and pressure; existing muscle activation still follows its native decay law.
- Cursor force uses the selected entity's current rigid and affine transform for its material grab offset. It sends `forces:[{id,point_m,force_n}]` to the unified runtime. It is not a visual-only drag.

The default right stack contains current-body cards. Recorded signals, playback, experiment controls and previous inspector/evidence widgets remain independently addable from New pane; saved workspace layouts retain their chosen widgets. Three live signal cards are available in this version, rather than an unbounded number of subscriptions.

The frame scope reports source registration, controller and metabolic limitations. The native support law currently uses proxy contacts; there is not yet a validated single canonical-frame bed/floor mesh for the collection of segment registrations. The unified view therefore does not draw the reduced experiment's guessed support box as if it represented those contacts.

## Session lifecycle

Create immediately returns a retained actor ID with `initializing` or `ready` status. The viewer polls GET until a full `ihm.embodied-frame.v1` exists. Pending initialization disables recorded playback ownership. Reset requests native cleanup and polls until `closed:true`; source/region/history handoffs wait for that confirmation. Errors preserve the last accepted state and show stopped/paused status.

An uncertain step is never retried: GET retrieves the committed sequence and the viewer pauses. After a lost create response, only an actor newly allocated since that attempt can be recovered automatically; a definite 4xx rejection does not adopt another actor. Reconnect existing body is an explicit read/attach action and leaves it paused. Concurrent controllers are still discouraged and backend sequence checks remain authoritative.

Changing source, starting a separate protocol or selecting a regional/systemic replay first awaits shared reset. Closing a pending actor does not free the resource slot until backend cleanup completes. No simulation occurs merely from opening the viewer.

## Resource budget and tests

Existing 30 FPS, one device-pixel ratio, hidden-tab pause, two concurrent geometry fetches and opt-in hair dynamics remain intact. Live material projection and monitor redraws are capped at 5 Hz; the backend's 20 ms exchange time is unchanged. Samples are taken from actual accepted states, with no wall-time catch-up. Skin normals are restored once when leaving a respiratory field, rather than recomputed forever after the field is gone.

Eleven focused CPU/DOM tests cover endpoint ownership, bounded histories, independent real monitor event handlers, selective drive/block/pressure/release, material-offset rotation/affine deformation, uncertain command recovery and pending cleanup. The lightweight DOM test uses no browser or graphics. The intercepted browser fixture has been updated for the pending embodied API but has **not** been executed during this resource-constrained work.

```
cd app
prlimit --as=1073741824 -- nice -n 10 node --jitless --max-old-space-size=128 \
  --test --test-concurrency=1 test/embodied-live.test.js test/embodied-panels.test.js \
  test/scene-interaction.test.js test/monitors.test.js
```

Observed peak RSS for these small test groups was about 50 MiB. No browser, native job or user-tab reload was launched. The production build must follow the root agent's native acceptance checkpoint; source UI checks alone do not establish integrated native behavior or successful default publication.
