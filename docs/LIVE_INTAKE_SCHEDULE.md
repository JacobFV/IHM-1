# Live intake UI contract

The scheduler and native ownership contract is documented in
[research/LIVE_INTAKE_SCHEDULE.md](research/LIVE_INTAKE_SCHEDULE.md).
`IntakeSchedule.add_events(events, current_tick)` validates an entire batch
before mutation, retains all IDs including completed events, and permits current
or future tick scheduling within its horizon. Pending or uncertain commands
prevent changes. `snapshot()` exposes queued, issued, accepted and uncertain
states as `{events: [...]}`. Snapshots are observations, not resumable checkpoints.

## Component integration

`app/src/intake-monitor.js` exports:

```js
const monitor = mountIntakeMonitor(host, {
  submit: async events => {
    // Use the owning controller's serialized queue request.
    // Return the authoritative IntakeSchedule.snapshot() object.
    return submitToCurrentOwner(events);
  },
});
monitor.update(authoritativeIntakeSnapshot);
monitor.setConnected(true);
// Before changing owner or replacing its DOM:
monitor.dispose();
```

The injected callback receives an array containing one event:
`{event_id, time_s, meal: {name, carbohydrate_g, protein_g, fat_g, sodium_g,
calcium_mg, water_ml}}`. The form validates the native Meal fields and 20 ms
schedule time alignment. All values are numeric, finite and bounded by the
native API. At least one amount must be positive. The server must validate again,
including owner clock, duplicate IDs, horizon and capacity.

Snapshot entries contain that event payload plus `state` and optional
`issued_tick`, `native_sequence`, `reason`. The UI reports queued, issued,
accepted or uncertain only from authoritative snapshots. A request in flight
is labeled as submitting; it does not fabricate an issued native command.
Acknowledgment does not imply absorption. Snapshot identity and state cannot
regress or disappear within one owner. The component starts disconnected and
requires explicit `setConnected(true)` from the controller.

Unknown callback failures lock further submission without a retry control.
The controller may perform read-only recovery after a lost response and return
an authoritative snapshot containing every submitted event. It must not repeat
the POST. For a proven nonmutating validation rejection only, it may throw an
Error with `definitelyRejected = true`; this lets the user correct the request.
A generic HTTP error alone does not establish rejection before mutation.
Any reported uncertain native event permanently locks the component instance.
Reconnecting does not reset uncertainty. Dispose and remount for a genuinely
new owner; do not remount to bypass uncertainty on the same native session.
`dispose()` disables old controls and ignores late asynchronous responses.

## Suggested monitor registration

The shared monitor owner should add a catalog item with ID `intake`, title
`Food & water`, tags `['live','controls','physiology']`, and a host div with ID
`intake-monitor` in its existing widget registry. The live controller owns the
mount call, serialized submit callback, owner changes and snapshot forwarding.
The component does not edit those shared modules or connect itself to an API.

Verification: `cd app && node --test test/intake-monitor.test.js` covers the
actual DOM handlers, double-submit prevention, native field validation, receipt
state rendering, terminal uncertainty, definite rejection, and late responses
after disposal. These bounded DOM tests launch no browser or native engine and
do not establish that a live integration is connected.
