import test from 'node:test';
import assert from 'node:assert/strict';
import {transportState, recordingCatalog, recordingHeader, SIMULATION, RECORDING} from '../src/transport.js';

test('the transport is the simulation until a recording is explicitly opened', () => {
 const stopped = transportState({});
 assert.equal(stopped.owner, SIMULATION);
 assert.equal(stopped.label, 'Simulation');
 assert.equal(stopped.action, 'start');
 assert.equal(stopped.scrub, false);
 assert.equal(stopped.speed, false);
 assert.equal(stopped.disabled, false);
 // The old control said "This does not run the simulation" in a tooltip. This
 // one runs the simulation, and nothing else claims to.
 assert.match(stopped.title, /Run the simulation/);
 assert.doesNotMatch(stopped.title, /does not run/i);

 const running = transportState({started: true, running: true});
 assert.equal(running.action, 'hold');
 assert.equal(running.symbol, '❚❚');
 const held = transportState({started: true, running: false});
 assert.equal(held.action, 'resume');
 assert.equal(held.symbol, '▶');
});

test('an open recording takes the transport and says so, without a mode toggle', () => {
 const open = transportState({open: true, recording: {label: 'Canonical trajectory'}, frames: 120, index: 4});
 assert.equal(open.owner, RECORDING);
 assert.equal(open.label, 'Recording');
 assert.equal(open.name, 'Canonical trajectory');
 assert.equal(open.scrub, true);
 assert.equal(open.speed, true);
 assert.equal(open.index, 4);
 assert.equal(open.disabled, false);
 // Running the simulation is not reachable from the recording's transport.
 assert.doesNotMatch(open.title, /Run the simulation/);
 // The index is clamped to the recording, never to a frame it does not hold.
 assert.equal(transportState({open: true, frames: 3, index: 99}).index, 2);
 assert.equal(transportState({open: true, frames: 0, index: 5}).index, 0);
});

test('a dead transport says why it is dead rather than looking broken', () => {
 const stale = transportState({open: true, frames: 0, recordingError: 'Canonical trajectory is stale'});
 assert.equal(stale.disabled, true);
 assert.equal(stale.title, 'Canonical trajectory is stale');
 const thin = transportState({open: true, frames: 1});
 assert.equal(thin.disabled, true);
 assert.match(thin.title, /fewer than two frames/);
 const volume = transportState({staticVolume: true});
 assert.equal(volume.owner, SIMULATION);
 assert.equal(volume.disabled, true);
 assert.match(volume.title, /static/);
});

test('only completed canonical runs are offered as recordings', () => {
 const entries = recordingCatalog({
  canonical: {available: true, detail: '120 computed frames'},
  scenarios: {runs: [
   {id: 'a', status: 'completed', canonical_body: true, summary: {canonical_body: {frames: 40}}},
   {id: 'b', status: 'failed', canonical_body: true},
   {id: 'c', status: 'completed', canonical_body: false},
   {id: 'd', status: 'running', canonical_body: true},
  ]},
 });
 assert.deepEqual(entries.map(e => e.id), ['canonical', 'a']);
 assert.equal(entries[1].detail, '40 computed frames');
 assert.match(entries[1].url, /run=a$/);
 // Nothing to open is not a recording with an empty name.
 assert.deepEqual(recordingCatalog({canonical: {available: false}}), []);
});

test('the header names the open recording and nothing else', () => {
 assert.equal(recordingHeader(null, {frames: 5}), null);
 const header = recordingHeader({label: 'Scenario 2026', detail: '40 computed frames'},
  {frames: 40, index: 3, time_s: 1.25});
 assert.equal(header.name, 'Scenario 2026');
 assert.match(header.detail, /frame 4 \/ 40/);
 assert.match(header.detail, /t 1\.250 s/);
 assert.equal(header.failed, false);
 // A recording that would not load still gets its header, saying why rather
 // than reciting a frame count the file does not deliver.
 const broken = recordingHeader({label: 'Scenario 2026', detail: '40 computed frames'},
  {frames: 0, error: 'Canonical trajectory is stale; rematerialize the native run'});
 assert.equal(broken.failed, true);
 assert.equal(broken.detail, 'Canonical trajectory is stale; rematerialize the native run');
 assert.doesNotMatch(broken.detail, /40 computed frames/);
});
