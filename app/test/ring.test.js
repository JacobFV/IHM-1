import test from 'node:test';
import assert from 'node:assert/strict';
import {arrowWidth, magnitudeText, contributors, fibreRows, placeLabels, HAIRLINE_PX, FULL_PX} from '../src/ring.js';

const graph = {
 systems: [{id: 'precentral', label: 'Precentral (motor)'}, {id: 'cord', label: 'Spinal cord'}],
 edges: [
  {id: 'cortex.postcentral_to_precentral', kind: 'materialization', label: 'Sensory sites → motor sites',
   system: 'postcentral', target: 'precentral', magnitude: 0.0005, magnitude_kind: 'measured'},
  {id: 'cortex.stance_correction', kind: 'materialization', label: 'Cortical sheet → 98 muscle commands',
   system: 'cortical_sheet', target: 'cord', magnitude: 1, magnitude_kind: 'measured'},
  {id: 'cord.arc.stretch', kind: 'arc', label: 'Monosynaptic stretch', system: 'cord', target: 'cord',
   magnitude: null, magnitude_kind: 'declared', declared_gain: 0.4, loop_delay_s: 0.03},
  {id: 'peripheral-nerve-left-median', kind: 'route', label: 'left median nerve', system: 'upper_limb',
   target: 'cord', magnitude: null, magnitude_kind: 'unmeasured', side: 'left',
   fibres: [{fibre_class: 'ia', velocity_m_s: 100, delay_s: 0.007},
            {fibre_class: 'c', velocity_m_s: 1, delay_s: 0.7}]},
  {id: 'peripheral-nerve-right-median', kind: 'route', label: 'right median nerve', system: 'upper_limb',
   target: 'cord', magnitude: null, magnitude_kind: 'unmeasured', side: 'right',
   fibres: [{fibre_class: 'ia', velocity_m_s: 100, delay_s: 0.0072}]},
 ],
};

test('arrow width IS the measurement, linearly, so a dead path looks dead', () => {
 // 0.03-0.07% of a sensory drive reaches the disjoint motor region. Drawn at
 // that fraction it is a hairline, which is the point: an arrow glowing from
 // sensory to motor cortex would assert something measured to be false.
 const dead = arrowWidth(0.0005, 'measured');
 const alive = arrowWidth(1, 'measured');
 assert.equal(alive, FULL_PX);
 assert.ok(dead < HAIRLINE_PX + 0.02, `a measured 0.05% must stay a hairline, got ${dead}`);
 assert.ok(dead / alive < 0.05);
 // Linear, not logarithmic: half the magnitude is half the width above the
 // hairline floor. A log scale would draw 0.0005 at a third of full width.
 const half = arrowWidth(0.5, 'measured');
 assert.ok(Math.abs((half - HAIRLINE_PX) - (alive - HAIRLINE_PX) / 2) < 1e-9);
 // Declared parameters and unmeasured edges make no width claim at all.
 assert.equal(arrowWidth(null, 'declared'), HAIRLINE_PX);
 assert.equal(arrowWidth(0.9, 'unmeasured'), HAIRLINE_PX);
 assert.equal(arrowWidth(NaN, 'measured'), HAIRLINE_PX);
});

test('an unmeasured edge is never printed as a measured zero', () => {
 assert.equal(magnitudeText({magnitude: null, magnitude_kind: 'unmeasured'}), 'not measured');
 assert.match(magnitudeText({magnitude: 0.0005, magnitude_kind: 'measured'}), /0\.050%.*measured/);
 assert.match(magnitudeText({magnitude: 1, magnitude_kind: 'measured'}), /100\.0%/);
 assert.match(magnitudeText({magnitude: 0.4, magnitude_kind: 'live'}), /live this tick/);
 assert.equal(magnitudeText({magnitude_kind: 'declared', declared_gain: 0.4}), 'gain +0.40 · declared');
});

test('contributors are everything touching the system, ranked by what they carry', () => {
 const into = contributors(graph, null, 'precentral');
 assert.deepEqual(into.map(e => e.id), ['cortex.postcentral_to_precentral']);
 assert.equal(into[0].direction, 'into');

 const cord = contributors(graph, null, 'cord');
 // The load-bearing one first: severing it drops the body at 2.90 s.
 assert.equal(cord[0].id, 'cortex.stance_correction');
 assert.equal(cord[0].magnitude, 1);
 // Both sides of one nerve are one ring item, with both sets of delays kept.
 const median = cord.find(e => e.key === 'route:median');
 assert.equal(median.members.length, 2);
 assert.equal(median.label, 'median nerve');
 assert.deepEqual(fibreRows(median).map(r => r.fibre_class), ['ia', 'ia', 'c']);
 assert.equal(fibreRows(median)[0].delay_s, 0.007);
 // Declared and unmeasured sort below anything measured, and never above it.
 assert.ok(cord.findIndex(e => e.magnitude_kind === 'measured') <
           cord.findIndex(e => e.magnitude_kind === 'unmeasured'));
});

test('a live magnitude overrides the stored one and says which it is', () => {
 const state = {live: true, edges: {'cortex.stance_correction': {magnitude: 0.37, magnitude_kind: 'live',
   measurement: 'Live motor excitation peak 0.37.'}}};
 const cord = contributors(graph, state, 'cord');
 const stance = cord.find(e => e.id === 'cortex.stance_correction');
 assert.equal(stance.magnitude, 0.37);
 assert.equal(stance.magnitude_kind, 'live');
 // With no state at all nothing is overridden: the measured values stand and
 // are still labelled measured rather than promoted to live.
 const stored = contributors(graph, null, 'cord').find(e => e.id === 'cortex.stance_correction');
 assert.equal(stored.magnitude_kind, 'measured');
});

// --- placement: a label stands where its own point projects ----------------
// The labels used to be dealt onto a screen-space ellipse, so a label's position
// said nothing about where its structure was and only one system's 26 strongest
// edges fitted. These cover the placement that replaced it.

test('a label is placed at its own projected point, hanging outward', () => {
 // Outward, away from the middle of the view: a label placed toward the centre
 // covers the structure it names, and a hundred of them bury the body.
 const [left] = placeLabels([{id: 'a', x: 200, y: 300, w: 120, h: 20}], 1000, 600);
 assert.equal(left.side, 'l');
 assert.ok(left.cx < 200, 'a left-half label hangs further left');
 assert.equal(left.cy, 300, 'and level with its point when nothing is in the way');
 const [right] = placeLabels([{id: 'b', x: 800, y: 300, w: 120, h: 20}], 1000, 600);
 assert.equal(right.side, 'r');
 assert.ok(right.cx > 800, 'a right-half label hangs further right');
});

test('a label never leaves the viewport; it takes another direction or none', () => {
 // Hard against the left edge, outward is off-screen, so it must come back
 // across rather than hang half off the view.
 const [edge] = placeLabels([{id: 'a', x: 5, y: 300, w: 200, h: 20}], 1000, 600);
 assert.ok(!edge.hidden && edge.cx - 100 >= 2, 'the box stays on screen');
 assert.equal(edge.x, 5, 'while the point it belongs to does not move');
 // A viewport too small for the box at any angle drops it rather than clipping.
 const [tight] = placeLabels([{id: 'b', x: 50, y: 50, w: 400, h: 20}], 300, 100);
 assert.equal(tight.hidden, true);
});

test('two labels on the same point are nudged apart, biggest keeps its place', () => {
 const placed = placeLabels([
  {id: 'big', x: 400, y: 300, w: 120, h: 20},
  {id: 'small', x: 400, y: 300, w: 120, h: 20},
 ], 1000, 600);
 assert.equal(placed[0].cy, 300, 'the first, which ranks highest, keeps the point');
 assert.notEqual(placed[1].cy, 300, 'the second moves');
 assert.ok(Math.abs(placed[1].cy - 300) >= 15, 'by at least one step');
 assert.equal(placed[1].x, 400, 'and still leads back to the same point');
});

test('a label with nowhere to go is dropped rather than stacked', () => {
 // Forty labels on one point in a small viewport: the candidate directions run
 // out and the rest are hidden and counted, not piled on top of each other.
 const many = Array.from({length: 40}, (_, i) => ({id: `e${i}`, x: 300, y: 120, w: 120, h: 20}));
 const placed = placeLabels(many, 600, 240);
 assert.ok(placed.some((p) => p.hidden), 'some are dropped');
 // How many fit is set by the directions and radii a label may take before it
 // stops being about its own point, not by how many were asked for. The ones
 // that fit are the ones that came first: the ranking decides who keeps a place.
 const kept = placed.filter((p) => !p.hidden);
 assert.ok(kept.length >= 3, `the strongest are still placed, got ${kept.length}`);
 assert.deepEqual(kept.map((p) => p.id), many.slice(0, kept.length).map((p) => p.id));
});

test('placement is pure: the same inputs give the same answer', () => {
 const input = [{id: 'a', x: 300, y: 200, w: 100, h: 18}, {id: 'b', x: 305, y: 205, w: 100, h: 18}];
 assert.deepEqual(placeLabels(input, 800, 600), placeLabels(input, 800, 600));
});
