import test from "node:test";
import assert from "node:assert/strict";
import { filterMonitors } from "../src/monitors.js";
test("monitor search combines words with overlapping category membership", () => {
  assert.deepEqual(
    filterMonitors("  finite LAPLACE ").map((m) => m.id),
    ["temporal-spectrum", "signals"],
  );
  assert.deepEqual(
    filterMonitors("contact", "mechanics").map((m) => m.id),
    ["contact"],
  );
  assert.deepEqual(
    filterMonitors("contact", "evidence").map((m) => m.id),
    ["contact"],
  );
  assert.deepEqual(filterMonitors("contact", "vascular"), []);
  assert.deepEqual(filterMonitors("invented medicine"), []);
});
test('food and drink schedule is discoverable in the live monitor catalog',()=>{
 assert.deepEqual(filterMonitors('food','live').map(m=>m.id),['intake']);
 assert.deepEqual(filterMonitors('drink','controls').map(m=>m.id),['intake']);
});

test('live finite-window spectrum has its own searchable monitor card',()=>{
 assert.deepEqual(filterMonitors('laplace','live').map(m=>m.id),['temporal-spectrum']);
 assert.deepEqual(filterMonitors('damping','spectra').map(m=>m.id),['temporal-spectrum']);
});
test('local microvascular geometry is discoverable as a source-conditioned anatomy monitor',()=>{
 assert.deepEqual(filterMonitors('local','vascular').map(m=>m.id),['microvascular']);
});
